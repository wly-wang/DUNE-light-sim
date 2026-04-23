#include "Rtypes.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLine.h"
#include "TMath.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "Math/SpecFuncMathMore.h"
#include "TLatex.h"
#include "TStyle.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "functions.h"

double solid(const acc& out, const TVector3 &v);
double solid_lateral(const acc& out, const TVector3 &v);
double omega(const double &a, const double &b, const double &d);
double interpolate(const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate);
Double_t GaisserHillas(double x, double *par);
Double_t Omega_Dome_Model(double distance, double theta);

// ============================================================
// Shared constants / configuration
// ============================================================
static const double pi = 3.1416;
static const int N = 9;                     // theta bins
static const int M = 8;                     // dT bins
static const double d_min = 0.0;
static const double d_max = 1000.0;
static const double step_d = 50.0;
static const double FIT_MAX = 1000.0;
static const double BORDER_FIT_MAX = 600.0;
static const double range_d = 600.0;
static const double range_d_array[M+1] = {0.,75.,150.,225.,300.,375.,450.,525.,600.};
static const double delta_angulo = 90.0 / N;
static const double L_abs = 2000.0;
static const double arapuca_w = 10.0;
static const double arapuca_h = 47.75;
static const double radius = 8*2.54/2.;
static const int OpDetType = 1;             // 0=disk, 1=rectangular, 2=dome
static const bool IsRectangular = true;
static const bool IsSphere = false;
static const bool fVerticalBorderCorrectionMode = true;
static const bool fUseQuadraticNmax = false;
static const bool isDouble = true;

static const double centerY = 0.0;         // fit-side vertical border center
static const double y_foils = 0.0;         // validation-side equivalent center
static const double z_foils = 2904.0;      // only used if radial mode is enabled
static const double x_foils = 0.0;

static const double opdet_x_target = 0.0;
static const double opdet_x_window = 1.0;
static const double z_min_keep = 1000.0;
static const double z_max_keep = 4700.0;
static const int kMaxDevices = 6000;


// ============================================================
// Fit results to be filled at runtime and then reused by validation
// ============================================================
std::vector<double> g_theta_centers;
std::vector<double> g_fit_p1;
std::vector<double> g_fit_p2;
std::vector<double> g_fit_p3;
std::vector<double> g_fit_p4;
std::vector<double> g_slopes1;
std::vector<double> g_slopes2;
std::vector<double> g_slopes3;

// ============================================================
// Output containers from raw pair building
// ============================================================
struct PairData {
  std::vector<double> v_distance;
  std::vector<double> v_hits;
  std::vector<double> v_rec_hits;
  std::vector<double> v_offset_angle;
  std::vector<double> v_d_center;
  std::vector<double> v_x;
  std::vector<double> v_devx;
};

// ============================================================
// Shared detector map loading
// ============================================================
struct DetMap {
  std::vector<double> id;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
};

DetMap LoadDetectorMap(const std::string& positions, bool restrict_x = false) {
  DetMap dm;
  std::ifstream in(positions.c_str());
  if (!in) {
    std::cerr << "WARNING: failed to open detector map: " << positions << std::endl;
    return dm;
  }

  std::string line;
  while (std::getline(in, line)) {
    double id, x, y, z;
    if (sscanf(line.c_str(), "%lf %lf %lf %lf", &id, &x, &y, &z) != 4) continue;
    if (restrict_x && std::fabs(x - opdet_x_target) > opdet_x_window) continue;
    dm.id.push_back(id);
    dm.x.push_back(x);
    dm.y.push_back(y);
    dm.z.push_back(z);
  }
  return dm;
}

// ============================================================
// Raw pair builder used by fit stage
// ============================================================
PairData BuildPairData(const std::string& positions, const std::string& input_file) {
  PairData out;

  DetMap dm = LoadDetectorMap(positions, false);
  std::cout << "Loaded detector positions: " << dm.x.size() << std::endl;

  TFile* f = new TFile(input_file.c_str());
  TTree* tree = (TTree*)f->Get("myTree");
  if (!tree) {
    std::cerr << "ERROR: tree 'myTree' not found in " << input_file << std::endl;
    return out;
  }

  int VUV_hits[kMaxDevices];
  int Vis_hits[kMaxDevices];
  double X, Y, Z;
  int genPhotons, numberDevices;
  tree->SetBranchAddress("numberDevices", &numberDevices);
  tree->SetBranchAddress("X", &X);
  tree->SetBranchAddress("Y", &Y);
  tree->SetBranchAddress("Z", &Z);
  tree->SetBranchAddress("VUV_hits", VUV_hits);
  tree->SetBranchAddress("Vis_hits", Vis_hits);
  tree->SetBranchAddress("genPhotons", &genPhotons);

  for (int n = 0; n < tree->GetEntries(); ++n) {
    tree->GetEntry(n);

    if (numberDevices > kMaxDevices) continue;
    if (numberDevices > (int)dm.x.size()) continue;

    double posSource[3] = {X, Y, Z};
    if (posSource[2] < z_min_keep || posSource[2] > z_max_keep) continue;

    double dT = std::abs(posSource[1] - centerY);
    if (dT <= 0) continue;

    for (int i = 0; i < numberDevices; ++i) {
      if (std::fabs(dm.x.at(i) - opdet_x_target) > opdet_x_window) continue;

      int entries = VUV_hits[i];
      double distance_to_pmt = std::sqrt(
        std::pow(posSource[0] - dm.x.at(i), 2) +
        std::pow(posSource[1] - dm.y.at(i), 2) +
        std::pow(posSource[2] - dm.z.at(i), 2)
      );
      if (distance_to_pmt <= 0) continue;

      double coseno = std::sqrt(std::pow(posSource[0] - dm.x.at(i), 2)) / distance_to_pmt;
      double theta = std::acos(coseno) * 180.0 / pi;
      if (theta >= 90 || theta < 0) continue;

      double geo_factor = -1.0;
      if (IsRectangular) {
        acc detPoint;
        detPoint.ax = dm.x.at(i);
        detPoint.ay = dm.y.at(i);
        detPoint.az = dm.z.at(i);
        detPoint.w = arapuca_w;
        detPoint.h = arapuca_h;
        TVector3 ScintPoint(posSource[0], posSource[1], posSource[2]);
        TVector3 OpDetPoint(dm.x.at(i), dm.y.at(i), dm.z.at(i));
        TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;
        geo_factor = solid(detPoint, ScintPoint_rel);
      } else if (IsSphere) {
        geo_factor = Omega_Dome_Model(distance_to_pmt, theta);
      } else {
        double d = std::sqrt(std::pow(posSource[2] - dm.z.at(i), 2) + std::pow(posSource[1] - dm.y.at(i), 2));
        double h = std::sqrt(std::pow(dm.x.at(i) - posSource[0], 2));
        (void)d;
        (void)h;
      }

      if (geo_factor <= 0) continue;
      double rec_N = std::exp(-distance_to_pmt / L_abs) * (genPhotons * geo_factor / (4.0 * TMath::Pi()));
      if (rec_N <= 0) continue;

      out.v_hits.push_back(entries);
      out.v_x.push_back(posSource[0]);
      out.v_devx.push_back(dm.x.at(i));
      out.v_distance.push_back(distance_to_pmt);
      out.v_rec_hits.push_back(rec_N);
      out.v_offset_angle.push_back(theta);
      out.v_d_center.push_back(dT);
    }
  }

  delete f;
  return out;
}

// ============================================================
// Validation predictor uses fit outputs from this same executable
// ============================================================
int VUVHitsPredicted(const int &Nphotons_created,
                     const TVector3 &ScintPoint,
                     const TVector3 &OpDetPoint,
                     const int &optical_detector_type,
                     const double &cosine,
                     const double &theta,
                     const double distance,
                     const int &j) {

  double solid_angle = 0.0;
  if (optical_detector_type == 1) {
    acc detPoint;
    detPoint.ax = OpDetPoint[0];
    detPoint.ay = OpDetPoint[1];
    detPoint.az = OpDetPoint[2];
    detPoint.w = arapuca_w;
    detPoint.h = arapuca_h;
    TVector3 rel = ScintPoint - OpDetPoint;
    solid_angle = solid(detPoint, rel);
  } else if (optical_detector_type == 2) {
    solid_angle = Omega_Dome_Model(distance, theta);
  } else {
    std::cerr << "Unsupported optical detector type in VUVHitsPredicted" << std::endl;
    return 0;
  }

  double hits_geo = std::exp(-distance / L_abs) * (solid_angle / (4 * pi)) * Nphotons_created;

  double r_distance = 0.0;
  if (fVerticalBorderCorrectionMode) {
    r_distance = std::abs(ScintPoint[1] - centerY);
  } else {
    r_distance = std::sqrt(std::pow(ScintPoint[1] - y_foils, 2) + std::pow(ScintPoint[2] - z_foils, 2));
  }

  double pars[4] = {g_fit_p1.at(j), g_fit_p2.at(j), g_fit_p3.at(j), g_fit_p4.at(j)};
  pars[0] += g_slopes1.at(j) * r_distance;
  pars[1] += g_slopes2.at(j) * r_distance;
  pars[2] += g_slopes3.at(j) * r_distance;

  // If your functions.h expects pointer-style GH, use:
  // double xgh = distance;
  double xgh = distance;
  double GH_correction = GaisserHillas(&xgh, pars);

  double hits_rec = GH_correction * hits_geo / cosine;
  return std::round(hits_rec);
}

// ============================================================
// Validation block
// ============================================================
void RunValidation(const std::string& inputfilename, const std::string& positions, const std::string& save_dir) {
  DetMap dm = LoadDetectorMap(positions, true);
  std::cout << "Validation loaded detector positions: " << dm.x.size() << std::endl;
  int numberPMTs = dm.x.size();

  std::vector<double> v_hits_sim;
  std::vector<double> v_hits_geo;
  std::vector<double> v_r;
  std::vector<double> v_prop_dist;

  TFile* f = new TFile(inputfilename.c_str());
  TTree *tree = (TTree *)f->Get("myTree");
  if (!tree) {
    std::cerr << "ERROR: validation tree 'myTree' not found" << std::endl;
    return;
  }

  int VUV_hits[kMaxDevices];
  int Vis_hits[kMaxDevices];
  double X, Y, Z;
  int genPhotons, numberDevices;
  tree->SetBranchAddress("numberDevices", &numberDevices);
  tree->SetBranchAddress("X", &X);
  tree->SetBranchAddress("Y", &Y);
  tree->SetBranchAddress("Z", &Z);
  tree->SetBranchAddress("VUV_hits", VUV_hits);
  tree->SetBranchAddress("Vis_hits", Vis_hits);
  tree->SetBranchAddress("genPhotons", &genPhotons);

  TH2F* h_truePE_vs_predPE = new TH2F("h_truePE_vs_predPE", "", 50, 0, 2.1e6, 50, 0, 2.1e6);

  TH1D* h_tot_truth = new TH1D(
  "h_tot_truth", "", 120, 0, 2.1e6
  );
  TH1D* h_tot_pred = new TH1D(
  "h_tot_pred", "", 120, 0, 2.1e6
  );

  for (int n = 0; n < tree->GetEntries(); ++n) {
    tree->GetEntry(n);
    if (Z < z_min_keep || Z > z_max_keep) continue;

    double posSource[3] = {X, Y, Z};
    double R = fVerticalBorderCorrectionMode
      ? std::abs(posSource[1] - centerY)
      : std::sqrt(std::pow(posSource[1] - y_foils,2) + std::pow(posSource[2] - z_foils,2));

    double total_pe_truth = 0.0;
    double total_pe_prediction = 0.0;

    for (int nPMT = 0; nPMT < numberPMTs; ++nPMT) {
      int detID = static_cast<int>(dm.id.at(nPMT));
      if (detID < 0 || detID >= numberDevices) continue;

      TVector3 ScintPoint(posSource[0], posSource[1], posSource[2]);
      TVector3 OpDetPoint(dm.x.at(nPMT), dm.y.at(nPMT), dm.z.at(nPMT));

      double distance = std::sqrt(
        std::pow(ScintPoint[0] - OpDetPoint[0], 2) +
        std::pow(ScintPoint[1] - OpDetPoint[1], 2) +
        std::pow(ScintPoint[2] - OpDetPoint[2], 2)
      );
      if (distance <= 0) continue;

      double cosine = std::abs(ScintPoint[0] - OpDetPoint[0]) / distance;
      if (cosine <= 0) continue;
      double theta = std::acos(cosine) * 180.0 / pi;
      int j = int(theta / delta_angulo);
      if (j < 0) j = 0;
      if (j >= N) j = N - 1;

      double nPhotons_solid = VUVHitsPredicted(genPhotons, ScintPoint, OpDetPoint, OpDetType, cosine, theta, distance, j);

      total_pe_prediction += nPhotons_solid;
      total_pe_truth += VUV_hits[detID];
      v_hits_sim.push_back(VUV_hits[detID]);
      v_hits_geo.push_back(nPhotons_solid);
      v_r.push_back(R);
      v_prop_dist.push_back(distance);
    }

    h_tot_truth->Fill(total_pe_truth);
    h_tot_pred->Fill(total_pe_prediction);
    std::vector<double> v_total_truth;
    std::vector<double> v_total_pred;
    v_total_truth.push_back(total_pe_truth);
    v_total_pred.push_back(total_pe_prediction);

    h_truePE_vs_predPE->Fill(total_pe_truth, total_pe_prediction);
  }

  std::cout << "\n=== totals between 1.3e6 and 1.8e6 ===\n";
  int n_truth_band = 0, n_pred_band = 0;

  for (size_t i = 0; i < v_total_truth.size(); ++i) {
    if (v_total_truth[i] >= 1.3e6 && v_total_truth[i] <= 1.8e6) {
      n_truth_band++;
      std::cout << "truth total[" << i << "] = " << v_total_truth[i] << "\n";
    }
  }

  for (size_t i = 0; i < v_total_pred.size(); ++i) {
    if (v_total_pred[i] >= 1.3e6 && v_total_pred[i] <= 1.8e6) {
      n_pred_band++;
      std::cout << "pred total[" << i << "] = " << v_total_pred[i] << "\n";
    }
  }

  std::cout << "N truth totals in band = " << n_truth_band << "\n";
  std::cout << "N pred totals in band  = " << n_pred_band << "\n";

  for (size_t i = 0; i < v_total_pred.size(); ++i) {
    if (v_total_pred[i] >= 1.3e6 && v_total_pred[i] <= 1.8e6) {
      n_pred_band++;
    std::cout << "pred total[" << i << "] = " << v_total_pred[i] << "\n";
    }
  }

  std::cout << "N truth totals in band = " << n_truth_band << "\n";
  std::cout << "N pred totals in band  = " << n_pred_band << "\n";

  if (total_pe_truth > 1.8e6) {
  std::cout << "high-total source: X=" << X
            << " Y=" << Y
            << " Z=" << Z
            << " truth=" << total_pe_truth
            << " pred=" << total_pe_prediction
            << "\n";
}


  gSystem->Exec(("mkdir -p " + save_dir).c_str());

  TCanvas* ctest = new TCanvas("c_validation_2d", "", 200, 10, 700, 500);
  gPad->SetLogz();
  h_truePE_vs_predPE->SetTitle(" ");
  h_truePE_vs_predPE->GetXaxis()->SetTitle("total PE truth");
  h_truePE_vs_predPE->GetYaxis()->SetTitle("total PE prediction");
  h_truePE_vs_predPE->SetStats(0);
  h_truePE_vs_predPE->Draw("colz");
  TLine* line1 = new TLine(0,0,2.1e6,2.1e6);
  line1->SetLineColor(kRed);
  line1->Draw("same");
  ctest->SaveAs((save_dir + "/validation_truth_vs_pred.pdf").c_str());

  TCanvas* c_occ = new TCanvas("c_occ", "", 200, 10, 900, 700);
  c_occ->Divide(1,2);
  c_occ->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.12);
  h_tot_truth->SetStats(0);
  h_tot_truth->SetLineWidth(2);
  h_tot_truth->SetTitle("Total-PE occupancy check");
  h_tot_truth->GetXaxis()->SetTitle("total PE truth");
  h_tot_truth->GetYaxis()->SetTitle("events");
  h_tot_truth->Draw("hist");
  c_occ->cd(2);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.12);
  h_tot_pred->SetStats(0);
  h_tot_pred->SetLineColor(kBlue+1);
  h_tot_pred->SetLineWidth(2);
  h_tot_pred->SetTitle("");
  h_tot_pred->GetXaxis()->SetTitle("total PE prediction");
  h_tot_pred->GetYaxis()->SetTitle("events");
  h_tot_pred->Draw("hist");
  c_occ->SaveAs((save_dir + "/validation_truth_pred_occupancy.pdf").c_str());

  const double VALID_MAX = 2000.0;
  // const Double_t edges[] = {
  //   0, 100, 200, 300, 400, 500, 600,
  //   800, 1000, 1200, 1400, 1700, 2000
  // };
  const Double_t edges[] = {
  0, 100, 200, 300, 400,
  500, 600, 800, 1000, 1200, 1400, 1700, 2000
};
  const int n_bins = sizeof(edges)/sizeof(edges[0]) - 1;
  const double DT_CUTOFF = 600.0;

  TProfile* profile = new TProfile("profile_validation", "", n_bins, edges, "s");
  TProfile* profile_core = new TProfile("profile_validation_core", "", n_bins, edges, "s");

  std::cout << "1 source + 1 opdet - Validation total = " << v_prop_dist.size() << std::endl;

  int n_under_1000 = 0;
  int n_under_1500 = 0;
  int n_under_2000 = 0;

  for (size_t i = 0; i < v_prop_dist.size(); ++i) {
    if (v_prop_dist[i] <= 1000.0) n_under_1000++;
    if (v_prop_dist[i] <= 1500.0) n_under_1500++;
    if (v_prop_dist[i] <= 2000.0) n_under_2000++;
  }

  std::cout << "combinations with distance <= 1000 cm: " << n_under_1000 << std::endl;
  std::cout << "combinations with distance <= 1500 cm: " << n_under_1500 << std::endl;
  std::cout << "combinations with distance <= 2000 cm: " << n_under_2000 << std::endl;

  for (size_t i = 0; i < v_prop_dist.size(); ++i) {
    if (v_prop_dist[i] > VALID_MAX) continue;
    if (v_hits_sim[i] <= 0) continue;
    if (v_prop_dist[i] < 1.0) continue;

    double discrepancy = (v_hits_geo[i] - v_hits_sim[i]) / v_hits_sim[i];
    double weight = v_hits_sim[i];
    profile->Fill(v_prop_dist[i], discrepancy, weight);
    if (v_r[i] <= DT_CUTOFF) profile_core->Fill(v_prop_dist[i], discrepancy, weight);
  }

  // ------------------------------------------------------------
  // DEBUG: print everything in the FIRST validation bin
  // ------------------------------------------------------------
  double first_bin_lo = edges[0];
  double first_bin_hi = edges[1];

  int first_bin_count = 0;
  double first_bin_sum_Nhit = 0.0;

  std::cout << "\n=== DEBUG: ultra-close entries in first bin (distance < 1 cm) ===" << std::endl;

  int n_ultraclose = 0;
  double sum_ultraclose_weight = 0.0;

  for (size_t i = 0; i < v_prop_dist.size(); ++i) {
    if (v_prop_dist[i] < edges[0] || v_prop_dist[i] >= edges[1]) continue;
    if (v_hits_sim[i] <= 0) continue;
    if (v_prop_dist[i] >= 1.0) continue;
    

    double discrepancy = (v_hits_geo[i] - v_hits_sim[i]) / v_hits_sim[i];

    n_ultraclose++;
    sum_ultraclose_weight += v_hits_sim[i];

    std::cout
      << " entry " << n_ultraclose
      << " : distance=" << v_prop_dist[i]
      << " cm"
      << " , N_hit=" << v_hits_sim[i]
      << " , N_omega=" << v_hits_geo[i]
      << " , bias=" << discrepancy
      << " , dT=" << v_r[i]
      << std::endl;
  }

  std::cout << "Number of entries (<1 cm) in first bin = "
            << n_ultraclose << std::endl;
  std::cout << "Total N_hit for those entries = "
            << sum_ultraclose_weight << std::endl;
  std::cout << "====================================================\n" << std::endl;

  auto MakeBiasPlot = [&](TProfile* prof, const std::string& outName, const std::string& title) {
    std::vector<double> xs, exs, means, rmss;
    for (int iBin = 1; iBin <= n_bins; ++iBin) {
      xs.push_back(prof->GetXaxis()->GetBinCenter(iBin));
      exs.push_back(0.5 * (prof->GetXaxis()->GetBinUpEdge(iBin) - prof->GetXaxis()->GetBinLowEdge(iBin)));
      means.push_back(prof->GetBinContent(iBin));
      rmss.push_back(prof->GetBinError(iBin));
    }

    std::cout << "\n" << outName << std::endl;
    for (int iBin = 1; iBin <= n_bins; ++iBin) {
      std::cout
        << "bin " << iBin
        << " x=" << prof->GetXaxis()->GetBinCenter(iBin)
        << " entries=" << prof->GetBinEntries(iBin)
        << " bias=" << prof->GetBinContent(iBin)
        << " rms=" << prof->GetBinError(iBin)
        << std::endl;
    }

    TCanvas* c = new TCanvas(("c_" + outName).c_str(), "", 200, 10, 1080, 1080);
    c->SetGrid();
    c->SetBottomMargin(0.105);
    c->SetLeftMargin(0.16);
    

    TGraphErrors* grRMS = new TGraphErrors(xs.size(), &xs[0], &rmss[0], &exs[0], 0);
    TGraphErrors* grBias = new TGraphErrors(xs.size(), &xs[0], &means[0], &exs[0], 0);

    grRMS->GetXaxis()->SetTitle("distance [cm]");
    grRMS->GetYaxis()->SetTitle("(N_{pred} - N_{Geant4}) / N_{Geant4}");
    grRMS->SetTitle(title.c_str());
    grRMS->GetXaxis()->SetRangeUser(0, 1000);
    grRMS->GetYaxis()->SetRangeUser(-0.5, 0.5);
    grRMS->SetMarkerStyle(kFullSquare);
    grRMS->SetMarkerSize(1.8);
    grRMS->SetMarkerColor(kBlack);

    grBias->SetMarkerStyle(kFullCircle);
    grBias->SetMarkerSize(1.8);
    grBias->SetMarkerColor(kRed);

    grRMS->Draw("AP");
    grBias->Draw("P SAME");

    TLegend* leg = new TLegend(0.60,0.70,0.90,0.90,NULL,"brNDC");
    leg->AddEntry(grRMS, "RMS", "P");
    leg->AddEntry(grBias, "Bias", "P");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw("same");

    c->SaveAs((save_dir + "/" + outName).c_str());
  };

  MakeBiasPlot(profile, "validation_bias_all.pdf", "DUNE FD-HD validation: all selected pairs");
  MakeBiasPlot(profile_core, "validation_bias_dTlt600.pdf", "DUNE FD-HD validation: d_{T} < 600 cm");

  delete f;
}

// ============================================================
// Main
// ============================================================
int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Please specify input file" << std::endl;
    return 1;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  std::string positions = "./protodune_optical_mapping.txt";
  std::string input_file = argv[1];
  std::string file_name = "output";
  std::string save_dir = "./plots";
  gSystem->Exec(("mkdir -p " + save_dir).c_str());

  for (int i = 0; i < N; ++i) g_theta_centers.push_back(delta_angulo/2.0 + i*delta_angulo);

  PairData pairs = BuildPairData(positions, input_file);

  TH1D* h = new TH1D("","", range_d, 0, range_d*1.05);
  TH1D* hd_centers[M];
  TLine* line[M+1];
  TF1* GH[N][M];
  std::string options[N];
  for (int k = 0; k < N; ++k) options[k] = "W0Q";

  for (int k = 0; k < M; ++k) {
    hd_centers[k] = new TH1D("","", M, range_d_array);
    line[k] = new TLine(range_d_array[k], 0., range_d_array[k], 2000);
    for (int j = 0; j < N; ++j) {
      GH[j][k] = new TF1("GH", GaisserHillas, 0., d_max, 4);
      double pars_ini[4] = {1., 128., 55., -2500.};
      GH[j][k]->SetParLimits(2, 10, 350);

      if      (j < 1) pars_ini[3] = -5000;
      else if (j < 2) pars_ini[3] = -3000;
      else if (j < 4) pars_ini[3] = -3000;
      else if (j < 5) pars_ini[3] = -3000;
      else if (j < 6) pars_ini[3] = -1000;
      else if (j < 8) pars_ini[3] = -500;
      else            pars_ini[3] = -100;

      GH[j][k]->FixParameter(3, pars_ini[3]);
      GH[j][k]->SetParameters(pars_ini);
    }
  }

  TProfile* pdiff_d[N][M];
  TProfile* pdiff_dx[N][M];
  for (int j = 0; j < N; ++j) {
    for (int k = 0; k < M; ++k) {
      pdiff_d[j][k]  = new TProfile("","", int((d_max-d_min)/step_d), d_min, d_max, "s");
      pdiff_dx[j][k] = new TProfile("","", int((d_max-d_min)/step_d), d_min, d_max, "s");
    }
  }

  for (size_t i = 0; i < pairs.v_distance.size(); ++i) {
    double costheta = std::cos(pi * pairs.v_offset_angle.at(i) / 180.0);
    int j = int(pairs.v_offset_angle.at(i) / delta_angulo);
    if (j < 0 || j >= N) continue;

    double temp = pairs.v_d_center.at(i);
    int k = std::lower_bound(range_d_array, range_d_array + M + 1, temp) - range_d_array - 1;
    if (k < 0 || k >= M) continue;
    if (pairs.v_rec_hits.at(i) <= 0) continue;
    if (pairs.v_distance.at(i) < 1.0) continue; 

    pdiff_d[j][k]->Fill(pairs.v_distance.at(i), pairs.v_hits.at(i)/pairs.v_rec_hits.at(i)*costheta);
    pdiff_dx[j][k]->Fill(pairs.v_distance.at(i), pairs.v_distance.at(i));
    hd_centers[k]->Fill(pairs.v_d_center.at(i));
    h->Fill(pairs.v_d_center.at(i));
    
  }

  double max_x[N][M], min_x[N][M], n_entries[N][M];
  for (int j = 0; j < N; ++j) {
    for (int k = 0; k < M; ++k) {
      n_entries[j][k] = 0;
      TAxis* xaxis = pdiff_d[j][k]->GetXaxis();
      int nbins = xaxis->GetNbins();
      double minv = d_max;
      double maxv = 0;
      for (int bin = 0; bin <= nbins; ++bin) {
        n_entries[j][k] += pdiff_d[j][k]->GetBinContent(bin);
        if (pdiff_d[j][k]->GetBinContent(bin) == 0) continue;
        if (minv > xaxis->GetBinCenter(bin)) minv = xaxis->GetBinCenter(bin);
        if (maxv < xaxis->GetBinCenter(bin)) maxv = xaxis->GetBinCenter(bin);
      }
      max_x[j][k] = maxv;
      min_x[j][k] = minv;
    }
  }

  std::vector<int> N_canvas;
  std::vector<double> d_center;
  std::string title[M];
  for (int k = 0; k < M; ++k) {
    for (int j = 0; j < N; ++j) {
      if (n_entries[j][k] > 0) {
        N_canvas.push_back(k);
        break;
      }
    }
    if (hd_centers[k]->GetEntries() > 0) {
      d_center.push_back(hd_centers[k]->GetMean());
      title[k] = "Distance to center " + std::to_string(int(hd_centers[k]->GetMean())) + "cm";
    }
  }

  const int dim = N_canvas.size();

  // crown plot
  TCanvas *canvas0 = new TCanvas("canvas0", "graph draw options", 200, 200, 500, 400);
  h->SetTitle("Help to choose the \"border-study\" distance bins");
  h->GetXaxis()->SetTitle("distance to centre in Y-Z plane [cm]");
  h->SetStats(0);
  h->Draw("hist");
  line[M] = new TLine(range_d, 0., range_d, 2000);
  for (int k = 0; k < M+1; k++) {
    line[k]->SetLineColor(2);
    line[k]->SetLineWidth(3);
    line[k]->SetLineStyle(kDashed);
    line[k]->Draw("same");
  }
  canvas0->SaveAs((save_dir + "/" + file_name + "_crown_plot_double_sided.pdf").c_str());

  TGraphErrors* gr[N][M];
  for (int jj = 0; jj < N; ++jj)
    for (int kk = 0; kk < M; ++kk)
      gr[jj][kk] = nullptr;

  for (int k = 0; k < M; ++k) {
    for (int j = 0; j < N; ++j) {
      if (n_entries[j][k] > 0) {
        gr[j][k] = Profile_to_Graph(pdiff_dx[j][k], pdiff_d[j][k], j);
        gr[j][k]->SetLineColor(1 + j%9);
        gr[j][k]->SetMarkerColor(1 + j%9);
        if (j == 4) {
          gr[j][k]->SetLineColor(kOrange+7);
          gr[j][k]->SetMarkerColor(kOrange+7);
        }
        gr[j][k]->SetMarkerStyle(20 + j);
        gr[j][k]->SetMarkerSize(0.6);
      }
    }
  }

  std::vector<double> p1[64], p2[64], p3[64], p4[64], ep1[64], ep2[64], ep3[64];
  std::vector<double> a1, a2, a3;

  // GH multi-pad canvas
  const int maxRows = 3;
  int nCols = (dim <= 3) ? dim : std::min(4, dim);
  int nRows = (dim + nCols - 1) / nCols;
  nRows = std::min(nRows, maxRows);
  const int padW = 900;
  const int padH = 700;

  TCanvas *canvas1 = new TCanvas("canvas1", "GH", 200, 200, nCols*padW, nRows*padH);
  canvas1->Divide(nCols, nRows, 0.001, 0.001);

  double x_0[2] = {0, d_max};
  double y_0[2] = {0, 2.};
  TGraph* gg0[64];
  TLegend *leg1 = new TLegend(0.62, 0.55, 0.95, 0.93, NULL, "brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.035);
  char label[N][40];

  for (int l = 0; l < dim; ++l) {
    int k = N_canvas.at(l);
    canvas1->cd(l + 1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.10);

    gg0[l] = new TGraph(2, x_0, y_0);
    gg0[l]->SetTitle(title[k].c_str());
    gg0[l]->GetXaxis()->SetTitle("distance [cm]");
    gg0[l]->GetYaxis()->SetTitle("N_{hit} / N_{#Omega} / cos(#theta)");

    double yMax = 0.0;
    for (int j = 0; j < N; ++j) {
      if (n_entries[j][k] <= 0 || gr[j][k] == nullptr) continue;
      int nPoints = gr[j][k]->GetN();
      for (int i = 0; i < nPoints; ++i) {
        double x, y;
        gr[j][k]->GetPoint(i, x, y);
        double ey = gr[j][k]->GetErrorY(i);
        yMax = std::max(yMax, y + ey);
      }
    }

    if (yMax <= 0) yMax = 1.0;
    gg0[l]->GetXaxis()->SetRangeUser(0, d_max);
    gg0[l]->GetYaxis()->SetRangeUser(0, yMax * 1.15);
    gg0[l]->GetXaxis()->SetLabelSize(0.045);
    gg0[l]->GetXaxis()->SetTitleSize(0.055);
    gg0[l]->GetXaxis()->SetTitleOffset(1.1);
    gg0[l]->GetYaxis()->SetLabelSize(0.045);
    gg0[l]->GetYaxis()->SetTitleSize(0.055);
    gg0[l]->GetYaxis()->SetTitleOffset(1.10);
    gg0[l]->Draw("AP");

    for (int j = 0; j < N; ++j) {
      double pars_GH[4]  = {-999, -999, -999, -999};
      double epars_GH[4] = {-999, -999, -999, -999};

      if (n_entries[j][k] > 0 && gr[j][k] != nullptr) {
        gr[j][k]->Draw("P");
        double fit_lo = min_x[j][k];
        double fit_hi = std::min(FIT_MAX, max_x[j][k]);
        if (fit_lo < fit_hi) gr[j][k]->Fit(GH[j][k], options[j].c_str(), "", fit_lo, fit_hi);
        GH[j][k]->GetParameters(pars_GH);
        GH[j][k]->SetLineColor(1 + j);
        if (j == 4) GH[j][k]->SetLineColor(kOrange + 7);
        GH[j][k]->SetLineStyle(kDashed);
        GH[j][k]->Draw("SAME");

        const double* epars_GH_ = GH[j][k]->GetParErrors();
        epars_GH[0] = epars_GH_[0];
        epars_GH[1] = epars_GH_[1];
        epars_GH[2] = epars_GH_[2];
      }

      p1[l].push_back(pars_GH[0]);
      p2[l].push_back(pars_GH[1]);
      p3[l].push_back(pars_GH[2]);
      p4[l].push_back(pars_GH[3]);
      ep1[l].push_back(epars_GH[0]);
      ep2[l].push_back(epars_GH[1]);
      ep3[l].push_back(epars_GH[2]);

      if (l == 0 && n_entries[j][k] > 0 && gr[j][k] != nullptr) {
        int a_min = j * delta_angulo;
        int a_max = (j + 1) * delta_angulo;
        sprintf(label[j], "#theta #in [%i, %i] deg", a_min, a_max);
        leg1->AddEntry(gr[j][k], label[j], "p");
      }
    }

    if (l == 0) leg1->Draw();
  }

  canvas1->SaveAs((save_dir + "/" + file_name + "_GH_double_sided.pdf").c_str());

  // overlay by angle canvases
  std::vector<int> chosenK;
  const int wantK = 6;
  if (dim <= wantK) {
    chosenK = N_canvas;
  } else {
    chosenK.reserve(wantK);
    for (int i = 0; i < wantK; ++i) {
      int idx = int(std::round(i * (dim - 1.0) / (wantK - 1.0)));
      chosenK.push_back(N_canvas.at(idx));
    }
  }

  const Color_t col6[6] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7};
  const Style_t mark6[6] = {20, 21, 22, 23, 33, 34};

  for (int j = 0; j < N; ++j) {
    bool hasAny = false;
    for (int ii = 0; ii < (int)chosenK.size(); ++ii) {
      int k = chosenK[ii];
      if (n_entries[j][k] > 0 && gr[j][k] != nullptr) { hasAny = true; break; }
    }
    if (!hasAny) continue;

    TString cname; cname.Form("c_angle_%d", j);
    TCanvas* cAng = new TCanvas(cname, cname, 900, 650);
    cAng->SetLeftMargin(0.14);
    cAng->SetRightMargin(0.05);
    cAng->SetBottomMargin(0.14);
    cAng->SetTopMargin(0.08);

    TH1F* frame = new TH1F("", "", 100, 0, d_max);
    frame->GetXaxis()->SetTitle("distance [cm]");
    frame->GetYaxis()->SetTitle("N_{hit} / N_{#Omega} / cos(#theta)");
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleOffset(1.10);

    double yMax = 0.0;
    for (int k = 0; k < M; ++k) {
      if (!gr[j][k] || n_entries[j][k] <= 0) continue;
      int nPoints = gr[j][k]->GetN();
      for (int i = 0; i < nPoints; ++i) {
        double x, y;
        gr[j][k]->GetPoint(i, x, y);
        double ey = gr[j][k]->GetErrorY(i);
        yMax = std::max(yMax, y + ey);
      }
    }

    if (yMax <= 0) yMax = 1.0;
    frame->GetYaxis()->SetRangeUser(0., yMax * 1.25);
    frame->Draw("AXIS");

    TString t; t.Form("#theta #in [%d, %d] deg", int(j*delta_angulo), int((j+1)*delta_angulo));
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.045);
    lat.DrawLatex(0.16, 0.93, t);

    TLegend* legA = new TLegend(0.60, 0.58, 0.94, 0.92);
    legA->SetFillStyle(0);
    legA->SetBorderSize(0);
    legA->SetTextSize(0.035);

    for (int ii = 0; ii < (int)chosenK.size(); ++ii) {
      int k = chosenK[ii];
      if (n_entries[j][k] <= 0 || gr[j][k] == nullptr) continue;

      gr[j][k]->SetLineColor(col6[ii]);
      gr[j][k]->SetMarkerColor(col6[ii]);
      gr[j][k]->SetMarkerStyle(mark6[ii]);
      gr[j][k]->SetMarkerSize(0.7);
      gr[j][k]->Draw("P SAME");

      GH[j][k]->SetLineColor(col6[ii]);
      GH[j][k]->SetLineStyle(kDashed);
      GH[j][k]->SetLineWidth(2);
      GH[j][k]->Draw("SAME");

      double dmean = (hd_centers[k]->GetEntries() > 0) ? hd_centers[k]->GetMean() : d_center[k];
      TString lab; lab.Form("d_{T} #approx %.0f cm", dmean);
      legA->AddEntry(gr[j][k], lab, "p");
    }

    legA->Draw();

    TString out;
    out.Form("%s/%s_overlayByAngle_theta%02d.pdf", save_dir.c_str(), file_name.c_str(), j);
    cAng->SaveAs(out);
  }

  std::vector<double> vd_center;
  std::vector<double> eslopes1, eslopes2, eslopes3;

  TGraphErrors* gNmax[N];
  TGraphErrors* gdmax[N];
  TGraphErrors* glambda[N];
  TF1* f1[N];
  TF1* f2[N];
  TF1* f3[N];

  for (int j = 0; j < N; ++j) {
    std::vector<double> vNmax, vdmax, vlambda;
    std::vector<double> veNmax, vedmax, velambda;

    if (j == 0) vd_center.clear();

    for (int k = 0; k < dim; ++k) {
      vNmax.push_back(p1[k].at(j));
      vdmax.push_back(p2[k].at(j));
      vlambda.push_back(p3[k].at(j));

      veNmax.push_back(ep1[k].at(j));
      vedmax.push_back(ep2[k].at(j));
      velambda.push_back(ep3[k].at(j));

      if (j == 0) vd_center.push_back(d_center.at(k));
    }

    gNmax[j]   = new TGraphErrors(vd_center.size(), &vd_center[0], &vNmax[0],   0, &veNmax[0]);
    gdmax[j]   = new TGraphErrors(vd_center.size(), &vd_center[0], &vdmax[0],   0, &vedmax[0]);
    glambda[j] = new TGraphErrors(vd_center.size(), &vd_center[0], &vlambda[0], 0, &velambda[0]);

    f1[j] = new TF1(Form("f1_%d", j), fUseQuadraticNmax ? "pol2" : "pol1", 0., range_d);
    f2[j] = new TF1(Form("f2_%d", j), "pol1", 0., range_d);
    f3[j] = new TF1(Form("f3_%d", j), "pol1", 0., range_d);

    f1[j]->SetLineColor(1 + j);
    f2[j]->SetLineColor(1 + j);
    f3[j]->SetLineColor(1 + j);

    if (j == 4) {
      f1[j]->SetLineColor(kOrange+7);
      f2[j]->SetLineColor(kOrange+7);
      f3[j]->SetLineColor(kOrange+7);
    }

    f1[j]->SetLineStyle(kDashed);
    f2[j]->SetLineStyle(kDashed);
    f3[j]->SetLineStyle(kDashed);

    double bf_hi = std::min(BORDER_FIT_MAX, range_d);
    gNmax[j]->Fit(f1[j], "Q", "", 0, bf_hi);
    gdmax[j]->Fit(f2[j], "Q", "", 0, bf_hi);
    glambda[j]->Fit(f3[j], "Q", "", 0, bf_hi);

    g_slopes1.push_back(f1[j]->GetParameter(1));
    g_slopes2.push_back(f2[j]->GetParameter(1));
    g_slopes3.push_back(f3[j]->GetParameter(1));

    eslopes1.push_back(f1[j]->GetParError(1));
    eslopes2.push_back(f2[j]->GetParError(1));
    eslopes3.push_back(f3[j]->GetParError(1));

    a1.push_back(f1[j]->GetParameter(0));
    a2.push_back(f2[j]->GetParameter(0));
    a3.push_back(f3[j]->GetParameter(0));

    gNmax[j]->SetLineColor(1 + j);
    gNmax[j]->SetMarkerColor(1 + j);
    gNmax[j]->SetMarkerStyle(20 + j);

    gdmax[j]->SetLineColor(1 + j);
    gdmax[j]->SetMarkerColor(1 + j);
    gdmax[j]->SetMarkerStyle(20 + j);

    glambda[j]->SetLineColor(1 + j);
    glambda[j]->SetMarkerColor(1 + j);
    glambda[j]->SetMarkerStyle(20 + j);

    if (j == 4) {
      gNmax[j]->SetLineColor(kOrange+7);
      gNmax[j]->SetMarkerColor(kOrange+7);
      gdmax[j]->SetLineColor(kOrange+7);
      gdmax[j]->SetMarkerColor(kOrange+7);
      glambda[j]->SetLineColor(kOrange+7);
      glambda[j]->SetMarkerColor(kOrange+7);
    }
  }

  for (int i = 0; i < N; ++i) {
    g_fit_p1.push_back(a1.at(i));
    g_fit_p2.push_back(a2.at(i));
    g_fit_p3.push_back(a3.at(i));
    g_fit_p4.push_back(p4[0].at(i));
  }

  std::cout << "\nCorrections to plug into LArSoft PhotonVisibilityServices" << std::endl;
  std::cout << "double p1[9] = {";
  for (int i = 0; i < N; ++i) std::cout << g_fit_p1[i] << (i+1<N ? ", " : "};\n");
  std::cout << "double p2[9] = {";
  for (int i = 0; i < N; ++i) std::cout << g_fit_p2[i] << (i+1<N ? ", " : "};\n");
  std::cout << "double p3[9] = {";
  for (int i = 0; i < N; ++i) std::cout << g_fit_p3[i] << (i+1<N ? ", " : "};\n");
  std::cout << "double p4[9] = {";
  for (int i = 0; i < N; ++i) std::cout << g_fit_p4[i] << (i+1<N ? ", " : "};\n");

  std::cout << "double slopes1[9] = {";
  for (int i = 0; i < N; ++i) std::cout << g_slopes1[i] << (i+1<N ? ", " : "};\n");
  std::cout << "double slopes2[9] = {";
  for (int i = 0; i < N; ++i) std::cout << g_slopes2[i] << (i+1<N ? ", " : "};\n");
  std::cout << "double slopes3[9] = {";
  for (int i = 0; i < N; ++i) std::cout << g_slopes3[i] << (i+1<N ? ", " : "};\n");

  // border slopes canvas
  double xx[3][2];
  double yy[3][2];
  TGraphErrors* gf[3];
  for (int j = 0; j < 3; j++) {
    xx[j][0] = 0;
    xx[j][1] = 1200;
    yy[j][0] = 0;
    yy[j][1] = 200;
    gf[j] = new TGraphErrors(2, xx[j], yy[j]);
    gf[j]->SetTitle("");
  }

  std::vector<double> slopes1_plot, slopes2_plot, slopes3_plot;
  std::vector<double> eslopes1_plot, eslopes2_plot, eslopes3_plot;

  for (int i = 0; i < N; ++i) {
    slopes1_plot.push_back(1000.0 * g_slopes1[i]);
    slopes2_plot.push_back(g_slopes2[i]);
    slopes3_plot.push_back(g_slopes3[i]);

    eslopes1_plot.push_back(1000.0 * eslopes1[i]);
    eslopes2_plot.push_back(eslopes2[i]);
    eslopes3_plot.push_back(eslopes3[i]);
  }

  TGraphErrors* g1 = new TGraphErrors(g_theta_centers.size(), &g_theta_centers[0], &slopes1_plot[0], 0, &eslopes1_plot[0]);
  TGraphErrors* g2 = new TGraphErrors(g_theta_centers.size(), &g_theta_centers[0], &slopes2_plot[0], 0, &eslopes2_plot[0]);
  TGraphErrors* g3 = new TGraphErrors(g_theta_centers.size(), &g_theta_centers[0], &slopes3_plot[0], 0, &eslopes3_plot[0]);
  g1->SetTitle("");
  g2->SetTitle("");
  g3->SetTitle("");
  TF1* p0_m1 = new TF1("p0_m1", "pol1", 0., 90.);
  TF1* p0_m2 = new TF1("p0_m2", "pol1", 0., 90.);
  TF1* p0_m3 = new TF1("p0_m3", "pol1", 0., 90.);
  g1->Fit(p0_m1, "W0Q", "", 0, 90);
  g2->Fit(p0_m2, "W0Q", "", 0, 90);
  g3->Fit(p0_m3, "W0Q", "", 0, 90);

  TCanvas *c = new TCanvas("c", "canvas", 900, 2000);
  c->Divide(1,3);

  c->cd(1);
  TPad *pad1 = new TPad("pad1", "pad1", 0., 0.37, 1, 1.0);
  pad1->SetBottomMargin(0.1);
  pad1->SetLeftMargin(0.16);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();
  gf[0]->GetYaxis()->SetRangeUser(0, 2.);
  gf[0]->GetXaxis()->SetRangeUser(0, range_d*1.05);
  gf[0]->GetXaxis()->SetTitle("d_{T} [cm]");
  gf[0]->GetYaxis()->SetTitle("N_{max}");
  gf[0]->GetXaxis()->SetLabelSize(0.05);
  gf[0]->GetXaxis()->SetTitleSize(0.06);
  gf[0]->GetXaxis()->SetTitleOffset(1.0);
  gf[0]->GetYaxis()->SetLabelSize(0.05);
  gf[0]->GetYaxis()->SetTitleSize(0.06);
  gf[0]->GetYaxis()->SetTitleOffset(1.15);
  gf[0]->Draw("AP");
  for (int j = 0; j < N; j++) {
    gNmax[j]->Draw("P SAME");
    f1[j]->Draw("L SAME");
  }

  c->cd(1);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
  pad2->SetTopMargin(0.00);
  pad2->SetBottomMargin(0.20);
  pad2->SetLeftMargin(0.16);
  pad2->SetRightMargin(0.04);
  pad2->Draw();
  pad2->cd();
  g1->GetXaxis()->SetTitle("#theta [deg]");
  g1->GetYaxis()->SetTitle("slope N_{max} [cm^{-1}] x10^{-3}");
  g1->GetXaxis()->SetLabelSize(0.08);
  g1->GetXaxis()->SetTitleSize(0.09);
  g1->GetXaxis()->SetTitleOffset(0.95);
  g1->GetYaxis()->SetLabelSize(0.08);
  g1->GetYaxis()->SetTitleSize(0.09);
  g1->GetYaxis()->SetTitleOffset(0.75);
  g1->SetMarkerStyle(20);
  g1->Draw("AP");
  p0_m1->SetLineColor(kGray+1);
  p0_m1->SetLineStyle(kDashed);
  p0_m1->Draw("SAME");

  c->cd(2);
  TPad *padc1 = new TPad("padc1", "padc1", 0., 0.37, 1, 1.0);
  padc1->SetBottomMargin(0.10);
  padc1->SetLeftMargin(0.16);
  padc1->SetRightMargin(0.04);
  padc1->Draw();
  padc1->cd();
  gf[1]->GetYaxis()->SetRangeUser(70., 300.);
  gf[1]->GetXaxis()->SetRangeUser(0, range_d*1.05);
  gf[1]->GetXaxis()->SetTitle("d_{T} [cm]");
  gf[1]->GetYaxis()->SetTitle("d_{max} [cm]");
  gf[1]->GetXaxis()->SetLabelSize(0.05);
  gf[1]->GetXaxis()->SetTitleSize(0.06);
  gf[1]->GetXaxis()->SetTitleOffset(1.0);
  gf[1]->GetYaxis()->SetLabelSize(0.05);
  gf[1]->GetYaxis()->SetTitleSize(0.06);
  gf[1]->GetYaxis()->SetTitleOffset(1.15);
  gf[1]->Draw("AP");
  for (int j = 0; j < N; j++) {
    gdmax[j]->Draw("P SAME");
    f2[j]->Draw("L SAME");
  }

  c->cd(2);
  TPad *padc2 = new TPad("padc2", "padc2", 0, 0, 1, 0.35);
  padc2->SetTopMargin(0.00);
  padc2->SetBottomMargin(0.20);
  padc2->SetLeftMargin(0.16);
  padc2->SetRightMargin(0.04);
  padc2->Draw();
  padc2->cd();
  g2->GetXaxis()->SetTitle("#theta [deg]");
  g2->GetYaxis()->SetTitle("slope d_{max}");
  g2->GetXaxis()->SetLabelSize(0.08);
  g2->GetXaxis()->SetTitleSize(0.09);
  g2->GetXaxis()->SetTitleOffset(0.95);
  g2->GetYaxis()->SetLabelSize(0.08);
  g2->GetYaxis()->SetTitleSize(0.09);
  g2->GetYaxis()->SetTitleOffset(0.75);
  g2->SetMarkerStyle(20);
  g2->Draw("AP");
  p0_m2->SetLineColor(kGray+1);
  p0_m2->SetLineStyle(kDashed);
  p0_m2->Draw("SAME");

  c->cd(3);
  TPad *padcc1 = new TPad("padcc1", "padcc1", 0., 0.37, 1, 1.0);
  padcc1->SetBottomMargin(0.10);
  padcc1->SetLeftMargin(0.16);
  padcc1->SetRightMargin(0.04);
  padcc1->Draw();
  padcc1->cd();
  gf[2]->GetYaxis()->SetRangeUser(-10, 210);
  gf[2]->GetXaxis()->SetRangeUser(0, range_d*1.05);
  gf[2]->GetXaxis()->SetTitle("d_{T} [cm]");
  gf[2]->GetYaxis()->SetTitle("#Lambda [cm]");
  gf[2]->Draw("AP");
  for (int j = 0; j < N; j++) {
    glambda[j]->Draw("P SAME");
    f3[j]->Draw("L SAME");
  }

  c->cd(3);
  TPad *padcc2 = new TPad("padcc2", "padcc2", 0, 0, 1, 0.35);
  padcc2->SetTopMargin(0.00);
  padcc2->SetBottomMargin(0.20);
  padcc2->SetLeftMargin(0.16);
  padcc2->SetRightMargin(0.04);
  padcc2->Draw();
  padcc2->cd();
  g3->GetXaxis()->SetTitle("#theta [deg]");
  g3->GetYaxis()->SetTitle("slope #Lambda");
  g3->SetMarkerStyle(20);
  g3->Draw("AP");
  p0_m3->SetLineColor(kGray+1);
  p0_m3->SetLineStyle(kDashed);
  p0_m3->Draw("SAME");

  c->SaveAs((save_dir + "/" + file_name + "_border_slopes_double_sided.pdf").c_str());

  RunValidation(input_file, positions, save_dir);
  return 0;
}

// ============================================================
// Shared math utilities
// ============================================================
Double_t GaisserHillas(double x, double *par) {
  Double_t X_mu_0 = par[3];
  Double_t Normalization = par[0];
  Double_t Diff = par[1] - X_mu_0;
  Double_t Term = pow((x - X_mu_0)/Diff, Diff/par[2]);
  Double_t Exponential = TMath::Exp((par[1] - x)/par[2]);
  return Normalization * Term * Exponential;
}

double omega(const double &a, const double &b, const double &d) {
  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));
}

double solid(const acc& out, const TVector3 &v) {
  if (v.Y()==0.0 && v.Z()==0.0) return omega(out.w,out.h,v.X());
  if ((std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)) {
    double A = std::abs(v.Y())-out.w/2.0, B = std::abs(v.Z())-out.h/2.0, d = std::abs(v.X());
    return (omega(2*(A+out.w),2*(B+out.h),d)-omega(2*A,2*(B+out.h),d)-omega(2*(A+out.w),2*B,d)+omega(2*A,2*B,d))/4.0;
  }
  if ((std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)) {
    double A = -std::abs(v.Y())+out.w/2.0, B = -std::abs(v.Z())+out.h/2.0, d = std::abs(v.X());
    return (omega(2*(out.w-A),2*(out.h-B),d)+omega(2*A,2*(out.h-B),d)+omega(2*(out.w-A),2*B,d)+omega(2*A,2*B,d))/4.0;
  }
  if ((std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)) {
    double A = std::abs(v.Y())-out.w/2.0, B = -std::abs(v.Z())+out.h/2.0, d = std::abs(v.X());
    return (omega(2*(A+out.w),2*(out.h-B),d)-omega(2*A,2*(out.h-B),d)+omega(2*(A+out.w),2*B,d)-omega(2*A,2*B,d))/4.0;
  }
  if ((std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)) {
    double A = -std::abs(v.Y())+out.w/2.0, B = std::abs(v.Z())-out.h/2.0, d = std::abs(v.X());
    return (omega(2*(out.w-A),2*(B+out.h),d)-omega(2*(out.w-A),2*B,d)+omega(2*A,2*(B+out.h),d)-omega(2*A,2*B,d))/4.0;
  }
  return 0.0;
}

double solid_lateral(const acc& out, const TVector3 &v) {
  if (v.X()==0.0 && v.Z()==0.0) return omega(out.w,out.h,v.Y());
  if ((std::abs(v.X()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)) {
    double A = std::abs(v.X())-out.w/2.0, B = std::abs(v.Z())-out.h/2.0, d = std::abs(v.Y());
    return (omega(2*(A+out.w),2*(B+out.h),d)-omega(2*A,2*(B+out.h),d)-omega(2*(A+out.w),2*B,d)+omega(2*A,2*B,d))/4.0;
  }
  if ((std::abs(v.X()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)) {
    double A = -std::abs(v.X())+out.w/2.0, B = -std::abs(v.Z())+out.h/2.0, d = std::abs(v.Y());
    return (omega(2*(out.w-A),2*(out.h-B),d)+omega(2*A,2*(out.h-B),d)+omega(2*(out.w-A),2*B,d)+omega(2*A,2*B,d))/4.0;
  }
  if ((std::abs(v.X()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)) {
    double A = std::abs(v.X())-out.w/2.0, B = -std::abs(v.Z())+out.h/2.0, d = std::abs(v.Y());
    return (omega(2*(A+out.w),2*(out.h-B),d)-omega(2*A,2*(out.h-B),d)+omega(2*(A+out.w),2*B,d)-omega(2*A,2*B,d))/4.0;
  }
  if ((std::abs(v.X()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)) {
    double A = -std::abs(v.X())+out.w/2.0, B = std::abs(v.Z())-out.h/2.0, d = std::abs(v.Y());
    return (omega(2*(out.w-A),2*(B+out.h),d)-omega(2*(out.w-A),2*B,d)+omega(2*A,2*(B+out.h),d)-omega(2*A,2*B,d))/4.0;
  }
  return 0.0;
}

double interpolate(const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate) {
  int size = xData.size();
  int i = 0;
  if (x >= xData[size - 2]) i = size - 2;
  else while (x > xData[i+1]) i++;
  double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];
  if (!extrapolate) {
    if (x < xL) yR = yL;
    if (x > xR) yL = yR;
  }
  double dydx = (yR - yL) / (xR - xL);
  return yL + dydx * (x - xL);
}


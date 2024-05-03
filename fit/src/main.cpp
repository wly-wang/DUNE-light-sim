// This macro calculates the corrections needed in the Semi-analytic
// scintillation light model. These corrections are parametrised as
// Gaisser-Hillas functions. These functions depend on: 1) distance,
// 2) offset-angle, and 3) distance to the center (border) of the
// LArTPC active volume. It allows to correct for the LAr light signals
// for optical detector with flat sensitive windows (circular and
// rectangular shepes) or semisphere (by an approach model better than
// few%).

// Author: Diego Garcia-Gamez (dgamez@fnal.gov)
// Created: 15.02.2021

#include "Rtypes.h"
#include "TF1.h"
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

#include <vector>
#include <fstream>
#include <iostream>

#include "functions.h"

void calcula(std::string positions, std::string input_file, std::vector<double> &v1, std::vector<double> &v2,
       std::vector<double> &v3, std::vector<double> &v4, std::vector<double> &v5, bool IsRectangular,
       bool IsSphere, std::vector<double> &v6, std::vector<double> &v7) {
  std::cout<<"calcula function ..."<<std::endl;

  double min_number_entries = 1000;
  // width and height in cm of single arapuca active window
  double arapuca_w = 10;
  double arapuca_h = 47.75;
  // 8" PMT radius
  double b = 8*2.54/2.;
  // Y-Z coordinates of the active volume center
  const double centerYZ[2] = {305.534, 231.266};
  //const double centerXZ[2] = {0., 1000.};
  //double x_anode = 324.97;//325.01;
  //double y_anode = 0;
  //double z_anode = 1000;

  // Modification requried from Users:
  // LAr absorption length in cm
  // This needs to match with the value used in the full-geant4 simulation
  //const double L_abs = 8000.; //xenon
  const double L_abs = 2000.; //argon

  gRandom->SetSeed(0);
  //getting the pmt positions (y and z)
  std::ifstream Traks_file1(positions.c_str());
  if(!Traks_file1) {
    std::cerr << "WARNING:  Failed to open file with optical detector positions"<< std::endl;
  }
  Traks_file1.seekg(0);
  std::vector<double> devx;
  std::vector<double> devy;
  std::vector<double> devz;

  double id, x, y, z;
  /*while(!(Traks_file1.eof())) {
    Traks_file1 >> id >> x >> y >> z;
    devx.push_back(x);
    devy.push_back(y);
    devz.push_back(z);
  }*/
  std::string line;
  while(std::getline(Traks_file1, line)) {
    sscanf(line.c_str(), "%lf %lf %lf %lf", &id, &x, &y, &z);
    devx.push_back(x);
    devy.push_back(y);
    devz.push_back(z);
  }

  std::cout<<devx.size()<<std::endl;

  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"Input File: "<<input_file<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;

  std::vector<double> v_hits, v_distance, v_rec_hits, v_offset_angle, v_d_center, v_x, v_devx;
  //loop over files
  TFile* f = new TFile(input_file.c_str());
  TTree *tree = (TTree *)f->Get("myTree");
  const Int_t kMaxDevices = 3000;
  int VUV_hits[kMaxDevices];
  int Vis_hits[kMaxDevices];
  double X, Y, Z;
  int genPhotons, numberDevices;
  tree->SetBranchAddress("numberDevices",&numberDevices);
  tree->SetBranchAddress("X", &X);
  tree->SetBranchAddress("Y", &Y);
  tree->SetBranchAddress("Z", &Z);
  tree->SetBranchAddress("VUV_hits",VUV_hits);
  tree->SetBranchAddress("Vis_hits",Vis_hits);
  tree->SetBranchAddress("genPhotons", &genPhotons);

  for(int n=0; n!=tree->GetEntries(); ++n){
    tree->GetEntry(n);

    //if (n != 3000) continue;
    //std::cout << X << ", " << Y << ", " << Z << std::endl;
    double posSource[3]={X, Y, Z};
    double num_phot_generated = genPhotons;

    double distance_to_center = GetDistanceCenter(centerYZ, posSource[2], posSource[1]);
    //double distance_to_center_lateral = sqrt( pow(posSource[0] - x_anode, 2) + pow(posSource[2] - z_anode, 2));

    //std::cout << distance_to_center << std::endl;
    if(distance_to_center<=0) {
      std::cout << "Warning: distance to center < 0" << std::endl;
      continue;
    }

    //loop over the channels
    for(int i=0; i<numberDevices; i++) {

      // Modification requried from Users:
      bool isDouble=false;
      if (devx.at(i) > 350.){
        if (devz.at(i) < 231.){
          if ( (devy.at(i) < 35.) || ((devy.at(i) > 210.) && (devy.at(i) < 340.)) || ((devy.at(i) > 450.) && (devy.at(i) < 510.)) || (devy.at(i) > 570.) ){
            isDouble = true;
          }
        }
        if(devz.at(i) > 231.){
          if ( ((devy.at(i) > 90.) && (devy.at(i) < 150.)) || ((devy.at(i) > 210.) && (devy.at(i) < 280.)) || ((devy.at(i) > 390.) && (devy.at(i) < 450.)) || ((devy.at(i) > 510.) && (devy.at(i) < 580.)) ) {
            isDouble = true;
          }
        }
      }
      // Reminder: continue: skip the rest of the code when meeting the if-condition.
      //if (devx.at(i) > 0.) continue; // for single sided GH curves, only select the one TPC
      if (!isDouble) continue;    // remove double-sided supercells.

      int entries = VUV_hits[i];
      if(entries < min_number_entries) continue;

      double distance_to_pmt = sqrt(pow(posSource[0] - devx.at(i),2) +
                                    pow(posSource[1] - devy.at(i),2) +
                                    pow(posSource[2] - devz.at(i),2));

      // calculate theta
      double theta = -1;
      //if (abs(devy.at(i)) > 730) {
        // laterals
        //double coseno = sqrt(pow(posSource[1] - devy.at(i),2))/distance_to_pmt;
        //theta = acos(coseno)*180/3.1416;
      //}else {
        // cathode
        double coseno = sqrt(pow(posSource[0] - devx.at(i),2))/distance_to_pmt;
        theta = acos(coseno)*180/3.1416;
      //}

      if(theta >= 90 || theta == -1) {
        std::cout << "Warning -- Theta: " << theta << std::endl;
        continue;
      }

      double geo_factor = -1;
      if(IsRectangular) {
        //------Rectangle---------
        acc detPoint;
        TVector3 ScintPoint_rel;

        // centre coordinates of optical detector
        detPoint.ax = devx.at(i); detPoint.ay = devy.at(i); detPoint.az = devz.at(i);
        // width and height in cm of arapuca active window
        detPoint.w = arapuca_w; detPoint.h = arapuca_h;
        // get scintillation point coordinates relative to arapuca window centre
        TVector3 ScintPoint(posSource[0], posSource[1], posSource[2]);
        TVector3 OpDetPoint(devx.at(i), devy.at(i), devz.at(i));
        ScintPoint_rel = ScintPoint - OpDetPoint;

        //if (abs(devy.at(i)) > 730) {
          // laterals
          // calculate solid angle
          //geo_factor  = Rectangle_SolidAngle_lateral(detPoint, ScintPoint_rel);

          // altenative -- project points?
          //TVector3 ScintPointNew(posSource[1], posSource[0], posSource[2]);
          //TVector3 OpDetPointNew(devy.at(i), devx.at(i), devz.at(i));
          //ScintPoint_rel = ScintPointNew - OpDetPointNew;
          //geo_factor  = Rectangle_SolidAngle(detPoint, ScintPoint_rel);
        //}else {
          // cathode
          // calculate solid angle
          //continue;
          //TVector3 ScintPointFlip(-posSource[0], posSource[1], posSource[2]);
          //TVector3 OpDetPointFlip(-devx.at(i), devy.at(i), devz.at(i));
          //ScintPoint_rel = ScintPointFlip - OpDetPointFlip;
          geo_factor  = Rectangle_SolidAngle(detPoint, ScintPoint_rel);
        //}
      }else {
        if(IsSphere) {
          //------Semi-sphere---------
          // calculate solid angle
          geo_factor = Omega_Dome_Model(distance_to_pmt, theta);
        }
        else {
          //------Disk---------
          double d,h;
          //offset in z-y plane
          d = sqrt((posSource[2] - devz.at(i))*(posSource[2] - devz.at(i)) + (posSource[1] - devy.at(i))*(posSource[1] - devy.at(i)));
          //drift distance (in x)
          h =  sqrt((devx.at(i) - posSource[0])*(devx.at(i) - posSource[0]));
          // calculate solid angle
          geo_factor = Omega(d, h, b);
        }
      }

      //pure geometric estimation of the number of arriving VUV photons
      double rec_N =  exp(-1.*distance_to_pmt/L_abs)*gRandom->Poisson(num_phot_generated*geo_factor/(4*3.1416));

      //if (entries > 1000) std::cout << "OpDet: " << i << ", distance: " << distance_to_pmt << ", theta: " << theta << ", G4: " << entries << ", S-A: " << rec_N << ", geo = " << geo_factor << std::endl;
      v_x.push_back(posSource[0]);
      v_devx.push_back(devx.at(i));
      v_hits.push_back(entries);
      v_distance.push_back(distance_to_pmt);
      v_rec_hits.push_back(rec_N);
      v_offset_angle.push_back(theta);

      //if (abs(devy.at(i)) > 730) {
        // laterals
        // switch between two methods.
        //v_d_center.push_back(distance_to_center_lateral); // used for GH curves fitting method.
        //v_d_center.push_back(abs(x_anode - posSource[0])); // used for interpolation method.
      //}else{
        v_d_center.push_back(distance_to_center);
      //}

      //if (devy.at(i) == Y && devz.at(i) == Z){
        //std::cout <<" PMT: " <<  i <<", dev x,y,z: "<<devx.at(i)<<", "<<devy.at(i)<<", "<<devz.at(i)<<", points: ("<<X<<", "<<Y<<", "<<Z<< "), Hits: " << entries << ", Geo: " << rec_N << ", ratio (Geant4/Geo/cos): "<<entries/rec_N*(cos(theta))<<" "<<exp(1.*distance_to_pmt/L_abs)<<", Fractional differences: " << (entries-rec_N)/entries << ", Distance: " << distance_to_pmt << ", Angle: " << theta << std::endl <<std::endl;
      //}


    }// channels
  }//files
  v1 = v_distance;
  v2 = v_hits;
  v3 = v_rec_hits;
  v4 = v_offset_angle;
  v5 = v_d_center;
  v6 = v_x;
  v7 = v_devx;
}

int main(int argc, char * argv[]) {

  if (argc < 2) {
    std::cerr << "Please specify input file" << std::endl;
    return 1;
  }
  // path to the directory where the files are
  // path + file containing the LDS id and positions (x, y, z)
  std::string positions = "./protodune_optical_mapping.txt";
  // path + file containing the Sim-Info
  std::string file_name = "output";
  //std::string input_file = file_name + ".root";
  std::string input_file = argv[1];
  std::string save_dir = "./plots/";

  // If PMTs (IsSphere), if (X-)ARAPUCAS (IsRectangular),
  // acrylic disk in front of the PMTs (both false)
  bool IsRectangular = true;
  bool IsSphere = false;
  //offset-angle theta binning
  const int N = 9;
  double delta_angulo = 90./N;
  std::vector<double> angulo;

  double theta[N];
  for(int i=0; i<N; i++){
    theta[i] =  delta_angulo/2 + i*delta_angulo;
  }

  //Distance range and step for the profile to fit with GH
  double d_min = 0;
  double d_max = 700.; // changable 1400 - argon; 2000 -xenon
  double step_d = 50;

  bool isDouble=true;
  //Center distance bins
  double range_d = 400;//1100; 300;

  // Modification requried from Users:
  const int M = 6;
  double range_d_array[M+1];
  if (!isDouble) {
    double range_d_array_temp[M+1] = {0.,100.,150.,200.,250.,300.,range_d};
    std::copy(std::begin(range_d_array_temp), std::end(range_d_array_temp), std::begin(range_d_array));
  }else {
    double range_d_array_temp[M+1] = {0.,100.,150.,200.,250.,280.,range_d};
    std::copy(std::begin(range_d_array_temp), std::end(range_d_array_temp), std::begin(range_d_array));
  }
  // End.

  //double delta_d = range_d/M;
  TH1D* h=new TH1D("","",range_d, 0, range_d*1.05);

  // Initializing parameters for GH fit
  // *** User *** must play and optimise the fit limits/options to guarantee
  // the correct convergency of the fit
  std::string options[N];
  for(int k=0; k<N; k++){
    options[k] = "W0Q"; //W -> set all weights to 1
  }
  //options[6] = "W0Q";
  //options[7] = "W0Q";
  //options[8] = "W0Q";

  TF1 *GH[N][M];
  TH1D* hd_centers[M];
  TLine* line[M+1];
  for(int k=0; k < M; k++) { // M -- the bin to use for border effect.
    hd_centers[k] = new TH1D("","",M,range_d_array);
    line[k] = new TLine(range_d_array[k], 0.,range_d_array[k],2000);
    //line[k] = new TLine(k*delta_d, 0.,k*delta_d,2000);
    for(int j=0; j < N; j++) { // N == 9

      GH[j][k] =  new TF1("GH",GaisserHillas,0.,d_max,4);

      double pars_ini[4] = {1., 128., 55, -2500};
      GH[j][k]->SetParLimits(2,10,350);
      //GH[j][k]->SetParLimits(1,0,0.8);

      // xenon
      /*if (j < 7) pars_ini[3] = -5000;
      else pars_ini[3] = -3000;*/

      // argon cathode
      /*if (j < 1) pars_ini[3] = -1800;
      else if (j < 4) pars_ini[3] = -1000;
      else if (j < 5) pars_ini[3] = -500;
      else if (j < 6) pars_ini[3] = -250;
      else if (j < 7) pars_ini[3] = -200;
      else pars_ini[3] = -75;*/

      // for cathode. the tune
      /*if (j == 1 && k < 5) {
        GH[j][k]->SetParLimits(1,125,300);
        if (k == 3 || k == 4)
          GH[j][k]->SetParLimits(0,1.15,1.25);
        if (k == 5){
          //GH[j][k]->SetParLimits(0,1.06,1.25);
          GH[j][k]->SetParLimits(1,125,300);
        }
      }

      if (k == 5 && j == 0) pars_ini[3] = -180;
      //if (k == 5 && j == 1) pars_ini[3] = -150;
      if (k == 6 && j == 0) pars_ini[3] = -150;
      if (k == 6 && j == 1) pars_ini[3] = -100;
      if (k == 6 && j == 2) pars_ini[3] = -500;
      if (k == 6 && j == 3) pars_ini[3] = -300;
      */

      // argon lateral
      if (j < 1) pars_ini[3] = -275;
      else if (j < 2) pars_ini[3] = -250;
      else if (j < 4) pars_ini[3] = -200;
      else if (j < 5) pars_ini[3] = -150;
      else if (j < 6) pars_ini[3] = -100;
      else if (j < 8) pars_ini[3] = -75;
      else pars_ini[3] = -50;

      GH[j][k]->FixParameter(3,pars_ini[3]);
      GH[j][k]->SetParameters(pars_ini);
    }
  }
  // end of Initialisation.

  // Plot correction curves.
  // 2d histo, the ratio between Geant 4 and geometric predition vs. drift distances / distances
  TH2F *hist = new TH2F("", "", 100, -300, 300, 100, 0, 1);
  //TH1F *hist = new TH1F("", "", 100, 0, 1);
  TProfile* pdiff_d[N][M]; // correction curves (vs. distances. )
  TProfile* pdiff_dx[N][M]; // vs. drift direction.
  for(int j=0; j < N; j++) {
    double theta_min = j*delta_angulo;
    double theta_max = theta_min + delta_angulo;
    angulo.push_back(theta_min + delta_angulo/2);
    for(int k=0; k < M; k++) {
      pdiff_d[j][k]  = new TProfile("","", int((d_max-d_min)/step_d), d_min, d_max, "s");
      pdiff_dx[j][k] = new TProfile("","", int((d_max-d_min)/step_d), d_min, d_max, "s");
    }
  }

  std::vector<double> v_distance, v_hits, v_rec_hits, v_offset_angle, v_d_center, v_x, v_devx;
  calcula(positions, input_file, v_distance, v_hits, v_rec_hits, v_offset_angle, v_d_center, IsRectangular, IsSphere, v_x, v_devx);

  for(int i=0; i<v_distance.size(); i++){
    double costheta = cos(3.1416*v_offset_angle.at(i)/180.);
    //which angulat bin
    int j = int(v_offset_angle.at(i)/delta_angulo);

    //which "center/crown" bin
    double temp = v_d_center.at(i);
    int k = std::lower_bound(range_d_array, range_d_array+M+1, temp) - range_d_array - 1;
    //int k = int(v_d_center.at(i)/delta_d);
    if(k>=M) continue;

    pdiff_d[j][k]->Fill(v_distance.at(i), v_hits.at(i)/v_rec_hits.at(i)*costheta);
    hd_centers[k]->Fill(v_d_center.at(i));
    h->Fill(v_d_center.at(i));
    pdiff_dx[j][k]->Fill(v_distance.at(i), v_distance.at(i));
  }


//------------------------------------------------------------------------------------------------------------------------------
// Canvas 0, the crowns plot, to study for border effect.
  TCanvas *canvas0 = new TCanvas("canvas0", "graph draw options",200,200,500,400);
  h->SetTitle("Help to choose the \"border-study\" distance bins");
  h->GetXaxis()->SetTitle("distance to centre in Y-Z plane [cm]");
  h->SetStats(0);
  h->Draw("hist");
  line[M] = new TLine(range_d, 0.,range_d,2000);
  for(int k=0; k < M+1; k++) {
    line[k]->SetLineColor(2);
    line[k]->SetLineWidth(3);
    line[k]->SetLineStyle(kDashed);
    line[k]->Draw("same");
  }
  canvas0->Update();
  canvas0->Modified();
  //canvas0->WaitPrimitive();

  TString fig_name;
  if (isDouble) fig_name = save_dir + file_name+"_crown_plot_double_sided.pdf";
  else fig_name = save_dir + file_name+"_crown_plot_single_sided.pdf";
  canvas0->SaveAs(fig_name);

  //------------------------------------------------------------------------------------------------------------------------------
  //Estimate range to fit in distance
  double max_x[N][M], min_x[N][M], n_entries[N][M];
  for(int j=0; j < N; j++) {
    for(int k=0; k < M; k++) {
      n_entries[j][k] = 0;
      TAxis *xaxis  =  pdiff_d[j][k]->GetXaxis();
      Int_t nbins  = xaxis->GetNbins();
      double min = d_max;
      double max = 0;
      for (Int_t bin=0;bin<=nbins;bin++) {
        n_entries[j][k] += pdiff_d[j][k]->GetBinContent(bin);
        if(pdiff_d[j][k]->GetBinContent(bin) == 0) continue;
        if(min > xaxis->GetBinCenter(bin)) min = xaxis->GetBinCenter(bin);
        if(max < xaxis->GetBinCenter(bin)) max = xaxis->GetBinCenter(bin);
      }
      max_x[j][k] = max;
      min_x[j][k] = min;
    }
  }

  std::vector<int> N_canvas;
  std::vector<double> d_center;
  std::string title[M];
  for(int k=0; k < M; k++) {
    for(int j=0; j < N; j++){
      if(n_entries[j][k]>0){
        N_canvas.push_back(k);
        break;
      }
    }
    if(hd_centers[k]->GetEntries()>0) {
      d_center.push_back(hd_centers[k]->GetMean());
      //title[k] =  "Distance to anode " +  std::to_string(int(hd_centers[k]->GetMean())) + "cm";
      title[k] =  "Distance to center " +  std::to_string(int(hd_centers[k]->GetMean())) + "cm";
    }
  }

  // Making plots and fits ...
  const int dim = N_canvas.size();
  TGraphErrors *gr[N][M];
  TGraphErrors *grx[N][M];
  for(int k=0; k < M; k++){
    for(int j=0; j < N; j++){
      if(n_entries[j][k] > 0){
        gr[j][k] = Profile_to_Graph(pdiff_dx[j][k], pdiff_d[j][k], j); // correction curves.
        gr[j][k]->SetLineColor(1 + j%9);
        gr[j][k]->SetMarkerColor(1 + j%9);
        if(j==4) {gr[j][k]->SetLineColor(kOrange+7);gr[j][k]->SetMarkerColor(kOrange+7);}
        gr[j][k]->SetMarkerStyle(20 + j);
        gr[j][k]->SetMarkerSize(0.6);
      }
    }
  }

  // Canvas to plot correction curves.
  const int dimension = ceil(double(N_canvas.size())/2.);
  TCanvas *canvas1 = new TCanvas("canvas1", "GH",200,200,dimension*450,800);
  if(N_canvas.size() > 1) {
    double nn = double(N_canvas.size())/2.;
    canvas1->Divide(ceil(nn),2);
  }
  double x_0[2] = {0, d_max};
  double y_0[2] = {0, 2.};
  TGraph* gg0[dim];
  std::vector<double> p1[dim], p2[dim], p3[dim], p4[dim], ep1[dim], ep2[dim], ep3[dim];
  TLegend *leg1=new TLegend(0.55, 0.4, 0.89, 0.89,NULL,"brNDC");
  leg1->SetHeader("");
  leg1->SetBorderSize(0);
  char label[N][20];

  for(int l=0; l < dim; l++){
    int k = N_canvas.at(l);
    canvas1->cd(1 + l);
    gg0[l] = new TGraph(2,x_0,y_0);
    gg0[l]->SetTitle(title[k].c_str());
    gg0[l]->GetYaxis()->SetTitleSize(0.05);
    gg0[l]->GetYaxis()->SetTitleOffset(1.);
    gg0[l]->GetXaxis()->SetTitleSize(0.05);
    gg0[l]->GetXaxis()->SetTitleOffset(1.);
    gg0[l]->GetXaxis()->SetRangeUser(0,d_max);
    gg0[l]->GetYaxis()->SetRangeUser(0, 1.5); // changable.
    gg0[l]->GetXaxis()->SetTitle("distance [cm]");
    gg0[l]->GetYaxis()->SetTitle("N_{hit} / N_{#Omega} / cos(#theta)");
    gg0[l]->Draw("ap");

//Info in <ROOT::Math::ParameterSettings>: lower/upper bounds outside current parameter value. The value will be set to (low+up)/2
    for(int j=0; j < N; j++) {
      double pars_GH[4] = {-999, -999, -999, -999};
      double epars_GH[4]= {-999, -999, -999, -999};

      if(n_entries[j][k]>0) {

        gr[j][k]->Draw("p");

        // Fitting the simulation data
        gr[j][k]->Fit(GH[j][k], options[j].c_str(),"", min_x[j][k],max_x[j][k]);
        //Loading parameters
        GH[j][k]->GetParameters(pars_GH);
        GH[j][k]->SetParameters(pars_GH);
        GH[j][k]->SetLineColor(1 + j);
        if(j==4) GH[j][k]->SetLineColor(kOrange+7);
        GH[j][k]->SetLineStyle(kDashed);
        GH[j][k]->Draw("same");
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
      if(l==0) {
        int a_min = j*delta_angulo;
        int a_max = (j+1)*delta_angulo;
        sprintf(label[j],"#theta #in [%i, %i] deg",a_min, a_max);
        leg1->AddEntry(gr[j][k],label[j],"p");
      }
    }
   //if(l==0)  leg1->Draw();
  }

  canvas1->Update();
  canvas1->Modified();
  //canvas1->WaitPrimitive();

  TString fig_name2;
  if (isDouble) fig_name2 = save_dir + file_name+"_GH_double_sided.pdf";
  else fig_name2 = save_dir + file_name+"_GH_single_sided.pdf";
  canvas1->SaveAs(fig_name2);

  std::cout<<"  GH p1: "<<std::endl <<"{";
  for(int k=0; k < dim; k++){
    std::cout<<"{";
    for(int j=0; j<N; j++)
      std::cout<<p1[k].at(j)<<", ";
    std::cout<<"},"<<std::endl;
  }
  /*  std::cout<<"  ERROR GH p1: "<<std::endl;
  for(int k=0; k < dim; k++){
    for(int j=0; j<N; j++)
      std::cout<<ep1[k].at(j)<<", ";
    std::cout<<""<<std::endl;
    }*/
  std::cout<<"  GH p2: "<<std::endl<<"{";
  for(int k=0; k < dim; k++){
    std::cout<<"{";
    for(int j=0; j<N; j++)
      std::cout<<p2[k].at(j)<<", ";
    std::cout<<"},"<<std::endl;
  }
  /* std::cout<<"  ERROR GH p2: "<<std::endl;
  for(int k=0; k < dim; k++){
    for(int j=0; j<N; j++)
      std::cout<<ep2[k].at(j)<<", ";
    std::cout<<""<<std::endl;
  }*/
  std::cout<<"  GH p3: "<<std::endl<<"{";
  for(int k=0; k < dim; k++){
    std::cout<<"{";
    for(int j=0; j<N; j++)
      std::cout<<p3[k].at(j)<<", ";
    std::cout<<"},"<<std::endl;
  }
  /*  std::cout<<"  ERROR GH p3: "<<std::endl;
  for(int k=0; k < dim; k++){
    for(int j=0; j<N; j++)
      std::cout<<ep3[k].at(j)<<", ";
    std::cout<<""<<std::endl;
    }*/
  std::cout<<"  GH p4: "<<std::endl<<"{";
  for(int k=0; k < dim; k++){
    std::cout<<"{";
    for(int j=0; j<N; j++)
      std::cout<<p4[k].at(j)<<", ";
    std::cout<<"},"<<std::endl;
  }

// The calculation for border effect.
//------------------------------------------------------------------------------------------------------------------------------
  std::vector<double> vNmax[N], vdmax[N], vlambda[N], veNmax[N], vedmax[N], velambda[N], vd_center; // N -- 9.
  TGraphErrors* gNmax[N];
  TGraphErrors* gdmax[N];
  TGraphErrors* glambda[N];
  std::vector<double> slopes1, slopes2, slopes3, eslopes1, eslopes2, eslopes3, angles;
  std::vector<double> a1, a2, a3;
  TF1* f1[N];
  TF1* f2[N];
  TF1* f3[N];

  for(int j=0; j < N; j++) { // N -- 9

    for(int k=0; k < dim; k++){ // dim -- 4
      vNmax[j].push_back(p1[k].at(j));
      vdmax[j].push_back(p2[k].at(j));
      vlambda[j].push_back(p3[k].at(j));
      if(j==0) {
        vd_center.push_back(d_center.at(j));
        //std::cout<< "d_center.at(j)"<< d_center.at(j)<< std::endl;
      }
    }

    gNmax[j] =   new TGraphErrors(vd_center.size(), &d_center[0], &vNmax[j][0], 0, &veNmax[j][0]);
    gdmax[j] =   new TGraphErrors(vd_center.size(), &d_center[0], &vdmax[j][0], 0, &vedmax[j][0]);
    glambda[j] = new TGraphErrors(vd_center.size(), &d_center[0], &vlambda[j][0], 0, &velambda[j][0]);

    f1[j] =  new TF1("f1",pol1,0.,range_d,2);
    f1[j]->SetLineColor(1 + j);
    if(j==4) {f1[j]->SetLineColor(kOrange+7);}
    f1[j]->SetLineStyle(kDashed);
    f2[j] =  new TF1("f2",pol1,0.,range_d,4);
    f2[j]->SetLineColor(1 + j);
    if(j==4) {f2[j]->SetLineColor(kOrange+7);}
    f2[j]->SetLineStyle(kDashed);
    f3[j] =  new TF1("f3",pol1,0.,range_d,2);
    f3[j]->SetLineColor(1 + j);
    if(j==4) {f3[j]->SetLineColor(kOrange+7);}
    f3[j]->SetLineStyle(kDashed);

    gNmax[j]->Fit(f1[j],"Q","",0,range_d);
    gdmax[j]->Fit(f2[j],"Q","",0,range_d);
    glambda[j]->Fit(f3[j],"Q","",0,range_d);

    slopes1.push_back(f1[j]->GetParameter(1));
    slopes2.push_back(f2[j]->GetParameter(1));
    slopes3.push_back(f3[j]->GetParameter(1));
    const double* tmp1 = f1[j]->GetParErrors();
    const double* tmp2 = f2[j]->GetParErrors();
    const double* tmp3 = f3[j]->GetParErrors();
    eslopes1.push_back(tmp1[1]);
    eslopes2.push_back(tmp2[1]);
    eslopes3.push_back(tmp3[1]);

    a1.push_back(f1[j]->GetParameter(0));
    a2.push_back(f2[j]->GetParameter(0));
    a3.push_back(f3[j]->GetParameter(0));

    angles.push_back(theta[j]);

    gNmax[j]->SetLineColor(1 + j);
    gNmax[j]->SetMarkerColor(1 + j);
    gNmax[j]->SetMarkerStyle(20 + j);
    if(j==4) {gNmax[j]->SetLineColor(kOrange+7);gNmax[j]->SetMarkerColor(kOrange+7);}
    gdmax[j]->SetLineColor(1 + j);
    gdmax[j]->SetMarkerColor(1 + j);
    gdmax[j]->SetMarkerStyle(20 + j);
    if(j==4) {gdmax[j]->SetLineColor(kOrange+7);gdmax[j]->SetMarkerColor(kOrange+7);}
    glambda[j]->SetLineColor(1 + j);
    glambda[j]->SetMarkerColor(1 + j);
    glambda[j]->SetMarkerStyle(20 + j);
    if(j==4) {glambda[j]->SetLineColor(kOrange+7);glambda[j]->SetMarkerColor(kOrange+7);}
    if(j==0) {
      gNmax[j]->SetTitle(0);
      gNmax[j]->GetXaxis()->SetTitle("d_{T} [cm]");
      gNmax[j]->GetYaxis()->SetTitle("N_{max}");
      gNmax[j]->GetYaxis()->SetRangeUser(-0.5,2.0);
      gdmax[j]->SetTitle(0);
      gdmax[j]->GetXaxis()->SetTitle("d_{T} [cm]");
      gdmax[j]->GetYaxis()->SetTitle("d_{max} [cm]");
      gdmax[j]->GetYaxis()->SetRangeUser(-100,250);
      glambda[j]->SetTitle(0);
      glambda[j]->GetXaxis()->SetTitle("d_{T} [cm]");
      glambda[j]->GetYaxis()->SetTitle("#Lambda [cm]");
      glambda[j]->GetYaxis()->SetRangeUser(0,400);
    }
  }

  TGraphErrors* g1= new TGraphErrors(angles.size(),&angles[0],&slopes1[0],0, &eslopes1[0]);
  TF1* p0_m1 =  new TF1("p0_m1",pol1,0.,90,2);
  TGraphErrors* g2= new TGraphErrors(angles.size(),&angles[0],&slopes2[0],0, &eslopes2[0]);
  TF1* p0_m2 =  new TF1("p0_m2",pol1,0.,90,2);
  TGraphErrors* g3= new TGraphErrors(angles.size(),&angles[0],&slopes3[0],0, &eslopes3[0]);
  TF1* p0_m3 =  new TF1("p0_m3",pol1,0.,90,2);

  g1->SetTitle(0);
  g1->GetXaxis()->SetTitle("#theta [deg]");
  g1->GetYaxis()->SetTitle("slope N_{max} [cm^{-1}] x10^{-3}");
  g1->SetMarkerStyle(20);
  p0_m1->SetParameters(0,1);
  g1->Fit(p0_m1,"W0Q","",0,90);
  double y1 = p0_m1->GetParameter(0);
  double m1 = p0_m1->GetParameter(1);
  g2->SetTitle(0);
  g2->GetXaxis()->SetTitle("#theta [deg]");
  g2->GetYaxis()->SetTitle("slope d_{max}");
  g2->SetMarkerStyle(20);
  g2->Fit(p0_m2,"W0Q","",0,90);
  double y2 = p0_m2->GetParameter(0);
  double m2 = p0_m2->GetParameter(1);
  g3->SetTitle(0);
  g3->GetXaxis()->SetTitle("#theta [deg]");
  g3->GetYaxis()->SetTitle("slope #Lambda");
  g3->SetMarkerStyle(20);
  g3->Fit(p0_m3,"W0Q","",0,90);
  double y3 = p0_m3->GetParameter(0);
  double m3 = p0_m3->GetParameter(1);

  // Bringing the corrections to the Y-Z center
  double pGH[9][4];
  for(int i = 0; i < 9; i++) {
    //Both A) and B) options give very similar results (user should chose)
    //(A) linear fit intersection with y-axix
    pGH[i][0] = a1.at(i);
    pGH[i][1] = a2.at(i);
    pGH[i][2] = a3.at(i);

    //(B) propagation of closer point to the d_center=0
    //pGH[i][0]=  p1[0].at(i) + slopes1.at(i) * (0 - d_center.at(0));
    //pGH[i][1]=  p2[0].at(i) + slopes2.at(i) * (0 - d_center.at(0));
    //pGH[i][2]=  p3[0].at(i) + slopes3.at(i) * (0 - d_center.at(0));

    pGH[i][3]= p4[0].at(i);
  }
  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"Corrections to plug into LArSoft PhotonVisibilityServices"<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;
  std::cout<<"-------------------------------------------------------"<<std::endl;

  std::cout<<"double p1[9] = {";
  for(int i = 0; i < 8; i++)
    std::cout<<pGH[i][0]<<", ";
  std::cout<<pGH[8][0]<<"}; ";
  std::cout<<""<<std::endl;
  std::cout<<"double p2[9] = {";
  for(int i = 0; i < 8; i++)
    std::cout<<pGH[i][1]<<", ";
  std::cout<<pGH[8][1]<<"}; ";
  std::cout<<""<<std::endl;
  std::cout<<"double p3[9] = {";
  for(int i = 0; i < 8; i++)
    std::cout<<pGH[i][2]<<", ";
  std::cout<<pGH[8][2]<<"}; ";
  std::cout<<""<<std::endl;
  std::cout<<"double p4[9] = {";
  for(int i = 0; i < 8; i++)
    std::cout<<pGH[i][3]<<", ";
  std::cout<<pGH[8][3]<<"}; ";
  std::cout<<""<<std::endl;


  // Define the Canvas
  double xx[3][2];
  double yy[3][2];
  TGraphErrors* gf[3];
  for(int j=0; j < 3; j++) {
    xx[j][0] = 0;
    xx[j][1] = 1200;
    yy[j][0] = 0;
    yy[j][1] = 200;
    gf[j] = new TGraphErrors(2, xx[j], yy[j]);
  }

  TCanvas *c = new TCanvas("c", "canvas", 500, 1700);
  c->Divide(1,3);
  c->cd(1);
  TPad *pad1 = new TPad("pad1", "pad1", 0., 0.37, 1, 1.0);
  pad1->SetBottomMargin(0.1);
  pad1->Draw();
  pad1->cd();
  gf[0]->GetYaxis()->SetLabelSize(0.05);
  gf[0]->GetYaxis()->SetTitleSize(0.06);
  gf[0]->GetYaxis()->SetTitleOffset(0.74);
  gf[0]->GetYaxis()->SetRangeUser(0,1.5);    // need to rewrite this
  gf[0]->GetXaxis()->SetRangeUser(0,range_d*1.05); // need to rewrite this to sset the range of border effect. 
  gf[0]->GetXaxis()->SetLabelSize(0.05);
  gf[0]->GetXaxis()->SetTitleSize(0.05);
  gf[0]->GetXaxis()->SetTitleOffset(1.0);
  gf[0]->GetXaxis()->SetTitle("d_{T} [cm]");
  gf[0]->GetYaxis()->SetTitle("N_{max}");
  gf[0]->Draw("ap");
  for(int j=0; j < N; j++){
    gNmax[j]->Draw("p same");
    f1[j]->Draw("l same");
  }

  c->cd(1);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();
  g1->GetYaxis()->SetLabelSize(0.07);
  g1->GetYaxis()->SetTitleSize(0.08);
  g1->GetYaxis()->SetTitleOffset(0.55);
  g1->GetXaxis()->SetLabelSize(0.07);
  g1->GetXaxis()->SetTitleSize(0.08);
  g1->GetXaxis()->SetTitleOffset(0.9);
  g1->Draw("ap");
  p0_m1->SetLineColor(kGray+1);
  p0_m1->SetLineStyle(kDashed);

  c->cd(2);
  TPad *padc1 = new TPad("padc1", "padc1", 0., 0.37, 1, 1.0);
  padc1->SetBottomMargin(0.1);
  padc1->Draw();
  padc1->cd();
  gf[1]->GetYaxis()->SetLabelSize(0.05);
  gf[1]->GetYaxis()->SetTitleSize(0.06);
  gf[1]->GetYaxis()->SetTitleOffset(0.74);
  gf[1]->GetYaxis()->SetRangeUser(0,300);  // need to rewrite this
  gf[1]->GetXaxis()->SetRangeUser(0,range_d*1.05); // need to rewrite this to sset the range of border effect. 
  gf[1]->GetXaxis()->SetLabelSize(0.05);
  gf[1]->GetXaxis()->SetTitleSize(0.05);
  gf[1]->GetXaxis()->SetTitleOffset(1.0);
  gf[1]->GetXaxis()->SetTitle("d_{T} [cm]");
  gf[1]->GetYaxis()->SetTitle("d_{max} [cm]");
  gf[1]->Draw("ap");
  for(int j=0; j < N; j++){
    gdmax[j]->Draw("p same");
    f2[j]->Draw("l same");
  }

  c->cd(2);
  TPad *padc2 = new TPad("padc2", "padc2", 0, 0, 1, 0.35);
  padc2->SetTopMargin(0);
  padc2->SetBottomMargin(0.2);
  padc2->Draw();
  padc2->cd();
  g2->GetYaxis()->SetLabelSize(0.07);
  g2->GetYaxis()->SetTitleSize(0.08);
  g2->GetYaxis()->SetTitleOffset(0.55);
  g2->GetXaxis()->SetLabelSize(0.07);
  g2->GetXaxis()->SetTitleSize(0.08);
  g2->GetXaxis()->SetTitleOffset(0.9);
  g2->Draw("ap");
  p0_m2->SetLineColor(kGray+1);
  p0_m2->SetLineStyle(kDashed);

  c->cd(3);
  TPad *padcc1 = new TPad("padcc1", "padcc1", 0., 0.37, 1, 1.0);
  padcc1->SetBottomMargin(0.1);
  padcc1->Draw();
  padcc1->cd();
  gf[2]->GetYaxis()->SetLabelSize(0.05);
  gf[2]->GetYaxis()->SetTitleSize(0.06);
  gf[2]->GetYaxis()->SetTitleOffset(0.74);
  gf[2]->GetYaxis()->SetRangeUser(0,300);  // need to rewrite this
  gf[2]->GetXaxis()->SetRangeUser(0,range_d*1.05); // need to rewrite this to sset the range of border effect.
  gf[2]->GetXaxis()->SetLabelSize(0.05);
  gf[2]->GetXaxis()->SetTitleSize(0.05);
  gf[2]->GetXaxis()->SetTitleOffset(1.0);
  gf[2]->GetXaxis()->SetTitle("d_{T} [cm]");
  gf[2]->GetYaxis()->SetTitle("#Lambda [cm]");
  gf[2]->Draw("ap");
  for(int j=0; j < N; j++){
    glambda[j]->Draw("p same");
    f3[j]->Draw("l same");
    }

  c->cd(3);
  TPad *padcc2 = new TPad("padcc2", "padcc2", 0, 0, 1, 0.35);
  padcc2->SetTopMargin(0);
  padcc2->SetBottomMargin(0.2);
  padcc2->Draw();
  padcc2->cd();       // pad2 becomes the current pad
  g3->GetYaxis()->SetLabelSize(0.07);
  g3->GetYaxis()->SetTitleSize(0.08);
  g3->GetYaxis()->SetTitleOffset(0.55);
  g3->GetXaxis()->SetLabelSize(0.07);
  g3->GetXaxis()->SetTitleSize(0.08);
  g3->GetXaxis()->SetTitleOffset(0.9);
  g3->Draw("ap");
  p0_m3->SetLineColor(kGray+1);
  p0_m3->SetLineStyle(kDashed);

  TString fig_name3;
  if (isDouble) fig_name3 = save_dir + file_name+"_border_slopes_double_sided.pdf";
  else fig_name3 = save_dir + file_name+"_border_slopes_single_sided.pdf";
  c->SaveAs(fig_name3);
  /////////////////////////////////////////////////////////////////////


  std::cout<<"double slopes1[9] = {";
  for(int j=0; j < N-1; j++)
    std::cout<<slopes1[j]<<", ";
  std::cout<<slopes1[N-1]<<"};";
  std::cout<<"  "<<std::endl;
  std::cout<<"double slopes2[9] = {";
  for(int j=0; j < N-1; j++)
    std::cout<<slopes2[j]<<", ";
  std::cout<<slopes2[N-1]<<"};";
  std::cout<<"  "<<std::endl;
  std::cout<<"double slopes3[9] = {";
  for(int j=0; j < N-1; j++)
    std::cout<<slopes3[j]<<", ";
  std::cout<<slopes3[N-1]<<"};";
  std::cout<<"  "<<std::endl;

}

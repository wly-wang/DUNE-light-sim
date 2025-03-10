#include "TVector3.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TLine.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

//
// calculated semi-analytic model border corrections for visible photons
//

struct acc{
  // ax,ay,az = centre of rectangle; w = y dimension; h = z dimension
  double ax, ay, az, w, h;
};

// function prototypes
int VisHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type);
int VUVHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type, const double &cosine, const double &theta, const double distance, const int &j);

Double_t GaisserHillas(double x,double *par);
double omega(const double &a, const double &b, const double &d);
double solid(const acc& out, const TVector3 &v);
double solid_lateral(const acc& out, const TVector3 &v);
//double Disk_SolidAngle(double *x, double *p);
//double Disk_SolidAngle(double d, double h, double b);
double interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate );

double GHslope(double theta, double p0, double p1) {
  double y = p1*theta + p0;
  return y;
}
Double_t Omega_Dome_Model(double distance,double theta)
{
  // this function calculates the solid angle of a semi-sphere of radius b,
  // as a correction to the analytic formula of the on-axix solid angle,
  // as we move off-axis an angle theta. We have used 9-angular bins
  // with delta_theta width.

  // par0 = Radius correction close
  // par1 = Radius correction far
  // par2 = breaking distance betwween "close" and "far"

  double par0[9] = {0., 0., 0., 0., 0., 0.597542, 1.00872, 1.46993, 2.04221};
  double par1[9] = {0, 0, 0.19569, 0.300449, 0.555598, 0.854939, 1.39166, 2.19141, 2.57732};
  const double delta_theta = 10.;
  int j = int(theta/delta_theta);
  // 8" PMT radius
  const double b = 8*2.54/2.;
  // distance form which the model parameters break (empirical value)
  const double d_break = 5*b;//par2

  if(distance >= d_break) {
    double R_apparent_far = b - par1[j];
    return  (2*3.1416 * (1 - sqrt(1 - pow(R_apparent_far/distance,2))));

  }
  else {
    double R_apparent_close = b - par0[j];
    return (2*3.1416 * (1 - sqrt(1 - pow(R_apparent_close/distance,2))));
  }

}


// required constants:
// optical detector information
const int OpDetType = 1;  // 0 = disk, 1 = rectangular
const double y_dimension_detector = 10;    // cm
const double z_dimension_detector = 47.75; // cm
const double radius = 8*2.54/2.;      // 8" diameter PMT   // cm

// cathode foils information - SBND
// tpc drift length
const double plane_depth = 0; // cm
// size of cathode covered by foils
const double y_dimension_foils = 400; // cm
const double z_dimension_foils = 500; // cm
// centre coordinates of foils
const double x_foils = 360.279; const double y_foils = 0; const double z_foils = 679.572; // cm

// LAr properties
const double L_abs = 2000;
const double pi = 3.1416;


// VUV Gaisser-Hillas parameters
// DUNE flat surface (for cathode)
// DUNE-VD Gaisser-Hillas
  // USERS: update here with the print out from the fit!
  // Argon, RS = 99.9cm, flat PDs (Arapucas/Supercells)
  const std::vector<double> angulo = {5, 15, 25, 35, 45, 55, 65, 75, 85};
  // for DUNE FD1 HD 1x2x6 modded geometry (CPA side PDs)
  const double fGHVUVPars[4][9] = { {2.83662, 2.78267, 2.6701, 2.45801, 2.20845, 1.90297, 1.57653, 1.19622, 0.911946},
                                    {147.632, 152.908, 150.752, 150.93, 161.787, 178.607, 185.14, 177.828, 140.332},
                                    {60.6497, 92.39, 153.349, 175.23, 248.877, 271.499, 300.503, 313.338, 323.576},
                                    {-1800, -1000, -500, -500, -250, -200, -150, -150, -50}};

                                    //{-275, -250, -200, -200, -150, -100, -75, -75, -50}
  const std::vector<double> slopes1 = {0.00136662, 0.00133651, 0.00124256, 0.00111942, 0.00101549, 0.000879445, 0.00072483, 0.000559757, 0.00042345};
  const std::vector<double> slopes2 = {0.0685304, 0.0573245, 0.0623299, 0.0733342, 0.051667, 0.0366176, 0.0319569, 0.0194517, 0.0370617};
  const std::vector<double> slopes3 = {-0.0259989, -0.031503, -0.054348, -0.0765074, -0.103697, -0.117586, -0.140132, -0.138248, -0.135238};
  // double_sided
  const double fGHVUVPars_double_sided[4][9] = { {0.96311, 0.941925, 0.901074, 0.827969, 0.74171, 0.639441, 0.529001, 0.401473, 0.3031},
                                    {146.011, 149.516, 148.276, 148.898, 158.811, 171.994, 177.966, 171.344, 129.448},
                                    {190.428, 193.476, 226.5, 259.491, 293.266, 329.594, 353.936, 368.606, 363.766},
                                    {-400, -350, -300, -300, -250, -150, -100, -100, -50}};
  const std::vector<double> slopes1_double_sided = {0.000477157, 0.00046539, 0.000426561, 0.000380635, 0.000341355, 0.000297682, 0.000245961, 0.000189099, 0.000139638};
  const std::vector<double> slopes2_double_sided = {0.0550501, 0.0493876, 0.0538689, 0.0617787, 0.0477559, 0.032564, 0.0286264, 0.0184602, 0.0453299};
  const std::vector<double> slopes3_double_sided = {-0.0745289, -0.0629413, -0.0754967, -0.104142, -0.117681, -0.135749, -0.158582, -0.157273, -0.167454};

  // laterals
  const double fGHVUVPars_lateral[4][9] = { {0.96311, 0.941925, 0.901074, 0.827969, 0.74171, 0.639441, 0.529001, 0.401473, 0.3031},
                                    {146.011, 149.516, 148.276, 148.898, 158.811, 171.994, 177.966, 171.344, 129.448},
                                    {190.428, 193.476, 226.5, 259.491, 293.266, 329.594, 353.936, 368.606, 363.766},
                                    {-400, -350, -300, -300, -250, -150, -100, -100, -50}};
  const std::vector<double> slopes1_lateral = {0.000477157, 0.00046539, 0.000426561, 0.000380635, 0.000341355, 0.000297682, 0.000245961, 0.000189099, 0.000139638};
  const std::vector<double> slopes2_lateral = {0.0550501, 0.0493876, 0.0538689, 0.0617787, 0.0477559, 0.032564, 0.0286264, 0.0184602, 0.0453299};
  const std::vector<double> slopes3_lateral = {-0.0745289, -0.0629413, -0.0754967, -0.104142, -0.117681, -0.135749, -0.158582, -0.157273, -0.167454};

  // interpolate.
  /*std::vector<double> fGH_d_anode = {26., 105., 180., 254.}; // larger distances use 260cm case
  std::vector<std::vector<std::vector<double>>> fGHVUVPars_lateral_interpolate = {
    { {0.994607, 0.888712, 0.706136, 0.65617, 0.519425, 0.416587, 0.303893, 0.236757, 0.244876},
      {0.826898, 0.842618, 0.740369, 0.669461, 0.552277, 0.449452, 0.350628, 0.305361, 0.113616},
      {0.827474, 0.807339, 0.771067, 0.663827, 0.534322, 0.449413, 0.3554, 0.311696, 0.174768},
      {0.833167, 0.812198, 0.762091, 0.679776, 0.539244, 0.465872, 0.327116, 0.267314, 0.125724}
    },
    { {-81.7538, -65.5438, 2.53732, 12.0577, 54.4923, 85.4044, 110.599, 95.4083, 62.6809},
      {67.9336, 15.3034, 59.6169, 65.7805, 103.086, 133.466, 147.983, 102.703, 237.976},
      {89.657, 94.8347, 64.7458, 88.1447, 136.339, 157.599, 170.612, 121.523, 182.139},
      {91.9275, 99.1419, 85.3368, 95.3437, 145.187, 154.419, 204.796, 162.026, 254.467}
    },
    { {177.699, 200, 200, 199.998, 199.977, 200, 199.644, 199.978, 199.929},
      {145.268, 199.997, 200, 199.526, 199.999, 199.999, 199.999, 200, 199.68},
      {170.489, 156.77, 199.999, 200, 199.983, 200, 199.957, 200, 191.15},
      {187.388, 172.357, 199.998, 200, 199.991, 200, 199.977, 200, 199.799}
        },
    { {-275, -250, -200, -200, -150, -100, -75, -75, -50},
      {-275, -250, -200, -200, -150, -100, -75, -75, -50},
      {-275, -250, -200, -200, -150, -100, -75, -75, -50},
      {-275, -250, -200, -200, -150, -100, -75, -75, -50}
    }
  };
  */

  /*
  // DUNE flat surface (for cathode)
  // DUNE-VD Gaisser-Hillas
  // Argon, RS = 99.9cm, flat PDs (Arapucas/Supercells)
  const double fGHVUVPars[4][9] = { {1.2343, 1.19482, 1.13939, 1.06642, 0.957333, 0.8411, 0.721859, 0.582155, 0.420655},
                  {160.789, 163.041, 163.416, 176.419, 190.937, 205.741, 239.029, 244.506, 255.243},
                  {90.6463, 92.4193, 94.1577, 139.613, 139.729, 188.447, 200, 200, 186.966},
                  {-1000, -1000, -1000, -500, -500, -250, -100, -100, -100} };
  const std::vector<double> angulo = {5, 15, 25, 35, 45, 55, 65, 75, 85};
  const std::vector<double> slopes1 = {-1.76996e-05, -1.53516e-05, -2.32847e-05, -2.15235e-05, -1.22799e-05, -2.12407e-05, -2.28983e-05, -1.17738e-05, -9.59168e-06};
  const std::vector<double> slopes2 = {-0.0342994, -0.0421957, -0.0253359, -0.0275521, -0.0354278, -0.0183091, -0.022059, -0.0228965, -0.0109844};
  const std::vector<double> slopes3 = {-0.0149828, -0.00686603, -0.0105492, -0.00897181, -0.00182785, -0.00985137, -3.68105e-07, 1.15155e-08, -0.00467523};
  */

  // dome GH
  const double fGHVUVPars_dome[4][9] = { {0.96311, 0.941925, 0.901074, 0.827969, 0.74171, 0.639441, 0.529001, 0.401473, 0.3031},
                                    {146.011, 149.516, 148.276, 148.898, 158.811, 171.994, 177.966, 171.344, 129.448},
                                    {190.428, 193.476, 226.5, 259.491, 293.266, 329.594, 353.936, 368.606, 363.766},
                                    {-400, -350, -300, -300, -250, -150, -100, -100, -50}};
  const std::vector<double> angulo_dome = {5, 15, 25, 35, 45, 55, 65, 75, 85};
  const std::vector<double> slopes1_dome = {0.000477157, 0.00046539, 0.000426561, 0.000380635, 0.000341355, 0.000297682, 0.000245961, 0.000189099, 0.000139638};
  const std::vector<double> slopes2_dome = {0.0550501, 0.0493876, 0.0538689, 0.0617787, 0.0477559, 0.032564, 0.0286264, 0.0184602, 0.0453299};
  const std::vector<double> slopes3_dome = {-0.0745289, -0.0629413, -0.0754967, -0.104142, -0.117681, -0.135749, -0.158582, -0.157273, -0.167454};



// Visible centre correction parameters
/*
// SBND RS99cm, FLAT SURFACE
// with border Corrections
const std::vector<double> vDistances_x = {10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 195};      // cm [11]
// final
const std::vector<double> vDistances_r = {10, 60.4619, 136.421, 186.119, 234.324, 280.4, 299.1};    // cm [6]

// new model
// final
const std::vector<std::vector<std::vector<double>>> VIS_RS100_SBND_Borders = {
{ {1.29341, 1.06689, 0.917514, 0.815145, 0.753364, 0.702065, 0.666918, 0.644465, 0.621166, 0.595002, 0.588316},
  {1.31296, 1.08502, 0.937976, 0.841677, 0.772627, 0.725873, 0.691891, 0.672458, 0.651384, 0.623332, 0.618748},
  {1.43405, 1.20473, 1.05372, 0.947644, 0.873164, 0.819357, 0.782934, 0.75594, 0.728914, 0.703558, 0.695423},
  {1.59301, 1.36987, 1.19281, 1.06331, 0.966613, 0.896812, 0.845181, 0.803511, 0.774326, 0.740068, 0.7286},
  {1.78278, 1.56682, 1.3711, 1.21558, 1.09602, 1.01029, 0.9482, 0.896445, 0.854848, 0.817491, 0.803319},
  {2.25679, 2.05199, 1.77843, 1.53501, 1.33505, 1.19041, 1.07851, 0.990501, 0.929601, 0.872777, 0.859695},
  {2.46366, 2.29233, 1.92881, 1.61222, 1.38884, 1.20583, 1.07518, 1.0032, 0.922995, 0.870282, 0.85447}
},
{ {1.26185, 1.05699, 0.919836, 0.82496, 0.759135, 0.716973, 0.684797, 0.662843, 0.639668, 0.615167, 0.604396},
  {1.29166, 1.07927, 0.940517, 0.847525, 0.785951, 0.74171, 0.710172, 0.690287, 0.674349, 0.647164, 0.635898},
  {1.40541, 1.19361, 1.05497, 0.955724, 0.886237, 0.836534, 0.80311, 0.772294, 0.749522, 0.722264, 0.710622},
  {1.54293, 1.33626, 1.18082, 1.0663, 0.984753, 0.921948, 0.877609, 0.842915, 0.813647, 0.783677, 0.770201},
  {1.71983, 1.5348, 1.36362, 1.22346, 1.11845, 1.03964, 0.975172, 0.931689, 0.894116, 0.851895, 0.837597},
  {1.97572, 1.8517, 1.64545, 1.44385, 1.29247, 1.17045, 1.08171, 1.00954, 0.962294, 0.916043, 0.893014},
  {2.124, 2.08021, 1.82302, 1.56518, 1.37235, 1.24251, 1.13266, 1.03669, 0.990119, 0.933394, 0.947035}
},
{ {1.21386, 1.04246, 0.922076, 0.843702, 0.785706, 0.7451, 0.721179, 0.699071, 0.682837, 0.659844, 0.643677},
  {1.24777, 1.06882, 0.9486, 0.86897, 0.812549, 0.77522, 0.749193, 0.730754, 0.714888, 0.689542, 0.680021},
  {1.33396, 1.16049, 1.0439, 0.963931, 0.906936, 0.865632, 0.838478, 0.817183, 0.799853, 0.777907, 0.766699},
  {1.45335, 1.29304, 1.16917, 1.07182, 1.00087, 0.949423, 0.911154, 0.882191, 0.858197, 0.829756, 0.815373},
  {1.59895, 1.4705, 1.33064, 1.21224, 1.12157, 1.05265, 1.00036, 0.965166, 0.936406, 0.903511, 0.896755},
  {1.84783, 1.78559, 1.62281, 1.46062, 1.32436, 1.21775, 1.13392, 1.07288, 1.02545, 0.987848, 0.96399},
  {2.15177, 2.10504, 1.82927, 1.58502, 1.40705, 1.242, 1.1487, 1.0641, 1.00004, 0.957495, 0.950732}
},
{ {1.1727, 1.03578, 0.940702, 0.873906, 0.827511, 0.797897, 0.775774, 0.761414, 0.749438, 0.72382, 0.711179},
  {1.2132, 1.06484, 0.966703, 0.901601, 0.855773, 0.825201, 0.808204, 0.79367, 0.782621, 0.758138, 0.746555},
  {1.25039, 1.11687, 1.03042, 0.970275, 0.930286, 0.903234, 0.885374, 0.873358, 0.863847, 0.843897, 0.833992},
  {1.31485, 1.21002, 1.12804, 1.06291, 1.01489, 0.979885, 0.957625, 0.941524, 0.926478, 0.908791, 0.894145},
  {1.43499, 1.36914, 1.27599, 1.19562, 1.13189, 1.08483, 1.05269, 1.02897, 1.01285, 0.9901, 0.979448},
  {1.59552, 1.60331, 1.5141, 1.4107, 1.31642, 1.23953, 1.18134, 1.13969, 1.11084, 1.07396, 1.0626},
  {1.79348, 1.87074, 1.70682, 1.54105, 1.4167, 1.29971, 1.22392, 1.15606, 1.13197, 1.10927, 1.08277}
},
{ {1.15202, 1.04494, 0.97332, 0.919725, 0.890852, 0.869741, 0.858668, 0.853368, 0.844606, 0.824155, 0.81142},
  {1.15499, 1.04465, 0.974062, 0.928791, 0.899677, 0.883993, 0.874839, 0.873, 0.869601, 0.853134, 0.843116},
  {1.13576, 1.05186, 1.00015, 0.968852, 0.948286, 0.939056, 0.938627, 0.939464, 0.942101, 0.93046, 0.920415},
  {1.1822, 1.12938, 1.08292, 1.04746, 1.02311, 1.00959, 1.00423, 1.00404, 1.00415, 0.994619, 0.982961},
  {1.25265, 1.23815, 1.19676, 1.15935, 1.13046, 1.11252, 1.10558, 1.10391, 1.10331, 1.09726, 1.0894},
  {1.35865, 1.41995, 1.38585, 1.33898, 1.29103, 1.2456, 1.22157, 1.20555, 1.19286, 1.18111, 1.16382},
  {1.49942, 1.62699, 1.55145, 1.45453, 1.3655, 1.30066, 1.25959, 1.22223, 1.20069, 1.19475, 1.19202}
},
{ {1.22454, 1.12676, 1.05752, 1.01766, 0.993549, 0.976985, 0.967248, 0.970452, 0.967913, 0.945087, 0.937532},
  {1.11336, 1.02869, 0.976547, 0.94902, 0.935041, 0.931572, 0.937711, 0.946622, 0.956564, 0.947466, 0.935481},
  {1.05085, 0.994929, 0.968169, 0.958327, 0.960132, 0.970809, 0.99016, 1.01019, 1.02706, 1.03047, 1.02387},
  {1.06564, 1.04399, 1.02524, 1.01623, 1.01677, 1.02371, 1.04, 1.06211, 1.0812, 1.08682, 1.08115},
  {1.1062, 1.11927, 1.10719, 1.09901, 1.0981, 1.10266, 1.11982, 1.1389, 1.16453, 1.17172, 1.16851},
  {1.16715, 1.24818, 1.25329, 1.23966, 1.22603, 1.2182, 1.2223, 1.23606, 1.24656, 1.25359, 1.24914},
  {1.26056, 1.41104, 1.3769, 1.33283, 1.28912, 1.26094, 1.25143, 1.24354, 1.25471, 1.27858, 1.28191}
},
{ {1.22454, 1.12676, 1.05752, 1.01766, 0.993549, 0.976985, 0.967248, 0.970452, 0.967913, 0.945087, 0.937532},
  {1.11336, 1.02869, 0.976547, 0.94902, 0.935041, 0.931572, 0.937711, 0.946622, 0.956564, 0.947466, 0.935481},
  {1.0788, 1.02398, 1.0005, 0.997336, 1.00241, 1.02527, 1.05482, 1.08624, 1.10587, 1.12084, 1.10934},
  {1.05183, 1.03192, 1.01945, 1.01935, 1.02909, 1.05124, 1.0776, 1.11282, 1.14536, 1.1616, 1.15326},
  {1.04848, 1.06353, 1.05832, 1.05676, 1.07002, 1.08549, 1.11282, 1.15114, 1.18238, 1.21159, 1.20885},
  {1.0723, 1.14584, 1.16101, 1.15326, 1.15396, 1.1605, 1.1807, 1.20181, 1.24273, 1.26176, 1.26905},
  {1.14896, 1.27505, 1.25874, 1.21844, 1.19258, 1.18584, 1.19203, 1.21234, 1.23643, 1.28044, 1.29723}
},
{ {1.22454, 1.12676, 1.05752, 1.01766, 0.993549, 0.976985, 0.967248, 0.970452, 0.967913, 0.945087, 0.937532},
  {1.11336, 1.02869, 0.976547, 0.94902, 0.935041, 0.931572, 0.937711, 0.946622, 0.956564, 0.947466, 0.935481},
  {1.0788, 1.02398, 1.0005, 0.997336, 1.00241, 1.02527, 1.05482, 1.08624, 1.10587, 1.12084, 1.10934},
  {1.05183, 1.03192, 1.01945, 1.01935, 1.02909, 1.05124, 1.0776, 1.11282, 1.14536, 1.1616, 1.15326},
  {1.04848, 1.06353, 1.05832, 1.05676, 1.07002, 1.08549, 1.11282, 1.15114, 1.18238, 1.21159, 1.20885},
  {1.0723, 1.14584, 1.16101, 1.15326, 1.15396, 1.1605, 1.1807, 1.20181, 1.24273, 1.26176, 1.26905},
  {1.27474, 1.44297, 1.3349, 1.24249, 1.27868, 1.31034, 1.26115, 1.29498, 1.35899, 1.32796, 1.39645}
},
{ {1.22454, 1.12676, 1.05752, 1.01766, 0.993549, 0.976985, 0.967248, 0.970452, 0.967913, 0.945087, 0.937532},
  {1.11336, 1.02869, 0.976547, 0.94902, 0.935041, 0.931572, 0.937711, 0.946622, 0.956564, 0.947466, 0.935481},
  {1.0788, 1.02398, 1.0005, 0.997336, 1.00241, 1.02527, 1.05482, 1.08624, 1.10587, 1.12084, 1.10934},
  {1.05183, 1.03192, 1.01945, 1.01935, 1.02909, 1.05124, 1.0776, 1.11282, 1.14536, 1.1616, 1.15326},
  {1.04848, 1.06353, 1.05832, 1.05676, 1.07002, 1.08549, 1.11282, 1.15114, 1.18238, 1.21159, 1.20885},
  {1.0723, 1.14584, 1.16101, 1.15326, 1.15396, 1.1605, 1.1807, 1.20181, 1.24273, 1.26176, 1.26905},
  {1.27474, 1.44297, 1.3349, 1.24249, 1.27868, 1.31034, 1.26115, 1.29498, 1.35899, 1.32796, 1.39645}
}
};
*/
// SBND DOME SURFACE
// with border Corrections
const std::vector<double> vDistances_x = {10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 195, 199};   // cm [12]
// final
const std::vector<double> vDistances_r = {61.2, 138.5, 186.9, 233.2, 278.8};    // cm [6]

// new model
// final
const std::vector<std::vector<std::vector<double>>> VIS_RS100_SBND_Borders = {
{ {1.8763, 1.52373, 1.3242, 1.19976, 1.1151, 1.05578, 1.00907, 0.982144, 0.937224, 0.865445, 0.839462, 0.79318},
  {1.93968, 1.67467, 1.48614, 1.35155, 1.25265, 1.18543, 1.13355, 1.092, 1.05169, 0.995562, 0.975762, 0.955078},
  {2.2355, 1.91302, 1.67642, 1.51114, 1.38772, 1.30335, 1.22847, 1.17965, 1.12251, 1.05193, 1.01545, 0.985131},
  {2.44016, 2.17117, 1.92324, 1.73196, 1.58602, 1.47247, 1.3801, 1.30958, 1.24798, 1.16101, 1.13561, 1.10305},
  {3.14347, 2.78681, 2.40504, 2.09776, 1.8262, 1.65162, 1.50293, 1.39996, 1.3093, 1.18397, 1.17889, 1.12609}
},
{ {1.80745, 1.48535, 1.3008, 1.1898, 1.1061, 1.0493, 1.01484, 0.976277, 0.944248, 0.882041, 0.845325, 0.803434},
  {1.91318, 1.65796, 1.47698, 1.35147, 1.25976, 1.19285, 1.1433, 1.10679, 1.05892, 0.993415, 0.97606, 0.94869},
  {2.08913, 1.84576, 1.64763, 1.50504, 1.39226, 1.31155, 1.25001, 1.19946, 1.14859, 1.08283, 1.0575, 1.0335},
  {2.34601, 2.09526, 1.8806, 1.71392, 1.57715, 1.47382, 1.3906, 1.32411, 1.26456, 1.17708, 1.1429, 1.10686},
  {2.79837, 2.53838, 2.24301, 1.97426, 1.76432, 1.60871, 1.48067, 1.38432, 1.30774, 1.20871, 1.1756, 1.15198}
},
{ {1.77248, 1.48424, 1.31618, 1.20816, 1.13979, 1.08681, 1.05965, 1.02808, 0.994862, 0.921168, 0.892417, 0.843765},
  {1.83835, 1.63711, 1.48044, 1.36951, 1.29096, 1.234, 1.19044, 1.15832, 1.1188, 1.06318, 1.04285, 1.01846},
  {2.02176, 1.8108, 1.63644, 1.50929, 1.41091, 1.33912, 1.28675, 1.24143, 1.19999, 1.12673, 1.10533, 1.0717},
  {2.14757, 1.9664, 1.79718, 1.66008, 1.55188, 1.4671, 1.40252, 1.34688, 1.29602, 1.219, 1.19166, 1.15076},
  {2.59044, 2.41116, 2.17231, 1.95319, 1.77225, 1.63334, 1.52876, 1.44672, 1.37776, 1.28448, 1.24216, 1.20147}
},
{ {1.69775, 1.44526, 1.29862, 1.20978, 1.1489, 1.10433, 1.08007, 1.05234, 1.01887, 0.955227, 0.92068, 0.874745},
  {1.68841, 1.52384, 1.40315, 1.31799, 1.25804, 1.21605, 1.18363, 1.16103, 1.12936, 1.07502, 1.05167, 1.02826},
  {1.81104, 1.64695, 1.51421, 1.4171, 1.34749, 1.29791, 1.2612, 1.22904, 1.19453, 1.1359, 1.11008, 1.07742},
  {1.87606, 1.76023, 1.64998, 1.55816, 1.48252, 1.42679, 1.38541, 1.34829, 1.31035, 1.24114, 1.21324, 1.17937},
  {2.1958, 2.09713, 1.95069, 1.79933, 1.66706, 1.56632, 1.49174, 1.42973, 1.37534, 1.29159, 1.25829, 1.22535}
},
{ {1.5648, 1.35949, 1.25196, 1.18347, 1.14639, 1.11847, 1.10388, 1.08959, 1.07376, 1.0065, 0.972823, 0.926972},
  {1.46672, 1.35521, 1.27626, 1.22487, 1.19383, 1.17346, 1.16221, 1.15325, 1.13735, 1.08996, 1.07493, 1.05066},
  {1.54505, 1.44907, 1.36379, 1.31018, 1.27189, 1.24588, 1.23105, 1.21974, 1.20072, 1.15323, 1.13002, 1.09795},
  {1.5994, 1.53857, 1.47581, 1.42281, 1.38406, 1.35484, 1.33593, 1.32112, 1.30291, 1.2447, 1.2173, 1.18612},
  {1.90737, 1.87168, 1.77188, 1.67611, 1.58554, 1.52335, 1.46293, 1.43504, 1.40166, 1.33204, 1.29907, 1.25828}
},
{ {1.5241, 1.33631, 1.24821, 1.19466, 1.16756, 1.15967, 1.15183, 1.1454, 1.13249, 1.07067, 1.03428, 0.98441},
  {1.31055, 1.22337, 1.17281, 1.14424, 1.13092, 1.13168, 1.1358, 1.14018, 1.13582, 1.09974, 1.08056, 1.0597},
  {1.31779, 1.25705, 1.20763, 1.18268, 1.16782, 1.1641, 1.16677, 1.17098, 1.16924, 1.13356, 1.11169, 1.08408},
  {1.34233, 1.31486, 1.28581, 1.26407, 1.25103, 1.24717, 1.24907, 1.25334, 1.25111, 1.21344, 1.19085, 1.15777},
  {1.48501, 1.48441, 1.44596, 1.3993, 1.36101, 1.32453, 1.31428, 1.30644, 1.29399, 1.25435, 1.22196, 1.19599}
},
{ {1.5241, 1.33631, 1.24821, 1.19466, 1.16756, 1.15967, 1.15183, 1.1454, 1.13249, 1.07067, 1.03428, 0.98441},
  {1.34662, 1.27868, 1.23111, 1.21245, 1.20235, 1.20926, 1.22302, 1.24028, 1.23331, 1.19001, 1.18624, 1.15856},
  {1.32141, 1.25591, 1.207, 1.195, 1.18634, 1.19783, 1.20765, 1.2252, 1.22652, 1.19729, 1.17014, 1.13967},
  {1.22616, 1.2038, 1.1849, 1.17186, 1.17342, 1.18063, 1.19334, 1.20904, 1.21127, 1.18663, 1.16648, 1.13347},
  {1.31127, 1.31864, 1.28143, 1.24918, 1.22613, 1.20665, 1.20743, 1.20787, 1.21325, 1.1803, 1.15383, 1.13478}
},
{ {1.5241, 1.33631, 1.24821, 1.19466, 1.16756, 1.15967, 1.15183, 1.1454, 1.13249, 1.07067, 1.03428, 0.98441},
  {1.34662, 1.27868, 1.23111, 1.21245, 1.20235, 1.20926, 1.22302, 1.24028, 1.23331, 1.19001, 1.18624, 1.15856},
  {1.32141, 1.25591, 1.207, 1.195, 1.18634, 1.19783, 1.20765, 1.2252, 1.22652, 1.19729, 1.17014, 1.13967},
  {1.22616, 1.2038, 1.1849, 1.17186, 1.17342, 1.18063, 1.19334, 1.20904, 1.21127, 1.18663, 1.16648, 1.13347},
  {1.31127, 1.31864, 1.28143, 1.24918, 1.22613, 1.20665, 1.20743, 1.20787, 1.21325, 1.1803, 1.15383, 1.13478}
},
{ {1.5241, 1.33631, 1.24821, 1.19466, 1.16756, 1.15967, 1.15183, 1.1454, 1.13249, 1.07067, 1.03428, 0.98441},
  {1.34662, 1.27868, 1.23111, 1.21245, 1.20235, 1.20926, 1.22302, 1.24028, 1.23331, 1.19001, 1.18624, 1.15856},
  {1.32141, 1.25591, 1.207, 1.195, 1.18634, 1.19783, 1.20765, 1.2252, 1.22652, 1.19729, 1.17014, 1.13967},
  {1.22616, 1.2038, 1.1849, 1.17186, 1.17342, 1.18063, 1.19334, 1.20904, 1.21127, 1.18663, 1.16648, 1.13347},
  {1.31127, 1.31864, 1.28143, 1.24918, 1.22613, 1.20665, 1.20743, 1.20787, 1.21325, 1.1803, 1.15383, 1.13478}
}
};

// fit settings
// angular bin size, deg
const int number_angle_bins = 9;
const double delta_angulo = 90. / number_angle_bins;

// range and step of profiling
const double d_min = 0;
const double d_max = 365.;
const double step_d = 1;
// optical detector list
const std::string pmtlistname = "./protodune_optical_mapping.txt"; // DUNE modified geometry

// input filename
std::string inputfilename = "./output.root";
//std::string inputfilename = "./data/semi_analytic_hits_vd_negative_x.root";
//semi_analytic_hits_vd_negative_x.root semi_analytic_hits_vd_protoDUNEhd semi_analytic_hits_vd_postive_x.root


// function to calculate corrections for borders
//void vuv_rec_profile_argon(){
int main(int argc, char * argv[]){

  if (argc < 2) {
    std::cerr << "Please specify input file" << std::endl;
    return 1;
  }
  inputfilename = argv[1];

  gRandom->SetSeed(0);

  // load optical detector positions from text file
  std::vector<std::vector<double>> optical_detector_positions;
  std::cout << "Loading Photon Detector positions...\n";
  std::ifstream myfile;
  myfile.open(pmtlistname);
  if(myfile.is_open()) {
    std::cout << "File opened successfully" << std::endl;
    while(!myfile.eof()) {
        double num_pmt, x_pmt, y_pmt, z_pmt;
        if(myfile >> num_pmt >> x_pmt >> y_pmt >> z_pmt) {
            //if (x_pmt == 211.465) { // keep only pmts in TPC 1 (positive x)
              std::vector<double> line_data({num_pmt, x_pmt, y_pmt, z_pmt});
              optical_detector_positions.push_back(line_data);
          //}
        }else{ break; }
    }
    myfile.close();
    std::cout << "Photon detector positions loaded: " << optical_detector_positions.size() << std::endl << std::endl;
  }else std::cout << "Optical detector list not found." << std::endl;

  int numberPMTs = optical_detector_positions.size();

  std::cout<<"----> number optical detectors: " << numberPMTs << std::endl;

  // vectors to store calulated values
  std::vector<double> v_distance; v_distance.reserve(1e6);
  std::vector<double> v_hits_sim; v_hits_sim.reserve(1e6);
  std::vector<double> v_hits_geo; v_hits_geo.reserve(1e6);
  std::vector<double> v_offset_angle; v_offset_angle.reserve(1e6);
  std::vector<double> v_r; v_r.reserve(1e6);
  std::vector<double> v_prop_dist; v_prop_dist.reserve(1e6);

  // open file
  TFile* f = new TFile(inputfilename.c_str());
  TTree *tree = (TTree *)f->Get("myTree");
  const Int_t kMaxDevices = 960;
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

  // loop through TTree
  TH2F* h_truePE_vs_predPE = new TH2F("h_truePE_vs_predPE", "", 50, 0, 1.6e6, 50, 0, 1.6e6);
  for(int n=0; n < tree->GetEntries(); n++) {

    tree->GetEntry(n);
    double posSource[3] = {X, Y, Z};
    // check if within range of chosen parameter set
    double R = sqrt( pow(posSource[1],2) + pow(posSource[2] - z_foils, 2) );
    bool debug = false;
    if (debug) std::cout << "Scintillation point --> (x, y, z) = " << posSource[0] << "  " << posSource[1] << "  " << posSource[2] << std::endl;

    // ****
    // Calculate amount of light:
    // ****

    // loop through optical channels
    //if (posSource[0] > 0) continue;  // for single sided GH curves, only select the one TPC
    double total_pe_truth(0), total_pe_prediction(0);
    // double total_pe_prediction = 0;
    // double total_pe_truth = 0;
    for (int nPMT = 0; nPMT < numberPMTs; nPMT++) {
      //if (nPMT>479) continue;

      //std::cout << "channel started" << std::endl;

      // get Scintpoint
      TVector3 ScintPoint(posSource[0],posSource[1],posSource[2]);
      // get OpDetPoint
      int pmt_index = optical_detector_positions.at(nPMT).at(0);
      TVector3 OpDetPoint(optical_detector_positions.at(nPMT).at(1),optical_detector_positions.at(nPMT).at(2),optical_detector_positions.at(nPMT).at(3));

      // orientation hack
      int op_channel_orientation = 0;
      // if (abs(OpDetPoint[1]) > 730) op_channel_orientation = 1; //lateral
      // else op_channel_orientation = 0;

      // double vs. single sided arapucas hack
      bool isDouble=false; //This gets rid of the double-sided XARPUCAS design, maybe? Uncomment below to reinclude
      /*if (OpDetPoint[0] > 350.){
        if (OpDetPoint[2] < 231.){
          if ( (OpDetPoint[1] < 35.) || ((OpDetPoint[1] > 210.) && (OpDetPoint[1] < 340.)) || ((OpDetPoint[1] > 450.) && (OpDetPoint[1] < 510.)) || (OpDetPoint[1] > 570.) ){
            isDouble = true;
          }
        }
        if(OpDetPoint[2] > 231.){
          if ( ((OpDetPoint[1] > 90.) && (OpDetPoint[1] < 150.)) || ((OpDetPoint[1] > 210.) && (OpDetPoint[1] < 280.)) || ((OpDetPoint[1] > 390.) && (OpDetPoint[1] < 450.)) || ((OpDetPoint[1] > 510.) && (OpDetPoint[1] < 580.)) ) {
            isDouble = true;
          }
        }
      }*/

      // for DUNE-VD
      //if (op_channel_orientation == 1) continue;  // remove lateral
      //if (op_channel_orientation == 0) continue;  // remove cathode.
      //if ( abs(posSource[1] - optical_detector_positions.at(nPMT).at(2) ) < 57)  continue;// remove the scintillation points outside the field cage.
      //if (abs(ScintPoint[2] - 1000) < 750 && abs(ScintPoint[2] - 1000) > 1000) continue;

      // for protoDUNE-hd
      //if ( OpDetPoint[0] > 0 ) continue;
      //if (!isDouble) continue; // remove non double shifted.

      //calculate cosine, and the angle, this is for us to be able to plot only specific offset bins, in the rms/bias plot.
      double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
      double cosine;
      cosine = std::abs(ScintPoint[0] - OpDetPoint[0]) / distance;
      double theta = acos(cosine)*180./pi;
      //now the offset angle bin
      int j = theta/delta_angulo;

      // solid angle prediction
      double nPhotons_solid = VUVHits(genPhotons,ScintPoint,OpDetPoint,OpDetType,cosine,theta,distance,j);
      double distance_vuv = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));

      if (VUV_hits[pmt_index] < 100) continue;

      //if (j!=1) continue;
      total_pe_truth += VUV_hits[nPMT];
      total_pe_prediction += nPhotons_solid;
      

    //   // distance and angle between ScintPoint and OpDetPoint
    //   double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
    //   double cosine;
    //   // anode/cathode detectors -- fixed in x dimension
    //   if (op_channel_orientation == 0) cosine = std::abs(ScintPoint[0] - OpDetPoint[0]) / distance;
    //   // laterals -- fixed in y dimension
    //   else if (op_channel_orientation == 1) cosine = std::abs(ScintPoint[1] - OpDetPoint[1]) / distance;
    //   else {
    //     std::cout << "Error: Invalid optical detector orientation." << std::endl;
    //     exit(1);
    //   }
    // double theta = acos(cosine)*180./pi;
    //if (pmt_index==28 && ScintPoint[0]>89 && ScintPoint[0]<90 && ScintPoint[1]>510  && ScintPoint[1]<520)
    /*if (abs(VUV_hits[pmt_index]-nPhotons_solid)/VUV_hits[pmt_index]>0.2 && pmt_index==1){
      std::cout<<pmt_index<<", "<<theta<<", "<<distance_vuv<<", "<<ScintPoint[0]<<", "<<ScintPoint[1]<<", "<<ScintPoint[2]<<", "<<OpDetPoint[0]<<", "<<OpDetPoint[1]<<", "<<OpDetPoint[2]<<", "<<VUV_hits[pmt_index]<<", "<<nPhotons_solid<<", "<<(VUV_hits[pmt_index]-nPhotons_solid)/VUV_hits[pmt_index]<<std::endl;
    }*/
    // save values
    //v_distance.push_back(x_distance);     // distance, x position
    //v_offset_angle.push_back(theta);      // offset angle
    v_hits_sim.push_back(VUV_hits[pmt_index]);  // full simulation visible hits
    v_hits_geo.push_back(nPhotons_solid);   // solid angle visible hits
    v_r.push_back(R);             // radial distance from centre
    // propagation_distance
    v_prop_dist.push_back(distance_vuv);
    //v_prop_dist.push_back(325.01 - ScintPoint[0]);
    } // end loop through optical channels
    //std::cout<< total_pe_prediction << " " << total_pe_truth << std::endl;
    h_truePE_vs_predPE->Fill(total_pe_truth, total_pe_prediction);
  } // end of loop over points

  // Validation plots!
  //TString save_dir = "./validation_plots/";
  //gSystem->Exec("mkdir -p " + save_dir);

  auto ctest = new TCanvas("c1","",200,10,700,500);

  //TString type = "lateral"; //cathode lateral
  double chosen_x = 946;
  auto line1 = new TLine(0,0,1.6e6, 1.6e6);
  gPad->SetLogz();
  h_truePE_vs_predPE->Draw("colz");
  line1->SetLineColor(kRed);
  line1->Draw("same");
  //h_truePE_vs_predPE->SetTitle("Argon "+type);
  h_truePE_vs_predPE->SetTitle(" ");
  h_truePE_vs_predPE->GetXaxis()->SetTitle("total PE truth");
  h_truePE_vs_predPE->GetYaxis()->SetTitle("total PE prediction");
  h_truePE_vs_predPE->SetStats(0);
  TString name=h_truePE_vs_predPE->GetName();
  //ctest->SaveAs(save_dir+"Argon_"+type+"_"+file_name+"2d.pdf");
  ctest->SaveAs("plots/Argon_lateral_protoDUNEhd_2d.pdf");
  // create profile
  // profile
  //int n_bins = 15;
  //int bin_min = 0;
  //int bin_max = 365;

  int n_bins = 20;
  int bin_min = 0;
  int bin_max = 1000;

  double binsize = (bin_max-bin_min)/n_bins;
  TProfile* profile = new TProfile("","", n_bins, bin_min, bin_max, "s");

  // populate profile
  for(int i=0; i < v_prop_dist.size(); i++) {
    double discrepancy = (v_hits_geo[i] - v_hits_sim[i]) / v_hits_sim[i];
    double weight = v_hits_sim[i];
    profile->Fill(v_prop_dist[i],discrepancy, weight);
    //profile->Fill(v_r[i],discrepancy, weight);
  }

  // extract points from profile
  std::vector<double> r_values, r_values_error, mean, rms;
  for(int i = 1; i <= n_bins; i++) {
      // skip empty bins
      if(!(profile->GetBinEntries(i) > 0)) {
        std::cout << "bin: " << i << " empty\n";
        continue;
      }
      std::cout << "bin: " << i << " " << profile->GetBinEntries(i) << std::endl;

      //if (i > 11) continue;     // cut out underpopulated last bins
      // extract values
      mean.push_back(profile->GetBinContent(i));
      rms.push_back(profile->GetBinError(i));
      r_values.push_back(((i-1)+0.5)*binsize);
      //std::cout << i << ", " << binsize << std::endl;
      r_values_error.push_back(0.5*binsize);
  }

  // print values
  std::cout << std::endl;
  std::cout << "R: ";
  for (int i = 0; i < r_values.size(); i++) {
    std::cout << r_values[i] << ", ";
  }
  std::cout << std::endl;
  std::cout << "eR: ";
  for (int i = 0; i < r_values.size(); i++) {
    std::cout << r_values_error[i] << ", ";
  }
  std::cout << std::endl;
  std::cout << "RMS: ";
  for (int i = 0; i < r_values.size(); i++) {
    std::cout << rms[i] << ", ";
  }
  std::cout << std::endl;
  std::cout << "Bias: ";
  for (int i = 0; i < r_values.size(); i++) {
    std::cout << mean[i] << ", ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  // plot
  TCanvas *c1 = new TCanvas("c1","",200,10,1080,1080);
  c1->SetGrid();
  c1->SetBottomMargin(0.105);
  c1->SetLeftMargin(0.105);


  //TGraph *gr1 = new TGraphErrors(r_values.size(), &(r_values[0]),&(rms[0]),0, 0);
  //TGraph *gr2 = new TGraphErrors(r_values.size(), &(r_values[0]),&(mean[0]),0, 0);
  TGraph *gr1 = new TGraphErrors(r_values.size(), &(r_values[0]),&(rms[0]),&(r_values_error[0]), 0);
  TGraph *gr2 = new TGraphErrors(r_values.size(), &(r_values[0]),&(mean[0]),&(r_values_error[0]), 0);

  //gr1->GetXaxis()->SetTitle("x [cm]");

  TString type = "single-sided supercells";//double-sided supercells; single-sided supercells
  //TString type = "X-ARAPUCA";//"double-shift WLS light guides"; X-ARAPUCA

  gr1->GetXaxis()->SetTitle("distance [cm]");
  //gr1->GetXaxis()->SetTitle("d_{T} [cm]");
  gr1->GetYaxis()->SetTitle("N_{#gamma} - N_{Geant4} / N_{Geant4}");
  gr1->SetTitle("protoDUNE-hd: "+type);

  gr1->GetXaxis()->SetRangeUser(0,1000);
  //gr1->GetXaxis()->SetRangeUser(0,365);
  gr1->GetYaxis()->SetRangeUser(-0.5,0.5);


  //gg->GetYaxis()->SetLabelSize(0.05);
  gr1->GetYaxis()->SetTitleSize(0.05);
  gr1->GetYaxis()->SetTitleOffset(1.);

  //gg->GetXaxis()->SetLabelSize(0.05);
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->GetXaxis()->SetTitleOffset(1.);



  gr1->SetMarkerStyle(kFullSquare);
  gr1->SetMarkerSize(1.8);
  gr1->SetMarkerColor(1);

  gr2->SetMarkerStyle(kFullCircle);
  gr2->SetMarkerSize(1.8);
  gr2->SetMarkerColor(2);
  c1->Update();
  gr1->Draw("AP");
  c1->Update();
  gr2->Draw("same P");

  /*
  gr1->SetMarkerStyle(0);
  gr1->SetMarkerColor(0);
  gr1->SetFillColor(1);
  gr1->SetFillStyle(3001);
  gr1->Draw("A F same");


  gr2->SetMarkerStyle(0);
  gr2->SetMarkerColor(0);
  gr2->SetFillColor(2);
  gr2->SetFillStyle(3001);
  gr2->Draw("F");
  */

  TLegend *legend=new TLegend(0.71,0.71,0.96,0.96,NULL,"brNDC");
  //legend->AddEntry(gr1,"RMS","F");
  //legend->AddEntry(gr2,"Bias","F");
  legend->AddEntry(gr1,"RMS","P");
  legend->AddEntry(gr2,"Bias","P");
  legend->SetHeader("");
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  c1->Update();
  legend->Draw("same");

  c1->SaveAs("plots/Bias_"+type+".pdf");
}


// VUV hits calculation
int VUVHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type, const double &cosine, const double &theta, const double distance, const int &j) {

  // orientation hack
  int op_channel_orientation = 0;
  // if (abs(OpDetPoint[1]) > 730) op_channel_orientation = 1; // lateral.
  // else op_channel_orientation = 0; // cathode

  // double vs. single sided arapucas hack
  bool isDouble=false; //This gets rid of the double-sided XARPUCAS design, maybe? Uncomment below to reinclude
  /*if (OpDetPoint[0] > 350.){
    if (OpDetPoint[2] < 231.){
      if ( (OpDetPoint[1] < 35.) || ((OpDetPoint[1] > 210.) && (OpDetPoint[1] < 340.)) || ((OpDetPoint[1] > 450.) && (OpDetPoint[1] < 510.)) || (OpDetPoint[1] > 570.) ){
        isDouble = true;
      }
    }
    if(OpDetPoint[2] > 231.){
      if ( ((OpDetPoint[1] > 90.) && (OpDetPoint[1] < 150.)) || ((OpDetPoint[1] > 210.) && (OpDetPoint[1] < 280.)) || ((OpDetPoint[1] > 390.) && (OpDetPoint[1] < 450.)) || ((OpDetPoint[1] > 510.) && (OpDetPoint[1] < 580.)) ) {
        isDouble = true;
      }
    }
  }*/

  // distance and angle between ScintPoint and OpDetPoint
  //double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
  //double cosine;
  // anode/cathode detectors -- fixed in x dimension
  //if (op_channel_orientation == 0) cosine = std::abs(ScintPoint[0] - OpDetPoint[0]) / distance;
  // laterals -- fixed in y dimension
  //else if (op_channel_orientation == 1) cosine = std::abs(ScintPoint[1] - OpDetPoint[1]) / distance;
  // else {
  //   std::cout << "Error: Invalid optical detector orientation." << std::endl;
  //   exit(1);
  // }
  //double theta = acos(cosine)*180./pi;

  // calculate solid angle:
  double solid_angle = 0;
  // rectangular aperture
  if (optical_detector_type == 1) {
    // set Arapuca geometry struct for solid angle function
    acc detPoint;
    detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
    detPoint.w = y_dimension_detector; detPoint.h = z_dimension_detector; // width and height in cm of arapuca active window

    // get scintillation point coordinates relative to arapuca window centre
    TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;

    // calculate solid angle
    // anode/cathode detectors
    if (op_channel_orientation == 0) solid_angle = solid(detPoint, ScintPoint_rel);
    else if (op_channel_orientation == 1) solid_angle = solid_lateral(detPoint, ScintPoint_rel);
  }
  // disk aperture
  else if (optical_detector_type == 0) {
    // offset in z-y plane
    double d = sqrt(pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
    // drift distance (in x)
    double h =  sqrt(pow(ScintPoint[0] - OpDetPoint[0],2));
    // Solid angle of a disk
    //solid_angle = Disk_SolidAngle(d, h, radius);
  }
  // dome aperture
  else if (optical_detector_type == 2){
    solid_angle = Omega_Dome_Model(distance, theta);
  }
  else {
    std::cout << "Error: Invalid optical detector type." << std::endl;
    exit(1);
  }

  // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
  double hits_geo = exp(-1.*distance/L_abs) * (solid_angle / (4*pi)) * Nphotons_created;

  // determine Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence, accounting for border effects
  // offset angle bin

  // int j = (theta/delta_angulo);
  // distance from center for border corrections (cathode), and from anode (laterals)
  double r_distance = 0;
  if (op_channel_orientation == 0)  r_distance = sqrt( pow(ScintPoint[1] - y_foils, 2) + pow(ScintPoint[2] - z_foils, 2)); 
  else if (op_channel_orientation == 1) r_distance = abs(325.01 - ScintPoint[0]);
  // GH parameters
  double pars_ini[4] = {0,0,0,0};
  double s1 = 0; double s2 = 0; double s3 = 0;
  if (op_channel_orientation == 0) {

    // determine initial parameters and border corrections by optical detector type
    // flat PDs
    // dome PDs
    if (optical_detector_type == 0) {
      pars_ini[0] = fGHVUVPars_dome[0][j];
      pars_ini[1] = fGHVUVPars_dome[1][j];
      pars_ini[2] = fGHVUVPars_dome[2][j];
      pars_ini[3] = fGHVUVPars_dome[3][j];
      s1 = interpolate( angulo_dome, slopes1_dome, theta, true);
      s2 = interpolate( angulo_dome, slopes2_dome, theta, true);
      s3 = interpolate( angulo_dome, slopes3_dome, theta, true);
    }else if (optical_detector_type == 1){
      if (isDouble){
        pars_ini[0] = fGHVUVPars_double_sided[0][j];
        pars_ini[1] = fGHVUVPars_double_sided[1][j];
        pars_ini[2] = fGHVUVPars_double_sided[2][j];
        pars_ini[3] = fGHVUVPars_double_sided[3][j];
        s1 = interpolate( angulo, slopes1_double_sided, theta, true);
        s2 = interpolate( angulo, slopes2_double_sided, theta, true);
        s3 = interpolate( angulo, slopes3_double_sided, theta, true);
        //std::cout<<"here2. "<<std::endl;
      }else{
        pars_ini[0] = fGHVUVPars[0][j];
        pars_ini[1] = fGHVUVPars[1][j];
        pars_ini[2] = fGHVUVPars[2][j];
        pars_ini[3] = fGHVUVPars[3][j];
        s1 = interpolate( angulo, slopes1, theta, true);
        s2 = interpolate( angulo, slopes2, theta, true);
        s3 = interpolate( angulo, slopes3, theta, true);
      }
    }else {
      std::cout << "Error: Invalid optical detector type." << std::endl;
      exit(1);
    }
  }else if (op_channel_orientation == 1) {
    pars_ini[0] = fGHVUVPars_lateral[0][j];
    pars_ini[1] = fGHVUVPars_lateral[1][j];
    pars_ini[2] = fGHVUVPars_lateral[2][j];
    pars_ini[3] = fGHVUVPars_lateral[3][j];
    s1 = interpolate( angulo, slopes1_lateral, theta, true);
    s2 = interpolate( angulo, slopes2_lateral, theta, true);
    s3 = interpolate( angulo, slopes3_lateral, theta, true);
  }
  else  {
    std::cout << "Error: Invalid optical detector orientation." << std::endl;
    exit(1);
  }

  // add border correction
  pars_ini[0] = pars_ini[0] + s1 * r_distance;
  pars_ini[1] = pars_ini[1] + s2 * r_distance;
  pars_ini[2] = pars_ini[2] + s3 * r_distance;
  pars_ini[3] = pars_ini[3];



  // calculate correction factor
  double GH_correction = GaisserHillas(distance, pars_ini);
  // apply correction
  double hits_rec = GH_correction*hits_geo/cosine;

  // round to integer value, cannot have non-integer number of hits
  int hits_vuv = std::round(hits_rec);

  return hits_vuv;
}



// visible number of hits calculation function
// does not apply any correction
// Visible hits calculation
int VisHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type) {

  // 1). calculate total number of hits of VUV photons on reflective foils via solid angle + Gaisser-Hillas corrections:

  // set cathode plane struct for solid angle function
  acc cathode_plane;
  cathode_plane.ax = x_foils; cathode_plane.ay = y_foils; cathode_plane.az = z_foils;       // centre coordinates of cathode plane
  cathode_plane.w = y_dimension_foils; cathode_plane.h = z_dimension_foils;                     // width and height in cm

  // get scintpoint coords relative to centre of cathode plane
  TVector3 cathodeCentrePoint(x_foils,y_foils,z_foils);
  TVector3 ScintPoint_relative = ScintPoint - cathodeCentrePoint;

  // calculate solid angle of cathode from the scintillation point
  double solid_angle_cathode = solid(cathode_plane, ScintPoint_relative);

  // calculate distance and angle between ScintPoint and hotspot
  // vast majority of hits in hotspot region directly infront of scintpoint,therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
  double distance_cathode = std::abs(x_foils - ScintPoint[0]);

  double cosine_cathode = 1;
  double theta_cathode = 0;

  // calculate hits on cathode plane via geometric acceptance
  double cathode_hits_geo = exp(-1.*distance_cathode/L_abs) * (solid_angle_cathode / (4.*pi)) * Nphotons_created;

  // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
  // offset angle bin
  int j = (theta_cathode/delta_angulo);
  // correction
  double pars_ini[4] = {fGHVUVPars[0][j], fGHVUVPars[1][j], fGHVUVPars[2][j], fGHVUVPars[3][j]};

  // gh border
  double r_distance = sqrt( pow(ScintPoint[1] - y_foils, 2) + pow(ScintPoint[2] - z_foils, 2));

  double s1 = interpolate( angulo, slopes1, theta_cathode, true);
  double s2 = interpolate( angulo, slopes2, theta_cathode, true);
  double s3 = interpolate( angulo, slopes3, theta_cathode, true);

  pars_ini[0] = pars_ini[0] + s1 * r_distance;
  pars_ini[1] = pars_ini[1] + s2 * r_distance;
  pars_ini[2] = pars_ini[2] + s3 * r_distance;
  pars_ini[3] = pars_ini[3];


  double GH_correction = GaisserHillas(distance_cathode, pars_ini);

  double cathode_hits_rec = GH_correction*cathode_hits_geo/cosine_cathode;


 // 2). calculate number of these hits which reach the optical channel from the hotspot via solid angle

  // calculate hotspot location
  TVector3 v_to_wall(x_foils - ScintPoint[0],0,0);
  TVector3 hotspot = ScintPoint + v_to_wall;

  // distance to hotspot
  double distance_vuv = sqrt(pow(ScintPoint[0] - hotspot[0],2) + pow(ScintPoint[1] - hotspot[1],2) + pow(ScintPoint[2] - hotspot[2],2));
  // distance from hotspot to arapuca
  double distance_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2) + pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
  // angle between hotspot and arapuca
  double cosine_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2)) / distance_vis;
  double theta_vis = acos(cosine_vis)*180./pi;

  // solid angle :
  double solid_angle_detector = 0;
  // rectangular aperture
  if (optical_detector_type == 1) {
    // set Arapuca geometry struct for solid angle function
    acc detPoint;
    detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
    detPoint.w = y_dimension_detector; detPoint.h = z_dimension_detector; // width and height in cm of arapuca active window

    // get hotspot coordinates relative to detpoint
    TVector3 emission_relative = hotspot - OpDetPoint;

    // calculate solid angle of optical channel
    solid_angle_detector = solid(detPoint, emission_relative);
  }
  // disk aperture
  else if (optical_detector_type == 0) {
    /*
    // offset in z-y plane
    double d = sqrt(pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
    // drift distance (in x)
    double h =  sqrt(pow(hotspot[0] - OpDetPoint[0],2));
    // Solid angle of a disk
    solid_angle_detector = Disk_SolidAngle(d, h, radius);
    */
    solid_angle_detector = Omega_Dome_Model(distance_vis, theta_vis);
  }
  else {
    std::cout << "Error: Invalid optical detector type." << std::endl;
    exit(1);
  }

  // calculate number of hits via geometeric acceptance
  double hits_geo = (solid_angle_detector / (2*pi)) * cathode_hits_rec;


  // apply correction curves, interpolation
  int k = (theta_vis/delta_angulo);
  //double x_position_corr = ScintPoint[0];
  //double corr = interpolate(fVIS_XPos, fVIS_Pars[k], x_position_corr, false);

  //double hits_rec = corr*hits_geo/cosine_vis
  double hits_rec = hits_geo;


  // apply border correction
  // calculate radial distance from centre
  //double r_distance = sqrt( pow(ScintPoint[1] - y_foils, 2) + pow(ScintPoint[2] - z_foils, 2));
  // interpolate in x for each r bin
  std::vector<double> interp_vals = {0,0,0,0,0};
  for (int i = 0; i < VIS_RS100_SBND_Borders[k].size(); i++){
    interp_vals[i] = interpolate(vDistances_x, VIS_RS100_SBND_Borders[k][i], std::abs(ScintPoint[0]), false);
  }
  // interpolate in r
  double border_correction = interpolate(vDistances_r, interp_vals, r_distance, false);
  // apply correction
  double hits_rec_borders = border_correction * hits_rec / cosine_vis;

  // round final result
  int hits_vis = std::round(hits_rec_borders);


//  int hits_vis = hits_rec;

  return hits_vis;
}

// semi-analytic model utility functions
// gaisser-hillas function definition
Double_t GaisserHillas(double x,double *par) {
  //This is the Gaisser-Hillas function
  Double_t X_mu_0=par[3];
  Double_t Normalization=par[0];
  Double_t Diff=par[1]-X_mu_0;
  Double_t Term=pow((x-X_mu_0)/Diff,Diff/par[2]);
  Double_t Exponential=TMath::Exp((par[1]-x)/par[2]);

  return ( Normalization*Term*Exponential);
}

// solid angle of rectanglular aperture calculation functions
double omega(const double &a, const double &b, const double &d){

  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));

}

double solid(const acc& out, const TVector3 &v){

  //v is the position of the track segment with respect to
  //the center position of the arapuca window

  // arapuca plane fixed in x direction

  if( v.Y()==0.0 && v.Z()==0.0){
    return omega(out.w,out.h,v.X());
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(b-B),d)+omega(2*A,2*(b-B),d)+omega(2*(a-A),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(b-B),d)-omega(2*A,2*(b-B),d)+omega(2*(A+a),2*B,d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(B+b),d)-omega(2*(a-A),2*B,d)+omega(2*A,2*(B+b),d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}

double solid_lateral(const acc& out, const TVector3 &v) {

  //v is the position of the track segment with respect to
  //the center position of the arapuca window

  // arapuca plane fixed in y direction

  if( v.X()==0.0 && v.Z()==0.0){
    return omega(out.w,out.h,v.Y());
  }

  if( (std::abs(v.X()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.X())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.Y());
    double to_return = (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.X()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.X())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.Y());
    double to_return = (omega(2*(a-A),2*(b-B),d)+omega(2*A,2*(b-B),d)+omega(2*(a-A),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.X()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.X())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.Y());
    double to_return = (omega(2*(A+a),2*(b-B),d)-omega(2*A,2*(b-B),d)+omega(2*(A+a),2*B,d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.X()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.X())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.Y());
    double to_return = (omega(2*(a-A),2*(B+b),d)-omega(2*(a-A),2*B,d)+omega(2*A,2*(B+b),d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}
/*
// solid angle of circular aperture
double Disk_SolidAngle(double *x, double *p) {
  const double d = x[0];
  const double h = x[1];
  const double b = p[0];
  if(b <= 0. || d < 0. || h <= 0.) return 0.;
  const double aa = TMath::Sqrt(h*h/(h*h+(b+d)*(b+d)));
  if(d == 0) {
    return 2.*TMath::Pi()*(1.-aa);
  }
  const double bb = TMath::Sqrt(4*b*d/(h*h+(b+d)*(b+d)));
  const double cc = 4*b*d/((b+d)*(b+d));

  if(TMath::Abs(ROOT::Math::comp_ellint_1(bb) - bb) < 1e-10 && TMath::Abs(ROOT::Math::comp_ellint_3(cc,bb) - cc) <1e-10) {
    throw(std::runtime_error("please do gSystem->Load(\"libMathMore.so\") before running Disk_SolidAngle for the first time!"));
  }
  if(d < b) {
    return 2.*TMath::Pi() - 2.*aa*(ROOT::Math::comp_ellint_1(bb) + TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb));
  }
  if(d == b) {
    return TMath::Pi() - 2.*aa*ROOT::Math::comp_ellint_1(bb);
  }
  if(d > b) {
    return 2.*aa*(TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb) - ROOT::Math::comp_ellint_1(bb));
  }

  return 0.;
}

double Disk_SolidAngle(double d, double h, double b) {
  double x[2] = { d, h };
  double p[1] = { b };
  return Disk_SolidAngle(x,p);
}
*/

double interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate ) {
  int size = xData.size();
  int i = 0;                                          // find left end of interval for interpolation
  if ( x >= xData[size - 2] )                         // special case: beyond right end
    {
      i = size - 2;
    }
  else
    {
      while ( x > xData[i+1] ) i++;
    }
  double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1]; // points on either side (unless beyond ends)
  if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
    {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
    }
  double dydx = ( yR - yL ) / ( xR - xL );            // gradient
  return yL + dydx * ( x - xL );                      // linear interpolation
}
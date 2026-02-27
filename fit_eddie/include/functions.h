Double_t pol1(double *x,double *par)
{
  return (*x*par[1] + par[0]);
}
Double_t pol2(double *x,double *par)
{
  return (pow(*x,2)*par[2] + *x*par[1] + par[0]);
}
Double_t pol3(double *x,double *par)
{
  return (pow(*x,3)*par[3] + pow(*x,2)*par[2] + *x*par[1] + par[0]);
}
// Above are three definitions of polynomical funtions.

TGraphErrors *Profile_to_Graph(TProfile *px, TProfile *py, int j) {

  TAxis *xaxis = py->GetXaxis();
  Int_t nbins = xaxis->GetNbins();
  std::vector<double> vx, vy, vex, vey;
  for (Int_t bin=0; bin <= nbins; bin++) {
    double x_value = px->GetBinContent(bin);
    double y_value = py->GetBinContent(bin);
    double x_error = px->GetBinError(bin);
    double y_error = py->GetBinError(bin);

    //cout << x_value << " " << y_value << endl;
    if(y_value == 0) continue;
    vx.push_back(x_value);
    vy.push_back(y_value);
    vex.push_back(x_error);
    vey.push_back(y_error);
  }

  TGraphErrors *gr = new TGraphErrors(vx.size(), &vx[0], &vy[0], 0, &vey[0]);
  return gr;
}


// TGraphErrors *Profile_to_Graph(TProfile *px, TProfile *py, int /*j*/) {
//   TAxis *ax = py->GetXaxis();
//   const int nbins = ax->GetNbins();
//   std::vector<double> vx, vy, vex, vey;
//   vx.reserve(nbins); vy.reserve(nbins); vex.reserve(nbins); vey.reserve(nbins);

//   for (int bin = 1; bin <= nbins; ++bin) {              // 1..nbins (skip underflow/overflow)
//     const double entries = py->GetBinEntries(bin);
//     if (entries <= 0) continue;                         // only drop truly empty bins

//     const double x  = ax->GetBinCenter(bin);            // robust x
//     const double ex = 0.5 * ax->GetBinWidth(bin);       // or 0 if you prefer
//     const double y  = py->GetBinContent(bin);           // mean ratio (can be 0)
//     const double ey = py->GetBinError(bin);             // error on the mean

//     vx.push_back(x); vy.push_back(y);
//     vex.push_back(ex); vey.push_back(ey);
//   }

//   auto *gr = new TGraphErrors((int)vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
//   return gr;
// }



///////////////////////////////////////////////////
Double_t GaisserHillas(double *x,double *par)
{
  //This is the Gaisser-Hillas function
  Double_t X_mu_0=par[3];
  Double_t Normalization=par[0];
  Double_t Diff=par[1]-X_mu_0;
  Double_t Term=pow((*x-X_mu_0)/Diff,Diff/par[2]);
  Double_t Exponential=TMath::Exp((par[1]-*x)/par[2]);

  return ( Normalization*Term*Exponential);
}

//Distance to the center in the Y-Z Plane
double GetDistanceCenter(const double center[2], double z, double y){
  //z = abs(z - center[1]) - center[1];
  //y = abs(y) - center[0];
  z -= center[1];
  y -= center[0];

  return  sqrt(y*y + z*z);
}

// solid angle of rectanglular aperture
// structure definition for solid angle of rectangle function
struct acc{
  // ax,ay,az = centre of rectangle; w = width; h = height
  double ax, ay, az, w, h;
};

double Rectangle_SolidAngle(double a, double b, double d){

  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));

}

double Rectangle_SolidAngle(acc& out, TVector3 v){

  //v is the position of the track segment with respect to
  //the center position of the arapuca window

  // arapuca plane fixed in x direction

  if( v.Y()==0.0 && v.Z()==0.0){
    return Rectangle_SolidAngle(out.w,out.h,v.X());
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.X());
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(B+b),d)-Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*(A+a),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.X());
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(b-B),d)+Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.X());
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(b-B),d)-Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(A+a),2*B,d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.X());
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(B+b),d)-Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}

double Rectangle_SolidAngle_lateral(acc& out, TVector3 v){

  //v is the position of the track segment with respect to
  //the center position of the arapuca window

  // arapuca plane fixed in y direction

  if( v.X()==0.0 && v.Z()==0.0){
    return Rectangle_SolidAngle(out.w,out.h,v.Y());
  }

  if( (std::abs(v.X()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.X())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.Y());
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(B+b),d)-Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*(A+a),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.X()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.X())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.Y());
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(b-B),d)+Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.X()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.X())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.Y());
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(b-B),d)-Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(A+a),2*B,d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.X()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.X())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.Y());
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(B+b),d)-Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}

bool _mathmore_loaded_ = false;

double Omega(double* x, double *p) {
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
  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  if(TMath::Abs(ROOT::Math::comp_ellint_1(bb) - bb) < 1e-10 && TMath::Abs(ROOT::Math::comp_ellint_3(cc,bb) - cc) <1e-10) {
    throw(std::runtime_error("please do gSystem->Load(\"libMathMore.so\") before running Omega for the first time!"));
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

double Omega(double d, double h, double b) {
  double x[2] = { d, h };
  double p[1] = { b };
  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  return Omega(x,p);
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

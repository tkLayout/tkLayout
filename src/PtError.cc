#include <cmath>
#include <iostream>

#include <PtError.hh>
#include <global_constants.hh>
#include "SimParms.hh"
#include "Units.hh"

using namespace std;

// Default values for gloabl parameters
double ptError::IP_length = 70;     // mm

void ptError::defaultParameters() {
  Module_pitch = defaultModulePitch;
  Module_strip_l = defaultStripLength;
  Module_z = defaultModuleZ;
  Module_r = defaultModuleR;
  Module_d = defaultModuleD;
  Module_h = defaultModuleHeight;
  zCorrelation = defaultZCorrelation;
  moduleType = ModuleSubdetector::BARREL; // Barrel is the default moduleType 
  endcapType = ModuleShape::RECTANGULAR; // Barrel is the default moduleType 
}

double ptError::computeErrorBE(double p) {
  double A = 0.3 * SimParms::getInstance().magField() * Module_r / 2.;
  double a = pow(p/A,2);
  double g = 1/a;
  //std::cout << "g = " << g << std::endl;
  //std::cout << "a = " << a << std::endl;

  double one_over_sin_cos = (Module_z/Module_r) + (Module_r/Module_z);

  // Barrel and Endcap coefficients
  double f1_sq = pow((1-g),2)*(a-1);
  double f2_sq = (a-1); // For material interaction

  // Endcap only coefficients
  double f3 = 2-g; // TODO: this is 1 for radial strips
  double f4 = 1-g;
  double f5 = (1-g)*one_over_sin_cos; // For material interaction
  //double f5 = (1-g)/sin(Module_theta)/cos(Module_theta); // For material interaction
  //std::cout << "f1 = " << f1 <<std::endl;
  //std::cout << "f2 = " << f2 <<std::endl;
  //std::cout << "f3 = " << f3 <<std::endl;
  //std::cout << "f4 = " << f4 <<std::endl;
  //std::cout << "f5 = " << f5 <<std::endl;

  double relativeError2;
  double deltaPhi;
  double deltaTheta;

  if (material.radiation==0) {
    deltaPhi = 0;
    deltaTheta = 0;
  } else {
    double scatteringAngle;
    // Barrel modules:
    // p represents the transverse momentum p_T in GeV
    // in this case theta represents the average scattering angle in the (r,z) plane
    // (That is: i should divide p by sin(theta) [bigger momentum, less scattering]
    // to obtain the scattering angle in the plane containing the z axis and the track
    // and then I should use the distance between two hits [larger lever arm]
    // to obtain the error on phi. But the two effects cancel, so I can directly use pT
    // and deltaR of the hits
    scatteringAngle = (13.6*Units::MeV) / (p/Units::MeV) * sqrt(material.radiation) * (1 + 0.038 * log(material.radiation));
    deltaPhi = scatteringAngle;

    // Endcap modules:
    // Here I need the actual \Delta theta value (TODO: verify)
    // So I must go from pT -> p, which means I have to divide p by sin(theta)
    double sinTh = Module_r / sqrt( pow(Module_r,2) + pow(Module_z,2) );
    deltaTheta = scatteringAngle * sinTh;
  }

  switch ( moduleType ) {
  case BARREL:
    relativeError2 = f1_sq * pow(Module_pitch / Module_d, 2) / 6;
    if (material.radiation!=0) relativeError2 += f2_sq * pow( deltaPhi ,2);
    //#ifdef ADD_ERR
    //relativeError = f1 * sqrt(2*pow(Module_pitch,2)+pow(Module_strip_l*0.01,2)) /sqrt(12) / Module_d;
    //#endif
    break;
  case ENDCAP:
      {
//        static std::ofstream ofs("pterr.txt");
//        ofs << Module_z << "\t" << Module_r << "\t" << Module_d << "\t" << Module_ed << "\t" << Module_d*Module_r/Module_z << "...\t"; 
    relativeError2 = f1_sq * pow(Module_pitch / (Module_d * Module_r / Module_z),2) / 6; // * 1/tg(theta)
//        ofs << relativeError2 << "\t";
    relativeError2 += f2_sq * pow( deltaPhi ,2);
//        ofs << relativeError2 << "\t";
    relativeError2 += pow( f3 * Module_strip_l / Module_r ,2) / 12;
//        ofs << relativeError2 << "\t";
    relativeError2 += pow( f4 * IP_length / Module_z ,2);
//        ofs << relativeError2 << "\t";
    if (material.radiation!=0) relativeError2 += pow( f5 * deltaTheta, 2);
//        ofs << relativeError2 << std::endl;
      }
    break;
  default:
    // This should never happen
    relativeError2 = 1;
    std::cerr << "Error: in TODO:addFunctionNameHere"
	      << "I just found a module which is neither Barrel nor Endcap" << std::endl;
  }

  return sqrt(relativeError2);
}

double ptError::computeError(double p) { // CUIDADO TODO this should be refactored (remove type checks!), but we have to refit the new occupancy data first
  if (moduleType == BARREL && fabs(Module_tilt) > 1e-3) {
    double errB = computeErrorBE(p);
    moduleType = ENDCAP;
    double errE = computeErrorBE(p);
    moduleType = BARREL;
    return errB*pow(cos(Module_tilt),2) + errE*pow(sin(Module_tilt),2);
  } else {
    return computeErrorBE(p);
  }
}

double ptError::oneSidedGaussInt(double x) {
 return (0.5*erf(double(x)/sqrt(2)));
}

double ptError::probabilityInside_norm(double cut, double value) {
  double p_under_lower, p_under_higher;
  p_under_lower  = oneSidedGaussInt(-cut-value) /* + 0.5 */;
  p_under_higher = oneSidedGaussInt(cut-value) /* + 0.5 */;

  return (p_under_higher - p_under_lower);
}

double ptError::probabilityInside(double cut, double value, double value_err) {
  double mycut = fabs(cut/value_err);
  double myvalue = fabs(value/value_err);
  return probabilityInside_norm(mycut, myvalue);
}

/* */



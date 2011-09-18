#include <cmath>
#include <iostream>

#include <ptError.h>
#include <module.hh>
#include <global_constants.h>

using namespace std;

// Default values for gloabl parameters
double ptError::IP_length = 70;     // mm
double ptError::B = insur::magnetic_field; // T

double ptError::minimumPt = 0.3;
double ptError::maximumPt = 30;

void ptError::defaultParameters() {
  Module_pitch = defaultModulePitch;
  Module_strip_l = defaultStripLength;
  Module_z = defaultModuleZ;
  Module_r = defaultModuleR;
  Module_d = defaultModuleD;
  moduleType = Module::Barrel; // Barrel is the default moduleType 
  endcapType = Module::Rectangular; // Barrel is the default moduleType 
}

double ptError::geometricEfficiency() {
  double inefficiency;
  switch ( moduleType ) {
  case Module::Barrel:
    inefficiency = Module_d * Module_z / Module_strip_l / Module_r ;
    break;
  case Module::Endcap:
    inefficiency = Module_d * Module_r / Module_strip_l / Module_z ;
    break;
  default:
    // This should never happen
    inefficiency = 0;
    std::cerr << "Error: in TODO:addFunctionNameHere"
              << "I just found a module which is neither Barrel nor Endcap" << std::endl;
  }
  return (1-inefficiency);
}

double ptError::stripsToP(double strips) {
  double A = 0.3 * B * Module_r / 1000. / 2.; // GeV
  double p;
  double x;
  double effective_d;

  switch ( moduleType ) {
  case Module::Barrel:
    effective_d = Module_d;
    break;
  case Module::Endcap:
    effective_d = Module_d * Module_r / Module_z;
    break;
  default:
    // This should never happen
    strips = 0;
    effective_d = 0;
    std::cerr << "Error: in TODO:addFunctionNameHere"
              << "I just found a module which is neither Barrel nor Endcap" << std::endl;
  }

  x = strips * Module_pitch;
  p = A * sqrt( pow(effective_d/x,2) + 1 );
  return p;
}

double ptError::pToStrips(double p) {
  double A = 0.3 * B * Module_r / 1000. / 2.; // GeV
  double a = pow(p/A,2);
  double strips;

  switch ( moduleType ) {
  case Module::Barrel:
    strips = Module_d / sqrt(a-1);
    break;
  case Module::Endcap:
    strips = Module_d / sqrt(a-1) * (Module_r / Module_z);
    break;
  default:
    // This should never happen
    strips = 0;
    std::cerr << "Error: in TODO:addFunctionNameHere"
              << "I just found a module which is neither Barrel nor Endcap" << std::endl;
  }

  strips /= Module_pitch;
  return strips;
} 


double ptError::computeError(double p) {
  double A = 0.3 * B * Module_r / 1000. / 2.; // GeV
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
    scatteringAngle = (13.6) / (1000 * p) * sqrt(material.radiation) * (1 + 0.038 * log(material.radiation));
    deltaPhi = scatteringAngle;

    // Endcap modules:
    // Here I need the actual \Delta theta value (TODO: verify)
    // So I must go from pT -> p, which means I have to divide p by sin(theta)
    double sinTh = Module_r / sqrt( pow(Module_r,2) + pow(Module_z,2) );
    deltaTheta = scatteringAngle * sinTh;
  }

  switch ( moduleType ) {
  case Module::Barrel:
    relativeError2 = f1_sq * pow(Module_pitch / Module_d, 2) / 6;
    if (material.radiation!=0) relativeError2 += f2_sq * pow( deltaPhi ,2);
    //#ifdef ADD_ERR
    //relativeError = f1 * sqrt(2*pow(Module_pitch,2)+pow(Module_strip_l*0.01,2)) /sqrt(12) / Module_d;
    //#endif
    break;
  case Module::Endcap:
    relativeError2 = f1_sq * pow(Module_pitch / (Module_d * Module_r / Module_z),2) / 6; // * 1/tg(theta)
    relativeError2 += f2_sq * pow( deltaPhi ,2);
    relativeError2 += pow( f3 * Module_strip_l / Module_r ,2) / 12;
    relativeError2 += pow( f4 * IP_length / Module_z ,2);
    if (material.radiation!=0) relativeError2 += pow( f5 * deltaTheta, 2);
    break;
  default:
    // This should never happen
    relativeError2 = 1;
    std::cerr << "Error: in TODO:addFunctionNameHere"
	      << "I just found a module which is neither Barrel nor Endcap" << std::endl;
  }

  return sqrt(relativeError2);
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

double ptError::find_probability(double target, double ptCut) {

  double lowerPt = minimumPt;
  double higherPt = maximumPt;

  double lowerProbability =  geometricEfficiency() * probabilityInside(1/ptCut, 1/lowerPt, computeError(lowerPt)/lowerPt);
  double higerProbability =  geometricEfficiency() * probabilityInside(1/ptCut, 1/higherPt, computeError(higherPt)/higherPt);

  double testPt;
  double testProbability;

  //std::cerr << "LO: prob("<<lowerPt<<") = " << lowerProbability << std::endl;
  //std::cerr << "HI: prob("<<higherPt<<") = " << higerProbability << std::endl;

  if ((target<lowerProbability)||(target>higerProbability)) return -1; // probability not reachable within desired pt range

  for (int i=0; i<100; ++i) {
    //std::cerr << std::endl;
    //std::cerr << std::endl;
    //std::cerr << std::endl;
    //std::cerr << "****** STEP # " << i << "*********" << std::endl;
    testPt = sqrt(lowerPt*higherPt);
    testProbability = geometricEfficiency() * probabilityInside(1/ptCut, 1/testPt, computeError(testPt)/testPt);
    //std::cerr << "Geometric efficiency is " << geometricEfficiency() << std::endl;
    //std::cerr << "LO: prob("<<lowerPt<<") = " << lowerProbability << std::endl;
    //std::cerr << "HI: prob("<<higherPt<<") = " << higerProbability << std::endl;
    //std::cerr << "XX: prob("<<testPt<<") = " << testProbability << std::endl;

    if (testProbability>target) {
      higherPt = testPt;
      // higerProbability = testProbability; // debug
    } else {
      lowerPt = testPt;
      // lowerProbability = testProbability; // debug
    }

    if (fabs(testProbability-target)<1E-5) return (testPt);
  }

  return -1; // could not converge
}



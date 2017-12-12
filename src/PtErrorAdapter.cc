#include "PtErrorAdapter.hh"
#include "SimParms.hh"
#include "Units.hh"

const double PtErrorAdapter::minimumPt = 0.3*Units::GeV;
const double PtErrorAdapter::maximumPt = 30*Units::GeV;

const double PtErrorAdapter::ptMinFit = 0.22*Units::GeV;
const double PtErrorAdapter::ptMaxFit = 10.*Units::GeV;
// log(pt/z) distribution parameters for 12000 events at 14 TeV
const double PtErrorAdapter::ptFitParamsLow[]  = { 2.52523e+01, -6.84183e+00, -1.20149e+01,  1.89314e+00}; // 0.22 GeV to 1 GeV  Chi^2 / dof = 22.0408/16 = 1.37755
const double PtErrorAdapter::ptFitParamsMid[]  = { 7.27638e-01, -1.04041e+00,  8.56495e+00,  6.52714e-03}; // 1 GeV to 4 GeV     Chi^2 / dof = 84.2299/71 = 1.18634
const double PtErrorAdapter::ptFitParamsHigh[] = { 4.66514e+01, -2.88910e+00, -3.78716e+01,  1.26635e-01}; // 4 GeV to 10 GeV    Chi^2 / dof = 102.736/135 = 0.76101

void PtErrorAdapter::setPterrorParameters() {
  myPtError.setDistance( mod_.dsDistance() );
  myPtError.setEffectiveDistance( mod_.effectiveDsDistance() );
  myPtError.setPitch(mod_.outerSensor().pitch());
  myPtError.setStripLength( mod_.outerSensor().stripLength() );
  XYZVector center = mod_.center();
  myPtError.setZ(fabs(center.Z()));
  myPtError.setR(center.Rho());
  myPtError.setHeight(mod_.length());
  myPtError.setZCorrelation(mod_.zCorrelation());
  myPtError.setModuleType(mod_.subdet());
  myPtError.setTilt(mod_.tiltAngle());
}


double PtErrorAdapter::getTriggerProbability(const double& trackPt, const double& stereoDistance /*= 0*/, const int& triggerWindow /* = 0 */ ) {
  setPterrorParameters();
  if (stereoDistance!=0) {
    myPtError.setDistance(stereoDistance);
    if (fabs(mod_.tiltAngle()) < 1e-3) myPtError.setEffectiveDistance(mod_.dsDistance()); // CUIDADO temporary fix!!! this belongs inside a function
    else myPtError.setEffectiveDistance(mod_.dsDistance()*sin(mod_.center().Theta())/sin(mod_.center().Theta()+mod_.tiltAngle()));

  }
  int thisTriggerWindow;
  if (triggerWindow!=0) thisTriggerWindow = triggerWindow;
  else thisTriggerWindow = mod_.triggerWindow();
  double pt_cut = stripsToP(thisTriggerWindow/2.);
  // Error on curvatre is the relative error of trackPt times the
  // curvature (cur = 1/pt)
  double cur_error = myPtError.computeError(trackPt) / trackPt; 
  double result;
  result = myPtError.probabilityInside(1/pt_cut, 1/trackPt, cur_error) * mod_.geometricEfficiency();
  // std::cerr << "trigger prob @ " << trackPt << " GeV/c is " << result <<std::endl; // debug
  return result;
}

// TODO: this is VERY ugly!!! :( sorry: hurry !

double PtErrorAdapter::getTriggerFrequencyTruePerEventAbove(const double& myCut) {
  return getTriggerFrequencyTruePerEventBetween(myCut, ptMaxFit);
}

double PtErrorAdapter::getParticleFrequencyPerEventAbove(const double& myCut) {
  return getParticleFrequencyPerEventBetween(myCut, ptMaxFit);
}

double PtErrorAdapter::getTriggerFrequencyTruePerEventBelow(const double& myCut) {
  return getTriggerFrequencyTruePerEventBetween(ptMinFit, myCut);
}

// TODO: this and the next one should be merged in the same code...
double PtErrorAdapter::getTriggerFrequencyTruePerEventBetween(double myLowCut, double myHighCut) { 
  if (myLowCut<ptMinFit) myLowCut=ptMinFit;
  if (myHighCut>ptMaxFit) myHighCut=ptMaxFit;
  double r        = mod_.center().Rho();
  double ptMin    = MAX(0.3 * SimParms::getInstance().magField() * r, myLowCut);
  double integral = 0.0;
  double dPt      = 0.05*Units::GeV;
  double etaphi   = mod_.phiAperture()/(2.0*M_PI)*fabs(mod_.etaAperture())/6.0;
  double nPt;
  const double* ptFitParams;
  for (double pt = ptMin; pt < myHighCut; pt += dPt) {
    if      (pt < 1*Units::GeV) ptFitParams = ptFitParamsLow;
    else if (pt < 6*Units::GeV) ptFitParams = ptFitParamsMid;
    else                        ptFitParams = ptFitParamsHigh;

    // TODO: // nPt was fitted originally in different units, GeV=1, not GeV=1000 -> hence pt/1000 & final multiplying by 1000
    nPt = exp(ptFitParams[0]
              + ptFitParams[1] * pt/1000.
              + ptFitParams[2] * pow(pt/1000.,-0.1)
              + ptFitParams[3] * pow(pt/1000.,2))/12000 * 1000.;

    integral += getTriggerProbability(pt) * (nPt/4e-2) * etaphi * dPt;
  }

  return integral;
}

// TODO: this and the previous one should be merged in the same code...
double PtErrorAdapter::getParticleFrequencyPerEventBetween(double myLowCut, double myHighCut) {
  if (myLowCut<ptMinFit) myLowCut=ptMinFit;
  if (myHighCut>ptMaxFit) myHighCut=ptMaxFit;
  double r        = mod_.center().Rho();
  double ptMin    = MAX(0.3 * SimParms::getInstance().magField() * r, myLowCut);
  double integral = 0.0;
  double dPt      = 0.05*Units::GeV;
  double etaphi = mod_.phiAperture()/(2.0*M_PI)*fabs(mod_.etaAperture())/6.0;
  double nPt;
  const double* ptFitParams;
  for (double pt = ptMin; pt < myHighCut; pt += dPt) {
    if      (pt < 1*Units::GeV) ptFitParams = ptFitParamsLow;
    else if (pt < 6*Units::GeV) ptFitParams = ptFitParamsMid;
    else                        ptFitParams = ptFitParamsHigh;

    // TODO: nPt was fitted originally in different units, GeV=1, not GeV=1000 -> hence pt/1000 & final multiplying by 1000
    nPt = exp(ptFitParams[0]
              + ptFitParams[1] * pt/1000.
              + ptFitParams[2] * pow(pt/1000.,-0.1)
              + ptFitParams[3] * pow(pt/1000.,2))/12000 * 1000.;

    integral += (nPt/4e-2) * etaphi * dPt;
  }

  return integral;
}

double PtErrorAdapter::getTriggerFrequencyFakePerEvent() {
  return pow(mod_.hitOccupancyPerEvent(),2) * mod_.triggerWindow() * mod_.sensors().back().numChannels();
}

double PtErrorAdapter::getPtThreshold(const double& myEfficiency) {
  double pt_cut = stripsToP(mod_.triggerWindow()/2.);
  return find_probability(myEfficiency, pt_cut);
}

double PtErrorAdapter::stripsToP(double strips) const {
  double A = 0.3 * SimParms::getInstance().magField() * mod_.center().Rho() / 2.;
  double p;
  double x;
  double effective_d = mod_.effectiveDsDistance();

  x = strips * mod_.outerSensor().pitch();
  p = A * sqrt( pow(effective_d/x,2) + 1 );
  return p;
}

double PtErrorAdapter::pToStrips(double p) const {
  double A = 0.3 * SimParms::getInstance().magField() * mod_.center().Rho() / 2.;
  double a = pow(p/A,2);
  double strips = mod_.effectiveDsDistance() / sqrt(a-1) / mod_.outerSensor().pitch();

  return strips;
} 

double PtErrorAdapter::find_probability(double target, double ptCut) {

  double lowerPt = minimumPt;
  double higherPt = maximumPt;

  double lowerProbability =  mod_.geometricEfficiency() * myPtError.probabilityInside(1/ptCut, 1/lowerPt, myPtError.computeError(lowerPt)/lowerPt);
  double higerProbability =  mod_.geometricEfficiency() * myPtError.probabilityInside(1/ptCut, 1/higherPt, myPtError.computeError(higherPt)/higherPt);

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
    testProbability = mod_.geometricEfficiency() * myPtError.probabilityInside(1/ptCut, 1/testPt, myPtError.computeError(testPt)/testPt);
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

double PtErrorAdapter::getPtCut() const {
  return stripsToP(mod_.triggerWindow()/2.);
}

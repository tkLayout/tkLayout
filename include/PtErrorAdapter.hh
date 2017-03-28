#ifndef PT_ERROR_ADAPTER_H
#define PT_ERROR_ADAPTER_H

#include "global_constants.hh"
#include "PtError.hh"
#include "Module.hh"


class PtErrorAdapter {
  static const double minimumPt;
  static const double maximumPt;
  static const double ptMinFit;
  static const double ptMaxFit;
  // log(pt/z) distribution parameters for 12000 events at 14 TeV
  static const double ptFitParamsLow[]; // 0.22 GeV to 1 GeV  Chi^2 / dof = 22.0408/16 = 1.37755
  static const double ptFitParamsMid[]; // 1 GeV to 4 GeV     Chi^2 / dof = 84.2299/71 = 1.18634
  static const double ptFitParamsHigh[]; // 4 GeV to 10 GeV    Chi^2 / dof = 102.736/135 = 0.76101

  ptError myPtError;
  const DetectorModule& mod_;

   void setPterrorParameters();
public:
   PtErrorAdapter(const DetectorModule& m) : mod_(m) { setPterrorParameters(); }
   double getTriggerProbability(const double& trackPt, const double& stereoDistance = 0, const int& triggerWindow = 0);
   double getTriggerFrequencyTruePerEventAbove(const double& myCut);
   double getParticleFrequencyPerEventAbove(const double& myCut);
   double getTriggerFrequencyTruePerEventBelow(const double& myCut);
   double getTriggerFrequencyTruePerEventBetween(double myLowCut, double myHighCut);
   double getParticleFrequencyPerEventBetween(double myLowCut, double myHighCut);
   double getTriggerFrequencyFakePerEvent();
   double getPtThreshold(const double& myEfficiency);
   double stripsToP(double strips) const;
   double pToStrips(double strips) const;
   double find_probability(double target, double originalcut);
   double getPtCut() const;

};

#endif

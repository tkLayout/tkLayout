#ifndef TRIGGERDISTANCETUNINGPLOTS_H
#define TRIGGERDISTANCETUNINGPLOTS_H

#include <string>
#include <map>
#include <vector>
#include <utility>

#include "TH1.h"
#include "TGraphErrors.h"

#include "Tracker.hh"
#include "SimParms.hh"
#include "PtErrorAdapter.hh"
#include "Visitor.hh"
#include "Bag.hh"
#include "SummaryTable.hh"


typedef std::map<const DetectorModule*, std::map<int, double> > ModuleOptimalSpacings;


class TriggerDistanceTuningPlotsVisitor : public ConstGeometryVisitor {
  typedef std::vector<const DetectorModule*> ModuleVector;
  std::map<std::string, ModuleVector> selectedModules_;
  std::set<double> foundSpacing_;
  const std::vector<double>& triggerMomenta_;
  const unsigned int nWindows_ = 5;

  profileBag& myProfileBag_;
  std::map<std::string, bool> preparedProfiles_;
  std::map<std::string, bool> preparedTurnOn_;

  


  double findXThreshold(const TProfile& aProfile, const double& yThreshold, const bool& goForward) {  
    // TODO: add a linear interpolation here
    if (goForward) {
      double xThreshold=0;
      double binContent;
      for (int i=1; i<=aProfile.GetNbinsX(); ++i) {
        binContent = aProfile.GetBinContent(i);
        if ((binContent!=0)&&(binContent<yThreshold)) {
          // TODO: add a linear interpolation here
          xThreshold = aProfile.GetBinCenter(i);
          return xThreshold;
        }
      }
      return 100;
    } else {
      double xThreshold=100;
      double binContent;
      for (int i=aProfile.GetNbinsX(); i>=1; --i) {
        binContent = aProfile.GetBinContent(i);
        if ((binContent!=0)&&(binContent>yThreshold)) {
          // TODO: add a linear interpolation here
          xThreshold = aProfile.GetBinCenter(i);
          return xThreshold;
        }
      }
      return 0;
    }
  }
 

public:
  TH1D optimalSpacingDistribution, optimalSpacingDistributionAW;
  TH1D spacingTuningFrame;
  ModuleOptimalSpacings moduleOptimalSpacings;
  std::map<std::string, double> triggerRangeLowLimit;
  std::map<std::string, double> triggerRangeHighLimit;
  std::map<int, TGraphErrors> spacingTuningGraphs; // TODO: find a way to communicate the limits, not their plots!
  std::map<int, TGraphErrors> spacingTuningGraphsBad; // TODO: find a way to communicate the limits, not their plots!

  TriggerDistanceTuningPlotsVisitor(profileBag& myProfileBag,
                                    const std::vector<double>& triggerMomenta) :
    myProfileBag_(myProfileBag),
    triggerMomenta_(triggerMomenta) 
  {

    /************************************************/

    // TODO: clear only the relevant ones?
    myProfileBag_.clearTriggerNamedProfiles();
    optimalSpacingDistribution.SetName("optimalSpacing");
    optimalSpacingDistribution.SetTitle("Optimal spacing [default window]");
    optimalSpacingDistribution.SetXTitle("Spacing [mm]");
    optimalSpacingDistribution.SetYTitle("# modules");
    optimalSpacingDistribution.SetBins(100, 0.5, 6);
    optimalSpacingDistribution.Reset();

    optimalSpacingDistributionAW.SetName("optimalSpacingAW");
    optimalSpacingDistributionAW.SetTitle("Optimal spacing [actual window]");
    optimalSpacingDistributionAW.SetXTitle("Spacing [mm]");
    optimalSpacingDistributionAW.SetYTitle("# modules");
    optimalSpacingDistributionAW.SetBins(100, 0.5, 6);
    optimalSpacingDistributionAW.Reset();



  }

  void visit(const DetectorModule& m) {
    if (m.sensorLayout() && m.dsDistance() > 0.) foundSpacing_.insert(m.dsDistance());
  }

  void visit(const BarrelModule& aModule) {

    std::string myName = aModule.cntName() + "_L" + any2str(aModule.layer()) + "R" + any2str(aModule.ring());

    std::vector<const DetectorModule*>& theseBarrelModules = selectedModules_[myName];

    if (aModule.dsDistance() <= 0 || aModule.triggerWindow() == 0) return;

    // Prepare the variables to hold the profiles
    std::map<double, TProfile>& tuningProfiles = myProfileBag_.getNamedProfiles(profileBag::TriggerProfileName + myName);
    // Prepare the variables to hold the turn-on curve profiles
    std::map<double, TProfile>& turnonProfiles = myProfileBag_.getNamedProfiles(profileBag::TurnOnCurveName + myName);

    //  Profiles
    if (!preparedProfiles_[myName]) {
      preparedProfiles_[myName] = true;
      for (std::vector<double>::const_iterator it=triggerMomenta_.begin(); it!=triggerMomenta_.end(); ++it) {
        std::ostringstream tempSS; tempSS << "Trigger efficiency for " << myName.c_str() << ";Sensor spacing [mm];Efficiency [%]";
        tuningProfiles[*it].SetTitle(tempSS.str().c_str());
        tempSS.str(""); tempSS << "TrigEff" << myName.c_str() << "_" << (*it) << "GeV";
        tuningProfiles[*it].SetName(tempSS.str().c_str());
        tuningProfiles[*it].SetBins(100, 0.5, 6); // TODO: these numbers should go into some kind of const
      }     
    }

    // Turn-on curve
    if (!preparedTurnOn_[myName]) {
      preparedTurnOn_[myName] = true;
      for (unsigned int iWindow=0; iWindow<nWindows_; ++iWindow) {
        double windowSize=iWindow*2+1;
        std::ostringstream tempSS; tempSS << "Trigger efficiency for " << myName.c_str() << ";p_{T} [GeV/c];Efficiency [%]";
        turnonProfiles[windowSize].SetTitle(tempSS.str().c_str());
        tempSS.str(""); tempSS << "TurnOn" << myName.c_str() << "_window" << int(windowSize);
        turnonProfiles[windowSize].SetName(tempSS.str().c_str());
        turnonProfiles[windowSize].SetBins(100, 0.5, 10); // TODO: these numbers should go into some kind of const
      }     
    }


    XYZVector center = aModule.center();
    if ((center.Z()<0) || (center.Phi()<0) || (center.Phi()>M_PI/2)) return;
    theseBarrelModules.push_back(&aModule);

    PtErrorAdapter pterr(aModule);

    // Fill the tuning profiles for the windows actually set
    for (double dist=0.5; dist<=6; dist+=0.02) {
      for (std::vector<double>::const_iterator it=triggerMomenta_.begin(); it!=triggerMomenta_.end(); ++it) {
        double myPt = (*it);
        double myValue = 100 * pterr.getTriggerProbability(myPt, dist);
        if ((myValue>=0) && (myValue<=100))
          tuningProfiles[myPt].Fill(dist, myValue);
      }
    }

    // Fill the turnon curves profiles for the distance actually set
    for (double myPt=0.5; myPt<=10; myPt+=0.02) {
      for (unsigned int iWindow=0; iWindow<nWindows_; ++iWindow) {
        double windowSize=iWindow*2+1;
        double distance = aModule.dsDistance();
        double myValue = 100 * pterr.getTriggerProbability(myPt, distance, int(windowSize));
        if ((myValue>=0) && (myValue<=100))
          turnonProfiles[windowSize].Fill(myPt, myValue);
      }
    }
  }


  void visit(const EndcapModule& aModule) {



    string myBaseName = any2str(aModule.cntName()) + "_D" + any2str(aModule.disk());

    if ((aModule.dsDistance()<=0) || (aModule.triggerWindow()==0)) return;

    XYZVector center = aModule.center();
    // std::cerr << myBaseName << " z=" << center.Z() << ", Phi=" << center.Phi() << ", Rho=" << center.Rho() << std::endl; // debug
    if ((center.Z()<0) || (center.Phi()<0) || (center.Phi()>M_PI/2)) return;


    string myName = myBaseName + "R" + any2str(aModule.ring());
    ModuleVector& theseEndcapModules = selectedModules_[myName];
    theseEndcapModules.push_back(&aModule);

              // Prepare the variables to hold the profiles
    std::map<double, TProfile>& tuningProfiles = myProfileBag_.getNamedProfiles(profileBag::TriggerProfileName + myName);
    // Prepare the variables to hold the turn-on curve profiles
    std::map<double, TProfile>& turnonProfiles = myProfileBag_.getNamedProfiles(profileBag::TurnOnCurveName + myName);

    // Tuning profile
    if (!preparedProfiles_[myName]) {
      preparedProfiles_[myName] = true;
      for (std::vector<double>::const_iterator it=triggerMomenta_.begin(); it!=triggerMomenta_.end(); ++it) {
        std::ostringstream tempSS; tempSS << "Trigger efficiency for " << myName.c_str() << " GeV;Sensor spacing [mm];Efficiency [%]";
        tuningProfiles[*it].SetTitle(tempSS.str().c_str());
        tempSS.str(""); tempSS << "TrigEff" << myName.c_str() << "_" << (*it);
        tuningProfiles[*it].SetName(tempSS.str().c_str());
        tuningProfiles[*it].SetBins(100, 0.5, 6); // TODO: these numbers should go into some kind of const
      }
    }

    // Turn-on curve
    if (!preparedTurnOn_[myName]) {
      preparedTurnOn_[myName] = true;
      for (unsigned int iWindow=0; iWindow<nWindows_; ++iWindow) {
        double windowSize=iWindow*2+1;
        std::ostringstream tempSS; tempSS << "Trigger efficiency for " << myName.c_str() << ";p_{T} [GeV/c];Efficiency [%]";
        turnonProfiles[windowSize].SetTitle(tempSS.str().c_str());
        tempSS.str(""); tempSS << "TurnOn" << myName.c_str() << "_window" << windowSize;
        turnonProfiles[windowSize].SetName(tempSS.str().c_str());
        turnonProfiles[windowSize].SetBins(100, 0.5, 10); // TODO: these numbers should go into some kind of const
      }     
    }

    PtErrorAdapter pterr(aModule);

    // Fill the tuning profiles for the windows actually set
    for (double dist=0.5; dist<=6; dist+=0.02) {
      for (std::vector<double>::const_iterator it=triggerMomenta_.begin(); it!=triggerMomenta_.end(); ++it) {
        double myPt = (*it);
        double myValue = 100 * pterr.getTriggerProbability(myPt, dist);
        if ((myValue>=0) && (myValue<=100))
          tuningProfiles[myPt].Fill(dist, myValue);
      }
    }

    // Fill the turnon curves profiles for the distance actually set
    for (double myPt=0.5; myPt<=10; myPt+=0.02) {
      for (unsigned int iWindow=0; iWindow<nWindows_; ++iWindow) {
        double windowSize=iWindow*2+1;
        double distance = aModule.dsDistance();
        double myValue = 100 * pterr.getTriggerProbability(myPt, distance, int(windowSize));
        if ((myValue>=0) && (myValue<=100))
          turnonProfiles[windowSize].Fill(myPt, myValue);
      }
    }
  }

  void postVisit() {
      // TODO: put also the limits into a configurable parameter
      // Scan again over the plots I just made in order to find the
      // interesting range, if we have 1 and 2 in the range (TODO: find
      // a better way to find the interesting margins)

      // Run once per possible position in the tracker
      
    std::vector<double> spacingOptions(foundSpacing_.begin(), foundSpacing_.end());
    foundSpacing_.clear();    
     
    unsigned int nSpacingOptions = spacingOptions.size();             // TODO: keep this here!!

    std::pair<double, double> spacingTuningMomenta;
    spacingTuningMomenta.first = 1.;
    spacingTuningMomenta.second = 2.5;

    std::vector<std::string> profileNames = myProfileBag_.getProfileNames(profileBag::TriggerProfileName);
    for (std::vector<std::string>::const_iterator itName=profileNames.begin(); itName!=profileNames.end(); ++itName) {
      std::map<double, TProfile>& tuningProfiles = myProfileBag_.getNamedProfiles(*itName);
      TProfile& lowTuningProfile = tuningProfiles[spacingTuningMomenta.first];
      TProfile& highTuningProfile = tuningProfiles[spacingTuningMomenta.second];
      triggerRangeLowLimit[*itName] = findXThreshold(lowTuningProfile, 1, true);
      triggerRangeHighLimit[*itName] = findXThreshold(highTuningProfile, 90, false);
    }


    // Now loop over the selected modules and build the curves for the
    // hi-lo thingy (sensor spacing tuning)
    int windowSize;
    XYZVector center;
    double myPt;
    TProfile tempProfileLow("tempProfileLow", "", 100, 0.5, 6); // TODO: these numbers should go into some kind of const
    TProfile tempProfileHigh("tempProfileHigh", "", 100, 0.5, 6); // TODO: these numbers should go into some kind of const

    // TODO: IMPORTANT!!!!!! clear the spacing tuning graphs and frame here
    spacingTuningFrame.SetBins(selectedModules_.size(), 0, selectedModules_.size());
    spacingTuningFrame.SetYTitle("Optimal distance range [mm]");
    spacingTuningFrame.SetMinimum(0);
    spacingTuningFrame.SetMaximum(6);
    TAxis* xAxis = spacingTuningFrame.GetXaxis();

    int iType=0;
    std::map<double, bool> availableThinkness;
    // Loop over the selected module types
    for(std::map<std::string, ModuleVector>::iterator itTypes = selectedModules_.begin();
        itTypes!=selectedModules_.end(); ++itTypes) {
      const std::string& myName = itTypes->first;
      const ModuleVector& myModules = itTypes->second;
      xAxis->SetBinLabel(iType+1, myName.c_str());

      // Loop over the possible search windows
      for (unsigned int iWindow = 0; iWindow<nWindows_; ++iWindow) {
        windowSize = 1 + iWindow * 2;
        // Loop over the modules of type myName
        for (ModuleVector::const_iterator itModule = myModules.begin(); itModule!=myModules.end(); ++itModule) {
          const DetectorModule* aModule = (*itModule);
          PtErrorAdapter pterr(*aModule);
          // Loop over the possible distances
          double minDistBelow = 0.;
          availableThinkness[aModule->dsDistance()] = true;
          for (double dist=0.5; dist<=6; dist+=0.02) { // TODO: constant here
            // First with the high momentum
            myPt = (spacingTuningMomenta.second);
            double myValue = 100 * pterr.getTriggerProbability(myPt, dist, windowSize);
            if ((myValue>=0)&&(myValue<=100))
              tempProfileHigh.Fill(dist, myValue);
            // Then with low momentum
            myPt = (spacingTuningMomenta.first);
            myValue = 100 * pterr.getTriggerProbability(myPt, dist, windowSize);
            if ((myValue>=0)&&(myValue<=100))
              tempProfileLow.Fill(dist, myValue);
            if (myValue>1) minDistBelow = dist;
          }
          if (minDistBelow>=0) {
            if (windowSize==5) optimalSpacingDistribution.Fill(minDistBelow);
            if (windowSize==aModule->triggerWindow()) optimalSpacingDistributionAW.Fill(minDistBelow);
          }

          /*if ((myName=="ENDCAP_D02R05")&&(windowSize==5)) { // debug
            std::cout << myName << " - " << aModule->getTag() << " - minDistBelow = " << minDistBelow;
            std::cout <<  ", so[0]=" << spacingOptions[0] << ", so[n-1]=" << spacingOptions[nSpacingOptions-1];
            }*/
          if (minDistBelow<spacingOptions[0]) minDistBelow=spacingOptions[0];
          else if (minDistBelow>spacingOptions[nSpacingOptions-1]) minDistBelow=spacingOptions[nSpacingOptions-1];
          else {
            for (unsigned int iSpacing = 0; iSpacing < nSpacingOptions-1; ++iSpacing) {
              /*if ((myName=="ENDCAP_D02R05")&&(windowSize==5)) {// debug
                std::cout << " spacingOptions[" << iSpacing << "] = " << spacingOptions[iSpacing];
                std::cout << " spacingOptions[" << iSpacing+1 << "] = " << spacingOptions[iSpacing+1];
                }*/
              if ((minDistBelow>=spacingOptions[iSpacing]) && (minDistBelow<spacingOptions[iSpacing+1])) {
                minDistBelow=spacingOptions[iSpacing+1];
                /*if ((myName=="ENDCAP_D02R05")&&(windowSize==5)) // debug
                  std::cout << " here it is: ";*/
                break;
              }
            }
          }
          /*if ((myName=="ENDCAP_D02R05")&&(windowSize==5)) // debug
            std::cout << " - approx to " << minDistBelow << " for a window of " << windowSize << std::endl;*/
          moduleOptimalSpacings[aModule][windowSize] = minDistBelow;
        }
        // Find the "high" and "low" points
        double lowEdge = findXThreshold(tempProfileLow, 1, true);
        double highEdge = findXThreshold(tempProfileHigh, 90, false);
        // std::cerr << myName << ": " << lowEdge << " -> " << highEdge << std::endl; // debug
        double centerX; double sizeX;
        centerX = iType+(double(iWindow)+0.5)/(double(nWindows_));
        sizeX = 1./ (double(nWindows_)) * 0.8; // 80% of available space, so that they will not touch
        if (lowEdge<highEdge) {
          spacingTuningGraphs[iWindow].SetPoint(iType, centerX, (highEdge+lowEdge)/2.);
          spacingTuningGraphs[iWindow].SetPointError(iType, sizeX/2., (highEdge-lowEdge)/2.);
        } else {
          spacingTuningGraphsBad[iWindow].SetPoint(iType, centerX, (highEdge+lowEdge)/2.);
          spacingTuningGraphsBad[iWindow].SetPointError(iType, sizeX/2., (highEdge-lowEdge)/2.);
        }
        tempProfileLow.Reset();
        tempProfileHigh.Reset();
      }
      iType++;
    }

    // TODO: properly reset this!
    TGraphErrors& antani = spacingTuningGraphs[-1];
    int iPoints=0;
    for (std::map<double, bool>::iterator it = availableThinkness.begin(); it!= availableThinkness.end(); ++it) {
      iPoints++;
      antani.SetPoint(iPoints, selectedModules_.size()/2., it->first);
      antani.SetPointError(iPoints, selectedModules_.size()/2., 0);
      iPoints++;
    }
  }
};


#endif

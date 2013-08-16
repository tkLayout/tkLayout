/**
 * @file Analyzer.cc
 * @brief This is the implementation of the class that analyses a material budget
 */
#include <TH1D.h>
#include <TH2D.h>
#include <Analyzer.h>
#include <TProfile.h>
#include <TLegend.h>
#include <Palette.h>
#include <TTree.h>

#undef MATERIAL_SHADOW


namespace insur {

  int Analyzer::bsCounter =0;

  template<> void SummaryTable::setCell<std::string>(const int row, const int column, const std::string& content) {
    if (column > 0 && !hasCell(0, column)) summaryTable[std::make_pair(0, column)] = any2str(column + columnOffset_);
    if (row > 0 && !hasCell(row, 0)) summaryTable[std::make_pair(row, 0)] = any2str(row + rowOffset_);
    summaryTable[std::make_pair(row, column)]=content;
    numRows_ = row+1 > numRows_ ? row+1 : numRows_;
    numColumns_ = column+1 > numColumns_ ? column+1 : numColumns_;
  }

  template<> void SummaryTable::setSummaryCell<std::string>(std::string label, const std::string& content) {
    if (!hasSummaryCell()) {
      if (numRows_ > 2 && numColumns_ > 2) {
        summaryLabelPosition_ = std::make_pair(numRows_, 0); 
        summaryCellPosition_ = std::make_pair(numRows_++, numColumns_++);
      } else if (numRows_ >= 2 && numColumns_ == 2) {
        summaryLabelPosition_ = std::make_pair(numRows_, 0); 
        summaryCellPosition_ = std::make_pair(numRows_++, 1);
      } else if (numRows_ == 2 && numColumns_ > 2) {
        summaryLabelPosition_ = std::make_pair(0, numColumns_); 
        summaryCellPosition_ = std::make_pair(1, numColumns_++);
      }
    }
    summaryTable[summaryLabelPosition_] = label;
    summaryTable[summaryCellPosition_] = content;
  }

  const double graphBag::Triggerable     = 0.;
  const int graphBag::RhoGraph         = 0x001;
  const int graphBag::PhiGraph         = 0x002;
  const int graphBag::DGraph           = 0x003;
  const int graphBag::CtgthetaGraph    = 0x004;
  const int graphBag::Z0Graph          = 0x005;
  const int graphBag::PGraph           = 0x006;
  const int graphBag::TriggeredGraph   = 0x007;
  const int graphBag::IdealGraph       = 0x010;
  const int graphBag::RealGraph        = 0x020;
  const int graphBag::TriggerGraph     = 0x040;
  const int graphBag::StandardGraph    = 0x080;
  //const int graphBag::TriggerCorrelationGraph = 0x100;

  const int mapBag::efficiencyMap         = 0x001;
  const int mapBag::thresholdMap          = 0x002;
  const int mapBag::thicknessMap          = 0x004;
  const int mapBag::windowMap             = 0x008;
  const int mapBag::suggestedSpacingMap   = 0x010;
  const int mapBag::suggestedSpacingMapAW = 0x020;
  const int mapBag::nominalCutMap      = 0x040;
  const int mapBag::irradiatedPowerConsumptionMap = 0x080;
  const int mapBag::totalPowerConsumptionMap = 0x100;
  const int mapBag::moduleConnectionEtaMap = 0x200;
  const int mapBag::moduleConnectionPhiMap = 0x400;
  const int mapBag::moduleConnectionEndcapPhiMap = 0x800;

  const double profileBag::Triggerable    = 0.;
  const int profileBag::TriggeredProfile  = 0x0000007;
  const int profileBag::TriggeredFractionProfile  = 0x0000008;
  const int profileBag::TriggerPurityProfile = 0x0000010;
  const int profileBag::TriggerProfile    = 0x0000040;

  // These strings should be different from one another
  // Also one should never be a substring of the other
  const std::string profileBag::TriggerProfileName = "trigger";
  const std::string profileBag::TriggerProfileNameWindow = "windowTrigger";
  const std::string profileBag::TurnOnCurveName = "turnOnCurveTrigger";

  const double mapBag::dummyMomentum = 0.;

  int graphBag::clearTriggerGraphs() {
    return clearGraphs(graphBag::TriggerGraph);
  }

  int graphBag::clearStandardGraphs() {
    return clearGraphs(graphBag::StandardGraph);
  }

  int graphBag::clearGraphs(const int& attributeMask) {
    std::map<int, std::map<double, TGraph> >::iterator it;
    std::map<int, std::map<double, TGraph> >::iterator nextIt;

    int deleteCounter = 0;
    int anAttribute;
    for (it=graphMap_.begin(); it!=graphMap_.end(); ) {
      anAttribute=it->first;
      if ((anAttribute&attributeMask)==attributeMask) {
        nextIt = ((++it)--);
        graphMap_.erase(it);
        it=nextIt;
        ++deleteCounter;
      } else {
        ++it;
      }
    }
    return deleteCounter;
  }

  std::map<double, TGraph>& graphBag::getGraphs(const int& attribute) {
    return graphMap_[attribute];
  }

  std::map<double, TH2D>& mapBag::getMaps(const int& attribute) {
    return mapMap_[attribute];
  }

  int mapBag::clearMaps(const int& attributeMask) {
    std::map<int, std::map<double, TH2D> >::iterator it;
    std::map<int, std::map<double, TH2D> >::iterator nextIt;

    int deleteCounter = 0;
    int anAttribute;
    for (it=mapMap_.begin(); it!=mapMap_.end(); ) {
      anAttribute=it->first;
      if ((anAttribute&attributeMask)==attributeMask) {
        nextIt = it;
        ++nextIt;
        mapMap_.erase(it);
        it=nextIt;
        ++deleteCounter;
      } else {
        ++it;
      }
    }
    return deleteCounter;
  }

  int profileBag::clearTriggerProfiles() {
    return clearProfiles(profileBag::TriggerProfile);
  }

  int profileBag::clearTriggerNamedProfiles() {
    return clearNamedProfiles(profileBag::TriggerProfileName);
  }


  std::map<double, TProfile>& profileBag::getProfiles(const int& attribute) {
    return profileMap_[attribute];
  }

  // TODO: this looks like an invitation to use template classes Will
  // do as soon as I have time :D (copy-paste worked till now...)
  int profileBag::clearProfiles(const int& attributeMask) {
    std::map<int, std::map<double, TProfile> >::iterator it;
    std::map<int, std::map<double, TProfile> >::iterator nextIt;

    int deleteCounter = 0;
    int anAttribute;
    for (it=profileMap_.begin(); it!=profileMap_.end(); ) {
      anAttribute=it->first;
      if ((anAttribute&attributeMask)==attributeMask) {
        nextIt = it;
        ++nextIt;
        profileMap_.erase(it);
        it=nextIt;
        ++deleteCounter;
      } else {
        ++it;
      }
    }
    return deleteCounter;
  }

  int profileBag::clearNamedProfiles(const std::string& name) {
    std::map<std::string, std::map<double, TProfile> >::iterator it;
    std::map<std::string, std::map<double, TProfile> >::iterator nextIt;

    int deleteCounter = 0;
    std::string anAttribute;
    for (it=namedProfileMap_.begin(); it!=namedProfileMap_.end(); ) {
      anAttribute=it->first;
      if (anAttribute.substr(0, name.size())==name) {
        nextIt = it;
        ++nextIt;
        namedProfileMap_.erase(it);
        it=nextIt;
        ++deleteCounter;
      } else {
        ++it;
      }
    }
    return deleteCounter;    
  }

  std::vector<std::string> profileBag::getProfileNames(const std::string& name) {
    std::vector<std::string> result;

    std::map<std::string, std::map<double, TProfile> >::iterator it;

    std::string anAttribute;
    for (it=namedProfileMap_.begin(); it!=namedProfileMap_.end(); ++it) {
      anAttribute=it->first;
      if (anAttribute.substr(0, name.size())==name)
        result.push_back(anAttribute);
    }
    return result;    
  }


  std::map<double, TProfile>& profileBag::getNamedProfiles(const std::string& name) {
    return namedProfileMap_[name];
  }

  const double Analyzer::ZeroHitsRequired = 0;
  const double Analyzer::OneHitRequired = 0.0001;





  // public
  /**
   * A comparison function for the first elements in two pairs of integers.
   * @param p The first pair
   * @param q The second pair
   * @return True if <i>p.first</i> is smaller than <i>q.first</i>, false otherwise
   */
  bool compareIntPairFirst(std::pair<int, int> p, std::pair<int, int> q) {
    return (p.first < q.first);
  }

  /**
   * A comparison function for the second elements in two pairs of integers.
   * @param p The first pair
   * @param q The second pair
   * @return True if <i>p.second</i> is smaller than <i>q.second</i>, false otherwise
   */
  bool compareIntPairSecond(std::pair<int, int> p, std::pair<int, int> q) {
    return (p.second < q.second);
  }

  /**
   * The constructor sets a number of internal constants.
   */
  Analyzer::Analyzer() {
    // Not strictly necessary, but it's useful to keep
    // the color the same for the most used module types
    lastPickedColor = 1;
    //colorPicker("ptOut");
    //colorPicker("rphi");
    //colorPicker("stereo");
    //colorPicker("ptIn");
    geomLite = NULL;   geomLiteCreated=false;
    geomLiteXY = NULL; geomLiteXYCreated=false;
    geomLiteYZ = NULL; geomLiteYZCreated=false;
    geomLiteEC = NULL; geomLiteECCreated=false;
    geometryTracksUsed = 0;
    materialTracksUsed = 0;

    //etaMaxMaterial = 3.1;
    etaMaxGeometry = 2.6;

    trackingCuts.push_back(0.01);
    triggerCuts.push_back(0.01);
    double detaTrack = 0.8;
    //double detaTrig = 0.8; // This used to be 0.7
    addCut("C", detaTrack, detaTrack);
    addCut("I", detaTrack*2, detaTrack*2);
    addCut("F", detaTrack*3, 2.1);
    addCut("VF",detaTrack*4, 2.25);
    addCut("VF",detaTrack*5, 2.4);
  }

  // public
  // TODO: add documentation for this function
  void Analyzer::addCut(const std::string& cutName, const double& trackingCut, const double& triggerCut) {
    cutNames.push_back(cutName);
    trackingCuts.push_back(trackingCut);
    triggerCuts.push_back(triggerCut);
  }

  double Analyzer::getEtaMaxTracking() {
    return trackingCuts[trackingCuts.size()-1];
  }

  double Analyzer::getEtaMaxTrigger() {
    return triggerCuts[triggerCuts.size()-1];
  }


  // private
  /* High-level function finding all hits for a given tracker (and pixel)
   * and adding them to the track. The total crossed material is returned.
   * @param mb A reference to the instance of <i>MaterialBudget</i> that is to be analysed
   * @param momenta A list of momentum values for which to perform the efficiency measurements
   * @param etaSteps The number of wedges in the fan of tracks covered by the eta scan
   * @param pm A pointer to a second material budget associated to a pixel detector; may be <i>NULL</i>
   * @return the total crossed material amount
   */
  Material Analyzer::findAllHits(MaterialBudget& mb, MaterialBudget* pm, 
                                 double& eta, double& theta, double& phi, Track& track) {
    Material totalMaterial;
    //      active volumes, barrel
    totalMaterial  = findHitsModules(mb.getBarrelModuleCaps(), eta, theta, phi, track);
    //      active volumes, endcap
    totalMaterial += findHitsModules(mb.getEndcapModuleCaps(), eta, theta, phi, track);
    //      services, barrel
    totalMaterial += findHitsInactiveSurfaces(mb.getInactiveSurfaces().getBarrelServices(), eta, theta, track);
    //      services, endcap
    totalMaterial += findHitsInactiveSurfaces(mb.getInactiveSurfaces().getEndcapServices(), eta, theta, track);
    //      supports
    totalMaterial += findHitsInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track);
    //      pixels, if they exist
    if (pm != NULL) {
      totalMaterial += findHitsModules(pm->getBarrelModuleCaps(), eta, theta, phi, track, true);
      totalMaterial += findHitsModules(pm->getEndcapModuleCaps(), eta, theta, phi, track, true);
      totalMaterial += findHitsInactiveSurfaces(pm->getInactiveSurfaces().getBarrelServices(), eta, theta, track, true);
      totalMaterial += findHitsInactiveSurfaces(pm->getInactiveSurfaces().getEndcapServices(), eta, theta, track, true);
      totalMaterial += findHitsInactiveSurfaces(pm->getInactiveSurfaces().getSupports(), eta, theta, track, true);
    }
    return totalMaterial;
  }


  /**
   * The analysis function performing all necessary the trigger resolution
   * @param mb A reference to the instance of <i>MaterialBudget</i> that is to be analysed
   * @param momenta A list of momentum values for which to perform the efficiency measurements
   * @param etaSteps The number of wedges in the fan of tracks covered by the eta scan
   * @param pm A pointer to a second material budget associated to a pixel detector; may be <i>NULL</i>
   */
  void Analyzer::analyzeTrigger(MaterialBudget& mb,
                                const std::vector<double>& momenta,
                                const std::vector<double>& triggerMomenta,
                                const std::vector<double>& thresholdProbabilities,
                                int etaSteps,
                                MaterialBudget* pm) {

    Tracker& tracker = mb.getTracker();
    double efficiency = tracker.getEfficiency();

    materialTracksUsed = etaSteps;

    int nTracks;
    double etaStep, eta, theta, phi;

    // prepare etaStep, phiStep, nTracks, nScans
    if (etaSteps > 1) etaStep = getEtaMaxTrigger() / (double)(etaSteps - 1);
    else etaStep = getEtaMaxTrigger();
    nTracks = etaSteps;

    // prepareTriggerPerformanceHistograms(nTracks, getEtaMaxTrigger(), triggerMomenta, thresholdProbabilities);

    // reset the list of tracks
    std::vector<Track> tv;
    std::vector<Track> tvIdeal;

    // used fixed phi
    //phi = PI / 2.0;

    // Loop over nTracks (eta range [0, getEtaMaxTrigger()])
    for (int i_eta = 0; i_eta < nTracks; i_eta++) {
      phi = myDice.Rndm() * PI * 2.0;
      Material tmp;
      Track track;
      eta = i_eta * etaStep;
      theta = 2 * atan(exp(-eta));
      track.setTheta(theta);      
      track.setPhi(phi);

      tmp = findAllHits(mb, pm, eta, theta, phi, track);

      // Debug: material amount
      //std::cerr << "eta = " << eta
      //        << ", material.radiation = " << tmp.radiation
      //        << ", material.interaction = " << tmp.interaction
      //        << std::endl;

      // TODO: add the beam pipe as a user material eveywhere!
      // in a coherent way
      // Add the hit on the beam pipe
      Hit* hit = new Hit(23./sin(theta));
      hit->setOrientation(Hit::Horizontal);
      hit->setObjectKind(Hit::Inactive);
      Material beamPipeMat;
      beamPipeMat.radiation = 0.0023 / sin(theta);
      beamPipeMat.interaction = 0.0019 / sin(theta);
      hit->setCorrectedMaterial(beamPipeMat);
      track.addHit(hit);

      if (!track.noHits()) {
        // Keep only triggering hits
        // std::cerr << "Material before = " << track.getCorrectedMaterial().radiation;
        track.keepTriggerOnly();
        // std::cerr << " material after = " << track.getCorrectedMaterial().radiation << std::endl;

        // Assume the IP constraint here
        // TODO: make this configurable
        if (tracker.getUseIPConstraint())
          track.addIPConstraint(tracker.getRError(),tracker.getZErrorCollider());
        track.sort();
        track.setTriggerResolution(true);

        if (efficiency!=1) track.addEfficiency(efficiency, false);
        if (track.nActiveHits(true)>2) { // At least 3 points are needed to measure the arrow
          track.computeErrors(momenta);
          tv.push_back(track);

          Track trackIdeal = track;
          trackIdeal.removeMaterial();
          trackIdeal.computeErrors(momenta);
          tvIdeal.push_back(trackIdeal);    
        }    
      }
    }

    // Compute the number of triggering points along the selected tracks
    // fillTriggerEfficiencyGraphs(triggerMomenta, tv);

    // Fill the trigger performance maps
    // fillTriggerPerformanceMaps(tracker);

    calculateGraphs(momenta, tv, graphBag::TriggerGraph | graphBag::RealGraph);
    calculateGraphs(momenta, tvIdeal, graphBag::TriggerGraph | graphBag::IdealGraph);
  }




  /**
   * The analysis function performing all necessary the trigger efficiencies
   * @param tracker 
   * @param thresholdProbabilities
   * @param etaSteps The number of wedges in the fan of tracks covered by the eta scan
   */
  void Analyzer::analyzeTriggerEfficiency(Tracker& tracker,
                                          const std::vector<double>& triggerMomenta,
                                          const std::vector<double>& thresholdProbabilities,
                                          int etaSteps) {

    double efficiency = tracker.getEfficiency();

    materialTracksUsed = etaSteps;

    int nTracks;
    double etaStep, z0, eta, theta, phi;
    double zErrorCollider = tracker.getZErrorCollider();

    // prepare etaStep, phiStep, nTracks, nScans
    if (etaSteps > 1) etaStep = getEtaMaxTrigger() / (double)(etaSteps - 1);
    else etaStep = getEtaMaxTrigger();
    nTracks = etaSteps;

    prepareTriggerPerformanceHistograms(nTracks, getEtaMaxTrigger(), triggerMomenta, thresholdProbabilities);

    // reset the list of tracks
    std::vector<Track> tv;

    // Loop over nTracks (eta range [0, getEtaMaxTrigger()])
    for (int i_eta = 0; i_eta < nTracks; i_eta++) {
      phi = myDice.Rndm() * PI * 2.0;
      z0 = myDice.Gaus(0, zErrorCollider);
      int nHits;
      Track track;
      eta = i_eta * etaStep;
      theta = 2 * atan(exp(-eta));
      track.setTheta(theta);      
      track.setPhi(phi);

      nHits = findHitsModules(tracker, z0, eta, theta, phi, track);

      if (nHits) {
        // Keep only triggering hits
        // std::cerr << "Material before = " << track.getCorrectedMaterial().radiation;
        track.keepTriggerOnly();
        track.sort();
        track.setTriggerResolution(true);

        // std::cerr << " material after = " << track.getCorrectedMaterial().radiation << std::endl;

        if (efficiency!=1) track.addEfficiency(efficiency, false);
        if (track.nActiveHits(true)>0) { // At least 3 points are needed to measure the arrow
          tv.push_back(track);
        }    
      }
    }

    // Compute the number of triggering points along the selected tracks
    fillTriggerEfficiencyGraphs(triggerMomenta, tv, tracker);

    // Fill the trigger performance maps
    fillTriggerPerformanceMaps(tracker);

  }

  void Analyzer::fillAvailableSpacing(Tracker& tracker, std::vector<double>& spacingOptions) {
    LayerVector& layerSet = tracker.getLayers();
    LayerVector::iterator layIt;
    ModuleVector* moduleSet;
    ModuleVector::iterator modIt; 
    Layer* aLayer;
    Module* aModule;
    double aSpacing;
    std::set<double> foundSpacing;

    // Loop over all the layers
    for(layIt = layerSet.begin(); layIt!= layerSet.end(); ++layIt) {
      aLayer = (*layIt);
      // Actually scan the modules
      moduleSet = aLayer->getModuleVector();
      for(modIt = moduleSet->begin(); modIt != moduleSet->end(); ++modIt) {
        aModule = (*modIt);
        if (aModule->getReadoutType()==Module::Pt) {
          aSpacing = aModule->getStereoDistance();
          if (aSpacing>0) {
            foundSpacing.insert(aSpacing);
          }
        }
      }
    }

    spacingOptions.clear();
    for(std::set<double>::iterator it=foundSpacing.begin(); it!=foundSpacing.end(); ++it) {
      spacingOptions.push_back(*it);
    }

  }

  void Analyzer::createTriggerDistanceTuningPlots(Tracker& tracker, const std::vector<double>& triggerMomenta) {

    // TODO: put these in a configuration file somewhere
    /************************************************/
    std::pair<double, double> spacingTuningMomenta;
    std::vector<double> spacingOptions;
    fillAvailableSpacing(tracker, spacingOptions);
    //spacingOptions.push_back(1.1);
    //spacingOptions.push_back(1.8);
    //spacingOptions.push_back(2.5);
    //spacingOptions.push_back(5); // TODO: make these configurable
    //spacingOptions.push_back(0.8);
    //spacingOptions.push_back(1.6);
    //spacingOptions.push_back(3);
    //std::sort(spacingOptions.begin(), spacingOptions.end()); // TODO: keep this here!!
    unsigned int nSpacingOptions = spacingOptions.size();             // TODO: keep this here!!
    spacingTuningMomenta.first = 1.;
    spacingTuningMomenta.second = 2.5;
    const unsigned int nWindows = 5;
    /************************************************/

    // TODO: clear only the relevant ones?
    myProfileBag.clearTriggerNamedProfiles();
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


    std::string myName;
    std::string myBaseName;
    std::map<std::string, bool> preparedProfiles;
    std::map<std::string, bool> preparedTurnOn;

    // Loop over all the tracker
    LayerVector& layerSet = tracker.getLayers();
    LayerVector::iterator layIt;
    Layer* aLayer;
    ModuleVector::iterator modIt;
    ModuleVector* moduleSet;
    Module* aModule;
    double myValue;

    BarrelLayer* aBarrelLayer;
    EndcapLayer* anEndcapLayer;
    EndcapModule* anEndcapModule;
    std::ostringstream tempSS;


    std::map<std::string, ModuleVector> selectedModules;

    // Loop over all the layers
    for(layIt = layerSet.begin(); layIt!= layerSet.end(); ++layIt) {
      aLayer = (*layIt);
      //std::cerr << "layer " << aLayer->getIndex() << std::endl; // debug

      aBarrelLayer = dynamic_cast<BarrelLayer*>(aLayer);
      anEndcapLayer = dynamic_cast<EndcapLayer*>(aLayer);

      /****************/
      /* BARREL LAYER */
      /****************/
      if (aBarrelLayer) { 
        // If it's a barrel layer, scan over all its modules

        // Sort the plot by layer
        myName = aLayer->getContainerName() + "_" + aLayer->getName();
        ModuleVector& theseBarrelModules = selectedModules[myName];

        // Actually scan the modules
        moduleSet = aLayer->getModuleVector();
        for(modIt = moduleSet->begin(); modIt != moduleSet->end(); ++modIt) {
          aModule = (*modIt);

          if ((aModule->getStereoDistance()<=0) || (aModule->getTriggerWindow()==0)) continue;
          if (aModule->getReadoutType() != Module::Pt) {
            std::cerr << "WARNING: a non-pT module has a non-zero trigger window!"  << std::endl; // TODO: put this in the logger as a warning
            continue;
          }

          // Prepare the variables to hold the profiles
          std::map<double, TProfile>& tuningProfiles = myProfileBag.getNamedProfiles(profileBag::TriggerProfileName + myName);
          // Prepare the variables to hold the turn-on curve profiles
          std::map<double, TProfile>& turnonProfiles = myProfileBag.getNamedProfiles(profileBag::TurnOnCurveName + myName);

          //  Profiles
          if (!preparedProfiles[myName]) {
            preparedProfiles[myName] = true;
            for (std::vector<double>::const_iterator it=triggerMomenta.begin(); it!=triggerMomenta.end(); ++it) {
              tempSS.str(""); tempSS << "Trigger efficiency for " << myName.c_str() << ";Sensor spacing [mm];Efficiency [%]";
              tuningProfiles[*it].SetTitle(tempSS.str().c_str());
              tempSS.str(""); tempSS << "TrigEff" << myName.c_str() << "_" << (*it) << "GeV";
              tuningProfiles[*it].SetName(tempSS.str().c_str());
              tuningProfiles[*it].SetBins(100, 0.5, 6); // TODO: these numbers should go into some kind of const
            }     
          }

          // Turn-on curve
          if (!preparedTurnOn[myName]) {
            preparedTurnOn[myName] = true;
            for (unsigned int iWindow=0; iWindow<nWindows; ++iWindow) {
              double windowSize=iWindow*2+1;
              tempSS.str(""); tempSS << "Trigger efficiency for " << myName.c_str() << ";p_{T} [GeV/c];Efficiency [%]";
              turnonProfiles[windowSize].SetTitle(tempSS.str().c_str());
              tempSS.str(""); tempSS << "TurnOn" << myName.c_str() << "_window" << int(windowSize);
              turnonProfiles[windowSize].SetName(tempSS.str().c_str());
              turnonProfiles[windowSize].SetBins(100, 0.5, 10); // TODO: these numbers should go into some kind of const
            }     
          }


          XYZVector center = aModule->getMeanPoint();
          if ((center.Z()<0) || (center.Phi()<0) || (center.Phi()>M_PI/2)) continue;
          theseBarrelModules.push_back(aModule);

          // Fill the tuning profiles for the windows actually set
          for (double dist=0.5; dist<=6; dist+=0.02) {
            for (std::vector<double>::const_iterator it=triggerMomenta.begin(); it!=triggerMomenta.end(); ++it) {
              double myPt = (*it);
              myValue = 100 * aModule->getTriggerProbability(myPt, dist);
              if ((myValue>=0) && (myValue<=100))
                tuningProfiles[myPt].Fill(dist, myValue);
            }
          }

          // Fill the turnon curves profiles for the distance actually set
          for (double myPt=0.5; myPt<=10; myPt+=0.02) {
            for (unsigned int iWindow=0; iWindow<nWindows; ++iWindow) {
              double windowSize=iWindow*2+1;
              double distance = aModule->getStereoDistance();
              myValue = 100 * aModule->getTriggerProbability(myPt, distance, int(windowSize));
              if ((myValue>=0) && (myValue<=100))
                turnonProfiles[windowSize].Fill(myPt, myValue);
            }
          }
        }

        /****************/
        /* ENDCAP LAYER */
        /****************/
      } else if (anEndcapLayer) {
        // If it's an endcap layer, scan over all disk 1 modules
        //if (anEndcapLayer->getIndex() == 1 && (aLayer->getContainerName()!="")) {
        if (aLayer->getContainerName()!="") {

          // Sort the plot by ring (using only disk 1)
          tempSS.str(""); 
          tempSS << "_D";
          tempSS.width(2); tempSS.fill('0'); tempSS << anEndcapLayer->getIndex();
          myBaseName = aLayer->getContainerName() + tempSS.str();

          // Actually scan the modules
          moduleSet = aLayer->getModuleVector();
          for(modIt = moduleSet->begin(); modIt != moduleSet->end(); ++modIt) {
            aModule = (*modIt);
            if ((aModule->getStereoDistance()<=0) || (aModule->getTriggerWindow()==0)) continue;
            if (aModule->getReadoutType() != Module::Pt) {
              std::cerr << "WARNING: a non-pT module has a non-zero trigger window!"  << std::endl; // TODO: put this in the logger as a warning
              continue;
            }

            anEndcapModule = dynamic_cast<EndcapModule*>(aModule);
            XYZVector center = aModule->getMeanPoint();
            // std::cerr << myBaseName << " z=" << center.Z() << ", Phi=" << center.Phi() << ", Rho=" << center.Rho() << std::endl; // debug
            if ((center.Z()<0) || (center.Phi()<0) || (center.Phi()>M_PI/2)) continue;


            if (anEndcapModule) {
              tempSS.str("");
              tempSS << "R";
              tempSS.width(2); tempSS.fill('0');
              tempSS << anEndcapModule->getRing();
              myName = myBaseName + tempSS.str();
              ModuleVector& theseEndcapModules = selectedModules[myName];
              theseEndcapModules.push_back(aModule);

              // Prepare the variables to hold the profiles
              std::map<double, TProfile>& tuningProfiles = myProfileBag.getNamedProfiles(profileBag::TriggerProfileName + myName);
              // Prepare the variables to hold the turn-on curve profiles
              std::map<double, TProfile>& turnonProfiles = myProfileBag.getNamedProfiles(profileBag::TurnOnCurveName + myName);

              // Tuning profile
              if (!preparedProfiles[myName]) {
                preparedProfiles[myName] = true;
                for (std::vector<double>::const_iterator it=triggerMomenta.begin(); it!=triggerMomenta.end(); ++it) {
                  tempSS.str(""); tempSS << "Trigger efficiency for " << myName.c_str() << " GeV;Sensor spacing [mm];Efficiency [%]";
                  tuningProfiles[*it].SetTitle(tempSS.str().c_str());
                  tempSS.str(""); tempSS << "TrigEff" << myName.c_str() << "_" << (*it);
                  tuningProfiles[*it].SetName(tempSS.str().c_str());
                  tuningProfiles[*it].SetBins(100, 0.5, 6); // TODO: these numbers should go into some kind of const
                }
              }

              // Turn-on curve
              if (!preparedTurnOn[myName]) {
                preparedTurnOn[myName] = true;
                for (unsigned int iWindow=0; iWindow<nWindows; ++iWindow) {
                  double windowSize=iWindow*2+1;
                  tempSS.str(""); tempSS << "Trigger efficiency for " << myName.c_str() << ";p_{T} [GeV/c];Efficiency [%]";
                  turnonProfiles[windowSize].SetTitle(tempSS.str().c_str());
                  tempSS.str(""); tempSS << "TurnOn" << myName.c_str() << "_window" << windowSize;
                  turnonProfiles[windowSize].SetName(tempSS.str().c_str());
                  turnonProfiles[windowSize].SetBins(100, 0.5, 10); // TODO: these numbers should go into some kind of const
                }     
              }



              // Fill the tuning profiles for the windows actually set
              for (double dist=0.5; dist<=6; dist+=0.02) {
                for (std::vector<double>::const_iterator it=triggerMomenta.begin(); it!=triggerMomenta.end(); ++it) {
                  double myPt = (*it);
                  myValue = 100 * aModule->getTriggerProbability(myPt, dist);
                  if ((myValue>=0) && (myValue<=100))
                    tuningProfiles[myPt].Fill(dist, myValue);
                }
              }

              // Fill the turnon curves profiles for the distance actually set
              for (double myPt=0.5; myPt<=10; myPt+=0.02) {
                for (unsigned int iWindow=0; iWindow<nWindows; ++iWindow) {
                  double windowSize=iWindow*2+1;
                  double distance = aModule->getStereoDistance();
                  myValue = 100 * aModule->getTriggerProbability(myPt, distance, int(windowSize));
                  if ((myValue>=0) && (myValue<=100))
                    turnonProfiles[windowSize].Fill(myPt, myValue);
                }
              }

            } else {
              std::cerr << "ERROR: this should not happen: a not-endcap module was found in an endcap layer! Contact the developers" << std::endl;
            }
          }
        }
      }   
      }

      // TODO: put also the limits into a configurable parameter
      // Scan again over the plots I just made in order to find the
      // interesting range, if we have 1 and 2 in the range (TODO: find
      // a better way to find the interesting margins)

      // Run once per possible position in the tracker
      std::vector<std::string> profileNames = myProfileBag.getProfileNames(profileBag::TriggerProfileName);
      for (std::vector<std::string>::const_iterator itName=profileNames.begin(); itName!=profileNames.end(); ++itName) {
        std::map<double, TProfile>& tuningProfiles = myProfileBag.getNamedProfiles(*itName);
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
      spacingTuningFrame.SetBins(selectedModules.size(), 0, selectedModules.size());
      spacingTuningFrame.SetYTitle("Optimal distance range [mm]");
      spacingTuningFrame.SetMinimum(0);
      spacingTuningFrame.SetMaximum(6);
      TAxis* xAxis = spacingTuningFrame.GetXaxis();

      int iType=0;
      std::map<double, bool> availableThinkness;
      // Loop over the selected module types
      for(std::map<std::string, ModuleVector>::iterator itTypes = selectedModules.begin();
          itTypes!=selectedModules.end(); ++itTypes) {
        const std::string& myName = itTypes->first;
        const ModuleVector& myModules = itTypes->second;
        xAxis->SetBinLabel(iType+1, myName.c_str());

        // Loop over the possible search windows
        for (unsigned int iWindow = 0; iWindow<nWindows; ++iWindow) {
          windowSize = 1 + iWindow * 2;
          // Loop over the modules of type myName
          for (ModuleVector::const_iterator itModule = myModules.begin(); itModule!=myModules.end(); ++itModule) {
            aModule = (*itModule);
            // Loop over the possible distances
            double minDistBelow = 0.;
            availableThinkness[aModule->getStereoDistance()] = true;
            for (double dist=0.5; dist<=6; dist+=0.02) { // TODO: constant here
              // First with the high momentum
              myPt = (spacingTuningMomenta.second);
              myValue = 100 * aModule->getTriggerProbability(myPt, dist, windowSize);
              if ((myValue>=0)&&(myValue<=100))
                tempProfileHigh.Fill(dist, myValue);
              // Then with low momentum
              myPt = (spacingTuningMomenta.first);
              myValue = 100 * aModule->getTriggerProbability(myPt, dist, windowSize);
              if ((myValue>=0)&&(myValue<=100))
                tempProfileLow.Fill(dist, myValue);
              if (myValue>1) minDistBelow = dist;
            }
            if (minDistBelow>=0) {
              if (windowSize==5) optimalSpacingDistribution.Fill(minDistBelow);
              if (windowSize==aModule->getTriggerWindow()) optimalSpacingDistributionAW.Fill(minDistBelow);
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
            aModule->setOptimalSpacing(windowSize, minDistBelow);
          }
          // Find the "high" and "low" points
          double lowEdge = findXThreshold(tempProfileLow, 1, true);
          double highEdge = findXThreshold(tempProfileHigh, 90, false);
          // std::cerr << myName << ": " << lowEdge << " -> " << highEdge << std::endl; // debug
          double centerX; double sizeX;
          centerX = iType+(double(iWindow)+0.5)/(double(nWindows));
          sizeX = 1./ (double(nWindows)) * 0.8; // 80% of available space, so that they will not touch
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
        antani.SetPoint(iPoints, selectedModules.size()/2., it->first);
        antani.SetPointError(iPoints, selectedModules.size()/2., 0);
        iPoints++;
      }
    }


    double Analyzer::findXThreshold(const TProfile& aProfile, const double& yThreshold, const bool& goForward) {  
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




    void Analyzer::fillTriggerEfficiencyGraphs(const std::vector<double>& triggerMomenta,
                                               const std::vector<Track>& trackVector,
                                               const Tracker& tracker) {

      // Prepare the graphs to record the number of triggered points
      //std::map<double, TGraph>& trigGraphs = myGraphBag.getGraphs(graphBag::TriggerGraph|graphBag::TriggeredGraph);
      std::map<double, TProfile>& trigProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggeredProfile);
      std::map<double, TProfile>& trigFractionProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggeredFractionProfile);
      std::map<double, TProfile>& trigPurityProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggerPurityProfile);

      TProfile& totalProfile = trigProfiles[profileBag::Triggerable];

      double eta;

      for (std::vector<Track>::const_iterator itTrack = trackVector.begin();
           itTrack != trackVector.end(); ++itTrack) {
        const Track& myTrack=(*itTrack);

        eta = - log(tan(myTrack.getTheta() / 2));
        int nHits = myTrack.nActiveHits(false, false);
        totalProfile.Fill(eta, nHits);
        std::vector<Module*> hitModules = myTrack.getHitModules();

        for(std::vector<double>::const_iterator itMomentum = triggerMomenta.begin();
            itMomentum!=triggerMomenta.end(); ++itMomentum) {
          TProfile& myProfile = trigProfiles[(*itMomentum)];
          TProfile& myFractionProfile = trigFractionProfiles[(*itMomentum)];
          TProfile& myPurityProfile = trigPurityProfiles[(*itMomentum)];
          double nExpectedTriggerPoints = myTrack.expectedTriggerPoints(*itMomentum);
          if (nExpectedTriggerPoints>=0) { // sanity check (! nan)
            myProfile.Fill(eta, nExpectedTriggerPoints);
            if (nHits>0) {
              myFractionProfile.Fill(eta, nExpectedTriggerPoints*100/double(nHits));

              double curAvgTrue=0;
              double curAvgInteresting=0;
              double curAvgFake=0;
              double bgReductionFactor; // Reduction of the combinatorial background for ptMixed modules by turning off the appropriate pixels
              for (std::vector<Module*>::iterator itModule = hitModules.begin(); itModule!= hitModules.end(); ++itModule) {
                Module* hitModule = (*itModule);
                // Hits that we would like to have from tracks above this threshold
                curAvgInteresting += hitModule->getParticleFrequencyPerEventAbove(*itMomentum);
                // ... out of which we only see these
                curAvgTrue  += hitModule->getTriggerFrequencyTruePerEventAbove(*itMomentum);

                // The background is given by the contamination from low pT tracks...
                curAvgFake += hitModule->getTriggerFrequencyTruePerEventBelow(*itMomentum);
                // ... plus the combinatorial background from occupancy (can be reduced using ptMixed modules)
                if (hitModule->getType()=="ptMixed") bgReductionFactor = hitModule->getGeometricEfficiency(); else bgReductionFactor=1;
                curAvgFake += hitModule->getTriggerFrequencyFakePerEvent()*tracker.getNMB() * bgReductionFactor;
              }
              myPurityProfile.Fill(eta, 100*curAvgTrue/(curAvgTrue+curAvgFake));
            }
          }


        }
      }
      if (totalProfile.GetMaximum() < maximum_n_planes) totalProfile.SetMaximum(maximum_n_planes);

    }



    /**
     * The main analysis function provides a frame for the scan in eta, defers summing up the radiation
     * and interaction lengths for each volume category to subfunctions, and sorts those results into the
     * correct histograms.
     * @param mb A reference to the instance of <i>MaterialBudget</i> that is to be analysed
     * @param momenta A list of momentum values for the tracks that are shot through the layout
     * @param etaSteps The number of wedges in the fan of tracks covered by the eta scan
     * @param A pointer to a second material budget associated to a pixel detector; may be <i>NULL</i>
     */
    void Analyzer::analyzeMaterialBudget(MaterialBudget& mb, const std::vector<double>& momenta, int etaSteps,
                                         MaterialBudget* pm , bool computeResolution) {

      Tracker& tracker = mb.getTracker();
      double efficiency = tracker.getEfficiency();
      double pixelEfficiency = tracker.getPixelEfficiency();
      materialTracksUsed = etaSteps;
      int nTracks;
      double etaStep, eta, theta, phi;
      clearMaterialBudgetHistograms();
      clearCells();
      // prepare etaStep, phiStep, nTracks, nScans
      if (etaSteps > 1) etaStep = getEtaMaxMaterial() / (double)(etaSteps - 1);
      else etaStep = getEtaMaxMaterial();
      nTracks = etaSteps;
      // reset the number of bins and the histogram boundaries (0.0 to getEtaMaxMaterial()) for all histograms, recalculate the cell boundaries
      setHistogramBinsBoundaries(nTracks, 0.0, getEtaMaxMaterial());
      setCellBoundaries(nTracks, 0.0, outer_radius + volume_width, 0.0, getEtaMaxMaterial());
      // reset the list of tracks
      std::vector<Track> tv;
      std::vector<Track> tvIdeal;

      for (int i_eta = 0; i_eta < nTracks; i_eta++) {
        phi = myDice.Rndm() * PI * 2.0;
        Material tmp;
        Track track;
        eta = i_eta * etaStep;
        theta = 2 * atan(pow(E, -1 * eta)); // TODO: switch to exp() here
        track.setTheta(theta);
        track.setPhi(phi);
        //      active volumes, barrel
        std::map<std::string, Material> sumComponentsRI;
        tmp = analyzeModules(mb.getBarrelModuleCaps(), eta, theta, phi, track, sumComponentsRI);
        ractivebarrel.Fill(eta, tmp.radiation);
        iactivebarrel.Fill(eta, tmp.interaction);
        rbarrelall.Fill(eta, tmp.radiation);
        ibarrelall.Fill(eta, tmp.interaction);
        ractiveall.Fill(eta, tmp.radiation);
        iactiveall.Fill(eta, tmp.interaction);
        rglobal.Fill(eta, tmp.radiation);
        iglobal.Fill(eta, tmp.interaction);

        //      active volumes, endcap
        tmp = analyzeModules(mb.getEndcapModuleCaps(), eta, theta, phi, track, sumComponentsRI);
        ractiveendcap.Fill(eta, tmp.radiation);
        iactiveendcap.Fill(eta, tmp.interaction);
        rendcapall.Fill(eta, tmp.radiation);
        iendcapall.Fill(eta, tmp.interaction);
        ractiveall.Fill(eta, tmp.radiation);
        iactiveall.Fill(eta, tmp.interaction);
        rglobal.Fill(eta, tmp.radiation);
        iglobal.Fill(eta, tmp.interaction);

        for (std::map<std::string, Material>::iterator it = sumComponentsRI.begin(); it != sumComponentsRI.end(); ++it) {
          if (rComponents[it->first]==NULL) { 
            rComponents[it->first] = new TH1D();
            rComponents[it->first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
          }
          rComponents[it->first]->Fill(eta, it->second.radiation);
          if (iComponents[it->first]==NULL) {
            iComponents[it->first] = new TH1D();
            iComponents[it->first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
          }
          iComponents[it->first]->Fill(eta, it->second.interaction);
        }


        if (rComponents["Services"]==NULL) { 
          rComponents["Services"] = new TH1D();
          rComponents["Services"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
        }
        if (iComponents["Services"]==NULL) { 
          iComponents["Services"] = new TH1D();
          iComponents["Services"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
        }
        if (rComponents["Supports"]==NULL) { 
          rComponents["Supports"] = new TH1D();
          rComponents["Supports"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
        }
        if (iComponents["Supports"]==NULL) { 
          iComponents["Supports"] = new TH1D();
          iComponents["Supports"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
        }
        //      services, barrel
        tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getBarrelServices(), eta, theta, track, MaterialProperties::no_cat);
        rserfbarrel.Fill(eta, tmp.radiation);
        iserfbarrel.Fill(eta, tmp.interaction);
        rbarrelall.Fill(eta, tmp.radiation);
        ibarrelall.Fill(eta, tmp.interaction);
        rserfall.Fill(eta, tmp.radiation);
        iserfall.Fill(eta, tmp.interaction);
        rglobal.Fill(eta, tmp.radiation);
        iglobal.Fill(eta, tmp.interaction);
        rComponents["Services"]->Fill(eta, tmp.radiation);
        iComponents["Services"]->Fill(eta, tmp.interaction);
        //      services, endcap
        tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getEndcapServices(), eta, theta, track, MaterialProperties::no_cat);
        rserfendcap.Fill(eta, tmp.radiation);
        iserfendcap.Fill(eta, tmp.interaction);
        rendcapall.Fill(eta, tmp.radiation);
        iendcapall.Fill(eta, tmp.interaction);
        rserfall.Fill(eta, tmp.radiation);
        iserfall.Fill(eta, tmp.interaction);
        rglobal.Fill(eta, tmp.radiation);
        iglobal.Fill(eta, tmp.interaction);
        rComponents["Services"]->Fill(eta, tmp.radiation);
        iComponents["Services"]->Fill(eta, tmp.interaction);
        //      supports, barrel
        tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, MaterialProperties::b_sup);
        rlazybarrel.Fill(eta, tmp.radiation);
        ilazybarrel.Fill(eta, tmp.interaction);
        rbarrelall.Fill(eta, tmp.radiation);
        ibarrelall.Fill(eta, tmp.interaction);
        rlazyall.Fill(eta, tmp.radiation);
        ilazyall.Fill(eta, tmp.interaction);
        rglobal.Fill(eta, tmp.radiation);
        iglobal.Fill(eta, tmp.interaction);
        rComponents["Supports"]->Fill(eta, tmp.radiation);
        iComponents["Supports"]->Fill(eta, tmp.interaction);
        //      supports, endcap
        tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, MaterialProperties::e_sup);
        rlazyendcap.Fill(eta, tmp.radiation);
        ilazyendcap.Fill(eta, tmp.interaction);
        rendcapall.Fill(eta, tmp.radiation);
        iendcapall.Fill(eta, tmp.interaction);
        rlazyall.Fill(eta, tmp.radiation);
        ilazyall.Fill(eta, tmp.interaction);
        rglobal.Fill(eta, tmp.radiation);
        iglobal.Fill(eta, tmp.interaction);
        rComponents["Supports"]->Fill(eta, tmp.radiation);
        iComponents["Supports"]->Fill(eta, tmp.interaction);
        //      supports, tubes
        tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, MaterialProperties::o_sup);
        rlazytube.Fill(eta, tmp.radiation);
        ilazytube.Fill(eta, tmp.interaction);
        rlazyall.Fill(eta, tmp.radiation);
        ilazyall.Fill(eta, tmp.interaction);
        rglobal.Fill(eta, tmp.radiation);
        iglobal.Fill(eta, tmp.interaction);
        rComponents["Supports"]->Fill(eta, tmp.radiation);
        iComponents["Supports"]->Fill(eta, tmp.interaction);
        //      supports, barrel tubes
        tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, MaterialProperties::t_sup);
        rlazybtube.Fill(eta, tmp.radiation);
        ilazybtube.Fill(eta, tmp.interaction);
        rlazyall.Fill(eta, tmp.radiation);
        ilazyall.Fill(eta, tmp.interaction);
        rglobal.Fill(eta, tmp.radiation);
        iglobal.Fill(eta, tmp.interaction);
        rComponents["Supports"]->Fill(eta, tmp.radiation);
        iComponents["Supports"]->Fill(eta, tmp.interaction);
        //      supports, user defined
        tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, MaterialProperties::u_sup);
        rlazyuserdef.Fill(eta, tmp.radiation);
        ilazyuserdef.Fill(eta, tmp.interaction);
        rlazyall.Fill(eta, tmp.radiation);
        ilazyall.Fill(eta, tmp.interaction);
        rglobal.Fill(eta, tmp.radiation);
        iglobal.Fill(eta, tmp.interaction);
        rComponents["Supports"]->Fill(eta, tmp.radiation);
        iComponents["Supports"]->Fill(eta, tmp.interaction);
        //      pixels, if they exist
        if (pm != NULL) {
          std::map<std::string, Material> ignoredPixelSumComponentsRI;
          analyzeModules(pm->getBarrelModuleCaps(), eta, theta, phi, track, ignoredPixelSumComponentsRI, true);
          analyzeModules(pm->getEndcapModuleCaps(), eta, theta, phi, track, ignoredPixelSumComponentsRI, true);
          analyzeInactiveSurfaces(pm->getInactiveSurfaces().getBarrelServices(), eta, theta, track, MaterialProperties::no_cat, true);
          analyzeInactiveSurfaces(pm->getInactiveSurfaces().getEndcapServices(), eta, theta, track, MaterialProperties::no_cat, true);
          analyzeInactiveSurfaces(pm->getInactiveSurfaces().getSupports(), eta, theta, track, MaterialProperties::no_cat, true);
        }

        // Add the hit on the beam pipe
        Hit* hit = new Hit(23./sin(theta));
        hit->setOrientation(Hit::Horizontal);
        hit->setObjectKind(Hit::Inactive);
        Material beamPipeMat;
        beamPipeMat.radiation = 0.0023 / sin(theta);
        beamPipeMat.interaction = 0.0019 / sin(theta);
        hit->setCorrectedMaterial(beamPipeMat);
        track.addHit(hit);
        if (!track.noHits()) {
          track.sort();
          if (efficiency!=1) track.addEfficiency(efficiency, false);
          if (pixelEfficiency!=1) track.addEfficiency(efficiency, true);
          if (track.nActiveHits(true)>2) { // At least 3 points needed
            //#ifdef HIT_DEBUG_RZ
            //      if ((eta >1.617) &&(theta>0.383)) {
            //        Track::debugRZCorrelationMatrix = true;
            //        Track::debugRZCovarianceMatrix = true;
            //        Track::debugRZErrorPropagation = true;
            //      }
            //#endif
            if (computeResolution) track.computeErrors(momenta);
            tv.push_back(track);

            if (computeResolution) {
              Track trackIdeal = track;
              trackIdeal.removeMaterial();
              trackIdeal.computeErrors(momenta);
              tvIdeal.push_back(trackIdeal);
            }
          }

          // @@ Hadrons
          int nActive = track.nActiveHits();
          if (nActive>0) {
            hadronTotalHitsGraph.SetPoint(hadronTotalHitsGraph.GetN(),
                                          eta,
                                          nActive);
            double probability;
            std::vector<double> probabilities = track.hadronActiveHitsProbability();

            double averageHits=0;
            //double averageSquaredHits=0;
            double exactProb=0;
            double moreThanProb = 0;
            for (int i=probabilities.size()-1;
                 i>=0;
                 --i) {
              //if (nActive==10) { // debug
              //  std::cerr << "probabilities.at(" 
              //  << i << ")=" << probabilities.at(i)
              //  << endl;
              //}
              exactProb=probabilities.at(i)-moreThanProb;
              averageHits+=(i+1)*exactProb;
              //averageSquaredHits+=((i+1)*(i+1))*exactProb;
              moreThanProb+=exactProb;
            }
            hadronAverageHitsGraph.SetPoint(hadronAverageHitsGraph.GetN(),
                                            eta,
                                            averageHits);
            //hadronAverageHitsGraph.SetPointError(hadronAverageHitsGraph.GetN()-1,
            //                       0,
            //                       sqrt( averageSquaredHits - averageHits*averageHits) );

            unsigned int requiredHits;
            for (unsigned int i = 0;
                 i<hadronNeededHitsFraction.size();
                 ++i) {
              requiredHits = int(ceil(double(nActive) * hadronNeededHitsFraction.at(i)));
              if (requiredHits==0)
                probability=1;
              else if (requiredHits>probabilities.size())
                probability = 0;
              else
                probability = probabilities.at(requiredHits-1);
              //if (probabilities.size()==10) { // debug
              //  std::cerr << "required " << requiredHits
              //              << " out of " << probabilities.size()
              //              << " == " << nActive
              //              << endl;
              // std::cerr << "      PROBABILITY = " << probability << endl << endl;
              //}
              hadronGoodTracksFraction.at(i).SetPoint(hadronGoodTracksFraction.at(i).GetN(),
                                                      eta,
                                                      probability);
            }
          }
        }
      }

#ifdef MATERIAL_SHADOW       
      // integration over eta
      for (unsigned int i = 0; i < cells.size(); i++) {
        for (unsigned int j = 1; j < cells.at(i).size(); j++) {
          cells.at(i).at(j).rlength = cells.at(i).at(j).rlength + cells.at(i).at(j - 1).rlength;
          cells.at(i).at(j).ilength = cells.at(i).at(j).ilength + cells.at(i).at(j - 1).ilength;
        }
      }
      // transformation from (eta, r) to (z, r) coordinates
      transformEtaToZ();
#endif // MATERIAL_SHADOW

      if (computeResolution) {
        // fill TGraph map
        calculateGraphs(momenta, tv, graphBag::StandardGraph | graphBag::RealGraph);
        calculateGraphs(momenta, tvIdeal, graphBag::StandardGraph | graphBag::IdealGraph);
      }
    }

    void Analyzer::analyzePower(Tracker& tracker) {
      computeIrradiatedPowerConsumption(tracker);
      preparePowerHistograms();
      fillPowerMap(tracker);
    }

    void Analyzer::computeTriggerFrequency(Tracker& tracker) {
      std::map<std::string, std::map<std::pair<int,int>, int> >   triggerFrequencyCounts;
      std::map<std::string, std::map<std::pair<int,int>, double> > triggerFrequencyAverageTrue,
        triggerFrequencyInterestingParticleTrue,
        triggerFrequencyAverageFake; // trigger frequency by module in Z and R, averaged over Phi
      LayerVector& layers = tracker.getLayers();

      triggerFrequencyTrueSummaries_.clear();
      triggerFrequencyFakeSummaries_.clear();
      triggerRateSummaries_.clear();
      triggerPuritySummaries_.clear();
      triggerDataBandwidthSummaries_.clear();
      triggerDataBandwidths_.clear();
      stripOccupancySummaries_.clear();
      hitOccupancySummaries_.clear();

      double nMB = tracker.getNMB();

      for(LayerVector::iterator layIt = layers.begin(); layIt != layers.end(); ++layIt) {
        Layer* layer = *layIt;
        ModuleVector* modules = layer->getModuleVector();
        if (!modules) {
          logERROR("Cannot retrieve moduleVector from the layer");
          return;
        }
        std::string cntName = layer->getContainerName();
        //  int layIndex = (*layIt)->getIndex();
        if (cntName == "") {
          ostringstream tempSS;
          tempSS << "Skipping layer with no container name (" << modules->size() << " modules)";
          logWARNING(tempSS);
          continue;
        }

       
        triggerFrequencyTrueSummaries_[cntName].setHeader("Layer", "Ring");
        triggerFrequencyFakeSummaries_[cntName].setHeader("Layer", "Ring");
        triggerRateSummaries_[cntName].setHeader("Layer", "Ring");
        triggerPuritySummaries_[cntName].setHeader("Layer", "Ring");
        triggerDataBandwidthSummaries_[cntName].setHeader("Layer", "Ring");
        stripOccupancySummaries_[cntName].setHeader("Layer", "Ring");
        hitOccupancySummaries_[cntName].setHeader("Layer", "Ring");
        triggerFrequencyTrueSummaries_[cntName].setPrecision(3);
        triggerFrequencyFakeSummaries_[cntName].setPrecision(3);
        triggerRateSummaries_[cntName].setPrecision(3);
        triggerPuritySummaries_[cntName].setPrecision(3);
        triggerDataBandwidthSummaries_[cntName].setPrecision(3);
        stripOccupancySummaries_[cntName].setPrecision(3);
        hitOccupancySummaries_[cntName].setPrecision(3);


        //int nbins = (dynamic_cast<BarrelLayer*>(layer)) ? ((BarrelLayer*)layer)->getModulesOnRod() : ((EndcapLayer*)layer)->getRings();
        int nbins = layer->getRings();

        for (ModuleVector::iterator modIt = modules->begin(); modIt != modules->end(); ++modIt) {
          Module* module = (*modIt);
          XYZVector center = module->getMeanPoint();
          if ((center.Z()<0) || module->getPhiIndex() > 1/*(center.Phi()<0) || (center.Phi()>M_PI/2)*/ || (module->getStereoDistance()==0.0)) continue;

          TH1D* currentTotalHisto;
          TH1D* currentTrueHisto;

          if (totalStubRateHistos_.count(std::make_pair(cntName, layer->getIndex())) == 0) {
            currentTotalHisto = new TH1D(("totalStubsPerEventHisto" + cntName + any2str(layer->getIndex())).c_str(), ";Modules;MHz/cm^2", nbins, 0.5, nbins+0.5);
            currentTrueHisto = new TH1D(("trueStubsPerEventHisto" + cntName + any2str(layer->getIndex())).c_str(), ";Modules;MHz/cm^2", nbins, 0.5, nbins+0.5); 
            totalStubRateHistos_[std::make_pair(cntName, layer->getIndex())] = currentTotalHisto; 
            trueStubRateHistos_[std::make_pair(cntName, layer->getIndex())] = currentTrueHisto; 
          } else {
            currentTotalHisto = totalStubRateHistos_[std::make_pair(cntName, layer->getIndex())]; 
            currentTrueHisto = trueStubRateHistos_[std::make_pair(cntName, layer->getIndex())]; 
          }


          int curCnt = triggerFrequencyCounts[cntName][make_pair(module->getLayer(), module->getRing())]++;
          double curAvgTrue = triggerFrequencyAverageTrue[cntName][make_pair(module->getLayer(),module->getRing())];
          double curAvgInteresting = triggerFrequencyInterestingParticleTrue[cntName][make_pair(module->getLayer(),module->getRing())];
          double curAvgFake = triggerFrequencyAverageFake[cntName][make_pair(module->getLayer(),module->getRing())];

          //curAvgTrue  = curAvgTrue + (module->getTriggerFrequencyTruePerEvent()*tracker.getNMB() - curAvgTrue)/(curCnt+1);
          //curAvgFake  = curAvgFake + (module->getTriggerFrequencyFakePerEvent()*pow(tracker.getNMB(),2) - curAvgFake)/(curCnt+1); // triggerFrequencyFake scales with the square of Nmb!

          // TODO! Important <- make this interestingPt cut configurable
          const double interestingPt = 2;
          curAvgTrue  = curAvgTrue + (module->getTriggerFrequencyTruePerEventAbove(interestingPt)*tracker.getNMB() - curAvgTrue)/(curCnt+1);
          curAvgInteresting += (module->getParticleFrequencyPerEventAbove(interestingPt)*tracker.getNMB() - curAvgInteresting)/(curCnt+1);
          curAvgFake  = curAvgFake + (
                                      (module->getTriggerFrequencyFakePerEvent()*tracker.getNMB()+module->getTriggerFrequencyTruePerEventBelow(interestingPt))
                                      *tracker.getNMB() - curAvgFake)/(curCnt+1); // triggerFrequencyFake scales with the square of Nmb!

          double curAvgTotal = curAvgTrue + curAvgFake;

          triggerFrequencyAverageTrue[cntName][make_pair(module->getLayer(), module->getRing())] = curAvgTrue;            
          triggerFrequencyInterestingParticleTrue[cntName][make_pair(module->getLayer(), module->getRing())] = curAvgInteresting;    
          triggerFrequencyAverageFake[cntName][make_pair(module->getLayer(), module->getRing())] = curAvgFake;    

          int triggerDataHeaderBits  = tracker.getModuleType(module->getType()).getTriggerDataHeaderBits();
          int triggerDataPayloadBits = tracker.getModuleType(module->getType()).getTriggerDataPayloadBits();
          double triggerDataBandwidth = (triggerDataHeaderBits + curAvgTotal*triggerDataPayloadBits) / (tracker.getBunchSpacingNs()); // GIGABIT/second
          triggerDataBandwidths_[cntName][make_pair(module->getLayer(), module->getRing())] = triggerDataBandwidth;
          triggerFrequenciesPerEvent_[cntName][make_pair(module->getLayer(), module->getRing())] = curAvgTotal;

          module->setProperty("triggerDataBandwidth", triggerDataBandwidth); // averaged over phi
          module->setProperty("triggerFrequencyPerEvent", curAvgTotal); // averaged over phi


          //                currentTotalGraph->SetPoint(module->getRing()-1, module->getRing(), curAvgTotal*(1000/tracker.getBunchSpacingNs())*(100/module->getArea()));
          //                currentTrueGraph->SetPoint(module->getRing()-1, module->getRing(), curAvgTrue*(1000/tracker.getBunchSpacingNs())*(100/module->getArea()));

          currentTotalHisto->SetBinContent(module->getRing(), curAvgTotal*(1000/tracker.getBunchSpacingNs())*(100/module->getArea()));
          currentTrueHisto->SetBinContent(module->getRing(), curAvgTrue*(1000/tracker.getBunchSpacingNs())*(100/module->getArea()));

          triggerFrequencyTrueSummaries_[cntName].setCell(module->getLayer(), module->getRing(), curAvgTrue);
          triggerFrequencyInterestingSummaries_[cntName].setCell(module->getLayer(), module->getRing(), curAvgInteresting);
          triggerFrequencyFakeSummaries_[cntName].setCell(module->getLayer(), module->getRing(), curAvgFake);
          triggerRateSummaries_[cntName].setCell(module->getLayer(), module->getRing(), curAvgTotal);             
          triggerEfficiencySummaries_[cntName].setCell(module->getLayer(), module->getRing(), curAvgTrue/curAvgInteresting);                
          triggerPuritySummaries_[cntName].setCell(module->getLayer(), module->getRing(), curAvgTrue/(curAvgTrue+curAvgFake));                
          triggerDataBandwidthSummaries_[cntName].setCell(module->getLayer(), module->getRing(), triggerDataBandwidth);

          stripOccupancySummaries_[cntName].setCell(module->getLayer(), module->getRing(), module->getStripOccupancyPerEvent()*nMB*100); 
          hitOccupancySummaries_[cntName].setCell(module->getLayer(), module->getRing(), module->getHitOccupancyPerEvent()*nMB*100); 

        }
      }

    }


    std::pair<Circle, Circle> findCirclesTwoPoints(const Point& p1, const Point& p2, double r) {
      double x1 = p1.x, y1 = p1.y, x2 = p2.x, y2 = p2.y;

      double q = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
      double x3 = (x1+x2)/2;
      double y3 = (y1+y2)/2;

      double xc1 = x3 + sqrt(r*r - pow(q/2,2))*(y1-y2)/q;
      double yc1 = y3 + sqrt(r*r - pow(q/2,2))*(x2-x1)/q;

      double xc2 = x3 - sqrt(r*r - pow(q/2,2))*(y1-y2)/q;
      double yc2 = y3 - sqrt(r*r - pow(q/2,2))*(x2-x1)/q;
      
      return std::make_pair((Circle){ xc1, yc1, r }, (Circle){ xc2, yc2, r });
    }

    bool isPointInCircle(const Point& p, const Circle& c) { return pow(p.x - c.x0, 2) + pow(p.y - c.y0, 2) <= c.r*c.r; }



    double Analyzer::calculatePetalAreaMC(const Tracker& tracker, double crossoverR) const {
      static TRandom3 die;

      double r = tracker.getTriggerPtCut()/(0.3*insur::magnetic_field) * 1e3; // curvature radius of particles with the minimum accepted pt

      std::pair<Circle, Circle> cc = findCirclesTwoPoints((Point){0, 0}, (Point){0, crossoverR}, r);

      // Monte Carlo area calculation
      int hits = 0;
      double maxR = tracker.getMaxR(); // points randomly generated in a 40 degrees circle slice
      double minR = tracker.getMinR();
      double aperture = 0.34906585 * 2; // 40 degrees
      //double maxPhi = M_PI/2 + aperture/2;
      double minPhi = M_PI/2 - aperture/2;
      for (int i = 0; i < 100000; i++) { 
        double rr  = minR + die.Rndm()*(maxR-minR);
        double phi = minPhi + die.Rndm()*aperture;
        Polar2DPoint rndp(rr, phi); 
        bool inFirstCircle  = isPointInCircle((Point){rndp.X(), rndp.Y()}, cc.first);
        bool inSecondCircle = isPointInCircle((Point){rndp.X(), rndp.Y()}, cc.second); 
        if ((inFirstCircle && inSecondCircle) || (!inFirstCircle && !inSecondCircle)) hits++; // if it's in both circles means it's in the lower part of the petal (before the crossover), if it's outside both it means it's upper part of the petal (after the crossover)
      }

      return hits;
    }

    double Analyzer::calculatePetalAreaModules(Tracker& tracker, double crossoverR) const {
      double curvatureR = tracker.getParticleCurvatureR(tracker.getTriggerPtCut()); // curvature radius of particles with the minimum accepted pt
      int hits = 0;
      LayerVector* layers = tracker.getEndcapLayers();
      for(LayerVector::iterator layIt = layers->begin(); layIt != layers->end(); ++layIt) { // loop over layers
        ModuleVector* modules = (*layIt)->getModuleVector();
        std::string cntName = (*layIt)->getContainerName();
        Layer* layer = (*layIt);
        if (cntName == "" || !layer) continue;

        if ((*modules->begin())->getDisk() != 1 || (*modules->begin())->getZSide() < 0) continue;

        for (ModuleVector::iterator modIt = modules->begin(); modIt != modules->end(); ++modIt) { // loop over modules
          Module* module = (*modIt);

          if (isModuleInPetal(module, 0., curvatureR, crossoverR)) hits++;
        }
      }
      return hits;
    }

/*    double Analyzer::calculatePetalCrossoverRecursive(Tracker& tracker, double maxR, double minR, int recIndex) {
      double cr1 = (maxR + (maxR+minR)/2.)/2.;
      double cr2 = ((maxR+minR)/2. + minR)/2.;
      double area1 = calculatePetalArea(tracker, cr1);
      double area2 = calculatePetalArea(tracker, cr2);

      if (++recIndex == 100) { logERROR("Method didn't converge. Last calculated value: " + any2str(cr1)); return cr1; }
      double relDiff = fabs(area1 - area2)/area2; 
      if (relDiff < 1e-4) { logDEBUG("Method converged after " + any2str(recIndex) + " iterations. Value: " + any2str(cr1)); return cr1; }
      if (area1 < area2) { return calculatePetalCrossoverRecursive(tracker, maxR, (maxR+minR)/2., recIndex); }
      else { return calculatePetalCrossoverRecursive(tracker, (maxR+minR)/2., minR, recIndex); }
    }
*/
    double Analyzer::calculatePetalCrossover(Tracker& tracker) {
      class Trampoline : public ROOT::Math::IBaseFunctionOneDim {
        Analyzer& a_;
        Tracker& t_;
        double DoEval(double x) const { return a_.calculatePetalAreaModules(t_, x); }
      public:
        Trampoline(Analyzer& a, Tracker& t) : a_(a), t_(t) {}
        ROOT::Math::IBaseFunctionOneDim* Clone() const { return new Trampoline(a_, t_); }
      };

      Trampoline t(*this, tracker);

      ROOT::Math::BrentMinimizer1D minBrent;
      minBrent.SetFunction(t, 0., tracker.getMaxR());
      bool ok = minBrent.Minimize(100, 0.001, 0.001);


      if (ok) {
        logINFO("Multiple Trigger Towers: Searching for optimized petal crossover point");
        logINFO("  Method converged after " + any2str(minBrent.Iterations()) + " iterations.");
        logINFO("  Found minimum: crossover point = " + any2str(minBrent.XMinimum()) + "  petal area = " + any2str(minBrent.FValMinimum()));
      } else {
        logERROR("Multiple Trigger Towers: Search for optimized petal crossover point failed");
        logERROR("  Method did not converge after " + any2str(minBrent.Iterations()) + " iterations.");
        logERROR("  Last found value: crossover point = " + any2str(minBrent.XMinimum()) + "  petal area = " + any2str(minBrent.FValMinimum()));
      }

      return minBrent.XMinimum();
    }





    bool Analyzer::isModuleInEtaSector(const Tracker& tracker, const Module* module, int etaSector) const {
      int numProcEta = tracker.getTriggerProcessorsEta();
      double etaCut = tracker.getTriggerEtaCut();
      double etaSlice = etaCut*2 / numProcEta;
      double maxR = tracker.getMaxR();
      double zErrorCollider = tracker.getZErrorCollider();
      double eta = etaSlice*etaSector-etaCut;    

      double modMinZ = module->getMinZ();
      double modMaxZ = module->getMaxZ();
      double modMinR = module->getMinRho();                
      double modMaxR = module->getMaxRho();                

      double etaSliceZ1 = maxR/tan(2*atan(exp(-eta)));
      double etaSliceZ2 = maxR/tan(2*atan(exp(-eta-etaSlice)));

      double etaDist1 =  modMaxZ - ((etaSliceZ1 >= 0 ? modMinR : modMaxR)*(etaSliceZ1 + zErrorCollider)/maxR - zErrorCollider); // if etaDists are positive it means the module is in the slice
      double etaDist2 = -modMinZ + ((etaSliceZ2 >= 0 ? modMaxR : modMinR)*(etaSliceZ2 - zErrorCollider)/maxR + zErrorCollider); 

      return etaDist1 > 0 && etaDist2 > 0;
    }

    bool Analyzer::isModuleInPetal(const Module* module, double petalPhi, double curvatureR, double crossoverR) const {
      Polar2DPoint crossoverPoint(crossoverR, petalPhi);
      double proj = cos(module->getMeanPoint().Phi() - petalPhi); // check if module is in the same semi-plane as the petal by projecting its center on the petal symmetry line
      if (proj < 0.) return false;
      std::pair<Circle, Circle> cc = findCirclesTwoPoints((Point){0.,0.}, (Point){crossoverPoint.X(), crossoverPoint.Y()}, curvatureR);
      
      int inFirstCircle = 0, inSecondCircle = 0;
      for (int i = 0; i < 4; i++) {
        const XYZVector& corner = module->getCorner(i);
        inFirstCircle  |= (isPointInCircle((Point){corner.X(), corner.Y()}, cc.first) << i);
        inSecondCircle |= (isPointInCircle((Point){corner.X(), corner.Y()}, cc.second) << i);
      }
      return (inFirstCircle && inSecondCircle) || (inFirstCircle < 0xF && inSecondCircle < 0xF);
    }


    bool areClockwise(const Point& p1, const Point& p2) { return -p1.x*p2.y + p1.y*p2.x > 0; }
//#define OLD_PHI_SECTOR_CHECK

#ifndef OLD_PHI_SECTOR_CHECK
    bool Analyzer::isModuleInCircleSector(const Module* module, double startPhi, double endPhi) const {

      Point startArm = {cos(startPhi), sin(startPhi)};
      Point endArm   = {cos(endPhi), sin(endPhi)};

      for (int i = 0; i < 4; i++) {
        const XYZVector& corner = module->getCorner(i);
        if (!areClockwise(startArm, (Point){corner.X(), corner.Y()}) && areClockwise(endArm, (Point){corner.X(), corner.Y()})) return true;
      }
      return false;
    }
#else
    bool Analyzer::isModuleInCircleSector(const Module* module, double sliceMinPhi, double sliceMaxPhi) const {

      double modMinPhi = module->getMinPhi() >= 0 ? module->getMinPhi() : module->getMinPhi() + 2*M_PI;
      double modMaxPhi = module->getMaxPhi() >= 0 ? module->getMaxPhi() : module->getMaxPhi() + 2*M_PI;

      //double trajSlice = asin((modMaxR+modMinR)/2 * 0.0003 * magnetic_field / (2 * tracker.getTriggerPtCut())); // aka Alpha

      if (modMinPhi > modMaxPhi && sliceMaxPhi > 2*M_PI) modMaxPhi += 2*M_PI;      // this solves the issue with modules across the 2 PI line
      else if (modMinPhi > modMaxPhi && sliceMaxPhi < 2*M_PI) modMinPhi -= 2*M_PI; // 

      bool inSectorSlice = ((sliceMinPhi < modMaxPhi && modMinPhi < sliceMaxPhi) ||
                            (sliceMinPhi < modMaxPhi+2*M_PI && modMinPhi+2*M_PI < sliceMaxPhi) || // this catches the modules that are at a small angle but must be caught by a sweep crossing the 2 PI line
                            (sliceMinPhi < modMaxPhi-2*M_PI && modMinPhi-2*M_PI < sliceMaxPhi)); 

      return inSectorSlice;

    }
#endif

    bool Analyzer::isModuleInPhiSector(const Tracker& tracker, const Module* module, double crossoverR, int phiSector) const {
      static const double curvatureR = tracker.getTriggerPtCut()/(0.3*insur::magnetic_field) * 1e3; // curvature radius of particles with the minimum accepted pt

      double phiSlice = 2*M_PI / tracker.getTriggerProcessorsPhi();  // aka Psi
      double phi = phiSlice*phiSector;

      double sliceMinPhi = phi;
      double sliceMaxPhi = phi + phiSlice;

      bool inSectorSlice = isModuleInCircleSector(module, sliceMinPhi, sliceMaxPhi);
      bool inPetals = isModuleInPetal(module, sliceMinPhi, curvatureR, crossoverR) || isModuleInPetal(module, sliceMaxPhi, curvatureR, crossoverR);

      return inSectorSlice || inPetals;
    }


    int sebifyBarrelCoords(BarrelModule& m, BarrelLayer& l, Tracker& t) { // transform tkLayout module coordinates into the format used by Sebastian Viret (a.k.a. Sebification)
      int layoffset = 0;
      for (int i = 1; i < l.getContainerId(); i++) layoffset += t.getNumLayersInContainer(i);
      int layer = 4 + m.getLayer() + layoffset; 
      int ladder = m.getPhiIndex();
      int module = (m.getRing()*m.getZSide() + l.getModulesOnRodSide(-1) - (m.getZSide()>0));
      return layer*10000 + ladder*100 + module;
    }

    int sebifyEndcapCoords(EndcapModule& m, EndcapLayer&, Tracker&) {
      int disk = 10 + m.getDisk() + 7 * (m.getZSide()<0);
      int ring = m.getRing()-1;
      int module = m.getPhiIndex();
      return disk*10000 + ring*100 + module;
    }



    void Analyzer::computeTriggerProcessorsBandwidth(Tracker& tracker) {

      std::map<std::pair<int, int>, int> processorConnections;
      std::map<std::pair<int, int>, double> processorInboundBandwidths;
      std::map<std::pair<int, int>, double> processorInboundStubsPerEvent;

      processorConnectionSummary_.setHeader("Phi", "Eta");
      processorInboundBandwidthSummary_.setHeader("Phi", "Eta");
      processorInboundStubPerEventSummary_.setHeader("Phi", "Eta");

      processorCommonConnectionSummary_.setHeader("tEta,Phi", "tEta,Phi");

      processorInboundBandwidthSummary_.setPrecision(3);
      processorInboundStubPerEventSummary_.setPrecision(3);

      int numProcEta = tracker.getTriggerProcessorsEta();
      int numProcPhi = tracker.getTriggerProcessorsPhi();

      int totalProcs = numProcEta * numProcPhi;
      processorCommonConnectionMap_.SetBins(totalProcs, 0, totalProcs, totalProcs, 0, totalProcs);
      processorCommonConnectionMap_.SetXTitle("TT");
      processorCommonConnectionMap_.SetYTitle("TT");
      

      double crossoverR = calculatePetalCrossover(tracker);

      triggerPetalCrossoverR_ = crossoverR;
      sampleTriggerPetal_ = findCirclesTwoPoints((Point){0., 0.}, (Point){crossoverR, 0.}, tracker.getParticleCurvatureR(tracker.getTriggerPtCut()));

      LayerVector& layers = tracker.getLayers();
      for(LayerVector::iterator layIt = layers.begin(); layIt != layers.end(); ++layIt) { // loop over layers
        ModuleVector* modules = (*layIt)->getModuleVector();
        std::string cntName = (*layIt)->getContainerName();
        Layer* layer = (*layIt);
        if (cntName == "" || !layer) continue;

        for (ModuleVector::iterator modIt = modules->begin(); modIt != modules->end(); ++modIt) { // loop over modules
          Module* module = (*modIt);
          
          int etaConnections = 0, totalConnections = 0;

          for (int i=0; i < numProcEta; i++) {
            if (isModuleInEtaSector(tracker, module, i)) {
              etaConnections++;
              for (int j=0; j < numProcPhi; j++) {
                if (isModuleInPhiSector(tracker, module, crossoverR, j)) {
                  totalConnections++;

                  processorConnections[std::make_pair(j,i)] += 1;
                  processorConnectionSummary_.setCell(j+1, i+1, processorConnections[std::make_pair(j,i)]);

                  processorInboundBandwidths[std::make_pair(j,i)] += triggerDataBandwidths_[cntName][make_pair(module->getLayer(), module->getRing())]; // *2 takes into account negative Z's
                  processorInboundBandwidthSummary_.setCell(j+1, i+1, processorInboundBandwidths[std::make_pair(j,i)]);

                  processorInboundStubsPerEvent[std::make_pair(j,i)] += triggerFrequenciesPerEvent_[cntName][make_pair(module->getLayer(), module->getRing())];
                  processorInboundStubPerEventSummary_.setCell(j+1, i+1, processorInboundStubsPerEvent[std::make_pair(j,i)]);
                  
                  module->addConnectedProcessor(make_pair(i+1, j+1));

                  if (dynamic_cast<BarrelModule*>(module)) 
                    sectorMap_[make_pair(i+1, j+1)].insert(sebifyBarrelCoords(*(BarrelModule*)module, *(BarrelLayer*)layer, tracker));
                  else 
                    sectorMap_[make_pair(i+1, j+1)].insert(sebifyEndcapCoords(*(EndcapModule*)module, *(EndcapLayer*)layer, tracker));
                } 
              }
            }
          }
          module->setProcessorConnectionsEta(etaConnections);
          module->setProcessorConnectionsPhi(totalConnections > 0 ? totalConnections/etaConnections : 0);
        }
      }


      double inboundBandwidthTotal = 0.;
      int processorConnectionsTotal = 0;
      double inboundStubsPerEventTotal = 0.;
      for (std::map<std::pair<int, int>, double>::iterator it = processorInboundBandwidths.begin(); it != processorInboundBandwidths.end(); ++it)
        inboundBandwidthTotal += (*it).second;
      for (std::map<std::pair<int, int>, int>::iterator it = processorConnections.begin(); it != processorConnections.end(); ++it)
        processorConnectionsTotal += (*it).second;
      for (std::map<std::pair<int, int>, double>::iterator it = processorInboundStubsPerEvent.begin(); it != processorInboundStubsPerEvent.end(); ++it) 
        inboundStubsPerEventTotal += (*it).second;


      processorInboundBandwidthSummary_.setSummaryCell("Total", inboundBandwidthTotal);
      processorConnectionSummary_.setSummaryCell("Total", processorConnectionsTotal);
      processorInboundStubPerEventSummary_.setSummaryCell("Total", inboundStubsPerEventTotal);

      moduleConnectionsDistribution.Reset();
      moduleConnectionsDistribution.SetNameTitle("ModuleConnDist", "Number of connections to trigger processors;Connections;Modules");
      moduleConnectionsDistribution.SetBins(11, -.5, 10.5);

      std::map<std::pair<int, int>, int> processorCommonConnectionMatrix;

      for(LayerVector::iterator layIt = layers.begin(); layIt != layers.end(); ++layIt) { // loop over layers
        ModuleVector* modules = (*layIt)->getModuleVector();
        std::string cntName = (*layIt)->getContainerName();
        //BarrelLayer* layer = dynamic_cast<BarrelLayer*>(*layIt);
        if (cntName == "" /*|| !layer*/) continue;
        for (ModuleVector::iterator modIt = modules->begin(); modIt != modules->end(); ++modIt) { // loop over modules
          moduleConnectionsDistribution.Fill((*modIt)->getProcessorConnections(), 1);
          std::set<pair<int, int> > connectedProcessors = (*modIt)->getConnectedProcessors();
          if (connectedProcessors.size() == 1) {
            int ref = connectedProcessors.begin()->second + numProcPhi*(connectedProcessors.begin()->first-1);
            processorCommonConnectionMatrix[std::make_pair(ref, ref)] += 1;
          } else {
            while (!connectedProcessors.empty()) {
              pair<int, int> colRef = *connectedProcessors.begin();
              int col = colRef.second + numProcPhi*(colRef.first-1);
              connectedProcessors.erase(connectedProcessors.begin());
              for (std::set<pair<int, int> >::const_iterator pIt = connectedProcessors.begin(); pIt != connectedProcessors.end(); ++pIt) {
                int row = pIt->second + numProcPhi*(pIt->first-1);
                processorCommonConnectionMatrix[std::make_pair(row, col)] += 1;
              }
            }
          }
        }
      }

      //for (std::map<std::pair<int, int>, int>::const_iterator mit = processorCommonConnectionMatrix.begin(); mit != processorCommonConnectionMatrix.end(); ++mit) {
      //  std::cout << mit->first.first << ',' << mit->first.second << '=' << mit->second << std::endl;
      //}
      TAxis* xAxis = processorCommonConnectionMap_.GetXaxis();
      TAxis* yAxis = processorCommonConnectionMap_.GetYaxis();
      for (int i = 1; i <= numProcEta; i++) {
        for (int j = 1; j <= numProcPhi; j++) {
          processorCommonConnectionSummary_.setCell(0, j + (i-1)*numProcPhi, "t" + any2str(i) + "," + any2str(j));
          processorCommonConnectionSummary_.setCell(j + (i-1)*numProcPhi, 0, "t" + any2str(i) + "," + any2str(j));
          xAxis->SetBinLabel(j + (i-1)*numProcPhi, ("t" + any2str(i) + "," + any2str(j)).c_str());
          yAxis->SetBinLabel(j + (i-1)*numProcPhi, ("t" + any2str(i) + "," + any2str(j)).c_str());
        }
      }

      for (int col = 1; col <= numProcEta*numProcPhi; col++) {
        for (int row = col; row <= numProcEta*numProcPhi; row++) {
          if (processorCommonConnectionMatrix.count(std::make_pair(row, col))) {
            int val = processorCommonConnectionMatrix[std::make_pair(row, col)];
            processorCommonConnectionSummary_.setCell(row, col, val);
            processorCommonConnectionMap_.SetCellContent(row, col, val); 
          }
          //else processorCommonConnectionSummary_.setCell(row, col, "0");
        }
      }

    }
    // function that finds in which quadrant of a Fluka grid cell the module center is placed 
	int Analyzer::whichCellQuadrant(double z, double r){
		double zc = floor(z);
		double rc = floor(r);
		int quadrant(-1);
		if( z > (zc+0.5) && r > (rc+0.5) ){quadrant = 1;}
		if( z < (zc+0.5) && r > (rc+0.5) ){quadrant = 2;}
		if( z < (zc+0.5) && r < (rc+0.5) ){quadrant = 3;}
		if( z > (zc+0.5) && r < (rc+0.5) ){quadrant = 4;}
		return quadrant;
	}

	void Analyzer::computeIrradiatedPowerConsumption(Tracker& tracker) {

		//double numInvFemtobarns = tracker.getNumInvFemtobarns();
		double operatingTemp    = tracker.getOperatingTemp();
		double chargeDepletionVoltage    = tracker.getChargeDepletionVoltage();
		double alphaParam       = tracker.getAlphaParam();
		double referenceTemp    = tracker.getReferenceTemp();

		double irrStepZ  = 10.0*tracker.getIrradiationInfoZ().stepC;//in mm!
		double irrStepR  = 10.0*tracker.getIrradiationInfoR().stepC;//in mm!
		double irrMinZ = 10.0*tracker.getIrradiationInfoZ().minC; //in mm! 
		double irrMinR = 10.0*tracker.getIrradiationInfoR().minC;// in mm!


		// cout << "numInvFemtobarns = " << tracker.getNumInvFemtobarns() << endl;
		// cout << "operatingTemp    = " << tracker.getOperatingTemp() << endl;
		// cout << "chargeDepletionVoltage    = " << tracker.getChargeDepletionVoltage() << endl;
		// cout << "alphaParam       = " << tracker.getAlphaParam() << endl;
		// cout << "referenceTemp    = " << tracker.getReferenceTemp() << endl;
		irradiatedPowerConsumptionSummaries_.clear();   
		maxIrradiatedPowerConsumptionSummaries_.clear();   
		maxIrradiatedPowerConsumptionTableSummaries_.clear();   

		LayerVector& layers = tracker.getLayers();
		for(LayerVector::iterator layIt = layers.begin(); layIt != layers.end(); ++layIt) {
			ModuleVector* modules = (*layIt)->getModuleVector();
			if (!modules) {
				logERROR("cannot retrieve moduleVector from the layer");
				return;
			}
			std::string cntName = (*layIt)->getContainerName();
			if (cntName == "") {
				ostringstream tempSS;
				tempSS << "Skipping layer with no container name (" << modules->size() << " modules)";
				logWARNING(tempSS);
				continue;
			}

			irradiatedPowerConsumptionSummaries_[cntName].setHeader("layer", "ring");
			irradiatedPowerConsumptionSummaries_[cntName].setPrecision(3);        

			for (ModuleVector::iterator modIt = modules->begin(); modIt != modules->end(); ++modIt) {
				Module* module = (*modIt); 
				XYZVector center = module->getMeanPoint();
				// if (center.Z()<0) continue; // I want to assign the right value to all modules, for the totals

				std::string moduleName = module->getType();
				std::string stDistance = any2str<double>(module->getStereoDistance());
				std::string moduleDistance = moduleName + "_" + stDistance;
				maxIrradiatedPowerConsumptionTableSummaries_[moduleDistance].setHeader("layer", "ring");
				maxIrradiatedPowerConsumptionTableSummaries_[moduleDistance].setPrecision(3);        

				double volume  = tracker.getSensorThickness(module->getType()) * module->getArea() / 1000.0 * module->getNFaces(); // volume is in cm^3
				// in which bin is the center of the module ? 
				double x  = (center.Z() - irrMinZ)/irrStepZ;   
				double y  = (center.Rho() - irrMinR)/irrStepR;

				double x1(0.0), y1(0.0),x2(0.0), y2(0.0);
				x1 = floor(x); 
				y1 = floor(y); 
				x2 = ceil(x);  
				y2 = ceil(y);  
				if (x1==x2) x2++; // to avoid division by 0 if x==int(x)
				if (y1==y2) y2++; // to avoid division by 0 if y==int(y)
				double irr11 = tracker.getIrradiationMap()[make_pair(int(x1), int(y1))]; //this is ok 
				//double irr21 = tracker.getIrradiationMap()[make_pair(int(x2), int(y1))]; // TODO : cleanup unused vars?
				//double irr12 = tracker.getIrradiationMap()[make_pair(int(x1), int(y2))];
				//double irr22 = tracker.getIrradiationMap()[make_pair(int(x2), int(y2))];

				// no need for averaging (check with Stefano) - the FLUKA map already comes with averaged values for each cell
				//double irrxy = irr11/((x2-x1)*(y2-y1))*(x2-x)*(y2-y) + irr21/((x2-x1)*(y2-y1))*(x-x1)*(y2-y) + irr12/((x2-x1)*(y2-y1))*(x2-x)*(y-y1) + irr22/((x2-x1)*(y2-y1))*(x-x1)*(y-y1); // bilinear interpolation
				double irrxy = irr11;
				// 80 mb is the current estimate for pp cross-section @ 7 TeV (2013-05-10)
				// no need for normalisation - new FLUKA map is already normalized
				//double fluence = irrxy * numInvFemtobarns * 1e15 * 80 * 1e-3; // fluence is in 1MeV-equiv-neutrons/cm^2 
				double fluence = irrxy; // fluence is in 1MeV-equiv-neutrons/cm^2 
				double leakCurrentScaled = alphaParam * fluence * volume * pow((operatingTemp+273.15) / (referenceTemp+273.15), 2) * exp(-1.21/(2*8.617334e-5)*(1/(operatingTemp+273.15)-1/(referenceTemp+273.15))); 
				double irradiatedPowerConsumption = leakCurrentScaled * chargeDepletionVoltage;         
				module->setIrradiatedPowerConsumption(irradiatedPowerConsumption);
				module->setProperty("irradiatedPowerConsumption", irradiatedPowerConsumption);
				if (!(irradiatedPowerConsumption>0))
					std::cerr << module->getSensorGeoTag() << " => " << irradiatedPowerConsumption << std::endl;

				irradiatedPowerConsumptionSummaries_[cntName].setCell(module->getLayer(), module->getRing(), irradiatedPowerConsumption);
				maxIrradiatedPowerConsumptionTableSummaries_[moduleDistance].setCell(module->getLayer(), module->getRing(), irradiatedPowerConsumption);
				maxIrradiatedPowerConsumptionSummaries_[moduleDistance].push_back(irradiatedPowerConsumption);
			}
		}
	}


    // protected
    /**
     * Looping on the layers, and picking only modules on the YZ section
     * For all these modules prepare a different table with materials
     * @param tracker a reference to the <i>ModuleCap</i> vector of vectors that sits on top of the tracker modules
     * @param result a map of summary tables to be filled
     */
    void Analyzer::computeDetailedWeights(std::vector<std::vector<ModuleCap> >& tracker,std::map<std::string, SummaryTable>& result,
                                          bool byMaterial) {
      map<std::string, map<std::pair<int, int>, bool> > typeTaken;
      map<std::string, map<std::pair<int, int>, bool> > typeWritten;

      string tempString;
      ostringstream tempSS;

      std::vector<std::vector<ModuleCap> >::iterator layerIt;
      //std::vector<std::vector<ModuleCap> >::iterator layerGuard;
      std::vector<ModuleCap>::iterator moduleIt;
      //std::vector<ModuleCap>::iterator moduleGuard;

      ModuleCap* myModuleCap;
      Module* myModule;

      //  unsigned int nLocalMasses;
      //  unsigned int nExitingMasses;

      std::map<std::string, double>::const_iterator localmassesBegin;
      std::map<std::string, double>::const_iterator localmassesEnd;

      std::map<std::string, double>::const_iterator exitingmassesBegin;
      std::map<std::string, double>::const_iterator exitingmassesEnd;

      // First create a list of material used anywhere
      std::vector<std::string> materialTagV;
      std::vector<std::string>::iterator materialTagIt;

      double localMaterial;
      double exitingMaterial;
      string materialTag;

      // loop over layers
      for (layerIt = tracker.begin(); layerIt != tracker.end(); ++layerIt) {
        // Loop over modules
        for (moduleIt = layerIt->begin(); moduleIt != layerIt->end(); ++moduleIt) {
          // Check if the module is on the YZ section
          myModuleCap = &(*moduleIt);
          myModule = &(myModuleCap->getModule());
          if (myModule->getSection()==Layer::YZSection) {
            pair<int, int> myIndex = make_pair(myModule->getLayer()/*+myModule->getDisk()*/, myModule->getRing());
            tempString = myModule->getContainerName();
            if (!typeTaken[tempString][myIndex]) {
              typeTaken[tempString][myIndex]=true;
              if (tempString!="") {
                // TODO: put this in a better place
                // (and make a better module typing)
                tempSS.str("");
                if (myModule->getSubdetectorType()==Module::Barrel) {
                  tempSS << myModule->getLayer();
                  tempString+=" (L"+tempSS.str()+")";
                } else if (myModule->getSubdetectorType()==Module::Endcap) {
                  if (EndcapModule* myEcMod = dynamic_cast<EndcapModule*>(myModule))  {
                    tempSS << myEcMod->getLayer(); //getDisk();
                    tempString+=" (D"+tempSS.str()+")";
                  } else {
                    cerr << "ERROR in Analyzer::detailedWeights(): "
                        << "I found an endcap Module that cannot be cast to an EndcapModule!" << std::endl;        
                  }
                } else {
                  cerr << "ERROR in Analyzer::detailedWeights(): "
                      << "I found a module which is neither endcap nor barrel!" << std::endl;
                }
                result[tempString].setCell(0, myModule->getRing(), myModule->getTag());
              } else {
                cerr << "ERROR in Analyzer::detailedWeights(): "
                    << "I found a module with no reference to the container name." << endl;
              }
            }
            if (byMaterial) { // sort by Material tag
              //nLocalMasses = myModuleCap->localMassCount();
              //nExitingMasses = myModuleCap->exitingMassCount();
              localmassesBegin = myModuleCap->getLocalMasses().begin();
              localmassesEnd = myModuleCap->getLocalMasses().end();
              exitingmassesBegin = myModuleCap->getExitingMasses().begin();
              exitingmassesEnd = myModuleCap->getExitingMasses().end();
            } else { // sort by Component tag
              //nLocalMasses = myModuleCap->localMassCompCount();
              //nExitingMasses = myModuleCap->exitingMassCompCount();
              localmassesBegin = myModuleCap->getLocalMassesComp().begin();
              localmassesEnd = myModuleCap->getLocalMassesComp().end();
              exitingmassesBegin = myModuleCap->getExitingMassesComp().begin();
              exitingmassesEnd = myModuleCap->getExitingMassesComp().end();
            }
            for (std::map<std::string, double>::const_iterator it = localmassesBegin; it != localmassesEnd; ++it) {
              // if (byMaterial) materialTag = myModuleCap->getLocalTag(iLocalMasses); // sort by Material tag
              //  else materialTag = myModuleCap->getLocalTagComp(iLocalMasses);           // sort by Component tag
              materialTag = it->first;
              materialTagIt = find(materialTagV.begin(), materialTagV.end(), materialTag);
              if (materialTagIt==materialTagV.end()) {
                materialTagV.push_back(materialTag);
              }
            }
            for (std::map<std::string, double>::const_iterator it = exitingmassesBegin; it != exitingmassesEnd; ++it) {
              //  if (byMaterial) materialTag = myModuleCap->getExitingTag(iExitingMasses); // sort by Material tag
              //  else materialTag = myModuleCap->getExitingTagComp(iExitingMasses);           // sort by Component tag
              materialTag = it->first;
              materialTagIt = find(materialTagV.begin(), materialTagV.end(), materialTag);
              if (materialTagIt==materialTagV.end()) {
                materialTagV.push_back(materialTag);
              }
            }
          }
        }
      }

      // Alphabetically sort materials
      std::sort(materialTagV.begin(), materialTagV.end());

      // Prepare the columns of the tables
      for (map<string, SummaryTable>::iterator it=result.begin();
           it!=result.end(); ++it) {
        for (unsigned int materialTag_i=0; materialTag_i<materialTagV.size(); ++materialTag_i) {
          it->second.setCell(materialTag_i+1, 0, materialTagV[materialTag_i]);
        }
        it->second.setCell(materialTagV.size()+1, 0, "Total");
      }

      // Now fill the table
      // loop over layers
      for (layerIt = tracker.begin(); layerIt != tracker.end(); ++layerIt) { 
        // Loop over modules
        for (moduleIt = layerIt->begin(); moduleIt != layerIt->end(); ++moduleIt) { 
          // Check if the module is on the YZ section
          myModuleCap = &(*moduleIt);
          myModule = &(myModuleCap->getModule());

          if (byMaterial) {
            typeWeight[myModule->getTag()]+=myModuleCap->getLocalMass();
            typeWeight[myModule->getTag()]+=myModuleCap->getExitingMass();
            tagWeight[myModule->getSensorGeoTag()]+=myModuleCap->getLocalMass();
            tagWeight[myModule->getSensorGeoTag()]+=myModuleCap->getExitingMass();
          }
          if (myModule->getSection()==Layer::YZSection) {
            // If we did not write this module type yet
            pair<int, int> myIndex = make_pair(myModule->getLayer()/*+myModule->getDisk()*/, myModule->getRing());
            tempString = myModule->getContainerName();
            if (!typeWritten[tempString][myIndex]) {
              typeWritten[tempString][myIndex]=true;
              if (tempString!="") {
                // TODO: put this in a better place
                // (and make a better module typing)
                tempSS.str("");
                if (myModule->getSubdetectorType()==Module::Barrel) {
                  tempSS << myModule->getLayer();
                  tempString+=" (L"+tempSS.str()+")";
                } else if (myModule->getSubdetectorType()==Module::Endcap) {
                  tempSS << myModule->getLayer(); //getDisk();
                  tempString+=" (D"+tempSS.str()+")";
                } else {
                  cerr << "ERROR in Analyzer::detailedWeights(): "
                      << "I found a module which is neither endcap nor barrel!" << std::endl;
                }
                //      cout << "Here's a module: id = " << myModule->getId()
                //           << ", tag = " << myModule->getTag()
                //           << ", type = " << myModule->getType()
                //           << ", ring = " << myModule->getRing()
                //           << endl;
                // std::cout << "Material\tLocal\texiting\n";

                // Then we fill in the proper column
                // Prepare the columns of the table
                for (unsigned int materialTag_i=0; materialTag_i<materialTagV.size(); ++materialTag_i) {
                  materialTag = materialTagV[materialTag_i];
                  if (byMaterial) { // table by materials
                    try { localMaterial = myModuleCap->getLocalMass(materialTag); }
                    catch (exception e) { localMaterial = 0; }
                    try { exitingMaterial = myModuleCap->getExitingMass(materialTag); }
                    catch (exception e) { exitingMaterial = 0; }
                  } else { // table by components
                    try { localMaterial = myModuleCap->getLocalMassComp(materialTag); }
                    catch (exception e) { localMaterial = 0; }
                    try { exitingMaterial = myModuleCap->getExitingMassComp(materialTag); }
                    catch (exception e) { exitingMaterial = 0; }
                  }
                  //cout << materialTag << "\t"
                  //<< localMaterial << "\t"
                  //<< exitingMaterial << "\t" << endl;
                  tempSS.str("");
                  // TODO: move this to Vizard
                  tempSS << std::dec << std::fixed << std::setprecision(1) << localMaterial << "+"
                    << std::dec << std::fixed << std::setprecision(1) << exitingMaterial << "="
                    << std::dec << std::fixed << std::setprecision(1) << localMaterial+exitingMaterial;
                  result[tempString].setCell(materialTag_i+1, myModule->getRing(), tempSS.str());
                }
                localMaterial = myModuleCap->getLocalMass();
                exitingMaterial = myModuleCap->getExitingMass();
                tempSS.str("");
                tempSS << std::dec << std::fixed << std::setprecision(1) << localMaterial << "+"
                  << std::dec << std::fixed << std::setprecision(1) << exitingMaterial << "="
                  << std::dec << std::fixed << std::setprecision(1) << localMaterial+exitingMaterial;
                result[tempString].setCell(materialTagV.size()+1, myModule->getRing(), tempSS.str());
              } else {
                cerr << "ERROR in Analyzer::detailedWeights(): "
                    << "I found a module with no reference to the container name." << endl;
              }
            }
          }
        }
      }
    }

    // public
    /**
     * Produces a full material summary for the modules
     */
    void Analyzer::computeWeightSummary(MaterialBudget& mb) {

      typeWeight.clear();
      tagWeight.clear();
      barrelWeights.clear();
      computeDetailedWeights(mb.getBarrelModuleCaps(), barrelWeights, true);
      endcapWeights.clear();
      computeDetailedWeights(mb.getEndcapModuleCaps(), endcapWeights, true);

      barrelComponentWeights.clear();
      computeDetailedWeights(mb.getBarrelModuleCaps(), barrelComponentWeights, false);
      endcapComponentWeights.clear();
      computeDetailedWeights(mb.getEndcapModuleCaps(), endcapComponentWeights, false);
    }

    // protected
    /**
     * The layer-level analysis function for modules forms the frame for sending a single track through the active modules.
     * It loops through all layers and adds up the results returned from the analysis of each of them.
     * @param tr A reference to the <i>ModuleCap</i> vector of vectors that sits on top of the tracker modules
     * @param eta The pseudorapidity of the track
     * @param theta The track angle in the yz-plane
     * @param phi The track angle in the xy-plane
     * @param t A reference to the current track object
     * @param A boolean flag to indicate which set of active surfaces is analysed: true if the belong to a pixel detector, false if they belong to the tracker
     * @return The summed up radiation and interaction lengths for the given track, bundled into a <i>std::pair</i>
     */
    Material Analyzer::analyzeModules(std::vector<std::vector<ModuleCap> >& tr,
                                      double eta, double theta, double phi, Track& t, 
                                      std::map<std::string, Material>& sumComponentsRI,
                                      bool isPixel) {
      std::vector<std::vector<ModuleCap> >::iterator iter = tr.begin();
      std::vector<std::vector<ModuleCap> >::iterator guard = tr.end();
      Material res, tmp;
      res.radiation= 0.0;
      res.interaction = 0.0;
      while (iter != guard) {
        tmp = findModuleLayerRI(*iter, eta, theta, phi, t, sumComponentsRI, isPixel);
        res.radiation= res.radiation+ tmp.radiation;
        res.interaction= res.interaction + tmp.interaction;
        iter++;
      }
      return res;
    }


    /**
     * The module-level analysis function loops through all modules of a given layer, checking for collisions with the given track.
     * If one is found, the radiation and interaction lengths are scaled with respect to theta, then summed up into a grand total,
     * which is returned. As phi is fixed at the moment and the tracks hit the modules orthogonally with respect to it, it is so far
     * not used to scale the results further.
     * @param layer A reference to the <i>ModuleCap</i> vector linking the collection of material properties to the current layer
     * @param eta The pseudorapidity of the current track
     * @param theta The track angle in the yz-plane
     * @param phi The track angle in the xy-plane
     * @param t A reference to the current track object
     * @param A boolean flag to indicate which set of active surfaces is analysed: true if the belong to a pixel detector, false if they belong to the tracker
     * @return The scaled and summed up radiation and interaction lengths for the given layer and track, bundled into a <i>std::pair</i>
     */
    Material Analyzer::findModuleLayerRI(std::vector<ModuleCap>& layer,
                                         double eta, double theta, double phi, Track& t, 
                                         std::map<std::string, Material>& sumComponentsRI,
                                         bool isPixel) {
      std::vector<ModuleCap>::iterator iter = layer.begin();
      std::vector<ModuleCap>::iterator guard = layer.end();
      Material res, tmp;
      XYZVector origin, direction;
      Polar3DVector dir;
      double distance, r;
      int hits = 0;
      res.radiation = 0.0;
      res.interaction = 0.0;
      // set the track direction vector
      dir.SetCoordinates(1, theta, phi);
      direction = dir;
      while (iter != guard) {
        // collision detection: rays are in z+ only, so consider only modules that lie on that side
        // only consider modules that have type BarrelModule or EndcapModule
        if (iter->getModule().getMaxZ() > 0) {
          if ((iter->getModule().getSubdetectorType() == Module::Barrel) ||
              (iter->getModule().getSubdetectorType() == Module::Endcap)) {
            // same method as in Tracker, same function used
            // TODO: in case origin==0,0,0 and phi==0 just check if sectionYZ and minEta, maxEta
            distance = iter->getModule().trackCross(origin, direction);
            if (distance > 0) {
              // module was hit
              hits++;
              r = distance * sin(theta);
              tmp.radiation = iter->getRadiationLength();
              tmp.interaction = iter->getInteractionLength();
              // 2D material maps
              fillMapRT(r, theta, tmp);
              // radiation and interaction length scaling for barrels
              if (iter->getModule().getSubdetectorType() == Module::Barrel) {
                tmp.radiation = tmp.radiation / sin(theta);
                tmp.interaction = tmp.interaction / sin(theta);
              }
              // radiation and interaction length scaling for endcaps
              else {
                tmp.radiation = tmp.radiation / cos(theta);
                tmp.interaction = tmp.interaction / cos(theta);
              }

              double tmpr = 0., tmpi = 0.;

              std::map<std::string, Material> moduleComponentsRI = iter->getComponentsRI();
              for (std::map<std::string, Material>::iterator cit = moduleComponentsRI.begin(); cit != moduleComponentsRI.end(); ++cit) {
                sumComponentsRI[cit->first].radiation += cit->second.radiation / (iter->getModule().getSubdetectorType() == Module::Barrel ? sin(theta) : cos(theta));
                tmpr += sumComponentsRI[cit->first].radiation;
                sumComponentsRI[cit->first].interaction += cit->second.interaction / (iter->getModule().getSubdetectorType() == Module::Barrel ? sin(theta) : cos(theta));
                tmpi += sumComponentsRI[cit->first].interaction;
              }
              // 2D plot and eta plot results
              if (!isPixel) fillCell(r, eta, theta, tmp);
              res += tmp;
              // create Hit object with appropriate parameters, add to Track t
              Hit* hit = new Hit(distance, &(iter->getModule()));
              //if (iter->getModule().getSubdetectorType() == Module::Barrel) hit->setOrientation(Hit::Horizontal); // should not be necessary
              //else if(iter->getModule().getSubdetectorType() == Module::Endcap) hit->setOrientation(Hit::Vertical); // should not be necessary
              //hit->setObjectKind(Hit::Active); // should not be necessary
              hit->setCorrectedMaterial(tmp);
              hit->setPixel(isPixel);
              t.addHit(hit);
            }
          }
          else std::cout << msg_module_warning << std::endl;
        }
        iter++;
      }
      return res;
    }


    // protected
    /**
     * The layer-level function to find hits on modules and connect them to a track
     * It loops through all layers and adds up the results returned from the analysis of each of them.
     * @param tr A reference to the <i>ModuleCap</i> vector of vectors that sits on top of the tracker modules
     * @param eta The pseudorapidity of the track
     * @param theta The track angle in the yz-plane
     * @param phi The track angle in the xy-plane
     * @param t A reference to the current track object
     * @param A boolean flag to indicate which set of active surfaces is analysed: true if the belong to a pixel detector, false if they belong to the tracker
     * @return The summed up radiation and interaction lengths for the given track, bundled into a <i>std::pair</i>
     */
    Material Analyzer::findHitsModules(std::vector<std::vector<ModuleCap> >& tr,
                                       // TODO: add z0 here and in the hit finder for inactive surfaces
                                       double eta, double theta, double phi, Track& t, bool isPixel) {
      std::vector<std::vector<ModuleCap> >::iterator iter = tr.begin();
      std::vector<std::vector<ModuleCap> >::iterator guard = tr.end();
      Material res, tmp;
      res.radiation= 0.0;
      res.interaction = 0.0;
      while (iter != guard) {
        tmp = findHitsModuleLayer(*iter, eta, theta, phi, t, isPixel);
        res.radiation = res.radiation + tmp.radiation;
        res.interaction = res.interaction + tmp.interaction;
        iter++;
      }
      return res;
    }

    int Analyzer::findHitsModules(Tracker& tracker, double z0, double eta, double theta, double phi, Track& t) {

      Material emptyMaterial;
      XYZVector origin(0,0,z0);
      XYZVector direction;
      Polar3DVector dir;
      double distance;

      int hits = 0;
      emptyMaterial.radiation = 0.0;
      emptyMaterial.interaction = 0.0;

      // set the track direction vector
      dir.SetCoordinates(1, theta, phi);
      direction = dir;

      LayerVector allLayers = tracker.getLayers();
      LayerVector::iterator layIt;
      for (layIt=allLayers.begin(); layIt!=allLayers.end(); ++layIt) {
        Layer* aLayer = *layIt;
        ModuleVector* layerModules = aLayer->getModuleVector();
        ModuleVector::iterator modIt;
        for (modIt=layerModules->begin(); modIt!=layerModules->end(); ++modIt) {
          Module* aModule = *modIt;

          // collision detection: rays are in z+ only, so consider only modules that lie on that side
          // only consider modules that have type BarrelModule or EndcapModule
          if (aModule->getMaxZ() > 0) {

            // same method as in Tracker, same function used
            distance = aModule->trackCross(origin, direction);
            if (distance > 0) {
              // module was hit
              hits++;

              // create Hit object with appropriate parameters, add to Track t
              Hit* hit = new Hit(distance, aModule);
              hit->setCorrectedMaterial(emptyMaterial);
              t.addHit(hit);
            }
          }
        }
      }
      return hits;
    }


    /**
     * The module-level analysis function loops through all modules of a given layer, finds hits and connects them to the track.
     * If one is found, the radiation and interaction lengths are scaled with respect to theta, then summed up into a grand total,
     * which is returned. As phi is fixed at the moment and the tracks hit the modules orthogonally with respect to it, it is so far
     * not used to scale the results further.
     * @param layer A reference to the <i>ModuleCap</i> vector linking the collection of material properties to the current layer
     * @param eta The pseudorapidity of the current track
     * @param theta The track angle in the yz-plane
     * @param phi The track angle in the xy-plane
     * @param t A reference to the current track object
     * @param A boolean flag to indicate which set of active surfaces is analysed: true if the belong to a pixel detector, false if they belong to the tracker
     * @return The scaled and summed up radiation and interaction lengths for the given layer and track, bundled into a <i>std::pair</i>
     */
    Material Analyzer::findHitsModuleLayer(std::vector<ModuleCap>& layer,
                                           double eta, double theta, double phi, Track& t, bool isPixel) {
      std::vector<ModuleCap>::iterator iter = layer.begin();
      std::vector<ModuleCap>::iterator guard = layer.end();
      Material res, tmp;
      XYZVector origin, direction;
      Polar3DVector dir;
      double distance;
      //double r;
      int hits = 0;
      res.radiation = 0.0;
      res.interaction = 0.0;
      // set the track direction vector
      dir.SetCoordinates(1, theta, phi);
      direction = dir;
      while (iter != guard) {
        // collision detection: rays are in z+ only, so consider only modules that lie on that side
        // only consider modules that have type BarrelModule or EndcapModule
        if (iter->getModule().getMaxZ() > 0) {
          if ((iter->getModule().getSubdetectorType() == Module::Barrel) ||
              (iter->getModule().getSubdetectorType() == Module::Endcap)) {
            // same method as in Tracker, same function used
            // TODO: in case origin==0,0,0 and phi==0 just check if sectionYZ and minEta, maxEta
            distance = iter->getModule().trackCross(origin, direction);
            if (distance > 0) {
              // module was hit
              hits++;
              // r = distance * sin(theta);
              tmp.radiation = iter->getRadiationLength();
              tmp.interaction = iter->getInteractionLength();
              // radiation and interaction length scaling for barrels
              if (iter->getModule().getSubdetectorType() == Module::Barrel) {
                tmp.radiation = tmp.radiation / sin(theta);
                tmp.interaction = tmp.interaction / sin(theta);
              }
              // radiation and interaction length scaling for endcaps
              else {
                tmp.radiation = tmp.radiation / cos(theta);
                tmp.interaction = tmp.interaction / cos(theta);
              }
              res += tmp;
              // create Hit object with appropriate parameters, add to Track t
              Hit* hit = new Hit(distance, &(iter->getModule()));
              //if (iter->getModule().getSubdetectorType() == Module::Barrel) hit->setOrientation(Hit::Horizontal); // should not be necessary
              //else if(iter->getModule().getSubdetectorType() == Module::Endcap) hit->setOrientation(Hit::Vertical); // should not be necessary
              //hit->setObjectKind(Hit::Active); // should not be necessary
              hit->setCorrectedMaterial(tmp);
              hit->setPixel(isPixel);
              t.addHit(hit);
            }
          }
          else std::cout << msg_module_warning << std::endl;
        }
        iter++;
      }
      return res;
    }

    /**
     * The analysis function for inactive volumes loops through the given vector of elements, checking for collisions with
     * the given track. If one is found, the radiation and interaction lengths are scaled with respect to theta, then summed
     * up into a grand total, which is returned. As all inactive volumes are symmetric with respect to rotation around the
     * z-axis, the track angle phi is not necessary.
     * @param elements A reference to the collection of inactive surfaces that is to be checked for collisions with the track
     * @param eta The pseudorapidity of the current track
     * @param theta The track angle in the yz-plane
     * @param t A reference to the current track object
     * @param cat The category of inactive surfaces that need to be considered within the collection; none if the function is to look at all of them
     * @param A boolean flag to indicate which set of active surfaces is analysed: true if the belong to a pixel detector, false if they belong to the tracker
     * @return The scaled and summed up radiation and interaction lengths for the given collection of elements and track, bundled into a <i>std::pair</i>
     */

    Material Analyzer::analyzeInactiveSurfaces(std::vector<InactiveElement>& elements, double eta,
                                               double theta, Track& t, MaterialProperties::Category cat, bool isPixel) {
      std::vector<InactiveElement>::iterator iter = elements.begin();
      std::vector<InactiveElement>::iterator guard = elements.end();
      Material res, corr;
      std::pair<double, double> tmp;
      double s = 0.0;
      while (iter != guard) {
        // collision detection: rays are in z+ only, so only volumes in z+ need to be considered
        // only volumes of the requested category, or those without one (which should not exist) are examined
        if (((iter->getZOffset() + iter->getZLength()) > 0)
            && ((cat == MaterialProperties::no_cat) || (cat == iter->getCategory()))) {
          // collision detection: check eta range
          tmp = iter->getEtaMinMax();
          // volume was hit
          if ((tmp.first < eta) && (tmp.second > eta)) {
            double r, z;
            // radiation and interaction lenth scaling for vertical volumes
            if (iter->isVertical()) {
              z = iter->getZOffset() + iter->getZLength() / 2.0;
              r = z * tan(theta);
              // 2D maps for vertical surfaces
              fillMapRZ(r,z,iter->getMaterialLengths());
              // special treatment for user-defined supports as they can be very close to z=0
              if (cat == MaterialProperties::u_sup) {
                s = iter->getZLength() / cos(theta);
                if (s > (iter->getRWidth() / sin(theta))) s = iter->getRWidth() / sin(theta);
                // add the hit if it's declared as inside the tracking volume, add it to 'others' if not
                if (iter->track()) {
                  corr.radiation = iter->getRadiationLength() * s / iter->getZLength();
                  corr.interaction = iter->getInteractionLength() * s / iter->getZLength();
                  res += corr;
                  if (!isPixel) {
                    Material thisLength;
                    thisLength.radiation = iter->getRadiationLength() * s / iter->getZLength();
                    thisLength.interaction = iter->getInteractionLength() * s / iter->getZLength(); 
                    fillCell(r, eta, theta, thisLength); 
                  }
                }
                else {
                  if (!isPixel) {
                    rextrasupports.Fill(eta, iter->getRadiationLength() * s / iter->getZLength());
                    iextrasupports.Fill(eta, iter->getInteractionLength() * s / iter->getZLength());
                  }
                }
              }
              else {
                // add the hit if it's declared as inside the tracking volume, add it to 'others' if not
                if (iter->track()) {
                  corr.radiation = iter->getRadiationLength() / cos(theta);
                  corr.interaction = iter->getInteractionLength() / cos(theta);
                  res += corr;
                  if (!isPixel) {
                    Material thisLength;
                    thisLength.radiation = iter->getRadiationLength() / cos(theta); 
                    thisLength.interaction = iter->getInteractionLength() / cos(theta);
                    fillCell(r, eta, theta, thisLength);
                  }
                }
                else {
                  if (!isPixel) {
                    if ((iter->getCategory() == MaterialProperties::b_ser)
                        || (iter->getCategory() == MaterialProperties::e_ser)) {
                      rextraservices.Fill(eta, iter->getRadiationLength() / cos(theta));
                      iextraservices.Fill(eta, iter->getInteractionLength() / cos(theta));
                    }
                    else if ((iter->getCategory() == MaterialProperties::b_sup)
                             || (iter->getCategory() == MaterialProperties::e_sup)
                             || (iter->getCategory() == MaterialProperties::o_sup)
                             || (iter->getCategory() == MaterialProperties::t_sup)) {
                      rextrasupports.Fill(eta, iter->getRadiationLength() / cos(theta));
                      iextrasupports.Fill(eta, iter->getInteractionLength() / cos(theta));
                    }
                  }
                }
              }
            }
            // radiation and interaction length scaling for horizontal volumes
            else {
              r = iter->getInnerRadius() + iter->getRWidth() / 2.0;
              // 2D maps for horizontal surfaces
              fillMapRT(r,theta,iter->getMaterialLengths());
              // special treatment for user-defined supports; should not be necessary for now
              // as all user-defined supports are vertical, but just in case...
              if (cat == MaterialProperties::u_sup) {
                s = iter->getZLength() / sin(theta);
                if (s > (iter->getRWidth() / cos(theta))) s = iter->getRWidth() / cos(theta);
                // add the hit if it's declared as inside the tracking volume, add it to 'others' if not
                if (iter->track()) {
                  corr.radiation = iter->getRadiationLength() * s / iter->getZLength();
                  corr.interaction = iter->getInteractionLength() * s / iter->getZLength();
                  res += corr;
                  if (!isPixel) {
                    Material thisLength;
                    thisLength.radiation = iter->getRadiationLength() * s / iter->getZLength(); 
                    thisLength.interaction = iter->getInteractionLength() * s / iter->getZLength();
                    fillCell(r, eta, theta, thisLength);
                  }
                }
                else {
                  if (!isPixel) {
                    rextrasupports.Fill(eta, iter->getRadiationLength() * s / iter->getZLength());
                    iextrasupports.Fill(eta, iter->getInteractionLength() * s / iter->getZLength());
                  }
                }
              }
              else {
                // add the hit if it's declared as inside the tracking volume, add it to 'others' if not
                if (iter->track()) {
                  corr.radiation = iter->getRadiationLength() / sin(theta);
                  corr.interaction = iter->getInteractionLength() / sin(theta);
                  res += corr;
                  if (!isPixel) {
                    Material thisLength;
                    thisLength.radiation = iter->getRadiationLength() / sin(theta);
                    thisLength.interaction =  iter->getInteractionLength() / sin(theta);
                    fillCell(r, eta, theta, thisLength); 
                  }
                }
                else {
                  if (!isPixel) {
                    if ((iter->getCategory() == MaterialProperties::b_ser)
                        || (iter->getCategory() == MaterialProperties::e_ser)) {
                      rextraservices.Fill(eta, iter->getRadiationLength() / sin(theta));
                      iextraservices.Fill(eta, iter->getInteractionLength() / sin(theta));
                    }
                    else if ((iter->getCategory() == MaterialProperties::b_sup)
                             || (iter->getCategory() == MaterialProperties::e_sup)
                             || (iter->getCategory() == MaterialProperties::o_sup)
                             || (iter->getCategory() == MaterialProperties::t_sup)) {
                      rextrasupports.Fill(eta, iter->getRadiationLength() / sin(theta));
                      iextrasupports.Fill(eta, iter->getInteractionLength() / sin(theta));
                    }
                  }
                }
              }
            }
            // create Hit object with appropriate parameters, add to Track t
            Hit* hit = new Hit((theta == 0) ? r : (r / sin(theta)));
            if (iter->isVertical()) hit->setOrientation(Hit::Vertical);
            else hit->setOrientation(Hit::Horizontal);
            hit->setObjectKind(Hit::Inactive);
            hit->setCorrectedMaterial(corr);
            hit->setPixel(isPixel);
            t.addHit(hit);
          }
        }
        iter++;
      }
      return res;
    }

    /**
     * The analysis function for inactive volumes loops through the given vector of elements, checking for collisions with
     * the given track. If one is found, the radiation and interaction lengths are scaled with respect to theta.
     * All hits are added to the given track
     * The total crossed material is returned.
     * As all inactive volumes are symmetric with respect to rotation around the z-axis, the track angle phi is not necessary.
     * @param elements A reference to the collection of inactive surfaces that is to be checked for collisions with the track
     * @param eta The pseudorapidity of the current track
     * @param theta The track angle in the yz-plane
     * @param t A reference to the current track object
     * @return The scaled and summed up crossed material amount
     */
    Material Analyzer::findHitsInactiveSurfaces(std::vector<InactiveElement>& elements, double eta,
                                                double theta, Track& t, bool isPixel) {
      std::vector<InactiveElement>::iterator iter = elements.begin();
      std::vector<InactiveElement>::iterator guard = elements.end();
      Material res, corr;
      std::pair<double, double> tmp;
      double s_normal = 0;
      double s_alternate = 0;
      while (iter != guard) {
        // Collision detection: rays are in z+ only, so only volumes in z+ need to be considered
        // only volumes of the requested category, or those without one (which should not exist) are examined
        if ((iter->getZOffset() + iter->getZLength()) > 0) {
          // collision detection: check eta range
          tmp = iter->getEtaMinMax();
          // Volume was hit if:
          if ((tmp.first < eta) && (tmp.second > eta)) {
            double r, z;
            // radiation and interaction lenth scaling for vertical volumes
            if (iter->isVertical()) { // Element is vertical
              z = iter->getZOffset() + iter->getZLength() / 2.0;
              r = z * tan(theta);

              // In case we are crossing the material with a very shallow angle
              // we have to take into account its finite radial size
              s_normal = iter->getZLength() / cos(theta);
              s_alternate = iter->getRWidth() / sin(theta);
              if (s_normal > s_alternate) { 
                // Special case: it's easier to cross the material by going left-to-right than
                // by going bottom-to-top, so I have to rescale the material amount computation
                corr.radiation = iter->getRadiationLength() / iter->getZLength() * s_alternate;
                corr.interaction = iter->getInteractionLength() / iter->getZLength() * s_alternate;
                res += corr;
              } else {
                // Standard computing of the crossed material amount
                corr.radiation = iter->getRadiationLength() / cos(theta);
                corr.interaction = iter->getInteractionLength() / cos(theta);
                res += corr;
              }
            }
            // radiation and interaction length scaling for horizontal volumes
            else { // Element is horizontal
              r = iter->getInnerRadius() + iter->getRWidth() / 2.0;

              // In case we are crossing the material with a very shallow angle
              // we have to take into account its finite z length
              s_normal = iter->getRWidth() / sin(theta);
              s_alternate = iter->getZLength() / cos(theta);
              if (s_normal > s_alternate) { 
                // Special case: it's easier to cross the material by going left-to-right than
                // by going bottom-to-top, so I have to rescale the material amount computation
                corr.radiation = iter->getRadiationLength() / iter->getRWidth() * s_alternate;
                corr.interaction = iter->getInteractionLength() / iter->getRWidth() * s_alternate;
                res += corr;
              } else {
                // Standard computing of the crossed material amount
                corr.radiation = iter->getRadiationLength() / sin(theta);
                corr.interaction = iter->getInteractionLength() / sin(theta);
                res += corr;
              }
            }
            // create Hit object with appropriate parameters, add to Track t
            Hit* hit = new Hit((theta == 0) ? r : (r / sin(theta)));
            if (iter->isVertical()) hit->setOrientation(Hit::Vertical);
            else hit->setOrientation(Hit::Horizontal);
            hit->setObjectKind(Hit::Inactive);
            hit->setCorrectedMaterial(corr);
            hit->setPixel(isPixel);
            t.addHit(hit);
          }
        }
        iter++;
      }
      return res;
    }

    /**
     * Calculate the error graphs for the radius curvature, the distance and the angle, for each momentum,
     * and store them internally for later visualisation.
     * @param p The list of different momenta that the error graphs are calculated for
     */
    void Analyzer::calculateGraphs(const std::vector<double>& p, 
                                   const std::vector<Track>& trackVector,
                                   int graphAttributes) {

      std::map<double, TGraph>& thisRhoGraphs = myGraphBag.getGraphs(graphAttributes | graphBag::RhoGraph );
      std::map<double, TGraph>& thisPhiGraphs = myGraphBag.getGraphs(graphAttributes | graphBag::PhiGraph );
      std::map<double, TGraph>& thisDGraphs = myGraphBag.getGraphs(graphAttributes | graphBag::DGraph );
      std::map<double, TGraph>& thisCtgThetaGraphs = myGraphBag.getGraphs(graphAttributes | graphBag::CtgthetaGraph );
      std::map<double, TGraph>& thisZ0Graphs = myGraphBag.getGraphs(graphAttributes | graphBag::Z0Graph );
      std::map<double, TGraph>& thisPGraphs = myGraphBag.getGraphs(graphAttributes | graphBag::PGraph );

      std::map<double, double>::const_iterator miter, mguard;
      std::vector<double>::const_iterator iter, guard = p.end();
      int n = trackVector.size();
      double eta, R;
      thisRhoGraphs.clear();
      thisPhiGraphs.clear();
      thisDGraphs.clear();
      thisCtgThetaGraphs.clear();
      thisZ0Graphs.clear();
      thisPGraphs.clear();
      // momentum loop
      std::ostringstream aName;
      for (iter = p.begin(); iter != guard; iter++) {
        std::pair<double, TGraph> elem;
        TGraph graph;
        elem.first = *iter;
        elem.second = graph;
        // Prepare plots: pT
        thisRhoGraphs.insert(elem);
        thisRhoGraphs[elem.first].SetTitle("Transverse momentum error;#eta;#sigma (#delta p_{T}/p_{T}) [%]");
        aName.str(""); aName << "pt_vs_eta" << *iter;
        thisRhoGraphs[elem.first].SetName(aName.str().c_str());
        // Prepare plots: phi
        thisPhiGraphs.insert(elem);
        thisPhiGraphs[elem.first].SetTitle("Track azimuthal angle error;#eta;#sigma (#delta #phi) [rad]");
        aName.str(""); aName << "phi_vs_eta" << *iter;
        thisPhiGraphs[elem.first].SetName(aName.str().c_str());
        // Prepare plots: d
        thisDGraphs.insert(elem);
        thisDGraphs[elem.first].SetTitle("Transverse impact parameter error;#eta;#sigma (#delta d_{0}) [cm]");
        aName.str(""); aName << "d_vs_eta" << *iter;
        thisDGraphs[elem.first].SetName(aName.str().c_str());
        // Prepare plots: ctg(theta)
        thisCtgThetaGraphs.insert(elem);
        thisCtgThetaGraphs[elem.first].SetTitle("Track polar angle error;#eta;#sigma (#delta ctg(#theta))");
        aName.str(""); aName << "ctgTheta_vs_eta" << *iter;
        thisCtgThetaGraphs[elem.first].SetName(aName.str().c_str());
        // Prepare plots: z0
        thisZ0Graphs.insert(elem);
        thisZ0Graphs[elem.first].SetTitle("Longitudinal impact parameter error;#eta;#sigma (#delta z_{0}) [cm]");
        aName.str(""); aName << "z_vs_eta" << *iter;
        thisZ0Graphs[elem.first].SetName(aName.str().c_str());
        // Prepare plots: p
        thisPGraphs.insert(elem);
        thisPGraphs[elem.first].SetTitle("Momentum error;#eta;#sigma (#delta p/p) [%]");
        aName.str(""); aName << "p_vs_eta" << *iter;
        thisPGraphs[elem.first].SetName(aName.str().c_str());
      }
      // track loop
      std::map<double,int> rhoPointCount;
      std::map<double,int> phiPointCount;
      std::map<double,int> dPointCount;
      std::map<double,int> ctgPointCount;
      std::map<double,int> z0PointCount;
      std::map<double,int> pPointCount;
      double graphValue;
      for (int i = 0; i < n; i++) {
        const std::map<double, double>& drho = trackVector.at(i).getDeltaRho();
        const std::map<double, double>& dphi = trackVector.at(i).getDeltaPhi();
        const std::map<double, double>& dd = trackVector.at(i).getDeltaD();
        const std::map<double, double>& dctg = trackVector.at(i).getDeltaCtgTheta();
        const std::map<double, double>& dz0 = trackVector.at(i).getDeltaZ0();
        const std::map<double, double>& dp = trackVector.at(i).getDeltaP();
        eta = - log(tan(trackVector.at(i).getTheta() / 2));
        mguard = drho.end();
        // error by momentum loop
        for (miter = drho.begin(); miter != mguard; miter++) {
          if (thisRhoGraphs.find(miter->first) != thisRhoGraphs.end()) {
            R = miter->first / magnetic_field / 0.3 * 1E3; // radius in mm
            if ((miter->second)>0) {
              // deltaRho / rho = deltaRho * R
              graphValue = (miter->second * R) * 100; // in percent
              //if (miter->first>50) std::cout << "eta = " << eta << ", error = " << graphValue << std::endl;
              thisRhoGraphs[miter->first].SetPoint(rhoPointCount[miter->first]++, eta, graphValue);
            }
          }
        }
        mguard = dphi.end();
        for (miter = dphi.begin(); miter != mguard; miter++) {
          if (thisPhiGraphs.find(miter->first) != thisPhiGraphs.end())
            if ((miter->second)>0) {
              graphValue = miter->second; // radians is ok
              thisPhiGraphs[miter->first].SetPoint(phiPointCount[miter->first]++, eta, graphValue);
            }
        }
        mguard = dd.end();
        for (miter = dd.begin(); miter != mguard; miter++) {
          if (thisDGraphs.find(miter->first) != thisDGraphs.end())
            if ((miter->second)>0) {
              graphValue =  (miter->second) / 10.; // in cm
              thisDGraphs[miter->first].SetPoint(dPointCount[miter->first]++, eta, graphValue );
            }
        }
        mguard = dctg.end();
        for (miter = dctg.begin(); miter != mguard; miter++) {
          // Ctg theta (absolute number)
          if (thisCtgThetaGraphs.find(miter->first) != thisCtgThetaGraphs.end()) {
            graphValue = miter->second; // An absolute number
            thisCtgThetaGraphs[miter->first].SetPoint(ctgPointCount[miter->first]++, eta, graphValue);
          }
        }
        mguard = dz0.end();
        for (miter = dz0.begin(); miter != mguard; miter++) {
          if (thisZ0Graphs.find(miter->first) != thisZ0Graphs.end()) {
            graphValue =  (miter->second) / 10.; // in cm
            thisZ0Graphs[miter->first].SetPoint(z0PointCount[miter->first]++, eta, graphValue);
          }
        }       
        mguard = dp.end();
        for (miter = dp.begin(); miter != mguard; miter++) {
          if (thisPGraphs.find(miter->first) != thisPGraphs.end()) {
            graphValue =  (miter->second) * 100.; // in percent 
            thisPGraphs[miter->first].SetPoint(pPointCount[miter->first]++, eta, graphValue);
          }
        }
      }
    }

    /**
     * This convenience function resets and empties all histograms for the
     * material budget, so they are ready for a new round of analysis.
     */
    void Analyzer::clearMaterialBudgetHistograms() {
      // single category
      ractivebarrel.Reset();
      ractivebarrel.SetNameTitle("ractivebarrels", "Barrel Modules Radiation Length");
      ractiveendcap.Reset();
      ractiveendcap.SetNameTitle("ractiveendcap", "Endcap Modules Radiation Length");
      rserfbarrel.Reset();
      rserfbarrel.SetNameTitle("rserfbarrel", "Barrel Services Radiation Length");
      rserfendcap.Reset();
      rserfendcap.SetNameTitle("rserfendcap", "Endcap Services Radiation Length");
      rlazybarrel.Reset();
      rlazybarrel.SetNameTitle("rlazybarrel", "Barrel Supports Radiation Length");
      rlazyendcap.Reset();
      rlazyendcap.SetNameTitle("rlazyendcap", "Endcap Supports Radiation Length");
      rlazytube.Reset();
      rlazytube.SetNameTitle("rlazytube", "Support Tubes Radiation Length");
      rlazyuserdef.Reset();
      rlazyuserdef.SetNameTitle("rlazyuserdef", "Userdefined Supports Radiation Length");
      iactivebarrel.Reset();
      iactivebarrel.SetNameTitle("iactivebarrel", "Barrel Modules Interaction Length");
      iactiveendcap.Reset();
      iactiveendcap.SetNameTitle("iactiveendcap", "Endcap Modules Interaction Length");
      iserfbarrel.Reset();
      iserfbarrel.SetNameTitle("iserfbarrel", "Barrel Services Interaction Length");
      iserfendcap.Reset();
      iserfendcap.SetNameTitle("iserfendcap", "Endcap Services Interaction Length");
      ilazybarrel.Reset();
      ilazybarrel.SetNameTitle("ilazybarrel", "Barrel Supports Interaction Length");
      ilazyendcap.Reset();
      ilazyendcap.SetNameTitle("ilazyendcap", "Endcap Supports Interaction Length");
      ilazytube.Reset();
      ilazytube.SetNameTitle("ilazytube", "Support Tubes Interaction Length");
      ilazyuserdef.Reset();
      ilazyuserdef.SetNameTitle("ilazyuserdef", "Userdefined Supports Interaction Length");
      // composite
      rbarrelall.Reset();
      rbarrelall.SetNameTitle("rbarrelall", "Barrel Radiation Length");
      rendcapall.Reset();
      rendcapall.SetNameTitle("rendcapall", "Endcap Radiation Length");
      ractiveall.Reset();
      ractiveall.SetNameTitle("ractiveall", "Modules Radiation Length");
      rserfall.Reset();
      rserfall.SetNameTitle("rserfall", "Services Radiation Length");
      rlazyall.Reset();
      rlazyall.SetNameTitle("rlazyall", "Supports Radiation Length");
      ibarrelall.Reset();
      ibarrelall.SetNameTitle("ibarrelall", "Barrel Interaction Length");
      iendcapall.Reset();
      iendcapall.SetNameTitle("iendcapall", "Endcap Interaction Length");
      iactiveall.Reset();
      iactiveall.SetNameTitle("iactiveall", "Modules Interaction Length");
      iserfall.Reset();
      iserfall.SetNameTitle("iserfall", "Services Interaction Length");
      ilazyall.Reset();
      ilazyall.SetNameTitle("ilazyall", "Supports Interaction Length");
      // outside tracking volume
      rextraservices.Reset();
      rextraservices.SetNameTitle("rextraservices", "Services Outside Tracking Volume: Radiation Length");
      rextrasupports.Reset();
      rextrasupports.SetNameTitle("rextrasupports", "Supports Outside Tracking Volume: Radiation Length");
      iextraservices.Reset();
      iextraservices.SetNameTitle("iextraservices", "Services Outside Tracking Volume: Interaction Length");
      iextrasupports.Reset();
      iextrasupports.SetNameTitle("iextrasupports", "Supports Outside Tracking Volume: Interaction Length");
      // global
      rglobal.Reset();
      rglobal.SetNameTitle("rglobal", "Overall Radiation Length");
      iglobal.Reset();
      iglobal.SetNameTitle("iglobal", "Overall Interaction Length");
      // isolines
      isor.Reset();
      isor.SetNameTitle("isor", "Radiation Length Contours");
      isoi.Reset();
      isoi.SetNameTitle("isoi", "Interaction Length Contours");
      mapRadiation.Reset();
      mapRadiation.SetName("mapRadiation");
      mapRadiation.SetTitle("Radiation length map (raw);z [mm];r [mm]");
      mapInteraction.Reset();
      mapInteraction.SetName("mapInteraction");
      mapInteraction.SetTitle("Interaction length map (raw);z [mm];r [mm]");
      mapRadiationCount.Reset();
      mapRadiationCount.SetName("mapRadiationCount");
      mapRadiationCount.SetTitle("Radiation length hit count map;z [mm];r [mm]");
      mapInteractionCount.Reset();
      mapInteractionCount.SetName("mapInteractionCount");
      mapInteractionCount.SetTitle("Interaction length hit count map;z [mm];r [mm]");
      mapRadiationCalib.Reset();
      mapRadiationCalib.SetName("mapRadiationCalib");
      mapRadiationCalib.SetTitle("Radiation length map;z [mm];r [mm]");
      mapInteractionCalib.Reset();
      mapInteractionCalib.SetName("mapInteractionCalib");
      mapInteractionCalib.SetTitle("Interaction length map;z [mm];r [mm]");

      // Nuclear interactions
      while (hadronTotalHitsGraph.GetN()) hadronTotalHitsGraph.RemovePoint(0);
      hadronTotalHitsGraph.SetName("hadronTotalHitsGraph");
      while (hadronAverageHitsGraph.GetN()) hadronAverageHitsGraph.RemovePoint(0);
      hadronAverageHitsGraph.SetName("hadronAverageHitsGraph");

      // Clear the list of requested good hadron hits
      hadronNeededHitsFraction.clear();
      hadronGoodTracksFraction.clear();

      hadronNeededHitsFraction.push_back(ZeroHitsRequired);
      hadronNeededHitsFraction.push_back(OneHitRequired);
      //hadronNeededHitsFraction.push_back(.33);
      hadronNeededHitsFraction.push_back(.66);
      hadronNeededHitsFraction.push_back(1);

      std::sort(hadronNeededHitsFraction.begin(),
                hadronNeededHitsFraction.end());

      // Prepare the plots for the track survival fraction
      ostringstream tempSS;
      string tempString;
      TGraph* myGraph;
      for (unsigned int i=0;
           i<hadronNeededHitsFraction.size();
           ++i) {
        tempSS.str("");
        tempSS << "hadronGoodTracksFraction_at"
          << hadronNeededHitsFraction.at(i);
        tempString = tempSS.str();
        myGraph = new TGraph();
        myGraph->SetName(tempString.c_str());
        hadronGoodTracksFraction.push_back(*myGraph);
      }
    }


    /**
     * This convenience function computes the trigger performance maps
     * the trigger performance, so they are ready for a new round of
     * analysis.
     */  
    void Analyzer::fillTriggerPerformanceMaps(Tracker& tracker) {
      std::map<double, TH2D>& efficiencyMaps = myMapBag.getMaps(mapBag::efficiencyMap);
      std::map<double, TH2D>& thresholdMaps = myMapBag.getMaps(mapBag::thresholdMap);

      TH2D& thicknessMap = myMapBag.getMaps(mapBag::thicknessMap)[mapBag::dummyMomentum];
      TH2D& windowMap = myMapBag.getMaps(mapBag::windowMap)[mapBag::dummyMomentum];
      TH2D& suggestedSpacingMap = myMapBag.getMaps(mapBag::suggestedSpacingMap)[mapBag::dummyMomentum];
      TH2D& suggestedSpacingMapAW = myMapBag.getMaps(mapBag::suggestedSpacingMapAW)[mapBag::dummyMomentum];
      TH2D& nominalCutMap = myMapBag.getMaps(mapBag::nominalCutMap)[mapBag::dummyMomentum];

      LayerVector& layerSet = tracker.getLayers();
      LayerVector::iterator layIt;
      Layer* aLayer;
      ModuleVector::iterator modIt;
      ModuleVector* moduleSet;
      Module* aModule;
      double myValue;

      // Compute the trigger efficiency maps. The pT values are given by
      // the maps
      for (std::map<double, TH2D>::iterator it = efficiencyMaps.begin();
           it!=efficiencyMaps.end(); ++it) {
        double myPt = it->first;
        TH2D& myMap = it->second;

        // Reset the our map, in case it is not empty
        for (int i=1; i<=myMap.GetNbinsX(); ++i)
          for (int j=1; j<=myMap.GetNbinsY(); ++j)
            myMap.SetBinContent(i,j,0);

        // Create a map for the counter
        TH2D* counter = (TH2D*)myMap.Clone();

        // Loop over all the modules
        for(layIt = layerSet.begin(); layIt!= layerSet.end(); ++layIt) {
          aLayer = (*layIt);
          moduleSet = aLayer->getModuleVector();
          for(modIt = moduleSet->begin(); modIt != moduleSet->end(); ++modIt) {
            aModule = (*modIt);
            myValue = aModule->getTriggerProbability(myPt);

            if (myValue>=0) {
              // Draw the module
              XYZVector start = (aModule->getCorner(0)+aModule->getCorner(1))/2;
              XYZVector end = (aModule->getCorner(2)+aModule->getCorner(3))/2;
              XYZVector diff = end-start;
              XYZVector point;
              for (double l=0; l<=1; l+=0.1) {
                point = start + l * diff;
                myMap.Fill(point.Z(), point.Rho(), myValue);
                counter->Fill(point.Z(), point.Rho(), 1);
              }
            }
          }
        }

        // Normalize the counts to the number of hits per bin ...
        for (int i=1; i<=myMap.GetNbinsX(); ++i)
          for (int j=1; j<=myMap.GetNbinsY(); ++j)
            if (counter->GetBinContent(i,j)!=0)
              myMap.SetBinContent(i,j, myMap.GetBinContent(i,j) / counter->GetBinContent(i,j));
        // ... and get rid of the counter
        delete counter;
      }

      // Compute the trigger threshold maps. The efficiency values are
      // given by the maps
      for (std::map<double, TH2D>::iterator it = thresholdMaps.begin();
           it!=thresholdMaps.end(); ++it) {
        double myEfficiency = it->first;
        TH2D& myMap = it->second;

        // Reset the our map, in case it is not empty
        for (int i=1; i<=myMap.GetNbinsX(); ++i)
          for (int j=1; j<=myMap.GetNbinsY(); ++j)
            myMap.SetBinContent(i,j,0);

        // Create a map for the counter
        TH2D* counter = (TH2D*)myMap.Clone();

        // Loop over all the modules
        for(layIt = layerSet.begin(); layIt!= layerSet.end(); ++layIt) {
          aLayer = (*layIt);
          moduleSet = aLayer->getModuleVector();
          for(modIt = moduleSet->begin(); modIt != moduleSet->end(); ++modIt) {
            aModule = (*modIt);
            myValue = aModule->getPtThreshold(myEfficiency);

            if (myValue>=0) {
              // Draw the module
              XYZVector start = (aModule->getCorner(0)+aModule->getCorner(1))/2;
              XYZVector end = (aModule->getCorner(2)+aModule->getCorner(3))/2;
              XYZVector diff = end-start;
              XYZVector point;
              for (double l=0; l<1; l+=0.1) {
                point = start + l * diff;
                myMap.Fill(point.Z(), point.Rho(), myValue);
                counter->Fill(point.Z(), point.Rho(), 1);
              }
            }
          }
        }

        // Normalize the counts to the number of hits per bin ...
        for (int i=1; i<=myMap.GetNbinsX(); ++i)
          for (int j=1; j<=myMap.GetNbinsY(); ++j)
            if (counter->GetBinContent(i,j)!=0)
              myMap.SetBinContent(i,j, myMap.GetBinContent(i,j) / counter->GetBinContent(i,j));
        // ... and get rid of the counter
        delete counter;
      }

      // Then: single maps

      // Reset the our map, in case it is not empty
      //thicknessMap.Reset();
      //windowMap.Reset();
      //suggestedSpacingMap.Reset();
      //suggestedSpacingMapAW.Reset();
      for (int i=1; i<=thicknessMap.GetNbinsX(); ++i) {
        for (int j=1; j<=thicknessMap.GetNbinsY(); ++j) {
          thicknessMap.SetBinContent(i,j,0);
          windowMap.SetBinContent(i,j,0);
          suggestedSpacingMap.SetBinContent(i,j,0);
          suggestedSpacingMapAW.SetBinContent(i,j,0);
          nominalCutMap.SetBinContent(i,j,0);
        }
      }

      // Create a map for the counter
      TH2D* counter = (TH2D*)thicknessMap.Clone();
      TH2D* counterSpacing = (TH2D*)suggestedSpacingMap.Clone();
      TH2D* counterSpacingAW = (TH2D*)suggestedSpacingMapAW.Clone();

      double myThickness;
      double myWindow;
      double mySuggestedSpacing;
      double mySuggestedSpacingAW;
      //double myPitch;
      //double myWindowmm;
      // Loop over all the modules
      for(layIt = layerSet.begin(); layIt!= layerSet.end(); ++layIt) {
        aLayer = (*layIt);
        moduleSet = aLayer->getModuleVector();
        for(modIt = moduleSet->begin(); modIt != moduleSet->end(); ++modIt) {
          aModule = (*modIt);
          if (aModule->getReadoutType()!=Module::Pt) continue;
          myThickness = aModule->getStereoDistance();
          myWindow = aModule->getTriggerWindow();
          //myWindowmm = myWindow * (aModule->getLowPitch() + aModule->getHighPitch())/2.;
          mySuggestedSpacing = aModule->getOptimalSpacing(5); // TODO: put this 5 in a configuration of some sort
          mySuggestedSpacingAW = aModule->getOptimalSpacing(aModule->getTriggerWindow());
          double nominalCut = aModule->getPtCut();
          /*        XYZVector modCenter = aModule->getMeanPoint();
                    if (aModule->getSubdetectorType()==Module::Endcap) {
                    corrFactor = pow(modCenter.Rho(),2) / modCenter.Z();
                    } else {
                    corrFactor = modCenter.Rho();
                    }*/

          // Draw the module
          XYZVector start = (aModule->getCorner(0)+aModule->getCorner(1))/2;
          XYZVector end = (aModule->getCorner(2)+aModule->getCorner(3))/2;
          XYZVector diff = end-start;
          XYZVector point;
          double myZ, myRho;
          for (double l=0; l<1; l+=0.1) {
            point = start + l * diff;
            myZ=point.Z();
            myRho=point.Rho();
            thicknessMap.Fill(myZ, myRho, myThickness);
            windowMap.Fill(myZ, myRho, myWindow);
            nominalCutMap.Fill(myZ, myRho, nominalCut);
            counter->Fill(myZ, myRho, 1);
            if (mySuggestedSpacing!=0) {
              suggestedSpacingMap.Fill(myZ, myRho, mySuggestedSpacing);
              counterSpacing->Fill(myZ, myRho, 1);
            }
            if (mySuggestedSpacingAW!=0) {
              suggestedSpacingMapAW.Fill(myZ, myRho, mySuggestedSpacingAW);
              counterSpacingAW->Fill(myZ, myRho, 1);
            }
          }
        }
      }

      // Normalize the counts to the number of hits per bin ...
      for (int i=1; i<=thicknessMap.GetNbinsX(); ++i) {
        for (int j=1; j<=thicknessMap.GetNbinsY(); ++j) {
          if (counter->GetBinContent(i,j)!=0) {
            thicknessMap.SetBinContent(i,j, thicknessMap.GetBinContent(i,j) / counter->GetBinContent(i,j));
            windowMap.SetBinContent(i,j, windowMap.GetBinContent(i,j) / counter->GetBinContent(i,j));
            nominalCutMap.SetBinContent(i,j, nominalCutMap.GetBinContent(i,j) / counter->GetBinContent(i,j));
            if ((suggestedSpacingMap.GetBinContent(i,j)/counterSpacing->GetBinContent(i,j))>50) {
              std::cout << "debug: for bin " << i << ", " << j << " suggestedSpacing is " << suggestedSpacingMap.GetBinContent(i,j)
                << " and counter is " << counterSpacing->GetBinContent(i,j) << std::endl;
            }
            suggestedSpacingMap.SetBinContent(i,j, suggestedSpacingMap.GetBinContent(i,j) / counterSpacing->GetBinContent(i,j));
            suggestedSpacingMapAW.SetBinContent(i,j, suggestedSpacingMapAW.GetBinContent(i,j) / counterSpacingAW->GetBinContent(i,j));
          }
        }
      }
      // ... and get rid of the counter
      delete counter;
      delete counterSpacing;
      delete counterSpacingAW;
    }


    void Analyzer::fillPowerMap(Tracker& tracker) {
      TH2D& irradiatedPowerConsumptionMap = myMapBag.getMaps(mapBag::irradiatedPowerConsumptionMap)[mapBag::dummyMomentum];
      TH2D& totalPowerConsumptionMap = myMapBag.getMaps(mapBag::totalPowerConsumptionMap)[mapBag::dummyMomentum];

      LayerVector& layerSet = tracker.getLayers();
      LayerVector::iterator layIt;
      Layer* aLayer;
      ModuleVector::iterator modIt;
      ModuleVector* moduleSet;
      Module* aModule;
      string cntName;

      // Then: single maps

      // Reset the our map, in case it is not empty
      //thicknessMap.Reset();
      //windowMap.Reset();
      //suggestedSpacingMap.Reset();
      //suggestedSpacingMapAW.Reset();
      for (int i=1; i<=irradiatedPowerConsumptionMap.GetNbinsX(); ++i) {
        for (int j=1; j<=irradiatedPowerConsumptionMap.GetNbinsY(); ++j) {
          irradiatedPowerConsumptionMap.SetBinContent(i,j,0);
          totalPowerConsumptionMap.SetBinContent(i,j,0);
        }
      }

      // Create a map for the counter
      TH2D* counter = (TH2D*)irradiatedPowerConsumptionMap.Clone();

      // Loop over all the modules
      for(layIt = layerSet.begin(); layIt!= layerSet.end(); ++layIt) {
        aLayer = (*layIt);
        moduleSet = aLayer->getModuleVector();

        cntName = aLayer->getContainerName();
        if (cntName == "") continue;      

        for(modIt = moduleSet->begin(); modIt != moduleSet->end(); ++modIt) {
          aModule = (*modIt);

          if ((aModule->getMeanPoint().Z()<0) || (aModule->getMeanPoint().Phi()<0) || (aModule->getMeanPoint().Phi()>M_PI/2)) continue;
          double myPower = aModule->getIrradiatedPowerConsumption();
          ModuleType& myType = tracker.getModuleType( aModule->getType() );
          double myPowerChip = myType.getPower( aModule->getNChannels() );

          // Draw the module
          XYZVector start = (aModule->getCorner(0)+aModule->getCorner(1))/2;
          XYZVector end = (aModule->getCorner(2)+aModule->getCorner(3))/2;
          XYZVector diff = end-start;
          XYZVector point;
          double myZ, myRho;
          for (double l=0; l<1; l+=0.1) {
            point = start + l * diff;
            myZ=point.Z();
            myRho=point.Rho();
            irradiatedPowerConsumptionMap.Fill(myZ, myRho, myPower);
            totalPowerConsumptionMap.Fill(myZ, myRho, myPowerChip+myPower);
            counter->Fill(myZ, myRho, 1);
          }
        }
      }

      // Normalize the counts to the number of hits per bin ...
      for (int i=1; i<=irradiatedPowerConsumptionMap.GetNbinsX(); ++i) {
        for (int j=1; j<=irradiatedPowerConsumptionMap.GetNbinsY(); ++j) {
          if (counter->GetBinContent(i,j)!=0) {
            irradiatedPowerConsumptionMap.SetBinContent(i,j, irradiatedPowerConsumptionMap.GetBinContent(i,j) / counter->GetBinContent(i,j));
            totalPowerConsumptionMap.SetBinContent(i,j, totalPowerConsumptionMap.GetBinContent(i,j) / counter->GetBinContent(i,j));
          }
        }
      }
      // ... and get rid of the counter
      delete counter;

    }

    void Analyzer::prepareTrackerMap(TH2D& myMap, const std::string& name, const std::string& title) { 
      int mapBinsY = int( (outer_radius + volume_width) * 1.1 / 10.); // every cm
      int mapBinsX = int( (max_length) * 1.1 / 10.); // every cm
      myMap.SetName(name.c_str());
      myMap.SetTitle(title.c_str());
      myMap.SetXTitle("z [mm]");
      myMap.SetYTitle("r [mm]");
      myMap.SetBins(mapBinsX, 0.0, max_length*1.1, mapBinsY, 0.0, (outer_radius + volume_width) * 1.1);
      myMap.Reset();
    }

    void Analyzer::prepareRadialTrackerMap(TH2D& myMap, const std::string& name, const std::string& title) { 
      int mapBinsY = int( (2*outer_radius) * 1.1 / 10.); // every cm
      int mapBinsX = int( (2*outer_radius) * 1.1 / 10.); // every cm
      myMap.SetName(name.c_str());
      myMap.SetTitle(title.c_str());
      myMap.SetXTitle("x [mm]");
      myMap.SetYTitle("y [mm]");
      myMap.SetBins(mapBinsX, -outer_radius*1.1, outer_radius*1.1, mapBinsY, -outer_radius*1.1, outer_radius*1.1);
      myMap.Reset();
    }


    /**
     * This convenience function resets and empties all histograms for
     * the trigger performance, so they are ready for a new round of
     * analysis.
     * @parameter triggerMomenta the vector of pt to test the trigger
     */

    void Analyzer::prepareTriggerPerformanceHistograms(const int& nTracks, const double& myEtaMax, const std::vector<double>& triggerMomenta, const std::vector<double>& thresholdProbabilities) {
      // Clean-up and prepare the trigger performance maps
      myMapBag.clearMaps(mapBag::efficiencyMap);
      myMapBag.clearMaps(mapBag::thresholdMap);
      myMapBag.clearMaps(mapBag::thicknessMap);
      myMapBag.clearMaps(mapBag::windowMap);
      myMapBag.clearMaps(mapBag::suggestedSpacingMap);
      myMapBag.clearMaps(mapBag::suggestedSpacingMapAW);
      myMapBag.clearMaps(mapBag::nominalCutMap);

      std::map<double, TH2D>& thresholdMaps = myMapBag.getMaps(mapBag::thresholdMap);
      std::map<double, TH2D>& efficiencyMaps = myMapBag.getMaps(mapBag::efficiencyMap);

      // PT Threshold maps here
      ostringstream tempSS;
      // string tempString;
      for (std::vector<double>::const_iterator it = thresholdProbabilities.begin();
           it!=thresholdProbabilities.end(); ++it) {
        TH2D& myMap = thresholdMaps[(*it)];
        tempSS.str("");
        tempSS << "ptThresholdMap_" << std::dec << (*it) * 100 << "perc";
        prepareTrackerMap(myMap, tempSS.str(), tempSS.str());
      }

      // Efficiency maps here
      for (std::vector<double>::const_iterator it = triggerMomenta.begin();
           it!=triggerMomenta.end(); ++it) {
        TH2D& myMap = efficiencyMaps[(*it)];
        tempSS.str("");
        tempSS << "triggerEfficiencyMap_" << std::dec << (*it) << "GeV";
        prepareTrackerMap(myMap, tempSS.str(), tempSS.str());
      }

      // Single maps here
      TH2D& thicknessMap = myMapBag.getMaps(mapBag::thicknessMap)[mapBag::dummyMomentum];
      TH2D& windowMap = myMapBag.getMaps(mapBag::windowMap)[mapBag::dummyMomentum];
      TH2D& suggestedSpacingMap = myMapBag.getMaps(mapBag::suggestedSpacingMap)[mapBag::dummyMomentum];
      TH2D& suggestedSpacingMapAW = myMapBag.getMaps(mapBag::suggestedSpacingMapAW)[mapBag::dummyMomentum];
      TH2D& nominalCutMap = myMapBag.getMaps(mapBag::nominalCutMap)[mapBag::dummyMomentum];

      prepareTrackerMap(thicknessMap, "thicknessMap", "Module thickness map");
      //thicknessMap.SetMinimum(0);
      //thicknessMap.SetMaximum(6);
      prepareTrackerMap(windowMap, "windowMap", "Module window map");
      prepareTrackerMap(suggestedSpacingMap, "suggestedSpacingMap", "Map of nearest available spacing [using standard window]");
      prepareTrackerMap(suggestedSpacingMapAW, "suggestedSpacingMapAW", "Map of nearest available spacing [using selected windows]");
      prepareTrackerMap(nominalCutMap, "nominalCutMap", "Map of nominal pT cut");


      // Clear the graph (and profile) list
      myGraphBag.clearTriggerGraphs();
      myProfileBag.clearTriggerProfiles();

      // Prepare the graphs to record the number of triggered points
      // std::map<double, TGraph>& trigGraphs = myGraphBag.getGraphs(graphBag::TriggerGraph|graphBag::TriggeredGraph);
      std::map<double, TProfile>& trigProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggeredProfile);
      std::map<double, TProfile>& trigFractionProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggeredFractionProfile);
      std::map<double, TProfile>& trigPurityProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggerPurityProfile);

      // Prepare the graphs for the trigger performace
      std::ostringstream aName;

      // TODO: tune this parameter
      int nbins = int(nTracks/10.);
      for (vector<double>::const_iterator iter = triggerMomenta.begin();
           iter != triggerMomenta.end();
           ++iter) {
        //std::pair<double, TGraph> elemGraph;
        std::pair<double, TProfile> elemProfile;
        std::pair<double, TProfile> elemFractionProfile;
        std::pair<double, TProfile> elemPurityProfile;
        //TGraph graph;
        TProfile profile("dummyName", "dummyTitle", nbins, 0, myEtaMax);
        //elemGraph.first = *iter;
        //elemGraph.second = graph;
        elemProfile.first = *iter;
        elemProfile.second = profile;
        elemFractionProfile.first = *iter;
        elemFractionProfile.second = profile;
        elemPurityProfile.first = *iter;
        elemPurityProfile.second = profile;
        // Prepare plots: triggered graphs
        // trigGraphs.insert(elemGraph);
        trigProfiles.insert(elemProfile);
        trigFractionProfiles.insert(elemFractionProfile);
        trigPurityProfiles.insert(elemPurityProfile);

        // trigGraphs[*iter].SetTitle("Average triggered points;#eta;Triggered points <N>");
        trigProfiles[*iter].SetTitle("Average triggered points;#eta;Triggered points <N>");
        trigFractionProfiles[*iter].SetTitle("Average trigger efficiency;#eta;Efficiency [%]");
        trigPurityProfiles[*iter].SetTitle("Average stub purity;#eta;Purity [%]");
        aName.str(""); aName << "triggered_vs_eta" << *iter << "_profile";      
        trigProfiles[*iter].SetName(aName.str().c_str());
        aName.str(""); aName << "triggered_vs_eta" << *iter << "_fractionProfile";
        trigFractionProfiles[*iter].SetName(aName.str().c_str());
        aName.str(""); aName << "triggered_vs_eta" << *iter << "_purityProfile";
        trigPurityProfiles[*iter].SetName(aName.str().c_str());
      }

      //std::pair<double, TGraph> elemTotalGraph;
      std::pair<double, TProfile> elemTotalProfile;
      //TGraph totalGraph;
      char dummyName[256]; sprintf(dummyName, "dummyName%d", bsCounter++);
      TProfile totalProfile(dummyName, dummyName, nbins, 0, myEtaMax);
      //elemTotalGraph.first = graphBag::Triggerable;
      //elemTotalGraph.second = totalGraph;
      elemTotalProfile.first = profileBag::Triggerable;
      elemTotalProfile.second = totalProfile;
      // Prepare plot: total trigger points
      //trigGraphs.insert(elemTotalGraph);
      trigProfiles.insert(elemTotalProfile);
      //trigGraphs[graphBag::Triggerable].SetTitle("Average triggered points;#eta;Triggered points <N>");
      trigProfiles[profileBag::Triggerable].SetTitle("Average triggered points;#eta;Triggered points <N>");
      //aName.str(""); aName << "triggerable_vs_eta_graph";
      //trigGraphs[graphBag::Triggerable].SetName(aName.str().c_str());
      aName.str(""); aName << "triggerable_vs_eta_profile";
      trigProfiles[profileBag::Triggerable].SetName(aName.str().c_str());

    }


    void Analyzer::preparePowerHistograms() {
      myMapBag.clearMaps(mapBag::irradiatedPowerConsumptionMap);
      myMapBag.clearMaps(mapBag::totalPowerConsumptionMap);
      TH2D& irradiatedPowerConsumptionMap = myMapBag.getMaps(mapBag::irradiatedPowerConsumptionMap)[mapBag::dummyMomentum]; // dummyMomentum is supplied because it is a single map. Multiple maps are indexed like arrays (see above efficiency maps)
      TH2D& totalPowerConsumptionMap = myMapBag.getMaps(mapBag::totalPowerConsumptionMap)[mapBag::dummyMomentum]; // dummyMomentum is supplied because it is a single map. Multiple maps are indexed like arrays (see above efficiency maps)
      prepareTrackerMap(irradiatedPowerConsumptionMap, "irradiatedPowerConsumptionMap", "Map of power dissipation in sensors (after irradiation)");
      prepareTrackerMap(totalPowerConsumptionMap, "irradiatedPowerConsumptionMap", "Map of power dissipation in modules (after irradiation)");
    }

    void Analyzer::prepareTriggerProcessorHistograms() {
      myMapBag.clearMaps(mapBag::moduleConnectionEtaMap);
      TH2D& moduleConnectionEtaMap = myMapBag.getMaps(mapBag::moduleConnectionEtaMap)[mapBag::dummyMomentum];
      prepareTrackerMap(moduleConnectionEtaMap, "moduleConnectionEtaMap", "Map");

      myMapBag.clearMaps(mapBag::moduleConnectionPhiMap);
      TH2D& moduleConnectionPhiMap = myMapBag.getMaps(mapBag::moduleConnectionPhiMap)[mapBag::dummyMomentum];
      prepareRadialTrackerMap(moduleConnectionPhiMap, "moduleConnectionPhiMap", "Map");

      myMapBag.clearMaps(mapBag::moduleConnectionEndcapPhiMap);
      TH2D& moduleConnectionEndcapPhiMap = myMapBag.getMaps(mapBag::moduleConnectionEndcapPhiMap)[mapBag::dummyMomentum];
      prepareRadialTrackerMap(moduleConnectionEndcapPhiMap, "moduleConnectionEndcapPhiMap", "Map");

      for (int i=1; i<=moduleConnectionEtaMap.GetNbinsX(); ++i) {
        for (int j=1; j<=moduleConnectionEtaMap.GetNbinsY(); ++j) {
          moduleConnectionEtaMap.SetBinContent(i,j,0);
        }
      }

      for (int i=1; i<=moduleConnectionPhiMap.GetNbinsX(); ++i) {
        for (int j=1; j<=moduleConnectionPhiMap.GetNbinsY(); ++j) {
          moduleConnectionPhiMap.SetBinContent(i,j,0);
        }
      }

      for (int i=1; i<=moduleConnectionEndcapPhiMap.GetNbinsX(); ++i) {
        for (int j=1; j<=moduleConnectionEndcapPhiMap.GetNbinsY(); ++j) {
          moduleConnectionEndcapPhiMap.SetBinContent(i,j,0);
        }
      }
    }

    /**
     * This convenience function resets and empties all histograms for the
     * geometry so they are ready for a new round of analysis.
     */
    void Analyzer::clearGeometryHistograms() {
      // geometry analysis
      mapPhiEta.Reset();
      mapPhiEta.SetNameTitle("mapPhiEta", "Number of hits;#phi;#eta");
      etaProfileCanvas.SetName("etaProfileCanvas"); etaProfileCanvas.SetTitle("Eta Profiles");
      hitDistribution.Reset();
      hitDistribution.SetNameTitle("hitDistribution", "Hit distribution");
      //geomLite->SetName("geometryLite");   geomLite->SetTitle("Modules geometry");
      //geomLiteXY->SetName("geometryLiteXY"); geomLiteXY->SetTitle("Modules geometry (XY Section)");
      //geomLiteYZ->SetName("geometryLiteYZ"); geomLiteYZ->SetTitle("Modules geometry (EC Section)");
      //geomLiteEC->SetName("geometryLiteEC"); geomLiteEC->SetTitle("Modules geometry (Endcap)");

      // Power density
      while (powerDensity.GetN()) powerDensity.RemovePoint(0);
      powerDensity.SetName("powerdensity");
      powerDensity.SetTitle("Power density;Total area [m^{2}];Power density [kW / m^{2}]");
    }

    /**
     * This convenience function sets all values in the internal array <i>cells</i> to zero.
     */
    void Analyzer::clearCells() {
      for (unsigned int i = 0; i < cells.size(); i++) {
        cells.at(i).clear();
      }
      cells.clear();
    }

    /**
     * This convenience function sets the number of bins and the lower and upper range for their contents for
     * each of the available histograms.
     * @param bins The number of bins in each 1D histogram
     * @param min The minimal eta value that should be plotted
     * @param max the maximal eta value that should be plotted
     */
    void Analyzer::setHistogramBinsBoundaries(int bins, double min, double max) {
      // single category
      ractivebarrel.SetBins(bins, min, max);
      ractiveendcap.SetBins(bins, min, max);
      rserfbarrel.SetBins(bins, min, max);
      rserfendcap.SetBins(bins, min, max);
      rlazybarrel.SetBins(bins, min, max);
      rlazyendcap.SetBins(bins, min, max);
      rlazytube.SetBins(bins, min, max);
      rlazyuserdef.SetBins(bins, min, max);
      iactivebarrel.SetBins(bins, min, max);
      iactiveendcap.SetBins(bins, min, max);
      iserfbarrel.SetBins(bins, min, max);
      iserfendcap.SetBins(bins, min, max);
      ilazybarrel.SetBins(bins, min, max);
      ilazyendcap.SetBins(bins, min, max);
      ilazytube.SetBins(bins, min, max);
      ilazyuserdef.SetBins(bins, min, max);
      // composite
      rbarrelall.SetBins(bins, min, max);
      rendcapall.SetBins(bins, min, max);
      ractiveall.SetBins(bins, min, max);
      rserfall.SetBins(bins, min, max);
      rlazyall.SetBins(bins, min, max);
      ibarrelall.SetBins(bins, min, max);
      iendcapall.SetBins(bins, min, max);
      iactiveall.SetBins(bins, min, max);
      iserfall.SetBins(bins, min, max);
      ilazyall.SetBins(bins, min, max);
      // outside tracking volume
      rextraservices.SetBins(bins, min, max);
      rextrasupports.SetBins(bins, min, max);
      iextraservices.SetBins(bins, min, max);
      iextrasupports.SetBins(bins, min, max);
      // global
      rglobal.SetBins(bins, min, max);
      iglobal.SetBins(bins, min, max);
      // isolines
      isor.SetBins(bins, 0.0, max_length, bins / 2, 0.0, outer_radius + volume_width);
      isoi.SetBins(bins, 0.0, max_length, bins / 2, 0.0, outer_radius + volume_width);
      // Material distribution maps
      int materialMapBinsY = int( (outer_radius + volume_width) * 1.1 / 5.); // every half a cm
      int materialMapBinsX = int( (max_length) * 1.1 / 5.); // every half a cm
      mapRadiation.SetBins(materialMapBinsX, 0.0, max_length*1.1, materialMapBinsY, 0.0, (outer_radius + volume_width) * 1.1);
      mapInteraction.SetBins(materialMapBinsX, 0.0, max_length*1.1, materialMapBinsY, 0.0, (outer_radius + volume_width) * 1.1);
      mapRadiationCount.SetBins(materialMapBinsX, 0.0, max_length*1.1, materialMapBinsY, 0.0, (outer_radius + volume_width) * 1.1);
      mapInteractionCount.SetBins(materialMapBinsX, 0.0, max_length*1.1, materialMapBinsY, 0.0, (outer_radius + volume_width) * 1.1);
      mapRadiationCalib.SetBins(materialMapBinsX, 0.0, max_length*1.1, materialMapBinsY, 0.0, (outer_radius + volume_width) * 1.1);
      mapInteractionCalib.SetBins(materialMapBinsX, 0.0, max_length*1.1, materialMapBinsY, 0.0, (outer_radius + volume_width) * 1.1);
    }

    /**
     * This convenience function sets the number of bins and the lower and upper range for their contents for
     * each of the cells that make up the basis for the 2D histograms.
     * @param bins The number of bins in eta; the number of bins in r will be half that
     * @param minr The minimum radius that will be considered for tracking
     * @param maxr The maximum radius that will be considered for tracking
     * @param mineta The minimum value of eta that will be considered for tracking
     * @param maxeta The maximum value of eta that will be considered for tracking
     */
    void Analyzer::setCellBoundaries(int bins, double minr, double maxr, double mineta, double maxeta) {
      double rstep, etastep;
      rstep = 2 * (maxr - minr) / bins;
      etastep = (maxeta - mineta) / bins;
      Cell c;
      c.rmin = 0.0;  // TODO: is this right?
      c.rmax = 0.0;
      c.etamin = 0.0;
      c.etamax = 0.0;
      c.rlength = 0.0;
      c.ilength = 0.0;
      cells.resize(bins);
      for (unsigned int i = 0; i < cells.size(); i++) if (bins != 1) cells.at(i).resize(bins / 2, c);
      for (unsigned int i = 0; i < cells.size(); i++) {
        for (unsigned int j = 0; j < cells.at(i).size(); j++) {
          cells.at(i).at(j).rmin = minr + j * rstep;
          cells.at(i).at(j).rmax = minr + (j+1) * rstep;
          cells.at(i).at(j).etamin = mineta + i * etastep;;
          cells.at(i).at(j).etamax = mineta + (i+1) * etastep;
        }
      }
    }


    /**
     * Fills the material distribution maps
     * @param r The radius at which the hit was detected
     * @param theta The angle of the track used for meterial detection
     * @param rl The local radiation length
     * @param il The local interaction length
     */
    void Analyzer::fillMapRT(const double& r, const double& theta, const Material& mat) {
      double z = r /tan(theta);
      if (mat.radiation>0){
        mapRadiation.Fill(z,r,mat.radiation);
        mapRadiationCount.Fill(z,r);
      } 
      if (mat.interaction>0) {
        mapInteraction.Fill(z,r,mat.interaction);
        mapInteractionCount.Fill(z,r);
      }
    }

    /**
     * Fills the material distribution maps
     * @param r The radius at which the hit was detected
     * @param z The z coordinate of the hit
     * @param rl The local radiation length
     * @param il The local interaction length
     */
    void Analyzer::fillMapRZ(const double& r, const double& z, const Material& mat) {
      if (mat.radiation>0){
        mapRadiation.Fill(z,r,mat.radiation);
        mapRadiationCount.Fill(z,r);
      } 
      if (mat.interaction>0) {
        mapInteraction.Fill(z,r,mat.interaction);
        mapInteractionCount.Fill(z,r);
      }
    }

    /**
     * @return a (hit-scaled) map of radiation length
     */
    TH2D& Analyzer::getHistoMapRadiation() {
      int nBins = mapRadiation.GetNbinsX()*mapRadiation.GetNbinsY();
      double content;
      double count;
      for (int iBin=1; iBin<=nBins; ++iBin) {
        content = mapRadiation.GetBinContent(iBin);
        count = mapRadiationCount.GetBinContent(iBin);
        //mapRadiationCalib.SetBinContent(iBin,content);
        if (count==1) mapRadiationCalib.SetBinContent(iBin,content);
        else if (count>1) mapRadiationCalib.SetBinContent(iBin,content/double(count));
        //else if (count==0) mapRadiationCalib.SetBinContent(iBin, 0.);
      }
      return mapRadiationCalib;
    }

    /**
     * @return a (hit-scaled) map of interaction length
     */
    TH2D& Analyzer::getHistoMapInteraction() {
      int nBins = mapInteraction.GetNbinsX()*mapInteraction.GetNbinsY();
      double content;
      double count;
      for (int iBin=1; iBin<=nBins; ++iBin) {
        content = mapInteraction.GetBinContent(iBin);
        count = mapInteractionCount.GetBinContent(iBin);
        //mapInteractionCalib.SetBinContent(iBin,content);
        if (count==1) mapInteractionCalib.SetBinContent(iBin,content);
        else if (count>1) mapInteractionCalib.SetBinContent(iBin,content/double(count));
        //else if (count==0) mapInteractionCalib.SetBinContent(iBin, 0.);
      }
      return mapInteractionCalib;
    }


    /**
     * This function assigns the local radiation and interaction lengths of a detected hit to their position in the
     * (eta, r) space.
     * @param r The radius at which the hit was detected
     * @param eta The eta value of the current track
     * @param rl The local radiation length
     * @param il The local interaction length
     */
    void Analyzer::fillCell(double r, double eta, double theta, Material mat) {
      double rl = mat.radiation;
      double il = mat.interaction;
      int rindex, etaindex;
      if (cells.size() > 0) {
        for (rindex = 0; (unsigned int) rindex < cells.at(0).size(); rindex++) {
          if ((cells.at(0).at(rindex).rmin <= r) && (cells.at(0).at(rindex).rmax > r)) break;
        }
        if ((unsigned int) rindex < cells.at(0).size()) {
          for (etaindex = 0; (unsigned int) etaindex < cells.size(); etaindex++) {
            if ((cells.at(etaindex).at(rindex).etamin <= eta) && (cells.at(etaindex).at(rindex).etamax > eta)) break;
          }
          if ((unsigned int) etaindex < cells.size()) {
            cells.at(etaindex).at(rindex).rlength = cells.at(etaindex).at(rindex).rlength + rl;
            cells.at(etaindex).at(rindex).ilength = cells.at(etaindex).at(rindex).ilength + il;
          }
        }
      }
    }

    /**
     * The integrated radiation and interaction lengths in (eta, r)
     * are converted to (z, r) coordinates and stored in <i>isor</i>
     * and <i>isoi</i> in this function.
     */
    void Analyzer::transformEtaToZ() {
      int size_z, size_r, rindex, etaindex;
      double z, r, eta, z_max, z_min, r_max, r_min, z_c, r_c;
      // init: sizes and boundaries
      size_z = isor.GetNbinsX();
      size_r = isor.GetNbinsY();
      z_max = isor.GetXaxis()->GetXmax();
      z_min = isor.GetXaxis()->GetXmin();
      r_max = isor.GetYaxis()->GetXmax();
      r_min = isor.GetYaxis()->GetXmin();
      // new number of bins in z and r
      z_c = (z_max - z_min) / (2 * size_z);
      r_c = (r_max - r_min) / (2 * size_r);
      // radiation length loop
      for (int i = 1; i <= size_z; i++) {
        // calculate current z bin
        z = z_min + 2 * (i - 1) * z_c + z_c;
        for (int j = 1; j <= size_r; j++) {
          // calculate current r bin
          r = r_min + 2 * (j - 1) * r_c + r_c;
          eta = -log(tan(atan(r / z) / 2.0));
          // find corresponding r and eta positions
          etaindex = findCellIndexEta(eta);
          rindex = findCellIndexR(r);
          // fill in radiation length in r and z
          if ((etaindex >= 0) && (rindex >= 0)) isor.Fill(z, r, cells.at(etaindex).at(rindex).rlength);
        }
      }
      // init: sizes and boundaries
      size_z = isoi.GetNbinsX();
      size_r = isoi.GetNbinsY();
      z_max = isoi.GetXaxis()->GetXmax();
      z_min = isoi.GetXaxis()->GetXmin();
      r_max = isoi.GetYaxis()->GetXmax();
      r_min = isoi.GetYaxis()->GetXmin();
      // new number of bins in z and r
      z_c = (z_max - z_min) / (2 * size_z);
      r_c = (r_max - r_min) / (2 * size_r);
      // interaction length loop
      for (int i = 0; i < size_z; i++) {
        // calculate current z bin
        z = z_min + 2 * i * z_c + z_c;
        for (int j = 0; j < size_r; j++) {
          // calculate current r bin
          r = r_min + 2 * j * r_c + r_c;
          eta = -log(tan(atan(r / z) / 2.0));
          // find corresponding r and eta positions
          etaindex = findCellIndexEta(eta);
          rindex = findCellIndexR(r);
          // fill in interaction length in r and z
          if ((etaindex >= 0) && (rindex >= 0)) isoi.Fill(z, r, cells.at(etaindex).at(rindex).ilength);
        }
      }
    }

    // private
    /**
     * This is a convenience function that converts a radius to an index into the second dimension of <i>cells</i>.
     * @param r The given radius value
     * @return The corresponding index into the vector
     */
    int Analyzer::findCellIndexR(double r) {
      int index = -1;
      if (r >= 0) {
        if (cells.size() > 0) {
          index = 0;
          while (index < (int)cells.at(0).size()) {
            if (cells.at(0).at(index).rmax < r) index++;
            else break;
          }
          if (index == (int)cells.at(0).size()) index = -1;
          return index;
        }
      }
      return index;
    }

    /**
     * This is a convenience function that converts an eta value to an index into the first dimension of <i>cells</i>.
     * @param eta The given eta value
     * @return the corresponding index into the vector
     */
    int Analyzer::findCellIndexEta(double eta) {
      int index = -1;
      if (eta >= 0) {
        index = 0;
        while (index < (int)cells.size()) {
          if (cells.at(index).size() > 0) {
            if (cells.at(index).at(0).etamax < eta) index++;
            else break;
          }
          else {
            index = -1;
            break;
          }
        }
        if (index == (int)cells.size()) index = -1;
        return index;
      }
      return index;
    }

    // public
    /**
     * Creates the histograms to analyze the tracker coverage
     * @param tracker the tracker to be analyzed
     * @param nTracker the number of tracks to be used to analyze the coverage (defaults to 1000)
     */
    void Analyzer::analyzeGeometry(Tracker& tracker, int nTracks /*=1000*/ ) {
      geometryTracksUsed = nTracks;
      savingGeometryV.clear();
      clearGeometryHistograms();


      // A bunch of pointers
      std::map <std::string, int> moduleTypeCount;
      std::map<PosRef, std::string> layerNames;
      std::map <std::string, TProfile> etaProfileByType;
      //TH2D* aPlot;
      std::string aType;


      // Optimize the track creation on the real tracker
      std::pair <double, double> etaMinMax = tracker.getEtaMinMax();
      double absMinEta = fabs(etaMinMax.first);
      double absMaxEta = fabs(etaMinMax.second);
      double maxEta = (absMinEta>absMaxEta) ? absMinEta : absMaxEta;

      // Computing the margin of the tracks to shoot
      double randomPercentMargin = 0.04;
      double randomSpan = (etaMinMax.second - etaMinMax.first)*(1. + randomPercentMargin);
      double randomBase = etaMinMax.first - (etaMinMax.second - etaMinMax.first)*(randomPercentMargin)/2.;

      maxEta *= (1 + randomPercentMargin);
      if (maxEta<3) maxEta=3; // TODO: make this configurable

      // Initialize random number generator, counters and histograms
      myDice.SetSeed(MY_RANDOM_SEED);
      createResetCounters(tracker, moduleTypeCount);
      createLayerNames(tracker, layerNames);
      

      // Creating the layer hit coverage profiles
      layerEtaCoverageProfile.clear();
      for (std::map<PosRef, std::string>::iterator it = layerNames.begin();
           it!=layerNames.end(); ++it) {
        TProfile* aProfile = new TProfile(Form("layerEtaCoverageProfile%s", (it->second).c_str()), (it->second).c_str(), 200, maxEta, maxEta);
        layerEtaCoverageProfile[it->second] = (*aProfile);
        delete aProfile;
      }

      etaProfileByType.clear();

      for (std::map <std::string, int>::iterator it = moduleTypeCount.begin();
           it!=moduleTypeCount.end(); it++) {
        TProfile& aProfile = etaProfileByType[(*it).first];
        aProfile.SetBins(100, 0, maxEta);
        aProfile.SetName((*it).first.c_str());
        aProfile.SetTitle((*it).first.c_str());
      }

      LayerVector::iterator layIt;
      ModuleVector* moduleV = NULL;
      ModuleVector::iterator modIt;
      ModuleVector allModules;
      LayerVector& layerSet = tracker.getLayers();
      double zErrorCollider = tracker.getZErrorCollider();

      // Build the proper list of modules
      for (layIt=layerSet.begin(); layIt!=layerSet.end(); layIt++) {
        moduleV = (*layIt)->getModuleVector();
        if (!moduleV) {
          std::cerr << "ERROR in Analyzer::analyzeGeometry: cannot retrieve moduleVector from the layer" << std::endl;
          return;
        }
        for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
          // I pre-compute the boxes to reduce the calculations
          (*modIt)->computeBoundaries(zErrorCollider);
          allModules.push_back(*modIt);
        }
      }

      // The real simulation
      std::pair <XYZVector, double> aLine;
      ModuleVector hitModules;


      int nTrackHits;
      int nTracksPerSide = int(pow(nTracks, 0.5));
      int nBlocks = int(nTracksPerSide/2.);
      nTracks = nTracksPerSide*nTracksPerSide;
      mapPhiEta.SetBins(nBlocks, -1*M_PI, M_PI, nBlocks, -maxEta, maxEta);
      TH2I mapPhiEtaCount("mapPhiEtaCount ", "phi Eta hit count", nBlocks, -1*M_PI, M_PI, nBlocks, -maxEta, maxEta);
      totalEtaProfile.Reset();
      totalEtaProfile.SetName("totalEtaProfile");
      totalEtaProfile.SetMarkerStyle(8);
      totalEtaProfile.SetMarkerColor(1);
      totalEtaProfile.SetMarkerSize(1.5);
      totalEtaProfile.SetTitle("Number of hit modules;#eta;Number of hits");
      totalEtaProfile.SetBins(100, 0, maxEta);

      double layerHit;

      // Shoot nTracksPerSide^2 tracks
      for (int i=0; i<nTracksPerSide; i++) {
        for (int j=0; j<nTracksPerSide; j++) {
          // Reset the hit counter
          nTrackHits=0;
          // Generate a straight track and collect the list of hit modules
          aLine = shootDirection(randomBase, randomSpan);
          hitModules = trackHit( XYZVector(0, 0, myDice.Gaus(0, zErrorCollider)), aLine.first, &allModules);
          // Reset the per-type hit counter and fill it
          resetTypeCounter(moduleTypeCount);
          for (ModuleVector::iterator it = hitModules.begin(); it!=hitModules.end(); it++) {
            moduleTypeCount[(*it)->getType()]++;
            nTrackHits++;
          }
          // Fill the module type hit plot
          for (std::map <std::string, int>::iterator it = moduleTypeCount.begin(); it!=moduleTypeCount.end(); it++) {
            etaProfileByType[(*it).first].Fill(fabs(aLine.second), (*it).second);
          }
          // Fill other plots
          totalEtaProfile.Fill(fabs(aLine.second), hitModules.size());                // Total number of hits
          mapPhiEta.Fill(aLine.first.Phi(), aLine.second, hitModules.size()); // phi, eta 2d plot
          mapPhiEtaCount.Fill(aLine.first.Phi(), aLine.second);               // Number of shot tracks

          for (std::map<PosRef, std::string>::iterator it = layerNames.begin(); it!=layerNames.end(); ++it) {
            layerHit = 0;
            const PosRef& aLayerPosRef = it->first;
            for (ModuleVector::iterator moduleIt = hitModules.begin(); moduleIt!=hitModules.end(); moduleIt++) {
              if ((*moduleIt)->getLayerPositionalReference()==aLayerPosRef) {
                layerHit=1;
                continue;
              }
            }
            layerEtaCoverageProfile[it->second].Fill(aLine.second, layerHit);

          }

        }
      }
      
      // Create and archive for saving our 2D map of hits
      double hitCount;
      double trackCount;
      for (int nx=0; nx<=mapPhiEtaCount.GetNbinsX()+1; nx++) {
        for (int ny=0; ny<=mapPhiEtaCount.GetNbinsY()+1; ny++) {
          trackCount=mapPhiEtaCount.GetBinContent(nx, ny);
          if (trackCount>0) {
            hitCount=mapPhiEta.GetBinContent(nx, ny);
            mapPhiEta.SetBinContent(nx, ny, hitCount/trackCount);
            //hitCount=mapPhiEta.GetBinContent(nx, ny); //debug
            //std::cerr << "hitCount " << hitCount << std::endl; //debug
          }
        }
      }

      savingGeometryV.push_back(mapPhiEta);

      // Eta profile compute
      //TProfile *myProfile;

      etaProfileCanvas.cd();
      savingGeometryV.push_back(etaProfileCanvas);
      int plotCount=0;

      //TProfile* total = total2D.ProfileX("etaProfileTotal");
      char profileName_[256];
      sprintf(profileName_, "etaProfileTotal%d", bsCounter++);
      // totalEtaProfile = TProfile(*total2D.ProfileX(profileName_));
      savingGeometryV.push_back(totalEtaProfile);
      if (totalEtaProfile.GetMaximum()<maximum_n_planes) totalEtaProfile.SetMaximum(maximum_n_planes);
      totalEtaProfile.Draw();
      std::string profileName;
      TProfile* myProfile;
      for (std::map <std::string, TProfile>::iterator it = etaProfileByType.begin();
           it!=etaProfileByType.end(); it++) {
        plotCount++;
        myProfile=(TProfile*)it->second.Clone();
        savingGeometryV.push_back(*myProfile); // TODO: remove savingGeometryV everywhere :-) [VERY obsolete...]
        myProfile->SetMarkerStyle(8);
        myProfile->SetMarkerColor(Palette::color((*it).first));
        myProfile->SetMarkerSize(1);
        profileName = "etaProfile"+(*it).first;
        myProfile->SetName(profileName.c_str());
        myProfile->SetTitle((*it).first.c_str());
        myProfile->GetXaxis()->SetTitle("eta");
        myProfile->GetYaxis()->SetTitle("Number of hits");
        myProfile->Draw("same");
        typeEtaProfile.push_back(*myProfile);
      }

      // Record the fraction of hits per module
      hitDistribution.SetBins(nTracks, 0 , 1);
      savingGeometryV.push_back(hitDistribution);
      for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
        hitDistribution.Fill((*modIt)->getNHits()/double(nTracks));
      }


      std::map<std::string, ModuleType>& typeMap = tracker.getTypes();
      std::map<std::string, int> typeToCount;
      std::map<std::string, double> typeToPower;
      std::map<std::string, double> typeToSurface;
      std::string aSensorType;
      ModuleType aModuleType;

      for (ModuleVector::iterator aModule = allModules.begin();
           aModule != allModules.end(); ++aModule) {
        aSensorType = (*aModule)->getType();
        aModuleType = typeMap[aSensorType];
        typeToCount[aSensorType] ++;
        typeToSurface[aSensorType] += (*aModule)->getArea() / 1e6; // in mq
        typeToPower[aSensorType] += aModuleType.getPower( (*aModule)->getNChannels() ) / 1e3; // in kW
      }

      int iPoints = 0;
      for (std::map<std::string, int>::iterator typesIt = typeToCount.begin();
           typesIt != typeToCount.end(); ++typesIt ) {
        aSensorType = typesIt->first;
        powerDensity.SetPoint(iPoints++,
                              typeToSurface[aSensorType],
                              typeToPower[aSensorType] / typeToSurface[aSensorType] );
      }
      return;
    }

    // public
    // TODO!!!
    // Creates the geometry objects geomLite
    // (now they are canvases, but in the future
    // they could be just lists of 3d poly objects)
    // Moreover now we create them by calling the tracker object, which seems improper
    void Analyzer::createGeometryLite(Tracker& tracker) {
      if (!(geomLiteCreated &&
            geomLiteXYCreated &&
            geomLiteYZCreated &&
            geomLiteECCreated)) {
        tracker.createGeometry(true);
        geomLiteCreated=true;
        geomLiteXYCreated=true;
        geomLiteYZCreated=true;
        geomLiteECCreated=true;
      }
      geomLite = tracker.getGeomLite();
      geomLiteXY = tracker.getGeomLiteXY();
      geomLiteYZ = tracker.getGeomLiteYZ();
      geomLiteEC = tracker.getGeomLiteEC();
    }

    // private
    /**
     * Creates a module type map
     * It sets a different integer for each one
     * @param tracker the tracker to be analyzed
     * @param moduleTypeCount the map to count the different module types
     * @return the total number of module types
     */
    int Analyzer::createResetCounters(Tracker& tracker, std::map <std::string, int> &moduleTypeCount) {
      ModuleVector result;
      LayerVector::iterator layIt;
      ModuleVector* moduleV;
      ModuleVector::iterator modIt;

      std::string aType;
      int typeCounter=0;

      LayerVector& layerSet = tracker.getLayers();
      for (layIt=layerSet.begin(); layIt!=layerSet.end(); layIt++) {
        moduleV = (*layIt)->getModuleVector();
        for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
          aType = (*modIt)->getType();
          (*modIt)->resetNHits();
          if (moduleTypeCount.find(aType)==moduleTypeCount.end()) {
            moduleTypeCount[aType]=typeCounter++;
          }
        }
      }

      return(typeCounter);
    }

  // private
  /**
   * Creates a map of layer positional references to their readable names (string)
   */
  int Analyzer::createLayerNames(Tracker& tracker, std::map<PosRef, std::string>& layerNames ) {
    layerNames.clear();
    LayerVector::iterator layIt;
    ModuleVector* moduleV;
    ModuleVector::iterator modIt;
    Module* aModule;

    std::string aLayerName;
    std::ostringstream oss;
    PosRef aModuleLayerPosRef;

    LayerVector& layerSet = tracker.getLayers();
    for (layIt=layerSet.begin(); layIt!=layerSet.end(); layIt++) {
      moduleV = (*layIt)->getModuleVector();
      for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
        aModule = (*modIt);
        // if (!aModule->getReadoutType()) continue; // TODO BAD: this is just a patch to a misbehaviour!
        // XYZVector cen = aModule->getMeanPoint();
        // std::cerr << aModule->getReadoutType() << ", " << cen.Rho() << ", " << cen.Z() << ", " << cen.Phi() << ", \"" << aModule->getContainerName() << "\"" << std::endl;
        aModuleLayerPosRef = aModule->getLayerPositionalReference();
        if (layerNames.find(aModuleLayerPosRef)==layerNames.end()) {
          oss.str("");
          oss << aModule->getContainerName() << " " << aModule->getLayer();
          layerNames[aModuleLayerPosRef]=oss.str();
          // std::cerr << "PosRef = " << aModuleLayerPosRef << ", name = " << oss.str() << std::endl;
        }
      }
    }
    
    return layerNames.size();
  }

    // private
    /**
     * Shoots directions with random (flat) phi, random (flat) pseudorapidity
     * gives also the direction's eta
     * @param minEta minimum eta to shoot tracks
     * @param spanEta difference between minimum and maximum eta
     * @return the pair of value: pointing XYZVector and eta of the track
     */
    std::pair <XYZVector, double > Analyzer::shootDirection(double minEta, double spanEta) {
      std::pair <XYZVector, double> result;

      double eta;
      double phi;
      double theta;

      // phi is random [0, 2pi)
      phi = myDice.Rndm() * 2 * M_PI; // debug

      // eta is random (-4, 4]
      eta = myDice.Rndm() * spanEta + minEta;
      theta=2*atan(exp(-1*eta));

      // Direction
      result.first  = XYZVector(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
      result.second = eta;
      return result;
    }

    // private
    /**
     * Checks whether a track would hit a module
     * @param origin XYZVector of origin of the track
     * @param direction pointing XYZVector of the track
     * @param moduleV vector of modules to be checked
     * @return the vector of hit modules
     */
    ModuleVector Analyzer::trackHit(const XYZVector& origin, const XYZVector& direction, ModuleVector* moduleV) {
      ModuleVector result;
      ModuleVector::iterator modIt;
      double distance;

      for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
        // A module can be hit if it fits the phi (precise) contraints
        // and the eta constaints (taken assuming origin within 5 sigma)
        if ((*modIt)->couldHit(direction.Eta(), direction.Phi())) {
          distance=(*modIt)->trackCross(origin, direction);
          if (distance>0) {
            result.push_back(*modIt);
          }
        }
      }
      return result;
    }

    // Resets a module type counter
    void Analyzer::resetTypeCounter(std::map <std::string, int> &modTypes) {
      for (std::map <std::string, int>::iterator it = modTypes.begin();
           it!=modTypes.end(); it++) {
        (*it).second = 0;
      }
    }

    double Analyzer::diffclock(clock_t clock1, clock_t clock2) {
      double diffticks=clock1-clock2;
      double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
      return diffms;
    }

    std::vector<TObject> Analyzer::getSavingVector() {
      std::vector<TObject> result;
      std::vector<TObject>::iterator it;

      for (it=savingGeometryV.begin(); it!=savingGeometryV.end(); ++it) {
        result.push_back(*it);
      }
      for (it=savingMaterialV.begin(); it!=savingMaterialV.end(); ++it) {
        result.push_back(*it);
      }
      return result;

    }


    void Analyzer::computeBandwidth(Tracker& tracker) {
      LayerVector::iterator layIt;
      ModuleVector::iterator modIt;
      ModuleVector* aLay;
      double hitChannels;

      // Clear and reset the histograms
      chanHitDistribution.Reset();
      bandwidthDistribution.Reset();
      bandwidthDistributionSparsified.Reset();
      chanHitDistribution.SetNameTitle("NHitChannels", "Number of hit channels;Hit Channels;Modules");
      bandwidthDistribution.SetNameTitle("BandWidthDist", "Module Needed Bandwidth;Bandwidth (bps);Modules");
      bandwidthDistributionSparsified.SetNameTitle("BandWidthDistSp", "Module Needed Bandwidth (sparsified);Bandwidth (bps);Modules");
      chanHitDistribution.SetBins(200, 0., 400);
      bandwidthDistribution.SetBins(100, 0., 6E+8);
      bandwidthDistributionSparsified.SetBins(100, 0., 6E+8);
      bandwidthDistribution.SetLineColor(kBlack);
      bandwidthDistributionSparsified.SetLineColor(kRed);

      int nChips;
      LayerVector layerSet = tracker.getLayers();
      double nMB = tracker.getNMB();
      for (layIt=layerSet.begin(); layIt!=layerSet.end(); layIt++) {
        aLay = (*layIt)->getModuleVector();
        for (modIt=aLay->begin(); modIt!=aLay->end(); modIt++) {
          if ((*modIt)->getReadoutType()==Module::Strip) {
            for (int nFace=1; nFace<=(*modIt)->getNFaces() ; nFace++) {
              hitChannels = (*modIt)->getHitOccupancyPerEvent()*nMB*((*modIt)->getNChannelsFace(nFace));
              chanHitDistribution.Fill(hitChannels);
              nChips=int(ceil((*modIt)->getNChannelsFace(nFace)/128.));

              // TODO: place the computing model choice here

              // ACHTUNG!!!! whenever you change the numbers here, you have to change
              // also the numbers in the summary

              // Binary unsparsified (bps)
              bandwidthDistribution.Fill((16*nChips+(*modIt)->getNChannelsFace(nFace))*100E3);

              int spHdr = tracker.getSparsifiedHeaderBits((*modIt)->getType());
              int spPay = tracker.getSparsifiedPayloadBits((*modIt)->getType());      

              //cout << "sparsified header: " << spHdr << " payload: " << spPay << endl;
              // Binary sparsified
              bandwidthDistributionSparsified.Fill(((spHdr*nChips)+(hitChannels*spPay))*100E3);
            }
          }
        }
      }

      savingGeometryV.push_back(chanHitDistribution);
      savingGeometryV.push_back(bandwidthDistribution);
      savingGeometryV.push_back(bandwidthDistributionSparsified);
    }


    std::vector<double> Analyzer::average(TGraph& myGraph, std::vector<double> cuts) {
      std::vector<double> averages;
      if (cuts.size()<2) return averages;
      std::sort(cuts.begin(), cuts.end());
      int iBorder;
      int nBoarders=cuts.size();
      double valuesCount, valuesSum;
      double vx, vy;
      for (iBorder=0; iBorder<nBoarders-1; ++iBorder) {
        // Here we have a cut between
        // cuts[iBorder] and cuts[iBorder+1]

        //std::cerr << cuts[iBorder] << "< x <= "
        //<< cuts[iBorder+1] << std::endl;

        // Average on points within the cut
        valuesSum=0;
        valuesCount=0;
        for (int iPoint=0; iPoint<myGraph.GetN(); ++iPoint) {
          myGraph.GetPoint(iPoint, vx, vy);
          if ((vx>=cuts[iBorder])
              && (vx<cuts[iBorder+1])) {
            valuesCount++;
            valuesSum+=vy;
          }
        }
        averages.push_back(valuesSum/valuesCount);
      }

      return averages;
    }


    int graphBag::buildAttribute(bool ideal, bool isTrigger) {
      int result;
      if (ideal) result = IdealGraph;
      else result = RealGraph;

      if (isTrigger) result |= TriggerGraph;
      else result |= StandardGraph;

      return result;
    }

    std::map<double, TGraph>& Analyzer::getRhoGraphs(bool ideal, bool isTrigger) {
      int attribute = graphBag::buildAttribute(ideal, isTrigger);
      attribute |= graphBag::RhoGraph;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<double, TGraph>& Analyzer::getPhiGraphs(bool ideal, bool isTrigger) {
      int attribute = graphBag::buildAttribute(ideal, isTrigger);
      attribute |= graphBag::PhiGraph;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<double, TGraph>& Analyzer::getDGraphs(bool ideal, bool isTrigger) {
      int attribute = graphBag::buildAttribute(ideal, isTrigger);
      attribute |= graphBag::DGraph;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<double, TGraph>& Analyzer::getCtgThetaGraphs(bool ideal, bool isTrigger) {
      int attribute = graphBag::buildAttribute(ideal, isTrigger);
      attribute |= graphBag::CtgthetaGraph;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<double, TGraph>& Analyzer::getZ0Graphs(bool ideal, bool isTrigger) {
      int attribute = graphBag::buildAttribute(ideal, isTrigger);
      attribute |= graphBag::Z0Graph;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<double, TGraph>& Analyzer::getPGraphs(bool ideal, bool isTrigger) {
      int attribute = graphBag::buildAttribute(ideal, isTrigger);
      attribute |= graphBag::PGraph;
      return myGraphBag.getGraphs(attribute);
    }

    TH1D& Analyzer::getHistoOptimalSpacing(bool actualWindow) {
      if (actualWindow) return optimalSpacingDistributionAW;
      else return optimalSpacingDistribution;
    }

  }


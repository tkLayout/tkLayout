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

#undef MATERIAL_SHADOW

namespace insur {

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

  const double profileBag::Triggerable    = 0.;
  const int profileBag::TriggeredProfile  = 0x0000007;
  const int profileBag::TriggeredFractionProfile  = 0x0000008;
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
    
#ifdef DEBUG_PERFORMANCE
    clock_t starttime = clock();
#endif
    int nTracks;
    double etaStep, eta, theta, phi;
    
    // prepare etaStep, phiStep, nTracks, nScans
    if (etaSteps > 1) etaStep = etaMax / (double)(etaSteps - 1);
    else etaStep = etaMax;
    nTracks = etaSteps;

    prepareTriggerPerformanceHistograms(nTracks, etaMax, triggerMomenta, thresholdProbabilities);
    
    // reset the list of tracks
    std::vector<Track> tv;
    std::vector<Track> tvIdeal;
    
    // used fixed phi
    //phi = PI / 2.0;
    
    // Loop over nTracks (eta range [0, etaMax])
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
      //		<< ", material.radiation = " << tmp.radiation
      //		<< ", material.interaction = " << tmp.interaction
      //		<< std::endl;
      
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
	  track.addIPConstraint(tracker.getRError(),tracker.getZError());
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

#ifdef DEBUG_PERFORMANCE
    std::cerr << "DEBUG_PERFORMANCE: material summary by analyzeTrigger(): ";
    clock_t endtime = clock();
    std::cerr << "elapsed time: " << diffclock(endtime, starttime)/1000. << "s" << std::endl;
#endif

    // Compute the number of triggering points along the selected tracks
    fillTriggerEfficiencyGraphs(triggerMomenta, tv);

    // Fill the trigger performance maps
    fillTriggerPerformanceMaps(tracker);

#ifdef DEBUG_PERFORMANCE
    starttime = clock();
#endif
    calculateGraphs(momenta, tv, graphBag::TriggerGraph | graphBag::RealGraph);
    calculateGraphs(momenta, tvIdeal, graphBag::TriggerGraph | graphBag::IdealGraph);
#ifdef DEBUG_PERFORMANCE
    std::cerr << "DEBUG_PERFORMANCE: tracking performance summary by analyzeTrigger(): ";
    endtime = clock();
    std::cerr << "elapsed time: " << diffclock(endtime, starttime)/1000. << "s" << std::endl;
#endif
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
					     const std::vector<Track>& trackVector) {
    
    // Prepare the graphs to record the number of triggered points
    //std::map<double, TGraph>& trigGraphs = myGraphBag.getGraphs(graphBag::TriggerGraph|graphBag::TriggeredGraph);
    std::map<double, TProfile>& trigProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggeredProfile);
    std::map<double, TProfile>& trigFractionProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggeredFractionProfile);

    TProfile& totalProfile = trigProfiles[profileBag::Triggerable];
    
    double eta;

    for (std::vector<Track>::const_iterator itTrack = trackVector.begin();
	 itTrack != trackVector.end(); ++itTrack) {
      const Track& myTrack=(*itTrack);

      eta = - log(tan(myTrack.getTheta() / 2));
      int nHits = myTrack.nActiveHits(false, false);
      totalProfile.Fill(eta, nHits);

      for(std::vector<double>::const_iterator itMomentum = triggerMomenta.begin();
	  itMomentum!=triggerMomenta.end(); ++itMomentum) {
	TProfile& myProfile = trigProfiles[(*itMomentum)];
	TProfile& myFractionProfile = trigFractionProfiles[(*itMomentum)];
	double nExpectedTriggerPoints = myTrack.expectedTriggerPoints(*itMomentum);
	if (nExpectedTriggerPoints>=0) { // sanity check (! nan)
	  myProfile.Fill(eta, nExpectedTriggerPoints);
	  if (nHits>0) {
	    myFractionProfile.Fill(eta, nExpectedTriggerPoints*100/double(nHits));
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
#ifdef DEBUG_PERFORMANCE
        clock_t starttime = clock();
#endif
        int nTracks;
        double etaStep, eta, theta, phi;
        clearMaterialBudgetHistograms();
        clearCells();
        // prepare etaStep, phiStep, nTracks, nScans
        if (etaSteps > 1) etaStep = etaMax / (double)(etaSteps - 1);
        else etaStep = etaMax;
        nTracks = etaSteps;
        // reset the number of bins and the histogram boundaries (0.0 to etaMax) for all histograms, recalculate the cell boundaries
        setHistogramBinsBoundaries(nTracks, 0.0, etaMax);
        setCellBoundaries(nTracks, 0.0, outer_radius + volume_width, 0.0, etaMax);
        // reset the list of tracks
	std::vector<Track> tv;
        std::vector<Track> tvIdeal;
        // tv.clear();
        // used fixed phi
        // phi = PI / 2.0;
        //      loop over nTracks (eta range [0, etaMax])
        for (int i_eta = 0; i_eta < nTracks; i_eta++) {
            phi = myDice.Rndm() * PI * 2.0;
            Material tmp;
            Track track;
            eta = i_eta * etaStep;
            theta = 2 * atan(pow(E, -1 * eta)); // TODO: switch to exp() here
            track.setTheta(theta);
            track.setPhi(phi);
            //      active volumes, barrel
            tmp = analyzeModules(mb.getBarrelModuleCaps(), eta, theta, phi, track);
            ractivebarrel.Fill(eta, tmp.radiation);
            iactivebarrel.Fill(eta, tmp.interaction);
            rbarrelall.Fill(eta, tmp.radiation);
            ibarrelall.Fill(eta, tmp.interaction);
            ractiveall.Fill(eta, tmp.radiation);
            iactiveall.Fill(eta, tmp.interaction);
            rglobal.Fill(eta, tmp.radiation);
            iglobal.Fill(eta, tmp.interaction);
            //      active volumes, endcap
            tmp = analyzeModules(mb.getEndcapModuleCaps(), eta, theta, phi, track);
            ractiveendcap.Fill(eta, tmp.radiation);
            iactiveendcap.Fill(eta, tmp.interaction);
            rendcapall.Fill(eta, tmp.radiation);
            iendcapall.Fill(eta, tmp.interaction);
            ractiveall.Fill(eta, tmp.radiation);
            iactiveall.Fill(eta, tmp.interaction);
            rglobal.Fill(eta, tmp.radiation);
            iglobal.Fill(eta, tmp.interaction);
            //      services, barrel
            tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getBarrelServices(), eta, theta, track);
            rserfbarrel.Fill(eta, tmp.radiation);
            iserfbarrel.Fill(eta, tmp.interaction);
            rbarrelall.Fill(eta, tmp.radiation);
            ibarrelall.Fill(eta, tmp.interaction);
            rserfall.Fill(eta, tmp.radiation);
            iserfall.Fill(eta, tmp.interaction);
            rglobal.Fill(eta, tmp.radiation);
            iglobal.Fill(eta, tmp.interaction);
            //      services, endcap
            tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getEndcapServices(), eta, theta, track);
            rserfendcap.Fill(eta, tmp.radiation);
            iserfendcap.Fill(eta, tmp.interaction);
            rendcapall.Fill(eta, tmp.radiation);
            iendcapall.Fill(eta, tmp.interaction);
            rserfall.Fill(eta, tmp.radiation);
            iserfall.Fill(eta, tmp.interaction);
            rglobal.Fill(eta, tmp.radiation);
            iglobal.Fill(eta, tmp.interaction);
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
            //      supports, tubes
            tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, MaterialProperties::o_sup);
            rlazytube.Fill(eta, tmp.radiation);
            ilazytube.Fill(eta, tmp.interaction);
            rlazyall.Fill(eta, tmp.radiation);
            ilazyall.Fill(eta, tmp.interaction);
            rglobal.Fill(eta, tmp.radiation);
            iglobal.Fill(eta, tmp.interaction);
            //      supports, barrel tubes
            tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, MaterialProperties::t_sup);
            rlazybtube.Fill(eta, tmp.radiation);
            ilazybtube.Fill(eta, tmp.interaction);
            rlazyall.Fill(eta, tmp.radiation);
            ilazyall.Fill(eta, tmp.interaction);
            rglobal.Fill(eta, tmp.radiation);
            iglobal.Fill(eta, tmp.interaction);
            //      supports, user defined
            tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, MaterialProperties::u_sup);
            rlazyuserdef.Fill(eta, tmp.radiation);
            ilazyuserdef.Fill(eta, tmp.interaction);
            rlazyall.Fill(eta, tmp.radiation);
            ilazyall.Fill(eta, tmp.interaction);
            rglobal.Fill(eta, tmp.radiation);
            iglobal.Fill(eta, tmp.interaction);
            //      pixels, if they exist
            if (pm != NULL) {
                analyzeModules(pm->getBarrelModuleCaps(), eta, theta, phi, track, true);
                analyzeModules(pm->getEndcapModuleCaps(), eta, theta, phi, track, true);
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
		//		if ((eta >1.617) &&(theta>0.383)) {
		//		  Track::debugRZCorrelationMatrix = true;
		//		  Track::debugRZCovarianceMatrix = true;
		//		  Track::debugRZErrorPropagation = true;
		//		}
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
		//					     0,
		//					     sqrt( averageSquaredHits - averageHits*averageHits) );
		
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
		  //		      << " out of " << probabilities.size()
		  //		      << " == " << nActive
		  //		      << endl;
		  // std::cerr << "      PROBABILITY = " << probability << endl << endl;
		  //}
		  hadronGoodTracksFraction.at(i).SetPoint(hadronGoodTracksFraction.at(i).GetN(),
							  eta,
							  probability);
		}
	      }
            }
        }

#ifdef DEBUG_PERFORMANCE
        std::cerr << "DEBUG_PERFORMANCE: material summary by analyzeMaterialBudget(): ";
        clock_t endtime = clock();
        std::cerr << "elapsed time: " << diffclock(endtime, starttime)/1000. << "s" << std::endl;
#endif
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
#ifdef DEBUG_PERFORMANCE
	starttime = clock();
#endif
	calculateGraphs(momenta, tv, graphBag::StandardGraph | graphBag::RealGraph);
	calculateGraphs(momenta, tvIdeal, graphBag::StandardGraph | graphBag::IdealGraph);
#ifdef DEBUG_PERFORMANCE
        std::cerr << "DEBUG_PERFORMANCE: tracking performance summary by analyzeMaterialBudget(): ";
        endtime = clock();
        std::cerr << "elapsed time: " << diffclock(endtime, starttime)/1000. << "s" << std::endl;
#endif
}
    }

	void Analyzer::analyzePower(Tracker& tracker) {
		computeIrradiatedPowerConsumption(tracker);
		preparePowerHistograms();
		fillPowerMap(tracker);
	}

	void Analyzer::computeTriggerFrequency(Tracker& tracker) {
		std::map<std::string, std::map<std::pair<int,int>, int> >   triggerFrequencyCounts;
		std::map<std::string, std::map<std::pair<int,int>, double> > triggerFrequencyAverageTrue, triggerFrequencyAverageFake; // trigger frequency by module in Z and R, averaged over Phi
		LayerVector& layers = tracker.getLayers();

        triggerFrequencyTrueSummaries_.clear();
        triggerFrequencyFakeSummaries_.clear();
		triggerRateSummaries_.clear();
		triggerPuritySummaries_.clear();
		triggerDataBandwidthSummaries_.clear();

		for(LayerVector::iterator layIt = layers.begin(); layIt != layers.end(); ++layIt) {
			ModuleVector* modules = (*layIt)->getModuleVector();
			if (!modules) {
				std::cerr << "ERROR in Analyzer::computeTriggerFrequency: cannot retrieve moduleVector from the layer\n";
				return;
			}
			std::string cntName = (*layIt)->getContainerName();
			if (cntName == "") {
				cout << "computeTriggerFrequency(): Skipping layer with no container name (" << modules->size() << " modules)." << endl;
				continue;
			}
			for (ModuleVector::iterator modIt = modules->begin(); modIt != modules->end(); ++modIt) {
				Module* module = (*modIt);
				XYZVector center = module->getMeanPoint();
				if ((center.Z()<0) || (center.Phi()<0) || (center.Phi()>M_PI/2) || (module->getStereoDistance()==0.0)) continue;
				int curCnt = triggerFrequencyCounts[cntName][make_pair(module->getLayer(), module->getRing())]++;
				double curAvgTrue = triggerFrequencyAverageTrue[cntName][make_pair(module->getLayer(),module->getRing())];
				double curAvgFake = triggerFrequencyAverageFake[cntName][make_pair(module->getLayer(),module->getRing())];

				curAvgTrue  = curAvgTrue + (module->getTriggerFrequencyTruePerEvent()*tracker.getNMB() - curAvgTrue)/(curCnt+1);
				curAvgFake  = curAvgFake + (module->getTriggerFrequencyFakePerEvent()*pow(tracker.getNMB(),2) - curAvgFake)/(curCnt+1); // triggerFrequencyFake scales with the square of Nmb!

				double curAvgTotal = curAvgTrue + curAvgFake;

				triggerFrequencyAverageTrue[cntName][make_pair(module->getLayer(), module->getRing())] = curAvgTrue;			
				triggerFrequencyAverageFake[cntName][make_pair(module->getLayer(), module->getRing())] = curAvgFake;	

				int triggerDataHeaderBits  = tracker.getModuleType(module->getType()).getTriggerDataHeaderBits();
				int triggerDataPayloadBits = tracker.getModuleType(module->getType()).getTriggerDataPayloadBits();


				std::stringstream ss1(""), ss2(""), ss3(""), ss4(""), ss5("");
				ss1.precision(6);
				ss2.precision(6);
				ss3.precision(6);
				ss4.precision(6);
				ss5.precision(6);
				ss1.setf(ios::fixed, ios::floatfield);
				ss2.setf(ios::fixed, ios::floatfield);
				ss3.setf(ios::fixed, ios::floatfield);
				ss4.setf(ios::fixed, ios::floatfield);
				ss5.setf(ios::fixed, ios::floatfield);
				ss1 << curAvgTrue;
				ss2 << curAvgFake;
				ss3 << curAvgTotal;
				ss4 << curAvgTrue/(curAvgTrue+curAvgFake);
				ss5 << triggerDataHeaderBits + curAvgTotal*triggerDataPayloadBits;
				triggerFrequencyTrueSummaries_[cntName].setCell(module->getLayer(), module->getRing(), ss1.str());
				triggerFrequencyFakeSummaries_[cntName].setCell(module->getLayer(), module->getRing(), ss2.str());
				triggerRateSummaries_[cntName].setCell(module->getLayer(), module->getRing(), ss3.str());				
				triggerPuritySummaries_[cntName].setCell(module->getLayer(), module->getRing(), ss4.str());				
				triggerDataBandwidthSummaries_[cntName].setCell(module->getLayer(), module->getRing(), ss5.str());

				if (!triggerFrequencyTrueSummaries_[cntName].hasCell(module->getLayer(), 0)) {
					std::stringstream ss("");
					ss << module->getLayer();
					triggerFrequencyTrueSummaries_[cntName].setCell(module->getLayer(), 0, ss.str());
					triggerFrequencyFakeSummaries_[cntName].setCell(module->getLayer(), 0, ss.str());
				}
				if (!triggerFrequencyTrueSummaries_[cntName].hasCell(0, module->getRing())) {
					std::stringstream ss("");
					ss << module->getRing();
					triggerFrequencyTrueSummaries_[cntName].setCell(0, module->getRing(), ss.str());
					triggerFrequencyFakeSummaries_[cntName].setCell(0, module->getRing(), ss.str());
				}
				if (!triggerRateSummaries_[cntName].hasCell(module->getLayer(), 0)) {
					std::stringstream ss("");
					ss << module->getLayer();
					triggerRateSummaries_[cntName].setCell(module->getLayer(), 0, ss.str());
					triggerPuritySummaries_[cntName].setCell(module->getLayer(), 0, ss.str());
					triggerDataBandwidthSummaries_[cntName].setCell(module->getLayer(), 0, ss.str());
				}
				if (!triggerRateSummaries_[cntName].hasCell(0, module->getRing())) {
					std::stringstream ss("");
					ss << module->getRing();
					triggerRateSummaries_[cntName].setCell(0, module->getRing(), ss.str());
					triggerPuritySummaries_[cntName].setCell(0, module->getRing(), ss.str());
					triggerDataBandwidthSummaries_[cntName].setCell(0, module->getRing(), ss.str());
				}
			}
			triggerFrequencyTrueSummaries_[cntName].setCell(0, 0, "Ring/<br>Layer");
			triggerFrequencyFakeSummaries_[cntName].setCell(0, 0, "Ring/<br>Layer");
			triggerRateSummaries_[cntName].setCell(0, 0, "Ring/<br>Layer");
			triggerPuritySummaries_[cntName].setCell(0, 0, "Ring/<br>Layer");
			triggerDataBandwidthSummaries_[cntName].setCell(0, 0, "Ring/<br>Layer");
		}
	}
	
void Analyzer::computeIrradiatedPowerConsumption(Tracker& tracker) {
	// loop over layers
	// 	  loop over modules
	// 	     compute the power formula for the given module P(V*, I(Fluence(Z, r), T*, alpha*))
	//		 populate R-Z table
	
	double numInvFemtobarns = tracker.getNumInvFemtobarns();
	double operatingTemp    = tracker.getOperatingTemp();
	double chargeDepletionVoltage    = tracker.getChargeDepletionVoltage();
	double alphaParam       = tracker.getAlphaParam();
    double referenceTemp    = tracker.getReferenceTemp();

	cout << "numInvFemtobarns = " << tracker.getNumInvFemtobarns() << endl;
	cout << "operatingTemp    = " << tracker.getOperatingTemp() << endl;
	cout << "chargeDepletionVoltage    = " << tracker.getChargeDepletionVoltage() << endl;
	cout << "alphaParam       = " << tracker.getAlphaParam() << endl;
    cout << "referenceTemp    = " << tracker.getReferenceTemp() << endl;
	irradiatedPowerConsumptionSummaries_.clear();
	

	LayerVector& layers = tracker.getLayers();
	for(LayerVector::iterator layIt = layers.begin(); layIt != layers.end(); ++layIt) {
		ModuleVector* modules = (*layIt)->getModuleVector();
		if (!modules) {
			std::cerr << "ERROR in Analyzer::computeTriggerFrequency: cannot retrieve moduleVector from the layer\n";
			return;
		}
		std::string cntName = (*layIt)->getContainerName();
		if (cntName == "") {
			cout << "computeIrradiatedPowerConsumption(): Skipping layer with no container name (" << modules->size() << " modules)." << endl;
			continue;
		}
		for (ModuleVector::iterator modIt = modules->begin(); modIt != modules->end(); ++modIt) {
			Module* module = (*modIt); 
			XYZVector center = module->getMeanPoint();
			if ((center.Z()<0) || (center.Phi()<0) || (center.Phi()>M_PI/2)) continue;
			double volume  = tracker.getSensorThickness(module->getType()) * module->getArea() / 1000.0 * module->getNFaces(); // volume is in cm^3
			double x  = center.Z()/25;
			double y  = center.Rho()/25;
			double x1 = floor(x);
			double x2 = ceil(x);
			double y1 = floor(y);
			double y2 = ceil(y);
			double irr11 = tracker.getIrradiationMap()[make_pair(int(x1), int(y1))]; 
			double irr21 = tracker.getIrradiationMap()[make_pair(int(x2), int(y1))];
			double irr12 = tracker.getIrradiationMap()[make_pair(int(x1), int(y2))];
			double irr22 = tracker.getIrradiationMap()[make_pair(int(x2), int(y2))];
			double irrxy = irr11/((x2-x1)*(y2-y1))*(x2-x)*(y2-y) + irr21/((x2-x1)*(y2-y1))*(x-x1)*(y2-y) + irr12/((x2-x1)*(y2-y1))*(x2-x)*(y-y1) + irr22/((x2-x1)*(y2-y1))*(x-x1)*(y-y1); // bilinear interpolation
			double fluence = irrxy * numInvFemtobarns * 1e15 * 77 * 1e-3; // fluence is in 1MeV-equiv-neutrons/cm^2 
			double leakCurrentScaled = alphaParam * fluence * volume * pow((operatingTemp+273.15) / (referenceTemp+273.15), 2) * exp(-1.21/(2*8.617334e-5)*(1/(operatingTemp+273.15)-1/(referenceTemp+273.15))); 
			double irradiatedPowerConsumption = leakCurrentScaled * chargeDepletionVoltage;			
			//cout << "mod irr: " << cntName << "," << module->getLayer() << "," << module->getRing() << ";  " << module->getThickness() << "," << center.Rho() << ";  " << volume << "," << fluence << "," << leakCurrentScaled << "," << irradiatedPowerConsumption << endl;
			module->setIrradiatedPowerConsumption(irradiatedPowerConsumption);
			std::stringstream ss("");
			ss.precision(6);
			ss.setf(ios::fixed, ios::floatfield);
			ss << irradiatedPowerConsumption;
			irradiatedPowerConsumptionSummaries_[cntName].setCell(module->getLayer(), module->getRing(), ss.str());
			if (!irradiatedPowerConsumptionSummaries_[cntName].hasCell(module->getLayer(), 0)) {
				std::stringstream ss("");
				ss << module->getLayer();
				irradiatedPowerConsumptionSummaries_[cntName].setCell(module->getLayer(), 0, ss.str());
				irradiatedPowerConsumptionSummaries_[cntName].setCell(module->getLayer(), 0, ss.str());
			}
			if (!irradiatedPowerConsumptionSummaries_[cntName].hasCell(0, module->getRing())) {
				std::stringstream ss("");
				ss << module->getRing();
				irradiatedPowerConsumptionSummaries_[cntName].setCell(0, module->getRing(), ss.str());
				irradiatedPowerConsumptionSummaries_[cntName].setCell(0, module->getRing(), ss.str());
			}
		}
		irradiatedPowerConsumptionSummaries_[cntName].setCell(0, 0, "Ring/<br>Layer");
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
    
    unsigned int nLocalMasses;
    unsigned int nExitingMasses;

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
	    nLocalMasses = myModuleCap->localMassCount();
	    nExitingMasses = myModuleCap->exitingMassCount();
	  } else { // sort by Component tag
	    nLocalMasses = myModuleCap->localMassCompCount();
	    nExitingMasses = myModuleCap->exitingMassCompCount();
	  }
	  for (unsigned int iLocalMasses=0; iLocalMasses < nLocalMasses; ++iLocalMasses) {
	    if (byMaterial) materialTag = myModuleCap->getLocalTag(iLocalMasses); // sort by Material tag
	    else materialTag = myModuleCap->getLocalTagComp(iLocalMasses);           // sort by Component tag
	    materialTagIt = find(materialTagV.begin(), materialTagV.end(), materialTag);
	    if (materialTagIt==materialTagV.end()) {
	      materialTagV.push_back(materialTag);
	    }
	  }
	  for (unsigned int iExitingMasses=0; iExitingMasses < nExitingMasses; ++iExitingMasses) {
	    if (byMaterial) materialTag = myModuleCap->getExitingTag(iExitingMasses); // sort by Material tag
	    else materialTag = myModuleCap->getExitingTagComp(iExitingMasses);           // sort by Component tag
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
	      //	  cout << "Here's a module: id = " << myModule->getId()
	      //	       << ", tag = " << myModule->getTag()
	      //	       << ", type = " << myModule->getType()
	      //	       << ", ring = " << myModule->getRing()
	      //	       << endl;
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
            double eta, double theta, double phi, Track& t, bool isPixel) {
        std::vector<std::vector<ModuleCap> >::iterator iter = tr.begin();
        std::vector<std::vector<ModuleCap> >::iterator guard = tr.end();
        Material res, tmp;
        res.radiation= 0.0;
        res.interaction = 0.0;
        while (iter != guard) {
            tmp = findModuleLayerRI(*iter, eta, theta, phi, t, isPixel);
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
            double eta, double theta, double phi, Track& t, bool isPixel) {
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

 

  /**
   * This convenience function resets and empties all histograms for
   * the trigger performance, so they are ready for a new round of
   * analysis.
   * @parameter triggerMomenta the vector of pt to test the trigger
   */
  
  void Analyzer::prepareTriggerPerformanceHistograms(const int& nTracks, const double& etaMax, const std::vector<double>& triggerMomenta, const std::vector<double>& thresholdProbabilities) {
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
      //TGraph graph;
      TProfile profile("dummyName", "dummyTitle", nbins, 0, 2.4);
      //elemGraph.first = *iter;
      //elemGraph.second = graph;
      elemProfile.first = *iter;
      elemProfile.second = profile;
      elemFractionProfile.first = *iter;
      elemFractionProfile.second = profile;
      // Prepare plots: triggered graphs
      // trigGraphs.insert(elemGraph);
      trigProfiles.insert(elemProfile);
      trigFractionProfiles.insert(elemFractionProfile);
      // trigGraphs[*iter].SetTitle("Average triggered points;#eta;Triggered points <N>");
      trigProfiles[*iter].SetTitle("Average triggered points;#eta;Triggered points <N>");
      trigFractionProfiles[*iter].SetTitle("Average trigger efficiency;#eta;Efficiency [%]");
      aName.str(""); aName << "triggered_vs_eta" << *iter << "_profile";      
      trigProfiles[*iter].SetName(aName.str().c_str());
      aName.str(""); aName << "triggered_vs_eta" << *iter << "_fractionProfile";
      trigFractionProfiles[*iter].SetName(aName.str().c_str());
    }

    //std::pair<double, TGraph> elemTotalGraph;
    std::pair<double, TProfile> elemTotalProfile;
    //TGraph totalGraph;
    TProfile totalProfile("dummyName", "dummyTitle", nbins, 0, 2.4); // where is this 2.4 defined normally? TODO: fix it
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
        std::map <std::string, TH2D> etaProfileByType;
        TH2D* aPlot;
        std::string aType;
        
        
        // Optimize the track creation on the real tracker
        std::pair <double, double> etaMinMax = tracker.getEtaMinMax();
        double absMinEta = fabs(etaMinMax.first);
        double absMaxEta = fabs(etaMinMax.second);
        double maxEta = (absMinEta>absMaxEta) ? absMinEta : absMaxEta;
        
        // Computing the margin of the tracks to shoot
        double randomPercentMargin = 0.1;
        double randomSpan = (etaMinMax.second - etaMinMax.first)*(1. + randomPercentMargin);
        double randomBase = etaMinMax.first - (etaMinMax.second - etaMinMax.first)*(randomPercentMargin)/2.;
        
        // Initialize random number generator, counters and histograms
        myDice.SetSeed(MY_RANDOM_SEED);
        createResetCounters(tracker, moduleTypeCount);
        
        /*for (std::map <std::string, TH2D*>::iterator it = etaProfileByType.begin();
         it!=etaProfileByType.end(); it++) {
      aPlot = (*it).second;
      if (aPlot) delete aPlot;
      }*/
        etaProfileByType.clear();
        
        for (std::map <std::string, int>::iterator it = moduleTypeCount.begin();
        it!=moduleTypeCount.end(); it++) {
            // std::cerr << "Creating plot for module type " << (*it).first << std::endl; //debug
            aPlot = new TH2D( (*it).first.c_str(), (*it).first.c_str(),
                    100, 0., maxEta*1.1,
                    1000, 0., 10.);
            etaProfileByType[(*it).first]=(*aPlot);
            delete aPlot;
        }
        
        LayerVector::iterator layIt;
        ModuleVector* moduleV = NULL;
        ModuleVector::iterator modIt;
        ModuleVector allModules;
        LayerVector& layerSet = tracker.getLayers();
        double zError = tracker.getZError();
        
        // Build the proper list of modules
        for (layIt=layerSet.begin(); layIt!=layerSet.end(); layIt++) {
            moduleV = (*layIt)->getModuleVector();
            if (!moduleV) {
              std::cerr << "ERROR in Analyzer::analyzeGeometry: cannot retrieve moduleVector from the layer" << std::endl;
              return;
            }
            for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
                // I pre-compute the boxes to reduce the calculations
                (*modIt)->computeBoundaries(zError);
                allModules.push_back(*modIt);
            }
        }
        
        // The real simulation
        std::pair <XYZVector, double> aLine;
        ModuleVector hitModules;
        
#ifdef DEBUG_PERFORMANCE
        clock_t starttime = clock();
#endif
        
        std::cout << "Shooting tracks..." << std::endl;
        int nTrackHits;
        int nTracksPerSide = int(pow(nTracks, 0.5));
        int nBlocks = int(nTracksPerSide/2.);
        nTracks = nTracksPerSide*nTracksPerSide;
        mapPhiEta.SetBins(nBlocks, -1*M_PI, M_PI, nBlocks, -3., 3.);
        TH2I mapPhiEtaCount("mapPhiEtaCount ", "phi Eta hit count", nBlocks, -1*M_PI, M_PI, nBlocks, -3., 3.);
        TH2D total2D("total2d", "Total 2D", 100, 0., maxEta*1.2, 4000 , 0., 40.);
        
        // Shoot nTracksPerSide^2 tracks
        for (int i=0; i<nTracksPerSide; i++) {
            for (int j=0; j<nTracksPerSide; j++) {
                // Reset the hit counter
                nTrackHits=0;
                // Generate a straight track and collect the list of hit modules
                aLine = shootDirection(randomBase, randomSpan);
                hitModules = trackHit( XYZVector(0, 0, myDice.Gaus(0, zError)), aLine.first, &allModules);
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
                total2D.Fill(fabs(aLine.second), hitModules.size());                // Total number of hits
                mapPhiEta.Fill(aLine.first.Phi(), aLine.second, hitModules.size()); // phi, eta 2d plot
                mapPhiEtaCount.Fill(aLine.first.Phi(), aLine.second);               // Number of shot tracks
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
#ifdef DEBUG_PERFORMANCE
        std::cerr << "DEBUG_PERFORMANCE: tracks for analyzeGeometry(): ";
        clock_t endtime = clock();
        std::cerr << "elapsed time: " << diffclock(endtime, starttime)/1000. << "s" << std::endl;
#endif

        // Eta profile compute
        //TProfile *myProfile;
        
        etaProfileCanvas.cd();
        savingGeometryV.push_back(etaProfileCanvas);
        int plotCount=0;

        //TProfile* total = total2D.ProfileX("etaProfileTotal");
        totalEtaProfile = TProfile(*total2D.ProfileX("etaProfileTotal"));
        savingGeometryV.push_back(totalEtaProfile);
        totalEtaProfile.SetMarkerStyle(8);
        totalEtaProfile.SetMarkerColor(1);
        totalEtaProfile.SetMarkerSize(1.5);
        totalEtaProfile.SetTitle("Number of hit modules;#eta;Number of hits");
        if (totalEtaProfile.GetMaximum()<maximum_n_planes) totalEtaProfile.SetMaximum(maximum_n_planes);
        totalEtaProfile.Draw();
        std::string profileName;
        TProfile* myProfile;
        for (std::map <std::string, TH2D>::iterator it = etaProfileByType.begin();
        it!=etaProfileByType.end(); it++) {
            plotCount++;
            myProfile=(*it).second.ProfileX();
            savingGeometryV.push_back(*myProfile);
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
		hitChannels = (*modIt)->getOccupancyPerEvent()*nMB*((*modIt)->getNChannelsFace(nFace));
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


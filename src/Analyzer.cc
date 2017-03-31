/**
 * @file Analyzer.cc
 * @brief This is the implementation of the class that analyses a material budget
 */
#include <TH1D.h>
#include <TH2D.h>
#include <Analyzer.hh>
#include "MainConfigHandler.hh"
#include <HitNew.hh>
#include <TProfile.h>
#include <TLegend.h>
#include <Palette.hh>
#include "SimParms.hh"
#include "AnalyzerVisitors/MaterialBillAnalyzer.hh"
#include <Units.hh>
#include "TProfile.h"
#include "TMath.h"

#undef MATERIAL_SHADOW

#include <TError.h>
Int_t gErrorIgnoreLevel = kError;

namespace insur {

  int Analyzer::bsCounter =0;

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
    //colorPicker("pt2S");
    //colorPicker("rphi");
    //colorPicker("stereo");
    //colorPicker("ptIn");
    geomLite           = nullptr; geomLiteCreated=false;
    geomLiteXY         = nullptr; geomLiteXYCreated=false;
    geomLiteYZ         = nullptr; geomLiteYZCreated=false;
    geomLiteEC         = nullptr; geomLiteECCreated=false;
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
  Material Analyzer::findAllHits(MaterialBudget& mb, MaterialBudget* pm, TrackNew& track) {
    Material totalMaterial;
    //      active volumes, barrel
    totalMaterial  = findHitsModules(mb.getBarrelModuleCaps(), track);
    //      active volumes, endcap
    totalMaterial += findHitsModules(mb.getEndcapModuleCaps(), track);
    //      services, barrel
    totalMaterial += findHitsInactiveSurfaces(mb.getInactiveSurfaces().getBarrelServices(), track);
    //      services, endcap
    totalMaterial += findHitsInactiveSurfaces(mb.getInactiveSurfaces().getEndcapServices(), track);
    //      supports
    totalMaterial += findHitsInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), track);
    //      pixels, if they exist
    if (pm != NULL) {
      totalMaterial += findHitsModules(pm->getBarrelModuleCaps(), track, true);
      totalMaterial += findHitsModules(pm->getEndcapModuleCaps(), track, true);
      totalMaterial += findHitsInactiveSurfaces(pm->getInactiveSurfaces().getBarrelServices(), track, true);
      totalMaterial += findHitsInactiveSurfaces(pm->getInactiveSurfaces().getEndcapServices(), track, true);
      totalMaterial += findHitsInactiveSurfaces(pm->getInactiveSurfaces().getSupports(), track, true);
    }
    return totalMaterial;
  }


  /* TODO: finish this :-)
void Analyzer::createTaggedTrackCollection(std::vector<MaterialBudget*> materialBudgets,
                                           int etaSteps,
                                           double maxEta,
                                           bool forceClean = false) {
  
  if (forceClean) taggedTrackCollectionMap_.clear();
  if (taggedTrackCollectionMap_.size()!=0) return;
    
  double etaStep, eta, theta, phi;
  
  // prepare etaStep, phiStep, nTracks, nScans
  if (etaSteps > 1) etaStep = maxEta / (double)(etaSteps - 1);
  else etaStep = maxEta;
  const int nTracks& = etaSteps;

  
        
}
    
  */                                

  void Analyzer::analyzeTaggedTracking(MaterialBudget& mb,
				       const std::vector<double>& momenta,
				       const std::vector<double>& triggerMomenta,
				       const std::vector<double>& thresholdProbabilities,
				       bool isPixel,
				       bool& debugResolution,
				       int etaSteps,
				       MaterialBudget* pm) {

  auto& simParms = SimParms::getInstance();

  materialTracksUsed = etaSteps;

  int nTracks;
  double etaStep, eta, theta, phi, zPos;

  // prepare etaStep, phiStep, nTracks, nScans
  if (etaSteps > 1) etaStep = getEtaMaxTrigger() / (double)(etaSteps - 1);
  else etaStep = getEtaMaxTrigger();
  nTracks = etaSteps;

  // prepareTriggerPerformanceHistograms(nTracks, getEtaMaxTrigger(), triggerMomenta, thresholdProbabilities);

  // reset the list of tracks
  //std::map<string, std::vector<Track>> tv;
  //std::map<string, std::vector<Track>> tvIdeal;
  std::map<std::string, TrackNewCollectionMap> taggedTrackPtCollectionMap;
  std::map<std::string, TrackNewCollectionMap> taggedTrackPCollectionMap;
  std::map<std::string, TrackNewCollectionMap> taggedTrackPtCollectionMapIdeal;
  std::map<std::string, TrackNewCollectionMap> taggedTrackPCollectionMapIdeal;


  for (int i_eta = 0; i_eta < nTracks; i_eta++) {
    phi = myDice.Rndm() * M_PI * 2.0;
    Material tmp;
    TrackNew track;
    eta = i_eta * etaStep;
    theta = 2 * atan(exp(-eta));
    //std::cout << " track's phi = " << phi << std::endl; 
    track.setThetaPhiPt(theta,phi,1*Units::TeV);
    track.setOrigin(0., 0., 0.); // TODO: Not assuming z-error when analyzing resolution (missing implementation of non-zero track starting point in inactive hits)
    //track.setTheta(theta);
    //track.setPhi(phi);

    // Assign material to the track
    tmp = findAllHits(mb, pm, track);

    // Debug: material amount
    // std::cerr << "eta = " << eta
    //           << ", material.radiation = " << tmp.radiation
    //           << ", material.interaction = " << tmp.interaction
    //           << std::endl;

    // TODO: add the beam pipe as a user material eveywhere!
    // in a coherent way
    // Add the hit on the beam pipe
    double rPos  = 23.*Units::mm;
    double zPos  = rPos/tan(theta);

    HitNewPtr hit(new HitNew(rPos, zPos));
    hit->setAsPassive();

    Material material;
    material.radiation   = 0.0022761 / sin(theta);  // was 0.0023, adapted to fit CMSSW 81X 2016/11/30
    material.interaction = 0.0020334 / sin(theta);  // was 0.0019, adapted to fit CMSSW 81X 2016/11/30

    hit->setCorrectedMaterial(material);
    hit->setBeamPipe(true);
    track.addHit(std::move(hit));

    if (!track.hasNoHits()) {
      for (string tag : track.getTags()) {

        bool useIPConstraint = false;
        if (SimParms::getInstance().useIPConstraint()) {

          useIPConstraint = true;
          track.addIPConstraint(SimParms::getInstance().rphiErrorCollider(), SimParms::getInstance().zErrorCollider());
        }
        track.keepTaggedHitsOnly(tag,useIPConstraint);
        //track.sort();
        //track.setTriggerResolution(true); // TODO: remove this (?)

        track.addEfficiency();
        // For each momentum/transverse momentum compute the tracks error
        for (const auto& pIter : momenta ) {
          int    parameter = pIter/Units::MeV; // Store p or pT in MeV as int (key to the map)
          double momentum  = pIter;

          // Case I) Initial momentum is equal to pT
          double pT = momentum;
  
          // Active+passive material
          TrackNewPtr trackPt(new TrackNew(track));
          trackPt->resetPt(pT);
          //trackPt.pruneHits();                // Remove hits from a track that is not able to reach a given radius due to its limited momentum
          if (trackPt->getNActiveHits(tag,useIPConstraint)>=3) { // Only keep tracks which have minimum 3 active hits
            //trackPt.computeErrors();
            TrackNewCollectionMap &myMap     = taggedTrackPtCollectionMap[tag];
            TrackNewCollection &myCollection = myMap[parameter];
            myCollection.push_back(std::move(trackPt));
          }

          // Ideal (no material)
          TrackNewPtr idealTrackPt(new TrackNew(track));
          idealTrackPt->resetPt(pT);
          idealTrackPt->removeMaterial();
          if (idealTrackPt->getNActiveHits(tag,useIPConstraint)>=3) { // Only keep tracks which have minimum 3 active hits
            //idealTrackPt.computeErrors();
            TrackNewCollectionMap &myMapIdeal     = taggedTrackPtCollectionMapIdeal[tag];
            TrackNewCollection &myCollectionIdeal = myMapIdeal[parameter];
            myCollectionIdeal.push_back(std::move(idealTrackPt));
          }

          // Case II) Initial momentum is equal to p
          pT = momentum*sin(theta);

          // Active+passive material
          TrackNewPtr trackP(new TrackNew(track));
          trackP->resetPt(pT);
          //trackP.pruneHits();                // Remove hits from a track that is not able to reach a given radius due to its limited momentum
          if (trackP->getNActiveHits(tag,useIPConstraint)>=3) { // Only keep tracks which have minimum 3 active hits
            //trackP.computeErrors();
            TrackNewCollectionMap &myMapII     = taggedTrackPCollectionMap[tag];
            TrackNewCollection &myCollectionII = myMapII[parameter];
            myCollectionII.push_back(std::move(trackP));
          }

          // Ideal (no material)
          TrackNewPtr idealTrackP(new TrackNew(track));
          idealTrackP->resetPt(pT);
          idealTrackP->removeMaterial();
          if (idealTrackP->getNActiveHits(tag,useIPConstraint)>=3) { // Only keep tracks which have minimum 3 active hits
            //idealTrackP.computeErrors();
            TrackNewCollectionMap &myMapIdealII     = taggedTrackPCollectionMapIdeal[tag];
            TrackNewCollection &myCollectionIdealII = myMapIdealII[parameter];
            myCollectionIdealII.push_back(std::move(idealTrackP));
          }
        }
      }
    }
  }

  if (!isPixel) {
    // Momentum = Pt
    for (/*const*/ auto& ttcmIt : taggedTrackPtCollectionMap) {
      const string& myTag = ttcmIt.first;
      clearGraphsPt(GraphBag::RealGraph, myTag);
      /*const*/ TrackNewCollectionMap& myTrackCollection = ttcmIt.second;
      for (const auto& tcmIt : myTrackCollection) {
      	const int &parameter = tcmIt.first;
	      const TrackNewCollection& myCollection = tcmIt.second;
	      //std::cout << myCollection.size() << std::endl;
	      calculateGraphsConstPt(parameter, myCollection, GraphBag::RealGraph, myTag);
      }
    }
    for (/*const*/ auto& ttcmIt : taggedTrackPtCollectionMapIdeal) {
      const string& myTag = ttcmIt.first;
      clearGraphsPt(GraphBag::IdealGraph, myTag);
      /*const*/ TrackNewCollectionMap& myTrackCollection = ttcmIt.second;
      for (const auto& tcmIt : myTrackCollection) {
	      const int &parameter = tcmIt.first;
	      const TrackNewCollection& myCollection = tcmIt.second;
	      calculateGraphsConstPt(parameter, myCollection, GraphBag::IdealGraph, myTag);
      }
    }

    // Momentum = P
    for (/*const*/ auto& ttcmIt : taggedTrackPCollectionMap) {
      const string& myTag = ttcmIt.first;
      clearGraphsP(GraphBag::RealGraph, myTag);
      /*const*/ TrackNewCollectionMap& myTrackCollection = ttcmIt.second;
      for (const auto& tcmIt : myTrackCollection) {
	      const int &parameter = tcmIt.first;
	      const TrackNewCollection& myCollection = tcmIt.second;
	      //std::cout << myCollection.size() << std::endl;
	      calculateGraphsConstP(parameter, myCollection, GraphBag::RealGraph, myTag);
      }
    }
    for (/*const*/ auto& ttcmIt : taggedTrackPCollectionMapIdeal) {
      const string& myTag = ttcmIt.first;
      clearGraphsP(GraphBag::IdealGraph, myTag);
      /*const*/ TrackNewCollectionMap& myTrackCollection = ttcmIt.second;
      for (const auto& tcmIt : myTrackCollection) {
	      const int &parameter = tcmIt.first;
	      const TrackNewCollection& myCollection = tcmIt.second;
	      calculateGraphsConstP(parameter, myCollection, GraphBag::IdealGraph, myTag);
      }
    }
  }
  if (debugResolution) calculateParametrizedResolutionPlots(taggedTrackPtCollectionMap);
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

    auto& simParms = SimParms::getInstance();

    materialTracksUsed = etaSteps;

    int nTracks;
    double etaStep, z0, eta, theta, phi;
    double zError = simParms.zErrorCollider();

    // prepare etaStep, phiStep, nTracks, nScans
    if (etaSteps > 1) etaStep = getEtaMaxTrigger() / (double)(etaSteps - 1);
    else etaStep = getEtaMaxTrigger();
    nTracks = etaSteps;

    prepareTriggerPerformanceHistograms(nTracks, getEtaMaxTrigger(), triggerMomenta, thresholdProbabilities);

    // reset the list of tracks
    std::vector<Track> tv;

    // Loop over nTracks (eta range [0, getEtaMaxTrigger()])
    for (int i_eta = 0; i_eta < nTracks; i_eta++) {
      phi = myDice.Rndm() * M_PI * 2.0;
      z0 = myDice.Gaus(0, zError);
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

        track.addEfficiency();
        if (track.nActiveHits(true)>0) { // At least 3 points are needed to measure the arrow
          tv.push_back(track);
        }    
      }
    }

    // Compute the number of triggering points along the selected tracks
    fillTriggerEfficiencyGraphs(tracker, triggerMomenta, tv);

    // Fill the trigger performance maps
    fillTriggerPerformanceMaps(tracker);

  }

//
// Check that a file can be opened
//
bool Analyzer::checkFile(const std::string& fileName, const std::string& filePath)
{
  fstream     file;
  std::string fullFileName(filePath+"/"+fileName);
  file.open(fullFileName);
  if (file.is_open()) {

    file.close();
    return true;
  }
  else {

    logERROR("Analyzer - failed opening file: " + fullFileName);
    return false;
  }
}

//
// Is starting triplet from different layers (avoid using overlapping modules in one layer)
//
bool Analyzer::isTripletFromDifLayers(TrackNew& track, int iHit, bool propagOutIn) {

  std::map<std::string, bool> hitIDs;

  for (int i=0; i<=iHit; i++) {

    std::string hitID = "";

    if (!propagOutIn) hitID = track.getMeasurableOrIPHit(i)->getDetName() +
                              std::string(track.getMeasurableOrIPHit(i)->isBarrel() ? "_L_" : "_D_") +
                              any2str(track.getMeasurableOrIPHit(i)->getLayerOrDiscID());
    else              hitID = track.getRMeasurableOrIPHit(i)->getDetName() +
                              std::string(track.getRMeasurableOrIPHit(i)->isBarrel() ? "_L_" : "_D_") +
                              any2str(track.getRMeasurableOrIPHit(i)->getLayerOrDiscID());

    hitIDs[hitID] = true;
  }

  if (hitIDs.size()>=3) return true;
  else                  return false;
}

//
// Analyze tracker pattern recognition capabilities. The quantities used to qualify the pattern reco capabilities are meant for primary tracks only
// rather than for secondary tracks, as one always starts in the innermost layer (+ IP contraint) or outermost layer in case of out-in approach.
// The extrapolation part of pattern recognition uses formulae calculated in parabolic approximation and particle fluences (occupancies) scaled to
// given pile-up scenario. The fluences are read in from an external file and are assumed to be simulated by Fluka for given p-p collision energies
// & rough detector (material) setup.
//
bool Analyzer::analyzePatterReco(MaterialBudget& mb, mainConfigHandler& mainConfig, int nTracks, MaterialBudget* pm) {

  bool isAnalysisOK = true;

  // Get fluence map
  const auto& directory     = mainConfig.getIrradiationDirectory();
  bool        fluenceMapOK  = false;
  IrradiationMap* fluenceMap= nullptr;
  int         nBins         = vis_n_bins*2;

  fluenceMapOK = checkFile(default_fluence_file, directory);
  std::cout << "Reading in: " << directory + "/" + default_fluence_file << std::endl;
  if (fluenceMapOK) fluenceMap  = new IrradiationMap(directory + "/" + default_fluence_file);

  // Set nTracks
  if (nTracks<=0) {

    std::ostringstream message;
    message << "PatternReco: Number of simulation tracks zero or negative: " << nTracks << " or zero number of trackers to be initialized";
    logERROR(message.str());
  }
  else {

    // Prepare histograms
    for (const auto& pIter : mainConfig.getMomenta()) {

      // In-Out
      std::string name = "hisBkgContInOut_pT"+any2str(pIter/Units::GeV);
      hisPatternRecoInOutPt.push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));
      name             = "hisBkgContInOut_p"+any2str(pIter/Units::GeV);
      hisPatternRecoInOutP.push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));

      // Out-In
      name             = "hisBkgContOutIn_pT"+any2str(pIter/Units::GeV);
      hisPatternRecoOutInPt.push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));
      name             = "hisBkgContOutIn_p"+any2str(pIter/Units::GeV);
      hisPatternRecoOutInP.push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));
    }
  }

  // Check that map available
  if (!fluenceMapOK) isAnalysisOK = false;
  else for (int iTrack = 0; iTrack <nTracks; iTrack++) {

    // Define track
    TrackNew matTrack;

    double eta   = 0.0 + geom_max_eta_coverage/nTracks*(iTrack+0.5);
    double theta = 2 * atan(exp(-eta));
    double phi   = myDice.Rndm() * M_PI * 2.0;
    double pT    = 100*Units::TeV; // Arbitrarily high number

    matTrack.setThetaPhiPt(theta, phi, pT);
    matTrack.setOrigin(0, 0, 0); // TODO: Not assuming z-error when analyzing resolution (missing implementation of non-zero track starting point in inactive hits)

    // Assign material to the track
    findAllHits(mb, pm, matTrack);

    // Add beam-pipe
    double rPos  = 23.*Units::mm;
    double zPos  = rPos/tan(theta);

    HitNewPtr hit(new HitNew(rPos, zPos));
    hit->setAsPassive();

    Material material;
    material.radiation   = 0.0022761 / sin(theta);  // was 0.0023, adapted to fit CMSSW 81X 2016/11/30
    material.interaction = 0.0020334 / sin(theta);  // was 0.0019, adapted to fit CMSSW 81X 2016/11/30

    hit->setCorrectedMaterial(material);
    hit->setBeamPipe(true);
    matTrack.addHit(std::move(hit));

    // For each momentum/transverse momentum compute
    int iMomentum = 0;
    for (const auto& pIter : mainConfig.getMomenta()) {

      // 2 options: pT & p
      for (int pOption=0;pOption<2;pOption++) {

        // InOut or OutIn approach
        for (int approachOption=0; approachOption<2; approachOption++) {

          bool propagOutIn = approachOption;

          //if (!propagOutIn) std::cout << "InOut approach" << std::endl;
          //else              std::cout << "OutIn approach" << std::endl;

          // Reset track total probability & calculate p/pT based on option
          double pNotContamTot      = 1;
          double pNotContamTotInner = 1;

          if (pOption==0) pT = pIter;             // pT option
          else            pT = pIter*sin(theta);  // p option

          // Set track & prune hits
          TrackNew track(matTrack);
          track.resetPt(pT);

          //
          // Pattern recognition procedure
          bool useIP        = true;
          int nMeasuredHits = 0;
          if (!propagOutIn) nMeasuredHits = track.getNMeasuredHits("all", !useIP);
          else              nMeasuredHits = track.getNMeasuredHits("all", !useIP);

          //track.printHits();
          //std::cout << ">> " << nMeasuredHits << " dD0=" << track.getDeltaD(0.0)/Units::um << " dZ0=" << track.getDeltaZ(0.0)/Units::um << std::endl;

          // Start with first 3 hits and end with N-1 hits (C counting from zero)
          bool testTriplet = true;
          for (int iHit=2; iHit<nMeasuredHits-1; iHit++) {

            // Start with triplet coming from 3 different layers -> avoid using 2 close measurements from overlapping modules in 1 layer
            if (testTriplet && !isTripletFromDifLayers(track, iHit, propagOutIn)) continue;
            else testTriplet = false;

            int nHitsUsed = iHit+1; // (C counting from zero)

            // Keep first/last N hits (based on approach) and find paramaters of the next hit
            double nextRPos    = 0;
            double nextZPos    = 0;
            double nextHitTilt = 0;
            double A           = 0; // r_i/2R
            std::string iHitID = "";

            if (!propagOutIn) {

              nextRPos    = track.getMeasurableOrIPHit(iHit+1)->getRPos();
              nextZPos    = track.getMeasurableOrIPHit(iHit+1)->getZPos();
              nextHitTilt = track.getMeasurableOrIPHit(iHit+1)->getTilt();
              if (SimParms::getInstance().isMagFieldConst()) A = track.getMeasurableOrIPHit(iHit+1)->getRPos()/2./track.getRadius(track.getMeasurableOrIPHit(iHit+1)->getZPos());  // r_i/2R

              iHitID      = track.getMeasurableOrIPHit(iHit+1)->getDetName() +
                            std::string(track.getMeasurableOrIPHit(iHit+1)->isBarrel() ? "_L_" : "_D_") +
                            any2str(track.getMeasurableOrIPHit(iHit+1)->getLayerOrDiscID());
            }
            else {

              nextRPos    = track.getRMeasurableOrIPHit(iHit+1)->getRPos();
              nextZPos    = track.getRMeasurableOrIPHit(iHit+1)->getZPos();
              nextHitTilt = track.getRMeasurableOrIPHit(iHit+1)->getTilt();
              if (SimParms::getInstance().isMagFieldConst()) A = track.getRMeasurableOrIPHit(iHit+1)->getRPos()/2./track.getRadius(track.getRMeasurableOrIPHit(iHit+1)->getZPos());  // r_i/2R
              iHitID      = track.getRMeasurableOrIPHit(iHit+1)->getDetName() +
                            std::string(track.getRMeasurableOrIPHit(iHit+1)->isBarrel() ? "_L_" : "_D_") +
                            any2str(track.getRMeasurableOrIPHit(iHit+1)->getLayerOrDiscID());
            }

            // Calculate dRPhi & dZ
            int    nPU    = SimParms::getInstance().numMinBiasEvents();
            // TODO: Assumed that map on ZxR grid defined in cm -> needs to be fixed
            double flux   = fluenceMap->calculateIrradiation(std::make_pair(nextZPos/Units::mm,nextRPos/Units::mm))*1/Units::cm2/Units::s * nPU;

            double corrFactor = cos(nextHitTilt) + track.getCotgTheta()/sqrt(1-A*A)*sin(nextHitTilt);
            double dZProj     = track.getDeltaZ(nextRPos,propagOutIn)/corrFactor;
            double dDProj     = track.getDeltaD(nextRPos,propagOutIn);

            // Print info
            //std::cout << ">> " << iHitID << " R=" << nextRPos/Units::mm << " dD=" << dDProj/Units::um << " " << " Z=" << nextZPos/Units::mm << " dZ=" << dZProj/Units::um << std::endl;

            // Calculate how many sigmas does one need to get in 2D Gauss. 5% coverge
            // F(mu + n*sigma) - F(mu - n*sigma) = erf(n/sqrt(2))
            // In 2D we assume independent measurement in r-phi & Z, hence 0.95 = erf(n/sqrt(2))*erf(n/sqrt(2)) assuming the same number of sigmas (n) in both r-phi & z
            // Hence n = InverseErf(sqrt(0.95))*sqrt(2)
            static double nSigmaFactor = TMath::ErfInverse(sqrt(0.95))*sqrt(2);

            // Calculate errorEllipse multiplied by nSigmaFactor(RPhi)*nSigmaFactor(Z)
            double errorEllipse = M_PI*dZProj*dDProj*nSigmaFactor*nSigmaFactor;

            // Accumulate probibilites accross all hits not to be contaminated anywhere -> final probability of being contaminated is 1 - pContam
            double pContam = flux*errorEllipse;
            if (pContam>1) pContam = 1.0;
            pNotContamTot *= 1-pContam;

            //
            // Fill d0proj, z0proj & pContam for individual layers/discs to histograms

            // Probability of inner tracker only
            if (iHitID.find("Inner")!=std::string::npos) pNotContamTotInner *= 1-pContam;

            // Create artificial map identifier with extra characters to have innermost first, then outermost & then fwd
            std::string iHitIDMap = "";
            if      (iHitID.find("Inner")!=std::string::npos) iHitIDMap = "A_" + iHitID;
            else if (iHitID.find("Outer")!=std::string::npos) iHitIDMap = "B_" + iHitID;
            else                                              iHitIDMap = "C_" + iHitID;

            // Pt & InOut
            if (pOption==0 && approachOption==0) {

              // Create profile histograms if don't exist yet
              if (hisPtHitDProjInOut.find(iHitIDMap)==hisPtHitDProjInOut.end()) {
                for (int iMom=0; iMom<mainConfig.getMomenta().size(); iMom++) {

                  std::string name = "hisPt_" + any2str(iMom) + "_Hit_" + iHitID + "_DProjInOut";
                  hisPtHitDProjInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));

                  name = "hisPt" + any2str(iMom) + "_Hit_" + iHitID + "_ZProjInOut";
                  hisPtHitZProjInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));

                  name = "hisPt" + any2str(iMom) + "_Hit_" + iHitID + "_ProbContamInOut";
                  hisPtHitProbContamInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));
                }
              }
              hisPtHitDProjInOut[iHitIDMap][iMomentum]->Fill(eta, dDProj/Units::um);
              hisPtHitZProjInOut[iHitIDMap][iMomentum]->Fill(eta, dZProj/Units::um);
              hisPtHitProbContamInOut[iHitIDMap][iMomentum]->Fill(eta, pContam);
            }
            // P & InOut
            if (pOption==1 && approachOption==0) {

              // Create profile histograms if don't exist yet
              if (hisPHitDProjInOut.find(iHitIDMap)==hisPHitDProjInOut.end()) {
                for (int iMom=0; iMom<mainConfig.getMomenta().size(); iMom++) {

                  std::string name = "hisP_" + any2str(iMom) + "_Hit_" + iHitID + "_DProjInOut";
                  hisPHitDProjInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));

                  name = "hisP" + any2str(iMom) + "_Hit_" + iHitID + "_ZProjInOut";
                  hisPHitZProjInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));

                  name = "hisP" + any2str(iMom) + "_Hit_" + iHitID + "_ProbContamInOut";
                  hisPHitProbContamInOut[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));
                }
              }
              // TODO: Check that above zero!!!
              hisPHitDProjInOut[iHitIDMap][iMomentum]->Fill(eta, dDProj/Units::um);
              hisPHitZProjInOut[iHitIDMap][iMomentum]->Fill(eta, dZProj/Units::um);
              hisPHitProbContamInOut[iHitIDMap][iMomentum]->Fill(eta, pContam);
            }
            // Pt & OutIn
            if (pOption==0 && approachOption==1) {

              // Create profile histograms if don't exist yet
              if (hisPtHitDProjOutIn.find(iHitIDMap)==hisPtHitDProjOutIn.end()) {
                for (int iMom=0; iMom<mainConfig.getMomenta().size(); iMom++) {

                  std::string name = "hisPt_" + any2str(iMom) + "_Hit_" + iHitID + "_DProjOutIn";
                  hisPtHitDProjOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));

                  name = "hisPt" + any2str(iMom) + "_Hit_" + iHitID + "_ZProjOutIn";
                  hisPtHitZProjOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));

                  name = "hisPt" + any2str(iMom) + "_Hit_" + iHitID + "_ProbContamOutIn";
                  hisPtHitProbContamOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));
                }
              }
              hisPtHitDProjOutIn[iHitIDMap][iMomentum]->Fill(eta, dDProj/Units::um);
              hisPtHitZProjOutIn[iHitIDMap][iMomentum]->Fill(eta, dZProj/Units::um);
              hisPtHitProbContamOutIn[iHitIDMap][iMomentum]->Fill(eta, pContam);
            }
            // P & OutIn
            if (pOption==1 && approachOption==1) {

              // Create profile histograms if don't exist yet
              if (hisPHitDProjOutIn.find(iHitIDMap)==hisPHitDProjOutIn.end()) {
                for (int iMom=0; iMom<mainConfig.getMomenta().size(); iMom++) {

                  std::string name = "hisP_" + any2str(iMom) + "_Hit_" + iHitID + "_DProjOutIn";
                  hisPHitDProjOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));

                  name = "hisP" + any2str(iMom) + "_Hit_" + iHitID + "_ZProjOutIn";
                  hisPHitZProjOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));

                  name = "hisP" + any2str(iMom) + "_Hit_" + iHitID + "_ProbContamOutIn";
                  hisPHitProbContamOutIn[iHitIDMap].push_back(new TProfile(name.c_str(),name.c_str(),nBins, 0, geom_max_eta_coverage));
                }
              }
              // TODO: Check that above zero!!!
              hisPHitDProjOutIn[iHitIDMap][iMomentum]->Fill(eta, dDProj/Units::um);
              hisPHitZProjOutIn[iHitIDMap][iMomentum]->Fill(eta, dZProj/Units::um);
              hisPHitProbContamOutIn[iHitIDMap][iMomentum]->Fill(eta, pContam);
            }
          } // Pattern reco loop

          //
          // Calculate total contamination based on different options: p/pT, in-out/out-in
          double pContamTot      = 1-pNotContamTot;
          double pContamInnerTot = 1-pNotContamTotInner;
          if (pOption==0) {

            if (approachOption==0) {
              hisPatternRecoInOutPt[iMomentum]->Fill(eta,pContamTot);
            }
            else {
              hisPatternRecoOutInPt[iMomentum]->Fill(eta,pContamTot);
            }
          }
          else  {

            if (approachOption==0) {
              hisPatternRecoInOutP[iMomentum]->Fill(eta,pContamTot);
            }
            else {
              hisPatternRecoOutInP[iMomentum]->Fill(eta,pContamTot);
            }
          }

        } //In-Out, Out-In approach - options loop
      } //pT/p option loop

      iMomentum++;
    } // Momenta loop

  } // Tracks loop

  return isAnalysisOK;
}

void Analyzer::fillAvailableSpacing(Tracker& tracker, std::vector<double>& spacingOptions) {
  double aSpacing;
  std::set<double> foundSpacing;

  // Loop over all the layers
  for (auto& aModule : tracker.modules()) {
    if (aModule->sensorLayout() == PT) {
      aSpacing = aModule->dsDistance();
      if (aSpacing>0) {
        foundSpacing.insert(aSpacing);
      }
    }
  }

  spacingOptions.clear();
  for(std::set<double>::iterator it=foundSpacing.begin(); it!=foundSpacing.end(); ++it) {
    spacingOptions.push_back(*it);
  }

}

void Analyzer::createTriggerDistanceTuningPlots(Tracker& tracker, const std::vector<double>& triggerMomenta) {
  TriggerDistanceTuningPlotsVisitor v(myProfileBag, 
                                      triggerMomenta);
  SimParms::getInstance().accept(v);
  tracker.accept(v);
  v.postVisit();
  optimalSpacingDistribution = v.optimalSpacingDistribution;
  optimalSpacingDistributionAW = v.optimalSpacingDistributionAW; 
  triggerRangeLowLimit = v.triggerRangeLowLimit; 
  triggerRangeHighLimit = v.triggerRangeHighLimit; 
  spacingTuningFrame = v.spacingTuningFrame;
  spacingTuningGraphs = v.spacingTuningGraphs;
  spacingTuningGraphsBad = v.spacingTuningGraphsBad;
  moduleOptimalSpacings = v.moduleOptimalSpacings;
}






void Analyzer::fillTriggerEfficiencyGraphs(const Tracker& tracker,
                                           const std::vector<double>& triggerMomenta,
                                           const std::vector<Track>& trackVector) {

  // Prepare the graphs to record the number of triggered points
  //std::map<double, TGraph>& trigGraphs = myGraphBag.getGraphs(GraphBag::TriggerGraph|GraphBag::TriggeredGraph);
  std::map<double, TProfile>& trigProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggeredProfile);
  std::map<double, TProfile>& trigFractionProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggeredFractionProfile);
  std::map<double, TProfile>& trigPurityProfiles = myProfileBag.getProfiles(profileBag::TriggerProfile|profileBag::TriggerPurityProfile);

  TProfile& totalProfile = trigProfiles[profileBag::Triggerable];

  std::map<std::string, std::map<std::string, TH1I*>>& stubEfficiencyCoverageProfiles = getStubEfficiencyCoverageProfiles();

  double maxEta = 4.0; //getEtaMaxTrigger();

  for (std::vector<Track>::const_iterator itTrack = trackVector.begin();
       itTrack != trackVector.end(); ++itTrack) {
    const Track& myTrack=(*itTrack);

    double eta = myTrack.getEta();
    int nHits = myTrack.nActiveHits(false, false);
    totalProfile.Fill(eta, nHits);
    std::vector<std::pair<Module*,HitType>> hitModules = myTrack.getHitModules();

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
           double bgReductionFactor; // Reduction of the combinatorial background for ptPS modules by turning off the appropriate pixels
           for (const auto& modAndType : hitModules) {
             Module* hitModule = modAndType.first;
             PtErrorAdapter pterr(*hitModule);
             // Hits that we would like to have from tracks above this threshold
             curAvgInteresting += pterr.getParticleFrequencyPerEventAbove(*itMomentum);
             // ... out of which we only see these
             curAvgTrue += pterr.getTriggerFrequencyTruePerEventAbove(*itMomentum);
               
             // The background is given by the contamination from low pT tracks...
             curAvgFake += pterr.getTriggerFrequencyTruePerEventBelow(*itMomentum);
             // ... plus the combinatorial background from occupancy (can be reduced using ptPS modules)
             if (hitModule->reduceCombinatorialBackground()) bgReductionFactor = hitModule->geometricEfficiency(); else bgReductionFactor=1;
             curAvgFake += pterr.getTriggerFrequencyFakePerEvent()*SimParms::getInstance().numMinBiasEvents() * bgReductionFactor;

             std::string layerName = hitModule->uniRef().cnt + "_" + any2str(hitModule->uniRef().layer);
             if (modAndType.second == HitType::STUB) {
               std::string momentumString = any2str(*itMomentum, 2);
               if (stubEfficiencyCoverageProfiles[layerName].count(momentumString) == 0) {
                 stubEfficiencyCoverageProfiles[layerName][momentumString] = new TH1I(Form("stubEfficiencyCoverageProfile%s%s", layerName.c_str(), momentumString.c_str()), (layerName + ";#eta;Stubs").c_str(), trackVector.size(), 0.0, maxEta); 
               }
               stubEfficiencyCoverageProfiles[layerName][momentumString]->Fill(myTrack.getEta(), 1);
             } 
           }
           myPurityProfile.Fill(eta, 100*curAvgTrue/(curAvgTrue+curAvgFake));
        }
      }
    }
  }
//  for (auto i : stubEfficiencyCoverageProfiles) {
//    std::cout << "--------------------- " << i.first << " ------------------ " << std::endl;
//    for (auto j : i.second) {
//      std::cout << j.first << std::endl;
//      j.second->Print("all");
//    }
//  }
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
                                     MaterialBudget* pm) {

  Tracker& tracker = mb.getTracker();
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
  setCellBoundaries(nTracks, 0.0, geom_max_radius + geom_inactive_volume_width, 0.0, getEtaMaxMaterial());

  // reset the list of tracks
  // std::vector<Track> tv;
  // std::vector<Track> tvIdeal;

  for (int i_eta = 0; i_eta < nTracks; i_eta++) {
    phi = myDice.Rndm() * M_PI * 2.0;
    Material tmp;
    Track track;
    eta = i_eta * etaStep;
    theta = 2 * atan(exp(-eta)); // TODO: switch to exp() here
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

    std::map<std::string, Material> sumServicesComponentsRI;

    //      services, barrel
    tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getBarrelServices(), eta, theta, track, sumServicesComponentsRI, MaterialProperties::no_cat);
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
    tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getEndcapServices(), eta, theta, track, sumServicesComponentsRI, MaterialProperties::no_cat);
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


    /*for (std::map<std::string, Material>::iterator it = sumServicesComponentsRI.begin(); it != sumServicesComponentsRI.end(); ++it) {
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
      }*/



    //      supports, barrel
    tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, sumServicesComponentsRI, MaterialProperties::b_sup);
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
    tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, sumServicesComponentsRI, MaterialProperties::e_sup);
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
    tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, sumServicesComponentsRI, MaterialProperties::o_sup);
    rlazytube.Fill(eta, tmp.radiation);
    ilazytube.Fill(eta, tmp.interaction);
    rlazyall.Fill(eta, tmp.radiation);
    ilazyall.Fill(eta, tmp.interaction);
    rglobal.Fill(eta, tmp.radiation);
    iglobal.Fill(eta, tmp.interaction);
    rComponents["Supports"]->Fill(eta, tmp.radiation);
    iComponents["Supports"]->Fill(eta, tmp.interaction);
    //      supports, barrel tubes
    tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, sumServicesComponentsRI, MaterialProperties::t_sup);
    rlazybtube.Fill(eta, tmp.radiation);
    ilazybtube.Fill(eta, tmp.interaction);
    rlazyall.Fill(eta, tmp.radiation);
    ilazyall.Fill(eta, tmp.interaction);
    rglobal.Fill(eta, tmp.radiation);
    iglobal.Fill(eta, tmp.interaction);
    rComponents["Supports"]->Fill(eta, tmp.radiation);
    iComponents["Supports"]->Fill(eta, tmp.interaction);
    //      supports, user defined
    tmp = analyzeInactiveSurfaces(mb.getInactiveSurfaces().getSupports(), eta, theta, track, sumServicesComponentsRI, MaterialProperties::u_sup);
    rlazyuserdef.Fill(eta, tmp.radiation);
    ilazyuserdef.Fill(eta, tmp.interaction);
    rlazyall.Fill(eta, tmp.radiation);
    ilazyall.Fill(eta, tmp.interaction);
    rglobal.Fill(eta, tmp.radiation);
    iglobal.Fill(eta, tmp.interaction);
    rComponents["Supports"]->Fill(eta, tmp.radiation);
    iComponents["Supports"]->Fill(eta, tmp.interaction);
    //      pixels, if they exist
    std::map<std::string, Material> ignoredPixelSumComponentsRI;
    std::map<std::string, Material> ignoredPixelSumServicesComponentsRI;
    if (pm != NULL) {
      analyzeModules(pm->getBarrelModuleCaps(), eta, theta, phi, track, ignoredPixelSumComponentsRI, true);
      analyzeModules(pm->getEndcapModuleCaps(), eta, theta, phi, track, ignoredPixelSumComponentsRI, true);
      analyzeInactiveSurfaces(pm->getInactiveSurfaces().getBarrelServices(), eta, theta, track, ignoredPixelSumServicesComponentsRI, MaterialProperties::no_cat, true);
      analyzeInactiveSurfaces(pm->getInactiveSurfaces().getEndcapServices(), eta, theta, track, ignoredPixelSumServicesComponentsRI, MaterialProperties::no_cat, true);
      analyzeInactiveSurfaces(pm->getInactiveSurfaces().getSupports(), eta, theta, track, ignoredPixelSumServicesComponentsRI, MaterialProperties::b_sup, true);
      }

    // Add the hit on the beam pipe
    Hit* hit = new Hit(23./sin(theta));
    hit->setOrientation(Hit::Horizontal);
    hit->setObjectKind(Hit::Inactive);
    hit->setObjectCategory(Hit::BeamPipe);
    Material beamPipeMat;
    beamPipeMat.radiation = 0.0022761 / sin(theta);  // was 0.0023, adapted to fit CMSSW 81X 2016/11/30
    beamPipeMat.interaction = 0.0020334 / sin(theta);  // was 0.0019, adapted to fit CMSSW 81X 2016/11/30
    hit->setCorrectedMaterial(beamPipeMat);
    track.addHit(hit);
    if (!track.noHits()) {
      track.sort();
      track.addEfficiency();
      track.addEfficiency();

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






    if (eta >= 0.) {

      for (const auto& hit : track.getHitV()) {
	// SERVICES DETAILS
	if (!hit->isPixel() && hit->getObjectCategory() == Hit::Service) {
	  
	  InactiveElement* inactive = hit->getHitInactiveElement();
	  std::map<std::string, Material> servicesComponentsRI = inactive->getComponentsRI();
	  for (const auto& it : servicesComponentsRI) {
	    Material res;
	    res.radiation = it.second.radiation / (inactive->isVertical() ? cos(theta) : sin(theta));  
	    res.interaction = it.second.interaction / (inactive->isVertical() ? cos(theta) : sin(theta));
	    if (rComponentsServicesDetails[it.first]==NULL) {
	      rComponentsServicesDetails[it.first] = new TH1D();
	      rComponentsServicesDetails[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	    }
	    rComponentsServicesDetails[it.first]->Fill(eta, res.radiation);
	    if (iComponentsServicesDetails[it.first]==NULL) { 
	      iComponentsServicesDetails[it.first] = new TH1D();
	      iComponentsServicesDetails[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	    }
	    iComponentsServicesDetails[it.first]->Fill(eta, res.interaction);
	  }

	  /*if (servicesComponentsRI.size() == 0) {
	    if (rComponentsServicesDetails["Services : others"]==NULL) {
	      rComponentsServicesDetails["Services : others"] = new TH1D();
	      rComponentsServicesDetails["Services : others"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	    }
	    rComponentsServicesDetails["Services : others"]->Fill(eta, hit->getCorrectedMaterial().radiation);
	    if (iComponentsServicesDetails["Services : others"]==NULL) {
	      iComponentsServicesDetails["Services : others"] = new TH1D();
	      iComponentsServicesDetails["Services : others"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	    }
	    iComponentsServicesDetails["Services : others"]->Fill(eta, hit->getCorrectedMaterial().interaction);
	    }*/
	}
      }

	  }





























    
    if (mb.getTracker().myid() == "Outer" && eta >= 0.) {
      track.assignTrackingVolumesToHits();


      
      for (const auto& hit : track.getHitV()) {
	if (hit->isTotalTrackingVolume() && hit->getObjectCategory() == Hit::BeamPipe) {
	  if (rComponentsBeamPipe["Beam Pipe"]==NULL) {
	    rComponentsBeamPipe["Beam Pipe"] = new TH1D();
	    rComponentsBeamPipe["Beam Pipe"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  rComponentsBeamPipe["Beam Pipe"]->Fill(eta, hit->getCorrectedMaterial().radiation);
	  if (iComponentsBeamPipe["Beam Pipe"]==NULL) { 
	    iComponentsBeamPipe["Beam Pipe"] = new TH1D();
	    iComponentsBeamPipe["Beam Pipe"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  iComponentsBeamPipe["Beam Pipe"]->Fill(eta, hit->getCorrectedMaterial().interaction);
	}

	if (hit->isPixelIntersticeVolume() && hit->getObjectCategory() != Hit::BeamPipe) {
	  if (rComponentsPixelInterstice["Services under Pixel Tracking Volume"]==NULL) {
	    rComponentsPixelInterstice["Services under Pixel Tracking Volume"] = new TH1D();
	    rComponentsPixelInterstice["Services under Pixel Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  rComponentsPixelInterstice["Services under Pixel Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().radiation);
	  if (iComponentsPixelInterstice["Services under Pixel Tracking Volume"]==NULL) { 
	    iComponentsPixelInterstice["Services under Pixel Tracking Volume"] = new TH1D();
	    iComponentsPixelInterstice["Services under Pixel Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  iComponentsPixelInterstice["Services under Pixel Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().interaction);
	}
	}
    



     
      for (const auto& it : ignoredPixelSumComponentsRI) {
	if (rComponentsPixelTrackingVolume[it.first]==NULL) {
	  rComponentsPixelTrackingVolume[it.first] = new TH1D();
	  rComponentsPixelTrackingVolume[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	}
	rComponentsPixelTrackingVolume[it.first]->Fill(eta, it.second.radiation);
	if (iComponentsPixelTrackingVolume[it.first]==NULL) {
	  iComponentsPixelTrackingVolume[it.first] = new TH1D();
	  iComponentsPixelTrackingVolume[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	}
	iComponentsPixelTrackingVolume[it.first]->Fill(eta, it.second.interaction);
      }









      // WRONG
      /*for (const auto& it : ignoredPixelSumServicesComponentsRI) {
	if (rComponentsPixelTrackingVolume[it.first]==NULL) {
	  rComponentsPixelTrackingVolume[it.first] = new TH1D();
	  rComponentsPixelTrackingVolume[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	}
	rComponentsPixelTrackingVolume[it.first]->Fill(eta, it.second.radiation);
	if (iComponentsPixelTrackingVolume[it.first]==NULL) {
	  iComponentsPixelTrackingVolume[it.first] = new TH1D();
	  iComponentsPixelTrackingVolume[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	}
	iComponentsPixelTrackingVolume[it.first]->Fill(eta, it.second.interaction);
	}*/


      for (const auto& hit : track.getHitV()) {

	if (hit->isPixelTrackingVolume() && hit->getObjectCategory() == Hit::Service) {
	  if (rComponentsPixelTrackingVolume["Services in Pixel Tracking Volume"]==NULL) {	   
	    rComponentsPixelTrackingVolume["Services in Pixel Tracking Volume"] = new TH1D();
	    rComponentsPixelTrackingVolume["Services in Pixel Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  rComponentsPixelTrackingVolume["Services in Pixel Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().radiation);
	  if (iComponentsPixelTrackingVolume["Services in Pixel Tracking Volume"]==NULL) { 
	    iComponentsPixelTrackingVolume["Services in Pixel Tracking Volume"] = new TH1D();
	    iComponentsPixelTrackingVolume["Services in Pixel Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  iComponentsPixelTrackingVolume["Services in Pixel Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().interaction);
	  }



	if (hit->isPixelTrackingVolume() && hit->getObjectCategory() == Hit::Support) {
	  if (rComponentsPixelTrackingVolume["Supports in Pixel Tracking Volume"]==NULL) {
	    rComponentsPixelTrackingVolume["Supports in Pixel Tracking Volume"] = new TH1D();
	    rComponentsPixelTrackingVolume["Supports in Pixel Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  rComponentsPixelTrackingVolume["Supports in Pixel Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().radiation);
	  if (iComponentsPixelTrackingVolume["Supports in Pixel Tracking Volume"]==NULL) { 
	    iComponentsPixelTrackingVolume["Supports in Pixel Tracking Volume"] = new TH1D();
	    iComponentsPixelTrackingVolume["Supports in Pixel Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  iComponentsPixelTrackingVolume["Supports in Pixel Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().interaction);
	}


	if (hit->isIntersticeVolume()) {
	  if (rComponentsInterstice["Services and supports in interstice"]==NULL) {
	    rComponentsInterstice["Services and supports in interstice"] = new TH1D();
	    rComponentsInterstice["Services and supports in interstice"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  rComponentsInterstice["Services and supports in interstice"]->Fill(eta, hit->getCorrectedMaterial().radiation);
	  if (iComponentsInterstice["Services and supports in interstice"]==NULL) { 
	    iComponentsInterstice["Services and supports in interstice"] = new TH1D();
	    iComponentsInterstice["Services and supports in interstice"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  iComponentsInterstice["Services and supports in interstice"]->Fill(eta, hit->getCorrectedMaterial().interaction);
	}
      }




      for (const auto& it : sumComponentsRI) {
	if (rComponentsOuterTrackingVolume[it.first]==NULL) {
	  rComponentsOuterTrackingVolume[it.first] = new TH1D();
	  rComponentsOuterTrackingVolume[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	}
	rComponentsOuterTrackingVolume[it.first]->Fill(eta, it.second.radiation);
	if (iComponentsOuterTrackingVolume[it.first]==NULL) {
	  iComponentsOuterTrackingVolume[it.first] = new TH1D();
	  iComponentsOuterTrackingVolume[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	}
	iComponentsOuterTrackingVolume[it.first]->Fill(eta, it.second.interaction);
	}

      // WRONG
      /*for (const auto& it : sumServicesComponentsRI) {
	if (rComponentsOuterTrackingVolume[it.first]==NULL) {
	  rComponentsOuterTrackingVolume[it.first] = new TH1D();
	  rComponentsOuterTrackingVolume[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	}
	rComponentsOuterTrackingVolume[it.first]->Fill(eta, it.second.radiation);
	if (iComponentsOuterTrackingVolume[it.first]==NULL) {
	  iComponentsOuterTrackingVolume[it.first] = new TH1D();
	  iComponentsOuterTrackingVolume[it.first]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	}
	iComponentsOuterTrackingVolume[it.first]->Fill(eta, it.second.interaction);
	}*/
 

      for (const auto& hit : track.getHitV()) {

	if (hit->isOuterTrackingVolume() && hit->getObjectCategory() == Hit::Service) {
	  if (rComponentsOuterTrackingVolume["Services in Outer Tracking Volume"]==NULL) {	   
	    rComponentsOuterTrackingVolume["Services in Outer Tracking Volume"] = new TH1D();
	    rComponentsOuterTrackingVolume["Services in Outer Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  rComponentsOuterTrackingVolume["Services in Outer Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().radiation);
	  if (iComponentsOuterTrackingVolume["Services in Outer Tracking Volume"]==NULL) { 
	    iComponentsOuterTrackingVolume["Services in Outer Tracking Volume"] = new TH1D();
	    iComponentsOuterTrackingVolume["Services in Outer Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  iComponentsOuterTrackingVolume["Services in Outer Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().interaction);
	  }


	if (hit->isOuterTrackingVolume() && hit->getObjectCategory() == Hit::Support) {
	  if (rComponentsOuterTrackingVolume["Supports in Outer Tracking Volume"]==NULL) {	    
	    rComponentsOuterTrackingVolume["Supports in Outer Tracking Volume"] = new TH1D();
	    rComponentsOuterTrackingVolume["Supports in Outer Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  rComponentsOuterTrackingVolume["Supports in Outer Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().radiation);
	  if (iComponentsOuterTrackingVolume["Supports in Outer Tracking Volume"]==NULL) { 
	    iComponentsOuterTrackingVolume["Supports in Outer Tracking Volume"] = new TH1D();
	    iComponentsOuterTrackingVolume["Supports in Outer Tracking Volume"]->SetBins(nTracks, 0.0, getEtaMaxMaterial()); 
	  }
	  iComponentsOuterTrackingVolume["Supports in Outer Tracking Volume"]->Fill(eta, hit->getCorrectedMaterial().interaction);
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

}





void Analyzer::computeTriggerFrequency(Tracker& tracker) {
  TriggerFrequencyVisitor v; 
  SimParms::getInstance().accept(v);
  tracker.accept(v);

  triggerFrequencyTrueSummaries_ = v.triggerFrequencyTrueSummaries;
  triggerFrequencyFakeSummaries_ = v.triggerFrequencyFakeSummaries;
  triggerFrequencyMisfilteredSummaries_ = v.triggerFrequencyMisfilteredSummaries;
  triggerFrequencyCombinatorialSummaries_ = v.triggerFrequencyCombinatorialSummaries;
  triggerFrequencyInterestingSummaries_ = v.triggerFrequencyInterestingSummaries;
  triggerRateSummaries_ = v.triggerRateSummaries; 
  triggerEfficiencySummaries_ = v.triggerEfficiencySummaries; 
  triggerPuritySummaries_ = v.triggerPuritySummaries; 
  triggerDataBandwidthSummaries_ = v.triggerDataBandwidthSummaries;
  triggerFrequenciesPerEvent_ = v.triggerFrequenciesPerEvent;
  stripOccupancySummaries_ = v.stripOccupancySummaries;
  hitOccupancySummaries_ = v.hitOccupancySummaries;
}


  
    
void Analyzer::computeTriggerProcessorsBandwidth(Tracker& tracker) {
  TriggerProcessorBandwidthVisitor v(triggerDataBandwidths_, triggerFrequenciesPerEvent_);
  v.preVisit();
  SimParms::getInstance().accept(v);
  tracker.accept(v);
  v.postVisit();

  processorConnectionSummary_ = v.processorConnectionSummary; 
  processorCommonConnectionSummary_ = v.processorCommonConnectionSummary;
  processorInboundBandwidthSummary_ = v.processorInboundBandwidthSummary; 
  processorInboundStubPerEventSummary_ = v.processorInboundStubPerEventSummary; 
  moduleConnectionsDistribution = v.moduleConnectionsDistribution; 
  moduleConnections_ = v.moduleConnections;
  triggerSectorMap_ = v.sectorMap;
  processorCommonConnectionMap_ = v.processorCommonConnectionMap;
  sampleTriggerPetal_ = v.sampleTriggerPetal;
  triggerPetalCrossoverR_ = v.crossoverR;
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

  std::map<std::string, double>::const_iterator localmassesBegin;
  std::map<std::string, double>::const_iterator localmassesEnd;

  // First create a list of material used anywhere
  std::vector<std::string> materialTagV;
  std::vector<std::string>::iterator materialTagIt;

  double localMaterial;
  string materialTag;

  // loop over layers
  for (layerIt = tracker.begin(); layerIt != tracker.end(); ++layerIt) {
    // Loop over modules
    for (moduleIt = layerIt->begin(); moduleIt != layerIt->end(); ++moduleIt) {
      // Check if the module is on the YZ section
      myModuleCap = &(*moduleIt);
      myModule = &(myModuleCap->getModule());
      if (myModule->posRef().phi==1) {
        pair<int, int> myIndex = make_pair(myModule->tableRef().row, myModule->tableRef().col);
        tempString = myModule->cntName();
        if (!typeTaken[tempString][myIndex]) {
          typeTaken[tempString][myIndex]=true;
          // TODO: put this in a better place
          // (and make a better module typing)
          struct Visitor : public ConstGeometryVisitor {
            std::map<string, SummaryTable>& result;

            Visitor(std::map<std::string, SummaryTable>& result_) : result(result_) {}
            void visit(const BarrelModule& m) {
              string s = m.cntName() + " (L" + any2str(m.layer()) + ")";
              result[s].setCell(0, m.ring(), TagMaker::makePosTag(m));
            }
            void visit(const EndcapModule& m) {
              string s = m.cntName() + " (D" + any2str(m.disk()) + ")";
              result[s].setCell(0, m.ring(), TagMaker::makePosTag(m));
            }
          };
          Visitor v(result);
          myModule->accept(v);
        }
        if (byMaterial) { // sort by Material tag
          //nLocalMasses = myModuleCap->localMassCount();
          localmassesBegin = myModuleCap->getLocalMasses().begin();
          localmassesEnd = myModuleCap->getLocalMasses().end();
        } else { // sort by Component tag
          //nLocalMasses = myModuleCap->localMassCompCount();
          localmassesBegin = myModuleCap->getLocalMassesComp().begin();
          localmassesEnd = myModuleCap->getLocalMassesComp().end();
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
      TagMaker tmak(*myModule);

      if (byMaterial) {
        typeWeight[tmak.posTag]+=myModuleCap->getLocalMass();
        tagWeight[tmak.sensorGeoTag]+=myModuleCap->getLocalMass();
      }
      if (myModule->posRef().phi == 1) {
        // If we did not write this module type yet
        pair<int, int> myIndex = make_pair(myModule->tableRef().row/*+myModule->getDisk()*/, myModule->tableRef().col);
        tempString = myModule->cntName();
        if (!typeWritten[tempString][myIndex]) {
          typeWritten[tempString][myIndex]=true;
          if (tempString!="") {
            // TODO: put this in a better place
            // (and make a better module typing)
            tempSS.str("");
            if (myModule->subdet() == BARREL) {
              tempSS << ((BarrelModule*)myModule)->layer();
              tempString+=" (L"+tempSS.str()+")";
            } else if (myModule->subdet() == ENDCAP) {
              tempSS << ((EndcapModule*)myModule)->disk(); //getDisk();
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
            // std::cout << "Material\tLocal\n";

            // Then we fill in the proper column
            // Prepare the columns of the table
            for (unsigned int materialTag_i=0; materialTag_i<materialTagV.size(); ++materialTag_i) {
              materialTag = materialTagV[materialTag_i];
              if (byMaterial) { // table by materials
                try { localMaterial = myModuleCap->getLocalMass(materialTag); }
                catch (exception e) { localMaterial = 0; }
              } else { // table by components
                try { localMaterial = myModuleCap->getLocalMassComp(materialTag); }
                catch (exception e) { localMaterial = 0; }
              }
              //cout << materialTag << "\t"
              //<< localMaterial << "\t"
              //<< endl;
              tempSS.str("");
              // TODO: move this to Vizard
              tempSS << std::dec << std::fixed << std::setprecision(1) << localMaterial;
              result[tempString].setCell(materialTag_i+1, myModule->tableRef().col, tempSS.str());
            }
            localMaterial = myModuleCap->getLocalMass();
            tempSS.str("");
            tempSS << std::dec << std::fixed << std::setprecision(1) << localMaterial;
            result[tempString].setCell(materialTagV.size()+1, myModule->tableRef().col, tempSS.str());
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

  MaterialBillAnalyzer v;
  v.inspectTracker(mb);
  billOfMaterials_ = v.outputTable;
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

void printPosRefString(std::ostream& os, const Module& m, const string& delim = " ") {
  os << "cnt=" << m.posRef().cnt << delim << "z=" << m.posRef().z << delim << "rho=" << m.posRef().rho << " (" << m.center().Rho() << ")" << delim << "phi=" << m.posRef().phi << delim << "side=" << m.side();
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
    if (iter->getModule().maxZ() > 0) {
        // same method as in Tracker, same function used
        // TODO: in case origin==0,0,0 and phi==0 just check if sectionYZ and minEta, maxEta
        //distance = iter->getModule().trackCross(origin, direction);
        auto h = iter->getModule().checkTrackHits(origin, direction);
        if (h.second != HitType::NONE) {
          distance = h.first.R();
          HitType type = h.second;
          // module was hit
          hits++;
          r = distance * sin(theta);
          tmp.radiation = iter->getRadiationLength();
          tmp.interaction = iter->getInteractionLength();

          Module& m = iter->getModule();
          double tiltAngle = m.tiltAngle();
          // 2D material maps
          fillMapRT(r, theta, tmp);
          // radiation and interaction length scaling for barrels
          if (iter->getModule().subdet() == BARREL) {
            tmp.radiation = tmp.radiation / sin(theta + tiltAngle);
            tmp.interaction = tmp.interaction / sin(theta + tiltAngle);
          }
          // radiation and interaction length scaling for endcaps
          else {
            tmp.radiation = tmp.radiation / cos(theta + tiltAngle - M_PI/2);
            tmp.interaction = tmp.interaction / cos(theta + tiltAngle - M_PI/2);
          }

          double tmpr = 0., tmpi = 0.;

          std::map<std::string, Material> moduleComponentsRI = iter->getComponentsRI();
          for (std::map<std::string, Material>::iterator cit = moduleComponentsRI.begin(); cit != moduleComponentsRI.end(); ++cit) {
            sumComponentsRI[cit->first].radiation += cit->second.radiation / (iter->getModule().subdet() == BARREL ? sin(theta + tiltAngle) : cos(theta + tiltAngle - M_PI/2));
            //if (cit->first == "SupportMechanics") std::cout << eta << " " << distance << " " << cit->second.radiation / sin(theta + tiltAngle) << " " << cit->second.radiation << std::endl;
            tmpr += sumComponentsRI[cit->first].radiation;
            sumComponentsRI[cit->first].interaction += cit->second.interaction / (iter->getModule().subdet() == BARREL ? sin(theta + tiltAngle) : cos(theta + tiltAngle - M_PI/2));
            tmpi += sumComponentsRI[cit->first].interaction;
          }
          // 2D plot and eta plot results
          if (!isPixel) fillCell(r, eta, theta, tmp);
          res += tmp;
          // create Hit object with appropriate parameters, add to Track t
          Hit* hit = new Hit(distance, &(iter->getModule()), type);
          //if (iter->getModule().getSubdetectorType() == Module::Barrel) hit->setOrientation(Hit::Horizontal); // should not be necessary
          //else if(iter->getModule().getSubdetectorType() == Module::Endcap) hit->setOrientation(Hit::Vertical); // should not be necessary
          //hit->setObjectKind(Hit::Active); // should not be necessary
          hit->setCorrectedMaterial(tmp);
          hit->setPixel(isPixel);
          t.addHit(hit);
        }
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
                                   TrackNew& t, bool isPixel) {
  std::vector<std::vector<ModuleCap> >::iterator iter = tr.begin();
  std::vector<std::vector<ModuleCap> >::iterator guard = tr.end();
  Material res, tmp;
  res.radiation= 0.0;
  res.interaction = 0.0;
  while (iter != guard) {
    tmp = findHitsModuleLayer(*iter, t, isPixel);
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

  for (auto aModule : tracker.modules()) {

    // collision detection: rays are in z+ only, so consider only modules that lie on that side
    if (aModule->maxZ() > 0) {

      // same method as in Tracker, same function used
      //distance = aModule->trackCross(origin, direction);
      auto ht = aModule->checkTrackHits(origin, direction);
      //if (distance > 0) {
      if (ht.second != HitType::NONE) {
        double distance = ht.first.R();
        // module was hit
        hits++;

        // create Hit object with appropriate parameters, add to Track t
        Hit* hit = new Hit(distance, aModule, ht.second);
        hit->setCorrectedMaterial(emptyMaterial);
        t.addHit(hit);
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
Material Analyzer::findHitsModuleLayer(std::vector<ModuleCap>& layer, TrackNew& t, bool isPixel) {
  std::vector<ModuleCap>::iterator iter = layer.begin();
  std::vector<ModuleCap>::iterator guard = layer.end();
  Material res, tmp;
  XYZVector origin, direction;
  origin    = t.getOrigin();
  direction = t.getDirection();
  double distance;
  //double r;
  int hits = 0;
  res.radiation = 0.0;
  res.interaction = 0.0;
  // set the track direction vector
  while (iter != guard) {
    // collision detection: rays are in z+ only, so consider only modules that lie on that side
    if (iter->getModule().maxZ() > 0) {
        // same method as in Tracker, same function used
        // TODO: in case origin==0,0,0 and phi==0 just check if sectionYZ and minEta, maxEta
        //distance = iter->getModule().trackCross(origin, direction);
        auto h = iter->getModule().checkTrackHits(origin, direction); 
        if (h.second != HitType::NONE) {
        //if (distance > 0) {
          double distance = h.first.R();
          // module was hit
          hits++;
          // r = distance * sin(theta);
          tmp.radiation = iter->getRadiationLength();
          tmp.interaction = iter->getInteractionLength();
          // radiation and interaction length scaling for barrels
          if (iter->getModule().subdet() == BARREL) {
            tmp.radiation = tmp.radiation / sin(t.getTheta());
            tmp.interaction = tmp.interaction / sin(t.getTheta());
          }
          // radiation and interaction length scaling for endcaps
          else {
            tmp.radiation = tmp.radiation / cos(t.getTheta());
            tmp.interaction = tmp.interaction / cos(t.getTheta());
          }
          res += tmp;
          // Create Hit object with appropriate parameters, add to Track t
          auto hitRPos = h.first.rho();
          auto hitZPos = h.first.z();
          auto hitType = h.second;

          HitNewPtr hit(new HitNew(hitRPos, hitZPos, &(iter->getModule()), hitType));
          hit->setCorrectedMaterial(tmp);
          if (isPixel) hit->setAsPixel();
          t.addHit(std::move(hit));
        }
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
 * @param A boolean flag to indicate which set of active surfaces isa nalysed: true if the belong to a pixel detector, false if they belong to the tracker
 * @return The scaled and summed up radiation and interaction lengths for the given collection of elements and track, bundled into a <i>std::pair</i>
 */

Material Analyzer::analyzeInactiveSurfaces(std::vector<InactiveElement>& elements, double eta,
                                           double theta, Track& t, std::map<std::string, Material>& sumServicesComponentsRI, MaterialProperties::Category cat, bool isPixel) {

  /*
  for (InactiveElement& currElem : elements) {
    currElem.calculateTotalMass();
    currElem.calculateRadiationLength();
    currElem.calculateInteractionLength();
  }
  */
  
  std::vector<InactiveElement>::iterator iter = elements.begin();
  std::vector<InactiveElement>::iterator guard = elements.end();
  Material res, corr;
  std::pair<double, double> tmp;
  double s = 0.0;
  while ((iter != guard)) {
    //if  ((iter->getInteractionLength() > 0) && (iter->getRadiationLength() > 0)) {
    // collision detection: rays are in z+ only, so only volumes in z+ need to be considered
    // only volumes of the requested category, or those without one (which should not exist) are examined
    if (((iter->getZOffset() + iter->getZLength()) > 0)
        && ((cat == MaterialProperties::no_cat) || (cat == iter->getCategory()))) {
      // collision detection: check eta range
      tmp = iter->getEtaMinMax();
      // volume was hit
      if ((tmp.first < eta) && (tmp.second > eta)) {
        double r, z;
        /*
        if (eta<0.01) {
          std::cout << "Hitting an inactive surface at z=("
                    << iter->getZOffset() << " to " << iter->getZOffset()+iter->getZLength()
                    << ") r=(" << iter->getInnerRadius() << " to " << iter->getInnerRadius()+iter->getRWidth() << ")" << std::endl;
          const std::map<std::string, double>& localMasses = iter->getLocalMasses();
          for (auto massIt : localMasses) std::cerr   << "       localMass" <<  massIt.first << " = " << any2str(massIt.second) << " g" << std::endl;
        }
        */
        // radiation and interaction lenth scaling for vertical volumes
        if (iter->isVertical()) {
          z = iter->getZOffset() + iter->getZLength() / 2.0;
          r = z * tan(theta);
          // 2D maps for vertical surfaces
          fillMapRZ(r,z,iter->getMaterialLengths());        
         
	  corr.radiation = iter->getRadiationLength() / cos(theta);
	  corr.interaction = iter->getInteractionLength() / cos(theta);
	  res += corr;
	  if (!isPixel) {
	    Material thisLength;
	    thisLength.radiation = corr.radiation;
	    thisLength.interaction = corr.interaction;
	    fillCell(r, eta, theta, thisLength);
	  }        
        }
        // radiation and interaction length scaling for horizontal volumes
        else {
          r = iter->getInnerRadius() + iter->getRWidth() / 2.0;
          // 2D maps for horizontal surfaces
          fillMapRT(r,theta,iter->getMaterialLengths());
          // special treatment for user-defined supports; should not be necessary for now
          // as all user-defined supports are vertical, but just in case...

	  corr.radiation = iter->getRadiationLength() / sin(theta);
	  corr.interaction = iter->getInteractionLength() / sin(theta);
	  res += corr;
	  if (!isPixel) {
	    Material thisLength;
	    thisLength.radiation = corr.radiation;
	    thisLength.interaction =  corr.interaction;
	    fillCell(r, eta, theta, thisLength); 
	  }             
        }

        // create Hit object with appropriate parameters, add to Track t
        Hit* hit = new Hit((theta == 0) ? r : (r / sin(theta)));
	hit->setHitInactiveElement(&(*iter));
        if (iter->isVertical()) hit->setOrientation(Hit::Vertical);
        else hit->setOrientation(Hit::Horizontal);
        hit->setObjectKind(Hit::Inactive);
        hit->setCorrectedMaterial(corr);
        hit->setPixel(isPixel);


	if ((iter->getCategory() != MaterialProperties::b_sup)
	    && (iter->getCategory() != MaterialProperties::e_sup)
	    && (iter->getCategory() != MaterialProperties::o_sup)
	    && (iter->getCategory() != MaterialProperties::u_sup)
	    && (iter->getCategory() != MaterialProperties::t_sup)) {
	  /*std::map<std::string, Material> servicesComponentsRI = iter->getComponentsRI();
          for (const auto& it : servicesComponentsRI) {
            sumServicesComponentsRI[it.first].radiation += it.second.radiation / (iter->isVertical() ? cos(theta) : sin(theta));  
            sumServicesComponentsRI[it.first].interaction += it.second.interaction / (iter->isVertical() ? cos(theta) : sin(theta));
	    }*/
	  sumServicesComponentsRI["Services : others"].radiation += corr.radiation;
	  sumServicesComponentsRI["Services : others"].interaction += corr.interaction;

	}


	if ((iter->getCategory() == MaterialProperties::b_ser)
	    || (iter->getCategory() == MaterialProperties::e_ser)) {
	  hit->setObjectCategory(Hit::Service);
	}
	else if ((iter->getCategory() == MaterialProperties::b_sup)
		 || (iter->getCategory() == MaterialProperties::e_sup)
		 || (iter->getCategory() == MaterialProperties::o_sup)
		 || (iter->getCategory() == MaterialProperties::t_sup)) {
	  hit->setObjectCategory(Hit::Support);
	}
	else if (iter->getCategory() == MaterialProperties::no_cat) {
	  hit->setObjectCategory(Hit::Service);
	}
	else {
	  hit->setObjectCategory(Hit::Unknown);
	}


	t.addHit(hit);
      }
    }
    iter++;
    //}
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
Material Analyzer::findHitsInactiveSurfaces(std::vector<InactiveElement>& elements, TrackNew& t, bool isPixel) {
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
      if ((tmp.first < t.getEta()) && (tmp.second > t.getEta())) {
        double r, z;
        // radiation and interaction lenth scaling for vertical volumes
        if (iter->isVertical()) { // Element is vertical
          z = iter->getZOffset() + iter->getZLength() / 2.0;
          r = z * tan(t.getTheta());

          // In case we are crossing the material with a very shallow angle
          // we have to take into account its finite radial size
          s_normal = iter->getZLength() / cos(t.getTheta());
          s_alternate = iter->getRWidth() / sin(t.getTheta());
          if (s_normal > s_alternate) { 
            // Special case: it's easier to cross the material by going left-to-right than
            // by going bottom-to-top, so I have to rescale the material amount computation
            corr.radiation = iter->getRadiationLength() / iter->getZLength() * s_alternate;
            corr.interaction = iter->getInteractionLength() / iter->getZLength() * s_alternate;
            res += corr;
          } else {
            // Standard computing of the crossed material amount
            corr.radiation = iter->getRadiationLength() / cos(t.getTheta());
            corr.interaction = iter->getInteractionLength() / cos(t.getTheta());
            res += corr;
          }
        }
        // radiation and interaction length scaling for horizontal volumes
        else { // Element is horizontal
          r = iter->getInnerRadius() + iter->getRWidth() / 2.0;

          // In case we are crossing the material with a very shallow angle
          // we have to take into account its finite z length
          s_normal = iter->getRWidth() / sin(t.getTheta());
          s_alternate = iter->getZLength() / cos(t.getTheta());
          if (s_normal > s_alternate) { 
            // Special case: it's easier to cross the material by going left-to-right than
            // by going bottom-to-top, so I have to rescale the material amount computation
            corr.radiation = iter->getRadiationLength() / iter->getRWidth() * s_alternate;
            corr.interaction = iter->getInteractionLength() / iter->getRWidth() * s_alternate;
            res += corr;
          } else {
            // Standard computing of the crossed material amount
            corr.radiation = iter->getRadiationLength() / sin(t.getTheta());
            corr.interaction = iter->getInteractionLength() / sin(t.getTheta());
            res += corr;
          }
        }
        // Create Hit object with appropriate parameters, add to Track t
        double rPos = r;
        double zPos = r/tan(t.getTheta());

        HitNewPtr hit(new HitNew(rPos, zPos));
        hit->setAsPassive();
        hit->setCorrectedMaterial(corr);
        t.addHit(std::move(hit));

//        Hit* hit = new Hit((theta == 0) ? r : (r / sin(theta)));
//        if (iter->isVertical()) hit->setOrientation(Hit::Vertical);
//        else hit->setOrientation(Hit::Horizontal);
//        hit->setObjectKind(Hit::Inactive);
//        hit->setCorrectedMaterial(corr);
//        hit->setPixel(isPixel);
//        t.addHit(hit);
	// std::cout << "OLD USED" << std::endl;
      }
    }
    iter++;
  }
  return res;
}
  
void Analyzer::clearGraphsPt(int graphAttributes, const std::string& graphTag) {
  std::map<int, TGraph>& thisRhoGraphs_Pt      = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::RhoGraph_Pt      ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::RhoGraph_Pt     , graphTag);
  std::map<int, TGraph>& thisPhiGraphs_Pt      = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::PhiGraph_Pt      ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::PhiGraph_Pt     , graphTag);
  std::map<int, TGraph>& thisDGraphs_Pt        = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::DGraph_Pt        ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::DGraph_Pt       , graphTag);
  std::map<int, TGraph>& thisCtgThetaGraphs_Pt = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::CtgthetaGraph_Pt ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::CtgthetaGraph_Pt, graphTag);
  std::map<int, TGraph>& thisZ0Graphs_Pt       = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::Z0Graph_Pt       ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::Z0Graph_Pt      , graphTag);
  std::map<int, TGraph>& thisPGraphs_Pt        = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::PGraph_Pt        ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::PGraph_Pt       , graphTag);
  std::map<int, TGraph>& thisLGraphs_Pt        = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::LGraph_Pt        ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::LGraph_Pt       , graphTag);
  std::map<int, TGraph>& thisBetaGraphs_Pt     = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::BetaGraph_Pt     ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::BetaGraph_Pt    , graphTag);
  std::map<int, TGraph>& thisOmegaGraphs_Pt    = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::OmegaGraph_Pt    ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::OmegaGraph_Pt   , graphTag);

  thisRhoGraphs_Pt.clear();
  thisPhiGraphs_Pt.clear();
  thisDGraphs_Pt.clear();
  thisCtgThetaGraphs_Pt.clear();
  thisZ0Graphs_Pt.clear();
  thisPGraphs_Pt.clear();
  thisLGraphs_Pt.clear();
  thisBetaGraphs_Pt.clear();
  thisOmegaGraphs_Pt.clear();

}

void Analyzer::clearGraphsP(int graphAttributes, const std::string& graphTag) {
  std::map<int, TGraph>& thisRhoGraphs_P      = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::RhoGraph_P      ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::RhoGraph_P     , graphTag);
  std::map<int, TGraph>& thisPhiGraphs_P      = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::PhiGraph_P      ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::PhiGraph_P     , graphTag);
  std::map<int, TGraph>& thisDGraphs_P        = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::DGraph_P        ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::DGraph_P       , graphTag);
  std::map<int, TGraph>& thisCtgThetaGraphs_P = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::CtgthetaGraph_P ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::CtgthetaGraph_P, graphTag);
  std::map<int, TGraph>& thisZ0Graphs_P       = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::Z0Graph_P       ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::Z0Graph_P      , graphTag);
  std::map<int, TGraph>& thisPGraphs_P        = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::PGraph_P        ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::PGraph_P       , graphTag);
  std::map<int, TGraph>& thisLGraphs_P        = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::LGraph_P        ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::LGraph_P       , graphTag);
  std::map<int, TGraph>& thisBetaGraphs_P     = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::BetaGraph_P     ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::BetaGraph_P    , graphTag);
  std::map<int, TGraph>& thisOmegaGraphs_P    = graphTag.empty() ? myGraphBag.getGraphs(graphAttributes | GraphBag::OmegaGraph_P    ) : myGraphBag.getTaggedGraphs(graphAttributes | GraphBag::OmegaGraph_P   , graphTag);

  thisRhoGraphs_P.clear();
  thisPhiGraphs_P.clear();
  thisDGraphs_P.clear();
  thisCtgThetaGraphs_P.clear();
  thisZ0Graphs_P.clear();
  thisPGraphs_P.clear();
  thisLGraphs_P.clear();
  thisBetaGraphs_P.clear();
  thisOmegaGraphs_P.clear();

}
 
/**
 * Calculate the error graphs for the radius curvature, the distance and the angle, for each momentum,
 * and store them internally for later visualisation.
 * @param parameter The list of different momenta that the error graphs are calculated for
 */
void Analyzer::calculateGraphsConstPt(const int& parameter,
                                      const TrackNewCollection& aTrackCollection,
                                      int graphAttributes,
                                      const string& graphTag) {

  // Get graphs from graphBag
  TGraph& thisRhoGraph_Pt       = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::RhoGraph_Pt     , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::RhoGraph_Pt       , graphTag, parameter);
  TGraph& thisPhiGraph_Pt       = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::PhiGraph_Pt     , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::PhiGraph_Pt       , graphTag, parameter);
  TGraph& thisDGraph_Pt         = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::DGraph_Pt       , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::DGraph_Pt         , graphTag, parameter);
  TGraph& thisCtgThetaGraph_Pt  = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::CtgthetaGraph_Pt, parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::CtgthetaGraph_Pt  , graphTag, parameter);
  TGraph& thisZ0Graph_Pt        = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::Z0Graph_Pt      , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::Z0Graph_Pt        , graphTag, parameter);
  TGraph& thisPGraph_Pt         = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::PGraph_Pt       , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::PGraph_Pt         , graphTag, parameter);
  TGraph& thisLGraph_Pt         = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::LGraph_Pt       , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::LGraph_Pt         , graphTag, parameter);
  TGraph& thisBetaGraph_Pt      = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::BetaGraph_Pt    , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::BetaGraph_Pt      , graphTag, parameter);
  TGraph& thisOmegaGraph_Pt     = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::OmegaGraph_Pt   , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::OmegaGraph_Pt     , graphTag, parameter);

  // Variables
  double eta, R;
  std::ostringstream aName;

  // Prepare plot for const pt across eta
  double momentum = double(parameter)*Units::MeV/Units::GeV;

  // Prepare plots: pt
  thisRhoGraph_Pt.SetTitle("p_{T} resolution versus #eta - const P_{T} across #eta;#eta;#delta p_{T}/p_{T} [%]");
  aName.str(""); aName << "pt_vs_eta" << momentum << graphTag;
  thisRhoGraph_Pt.SetName(aName.str().c_str());
  // Prepare plots: phi
  thisPhiGraph_Pt.SetTitle("Track azimuthal angle error - const P_{T} across #eta;#eta;#delta #phi [deg]");
  aName.str(""); aName << "phi_vs_eta" << momentum << graphTag;
  thisPhiGraph_Pt.SetName(aName.str().c_str());
  // Prepare plots: d
  thisDGraph_Pt.SetTitle("Transverse impact parameter error - const P_{T} across #eta;#eta;#delta d_{0} [#mum]");
  aName.str(""); aName << "d_vs_eta" << momentum << graphTag;
  thisDGraph_Pt.SetName(aName.str().c_str());
  // Prepare plots: ctg(theta)
  thisCtgThetaGraph_Pt.SetTitle("Track polar angle error - const P_{T} across #eta;#eta;#delta ctg(#theta)");
  aName.str(""); aName << "ctgTheta_vs_eta" << momentum << graphTag;
  thisCtgThetaGraph_Pt.SetName(aName.str().c_str());
  // Prepare plots: z0
  thisZ0Graph_Pt.SetTitle("Longitudinal impact parameter error - const P_{T} across #eta;#eta;#delta z_{0} [#mum]");
  aName.str(""); aName << "z_vs_eta" << momentum << graphTag;
  thisZ0Graph_Pt.SetName(aName.str().c_str());
  // Prepare plots: p
  thisPGraph_Pt.SetTitle("p resolution versus #eta - const P_{T} across #eta;#eta;#delta p/p [%]");
  aName.str(""); aName << "p_vs_eta" << momentum << graphTag;
  thisPGraph_Pt.SetName(aName.str().c_str());
  // Prepare plots: L
  thisLGraph_Pt.SetTitle("c#tau resolution versus #eta - const P_{T} across #eta;#eta;#delta c#tau [#mum]");
  aName.str(""); aName << "l_vs_eta" << momentum << graphTag;
  thisLGraph_Pt.SetName(aName.str().c_str());
  // Prepare plots: Beta
  thisBetaGraph_Pt.SetTitle("#beta - const P_{T} across #eta;#eta;#beta [rad]");
  aName.str(""); aName << "beta_vs_eta" << momentum << graphTag;
  thisBetaGraph_Pt.SetName(aName.str().c_str());
  // Prepare plots: Omega
  thisOmegaGraph_Pt.SetTitle("#omega - const P_{T} across #eta;#eta;#omega [rad]");
  aName.str(""); aName << "omega_vs_eta" << momentum << graphTag;
  thisOmegaGraph_Pt.SetName(aName.str().c_str());

  // track loop
  double graphValue;
  double rPos = 0.0; // At [0,0] point
  for ( const auto& myTrack : aTrackCollection ) {
    const double& dpt  = myTrack->getDeltaPtOverPt(rPos);
    const double& dphi0= myTrack->getDeltaPhi0();
    const double& dd0  = myTrack->getDeltaD0();
    const double& dctg = myTrack->getDeltaCtgTheta(rPos);
    const double& dz0  = myTrack->getDeltaZ0();
    const double& dp   = myTrack->getDeltaPOverP(rPos);

    /*std::vector<std::pair<Module*,HitType>> hitModules = myTrack.getHitModules();
    if ( hitModules.at(0).first->getObjectKind() == Active) {
    std::cout << "hitModules.at(0).first->getResolutionLocalX() = " << hitModules.at(0).first->getResolutionLocalX() << std::endl;
    }*/

//    std::vector<Hit*> hitModules = myTrack.getHitV();
//    //std::cout << "hitModules.at(0)->getObjectKind() = " << hitModules.at(0)->getObjectKind() << std::endl;
//    //std::cout << "Hit::Inactive = " << Hit::Inactive << std::endl;
//    for (auto& mh : hitModules) {
//    if ( mh->getObjectKind() == Hit::Active) {
//      if (mh->getHitModule()) {
//	//std::cout << "mh->getResolutionLocalX() = " << mh->getResolutionLocalX() << std::endl;
//      }
//    }
//    }


    eta = myTrack->getEta();
    double theta = myTrack->getTheta();
    if (dpt>0) {
      // deltaRho / rho = deltaRho * R
      graphValue = dpt*100; // in percent
      thisRhoGraph_Pt.SetPoint(thisRhoGraph_Pt.GetN(), eta, graphValue);
    }
    if (dphi0>0) {
      graphValue = dphi0 / Units::deg; // in degrees
      thisPhiGraph_Pt.SetPoint(thisPhiGraph_Pt.GetN(), eta, graphValue);
    }
    if (dd0>0) {
      graphValue = dd0 / Units::um ; // in um
      thisDGraph_Pt.SetPoint(thisDGraph_Pt.GetN(), eta, graphValue );
    }
    if (dctg>0) {
      graphValue = dctg; // An absolute number
      thisCtgThetaGraph_Pt.SetPoint(thisCtgThetaGraph_Pt.GetN(), eta, graphValue);
    }
    if (dz0>0) {
      graphValue =  (dz0) / Units::um ; // in um
      thisZ0Graph_Pt.SetPoint(thisZ0Graph_Pt.GetN(), eta, graphValue);
    }
    if ((dp>0)||true) {
      graphValue = dp * 100.; // in percent
      thisPGraph_Pt.SetPoint(thisPGraph_Pt.GetN(), eta, graphValue);
    }
    if ((dd0>0)&&(dz0>0)) {
      double resolution = (cos(theta)*cos(theta)+1)/dd0/dd0 + sin(theta)*sin(theta)/dz0/dz0;
      resolution = sqrt(2/resolution) / Units::um ; // resolution in um
      thisLGraph_Pt.SetPoint(thisLGraph_Pt.GetN(), eta, resolution);

      double beta = (cos(theta)*cos(theta)+1)*dz0*dz0*dz0/(sin(theta)*sin(theta)*dd0*dd0*dd0);
      beta = atan(beta);
      thisBetaGraph_Pt.SetPoint(thisBetaGraph_Pt.GetN(), eta, beta);

      double omega = ((cos(theta)*cos(theta)+1)-sin(theta)*sin(theta)*dz0*dz0/dd0/dd0)/dz0/dd0;
      omega = atan(omega);
      thisOmegaGraph_Pt.SetPoint(thisOmegaGraph_Pt.GetN(), eta, omega);
    }
  }
}

/**
 * Calculate the error graphs for case II with const P across eta: the radius curvature, the distance and the angle, for each momentum,
 * and store them internally for later visualisation.
 * @param parameter The list of different momenta that the error graphs are calculated for
 */
void Analyzer::calculateGraphsConstP(const int& parameter,
                                     const TrackNewCollection& aTrackCollection,
                                     int graphAttributes,
                                     const string& graphTag) {

  // Get graphs from graphBag
  TGraph& thisRhoGraph_P       = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::RhoGraph_P     , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::RhoGraph_P       , graphTag, parameter);
  TGraph& thisPhiGraph_P       = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::PhiGraph_P     , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::PhiGraph_P       , graphTag, parameter);
  TGraph& thisDGraph_P         = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::DGraph_P       , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::DGraph_P         , graphTag, parameter);
  TGraph& thisCtgThetaGraph_P  = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::CtgthetaGraph_P, parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::CtgthetaGraph_P  , graphTag, parameter);
  TGraph& thisZ0Graph_P        = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::Z0Graph_P      , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::Z0Graph_P        , graphTag, parameter);
  TGraph& thisPGraph_P         = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::PGraph_P       , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::PGraph_P         , graphTag, parameter);
  TGraph& thisLGraph_P         = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::LGraph_P       , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::LGraph_P         , graphTag, parameter);
  TGraph& thisBetaGraph_P      = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::BetaGraph_P    , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::BetaGraph_P      , graphTag, parameter);
  TGraph& thisOmegaGraph_P     = graphTag.empty() ? myGraphBag.getGraph(graphAttributes | GraphBag::OmegaGraph_P   , parameter ) : myGraphBag.getTaggedGraph(graphAttributes | GraphBag::OmegaGraph_P     , graphTag, parameter);

  // Variables
  double eta, R;
  std::ostringstream aName;

  // Prepare plot for const pt across eta
  double momentum = double(parameter)*Units::MeV/Units::GeV;

  // Prepare plots: pt
  thisRhoGraph_P.SetTitle("p_{T} resolution versus #eta - const P across #eta;#eta;#delta p_{T}/p_{T} [%]");
  aName.str(""); aName << "pt_vs_eta" << momentum << graphTag;
  thisRhoGraph_P.SetName(aName.str().c_str());
  // Prepare plots: phi
  thisPhiGraph_P.SetTitle("Track azimuthal angle error - const P across #eta;#eta;#delta #phi [deg]");
  aName.str(""); aName << "phi_vs_eta" << momentum << graphTag;
  thisPhiGraph_P.SetName(aName.str().c_str());
  // Prepare plots: d
  thisDGraph_P.SetTitle("Transverse impact parameter error - const P across #eta;#eta;#delta d_{0} [#mum]");
  aName.str(""); aName << "d_vs_eta" << momentum << graphTag;
  thisDGraph_P.SetName(aName.str().c_str());
  // Prepare plots: ctg(theta)
  thisCtgThetaGraph_P.SetTitle("Track polar angle error - const P across #eta;#eta;#delta ctg(#theta)");
  aName.str(""); aName << "ctgTheta_vs_eta" << momentum << graphTag;
  thisCtgThetaGraph_P.SetName(aName.str().c_str());
  // Prepare plots: z0
  thisZ0Graph_P.SetTitle("Longitudinal impact parameter error - const P across #eta;#eta;#delta z_{0} [#mum]");
  aName.str(""); aName << "z_vs_eta" << momentum << graphTag;
  thisZ0Graph_P.SetName(aName.str().c_str());
  // Prepare plots: p
  thisPGraph_P.SetTitle("p resolution versus #eta - const P across #eta;#eta;#delta p/p [%]");
  aName.str(""); aName << "p_vs_eta" << momentum << graphTag;
  thisPGraph_P.SetName(aName.str().c_str());
  // Prepare plots: L
  thisLGraph_P.SetTitle("c#tau resolution versus #eta - const P across #eta;#eta;#delta c#tau [#mum]");
  aName.str(""); aName << "l_vs_eta" << momentum << graphTag;
  thisLGraph_P.SetName(aName.str().c_str());
  // Prepare plots: Beta
  thisBetaGraph_P.SetTitle("#beta - const P across #eta;#eta;#beta [rad]");
  aName.str(""); aName << "beta_vs_eta" << momentum << graphTag;
  thisBetaGraph_P.SetName(aName.str().c_str());
  // Prepare plots: Omega
  thisOmegaGraph_P.SetTitle("#omega - const P across #eta;#eta;#omega [rad]");
  aName.str(""); aName << "omega_vs_eta" << momentum << graphTag;
  thisOmegaGraph_P.SetName(aName.str().c_str());

  // track loop
  double graphValue;
  double rPos; // At [0,0] point
  for ( const auto& myTrack : aTrackCollection ) {
    const double& dpt  = myTrack->getDeltaPtOverPt(rPos);
    const double& dphi0= myTrack->getDeltaPhi0();
    const double& dd0  = myTrack->getDeltaD0();
    const double& dctg = myTrack->getDeltaCtgTheta(rPos);
    const double& dz0  = myTrack->getDeltaZ0();
    const double& dp   = myTrack->getDeltaPOverP(rPos);

    eta = myTrack->getEta();
    double theta = myTrack->getTheta();
    if (dpt>0) {
      // deltaRho / rho = deltaRho * R
      graphValue = dpt* 100; // in percent
      thisRhoGraph_P.SetPoint(thisRhoGraph_P.GetN(), eta, graphValue);
    }
    if (dphi0>0) {
      graphValue = dphi0/Units::deg; // in degrees
      thisPhiGraph_P.SetPoint(thisPhiGraph_P.GetN(), eta, graphValue);
    }
    if (dd0>0) {
      graphValue = dd0 / Units::um;
      thisDGraph_P.SetPoint(thisDGraph_P.GetN(), eta, graphValue );
    }
    if (dctg>0) {
      graphValue = dctg; // An absolute number
      thisCtgThetaGraph_P.SetPoint(thisCtgThetaGraph_P.GetN(), eta, graphValue);
    }
    if (dz0>0) {
      graphValue =  (dz0) / Units::um;
      thisZ0Graph_P.SetPoint(thisZ0Graph_P.GetN(), eta, graphValue);
    }
    if ((dp>0)||true) {
      graphValue = dp * 100.; // in percent
      thisPGraph_P.SetPoint(thisPGraph_P.GetN(), eta, graphValue);
    }
    if ((dd0>0)&&(dz0>0)) {
      double resolution = (cos(theta)*cos(theta)+1)/dd0/dd0 + sin(theta)*sin(theta)/dz0/dz0;
      resolution = sqrt(2/resolution) / Units::um ; // resolution in um
      thisLGraph_P.SetPoint(thisLGraph_P.GetN(), eta, resolution);

      double beta = (cos(theta)*cos(theta)+1)*dz0*dz0*dz0/(sin(theta)*sin(theta)*dd0*dd0*dd0);
      beta = atan(beta);
      thisBetaGraph_P.SetPoint(thisBetaGraph_P.GetN(), eta, beta);

      double omega = ((cos(theta)*cos(theta)+1)-sin(theta)*sin(theta)*dz0*dz0/dd0/dd0)/dz0/dd0;
      omega = atan(omega);
      thisOmegaGraph_P.SetPoint(thisOmegaGraph_P.GetN(), eta, omega);
    }

  }
}

  /**
   * Creates the modules' parametrized spatial resolution profiles and distributions
   * @param taggedTrackPtCollectionMap Tagged collections of tracks
   */
  void Analyzer::calculateParametrizedResolutionPlots(std::map<std::string, TrackNewCollectionMap>& taggedTrackPtCollectionMap) {

    for (auto& ttcmIt : taggedTrackPtCollectionMap) {
      const string& myTag = ttcmIt.first;
      TrackNewCollectionMap& myTrackCollection = ttcmIt.second;

      // Modules' parametrized spatial resolution maps
      parametrizedResolutionLocalXBarrelMap[myTag].Reset();
      parametrizedResolutionLocalXBarrelMap[myTag].SetNameTitle("resoXBarMap","Resolution on local X coordinate vs cotg(alpha) (barrel modules)");
      parametrizedResolutionLocalXBarrelMap[myTag].SetBins(500, -0.3, 0.3, 500, 0, 30);
      parametrizedResolutionLocalXBarrelMap[myTag].GetXaxis()->SetTitle("cotg(alpha)");
      parametrizedResolutionLocalXBarrelMap[myTag].GetYaxis()->SetTitle("resolutionLocalX [um]");
      parametrizedResolutionLocalXBarrelMap[myTag].GetXaxis()->CenterTitle();
      parametrizedResolutionLocalXBarrelMap[myTag].GetYaxis()->CenterTitle();
      parametrizedResolutionLocalXBarrelMap[myTag].SetMarkerColor(kBlue + 2);

      parametrizedResolutionLocalXEndcapsMap[myTag].Reset();
      parametrizedResolutionLocalXEndcapsMap[myTag].SetNameTitle("resoXEndMap","Resolution on local X coordinate vs cotg(alpha) (endcaps modules)");
      parametrizedResolutionLocalXEndcapsMap[myTag].SetBins(500, -0.3, 0.3, 500, 0, 30);
      parametrizedResolutionLocalXEndcapsMap[myTag].GetXaxis()->SetTitle("cotg(alpha)");
      parametrizedResolutionLocalXEndcapsMap[myTag].GetYaxis()->SetTitle("resolutionLocalX [um]");
      parametrizedResolutionLocalXEndcapsMap[myTag].GetXaxis()->CenterTitle();
      parametrizedResolutionLocalXEndcapsMap[myTag].GetYaxis()->CenterTitle();
      parametrizedResolutionLocalXEndcapsMap[myTag].SetMarkerColor(kBlue + 2);

      parametrizedResolutionLocalYBarrelMap[myTag].Reset();
      parametrizedResolutionLocalYBarrelMap[myTag].SetNameTitle("resoYBarMap","Resolution on local Y coordinate vs |cotg(beta)| (barrel modules)");
      parametrizedResolutionLocalYBarrelMap[myTag].SetBins(500, 0, 5, 500, 0, 60);
      parametrizedResolutionLocalYBarrelMap[myTag].GetXaxis()->SetTitle("|cotg(beta)|");
      parametrizedResolutionLocalYBarrelMap[myTag].GetYaxis()->SetTitle("resolutionLocalY [um]");
      parametrizedResolutionLocalYBarrelMap[myTag].GetXaxis()->CenterTitle();
      parametrizedResolutionLocalYBarrelMap[myTag].GetYaxis()->CenterTitle();
      parametrizedResolutionLocalYBarrelMap[myTag].SetMarkerColor(kBlue + 2);

      parametrizedResolutionLocalYEndcapsMap[myTag].Reset();
      parametrizedResolutionLocalYEndcapsMap[myTag].SetNameTitle("resoYEndMap","Resolution on local Y coordinate vs |cotg(beta)| (endcaps modules)");
      parametrizedResolutionLocalYEndcapsMap[myTag].SetBins(500, 0, 0.8, 500, 0, 60);
      parametrizedResolutionLocalYEndcapsMap[myTag].GetXaxis()->SetTitle("|cotg(beta)|");
      parametrizedResolutionLocalYEndcapsMap[myTag].GetYaxis()->SetTitle("resolutionLocalY [um]");
      parametrizedResolutionLocalYEndcapsMap[myTag].GetXaxis()->CenterTitle();
      parametrizedResolutionLocalYEndcapsMap[myTag].GetYaxis()->CenterTitle();
      parametrizedResolutionLocalYEndcapsMap[myTag].SetMarkerColor(kBlue + 2);
 
      // Modules' parametrized spatial resolution distributions
      parametrizedResolutionLocalXBarrelDistribution[myTag].Reset();
      parametrizedResolutionLocalXBarrelDistribution[myTag].SetNameTitle("resoXBarDistr","Distribution of the resolution on local X coordinate (barrel modules)");
      parametrizedResolutionLocalXBarrelDistribution[myTag].SetBins(500, 0, 30);
      parametrizedResolutionLocalXBarrelDistribution[myTag].GetXaxis()->SetTitle("resolutionLocalX [um]");
      parametrizedResolutionLocalXBarrelDistribution[myTag].GetXaxis()->CenterTitle();

      parametrizedResolutionLocalXEndcapsDistribution[myTag].Reset();
      parametrizedResolutionLocalXEndcapsDistribution[myTag].SetNameTitle("resoXEndDistr","Distribution of the resolution on local X coordinate (endcaps modules)");
      parametrizedResolutionLocalXEndcapsDistribution[myTag].SetBins(500, 0, 30);
      parametrizedResolutionLocalXEndcapsDistribution[myTag].GetXaxis()->SetTitle("resolutionLocalX [um]");
      parametrizedResolutionLocalXEndcapsDistribution[myTag].GetXaxis()->CenterTitle();

      parametrizedResolutionLocalYBarrelDistribution[myTag].Reset();
      parametrizedResolutionLocalYBarrelDistribution[myTag].SetNameTitle("resoYBarDistr","Distribution of the resolution on local Y coordinate (barrel modules)");
      parametrizedResolutionLocalYBarrelDistribution[myTag].SetBins(500, 0, 60);
      parametrizedResolutionLocalYBarrelDistribution[myTag].GetXaxis()->SetTitle("resolutionLocalY [um]");
      parametrizedResolutionLocalYBarrelDistribution[myTag].GetXaxis()->CenterTitle();

      parametrizedResolutionLocalYEndcapsDistribution[myTag].Reset();
      parametrizedResolutionLocalYEndcapsDistribution[myTag].SetNameTitle("resoYEndDistr","Distribution of the resolution on local Y coordinate (endcaps modules)");
      parametrizedResolutionLocalYEndcapsDistribution[myTag].SetBins(500, 0, 60);
      parametrizedResolutionLocalYEndcapsDistribution[myTag].GetXaxis()->SetTitle("resolutionLocalY [um]");
      parametrizedResolutionLocalYEndcapsDistribution[myTag].GetXaxis()->CenterTitle();
 
      for (const auto& tcmIt : myTrackCollection) {
	//const int &parameter = tcmIt.first;
	const TrackNewCollection& myCollection = tcmIt.second;

 	// track loop
	for ( const auto& myTrack : myCollection ) {

	  // hit loop
	  for (std::vector<std::unique_ptr<HitNew>>::const_iterator iHit=myTrack->getBeginHits(); iHit!=myTrack->getEndHits(); iHit++) {

	    // In case the tag is "tracker", takes only the outer tracker
	    if (myTag != "tracker" || (myTag == "tracker" && !(*iHit)->isPixel())) {
	      // Consider hit modules	
	      if ((*iHit)->isActive() && (*iHit)->getHitModule()) {
		
		    const auto& hitModule = (*iHit)->getHitModule();
		// If any parameter for resolution on local X coordinate specified for hitModule, fill maps and distributions
		if (hitModule->hasAnyResolutionLocalXParam()) {
		  double cotAlpha = 1./tan(hitModule->alpha(myTrack->getPhi()));
		  double resolutionLocalX = hitModule->resolutionLocalX(myTrack->getPhi())/Units::um; // um
		  if ( hitModule->subdet() == BARREL ) {
		    parametrizedResolutionLocalXBarrelMap[myTag].Fill(cotAlpha, resolutionLocalX);
		    parametrizedResolutionLocalXBarrelDistribution[myTag].Fill(resolutionLocalX);
		  }
		  if ( hitModule->subdet() == ENDCAP ) {
		    parametrizedResolutionLocalXEndcapsMap[myTag].Fill(cotAlpha, resolutionLocalX);
		    parametrizedResolutionLocalXEndcapsDistribution[myTag].Fill(resolutionLocalX);
		  }
		}
		// If any parameter for resolution on local Y coordinate specified for hitModule, fill maps and distributions
		if (hitModule->hasAnyResolutionLocalYParam()) {
		  double absCotBeta = fabs(1./tan(hitModule->beta(myTrack->getTheta())));
		  double resolutionLocalY = hitModule->resolutionLocalY(myTrack->getTheta())/Units::um; // um
		  if ( hitModule->subdet() == BARREL ) {
		    parametrizedResolutionLocalYBarrelMap[myTag].Fill(absCotBeta, resolutionLocalY);
		    parametrizedResolutionLocalYBarrelDistribution[myTag].Fill(resolutionLocalY);
		  }
		  if (hitModule->subdet() == ENDCAP ) {
		    parametrizedResolutionLocalYEndcapsMap[myTag].Fill(absCotBeta, resolutionLocalY);
		    parametrizedResolutionLocalYEndcapsDistribution[myTag].Fill(resolutionLocalY);
		  }
		}
	      }
	    }

	  } // hit loop
	} // track loop
      } // collection loop
    } // tag loop
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


    TriggerEfficiencyMapVisitor temv(myMap, myPt);
    SimParms::getInstance().accept(temv);
    tracker.accept(temv);
    temv.postVisit();
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


    PtThresholdMapVisitor ptmv(myMap, myEfficiency);
    SimParms::getInstance().accept(ptmv);
    tracker.accept(ptmv);
    ptmv.postVisit();

  }

  // Then: single maps

  for (int i=1; i<=suggestedSpacingMap.GetNbinsX(); ++i) {
    for (int j=1; j<=suggestedSpacingMap.GetNbinsY(); ++j) {
      suggestedSpacingMap.SetBinContent(i,j,0);
      suggestedSpacingMapAW.SetBinContent(i,j,0);
      nominalCutMap.SetBinContent(i,j,0);
    }
  }

  SpacingCutVisitor scv(suggestedSpacingMap, suggestedSpacingMapAW, nominalCutMap, moduleOptimalSpacings);
  SimParms::getInstance().accept(scv);
  tracker.accept(scv);
  scv.postVisit();
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
  // std::map<double, TGraph>& trigGraphs = myGraphBag.getGraphs(GraphBag::TriggerGraph|GraphBag::TriggeredGraph);
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
  //elemTotalGraph.first = GraphBag::Triggerable;
  //elemTotalGraph.second = totalGraph;
  elemTotalProfile.first = profileBag::Triggerable;
  elemTotalProfile.second = totalProfile;
  // Prepare plot: total trigger points
  //trigGraphs.insert(elemTotalGraph);
  trigProfiles.insert(elemTotalProfile);
  //trigGraphs[GraphBag::Triggerable].SetTitle("Average triggered points;#eta;Triggered points <N>");
  trigProfiles[profileBag::Triggerable].SetTitle("Average triggered points;#eta;Triggered points <N>");
  //aName.str(""); aName << "triggerable_vs_eta_graph";
  //trigGraphs[GraphBag::Triggerable].SetName(aName.str().c_str());
  aName.str(""); aName << "triggerable_vs_eta_profile";
  trigProfiles[profileBag::Triggerable].SetName(aName.str().c_str());

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
  // global
  rglobal.SetBins(bins, min, max);
  iglobal.SetBins(bins, min, max);
  // isolines
  isor.SetBins(bins, 0.0, geom_max_length, bins / 2, 0.0, geom_max_radius + geom_inactive_volume_width);
  isoi.SetBins(bins, 0.0, geom_max_length, bins / 2, 0.0, geom_max_radius + geom_inactive_volume_width);
  // Material distribution maps
  int materialMapBinsY = int( (geom_max_radius + geom_inactive_volume_width) * geom_safety_factor / 5.); // every half a cm
  int materialMapBinsX = int( (geom_max_length) * geom_safety_factor / 5.); // every half a cm
  mapRadiation.SetBins(       materialMapBinsX, 0.0, geom_max_length*geom_safety_factor, materialMapBinsY, 0.0, (geom_max_radius + geom_inactive_volume_width)*geom_safety_factor);
  mapInteraction.SetBins(     materialMapBinsX, 0.0, geom_max_length*geom_safety_factor, materialMapBinsY, 0.0, (geom_max_radius + geom_inactive_volume_width)*geom_safety_factor);
  mapRadiationCount.SetBins(  materialMapBinsX, 0.0, geom_max_length*geom_safety_factor, materialMapBinsY, 0.0, (geom_max_radius + geom_inactive_volume_width)*geom_safety_factor);
  mapInteractionCount.SetBins(materialMapBinsX, 0.0, geom_max_length*geom_safety_factor, materialMapBinsY, 0.0, (geom_max_radius + geom_inactive_volume_width)*geom_safety_factor);
  mapRadiationCalib.SetBins(  materialMapBinsX, 0.0, geom_max_length*geom_safety_factor, materialMapBinsY, 0.0, (geom_max_radius + geom_inactive_volume_width)*geom_safety_factor);
  mapInteractionCalib.SetBins(materialMapBinsX, 0.0, geom_max_length*geom_safety_factor, materialMapBinsY, 0.0, (geom_max_radius + geom_inactive_volume_width)*geom_safety_factor);
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


std::pair<double, double> Analyzer::computeMinMaxTracksEta(const Tracker& t) const {
  std::pair <double, double> etaMinMax = t.computeMinMaxEta();
  if (SimParms::getInstance().maxTracksEta.state()) etaMinMax.second = SimParms::getInstance().maxTracksEta();
  if (SimParms::getInstance().minTracksEta.state()) etaMinMax.first = SimParms::getInstance().minTracksEta();
  return etaMinMax;
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
  std::map <std::string, int> moduleTypeCount; // counts hit per module type -- if any of the sensors (or both) are hit, it counts 1
  std::map <std::string, int> sensorTypeCount; // counts hit per sensor on module type -- if one sensor is hit, it counts 1, if both sensors are hit, it counts 2
  std::map <std::string, int> moduleTypeCountStubs; // counts stubs per module type -- if both sensors are hit and a stub is formed, it counts 1
  std::map <std::string, TProfile> etaProfileByType;
  std::map <std::string, TProfile> etaProfileByTypeSensors;
  std::map <std::string, TProfile> etaProfileByTypeStubs;
//  TH2D* aPlot;
  std::string aType;


  double randomPercentMargin = 0.04;
  // Optimize the track creation on the real tracker
  auto etaMinMax = computeMinMaxTracksEta(tracker);
  // Computing the margin of the tracks to shoot
  double randomSpan = (etaMinMax.second - etaMinMax.first)*(1. + randomPercentMargin);
  double randomBase = etaMinMax.first - (etaMinMax.second - etaMinMax.first)*(randomPercentMargin)/2.;
  double maxEta = etaMinMax.second *= (1 + randomPercentMargin);


  // Initialize random number generator, counters and histograms
  myDice.SetSeed(MY_RANDOM_SEED);
  createResetCounters(tracker, moduleTypeCount);
  createResetCounters(tracker, sensorTypeCount);
  createResetCounters(tracker, moduleTypeCountStubs);

  class LayerNameVisitor : public ConstGeometryVisitor {
    string id_;
  public:
    std::set<string> data;
    LayerNameVisitor(const Tracker& t) { t.accept(*this); }
    void visit(const Barrel& b) { id_ = b.myid(); }
    void visit(const Endcap& e) { id_ = e.myid(); }
    void visit(const Layer& l) { data.insert(id_ + " " + any2str(l.myid())); }
    void visit(const Disk& d) { data.insert(id_ + " " + any2str(d.myid())); }
  };

  LayerNameVisitor layerNames(tracker);

  // Creating the layer hit coverage profiles
  layerEtaCoverageProfile.clear();
  layerEtaCoverageProfileStubs.clear();
  for (auto it = layerNames.data.begin(); it!=layerNames.data.end(); ++it) { // CUIDADO this is horrible code (not mine). refactor!
    TProfile* aProfile = new TProfile(Form("layerEtaCoverageProfileHits%s", it->c_str()), it->c_str(), 200, maxEta, maxEta);
    TProfile* aProfileStubs = new TProfile(Form("layerEtaCoverageProfileStubs%s", it->c_str()), it->c_str(), 200, maxEta, maxEta);
    layerEtaCoverageProfile[*it] = (*aProfile);
    layerEtaCoverageProfileStubs[*it] = (*aProfileStubs);
    delete aProfile;
    delete aProfileStubs;
  }


  etaProfileByType.clear();
  etaProfileByTypeSensors.clear();
  etaProfileByTypeStubs.clear();

  for (std::map <std::string, int>::iterator it = moduleTypeCount.begin();
       it!=moduleTypeCount.end(); it++) {
    TProfile& aProfile = etaProfileByType[(*it).first];
    aProfile.SetBins(100, 0, maxEta);
    aProfile.SetName((*it).first.c_str());
    aProfile.SetTitle((*it).first.c_str());
  }

  for (auto mel : sensorTypeCount) {
    TProfile& aProfileStubs = etaProfileByTypeSensors[mel.first];
    aProfileStubs.SetBins(100, 0, maxEta);
    aProfileStubs.SetName(mel.first.c_str());
    aProfileStubs.SetTitle(mel.first.c_str());
  }

  for (auto mel : moduleTypeCountStubs) {
    TProfile& aProfileStubs = etaProfileByTypeStubs[mel.first];
    aProfileStubs.SetBins(100, 0, maxEta);
    aProfileStubs.SetName(mel.first.c_str());
    aProfileStubs.SetTitle(mel.first.c_str());
  }

  //  ModuleVector allModules;
  double zError = SimParms::getInstance().zErrorCollider();

  // The real simulation
  std::pair <XYZVector, double> aLine;


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
  totalEtaProfile.SetTitle("Number of modules with at least one hit;#eta;Number of hit modules");
  totalEtaProfile.SetBins(100, 0, maxEta);
  totalEtaProfile.SetStats(0);

  totalEtaProfileSensors.Reset();
  totalEtaProfileSensors.SetName("totalEtaProfileSensors");
  totalEtaProfileSensors.SetMarkerStyle(8);
  totalEtaProfileSensors.SetMarkerColor(1);
  totalEtaProfileSensors.SetMarkerSize(1.5);
  totalEtaProfileSensors.SetTitle("Number of hits;#eta;Number of hits");
  totalEtaProfileSensors.SetBins(100, 0, maxEta);
  totalEtaProfileSensors.SetStats(0);

  totalEtaProfileStubs.Reset();
  totalEtaProfileStubs.SetName("totalEtaProfileStubs");
  totalEtaProfileStubs.SetMarkerStyle(8);
  totalEtaProfileStubs.SetMarkerColor(1);
  totalEtaProfileStubs.SetMarkerSize(1.5);
  totalEtaProfileStubs.SetTitle("Number of modules with a stub;#eta;Number of stubs");
  totalEtaProfileStubs.SetBins(100, 0, maxEta);
  totalEtaProfileStubs.SetStats(0);

  totalEtaProfileLayers.Reset();
  totalEtaProfileLayers.SetName("totalEtaProfileLayers");
  totalEtaProfileLayers.SetMarkerStyle(8);
  totalEtaProfileLayers.SetMarkerColor(1);
  totalEtaProfileLayers.SetMarkerSize(1.5);
  totalEtaProfileLayers.SetTitle("Number of layers with at least a hit;#eta;Number of layers");
  totalEtaProfileLayers.SetBins(100, 0, maxEta);
  totalEtaProfileLayers.SetStats(0);

  std::map<std::string, int> modulePlotColors; // CUIDADO quick and dirty way of creating a map with all the module colors (a cleaner way would be to have the map already created somewhere else)

  //XYZVector dir(0, 1, 0);
  // Shoot nTracksPerSide^2 tracks
  double angle = M_PI/2/(double)nTracksPerSide;
  for (int i=0; i<nTracksPerSide; i++) {
    for (int j=0; j<nTracksPerSide; j++) {
      // Reset the hit counter
      // Generate a straight track and collect the list of hit modules
      aLine = shootDirection(randomBase, randomSpan);
      std::vector<std::pair<Module*, HitType>> hitModules = trackHit( XYZVector(0, 0, ((myDice.Rndm()*2)-1)* zError), aLine.first, tracker.modules());
      // Reset the per-type hit counter and fill it
      resetTypeCounter(moduleTypeCount);
      resetTypeCounter(sensorTypeCount);
      resetTypeCounter(moduleTypeCountStubs);
      int numStubs = 0;
      int numHits = 0;
      for (auto& mh : hitModules) {
        moduleTypeCount[mh.first->moduleType()]++;
        if (mh.second & HitType::INNER) {
          sensorTypeCount[mh.first->moduleType()]++;
          numHits++;
        }
        if (mh.second & HitType::OUTER) {
          sensorTypeCount[mh.first->moduleType()]++;
          numHits++;
        }
        if (mh.second == HitType::STUB) {
          moduleTypeCountStubs[mh.first->moduleType()]++;
          numStubs++;
        }
        modulePlotColors[mh.first->moduleType()] = mh.first->plotColor();
      }
      // Fill the module type hit plot
      for (std::map <std::string, int>::iterator it = moduleTypeCount.begin(); it!=moduleTypeCount.end(); it++) {
        etaProfileByType[(*it).first].Fill(fabs(aLine.second), (*it).second);
      }

      for (auto& mel : sensorTypeCount) {
        etaProfileByTypeSensors[mel.first].Fill(fabs(aLine.second), mel.second);
      }

      for (auto& mel : moduleTypeCountStubs) {
        etaProfileByTypeStubs[mel.first].Fill(fabs(aLine.second), mel.second);
      }
      // Fill other plots
      mapPhiEta.Fill(aLine.first.Phi(), aLine.second, hitModules.size()); // phi, eta 2d plot
      mapPhiEtaCount.Fill(aLine.first.Phi(), aLine.second);               // Number of shot tracks

      totalEtaProfile.Fill(fabs(aLine.second), hitModules.size());                // Total number of hits
      totalEtaProfileSensors.Fill(fabs(aLine.second), numHits);
      totalEtaProfileStubs.Fill(fabs(aLine.second), numStubs); 

      int nHitLayers = 0;
      for (auto layerName : layerNames.data) {
        int layerHit = 0;
        int layerStub = 0;
        for (auto mh : hitModules) {
          UniRef ur = mh.first->uniRef();
          if (layerName == (ur.cnt + " " + any2str(ur.layer))) {
            layerHit=1;
            if (mh.second == HitType::STUB) layerStub=1;
            if (layerHit && layerStub) break;
          }
        }
        layerEtaCoverageProfile[layerName].Fill(aLine.second, layerHit);
        layerEtaCoverageProfileStubs[layerName].Fill(aLine.second, layerStub);
	if (layerHit) nHitLayers++;
      }

      totalEtaProfileLayers.Fill(fabs(aLine.second), nHitLayers);
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
      }
    }
  }

  savingGeometryV.push_back(mapPhiEta);

  // Eta profile compute
  //TProfile *myProfile;

  etaProfileCanvas.cd();
  savingGeometryV.push_back(etaProfileCanvas);

  //TProfile* total = total2D.ProfileX("etaProfileTotal");
  char profileName_[256];
  sprintf(profileName_, "etaProfileTotal%d", bsCounter++);
  // totalEtaProfile = TProfile(*total2D.ProfileX(profileName_));
  savingGeometryV.push_back(totalEtaProfile);
  if (totalEtaProfile.GetMaximum()<maximum_n_planes) totalEtaProfile.SetMaximum(maximum_n_planes);
  if (totalEtaProfileSensors.GetMaximum()<maximum_n_planes) totalEtaProfileSensors.SetMaximum(maximum_n_planes);
  if (totalEtaProfileStubs.GetMaximum()<maximum_n_planes) totalEtaProfileStubs.SetMaximum(maximum_n_planes);
  totalEtaProfile.Draw();
  totalEtaProfileSensors.Draw();
  totalEtaProfileStubs.Draw();
  for (std::map <std::string, TProfile>::iterator it = etaProfileByType.begin();
       it!=etaProfileByType.end(); it++) {
    TProfile* myProfile=(TProfile*)it->second.Clone();
    savingGeometryV.push_back(*myProfile); // TODO: remove savingGeometryV everywhere :-) [VERY obsolete...]
    myProfile->SetMarkerStyle(8);
    myProfile->SetMarkerColor(Palette::color(modulePlotColors[it->first]));
    myProfile->SetMarkerSize(1);
    std::string profileName = "etaProfile"+(*it).first;
    myProfile->SetName(profileName.c_str());
    myProfile->SetTitle((*it).first.c_str());
    myProfile->GetXaxis()->SetTitle("eta");
    myProfile->GetYaxis()->SetTitle("Number of hit modules");
    myProfile->Draw("same");
    typeEtaProfile.push_back(*myProfile);
  }

  for (std::map <std::string, TProfile>::iterator it = etaProfileByTypeSensors.begin();
       it!=etaProfileByTypeSensors.end(); it++) {
    TProfile* myProfile=(TProfile*)it->second.Clone();
    myProfile->SetMarkerStyle(8);
    myProfile->SetMarkerColor(Palette::color(modulePlotColors[it->first]));
    myProfile->SetMarkerSize(1);
    std::string profileName = "etaProfileSensors"+(*it).first;
    myProfile->SetName(profileName.c_str());
    myProfile->SetTitle((*it).first.c_str());
    myProfile->GetXaxis()->SetTitle("eta");
    myProfile->GetYaxis()->SetTitle("Number of hits");
    myProfile->Draw("same");
    typeEtaProfileSensors.push_back(*myProfile);
  }

  for (std::map <std::string, TProfile>::iterator it = etaProfileByTypeStubs.begin();
       it!=etaProfileByTypeStubs.end(); it++) {
    TProfile* myProfile=(TProfile*)it->second.Clone();
    myProfile->SetMarkerStyle(8);
    myProfile->SetMarkerColor(Palette::color(modulePlotColors[it->first]));
    myProfile->SetMarkerSize(1);
    std::string profileName = "etaProfileStubs"+(*it).first;
    myProfile->SetName(profileName.c_str());
    myProfile->SetTitle((*it).first.c_str());
    myProfile->GetXaxis()->SetTitle("eta");
    myProfile->GetYaxis()->SetTitle("Number of stubs");
    myProfile->Draw("same");
    typeEtaProfileStubs.push_back(*myProfile);
  }


  // Record the fraction of hits per module
  hitDistribution.SetBins(nTracks, 0 , 1);
  savingGeometryV.push_back(hitDistribution);
  for (auto m : tracker.modules()) {
    hitDistribution.Fill(m->numHits()/double(nTracks));
  }


  std::map<std::string, int> typeToCount;
  std::map<std::string, double> typeToPower;
  std::map<std::string, double> typeToSurface;

  for (auto aModule : tracker.modules()) {
    string aSensorType = aModule->moduleType();
    typeToCount[aSensorType] ++;
    typeToSurface[aSensorType] += aModule->area() / 1e6; // in mq
    typeToPower[aSensorType] += aModule->sensorsIrradiationPowerMean() / 1e3; // in kW
  }

  int iPoints = 0;
  for (std::map<std::string, int>::iterator typesIt = typeToCount.begin();
       typesIt != typeToCount.end(); ++typesIt ) {
    string aSensorType = typesIt->first;
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

      for (auto m : tracker.modules()) {
          aType = m->moduleType();
          m->resetHits();
          if (moduleTypeCount.find(aType)==moduleTypeCount.end()) {
            moduleTypeCount[aType]=typeCounter++;
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
    std::vector<std::pair<Module*, HitType>> Analyzer::trackHit(const XYZVector& origin, const XYZVector& direction, Tracker::Modules& moduleV) {
      std::vector<std::pair<Module*, HitType>> result;
      double distance;
      static const double BoundaryEtaSafetyMargin = 5. ; // track origin shift in units of zError to compute boundaries

      //static std::ofstream ofs("hits.txt");
      for (auto& m : moduleV) {
        // A module can be hit if it fits the phi (precise) contraints
        // and the eta constaints (taken assuming origin within 5 sigma)
        if (m->couldHit(direction, SimParms::getInstance().zErrorCollider()*BoundaryEtaSafetyMargin)) {
          auto h = m->checkTrackHits(origin, direction); 
          if (h.second != HitType::NONE) {
            result.push_back(std::make_pair(m,h.second));
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
      BandwidthVisitor bv(chanHitDistribution, bandwidthDistribution, bandwidthDistributionSparsified);
      SimParms::getInstance().accept(bv);
      tracker.accept(bv);
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


    std::map<int, TGraph>& Analyzer::getRhoGraphs(bool ideal, bool isTrigger) {
      int attribute = GraphBag::buildAttribute(ideal, isTrigger);
      attribute |= GraphBag::RhoGraph_Pt;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<int, TGraph>& Analyzer::getPhiGraphs(bool ideal, bool isTrigger) {
      int attribute = GraphBag::buildAttribute(ideal, isTrigger);
      attribute |= GraphBag::PhiGraph_Pt;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<int, TGraph>& Analyzer::getDGraphs(bool ideal, bool isTrigger) {
      int attribute = GraphBag::buildAttribute(ideal, isTrigger);
      attribute |= GraphBag::DGraph_Pt;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<int, TGraph>& Analyzer::getCtgThetaGraphs(bool ideal, bool isTrigger) {
      int attribute = GraphBag::buildAttribute(ideal, isTrigger);
      attribute |= GraphBag::CtgthetaGraph_Pt;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<int, TGraph>& Analyzer::getZ0Graphs(bool ideal, bool isTrigger) {
      int attribute = GraphBag::buildAttribute(ideal, isTrigger);
      attribute |= GraphBag::Z0Graph_Pt;
     return myGraphBag.getGraphs(attribute);
    }

    std::map<int, TGraph>& Analyzer::getPGraphs(bool ideal, bool isTrigger) {
      int attribute = GraphBag::buildAttribute(ideal, isTrigger);
      attribute |= GraphBag::PGraph_Pt;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<int, TGraph>& Analyzer::getLGraphs(bool ideal, bool isTrigger) {
      int attribute = GraphBag::buildAttribute(ideal, isTrigger);
      attribute |= GraphBag::LGraph_Pt;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<int, TGraph>& Analyzer::getBetaGraphs(bool ideal, bool isTrigger) {
      int attribute = GraphBag::buildAttribute(ideal, isTrigger);
      attribute |= GraphBag::BetaGraph_Pt;
      return myGraphBag.getGraphs(attribute);
    }

    std::map<int, TGraph>& Analyzer::getOmegaGraphs(bool ideal, bool isTrigger) {
      int attribute = GraphBag::buildAttribute(ideal, isTrigger);
      attribute |= GraphBag::OmegaGraph_Pt;
      return myGraphBag.getGraphs(attribute);
    }

    TH1D& Analyzer::getHistoOptimalSpacing(bool actualWindow) {
      if (actualWindow) return optimalSpacingDistributionAW;
      else return optimalSpacingDistribution;
    }

  }


#include "InnerCabling/InnerCablingMap.hh"
#include <Tracker.hh>


/*
 * KEY POINT: CREATE THE INNER TRACKER CABLING MAP.
 */
InnerCablingMap::InnerCablingMap(Tracker* tracker) {
  try {
    // CONNECT MODULES TO SERIAL POWER CHAINS
    connectModulesToPowerChains(tracker);

    // CONNECT MODULES TO GBTS
    connectModulesToGBTs(powerChains_, GBTs_);
    
    // CONNECT GBTS TO BUNDLES
    connectGBTsToBundles(GBTs_, bundles_);

    // CONNECT BUNDLES TO DTCs
    connectBundlesToDTCs(bundles_, DTCs_);
  }

  catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }
}


/* MODULES TO SERIAL POWER CHAINS CONNECTIONS.
 */
void InnerCablingMap::connectModulesToPowerChains(Tracker* tracker) {
  // build power chains
  ModulesToPowerChainsConnector powerChainsBuilder;
  tracker->accept(powerChainsBuilder);
  powerChainsBuilder.postVisit();

  // store power chains
  for (const auto& powerChainIt : powerChainsBuilder.getPowerChains()) {
    std::unique_ptr<PowerChain> myPowerChain(powerChainIt.second);
    powerChains_.insert(std::make_pair(powerChainIt.first, std::move(myPowerChain)));
  }
}


/* MODULES TO GBTS */

/* MODULES TO GBTS CONNECTIONS.
 * Something important is that the connections from modules to GBTs are based on the power chains mapping.
 * Indeed, all modules of a given GBT must belong to the same power chain.
 * As a result, one can just 'split' the modules of a given power chain and assign them to GBTs.
 */
void InnerCablingMap::connectModulesToGBTs(std::map<int, std::unique_ptr<PowerChain> >& powerChains, std::map<std::string, std::unique_ptr<GBT> >& GBTs) {

  // Loops on all power chains
  for (auto& it : powerChains) {

    // COLLECT GENERAL INFORMATION NEEDED TO BUILD GBTS   
    PowerChain* myPowerChain = it.second.get();

    const bool isBarrel = myPowerChain->isBarrel();
    const std::string subDetectorName = myPowerChain->subDetectorName();
    const int layerNumber = myPowerChain->layerDiskNumber();
    const int ringNumber = myPowerChain->ringNumber();
    const int numModulesInPowerChain = myPowerChain->numModules();

    const int layerOrRingNumber = (isBarrel ? layerNumber : ringNumber);
    const int numELinksPerModule = inner_cabling_functions::computeNumELinksPerModule(subDetectorName, layerOrRingNumber);

    const std::pair<int, double> gbtsInPowerChain = computeNumGBTsInPowerChain(numELinksPerModule, numModulesInPowerChain, isBarrel);
    const int numGBTsInPowerChain = gbtsInPowerChain.first;
    const double numModulesPerGBTExact = gbtsInPowerChain.second;

    const int powerChainId = myPowerChain->myid();
    const bool isBarrelLong = myPowerChain->isBarrelLong();
    

    // Loops on all modules of the power chain
    for (auto& m : myPowerChain->modules()) {

      // SET NUMBER OF ELINKS PER MODULE
      m->setNumELinks(numELinksPerModule);

      // COLLECT MODULE INFORMATION NEEDED TO BUILD GBT
      const int ringRef = (isBarrelLong ? m->uniRef().ring - 1 : m->uniRef().ring - 2);
      const int phiRefInPowerChain = m->getPhiRefInPowerChain();
      
      const std::pair<int, int> myGBTIndexes = computeGBTPhiIndex(isBarrel, ringRef, phiRefInPowerChain, numModulesPerGBTExact, numGBTsInPowerChain);
      const int myGBTIndex = myGBTIndexes.first;
      const int myGBTIndexColor = myGBTIndexes.second;
      const std::string myGBTId = computeGBTId(powerChainId, myGBTIndex);

      // BUILD GBTS AND STORE THEM
      createAndStoreGBTs(myPowerChain, m, myGBTId, myGBTIndex, myGBTIndexColor, numELinksPerModule, GBTs);
    }
  }

  // CHECK GBTS
  checkModulesToGBTsCabling(GBTs);
}


/*
 * This is used to compute the number of GBTs in one power chain.
 * This also provides the exact average number of modules per GBT, in that power chain.
 * This is obviously based on the number of ELinks the modules are connected to, and the total number of modules in the power chain.
 */
const std::pair<int, double> InnerCablingMap::computeNumGBTsInPowerChain(const int numELinksPerModule, const int numModulesInPowerChain, const bool isBarrel) {

  int numModules = numModulesInPowerChain;
  if (isBarrel) {
    if (numModules % 2 == 1) logERROR(any2str("Found an odd number of modules in BPIX power chain, which is not supported."));
    else numModules /= 2;  // Divide by 2 because in BPIX, the GBTs assignment works by rod.
                           // This is becasue it makes the powering of the GBTs much easier.
  }

  const double maxNumModulesPerGBTExact = static_cast<double>(inner_cabling_maxNumELinksPerGBT) / numELinksPerModule;
  const int maxNumModulesPerGBT = (fabs(maxNumModulesPerGBTExact - round(maxNumModulesPerGBTExact)) < inner_cabling_roundingTolerance ? 
					       round(maxNumModulesPerGBTExact) 
					       : std::floor(maxNumModulesPerGBTExact)
					       );

  const double numGBTsExact = static_cast<double>(numModules) / maxNumModulesPerGBT;
  const int numGBTs = (fabs(numGBTsExact - round(numGBTsExact)) < inner_cabling_roundingTolerance ? 
		       round(numGBTsExact) 
		       : std::ceil(numGBTsExact)
		       );

  if (numGBTs == 0) logERROR(any2str("Power chain has ") + any2str(numModules) 
			     + any2str(" modules, but found numGBTs == ") +  any2str(numGBTs) + any2str(", that's not enough!!")
			     );

  const double numModulesPerGBTExact = static_cast<double>(numModules) / numGBTs;

  return std::make_pair(numGBTs, numModulesPerGBTExact);
}


/* Compute the phi index associated to each GBT.
 */
const std::pair<int, int> InnerCablingMap::computeGBTPhiIndex(const bool isBarrel, const int ringRef, const int phiRefInPowerChain, const double numModulesPerGBTExact, const int numGBTsInPowerChain) const {

  const int moduleRef = (isBarrel ? ringRef : phiRefInPowerChain);

  if (maxNumModulesPerGBTInPowerChain == 0) logERROR(any2str("Found maxNumModulesPerGBTInPowerChain == 0."));

  const double myGBTIndexExact = moduleRef / numModulesPerGBTExact;
  int myGBTIndex = (fabs(myGBTIndexExact - round(myGBTIndexExact)) < inner_cabling_roundingTolerance ? 
		    round(myGBTIndexExact) 
		    : std::floor(myGBTIndexExact)
		    );
  if (isBarrel && phiRefInPowerChain == 1) myGBTIndex += numGBTsInPowerChain;

  int myGBTIndexColor = myGBTIndex;
  if (isBarrel && phiRefInPowerChain == 1 && femod(numGBTsInPowerChain, 2) == 0) myGBTIndexColor += 1;
  myGBTIndexColor = femod(myGBTIndexColor, 2);

  return std::make_pair(myGBTIndex, myGBTIndexColor);
}


/* Compute the Id associated to each GBT.
 */
const std::string InnerCablingMap::computeGBTId(const int powerChainId, const int myGBTIndex) const {
  std::ostringstream GBTIdStream;
  GBTIdStream << powerChainId << "_" << myGBTIndex;
  const std::string GBTId = GBTIdStream.str();
  return GBTId;
}


/* Create a GBT, if does not exist yet.
 * Store it in the GBTs container.
 */
void InnerCablingMap::createAndStoreGBTs(PowerChain* myPowerChain, Module* m, const std::string myGBTId, const int myGBTIndex, const int myGBTIndexColor, const int numELinksPerModule, std::map<std::string, std::unique_ptr<GBT> >& GBTs) {

  auto found = GBTs.find(myGBTId);
  if (found == GBTs.end()) {
    std::unique_ptr<GBT> myGBT(new GBT(myPowerChain, myGBTId, myGBTIndex, myGBTIndexColor, numELinksPerModule));
    connectOneModuleToOneGBT(m, myGBT.get());
    GBTs.insert(std::make_pair(myGBTId, std::move(myGBT)));  
  }
  else {
    connectOneModuleToOneGBT(m, found->second.get());
  }
}


/* Connect module to GBT and vice-versa.
 */
void InnerCablingMap::connectOneModuleToOneGBT(Module* m, GBT* myGBT) const {
  myGBT->addModule(m);
  m->setGBT(myGBT);
}


/* Check bundles-cables connections.
 */

void InnerCablingMap::checkModulesToGBTsCabling(const std::map<std::string, std::unique_ptr<GBT> >& GBTs) const {
  for (const auto& it : GBTs) {
    const std::string myGBTId = it.first;
    const GBT* myGBT = it.second.get();

    // CHECK THE NUMBER OF ELINKS PER GBT.
    const int numELinks = myGBT->numELinks();

    if (numELinks > inner_cabling_maxNumELinksPerGBT) {
      logERROR(any2str("GBT ")  + any2str(myGBTId) + any2str(" is connected to ") + any2str(numELinks) + any2str(" ELinks. ")
	       + "Maximum number of ELinks per GBT allowed is " + any2str(inner_cabling_maxNumELinksPerGBT)
	       );
    }
  }
}


/* GBTs to BUNDLES */

/* GBTS TO BUNDLES CONNECTIONS.
 */
void InnerCablingMap::connectGBTsToBundles(std::map<std::string, std::unique_ptr<GBT> >& GBTs, std::map<int, std::unique_ptr<InnerBundle> >& bundles) {

  for (auto& it : GBTs) {
    // COLLECT ALL INFORMATION NEEDED TO BUILD BUNDLES   
    GBT* myGBT = it.second.get();

    const std::string subDetectorName = myGBT->subDetectorName();
    
    const int layerDiskNumber = myGBT->layerDiskNumber();
    const int powerChainPhiRef = myGBT->powerChainPhiRef();
    const int ringNumber = myGBT->ringNumber();
        
    const int myBundleIndex = computeBundleIndex(subDetectorName, layerDiskNumber, powerChainPhiRef, ringNumber);

    const bool isPositiveZEnd = myGBT->isPositiveZEnd();
    const bool isPositiveXSide = myGBT->isPositiveXSide();
    const int myBundleId = computeBundleId(isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber, myBundleIndex);

    // BUILD BUNDLES AND STORE THEM
    createAndStoreBundles(myGBT, bundles, myBundleId, isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber, myBundleIndex);    
  }

  // CHECK BUNDLES
  checkGBTsToBundlesCabling(bundles);
}


/*
 * Compute the index used to identify uniquely a Fiber Bundle.
 */
const int InnerCablingMap::computeBundleIndex(const std::string subDetectorName, const int layerNumber, const int powerChainPhiRef, const int ringNumber) const {
  int myBundleIndex = 0;

  if (subDetectorName == inner_cabling_tbpx) {
    int maxNumPowerChainsPerBundleBarrelLayer = 0;
    if (layerNumber == 1) maxNumPowerChainsPerBundleBarrelLayer = maxNumPowerChainsPerBundleBarrelLayer1;
    else if (layerNumber == 2) maxNumPowerChainsPerBundleBarrelLayer = maxNumPowerChainsPerBundleBarrelLayer2;
    else if (layerNumber == 3) maxNumPowerChainsPerBundleBarrelLayer = maxNumPowerChainsPerBundleBarrelLayer3;
    else if (layerNumber == 4) maxNumPowerChainsPerBundleBarrelLayer = maxNumPowerChainsPerBundleBarrelLayer4;
    else logERROR("Did not find supported layer number.");

    if (maxNumPowerChainsPerBundleBarrelLayer == 0) logERROR(any2str("Found maxNumPowerChainsPerBundleBarrelLayer == 0."));

    const double myBundleIndexExact = static_cast<double>(powerChainPhiRef) / maxNumPowerChainsPerBundleBarrelLayer;
    myBundleIndex = (fabs(myBundleIndexExact - round(myBundleIndexExact)) < inner_cabling_roundingTolerance ? 
		     round(myBundleIndexExact) 
		     : std::floor(myBundleIndexExact)
		     );
  }
  else if (subDetectorName == inner_cabling_tfpx || subDetectorName == inner_cabling_tepx) {
    myBundleIndex = (femod(ringNumber, 2) == 1 ? 0 : 1);
  }
  else logERROR("Unsupported detector name.");

  return myBundleIndex;
}


/* Compute the Id associated to each bundle.
 */
const int InnerCablingMap::computeBundleId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int myBundleIndex) const {

  const int innerTrackerQuarterIndex = inner_cabling_functions::computeInnerTrackerQuarterIndex(isPositiveZEnd, isPositiveXSide);
  const int subdetectorIndex = inner_cabling_functions::computeSubDetectorIndex(subDetectorName);

  const int bundleId = innerTrackerQuarterIndex * 1000 + subdetectorIndex * 100 + layerDiskNumber * 10 + myBundleIndex;
  return bundleId;
}


/* Create a Bundle, if does not exist yet.
 * Store it in the Bundles container.
 */
void InnerCablingMap::createAndStoreBundles(GBT* myGBT, std::map<int, std::unique_ptr<InnerBundle> >& bundles, const int bundleId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int myBundleIndex) {

  auto found = bundles.find(bundleId);
  if (found == bundles.end()) {
    std::unique_ptr<InnerBundle> myBundle(new InnerBundle(bundleId, isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber, myBundleIndex));
    connectOneGBTToOneBundle(myGBT, myBundle.get());
    bundles.insert(std::make_pair(bundleId, std::move(myBundle)));  
  }
  else {
    connectOneGBTToOneBundle(myGBT, found->second.get());
  }
}


/* Connect GBT to Bundle and vice-versa.
 */
void InnerCablingMap::connectOneGBTToOneBundle(GBT* myGBT, InnerBundle* myBundle) const {
  myBundle->addGBT(myGBT);
  myGBT->setBundle(myBundle);
}


/* Check GBTs-Bundle connections.
 */
void InnerCablingMap::checkGBTsToBundlesCabling(const std::map<int, std::unique_ptr<InnerBundle> >& bundles) const {
  for (const auto& it : bundles) {
    const int myBundleId = it.first;
    const InnerBundle* myBundle = it.second.get();

    // CHECK THE NUMBER OF GBTs per Bundle
    const int numGBTs = myBundle->numGBTs();

    if (numGBTs > inner_cabling_maxNumGBTsPerBundle) {
      logERROR(any2str("InnerBundle ")  + any2str(myBundleId) + any2str(" is connected to ") + any2str(numGBTs) + any2str(" GBTs. ")
	       + "Maximum number of GBTs per Bundle allowed is " + any2str(inner_cabling_maxNumGBTsPerBundle)
	       );
    }
  }
}


/* BUNDLES TO DTCS */

/* BUNDLES TO DTCS CONNECTIONS.
 */
void InnerCablingMap::connectBundlesToDTCs(std::map<int, std::unique_ptr<InnerBundle> >& bundles, std::map<int, std::unique_ptr<InnerDTC> >& DTCs) {

  for (auto& it : bundles) {
    // COLLECT ALL INFORMATION NEEDED TO BUILD DTCS   
    InnerBundle* myBundle = it.second.get();

    const std::string subDetectorName = myBundle->subDetectorName();   
    const int layerDiskNumber = myBundle->layerDiskNumber();
        
    const bool isPositiveZEnd = myBundle->isPositiveZEnd();
    const bool isPositiveXSide = myBundle->isPositiveXSide();
    const int myDTCId = computeDTCId(isPositiveZEnd, isPositiveXSide, subDetectorName, layerDiskNumber);

    // BUILD DTCS AND STORE THEM
    createAndStoreDTCs(myBundle, DTCs, myDTCId, isPositiveZEnd, isPositiveXSide);    
  }

  // CHECK DTCS
  checkBundlesToDTCsCabling(DTCs);
}


/* Compute the Id associated to each DTC.
 */
const int InnerCablingMap::computeDTCId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber) const {
  int myDTCId = 0;

  if (subDetectorName == inner_cabling_tbpx) myDTCId = layerDiskNumber;

  else if (subDetectorName == inner_cabling_tfpx) {
    if (layerDiskNumber == 1) myDTCId = 5;
    else if (layerDiskNumber == 2) myDTCId = 4;
    else if (layerDiskNumber == 3) myDTCId = 3;
    else if (layerDiskNumber == 4) myDTCId = 5;
    else if (layerDiskNumber == 5) myDTCId = 4;
    else if (layerDiskNumber == 6) myDTCId = 3;
    else if (layerDiskNumber == 7) myDTCId = 2;
    else if (layerDiskNumber == 8) myDTCId = 5;
    else logERROR(any2str("Unexpected diskNumber in FPX : ") + any2str(layerDiskNumber));
  }

  else if (subDetectorName == inner_cabling_tepx) {
    if (layerDiskNumber == 1) myDTCId = 6;
    else if (layerDiskNumber == 2) myDTCId = 6;
    else if (layerDiskNumber == 3) myDTCId = 7;
    else if (layerDiskNumber == 4) myDTCId = 7;
    else logERROR(any2str("Unexpected diskNumber in EPX : ") + any2str(layerDiskNumber));
  }

  const int innerTrackerQuarterIndex = inner_cabling_functions::computeInnerTrackerQuarterIndex(isPositiveZEnd, isPositiveXSide);
  myDTCId += innerTrackerQuarterIndex * 10;

  return myDTCId;
}


/* Create a DTC, if does not exist yet.
 * Store it in the DTC container.
 */
void InnerCablingMap::createAndStoreDTCs(InnerBundle* myBundle, std::map<int, std::unique_ptr<InnerDTC> >& DTCs, const int DTCId, const bool isPositiveZEnd, const bool isPositiveXSide) {

  auto found = DTCs.find(DTCId);
  if (found == DTCs.end()) {
    std::unique_ptr<InnerDTC> myDTC(new InnerDTC(DTCId, isPositiveZEnd, isPositiveXSide));
    connectOneBundleToOneDTC(myBundle, myDTC.get());
    DTCs.insert(std::make_pair(DTCId, std::move(myDTC)));  
  }
  else {
    connectOneBundleToOneDTC(myBundle, found->second.get());
  }
}


/* Connect Bundle to DTC and vice-versa.
 */
void InnerCablingMap::connectOneBundleToOneDTC(InnerBundle* myBundle, InnerDTC* myDTC) const {
  myDTC->addBundle(myBundle);
  myBundle->setDTC(myDTC);
}


/* Check Bundles-DTC connections.
 */

void InnerCablingMap::checkBundlesToDTCsCabling(const std::map<int, std::unique_ptr<InnerDTC> >& DTCs) const {
  for (const auto& it : DTCs) {
    const int myDTCId = it.first;
    const InnerDTC* myDTC = it.second.get();

    // CHECK THE NUMBER OF Bundles per DTC
    const int numBundles = myDTC->numBundles();

    if (numBundles > inner_cabling_maxNumBundlesPerCable) {
      logERROR(any2str("InnerDTC ")  + any2str(myDTCId) + any2str(" is connected to ") + any2str(numBundles) + any2str(" Bundles. ")
	       + "Maximum number of Bundles per DTC allowed is " + any2str(inner_cabling_maxNumBundlesPerCable)
	       );
    }
  }
}


#ifndef INNERCABLINGMAP_HH
#define INNERCABLINGMAP_HH

#include <global_constants.hh>
#include "global_funcs.hh"

#include "Property.hh"
#include "ITCabling/ModulesToPowerChainsConnector.hh"
#include "ITCabling/ModulesToELinksConnector.hh"
#include "ITCabling/InnerBundle.hh"
//#include "Cabling/DTC.hh"


/* Build the cabling map.
   The cabling map contains all the necessary bundles, cables and DTCs for a given tracker.
   All the relevant connections are made during the building process.
 */
class InnerCablingMap : public PropertyObject, public Buildable, public Identifiable<int> {
public:
  InnerCablingMap(Tracker* tracker);

  // positive cabling side
  const std::map<int, PowerChain*> getPowerChains() const { return powerChains_; }
  const std::map<std::string, ELink*> getELinks() const { return eLinks_; }
  const std::map<std::string, GBT*> getGBTs() const { return GBTs_; }
  const std::map<int, InnerBundle*>& getBundles() const { return bundles_; }
  //const std::map<int, Cable*>& getCables() const { return cables_; }
  //const std::map<const std::string, const DTC*>& getDTCs() const { return DTCs_; }

  // negative cabling side
  //const std::map<int, Bundle*>& getNegBundles() const { return negBundles_; }
  //const std::map<int, Cable*>& getNegCables() const { return negCables_; }
  //const std::map<const std::string, const DTC*>& getNegDTCs() const { return negDTCs_; }

private:
  // CONNECT MODULES TO POWER CHAINS
  void connectModulesToPowerChains(Tracker* tracker);
  void connectModulesToELinks(Tracker* tracker);
  void connectELinksToGBTs(std::map<int, PowerChain*>& powerChains, std::map<std::string, GBT*>& GBTs);
  void connectGBTsToBundles(std::map<std::string, GBT*>& GBTs, std::map<int, InnerBundle*>& bundles);



  

  // MODULES TO GBTS
  const std::pair<int, int> computeMaxNumModulesPerGBTInPowerChain(const int numELinksPerModule, const int numModulesInPowerChain, const bool isBarrel);
  const int computeGBTPhiIndex(const bool isBarrel, const int ringRef, const int phiRefInPowerChain, const int maxNumModulesPerGBTInPowerChain, const int numGBTsInPowerChain) const;
  const std::string computeGBTId(const int powerChainId, const int myGBTIndex) const;
  void createAndStoreGBTs(PowerChain* myPowerChain, Module& m, const std::string myGBTId, const int myGBTPhiIndex, const int numELinksPerModule, std::map<std::string, GBT*>& GBTs);
  void connectOneModuleToOneGBT(Module& m, GBT* GBT) const;
  void checkModulesToGBTsCabling(const std::map<std::string, GBT*>& GBTs) const;

  // GBTs to BUNDLES
  const int computeBundleIndex(const std::string subDetectorName, const int layerNumber, const int powerChainPhiRef, const bool isRingInnerEnd) const;
  const int computeBundleId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int myBundleIndex) const;
  void createAndStoreBundles(GBT* myGBT, std::map<int, InnerBundle*>& bundles, const int bundleId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int myBundleIndex);
  void connectOneGBTToOneBundle(GBT* myGBT, InnerBundle* myBundle) const;
  void checkGBTsToBundlesCabling(const std::map<int, InnerBundle*>& bundles) const;



  /*
  // CONNECT BUNDLES TO CABLES
  void connectBundlesToCables(std::map<int, Bundle*>& bundles, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs);
  const Category computeCableType(const Category& bundleType) const;
  const std::map<int, std::pair<int, int> > computeCablesPhiSectorRefAndSlot(const std::map<int, Bundle*>& bundles) const;
  const int computeCableTypeIndex(const Category& cableType) const;
  const int computeCableId(const int phiSectorRefCable, const int cableTypeIndex, const int slot, const bool isPositiveCablingSide) const;
  void createAndStoreCablesAndDTCs(Bundle* myBundle, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs, const int cableId, const double phiSectorWidth, const int phiSectorRefCable, const Category& type, const int slot, const bool isPositiveCablingSide); 
  void connectOneBundleToOneCable(Bundle* bundle, Cable* cable) const;
  void checkBundlesToCablesCabling(std::map<int, Cable*>& cables);  // check bundles to cables connections

  // COMPUTE SERVICES CHANNELS ASSIGNMENTS OF POWER CABLES
  void computePowerServicesChannels();
  void routeBarrelBundlesPoweringToSemiNonants(const bool isPositiveCablingSide);
  void checkBundlesToPowerServicesChannels(const std::map<int, Bundle*>& bundles);
  */

  // positive cabling side
  std::map<int, PowerChain*> powerChains_;
  std::map<std::string, ELink*> eLinks_;
  std::map<std::string, GBT*> GBTs_;
  std::map<int, InnerBundle*> bundles_;
  //std::map<const std::string, const DTC*> DTCs_;

  // negative cabling side
  //std::map<int, Bundle*> negBundles_;
  //std::map<int, Cable*> negCables_;
  //std::map<const std::string, const DTC*> negDTCs_;
};


#endif  // INNERCABLINGMAP_HH

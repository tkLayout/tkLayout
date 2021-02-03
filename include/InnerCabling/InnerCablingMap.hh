#ifndef INNERCABLINGMAP_HH
#define INNERCABLINGMAP_HH

#include <global_constants.hh>
#include "global_funcs.hh"

#include "Property.hh"

#include "InnerCabling/ModulesToPowerChainsConnector.hh"
#include "InnerCabling/InnerDTC.hh"


/* Build the IT cabling map.
   The cabling map contains all the mapping of the GBTs, Bundles, and DTCs with respect to each other.
   It also contains the Serial Power Chains mapping, as the readout mapping cannot be distangled from the powering.
   All the relevant connections are made during the building process.
 */
class InnerCablingMap : public PropertyObject, public Buildable, public Identifiable<int> {
public:
  // KEY POINT: CREATE THE INNER TRACKER CABLING MAP.
  InnerCablingMap(Tracker* tracker);

  const std::map<int, std::unique_ptr<PowerChain> >& getPowerChains() const { return powerChains_; }
  const std::map<std::string, std::unique_ptr<GBT> >& getGBTs() const { return GBTs_; }
  const std::map<int, std::unique_ptr<InnerBundle> >& getBundles() const { return bundles_; }
  const std::map<int, std::unique_ptr<InnerDTC> >& getDTCs() const { return DTCs_; }


private:
  // CONNECT MODULES TO POWER CHAINS
  void connectModulesToPowerChains(Tracker* tracker);
  void connectModulesToGBTs(std::map<int, std::unique_ptr<PowerChain> >& powerChains, std::map<std::string, std::unique_ptr<GBT> >& GBTs);
  void connectGBTsToBundles(std::map<std::string, std::unique_ptr<GBT> >& GBTs, std::map<int, std::unique_ptr<InnerBundle> >& bundles);
  void connectBundlesToDTCs(std::map<int, std::unique_ptr<InnerBundle> >& bundles, std::map<int, std::unique_ptr<InnerDTC> >& DTCs);

  // CONNECT MODULES TO GBTS
  const std::pair<int, double> computeNumGBTsInPowerChain(const int numELinksPerModule, const int numModulesInPowerChain, const bool isBarrel);
  const std::pair<int, int> computeGBTPhiIndex(const bool isBarrel, const int ringRef, const int phiRefInPowerChain, const double maxNumModulesPerGBTInPowerChain, const int numGBTsInPowerChain) const;
  const std::string computeGBTId(const int powerChainId, const int myGBTIndex) const;
  void createAndStoreGBTs(PowerChain* myPowerChain, Module* m, const std::string myGBTId, const int myGBTIndex, const int myGBTIndexColor, const int numELinksPerModule, std::map<std::string, std::unique_ptr<GBT> >& GBTs);
  void connectOneModuleToOneGBT(Module* m, GBT* GBT) const;
  void checkModulesToGBTsCabling(const std::map<std::string, std::unique_ptr<GBT> >& GBTs) const;

  // CONNECT GBTs TO BUNDLES
  const int computeBundleIndex(const std::string subDetectorName, const int layerNumber, const int powerChainPhiRef, const bool isAtSmallerAbsZDeeInDoubleDisk) const;
  const int computeBundleId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int myBundleIndex) const;
  void createAndStoreBundles(GBT* myGBT, std::map<int, std::unique_ptr<InnerBundle> >& bundles, const int bundleId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int myBundleIndex);
  void connectOneGBTToOneBundle(GBT* myGBT, InnerBundle* myBundle) const;
  void checkGBTsToBundlesCabling(const std::map<int, std::unique_ptr<InnerBundle> >& bundles) const;

  // CONNECT BUNDLES TO DTCS
  const int computeDTCId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber) const;
  void createAndStoreDTCs(InnerBundle* myBundle, std::map<int, std::unique_ptr<InnerDTC> >& DTCs, const int DTCId, const bool isPositiveZEnd, const bool isPositiveXSide);
  void connectOneBundleToOneDTC(InnerBundle* myBundle, InnerDTC* myDTC) const;
  void computeCMSSWIds(std::map<int, std::unique_ptr<InnerDTC> >& DTCs);
  void checkBundlesToDTCsCabling(const std::map<int, std::unique_ptr<InnerDTC> >& DTCs) const;


  std::map<int, std::unique_ptr<PowerChain> > powerChains_;
  std::map<std::string, std::unique_ptr<GBT> > GBTs_;
  std::map<int, std::unique_ptr<InnerBundle> > bundles_;
  std::map<int, std::unique_ptr<InnerDTC> > DTCs_;
  // All power chains, GBTs, bundles, and DTC are owned by the Cabling map, and the Cabling map only!!
};


#endif  // INNERCABLINGMAP_HH

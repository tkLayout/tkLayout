#ifndef CABLINGMAP_HH
#define CABLINGMAP_HH

#include <global_constants.hh>
#include "global_funcs.hh"

#include "Property.hh"
#include "OuterCabling/ModulesToBundlesConnector.hh"
#include "OuterCabling/DTC.hh"


/* Build the cabling map.
   The cabling map contains all the necessary bundles, cables and DTCs for a given tracker.
   All the relevant connections are made during the building process.
 */
class CablingMap : public PropertyObject, public Buildable, public Identifiable<int> {
public:
  CablingMap(Tracker* tracker);

  // positive cabling side
  const std::map<int, OuterBundle*>& getBundles() const { return bundles_; }
  const std::map<int, Cable*>& getCables() const { return cables_; }
  const std::map<const std::string, const DTC*>& getDTCs() const { return DTCs_; }

  // negative cabling side
  const std::map<int, OuterBundle*>& getNegBundles() const { return negBundles_; }
  const std::map<int, Cable*>& getNegCables() const { return negCables_; }
  const std::map<const std::string, const DTC*>& getNegDTCs() const { return negDTCs_; }

private:
  // CONNECT MODULES TO BUNDLES
  void connectModulesToBundles(Tracker* tracker);

  // CONNECT BUNDLES TO CABLES
  void connectBundlesToCables(std::map<int, OuterBundle*>& bundles, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs);
  const Category computeCableType(const Category& bundleType) const;
  const std::map<int, std::pair<int, int> > computeCablesPhiSectorRefAndSlot(const std::map<int, OuterBundle*>& bundles) const;
  const int computeCableTypeIndex(const Category& cableType) const;
  const int computeCableId(const int phiSectorRefCable, const int cableTypeIndex, const int slot, const bool isPositiveCablingSide) const;
  void createAndStoreCablesAndDTCs(OuterBundle* myBundle, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs, const int cableId, const double phiSectorWidth, const int phiSectorRefCable, const Category& type, const int slot, const bool isPositiveCablingSide); 
  void connectOneBundleToOneCable(OuterBundle* bundle, Cable* cable) const;
  void checkBundlesToCablesCabling(std::map<int, Cable*>& cables);  // check bundles to cables connections

  // COMPUTE SERVICES CHANNELS ASSIGNMENTS OF POWER CABLES
  void computePowerServicesChannels();
  void routeBarrelBundlesPoweringToSemiNonants(const bool isPositiveCablingSide);
  void checkBundlesToPowerServicesChannels(const std::map<int, OuterBundle*>& bundles);

  // positive cabling side
  std::map<int, OuterBundle*> bundles_;
  std::map<int, Cable*> cables_;
  std::map<const std::string, const DTC*> DTCs_;

  // negative cabling side
  std::map<int, OuterBundle*> negBundles_;
  std::map<int, Cable*> negCables_;
  std::map<const std::string, const DTC*> negDTCs_;
};


#endif  // CABLINGMAP_HH

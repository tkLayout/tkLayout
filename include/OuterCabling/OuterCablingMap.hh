#ifndef OUTERCABLINGMAP_HH
#define OUTERCABLINGMAP_HH

#include <global_constants.hh>
#include "global_funcs.hh"

#include "Property.hh"
#include "OuterCabling/ModulesToBundlesConnector.hh"
#include "OuterCabling/OuterDTC.hh"

namespace insur {

/* Build the cabling map.
   The cabling map contains all the necessary bundles, cables and DTCs for a given tracker.
   All the relevant connections are made during the building process.
 */
class OuterCablingMap : public PropertyObject, public Buildable, public Identifiable<int> {
public:
  OuterCablingMap(Tracker* tracker);

  // positive cabling side
  const std::map<const int, std::unique_ptr<OuterBundle> >& getBundles() const { return bundles_; }
  const std::map<const int, std::unique_ptr<OuterCable> >& getCables() const { return cables_; }
  const std::map<const std::string, std::unique_ptr<const OuterDTC> >& getDTCs() const { return DTCs_; }

  // negative cabling side
  const std::map<const int, std::unique_ptr<OuterBundle> >& getNegBundles() const { return negBundles_; }
  const std::map<const int, std::unique_ptr<OuterCable> >& getNegCables() const { return negCables_; }
  const std::map<const std::string, std::unique_ptr<const OuterDTC> >& getNegDTCs() const { return negDTCs_; }
  const std::map<const int, std::unique_ptr<OuterGBT> >& getGBTs() const { return gbts_; }
  

private:
  // CONNECT MODULES TO BUNDLES
  void connectModulesToBundles(Tracker* tracker);

  // CONNECT BUNDLES TO CABLES
  void connectBundlesToCables(std::map<const int, std::unique_ptr<OuterBundle> >& bundles, std::map<const int, std::unique_ptr<OuterCable> >& cables, std::map<const std::string, std::unique_ptr<const OuterDTC> >& DTCs);
  const Category computeCableType(const Category& bundleType) const;
  const std::map<int, std::pair<int, int> > computeCablesPhiSectorRefAndSlot(const std::map<const int, std::unique_ptr<OuterBundle> >& bundles) const;
  const int computeCableTypeIndex(const Category& cableType) const;
  const int computeCableId(const int phiSectorRefCable, const int cableTypeIndex, const int slot, const bool isPositiveCablingSide) const;
  void createAndStoreCablesAndDTCs(OuterBundle* myBundle, std::map<const int, std::unique_ptr<OuterCable> >& cables, std::map<const std::string, std::unique_ptr<const OuterDTC> >& DTCs, const int cableId, const double phiSectorWidth, const int phiSectorRefCable, const Category& type, const int slot, const bool isPositiveCablingSide); 
  void connectOneBundleToOneCable(OuterBundle* bundle, OuterCable* cable) const;
  void checkBundlesToCablesCabling(const std::map<const int, std::unique_ptr<OuterCable> >& cables);  // check bundles to cables connections
  void computeCMSSWIds(std::map<const std::string, std::unique_ptr<const OuterDTC> >& DTCs);

  // COMPUTE SERVICES CHANNELS ASSIGNMENTS OF POWER CABLES
  void computePowerServicesChannels();
  void routeBarrelBundlesPoweringToSemiNonants(const bool isPositiveCablingSide);
  void checkBundlesToPowerServicesChannels(const std::map<const int, std::unique_ptr<OuterBundle> >& bundles);

  // positive cabling side
  std::map<const int, std::unique_ptr<OuterBundle> > bundles_;
  std::map<const int, std::unique_ptr<OuterCable> > cables_;
  std::map<const std::string, std::unique_ptr<const OuterDTC> > DTCs_;

  // negative cabling side
  std::map<const int, std::unique_ptr<OuterBundle> > negBundles_;
  std::map<const int, std::unique_ptr<OuterCable> > negCables_;
  std::map<const std::string, std::unique_ptr<const OuterDTC> > negDTCs_;
  std::map<const int, std::unique_ptr<OuterGBT> > gbts_;
  // All bundles, cables, and DTC are owned by the Cabling map, and the Cabling map only!!
};

}

#endif  // OUTERCABLINGMAP_HH

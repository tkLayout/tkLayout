#ifndef OUTERCABLINGMAP_HH
#define OUTERCABLINGMAP_HH

#include <global_constants.hh>
#include "global_funcs.hh"

#include "Property.hh"
#include "OuterCabling/ModulesToBundlesConnector.hh"
#include "OuterCabling/OuterDTC.hh"


/* Build the cabling map.
   The cabling map contains all the necessary bundles, cables and DTCs for a given tracker.
   All the relevant connections are made during the building process.
 */
class OuterCablingMap : public PropertyObject, public Buildable, public Identifiable<int> {
public:
  OuterCablingMap(Tracker* tracker);

  // positive cabling side
  const std::map<const int, OuterBundle*> getBundles() const { return getBundles(bundles_); }
  const std::map<const int, OuterCable*> getCables() const { return getCables(cables_); }
  const std::map<const std::string, const OuterDTC*> getDTCs() const { return getDTCs(DTCs_); }

  // negative cabling side
  const std::map<const int, OuterBundle*> getNegBundles() const { return getNegBundles(negBundles_); }
  const std::map<const int, OuterCable*> getNegCables() const { return getNegCables(negCables_); }
  const std::map<const std::string, const OuterDTC*> getNegDTCs() const { return getNegDTCs(negDTCs_); }

  

private:
  // positive cabling side
  std::map<const int, OuterBundle*> getBundles(const std::map<const int, std::unique_ptr<OuterBundle> >& bundles) const { 
    std::map<const int, OuterBundle*> myBundles; 
    for (const auto& bundleIt : bundles) {
      myBundles.insert(std::make_pair(bundleIt.first, bundleIt.second.get()));
    }
    return myBundles; 
  }
  std::map<const int, OuterCable*> getCables(const std::map<const int, std::unique_ptr<OuterCable> >& cables) const { 
    std::map<const int, OuterCable*> myCables; 
    for (const auto& cableIt : cables) {
      myCables.insert(std::make_pair(cableIt.first, cableIt.second.get()));
    }
    return myCables; 
  }
  std::map<const std::string, const OuterDTC*> getDTCs(const std::map<const std::string, std::unique_ptr<const OuterDTC> >& DTCs) const { 
    std::map<const std::string, const OuterDTC*> myDTCs; 
    for (const auto& DTCIt : DTCs) {
      myDTCs.insert(std::make_pair(DTCIt.first, DTCIt.second.get()));
    }
    return myDTCs; 
  }

  // negative cabling side
  std::map<const int, OuterBundle*> getNegBundles(const std::map<const int, std::unique_ptr<OuterBundle> >& negBundles) const {
    std::map<const int, OuterBundle*> myNegBundles; 
    for (const auto& bundleIt : negBundles) {
      myNegBundles.insert(std::make_pair(bundleIt.first, bundleIt.second.get()));
    }
    return myNegBundles; 
  }
  std::map<const int, OuterCable*> getNegCables(const std::map<const int, std::unique_ptr<OuterCable> >& negCables) const { 
    std::map<const int, OuterCable*> myNegCables; 
    for (const auto& cableIt : negCables) {
      myNegCables.insert(std::make_pair(cableIt.first, cableIt.second.get()));
    }
    return myNegCables;
  }
  std::map<const std::string, const OuterDTC*> getNegDTCs(const std::map<const std::string, std::unique_ptr<const OuterDTC> >& negDTCs) const { 
    std::map<const std::string, const OuterDTC*> myNegDTCs; 
    for (const auto& DTCIt : negDTCs) {
      myNegDTCs.insert(std::make_pair(DTCIt.first, DTCIt.second.get()));
    }
    return myNegDTCs;
  }

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
};


#endif  // OUTERCABLINGMAP_HH

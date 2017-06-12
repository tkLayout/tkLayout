#ifndef CABLINGMAP_HH
#define CABLINGMAP_HH

#include <global_constants.hh>
#include "global_funcs.hh"







/*#include <vector>
#include <string>
#include <memory>
#include <set>
#include <limits.h>

#include <boost/ptr_container/ptr_vector.hpp>

#include <TCanvas.h>

#include "global_funcs.hh"
#include "Property.hh"
#include "Barrel.hh"
#include "Endcap.hh"
#include "SupportStructure.hh"
#include "Visitor.hh"
#include "Visitable.hh"
#include "Cabling/CablingMap.hh"
#include "Cabling/DTC.hh"
#include "Cabling/ModulesToBundlesConnector.hh"

#include "MainConfigHandler.hh"
#include "DetIdBuilder.hh"*/



















/*#include <vector>
#include <string>
#include <memory>
#include <set>
#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.hh"*/
#include "Property.hh"
//#include <Tracker.hh>
#include "Cabling/DTC.hh"
#include "Cabling/ModulesToBundlesConnector.hh"
//#include "Visitor.hh"
//#include "Visitable.hh"
//#include "Cabling/DTC.hh"






class CablingMap : public PropertyObject, public Buildable, public Identifiable<int> {
public:
  CablingMap(Tracker* tracker);
  //CablingMap(const CablingMap& otherMap);
  //~CablingMap();

  const std::map<int, Bundle*>& getBundles() const { return bundles_; }
  const std::map<int, Cable*>& getCables() const { return cables_; }
  const std::map<const std::string, const DTC*>& getDTCs() const { return DTCs_; }

  const std::map<int, Bundle*>& getNegBundles() const { return negBundles_; }
  const std::map<int, Cable*>& getNegCables() const { return negCables_; }
  const std::map<const std::string, const DTC*>& getNegDTCs() const { return negDTCs_; }

private:
  void connectModulesToBundles(Tracker* tracker);
  void connectBundlesToCables(std::map<int, Bundle*>& bundles, std::map<int, Cable*>& cables, std::map<const std::string, const DTC*>& DTCs);
  void checkBundlesToCablesCabling(std::map<int, Cable*>& cables);

  std::map<int, Bundle*> bundles_;
  std::map<int, Cable*> cables_;
  std::map<const std::string, const DTC*> DTCs_;

  std::map<int, Bundle*> negBundles_;
  std::map<int, Cable*> negCables_;
  std::map<const std::string, const DTC*> negDTCs_;
};












#endif  // CABLINGMAP_HH

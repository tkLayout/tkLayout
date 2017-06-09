#ifndef MODULESTOBUNDLESCONNECTOR_HH
#define MODULESTOBUNDLESCONNECTOR_HH

#include <global_constants.hh>
#include "global_funcs.hh"

#include "Cabling/DTC.hh"

//#include "Visitor.hh"
//#include "Visitable.hh"


//#include "Cabling/CablingMap"
//#include <Tracker.hh>




class ModulesToBundlesConnector : public GeometryVisitor {
  
public:
  std::map<int, Bundle*> getBundles() { return bundles_; }
  std::map<int, Bundle*> getNegBundles() { return negBundles_; }

  void visit(Barrel& b);
  void visit(Layer& l);
  void visit(RodPair& r);     
  void visit(BarrelModule& m);

  void visit(Endcap& e);
  void visit(Disk& d);
  void visit(Ring& r);
  void visit(EndcapModule& m);

  void postVisit();

private:
  int computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth) const;
  void staggerModules(std::map<int, Bundle*>& bundles);
  void checkModulesToBundlesCabling(std::map<int, Bundle*>& bundles);

  std::string barrelName;
  int layerNumber;
  int numRods;
  int totalNumFlatRings;   

  std::string endcapName;
  int diskNumber;
  int ringNumber;
  int numModulesInRing;

  std::string type;
  int typeIndex;
  bool side;
   
  int phiSegmentRef;
  int negPhiSegmentRef;
  double phiRegionWidth;
  const double phiSectorWidth = 40. * M_PI / 180.;
  
  int bundleId;
  int bundleFlatId;   
  int bundleFlatIdB;      
  int bundleTiltedId;

  int negBundleId;
  int negBundleFlatId;
  int negBundleFlatIdB;
  int negBundleTiltedId;

  Bundle* bundle = NULL;
  Bundle* bundleFlat = NULL;
  Bundle* bundleFlatB = NULL;
  Bundle* bundleTilted = NULL;

  Bundle* negBundle = NULL;
  Bundle* negBundleFlat = NULL;
  Bundle* negBundleFlatB = NULL;
  Bundle* negBundleTilted = NULL;

  std::map<int, Bundle*> bundles_;
  std::map<int, Bundle*> negBundles_;

  
};


#endif  // MODULESTOBUNDLESCONNECTOR_HH

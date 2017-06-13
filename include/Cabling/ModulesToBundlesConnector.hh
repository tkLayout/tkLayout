#ifndef MODULESTOBUNDLESCONNECTOR_HH
#define MODULESTOBUNDLESCONNECTOR_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "Cabling/DTC.hh"


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

  std::string barrelName_;
  int layerNumber_;
  int numRods_;
  int totalNumFlatRings_;   

  std::string endcapName_;
  int diskNumber_;
  int ringNumber_;
  int numModulesInRing_;

  std::string type_;
  int typeIndex_;
  bool side_;
   
  int phiSegmentRef_;
  int negPhiSegmentRef_;
  double phiRegionWidth_;
  const double phiSectorWidth_ = 40. * M_PI / 180.;

  const int maxNumModulesPerBundle_ = 12;
  
  int bundleId_;
  int bundleFlatId_;   
  int bundleFlatIdB_;      
  int bundleTiltedId_;

  int negBundleId_;
  int negBundleFlatId_;
  int negBundleFlatIdB_;
  int negBundleTiltedId_;

  Bundle* bundle_ = NULL;
  Bundle* bundleFlat_ = NULL;
  Bundle* bundleFlatB_ = NULL;
  Bundle* bundleTilted_ = NULL;

  Bundle* negBundle_ = NULL;
  Bundle* negBundleFlat_ = NULL;
  Bundle* negBundleFlatB_ = NULL;
  Bundle* negBundleTilted_ = NULL;

  std::map<int, Bundle*> bundles_;
  std::map<int, Bundle*> negBundles_;
};


#endif  // MODULESTOBUNDLESCONNECTOR_HH

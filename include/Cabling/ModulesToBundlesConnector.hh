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
  double computePhiSegmentStart(const double phi, const double phiSegmentWidth, const bool isPositiveCablingSide) const;
  int computePhiSegmentRef(const double phi, const double phiSegmentStart, const double phiSegmentWidth, const bool isPositiveCablingSide) const;
  int computePhiSliceRef(const double phi, const double phiSliceStart, const double phiSliceWidth, const bool isPositiveCablingSide) const;
  std::string computeBundleType(const bool isBarrel, const std::string subDetectorName, const int layerDiskNumber, const int ringNumber = 0) const;
  int computeBundleTypeIndex(const bool isBarrel, const std::string bundleType, const int totalNumFlatRings = 0, const int maxNumModulesPerBundle = 0, const bool isTilted = false, const bool isExtraFlatPart = false) const;
  int computeBundleId(const bool isBarrel, const bool isPositiveCablingSide, const int layerDiskNumber, const int phiRef, const int bundleTypeIndex) const;

  void createAndStoreBundle(std::map<int, Bundle*>& bundles, std::map<int, Bundle*>& negBundles, const int bundleId, const std::string bundleType, const std::string subDetectorName, const int layerDiskNumber, const double phiSegmentWidth, const int phiSegmentRef, const double phiRegionStart, const double phiRegionWidth, int phiRegionRef, const double phiSectorWidth, const int phiSectorRef, const bool isPositiveCablingSide, const bool isTiltedPart = false);
  //void connectModuleToBundle(Module& m, Bundle* bundle);

  void staggerModules(std::map<int, Bundle*>& bundles);
  void checkModulesToBundlesCabling(const std::map<int, Bundle*>& bundles) const;

  bool isBarrel_;
  std::string barrelName_;
  int layerNumber_;
  int numRods_;
  int totalNumFlatRings_;
  double rodPhi_; 

  std::string endcapName_;
  int diskNumber_;
  int ringNumber_;
  int numModulesInRing_;

  std::string bundleType_;
  int bundleTypeIndex_;
  bool side_;
   
  int phiSegmentRef_;
  int negPhiSegmentRef_;
  double phiRegionWidth_;
  const double phiSectorWidth_ = 40. * M_PI / 180.;

  const int maxNumModulesPerBundle_ = 12;
  
  int bundleFlatId_;   
  int bundleExtraFlatId_;      
  int bundleTiltedId_;

  int negBundleFlatId_;
  int negBundleExtraFlatId_;
  int negBundleTiltedId_;

  std::map<int, Bundle*> bundles_;
  std::map<int, Bundle*> negBundles_;
};


#endif  // MODULESTOBUNDLESCONNECTOR_HH

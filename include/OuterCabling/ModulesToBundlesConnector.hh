#ifndef MODULESTOBUNDLESCONNECTOR_HH
#define MODULESTOBUNDLESCONNECTOR_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "OuterCabling/PhiPosition.hh"
#include "OuterCabling/OuterDTC.hh"
#include "OuterCabling/OuterGBT.hh"


/* This class is used to CONNECT MODULES TO BUNDLES.
 * The idea is to compute the Phi position references asociated to each module.
 * Then, a bundle corresponding to that position is built or found.
 * Lastly, the module is connected to that bundle, and vice-versa.
 * At the end of the process, 2 maps containing all the bundles which have been built are returned.
 * As each module has its own lpGBT, we also create module<->lpGBT links here.
*/
class ModulesToBundlesConnector : public GeometryVisitor {
public:
  std::map<int, OuterBundle*> getBundles() { return bundles_; }        // positive cabling side
  std::map<int, OuterBundle*> getNegBundles() { return negBundles_; }  // negative cabling side
  std::map<int, OuterGBT*> getGBTs() { return gbts_; }

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
  // BUILDING
  const bool computeBarrelFlatPartRodCablingSide(const double rodPhi, const double phiSegmentWidth) const;

  const Category computeBundleType(const bool isBarrel, const std::string subDetectorName, const int layerDiskNumber, const int ringNumber = 0) const;
  void buildBundle(DetectorModule& m, std::map<int, OuterBundle*>& bundles, std::map<int, OuterBundle*>& negBundles, const Category& bundleType, const bool isBarrel, const std::string subDetectorName, const int layerDiskNumber, const PhiPosition& modulePhiPosition, const bool isPositiveCablingSide, const int totalNumFlatRings = 0, const bool isTiltedPart = false, const bool isExtraFlatPart = false);
  void buildGBT(DetectorModule& m, std::map<int, OuterGBT*>& gbts);
  const int computeBundleTypeIndex(const bool isBarrel, const Category& bundleType, const int totalNumFlatRings = 0, const bool isTilted = false, const bool isExtraFlatPart = false) const;
  const int computeBundleId(const bool isBarrel, const bool isPositiveCablingSide, const int layerDiskNumber, const int phiRef, const int bundleTypeIndex) const;
  const int computeStereoBundleId(const bool isBarrel, const bool isPositiveCablingSide, const int layerDiskNumber, const int phiRef, const int bundleTypeIndex) const;
  OuterBundle* createAndStoreBundle(std::map<int, OuterBundle*>& bundles, std::map<int, OuterBundle*>& negBundles, const int bundleId, const int stereoBundleId, const Category& bundleType, const std::string subDetectorName, const int layerDiskNumber, const PhiPosition& modulePhiPosition, const bool isPositiveCablingSide, const bool isTiltedPart = false);
  OuterGBT* createAndStoreGBT(std::map<int, OuterGBT*>&gbts, const int gbtId);
  void connectModuleToBundle(DetectorModule& m, OuterBundle* bundle) const;
  void connectModuleToGBT(DetectorModule& m, OuterGBT* bundle) const;

  // STAGERRING
  void staggerModules(std::map<int, OuterBundle*>& bundles);

  // CHECKING
  void checkModulesToBundlesCabling(const std::map<int, OuterBundle*>& bundles) const;

  // ENDCAPS ONLY: ASSIGN THE MFB FANOUTS BRANCHES TO THE MODULES
  void connectEndcapModulesToBundlesFanoutBranches(std::map<int, OuterBundle*>& bundles);
  void checkEndcapModulesToBundlesFanoutBranchesCabling(const std::map<int, OuterBundle*>& bundles) const;

  std::map<int, OuterBundle*> bundles_;     // positive cabling side bundles.
  std::map<int, OuterBundle*> negBundles_;  // negative cabling side bundles.
  std::map<int, OuterGBT*> gbts_;

  bool isBarrel_;
  std::string barrelName_;
  int layerNumber_;
  int numRods_;
  int totalNumFlatRings_;              // Total number of flat rings on both (+Z) side and (-Z) side
  double rodPhi_; 

  std::string endcapName_;
  int diskNumber_;
  int ringNumber_;
  int numModulesInRing_;

  Category bundleType_;
  bool side_;
};


#endif  // MODULESTOBUNDLESCONNECTOR_HH

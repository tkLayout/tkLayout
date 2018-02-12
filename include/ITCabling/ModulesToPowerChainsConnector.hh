#ifndef MODULESTOPOWERCHAINSCONNECTOR_HH
#define MODULESTOPOWERCHAINSCONNECTOR_HH

#include <global_constants.hh>
#include "global_funcs.hh"
//#include "ITCabling/PhiPosition.hh"
#include "ITCabling/DTC.hh"


/* This class is used to CONNECT MODULES TO BUNDLES.
 * The idea is to compute the Phi position references asociated to each module.
 * Then, a bundle corresponding to that position is built or found.
 * Lastly, the module is connected to that bundle, and vice-versa.
 * At the end of the process, 2 maps containing all the bundles which have been built are returned.
*/
class ModulesToPowerChainsConnector : public GeometryVisitor {
public:
  std::map<int, Bundle*> getPowerChains() { return bundles_; }        // positive cabling side
  std::map<int, Bundle*> getNegPowerChains() { return negPowerChains_; }  // negative cabling side

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
  void buildBundle(DetectorModule& m, std::map<int, Bundle*>& bundles, std::map<int, Bundle*>& negPowerChains, const Category& bundleType, const bool isBarrel, const std::string subDetectorName, const int layerDiskNumber, const PhiPosition& modulePhiPosition, const bool isPositiveCablingSide, const int totalNumFlatRings = 0, const bool isTiltedPart = false, const bool isExtraFlatPart = false);
  const int computeBundleTypeIndex(const bool isBarrel, const Category& bundleType, const int totalNumFlatRings = 0, const bool isTilted = false, const bool isExtraFlatPart = false) const;
  const int computeBundleId(const bool isBarrel, const bool isPositiveCablingSide, const int layerDiskNumber, const int phiRef, const int bundleTypeIndex) const;
  const int computeStereoBundleId(const bool isBarrel, const bool isPositiveCablingSide, const int layerDiskNumber, const int phiRef, const int bundleTypeIndex) const;
  Bundle* createAndStoreBundle(std::map<int, Bundle*>& bundles, std::map<int, Bundle*>& negPowerChains, const int bundleId, const int stereoBundleId, const Category& bundleType, const std::string subDetectorName, const int layerDiskNumber, const PhiPosition& modulePhiPosition, const bool isPositiveCablingSide, const bool isTiltedPart = false);
  void connectModuleToBundle(DetectorModule& m, Bundle* bundle) const;

  // STAGERRING
  void staggerModules(std::map<int, Bundle*>& bundles);

  // CHECKING
  void checkModulesToPowerChainsCabling(const std::map<int, Bundle*>& bundles) const;

  std::map<int, Bundle*> bundles_;     // positive cabling side bundles.
  std::map<int, Bundle*> negPowerChains_;  // negative cabling side bundles.

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

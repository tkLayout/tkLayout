#ifndef MODULESTOPOWERCHAINSCONNECTOR_HH
#define MODULESTOPOWERCHAINSCONNECTOR_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "ITCabling/inner_cabling_functions.hh"
//#include "ITCabling/InnerPhiPosition.hh"
#include "ITCabling/HvLine.hh"


/* This class is used to CONNECT MODULES TO BUNDLES.
 * The idea is to compute the Phi position references asociated to each module.
 * Then, a bundle corresponding to that position is built or found.
 * Lastly, the module is connected to that bundle, and vice-versa.
 * At the end of the process, 2 maps containing all the bundles which have been built are returned.
*/
class ModulesToPowerChainsConnector : public GeometryVisitor {
public:
  std::map<int, PowerChain*> getPowerChains() { return powerChains_; }
  //std::map<int, Bundle*> getNegPowerChains() { return negPowerChains_; }  // negative cabling side

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
  const bool computeXSide(const double modCenterX) const;
  const bool computeBarrelModuleZEnd(const int side, const int ring, const double rodPhi, const int numRods, const bool isPositiveXSide) const;
  const bool computeBarrelCentralModuleZEnd(const double rodPhi, const int numRods, const bool isPositiveXSide) const;
  const int computeForwardModulePhiPowerChain(const double modPhi, const int numModulesInRing, const bool isPositiveZEnd) const;

  //const Category computeBundleType(const bool isBarrel, const std::string subDetectorName, const int layerDiskNumber, const int ringNumber = 0) const;
  void buildPowerChain(DetectorModule& m, std::map<int, PowerChain*>& powerChains, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex = 0);
  //const int computeBundleTypeIndex(const bool isBarrel, const Category& bundleType, const int totalNumFlatRings = 0, const bool isTilted = false, const bool isExtraFlatPart = false) const;
  const int computePowerChainId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex) const;
  //const int computeStereoBundleId(const bool isBarrel, const bool isPositiveCablingSide, const int layerDiskNumber, const int phiRef, const int bundleTypeIndex) const;
  PowerChain* createAndStorePowerChain(std::map<int, PowerChain*>& powerChains, const int powerChainId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex);
  void connectModuleToPowerChain(DetectorModule& m, PowerChain* powerChain) const;

  // STAGERRING
  //void staggerModules(std::map<int, Bundle*>& bundles);

  // CHECKING
  void checkModulesToPowerChainsCabling(const std::map<int, PowerChain*>& powerChains) const;

  std::map<int, PowerChain*> powerChains_;
  //std::map<int, Bundle*> negPowerChains_;  // negative cabling side bundles.

  //bool isBarrel_;
  std::string barrelName_;
  int layerNumber_;
  int numRods_;
  double rodPhi_; 

  std::string endcapName_;
  bool endcapEnd_;
  int diskNumber_;
  int ringNumber_;
  int numModulesInRing_;

  //Category bundleType_;
  //bool isPositiveZEnd_;
  //bool isPositiveXSide_;
};


#endif  // MODULESTOBUNDLESCONNECTOR_HH

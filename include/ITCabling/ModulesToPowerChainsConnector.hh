#ifndef MODULESTOPOWERCHAINSCONNECTOR_HH
#define MODULESTOPOWERCHAINSCONNECTOR_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "ITCabling/inner_cabling_functions.hh"
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
  const bool computeBarrelModuleZEnd(const int side, const int ring, const int layerNumber) const;
  const bool computeBarrelCentralModuleZEnd(const int layerNumber) const;
  const std::pair<int, int> computeForwardModulePhiPowerChain(const double modPhi, const int numModulesInRing, const bool isPositiveZEnd) const;

  void buildPowerChain(DetectorModule& m, std::map<int, PowerChain*>& powerChains, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex = 0);
  const int computePowerChainId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex) const;
  PowerChain* createAndStorePowerChain(std::map<int, PowerChain*>& powerChains, const int powerChainId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringQuarterIndex);
  void connectModuleToPowerChain(DetectorModule& m, PowerChain* powerChain) const;

  // CHECKING
  void checkModulesToPowerChainsCabling(const std::map<int, PowerChain*>& powerChains) const;

  std::map<int, PowerChain*> powerChains_;

  std::string barrelName_;
  int layerNumber_;
  int numRods_;
  double rodPhi_; 

  std::string endcapName_;
  bool endcapEnd_;
  int diskNumber_;
  int ringNumber_;
  int numModulesInRing_;
};


#endif  // MODULESTOPOWERCHAINSCONNECTOR_HH

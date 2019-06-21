#ifndef MODULESTOPOWERCHAINSCONNECTOR_HH
#define MODULESTOPOWERCHAINSCONNECTOR_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "InnerCabling/inner_cabling_functions.hh"
#include "InnerCabling/HvLine.hh"


/* This class is used to CONNECT MODULES TO THE SERIAL POWER CHAINS.
 * The idea is to use generic information on the module location, and use it to assign the module to a given power chain.
 * Once the mapping is done, it is cross-checked that power chains are not assigned more modules than what they can afford.
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
  // COLLECTING MODULE INFO
  const bool computeXSide(const double modCenterX) const;
  const std::pair<bool, bool> computeBarrelModuleZEnd(const int side, const int ring, const int layerNumber) const;
  const bool computeBarrelCentralModuleZEnd(const int layerNumber) const;
  const std::pair<int, int> computeForwardModulePhiPowerChain(const double modPhi, const int numModulesInRing, const bool isPositiveZEnd) const;

  // BUILDING POWER CHAIN
  void buildPowerChain(DetectorModule& m, std::map<int, PowerChain*>& powerChains, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const bool isLongBarrel = false, const int halfRingIndex = 0);
  const int computePowerChainId(const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int halfRingIndex) const;
  PowerChain* createAndStorePowerChain(std::map<int, PowerChain*>& powerChains, const int powerChainId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const bool isLongBarrel, const int halfRingIndex);
  void connectModuleToPowerChain(DetectorModule& m, PowerChain* powerChain) const;

  // CHECKING POWER CHAIN
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

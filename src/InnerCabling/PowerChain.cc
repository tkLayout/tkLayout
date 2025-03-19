#include "InnerCabling/PowerChain.hh"
#include "InnerCabling/HvLine.hh"


PowerChain::PowerChain(const int powerChainId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const bool isLongBarrel, const int halfRingIndex, const bool isAtSmallerAbsZDeeInDoubleDisk, const bool isAtSmallerAbsZSideInDee, const bool isTEPXSpecialRing) :
  isPositiveZEnd_(isPositiveZEnd),
  isPositiveXSide_(isPositiveXSide),
  subDetectorName_(subDetectorName),
  layerDiskNumber_(layerDiskNumber),
  phiRef_(phiRef),
  isLongBarrel_(isLongBarrel),
  halfRingIndex_(halfRingIndex),
  isAtSmallerAbsZDeeInDoubleDisk_(isAtSmallerAbsZDeeInDoubleDisk),
  isAtSmallerAbsZSideInDee_(isAtSmallerAbsZSideInDee),
  isSplitOverRings_(isTEPXSpecialRing)
{
  myid(powerChainId);
  isBarrel_ = inner_cabling_functions::isBarrel(subDetectorName);
  ringNumber_ = inner_cabling_functions::computeRingNumber(halfRingIndex);
  isSmallerAbsZHalfRing_ = inner_cabling_functions::isSmallerAbsZHalfRing(halfRingIndex);

  plotColor_ = computePlotColor(isBarrel_, isPositiveZEnd, phiRef, halfRingIndex);

  // BUILD HVLINE, TO WHICH THE MODULES OF THE POWER CHAIN ARE ALL CONNECTED
  buildHvLine(powerChainId);
};


/*
 *  Assign a module to the power chain.
 */
void PowerChain::addModule(Module* m) { 
  modules_.push_back(m);
  hvLine_->addModule(m);
  m->setHvLine(hvLine_.get());
}


/*
 * Returns whether a power chain has 4 Ampere or 8 Ampere.
 * NB 1: Assumed all modules conected to the same power chain are of the same type.
 * NB 2: THIS DEPENDS ON THE MODULE TYPE: IT ONLY MAKES SENSE TO CALL THIS AFTER THE MODULES HAVE BEEN ASSIGNED TO THE POWER CHAIN.
 */
const PowerChainType PowerChain::powerChainType() const {
  if (modules().size() == 0) {
    logERROR("Tried to call PowerChain::powerChainType() on a power chain connected to 0 module.");
    return PowerChainType::IUNDEFINED;
  }
  else {
    const int numROCsPerModule = modules().front()->outerSensor().totalROCs();
    if (numROCsPerModule == 2) return PowerChainType::I4A;
    else if (numROCsPerModule == 4) return PowerChainType::I8A;
    else {
      logERROR(any2str("Found ")
	       + any2str(numROCsPerModule)
	       + any2str(" ROCs per module, which is not supported. If this is intended, tune PowerChain::powerChainType().")
	       );
      return PowerChainType::IUNDEFINED;
    }
  }
}


/*
 * Compute power chain color on website.
 * Power chains next to each other in space, must be of different colors.
 */
const int PowerChain::computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int phiRef, const int halfRingIndex) const {
  int plotColor = 0;

  const int plotPhi = femod(phiRef, 2);

  if (isBarrel) {
    const int plotZEnd = (isPositiveZEnd ? 0 : 1);
    plotColor = plotZEnd * 2 + plotPhi + 6;
  }
  else {
    const int plotRingQuarter = femod(halfRingIndex, 6);
    plotColor = plotRingQuarter * 2 + plotPhi + 1;
  }

  return plotColor;
}


/* Build HvLine asociated to the power chain.
 */
void PowerChain::buildHvLine(const int powerChainId) {
  std::string hvLineName = computeHvLineName(powerChainId);
  std::unique_ptr<HvLine> hvLine(new HvLine(hvLineName));
  hvLine->setPowerChain(this);
  hvLine_ = std::move(hvLine);
}


/* Compute HvLine name.
 */
const std::string PowerChain::computeHvLineName(const int powerChainId) const {
  std::ostringstream hvLineNameStream;
  hvLineNameStream << powerChainId << "_HV";
  const std::string hvLineName = hvLineNameStream.str();
  return hvLineName;
}

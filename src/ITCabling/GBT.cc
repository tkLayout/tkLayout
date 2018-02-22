#include "ITCabling/GBT.hh"
#include "ITCabling/InnerBundle.hh"


GBT::GBT(PowerChain* myPowerChain, const std::string GBTId, const int myGBTPhiIndex, const int numELinksPerModule) :
  myGBTId_(GBTId),
  myGBTPhiIndex_(myGBTPhiIndex),
  numELinksPerModule_(numELinksPerModule)
{
  myPowerChain_ = myPowerChain;

  //myPowerChain->setGBT(this);
  //plotColor_ = computePlotColor(isBarrel_, isPositiveZEnd, phiRef, ringQuarterIndex);
};


GBT::~GBT() {
  delete myPowerChain_;    // TO DO: switch to smart pointers and remove this!
  myPowerChain_ = nullptr;

  delete myBundle_;
  myBundle_ = nullptr;
}


void GBT::addModule(Module* m) { 
  modules_.push_back(m);
}



/*const int GBT::computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int phiRef, const int ringQuarterIndex) const {
  int plotColor = 0;

  const int plotPhi = femod(phiRef, 2);

  if (isBarrel) {
    const int plotZEnd = (isPositiveZEnd ? 0 : 1);
    plotColor = plotZEnd * 2 + plotPhi + 6;
  }
  else {
    const int plotRingQuarter = femod(ringQuarterIndex, 6);
    plotColor = plotRingQuarter * 2 + plotPhi + 1;
  }

  return plotColor;
  }*/

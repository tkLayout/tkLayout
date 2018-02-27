#include "ITCabling/GBT.hh"
#include "ITCabling/InnerBundle.hh"


GBT::GBT(PowerChain* myPowerChain, const std::string GBTId, const int myGBTPhiIndex, const int numELinksPerModule) :
  myGBTId_(GBTId),
  myGBTPhiIndex_(myGBTPhiIndex),
  numELinksPerModule_(numELinksPerModule)
{
  myPowerChain_ = myPowerChain;

  //myPowerChain->setGBT(this);
  plotColor_ = computePlotColor(myPowerChain);
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



const int GBT::computePlotColor(const PowerChain* myPowerChain) const {
  const int plotColor = myPowerChain->plotColor();
  return plotColor;
}

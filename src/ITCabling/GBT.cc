#include "ITCabling/GBT.hh"
#include "ITCabling/InnerBundle.hh"


GBT::GBT(PowerChain* myPowerChain, const std::string GBTId, const int myGBTIndex, const int myGBTIndexColor, const int numELinksPerModule) :
  myGBTId_(GBTId),
  myGBTIndex_(myGBTIndex),
  myGBTIndexColor_(myGBTIndexColor),
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

#include "InnerCabling/GBT.hh"
#include "InnerCabling/InnerBundle.hh"


GBT::GBT(PowerChain* myPowerChain, const std::string GBTId, const int myGBTIndexInPowerChain) :
  myGBTId_(GBTId),
  myGBTCMSSWId_(0), // Need to have consecutive integers, hence is done after full cabling map is created.
  myGBTIndexInPowerChain_(myGBTIndexInPowerChain),
  plotStyleGBTIndexInPowerChain_(computePlotStyleGBTIndexInPowerChain(myGBTIndexInPowerChain, myPowerChain))
{
  myPowerChain_ = myPowerChain;
  plotPowerChainColor_ = computePlotColor(myPowerChain);
};


/*
 *  Assign a module to the GBT.
 */
void GBT::addModule(Module* m, const int numELinks) {
  if (m == nullptr) throw PathfulException("GBT::addModule Error: Module pointer is nullptr");
  m->setGBT(this);
  m->setNumELinks(numELinks);
  modules_.push_back(m);
}


/*
 * Compute GBT plot style on website.
 * To distinguish different GBTs within the same power chain, alternation of fill, contour, dashed, and crosshatched styles is used.
 */
const int GBT::computePlotStyleGBTIndexInPowerChain(const int myGBTIndexInPowerChain, PowerChain* myPowerChain) const {

  const int myGBTIndexInPowerChainPlotStyle = femod(myGBTIndexInPowerChain, 4);
 
  return myGBTIndexInPowerChainPlotStyle;
}


/*
 * Compute GBT color on website.
 * GBT color is the same as the power chain its modules belong to.
 */
const int GBT::computePlotColor(const PowerChain* myPowerChain) const {
  const int plotPowerChainColor = myPowerChain->plotColor();
  return plotPowerChainColor;
}

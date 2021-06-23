#include "InnerCabling/GBT.hh"
#include "InnerCabling/InnerBundle.hh"


GBT::GBT(PowerChain* myPowerChain, const std::string GBTId, const int myGBTIndexInPowerChain, const int numELinksPerModule) :
  myGBTId_(GBTId),
  myGBTCMSSWId_(0), // Need to have consecutive integers, hence is done after full cabling map is created.
  myGBTIndexInPowerChain_(myGBTIndexInPowerChain),
  plotStyleGBTIndexInPowerChain_(computePlotStyleGBTIndexInPowerChain(myGBTIndexInPowerChain, myPowerChain)),
  numELinksPerModule_(numELinksPerModule)
{
  myPowerChain_ = myPowerChain;
  plotPowerChainColor_ = computePlotColor(myPowerChain);
};


/*
 *  Assign a module to the GBT.
 */
void GBT::addModule(Module* m) { 
  modules_.push_back(m);
}


/*
 * Compute GBT plot style on website.
 * To distinguish different GBTs within the same power chain, alternation of fill, contour and dashed styles is used.
 */
const int GBT::computePlotStyleGBTIndexInPowerChain(const int myGBTIndexInPowerChain, PowerChain* myPowerChain) const {
  const bool isBarrel = myPowerChain->isBarrel();
  const int numGBTsInPowerChain = myPowerChain->numGBTsInPowerChain();

  // Usually, only 2 plot styles are needed (for example, alternation of fill and empty plot styles) 
  // to distinguish the GBTs among a power chain.
  // Exception: case of an odd number of GBTs within the same power chain in the barrel:
  // a third plot style (for example, dashed) becomes necessary.
  const int myGBTIndexInPowerChainPlotStyle = ( (!isBarrel || femod(numGBTsInPowerChain, 2) == 0 ) 
                                                // 2 styles are enough (for example, full and contour)
                                                ? femod(myGBTIndexInPowerChain, 2)
                                                // 3 styles are necessary (for example, full, contour, and dashed)
                                                : femod(myGBTIndexInPowerChain, 3));
 
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

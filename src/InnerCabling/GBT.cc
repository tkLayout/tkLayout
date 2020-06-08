#include "InnerCabling/GBT.hh"
#include "InnerCabling/InnerBundle.hh"


GBT::GBT(PowerChain* myPowerChain, const std::string GBTId, const int myGBTIndex, const int myGBTIndexColor, const int numELinksPerModule) :
  myGBTId_(GBTId),
  myGBTCMSSWId_(0), // Need to have consecutive integers, hence is done after full cabling map is created.
  myGBTIndex_(myGBTIndex),
  myGBTIndexColor_(myGBTIndexColor),
  numELinksPerModule_(numELinksPerModule)
{
  myPowerChain_ = myPowerChain;
  plotColor_ = computePlotColor(myPowerChain);
};


/*
 *  Assign a module to the GBT.
 */
void GBT::addModule(Module* m) { 
  modules_.push_back(m);
}


/*
 * Compute GBT color on website.
 * GBT color is the same as the power chain its modules belong to.
 * To distinguish different GBTs on the same power chain, alternation of fill and contour style is used.
 */
const int GBT::computePlotColor(const PowerChain* myPowerChain) const {
  const int plotColor = myPowerChain->plotColor();
  return plotColor;
}

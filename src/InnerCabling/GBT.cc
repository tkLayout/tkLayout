#include "InnerCabling/GBT.hh"
#include "InnerCabling/InnerBundle.hh"


GBT::GBT(PowerChain* myPowerChain, const std::string GBTId, const int myGBTIndexInPowerChain, const int numELinksPerModule) :
  myGBTId_(GBTId),
  myGBTCMSSWId_(0), // Need to have consecutive integers, hence is done after full cabling map is created.
  myGBTIndexInPowerChain_(myGBTIndexInPowerChain),
  plotGBTIndexInPowerChain_(computePlotGBTIndexInPowerChain(myGBTIndexInPowerChain)),
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
 * To distinguish different GBTs within the same power chain, alternation of fill, contour and dashed styles is used
 */
const int GBT::computePlotGBTIndexInPowerChain(const int myGBTIndexInPowerChain) const {

  if (isBarrel) {
    std::cout << "ringRef = " << ringRef << std::endl;
    std::cout << "moduleRef = " << moduleRef << std::endl;
    std::cout << "myGBTIndexInPowerChainExact = " << myGBTIndexInPowerChainExact << std::endl;
    std::cout << "myGBTIndexInPowerChain = " << myGBTIndexInPowerChain << std::endl;
  }

  //int myGBTIndexInPowerChainPlotStyle = myGBTIndexInPowerChain;
  //if (isBarrel && phiRefInPowerChain == 1 && femod(numGBTsInPowerChain, 2) == 0) myGBTIndexInPowerChainPlotStyle += 1;

  //myGBTIndexInPowerChainPlotStyle = femod(myGBTIndexInPowerChainPlotStyle, 2);
  std::cout << "numGBTsInPowerChain = " << numGBTsInPowerChain << std::endl;
  int myGBTIndexInPowerChainPlotStyle = ( (!isBarrel || femod(numGBTsInPowerChain, 2) == 0 ) ? femod(myGBTIndexInPowerChain, 2) : femod(myGBTIndexInPowerChain, 3));
  //myGBTIndexInPowerChainPlotStyle = femod(myGBTIndexInPowerChainPlotStyle, 3);
  std::cout << "myGBTIndexInPowerChainPlotStyle = " << myGBTIndexInPowerChainPlotStyle << std::endl;

  //myGBTIndexInPowerChainPlotStyle = myGBTIndexInPowerChainPlotStyle;

}


/*
 * Compute GBT color on website.
 * GBT color is the same as the power chain its modules belong to.
 */
const int GBT::computePlotColor(const PowerChain* myPowerChain) const {
  const int plotPowerChainColor = myPowerChain->plotColor();
  return plotPowerChainColor;
}

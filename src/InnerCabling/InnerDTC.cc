#include "InnerCabling/InnerDTC.hh"


InnerDTC::InnerDTC(const int DTCId, const bool isPositiveZEnd, const bool isPositiveXSide) :
  isPositiveZEnd_(isPositiveZEnd),
  isPositiveXSide_(isPositiveXSide)
{
  myid(DTCId);

  plotColor_ = computePlotColor(DTCId);
};


/*
 *  Connect a Bundle to the DTC.
 */
void InnerDTC::addBundle(InnerBundle* bundle) { 
  bundles_.push_back(bundle);
}



/*
 * Compute DTC color on website.
 * DTCs next to each other in space, must be of different colors.
 */
const int InnerDTC::computePlotColor(const int DTCId) const {
  const int plotColor = DTCId % 10 + 4;

  return plotColor;
}

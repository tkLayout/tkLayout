#ifndef INNERDTC_HH
#define INNERDTC_HH

#include <vector>
#include <string>

#include "Property.hh"

#include "InnerCabling/InnerBundle.hh"
#include "InnerCabling/inner_cabling_functions.hh"


/*
 * Inner Tracker DTC class.
 * All Fiber Bundles connected to a given DTC can be accessed.
 * General info on the DTC is also provided.
 */
class InnerDTC : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef std::vector<InnerBundle*> Container; 

public:
  InnerDTC(const int DTCId, const bool isPositiveZEnd, const bool isPositiveXSide);

  // BUNDLES CONNECTED TO THE DTC
  const Container& bundles() const { return bundles_; }
  const int numBundles() const { return bundles_.size(); }
  void addBundle(InnerBundle* bundle);

  // GENERAL INFO
  void setCMSSWId(const int cmsswId) { myDTCCMSSWId_ = cmsswId; }
  const int getCMSSWId() const { return myDTCCMSSWId_; }

  const bool isPositiveZEnd() const { return isPositiveZEnd_; }
  const bool isPositiveXSide() const { return isPositiveXSide_; }
  const int plotColor() const { return plotColor_; }

private:
  const int computePlotColor(const int DTCId) const;

  Container bundles_;

  int myDTCCMSSWId_;

  bool isPositiveZEnd_;
  bool isPositiveXSide_;
  int plotColor_;
};



#endif

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


  const bool isPositiveZEnd() const { return isPositiveZEnd_; }
  const bool isPositiveXSide() const { return isPositiveXSide_; }

  // GENERAL INFO
  const int plotColor() const { return plotColor_; }

private:
  const int computePlotColor(const int DTCId) const;

  Container bundles_;

  bool isPositiveZEnd_;
  bool isPositiveXSide_;

  int plotColor_;
};



#endif

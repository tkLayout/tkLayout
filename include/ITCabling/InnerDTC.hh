#ifndef INNERDTC_HH
#define INNERDTC_HH

#include <vector>
#include <string>

#include "Property.hh"
//#include "Module.hh"
#include "ITCabling/InnerBundle.hh"
#include "ITCabling/inner_cabling_functions.hh"



class InnerDTC : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<InnerBundle> Container; 

public:
  InnerDTC(const int DTCId, const bool isPositiveZEnd, const bool isPositiveXSide);
  //~InnerDTC();

  // BUNDLES CONNECTED TO THE DTC
  const Container& bundles() const { return bundles_; }
  const int numBundles() const { return bundles_.size(); }
  void addBundle(InnerBundle* bundle);


  const bool isPositiveZEnd() const { return isPositiveZEnd_; }
  const bool isPositiveXSide() const { return isPositiveXSide_; }

  // GENERAL INFO
  //const InnerDTCType powerChainType() const { return powerChainType_; }
  const int plotColor() const { return plotColor_; }

private:
  const int computePlotColor(const int DTCId) const;

  Container bundles_;

  bool isPositiveZEnd_;
  bool isPositiveXSide_;

  //DTC* myDTC_ = nullptr;

  int plotColor_;
};



#endif

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
  InnerDTC(const int DTCId);
  //~InnerDTC();

  // BUNDLES CONNECTED TO THE DTC
  const Container& bundles() const { return bundles_; }
  const int numBundles() const { return bundles_.size(); }
  void addBundle(InnerBundle* bundle);


  // GENERAL INFO
  //const InnerDTCType powerChainType() const { return powerChainType_; }
  //const int plotColor() const { return plotColor_; }

private:
  //const int computePlotColor(const bool isBarrel, const bool isPositiveZEnd, const int phiRef, const int ringQuarterIndex) const;

  Container bundles_;

  //DTC* myDTC_ = nullptr;

  //int plotColor_;
};



#endif

#ifndef CABLE_HH
#define CABLE_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Bundle.hh"
//#include "CablingVisitable.h"

namespace insur { class DTC; }
using insur::DTC;


/*using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;*/

//class Cable : public PropertyObject, public Buildable, public Identifiable<int>, public CablingVisitable {
class Cable : public PropertyObject, public Buildable, public Identifiable<int> {
  double phiSectorWidth_;
  int phiSectorRef_;
  std::string type_;
  int slot_;

  int servicesChannel_;

  DTC* myDTC_ = NULL;

  typedef PtrVector<Bundle> Container;
  Container bundles_;

  int computeServicesChannel(int phiSectorRef, std::string type, int slot);

  //Property<int, Default> nBundlesPerCable;

public:
  Cable(int id, const double phiSectorWidth, int phiSectorRef, std::string type, int slot) ;

  const std::string type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int slot() const { return slot_; }

  const int servicesChannel() const { return servicesChannel_; }

  const DTC* getDTC() const { return myDTC_; }
  void setDTC(DTC* dtc) { myDTC_ = dtc; }


  void addBundle(Bundle* b) { bundles_.push_back(b); }
  const Container& bundles() const { return bundles_; }

  int numBundles() const { return bundles_.size(); }







  /*Cable() :
            nBundlesPerCable      ("nBundlesPerCable"      , parsedAndChecked(), 12)
  {}

  void setup() {}


  Container& bundles() { return bundles_; }
 
  int nBundles() const { return bundles_.size(); }
  int maxBundles() {return nBundlesPerCable(); }
   
  void check() override;
  void build();


  void accept(CablingVisitor& v) { 
    v.visit(*this); 
    for (auto& b : bundles_) { b.accept(v); }
  }
  void accept(ConstCablingVisitor& v) const { 
    v.visit(*this); 
    for (const auto& b :bundles_) { b.accept(v); }
    }*/

};



#endif

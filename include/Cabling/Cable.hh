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
  typedef PtrVector<Bundle> Container;
public:
  Cable(const int id, const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);
  ~Cable();

  const Container& bundles() const { return bundles_; }

  void addBundle(Bundle* b) { bundles_.push_back(b); }
  int numBundles() const { return bundles_.size(); }

  const DTC* getDTC() const { return myDTC_; }

  const Category& type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int slot() const { return slot_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }
  const int servicesChannel() const { return servicesChannel_; }


  /*Cable() :
            nBundlesPerCable      ("nBundlesPerCable"      , parsedAndChecked(), 12)
  {}

  void setup() {}
  Property<int, Default> nBundlesPerCable;

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

private:
  const int computeServicesChannel(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  void buildDTC(const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);
  const std::string computeDTCName(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  
  Container bundles_;

  DTC* myDTC_ = nullptr;

  double phiSectorWidth_;
  int phiSectorRef_;
  Category type_;
  int slot_;
  bool isPositiveCablingSide_;
  int servicesChannel_;
};



#endif

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

  //Property<int, Default> nBundlesPerCable;

public:

  Cable(int id, const double phiSectorWidth, int phiSectorRef, std::string type, int slot) {
    myid(id);
    phiSectorWidth_ = phiSectorWidth;
    phiSectorRef_ = phiSectorRef;
    type_ = type;
    slot_ = slot;

    // Assign servicesChannel (depends on cable type and phiSectorRef)
    servicesChannel_ = 0;
    if (type == "PS10G") {
      if (phiSectorRef == 0) servicesChannel_ = 1;
      else if (phiSectorRef == 1) servicesChannel_ = 3;
      else if (phiSectorRef == 2) servicesChannel_ = 4;
      else if (phiSectorRef == 3) servicesChannel_ = 5;
      else if (phiSectorRef == 4) servicesChannel_ = 6;
      else if (phiSectorRef == 5) servicesChannel_ = 7;
      else if (phiSectorRef == 6) servicesChannel_ = 9;
      else if (phiSectorRef == 7) servicesChannel_ = 10;
      else if (phiSectorRef == 8) servicesChannel_ = 12;   
    }
    else if (type == "PS5G") {
      if (slot != 3) {
	if (phiSectorRef == 0) servicesChannel_ = 1;
	else if (phiSectorRef == 1) servicesChannel_ = 3;
	else if (phiSectorRef == 2) servicesChannel_ = 4;
	else if (phiSectorRef == 3) servicesChannel_ = 6;
	else if (phiSectorRef == 4) servicesChannel_ = 7;
	else if (phiSectorRef == 5) servicesChannel_ = 8;
	else if (phiSectorRef == 6) servicesChannel_ = 10;
	else if (phiSectorRef == 7) servicesChannel_ = 11;
	else if (phiSectorRef == 8) servicesChannel_ = 12;
      }
      else {
	if (phiSectorRef == 0) servicesChannel_ = 1;
	else if (phiSectorRef == 1) servicesChannel_ = 3;
	else if (phiSectorRef == 2) servicesChannel_ = 4;
	else if (phiSectorRef == 3) servicesChannel_ = 5;
	else if (phiSectorRef == 4) servicesChannel_ = 6;
	else if (phiSectorRef == 5) servicesChannel_ = 7;
	else if (phiSectorRef == 6) servicesChannel_ = 9;
	else if (phiSectorRef == 7) servicesChannel_ = 10;
	else if (phiSectorRef == 8) servicesChannel_ = 12;
      }
    }
    else if (type == "2S") {
      if (slot != 1 && slot != 2) {
	if (phiSectorRef == 0) servicesChannel_ = 1;
	else if (phiSectorRef == 1) servicesChannel_ = 2;
	else if (phiSectorRef == 2) servicesChannel_ = 4;
	else if (phiSectorRef == 3) servicesChannel_ = 5;
	else if (phiSectorRef == 4) servicesChannel_ = 6;
	else if (phiSectorRef == 5) servicesChannel_ = 7;
	else if (phiSectorRef == 6) servicesChannel_ = 9;
	else if (phiSectorRef == 7) servicesChannel_ = 10;
	else if (phiSectorRef == 8) servicesChannel_ = 12;
      }
      if (slot == 1 || slot == 2) {
	if (phiSectorRef == 0) servicesChannel_ = 1;
	else if (phiSectorRef == 1) servicesChannel_ = 2;
	else if (phiSectorRef == 2) servicesChannel_ = 4;
	else if (phiSectorRef == 3) servicesChannel_ = 6;
	else if (phiSectorRef == 4) servicesChannel_ = 7;
	else if (phiSectorRef == 5) servicesChannel_ = 8;
	else if (phiSectorRef == 6) servicesChannel_ = 10;
	else if (phiSectorRef == 7) servicesChannel_ = 11;
	else if (phiSectorRef == 8) servicesChannel_ = 12;
      }
    }
  };


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

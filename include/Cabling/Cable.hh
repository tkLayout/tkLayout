#ifndef CABLE_HH
#define CABLE_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Ribbon.hh"
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
  int cableIndex_;

  DTC* myDTC_ = NULL;

  typedef PtrVector<Ribbon> Container;
  Container ribbons_;

  //Property<int, Default> nRibbonsPerCable;

public:

  Cable(int id, const double phiSectorWidth, int phiSectorRef, std::string type, int cableIndex) {
    myid(id);
    phiSectorWidth_ = phiSectorWidth;
    phiSectorRef_ = phiSectorRef;
    type_ = type;
    cableIndex_ = cableIndex;    
  };


  const std::string type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int cableIndex() const { return cableIndex_; }

  const DTC* getDTC() const { return myDTC_; }
  void setDTC(DTC* dtc) { myDTC_ = dtc; }


  void addRibbon(Ribbon* r) { ribbons_.push_back(r); }
  const Container& ribbons() const { return ribbons_; }

  int numRibbons() const { return ribbons_.size(); }







  /*Cable() :
            nRibbonsPerCable      ("nRibbonsPerCable"      , parsedAndChecked(), 12)
  {}

  void setup() {}


  Container& ribbons() { return ribbons_; }
 
  int nRibbons() const { return ribbons_.size(); }
  int maxRibbons() {return nRibbonsPerCable(); }
   
  void check() override;
  void build();


  void accept(CablingVisitor& v) { 
    v.visit(*this); 
    for (auto& r : ribbons_) { r.accept(v); }
  }
  void accept(ConstCablingVisitor& v) const { 
    v.visit(*this); 
    for (const auto& r :ribbons_) { r.accept(v); }
    }*/

};



#endif

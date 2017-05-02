#ifndef DTC_HH
#define DTC_HH

#include <vector>
#include <string>

#include "Cable.hh"
#include "Module.hh"


/*using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;*/

//class DTC : public PropertyObject, public Buildable, public Identifiable<int>, public CablingVisitable {
class DTC : public PropertyObject, public Buildable, public Identifiable<int> {
  std::string name_;

  double phiSectorWidth_;
  int phiSectorRef_;
  std::string type_;
  int cableIndex_;

  int plotColor_;


  typedef PtrVector<Cable> Container;
  Container cable_;


  //typedef PtrVector<Cable> Container;
  //Container cables_;

  //Property<int, Default> nCablesPerDTC;

public:

  DTC(std::string name, const double phiSectorWidth, int phiSectorRef, std::string type, int cableIndex) {
    name_ = name;
    phiSectorWidth_ = phiSectorWidth;
    phiSectorRef_ = phiSectorRef;
    type_ = type;
    cableIndex_ = cableIndex;

    if (type == "PS10G") plotColor_ = cableIndex;
    else if (type == "PS5G") plotColor_ = 1 + cableIndex;
    else if (type == "2S") plotColor_ = 6 + cableIndex;  
  };


  const std::string name() const { return name_; }
  const std::string type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int cableIndex() const { return cableIndex_; }


  void addCable(Cable* c) { cable_.push_back(c); }
  const Container& cable() const { return cable_; }

  const int plotColor() const { return plotColor_; }



  /*DTC() :
            nCablesPerDTC      ("nCablesPerDTC"      , parsedAndChecked(), 1)
  {}

  void setup() {
  }

  Container& cables() { return cables_; }
  const Container& cables() const { return cables_; }
  int nCables() const { return cables_.size(); }
  int maxCables() {return nCablesPerDTC(); }
  
  void check() override;
  void build();


  void accept(CablingVisitor& v) { 
    v.visit(*this); 
    for (auto& c : cables_) { c.accept(v); }
  }
  void accept(ConstCablingVisitor& v) const { 
    v.visit(*this); 
    for (const auto& c :cables_) { c.accept(v); }
    }*/

};



#endif

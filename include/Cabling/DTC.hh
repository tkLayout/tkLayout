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
  typedef PtrVector<Cable> Container;
  //typedef PtrVector<Cable> Container;
  //Container cables_;

public:
  DTC(const std::string name, const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);

  // CABLE CONNECTED TO THE DTC
  const Container& cable() const { return cable_; }

  void addCable(Cable* c) { cable_.push_back(c); }

  const std::string name() const { return name_; }
  const Category& type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int slot() const { return slot_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }

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
  Property<int, Default> nCablesPerDTC;

  void accept(CablingVisitor& v) { 
    v.visit(*this); 
    for (auto& c : cables_) { c.accept(v); }
  }
  void accept(ConstCablingVisitor& v) const { 
    v.visit(*this); 
    for (const auto& c :cables_) { c.accept(v); }
    }*/

private:
  const int computePlotColor(const int phiSectorRef, const Category& type, const int slot) const;

  Container cable_;

  std::string name_;

  double phiSectorWidth_;
  int phiSectorRef_;
  Category type_;
  int slot_;
  bool isPositiveCablingSide_;

  int plotColor_;
};



#endif

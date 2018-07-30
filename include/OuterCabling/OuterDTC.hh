#ifndef OuterDTC_HH
#define OuterDTC_HH

#include <vector>
#include <string>

#include "OuterCable.hh"
#include "Module.hh"


class OuterDTC : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<OuterCable> Container;

public:
  OuterDTC(const std::string name, const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);

  // CABLE CONNECTED TO THE DTC
  const Container& cable() const { return cable_; }

  void addCable(OuterCable* c) { cable_.push_back(c); }

  const std::string name() const { return name_; }
  const Category& type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int slot() const { return slot_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }

  const int plotColor() const { return plotColor_; }


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

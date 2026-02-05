#ifndef OuterDTC_HH
#define OuterDTC_HH

#include <vector>
#include <string>

#include "OuterCable.hh"
#include "OuterBundle.hh"
#include "Module.hh"

namespace insur {

class OuterDTC : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef std::vector<OuterCable*> Container;

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

  const int getCMSSWId() const { return myCMSSWId_; } 
  const int plotColor() const { return plotColor_; }


private:
  const int computeCMSSWId(const int oneCablingSideId, const bool isPositiveCablingSide) const;
  const int computePlotColor(const int oneCablingSideId) const;
  const int computeOneCablingSideId(const int phiSectorRef, const Category& type, const int slot) const;

  Container cable_;

  std::string name_;

  double phiSectorWidth_;
  int phiSectorRef_;
  Category type_;
  int slot_;
  bool isPositiveCablingSide_;

  int myCMSSWId_;
  int plotColor_;
};

}

#endif

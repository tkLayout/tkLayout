#ifndef HVLINE_HH
#define HVLINE_HH

#include <vector>
#include <string>

#include "ITCabling/PowerChain.hh"
#include "Module.hh"


class HvLine : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<Module> Container;

public:
  HvLine(const std::string name);
  ~HvLine();

  // CABLE CONNECTED TO THE HV LINE
  const Container& modules() const { return modules_; }
  const int numModules() const { return modules_.size(); }
  void addModule(Module* m) { modules_.push_back(m); }

  const std::string name() const { return name_; }
  /*
  const Category& type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int slot() const { return slot_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }

  const int plotColor() const { return plotColor_; }
  */

  // POWER CHAIN, TO WHICH THE MODULES OF THE HVLINE ARE ALL CONNECTED
  const PowerChain* getPowerChain() const {
    if (!powerChain_) throw PathfulException("powerChain_ is nullptr");
    return powerChain_;
  }
  void setPowerChain(PowerChain* powerChain) { powerChain_ = powerChain; }



private:
  //const int computePlotColor(const int phiSectorRef, const Category& type, const int slot) const;

  Container modules_;

  std::string name_;

  PowerChain* powerChain_ = nullptr;

  /*
  double phiSectorWidth_;
  int phiSectorRef_;
  Category type_;
  int slot_;
  bool isPositiveCablingSide_;

  int plotColor_;
  */
};



#endif

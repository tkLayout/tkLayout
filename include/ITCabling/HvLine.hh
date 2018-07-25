#ifndef HVLINE_HH
#define HVLINE_HH

#include <vector>
#include <string>

#include "ITCabling/PowerChain.hh"
#include "Module.hh"


/*
 * High Voltage line class.
 * All modules connected to a given HV line can be accessed, as well as the Power Chain they are connected to.
 */
class HvLine : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<Module> Container;

public:
  HvLine(const std::string name);
  ~HvLine();

  // MODULES CONNECTED TO THE HV LINE
  const Container& modules() const { return modules_; }
  const int numModules() const { return modules_.size(); }
  void addModule(Module* m) { modules_.push_back(m); }

  const std::string name() const { return name_; }


  // POWER CHAIN, TO WHICH THE MODULES OF THE HVLINE ARE ALL CONNECTED
  const PowerChain* getPowerChain() const {
    if (!powerChain_) throw PathfulException("powerChain_ is nullptr");
    return powerChain_;
  }
  void setPowerChain(PowerChain* powerChain) { powerChain_ = powerChain; }



private:
  Container modules_;

  std::string name_;

  PowerChain* powerChain_ = nullptr;
};



#endif

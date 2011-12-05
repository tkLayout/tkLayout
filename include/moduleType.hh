#ifndef _MODULETYPE_HH_
#define _MODULETYPE_HH_

// Standard stuff
#include <vector>
#include <map>
#include <string>
#include <iostream>

class ModuleType {
protected:
  std::map<int,double> perStripPower_;
  std::map<int,double> perModulePower_;
  double triggerErrorX_;
  double triggerErrorY_;
private:
  bool checkPowerType(int powerType);
public:
  ModuleType();
  enum { OpticalPower, ChipPower, PowerTypes };
  void setPowerPerStrip(double W, int powerType);
  void setPowerPerModule(double W, int powerType);
  double getPower(int powerType, int nStrips);
  double getPower(int nStrips);
  double getPowerPerStrip(int powerType);
  double getPowerPerModule(int powerType);
  void setTriggerErrorX(double newError);
  void setTriggerErrorY(double newError);
  double getTriggerErrorX() const ;
  double getTriggerErrorY() const ;
};

#endif

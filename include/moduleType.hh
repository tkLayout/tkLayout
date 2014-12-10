#ifndef _MODULETYPE_HH_
#define _MODULETYPE_HH_

// Standard stuff
#include <vector>
#include <map>
#include <string>
#include <iostream>

#include "ptError.h"

class ModuleType {
protected:
  std::map<int,double> perStripPower_;
  std::map<int,double> perModulePower_;
  double triggerErrorX_;
  double triggerErrorY_;
  int sparsifiedHeaderBits_, sparsifiedPayloadBits_;
  int triggerDataHeaderBits_, triggerDataPayloadBits_;
  double sensorThickness_;
  ptError::InefficiencyType inefficiencyType_;
private:
  bool checkPowerType(int powerType);
public:
  ModuleType();
  enum { OpticalPower, ChipPower, PowerTypes };
  void setPowerPerStrip(double W, int powerType);
  void setPowerPerModule(double W, int powerType);
  double getPower(int powerType, int nStrips);
  double getPower(int nStrips) const;
  double getPowerPerStrip(int powerType);
  double getPowerPerModule(int powerType);
  void setTriggerErrorX(double newError);
  void setTriggerErrorY(double newError);
  double getTriggerErrorX() const ;
  double getTriggerErrorY() const ;
  void setSparsifiedHeaderBits(int bits);
  void setSparsifiedPayloadBits(int bits);
  int  getSparsifiedHeaderBits()  const;
  int  getSparsifiedPayloadBits()  const;
  void setTriggerDataHeaderBits(int bits);
  void setTriggerDataPayloadBits(int bits);
  int  getTriggerDataHeaderBits()  const;
  int  getTriggerDataPayloadBits() const;
  void setSensorThickness(double thickness);
  double getSensorThickness() const;
  void setInefficiencyType(ptError::InefficiencyType type);
  ptError::InefficiencyType getInefficiencyType() const;
};

#endif

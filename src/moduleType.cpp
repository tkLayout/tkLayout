#include <moduleType.hh>

ModuleType::ModuleType() {
  triggerErrorX_=1;
  triggerErrorY_=1;
  sensorThickness_ = 0.2;
  sparsifiedHeaderBits_  = 13;
  sparsifiedPayloadBits_ = 9;  
  triggerDataHeaderBits_  = 20;
  triggerDataPayloadBits_ = 20;  
}

bool ModuleType::checkPowerType(int powerType) {
  if ((powerType<0) || (powerType>PowerTypes)) {
    std::cerr << "ERROR: unknown power type " << powerType << std::endl;
    return false;
  }  
  return true;
}

void ModuleType::setPowerPerStrip(double W, int powerType) {
  if (!checkPowerType(powerType)) return;
  perStripPower_[powerType] = W;
}

void ModuleType::setPowerPerModule(double W, int powerType) {
  if (!checkPowerType(powerType)) return;
  perModulePower_[powerType] = W;
}

double ModuleType::getPowerPerModule(int powerType) {
  if (!checkPowerType(powerType)) return 0;
  return perModulePower_[powerType];
}

double ModuleType::getPowerPerStrip(int powerType) {
  if (!checkPowerType(powerType)) return 0;
  return perStripPower_[powerType];
}


double ModuleType::getPower(int powerType, int nStrips) {
  return (perStripPower_[powerType] * nStrips +
	  perModulePower_[powerType]);
}

double ModuleType::getPower(int nStrips) {
  double perStripPower = 0;
  double perModulePower = 0;

  for (std::map<int,double>::iterator it = perStripPower_.begin();
       it != perStripPower_.end(); ++it) {
    perStripPower += it->second;
  }

  for (std::map<int,double>::iterator it = perModulePower_.begin();
       it != perModulePower_.end(); ++it) {
    perModulePower += it->second;
  }

  return (perStripPower * nStrips + perModulePower);
}

// Set and get triggerErrorx and triggerErrory
// error multipliers for the trigger
void ModuleType::setTriggerErrorX(double newError) { triggerErrorX_ = newError; }
void ModuleType::setTriggerErrorY(double newError) { triggerErrorY_ = newError; }
double ModuleType::getTriggerErrorX() const { return triggerErrorX_ ; }
double ModuleType::getTriggerErrorY() const { return triggerErrorY_ ; }

void ModuleType::setSparsifiedHeaderBits(int bits)  { sparsifiedHeaderBits_  = bits; }
void ModuleType::setSparsifiedPayloadBits(int bits) { sparsifiedPayloadBits_ = bits; }
int  ModuleType::getSparsifiedHeaderBits()  const { return sparsifiedHeaderBits_;  }
int  ModuleType::getSparsifiedPayloadBits() const { return sparsifiedPayloadBits_; }

void ModuleType::setTriggerDataHeaderBits(int bits)  { triggerDataHeaderBits_  = bits; }
void ModuleType::setTriggerDataPayloadBits(int bits) { triggerDataPayloadBits_ = bits; }
int  ModuleType::getTriggerDataHeaderBits()  const { return triggerDataHeaderBits_;  }
int  ModuleType::getTriggerDataPayloadBits() const { return triggerDataPayloadBits_; }

void ModuleType::setSensorThickness(double thickness) { sensorThickness_ = thickness; }
double ModuleType::getSensorThickness() const { return sensorThickness_; }

#include <moduleType.hh>

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


/*
 * IrradiationMapsManager.cpp
 *
 *  Created on: 19/feb/2014
 *      Author: stefano
 */

#include "IrradiationMapsManager.hh"

IrradiationMapsManager::IrradiationMapsManager() {
}

IrradiationMapsManager::~IrradiationMapsManager() {
  irradiationMaps.clear();
}

void IrradiationMapsManager::addIrradiationMap(const IrradiationMap& newIrradiationMap) {
  irradiationMaps.insert(newIrradiationMap);
}

void IrradiationMapsManager::addIrradiationMap(std::string newIrradiationMapFile) {
  IrradiationMap newIrradiationMap (newIrradiationMapFile);
  addIrradiationMap(newIrradiationMap);
}

double IrradiationMapsManager::calculateIrradiationPower(const std::pair<double,double>& coordinates) const{
  double irradiation = 0;
  bool mapFound = false;

  //Iterate through the ordered map set untill find a proper map
  for(std::set<IrradiationMap>::const_iterator iter = irradiationMaps.cbegin(); iter != irradiationMaps.cend(); ++ iter) {
    if (iter->isInRegion(coordinates)) {
      irradiation = iter->calculateIrradiation(coordinates);
      mapFound = true;
      break;
    }
  }
  if(!mapFound) {
    logERROR("Error while calculating irradiation, a proper irradiation map is not found");
  }

  return irradiation;
}

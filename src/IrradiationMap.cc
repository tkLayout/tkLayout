/*
 * IrradiationMap.cpp
 *
 *  Created on: 19/feb/2014
 *      Author: Stefano Martina
 */

#include"IrradiationMap.hh"
#include"Units.hh"

IrradiationMap::IrradiationMap(std::string irradiationMapFile) :
      rhoMin (0),
      rhoMax (0),
      rhoBinWidth (0),
      rhoBinNum (0),
      zMin (0),
      zMax (0),
      zBinWidth (0),
      zBinNum (0),
      invFemUnit (1)
{
  if (! irradiationMapFile.empty()) {
    ingest(irradiationMapFile);
  }
}

IrradiationMap::IrradiationMap() : IrradiationMap("")
{}

void IrradiationMap::ingest(std::string irradiationMapFile) {
  std::string line;
  bool found_rhoMin = false;
  bool found_rhoMax = false;
  bool found_rhoBinWidth = false;
  bool found_rhoBinNum = false;
  bool found_zMin = false;
  bool found_zMax = false;
  bool found_zBinWidth = false;
  bool found_zBinNum = false;
  bool found_invFemUnit = false;
  double irradiationValue = 0;
  std::ifstream filein(irradiationMapFile);

  if (!filein.is_open()) {
    logERROR("Failed opening irradiation map file!");
  }

  while(std::getline(filein, line)) {
    //find rhoMin
    if (line.find(comp_rhoMin) == 0) {
      line.erase(0,comp_rhoMin.length());
      rhoMin = strtod(line.c_str(),NULL);
      found_rhoMin = true;
      continue;
    }
    //find rhoMax
    if (line.find(comp_rhoMax) == 0) {
      line.erase(0,comp_rhoMax.length());
      rhoMax = strtod(line.c_str(),NULL);
      found_rhoMax = true;
      continue;
    }
    //find rhoBinWidth
    if (line.find(comp_rhoBinWidth) == 0) {
      line.erase(0,comp_rhoBinWidth.length());
      rhoBinWidth = strtod(line.c_str(),NULL);
      found_rhoBinWidth = true;
      continue;
    }
    //find rhoBinNum
    if (line.find(comp_rhoBinNum) == 0) {
      line.erase(0,comp_rhoBinNum.length());
      rhoBinNum = strtol(line.c_str(),NULL,10);
      found_rhoBinNum = true;
      continue;
    }
    //find zMin
    if (line.find(comp_zMin) == 0) {
      line.erase(0,comp_zMin.length());
      zMin = strtod(line.c_str(),NULL);
      found_zMin = true;
      continue;
    }
    //find zMax
    if (line.find(comp_zMax) == 0) {
      line.erase(0,comp_zMax.length());
      zMax = strtod(line.c_str(),NULL);
      found_zMax = true;
      continue;
    }
    //find zBinWidth
    if (line.find(comp_zBinWidth) == 0) {
      line.erase(0,comp_zBinWidth.length());
      zBinWidth = strtod(line.c_str(),NULL);
      found_zBinWidth = true;
      continue;
    }
    //find zBinNum
    if (line.find(comp_zBinNum) == 0) {
      line.erase(0,comp_zBinNum.length());
      zBinNum = strtol(line.c_str(),NULL,10);
      found_zBinNum = true;
      continue;
    }
    //find invFemUnit
    if (line.find(comp_invFemUnit) == 0) {
      line.erase(0,comp_invFemUnit.length());
      invFemUnit = strtod(line.c_str(),NULL);
      found_invFemUnit = true;
      continue;
    }

    //skip other comment or empty lines
    if (line.find_first_of("#//;")==0 || line=="") continue;

    //create a stream for reading values
    std::stringstream ss(line);
    //vector for storing a row
    std::vector<double> irradiationLine;
    while(!ss.eof()){
      //read a value
      ss >> irradiationValue;
      irradiationValue /= invFemUnit;
      if(!ss.eof()) {
        //add value to vector
        irradiationLine.push_back(irradiationValue);
      }
    }
    //add vector to matrix
    irradiation.push_back(irradiationLine);
  }

  //convert cm to mm
  zMin       *= Units::cm; //10;
  zMax       *= Units::cm; //10;
  zBinWidth  *= Units::cm; //10;
  rhoMin     *= Units::cm; //10;
  rhoMax     *= Units::cm; //10;
  rhoBinWidth*= Units::cm; //10;

  //change the min and max values because we are interested in the center points of the bins
  zMin += zBinWidth / 2;
  zMax -= zBinWidth / 2;
  rhoMin += rhoBinWidth / 2;
  rhoMax -= rhoBinWidth / 2;

  //control if found all values
  if (!found_rhoMin || !found_rhoMax || !found_rhoBinWidth ||
      !found_rhoBinNum || !found_zMin || !found_zMax ||
      !found_zBinWidth || !found_zBinNum || !found_invFemUnit) {
    std::ostringstream tempSS;
    tempSS << "Error while parsing irradiation map values: found_rMin " << found_rhoMin
        << "; found_rMax " << found_rhoMax << "; found_rBinWidth " << found_rhoBinWidth
        << "; found_rBinNum " << found_rhoBinNum << "; found_zMin " << found_zMin
        << "; found_zMax " << found_zMax << "; found_zBinWidth " << found_zBinWidth
        << "; found_zBinNum " << found_zBinNum << "; found_invFemUnit " << found_invFemUnit;
    logERROR(tempSS);
  }
}

double IrradiationMap::binArea() const{
  return zBinWidth*rhoBinWidth;
}

bool IrradiationMap::operator < (const IrradiationMap& confrontedMap) const{
  return (binArea() < confrontedMap.binArea());
}

bool IrradiationMap::isInRegion(const std::pair<double,double>& coordinates) const {
  return ((zMin <= coordinates.first) && (zMax >= coordinates.first) && (rhoMin <= coordinates.second) && (rhoMax >= coordinates.second));
}

double IrradiationMap::calculateIrradiation(const std::pair<double,double>& coordinates) const {
  double z = 0;
  double rho = 0;
  double z1 = 0;
  double z2 = 0;
  double rho1 = 0;
  double rho2 = 0;
  double irr1 = 0;
  double irr2 = 0;
  double irr11 = 0;
  double irr21 = 0;
  double irr12 = 0;
  double irr22 = 0;
  double irrxy = 0;

  if(isInRegion(coordinates)) {
    //correct the coordinates in the map matrix reference
    z = (coordinates.first - zMin) / zBinWidth;
    rho = (coordinates.second - rhoMin) / rhoBinWidth;

    //take the 4 nearest bin centers
    z1 = floor(z);
    z2 = ceil(z);
    rho1 = floor(rho);
    rho2 = ceil(rho);

    //if the point is in a intersection of the grid formed by the map bin centers
    if((z1 == z2) && (rho1 == rho2)) {
      //single value
      irrxy = irradiation[int(rho1)][int(z1)];
    }

    //if is in a z line
    else if (z1 == z2) {
      irr1 = irradiation[int(rho1)][int(z1)];
      irr2 = irradiation[int(rho2)][int(z1)];
      //linear interpolation in rho
      irrxy = irr1/(rho2-rho1) * (rho2-rho) + irr2/(rho2-rho1) * (rho-rho1);
    }

    //if is in a rho line
    else if (rho1 == rho2) {
      irr1 = irradiation[int(rho1)][int(z1)];
      irr2 = irradiation[int(rho1)][int(z2)];
      //linear interpolation in z
      irrxy = irr1/(z2-z1) * (z2-z) + irr2/(z2-z1) * (z-z1);
    }

    //if is in the middle
    else {
      irr11 = irradiation[int(rho1)][int(z1)];
      irr21 = irradiation[int(rho1)][int(z2)];
      irr12 = irradiation[int(rho2)][int(z1)];
      irr22 = irradiation[int(rho2)][int(z2)];
      //bilinear interpolation in z and rho
      irrxy = irr11/((z2-z1)*(rho2-rho1))*(z2-z)*(rho2-rho) + irr21/((z2-z1)*(rho2-rho1))*(z-z1)*(rho2-rho) + irr12/((z2-z1)*(rho2-rho1))*(z2-z)*(rho-rho1) + irr22/((z2-z1)*(rho2-rho1))*(z-z1)*(rho-rho1);
    }
  } else {
    logERROR("Error while calculating module irradiation, module out of region");
  }

  return irrxy;
}

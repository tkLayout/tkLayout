/*
 * IrradiationMap.cpp
 *
 *  Created on: 19/feb/2014
 *      Author: Stefano Martina, Zbynek Drasal
 */

#include"IrradiationMap.h"
#include <algorithm>
#include <Units.h>
#include <cmath>
#include "MessageLogger.h"

IrradiationMap::IrradiationMap(std::string irradiationMapFile) :
      m_fileName(irradiationMapFile),
      m_dataUnit("cm^-2"),
      m_dataFactor(1./Units::cm2),
      m_rUnit("cm"),
      m_rMin (0),
      m_rMax (0),
      m_rBinWidth (0),
      m_rBinNum (0),
      m_zUnit("cm"),
      m_zMin (0),
      m_zMax (0),
      m_zBinWidth (0),
      m_zBinNum (0),
      m_norm (1),
      m_typeMesh(false),
      m_typeHist(false)
{
  if (! m_fileName.empty()) {
    ingest(m_fileName);
  }
}

IrradiationMap::IrradiationMap() : IrradiationMap("")
{}

void IrradiationMap::ingest(std::string fileName) {
  std::string line;
  bool   found_rMin       = false;
  bool   found_rMax       = false;
  bool   found_rBinWidth  = false;
  bool   found_rBinNum    = false;
  bool   found_zMin       = false;
  bool   found_zMax       = false;
  bool   found_zBinWidth  = false;
  bool   found_zBinNum    = false;
  bool   found_norm       = false;
  double irradiationValue = 0;

  std::string m_fileName = fileName;
  std::ifstream filein(m_fileName);

  if (!filein.is_open()) {
    std::ostringstream message;
    message << "Failed opening map file: " << fileName << "!";
    logERROR(message.str());
  }

  // Check that all data correctly set to be able to read map
  bool init = false;

  while(std::getline(filein, line)) {
    //find fluxUnit
    if (line.find(comp_dataUnit) == 0) {
      line.erase(0,comp_dataUnit.length());
      m_dataUnit = line;
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), comp_EscLine.c_str() ), m_dataUnit.end());
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), '\r'), m_dataUnit.end());
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), '\n'), m_dataUnit.end());
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), ' ' ), m_dataUnit.end());
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), comp_EscValue.c_str()), m_dataUnit.end());
      continue;
    }
    //find rUnit
    if (line.find(comp_rUnit) == 0) {
      line.erase(0,comp_rUnit.length());
      m_rUnit = line;
      m_rUnit.erase(std::remove(m_rUnit.begin(), m_rUnit.end(), comp_EscLine.c_str() ), m_rUnit.end());
      m_rUnit.erase(std::remove(m_rUnit.begin(), m_rUnit.end(), '\r'), m_rUnit.end());
      m_rUnit.erase(std::remove(m_rUnit.begin(), m_rUnit.end(), '\n'), m_rUnit.end());
      m_rUnit.erase(std::remove(m_rUnit.begin(), m_rUnit.end(), ' ' ), m_rUnit.end());
      m_rUnit.erase(std::remove(m_rUnit.begin(), m_rUnit.end(), comp_EscValue.c_str()), m_rUnit.end());
      continue;
    }
    //find rMin
    if (line.find(comp_rMin) == 0) {
      line.erase(0,comp_rMin.length());
      m_rMin = strtod(line.c_str(),NULL);
      found_rMin = true;
      continue;
    }
    //find rMax
    if (line.find(comp_rMax) == 0) {
      line.erase(0,comp_rMax.length());
      m_rMax = strtod(line.c_str(),NULL);
      found_rMax = true;
      continue;
    }
    //find rBinWidth
    if (line.find(comp_rBinWidth) == 0) {
      line.erase(0,comp_rBinWidth.length());
      m_rBinWidth = strtod(line.c_str(),NULL);
      found_rBinWidth = true;
      continue;
    }
    //find rBinNum
    if (line.find(comp_rBinNum) == 0) {
      line.erase(0,comp_rBinNum.length());
      m_rBinNum = strtol(line.c_str(),NULL,10);
      found_rBinNum = true;
      continue;
    }
    //find zUnit
    if (line.find(comp_zUnit) == 0) {
      line.erase(0,comp_zUnit.length());
      m_zUnit = line;
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), comp_EscLine.c_str() ), m_zUnit.end());
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), '\r'), m_zUnit.end());
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), '\n'), m_zUnit.end());
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), ' ' ), m_zUnit.end());
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), comp_EscValue.c_str()), m_zUnit.end());
      continue;
    }
    //find zMin
    if (line.find(comp_zMin) == 0) {
      line.erase(0,comp_zMin.length());
      m_zMin = strtod(line.c_str(),NULL);
      found_zMin = true;
      continue;
    }
    //find zMax
    if (line.find(comp_zMax) == 0) {
      line.erase(0,comp_zMax.length());
      m_zMax = strtod(line.c_str(),NULL);
      found_zMax = true;
      continue;
    }
    //find zBinWidth
    if (line.find(comp_zBinWidth) == 0) {
      line.erase(0,comp_zBinWidth.length());
      m_zBinWidth = strtod(line.c_str(),NULL);
      found_zBinWidth = true;
      continue;
    }
    //find zBinNum
    if (line.find(comp_zBinNum) == 0) {
      line.erase(0,comp_zBinNum.length());
      m_zBinNum = strtol(line.c_str(),NULL,10);
      found_zBinNum = true;
      continue;
    }
    //find normalization (invFemUnit)
    if (line.find(comp_norm) == 0) {
      line.erase(0,comp_norm.length());
      m_norm = strtol(line.c_str(),NULL,10);
      found_norm = true;
      continue;
    }

    //skip other comment or empty lines
    if (line.find_first_of(comp_EscComment)==0 || line=="") continue;

    // Check that all values defined
    if (!init) {

      init = true;

      // Check if found all values
      if (!found_rMin      || !found_rMax    || !found_rBinWidth ||
          !found_rBinNum   || !found_zMin    || !found_zMax ||
          !found_zBinWidth || !found_zBinNum) {
        std::ostringstream tempSS;
        tempSS << "Error while parsing irradiation map values: found_rMin "  << found_rMin
            << "; found_rMax "      << found_rMax    << "; found_rBinWidth " << found_rBinWidth
            << "; found_rBinNum "   << found_rBinNum << "; found_zMin "      << found_zMin
            << "; found_zMax "      << found_zMax    << "; found_zBinWidth " << found_zBinWidth
            << "; found_zBinNum "   << found_zBinNum;
        logERROR(tempSS);

        // Irradiation field will be empty
        break;
      }
    }

    //create a stream for reading values
    std::stringstream ss(line);
    //vector for storing a row
    std::vector<double> irradiationLine;
    while(!ss.eof()){
      //read a value
      ss >> irradiationValue;
      irradiationValue /= m_norm;
      if(!ss.eof()) {
        //add value to vector
        irradiationLine.push_back(irradiationValue);
      }
    }
    //add vector to matrix
    m_irradiation.push_back(irradiationLine);
  }

  //conversion (if units not defined in a file, default cm assumed -> conversion to software default mm done)
  if      (m_dataUnit=="mm^-2") m_dataFactor = 1./Units::mm2;
  else if (m_dataUnit=="cm^-2") m_dataFactor = 1./Units::cm2;
  else if (m_dataUnit=="m^-2")  m_dataFactor = 1./Units::m2;
  else {
    std::ostringstream message;
    message << "IrradiationMap: Flux units " << m_dataUnit << " not recognized. Used cm^-2 instead!";
    logWARNING(message.str());
  }

  double rFactor = 1;
  if      (m_rUnit=="mm") rFactor = Units::mm;
  else if (m_rUnit=="cm") rFactor = Units::cm;
  else if (m_rUnit=="m" ) rFactor = Units::m;
  else {
    std::ostringstream message;
    message << "IrradiationMap: Radius units " << m_rUnit << " not recognized. Used cm instead!";
    logWARNING(message.str());
  }
  m_rMin      *= rFactor;
  m_rMax      *= rFactor;
  m_rBinWidth *= rFactor;

  double zFactor = 1;
  if      (m_zUnit=="mm") zFactor = Units::mm;
  else if (m_zUnit=="cm") zFactor = Units::cm;
  else if (m_zUnit=="m" ) zFactor = Units::m;
  else {
    std::ostringstream message;
    message << "IrradiationMap: Z units " << m_zUnit << " not recognized. Used cm instead!";
    logWARNING(message.str());
  }
  m_zMin      *= zFactor;
  m_zMax      *= zFactor;
  m_zBinWidth *= zFactor;

  // Check binning
  m_typeMesh = false;
  m_typeHist = false;
  if (long((m_zMax-m_zMin)/m_zBinWidth)==m_zBinNum    ) m_typeHist = true;
  if (long((m_zMax-m_zMin)/m_zBinWidth)==(m_zBinNum-1)) m_typeMesh = true;

  std::ostringstream message;
  if (m_typeMesh) {
    message << "Irradiation map of type mesh: " << m_fileName << " read in!";
    logINFO(message.str());
  }
  if (m_typeHist) {
    message << "Irradiation map of type histogram: " << m_fileName << " read in!";
    logINFO(message.str());
  }
  if (!m_typeMesh && !m_typeHist) {
    message << "Irradiation map: " << m_fileName << " - binning doesn't correspond to the defined range and the bin size!";
    logERROR(message.str());
  }
}

double IrradiationMap::binArea() const{
  return m_zBinWidth*m_rBinWidth;
}

bool IrradiationMap::operator < (const IrradiationMap& confrontedMap) const{
  return (binArea() < confrontedMap.binArea());
}

bool IrradiationMap::isInRegionZR(double zPos, double rPos) const {
  return ((m_zMin <= zPos) && (m_zMax >= zPos) && (m_rMin <= rPos) && (m_rMax >= rPos));
}

bool IrradiationMap::isInBinRegionZR(int& zPos, int& rPos) const {
  return ((zPos>=0) && (zPos<=m_zBinNum) && (rPos>=0) && (rPos<=m_rBinNum));
}

double IrradiationMap::calculateIrradiationZR(double zPos, double rPos) const {
  double z     = 0;
  double r     = 0;
  int    z1    = 0;
  int    z2    = 0;
  int    r1    = 0;
  int    r2    = 0;
  double irr1  = 0;
  double irr2  = 0;
  double irr11 = 0;
  double irr21 = 0;
  double irr12 = 0;
  double irr22 = 0;
  double irrxy = 0;

  // Check that irradiation map was correctly filled
  if (m_irradiation.size()==0) return 0;

  // Mesh type - refer to mesh points
  if (m_typeMesh) {
    z = (zPos - m_zMin) / m_zBinWidth;
    r = (rPos - m_rMin) / m_rBinWidth;
  }
  // Histogram type - always refer to the bin center
  if (m_typeHist) {
    z = (zPos - (m_zMin + m_zBinWidth/2.)) / m_zBinWidth;
    r = (rPos - (m_rMin + m_rBinWidth/2.)) / m_rBinWidth;
  }
  // Take the 4 nearest bin centers & interpolate
  z1 = int(floor(z));
  z2 = int(ceil(z));
  r1 = int(floor(r));
  r2 = int(ceil(r));

  // Take care about first and list bins in histogram -> no interpolation
  if (m_typeHist) {

    if ( z>=-1/2. && z<=0 ) z1 = z2;
    if ( r>=-1/2. && r<=0 ) r1 = r2;

    if ( z>=(m_zBinNum-1) && z<=((m_zBinNum-1)+1/2.) ) z2 = z1;
    if ( r>=(m_rBinNum-1) && r<=((m_rBinNum-1)+1/2.) ) r2 = r1;
  }

  // Avoid running out of range for max values -> use limit values instead and warn user!
  if (z2>=m_zBinNum) {
    z2 = m_zBinNum-1;
    z1 = m_zBinNum-1;
    std::ostringstream message;
    message << "Value out of range when using: " << m_fileName << "Irradiation calculated for Z= " << (z2+1)*(m_zBinWidth)+m_zMin << " instead for Z= " << zPos;
    logWARNING(message.str());
    zPos = (z2+1)*(m_zBinWidth)+m_zMin;
  }
  if (r2>=m_rBinNum) {
    r2 = m_rBinNum-1;
    r1 = m_rBinNum-1;
    std::ostringstream message;
    message << "Value out of range when using: " << m_fileName << "Irradiation calculated for R= " << (r2+1)*(m_rBinWidth)+m_rMin << " instead for R= " << rPos;
    logWARNING(message.str());
    rPos = (r2+1)*(m_rBinWidth)+m_rMin;
  }

  if (isInRegionZR(zPos, rPos)) {
    if (isInBinRegionZR(z1,r1) && isInBinRegionZR(z2,r2)) {

      //if the point is in a intersection of the grid formed by the map bin centers
      if((z1 == z2) && (r1 == r2)) {
        //single value
        irrxy = m_irradiation[z1][r1];
      }

      //if is in a z line
      else if (z1 == z2) {
        irr1 = m_irradiation[z1][r1];
        irr2 = m_irradiation[z1][r2];
        //linear interpolation in r
        irrxy = irr2/(r2-r1)*(r-r1) + irr1/(r2-r1)*(r2-r);
      }

      //if is in a r line
      else if (r1 == r2) {
        irr1 = m_irradiation[z1][r1];
        irr2 = m_irradiation[z2][r1];
        //linear interpolation in z
        irrxy = irr2/(z2-z1)*(z-z1) + irr1/(z2-z1)*(z2-z);
      }

      //if is in the middle
      else {
        irr11 = m_irradiation[z1][r1];
        irr21 = m_irradiation[z2][r1];
        irr12 = m_irradiation[z1][r2];
        irr22 = m_irradiation[z2][r2];
        //bilinear interpolation in z and r
        irrxy = irr11/((z2-z1)*(r2-r1))*(z2-z)*(r2-r) + irr21/((z2-z1)*(r2-r1))*(z-z1)*(r2-r) + irr12/((z2-z1)*(r2-r1))*(z2-z)*(r-r1) + irr22/((z2-z1)*(r2-r1))*(z-z1)*(r-r1);
      }
    }
    else logERROR("Irradiation map: (Z,R) bin out of region");
  }
  else logERROR("Irradiation map: (Z,R) position out of region");

  // Convert output to
  return irrxy*m_dataFactor;
}

double IrradiationMap::calculateIrradiationRZ(double rPos, double zPos) const {
  double z     = 0;
  double r     = 0;
  int    z1    = 0;
  int    z2    = 0;
  int    r1    = 0;
  int    r2    = 0;
  double irr1  = 0;
  double irr2  = 0;
  double irr11 = 0;
  double irr21 = 0;
  double irr12 = 0;
  double irr22 = 0;
  double irrxy = 0;

  // Check that irradiation map was correctly filled
  if (m_irradiation.size()==0) return 0;

  // Mesh type - refer to mesh points
  if (m_typeMesh) {
    z = (zPos - m_zMin) / m_zBinWidth;
    r = (rPos - m_rMin) / m_rBinWidth;
  }
  // Histogram type - always refer to the bin center
  if (m_typeHist) {
    z = (zPos - (m_zMin + m_zBinWidth/2.)) / m_zBinWidth;
    r = (rPos - (m_rMin + m_rBinWidth/2.)) / m_rBinWidth;
  }
  // Take the 4 nearest bin centers & interpolate
  z1 = int(floor(z));
  z2 = int(ceil(z));
  r1 = int(floor(r));
  r2 = int(ceil(r));

  // Take care about first and list bins in histogram -> no interpolation
  if (m_typeHist) {
    if (z>=-1/2. && z<=0) z1 = z2;
    if (r>=-1/2. && z<=0) r1 = r2;

    if ( z>=(m_zBinNum-1) && z<=((m_zBinNum-1)+1/2.) ) z2 = z1;
    if ( r>=(m_rBinNum-1) && r<=((m_rBinNum-1)+1/2.) ) r2 = r1;
  }

  // Avoid running out of range for max values -> use limit values instead and warn user!
  if (z2>=m_zBinNum) {
    z2 = m_zBinNum-1;
    z1 = m_zBinNum-1;
    std::ostringstream message;
    message << "Value out of range when using: " << m_fileName << "Irradiation calculated for Z= " << (z2+1)*(m_zBinWidth)+m_zMin << " instead for Z= " << zPos;
    logWARNING(message.str());
    zPos = (z2+1)*(m_zBinWidth)+m_zMin;
  }
  if (r2>=m_rBinNum) {
    r2 = m_rBinNum-1;
    r1 = m_rBinNum-1;
    std::ostringstream message;
    message << "Value out of range when using: " << m_fileName << "Irradiation calculated for R= " << (r2+1)*(m_rBinWidth)+m_rMin << " instead for R= " << rPos;
    logWARNING(message.str());
    rPos = (r2+1)*(m_rBinWidth)+m_rMin;
  }

  if(isInRegionZR(zPos, rPos)) {
    if (isInBinRegionZR(z1,r1) && isInBinRegionZR(z2,r2)) {

      //if the point is in a intersection of the grid formed by the map bin centers
      if((z1 == z2) && (r1 == r2)) {
        //single value
        irrxy = m_irradiation[r1][z1];
      }

      //if is in a z line
      else if (z1 == z2) {
        irr1 = m_irradiation[r1][z1];
        irr2 = m_irradiation[r2][z1];
        //linear interpolation in r
        irrxy = irr2/(r2-r1)*(r-r1) + irr1/(r2-r1)*(r2-r);
      }

      //if is in a r line
      else if (r1 == r2) {
        irr1 = m_irradiation[r1][z1];
        irr2 = m_irradiation[r1][z2];
        //linear interpolation in z
        irrxy = irr2/(z2-z1)*(z-z1) + irr1/(z2-z1)*(z2-z);
      }

      //if is in the middle
      else {
        irr11 = m_irradiation[r1][z1];
        irr21 = m_irradiation[r1][z2];
        irr12 = m_irradiation[r2][z1];
        irr22 = m_irradiation[r2][z2];
        //bilinear interpolation in z and r
        irrxy = irr11/((z2-z1)*(r2-r1))*(z2-z)*(r2-r) + irr21/((z2-z1)*(r2-r1))*(z-z1)*(r2-r) + irr12/((z2-z1)*(r2-r1))*(z2-z)*(r-r1) + irr22/((z2-z1)*(r2-r1))*(z-z1)*(r-r1);
      }
    }
    else logERROR("Iradiation map: (R,Z) bin out of region");
  }
  else logWARNING("Irradiation map: (R,Z) position out of region.");

  return irrxy*m_dataFactor;
}


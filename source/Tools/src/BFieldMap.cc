/*
 * BFieldMap.cpp
 *
 *  Created on: 03/Dec/2015
 *      Author: Zbynek Drasal
 */

#include"BFieldMap.h"
#include <algorithm>
#include <Units.h>
#include "MessageLogger.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TList.h"

BFieldMap::BFieldMap(std::string bFieldMapFile) :
      m_fileName(bFieldMapFile),
      m_dataUnit("T"),
      m_dataFactor(1./Units::T),
      m_xUnit("cm"),
      m_xMin (0),
      m_xMax (0),
      m_xBinWidth (0),
      m_xBinNum (0),
      m_yUnit("cm"),
      m_yMin (0),
      m_yMax (0),
      m_yBinWidth (0),
      m_yBinNum (0),
      m_zUnit("cm"),
      m_zMin (0),
      m_zMax (0),
      m_zBinWidth (0),
      m_zBinNum (0),
      m_typeMesh(false),
      m_typeHist(false),
      m_bFieldOK(false)
{
  if (!m_fileName.empty()) {
    ingest(m_fileName);
  }
}

void BFieldMap::ingest(std::string fileName) {
  std::string line;
  bool    found_xMin       = false;
  bool    found_xMax       = false;
  bool    found_xBinWidth  = false;
  bool    found_xBinNum    = false;
  bool    found_yMin       = false;
  bool    found_yMax       = false;
  bool    found_yBinWidth  = false;
  bool    found_yBinNum    = false;
  bool    found_zMin       = false;
  bool    found_zMax       = false;
  bool    found_zBinWidth  = false;
  bool    found_zBinNum    = false;

  std::vector<double> bFieldVector;
  bFieldVector.resize(3);

  m_fileName = fileName;
  std::ifstream filein(m_fileName);

  if (!filein.is_open()) {
    std::ostringstream message;
    message << "Failed opening B field file: " << fileName << "!";
    logERROR(message.str());
  }

  // Initialize bField variable
  bool init   = false;
  // Before reading file, B field indicis
  double prevXValue = -std::numeric_limits<double>::max();
  double prevYValue = -std::numeric_limits<double>::max();
  int    indexX = -1;
  int    indexY = -1;
  int    indexZ = -1;

  while(std::getline(filein, line)) {
    //find bFieldUnit
    if (line.find(c_fileDataUnit) == 0) {
      line.erase(0,c_fileDataUnit.length());
      m_dataUnit = line;
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), c_fileEscLine.c_str() ), m_dataUnit.end());
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), '\r'), m_dataUnit.end());
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), '\n'), m_dataUnit.end());
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), ' ' ), m_dataUnit.end());
      m_dataUnit.erase(std::remove(m_dataUnit.begin(), m_dataUnit.end(), c_fileEscValue.c_str()), m_dataUnit.end());
      continue;
    }
    //find xUnit
    if (line.find(c_fileXUnit) == 0) {
      line.erase(0,c_fileXUnit.length());
      m_xUnit = line;
      m_xUnit.erase(std::remove(m_xUnit.begin(), m_xUnit.end(), c_fileEscLine.c_str() ), m_xUnit.end());
      m_xUnit.erase(std::remove(m_xUnit.begin(), m_xUnit.end(), '\r'), m_xUnit.end());
      m_xUnit.erase(std::remove(m_xUnit.begin(), m_xUnit.end(), '\n'), m_xUnit.end());
      m_xUnit.erase(std::remove(m_xUnit.begin(), m_xUnit.end(), ' ' ), m_xUnit.end());
      m_xUnit.erase(std::remove(m_xUnit.begin(), m_xUnit.end(), c_fileEscValue.c_str()), m_xUnit.end());
      continue;
    }
    //find xMin
    if (line.find(c_fileXMin) == 0) {
      line.erase(0,c_fileXMin.length());
      m_xMin = strtod(line.c_str(),NULL);
      found_xMin = true;
      continue;
    }
    //find xMax
    if (line.find(c_fileXMax) == 0) {
      line.erase(0,c_fileXMax.length());
      m_xMax = strtod(line.c_str(),NULL);
      found_xMax = true;
      continue;
    }
    //find xBinWidth
    if (line.find(c_fileXBinWidth) == 0) {
      line.erase(0,c_fileXBinWidth.length());
      m_xBinWidth = strtod(line.c_str(),NULL);
      found_xBinWidth = true;
      continue;
    }
    //find xBinNum
    if (line.find(c_fileXBinNum) == 0) {
      line.erase(0,c_fileXBinNum.length());
      m_xBinNum = strtol(line.c_str(),NULL,10);
      found_xBinNum = true;
      continue;
    }
    //find yUnit
    if (line.find(c_fileYUnit) == 0) {
      line.erase(0,c_fileYUnit.length());
      m_yUnit = line;
      m_yUnit.erase(std::remove(m_yUnit.begin(), m_yUnit.end(), c_fileEscLine.c_str() ), m_yUnit.end());
      m_yUnit.erase(std::remove(m_yUnit.begin(), m_yUnit.end(), '\r'), m_yUnit.end());
      m_yUnit.erase(std::remove(m_yUnit.begin(), m_yUnit.end(), '\n'), m_yUnit.end());
      m_yUnit.erase(std::remove(m_yUnit.begin(), m_yUnit.end(), ' ' ), m_yUnit.end());
      m_yUnit.erase(std::remove(m_yUnit.begin(), m_yUnit.end(), c_fileEscValue.c_str()), m_yUnit.end());
      continue;
    }
    //find xMin
    if (line.find(c_fileYMin) == 0) {
      line.erase(0,c_fileYMin.length());
      m_yMin = strtod(line.c_str(),NULL);
      found_yMin = true;
      continue;
    }
    //find xMax
    if (line.find(c_fileYMax) == 0) {
      line.erase(0,c_fileYMax.length());
      m_yMax = strtod(line.c_str(),NULL);
      found_yMax = true;
      continue;
    }
    //find xBinWidth
    if (line.find(c_fileYBinWidth) == 0) {
      line.erase(0,c_fileYBinWidth.length());
      m_yBinWidth = strtod(line.c_str(),NULL);
      found_yBinWidth = true;
      continue;
    }
    //find xBinNum
    if (line.find(c_fileYBinNum) == 0) {
      line.erase(0,c_fileYBinNum.length());
      m_yBinNum = strtol(line.c_str(),NULL,10);
      found_yBinNum = true;
      continue;
    }
    //find m_zUnit
    if (line.find(c_fileZUnit) == 0) {
      line.erase(0,c_fileZUnit.length());
      m_zUnit = line;
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), c_fileEscLine.c_str() ), m_zUnit.end());
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), '\r'), m_zUnit.end());
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), '\n'), m_zUnit.end());
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), ' ' ), m_zUnit.end());
      m_zUnit.erase(std::remove(m_zUnit.begin(), m_zUnit.end(), c_fileEscValue.c_str()), m_zUnit.end());
      continue;
    }
    //find m_zMin
    if (line.find(c_fileZMin) == 0) {
      line.erase(0,c_fileZMin.length());
      m_zMin = strtod(line.c_str(),NULL);
      found_zMin = true;
      continue;
    }
    //find m_zMax
    if (line.find(c_fileZMax) == 0) {
      line.erase(0,c_fileZMax.length());
      m_zMax = strtod(line.c_str(),NULL);
      found_zMax = true;
      continue;
    }
    //find m_zBinWidth
    if (line.find(c_fileZBinWidth) == 0) {
      line.erase(0,c_fileZBinWidth.length());
      m_zBinWidth = strtod(line.c_str(),NULL);
      found_zBinWidth = true;
      continue;
    }
    //find m_zBinNum
    if (line.find(c_fileZBinNum) == 0) {
      line.erase(0,c_fileZBinNum.length());
      m_zBinNum = strtol(line.c_str(),NULL,10);
      found_zBinNum = true;
      continue;
    }

    //skip other comment or empty lines
    if (line.find_first_of(c_fileEscComment)==0 || line=="") continue;

    // Check that all values properly define & init
    if (!init) {

      init = true;
      // Check if found all values
      if (!found_xMin      || !found_xMax    || !found_xBinWidth || !found_xBinNum ||
          !found_yMin      || !found_yMax    || !found_yBinWidth || !found_yBinNum ||
          !found_zMin      || !found_zMax    || !found_zBinWidth || !found_zBinNum) {

        std::ostringstream tempSS;
        tempSS << "Error while parsing B field map values: found_xMin "
            << found_xMin << "; found_xMax " << found_xMax << "; found_xBinWidth " << found_xBinWidth << "; found_xBinNum " << found_xBinNum << "; found_yMin "
            << found_yMin << "; found_yMax " << found_yMax << "; found_yBinWidth " << found_yBinWidth << "; found_yBinNum " << found_yBinNum << "; found_xMin "
            << found_zMin << "; found_zMax " << found_zMax << "; found_zBinWidth " << found_zBinWidth << "; found_zBinNum " << found_zBinNum;
        logERROR(tempSS);

        // Irradiation field will be empty
        m_bFieldOK = false;
        break;
      }
      else {

        prevXValue = -std::numeric_limits<double>::max();
        prevYValue = -std::numeric_limits<double>::max();
        indexX = -1;
        indexY = -1;
        indexZ = -1;

        m_bField.resize(m_xBinNum);
        for (int x=0; x<m_xBinNum; x++) {
          m_bField[x].resize(m_yBinNum);
          for (int y=0; y<m_yBinNum; y++) {
            m_bField[x][y].resize(m_zBinNum);
            for (int z=0; z<m_zBinNum; z++) m_bField[x][y][z].resize(3);
          }
        }

        m_bFieldOK = true;
      }
    }

    // Read values
    double valueX;
    double valueY;
    double valueZ;
    std::vector<double> bField;
    bField.resize(3);

    std::stringstream ss(line);
    ss >> valueX;
    ss >> valueY;
    ss >> valueZ;
    ss >> bField[0];
    ss >> bField[1];
    ss >> bField[2];

    if (valueX!=prevXValue) {

      indexX++;
      indexY = 0;
      indexZ = 0;

      if (indexX>=m_xBinNum || indexY>=m_yBinNum || indexZ>=m_zBinNum) {
        logERROR("Error while parsing B field map - number of values in file not compatible with given number of bins in header");
        m_bFieldOK = false;
        break;
      }
      else {
        m_bField[indexX][indexY][indexZ][0] = bField[0];
        m_bField[indexX][indexY][indexZ][1] = bField[1];
        m_bField[indexX][indexY][indexZ][2] = bField[2];
        prevXValue = valueX;
        prevYValue = valueY;
      }
    }
    else if (valueY!=prevYValue) {

      indexY++;
      indexZ = 0;

      if (indexY>=m_yBinNum || indexZ>=m_zBinNum) {
        logERROR("Error while parsing B field map - number of values in file not compatible with given number of bins in header");
        m_bFieldOK = false;
        break;
      }
      else {
        m_bField[indexX][indexY][indexZ][0] = bField[0];
        m_bField[indexX][indexY][indexZ][1] = bField[1];
        m_bField[indexX][indexY][indexZ][2] = bField[2];
        prevYValue = valueY;
      }
    }
    else {
      indexZ++;

      if (indexZ>=m_zBinNum) {
        logERROR("Error while parsing B field map - number of values in file not compatible with given number of bins in header");
        m_bFieldOK = false;
        break;
      }
      else {
        m_bField[indexX][indexY][indexZ][0] = bField[0];
        m_bField[indexX][indexY][indexZ][1] = bField[1];
        m_bField[indexX][indexY][indexZ][2] = bField[2];
      }
    }
  }

  //conversion (if units not defined in a file, default cm assumed -> conversion to software default mm done)
  if      (m_dataUnit=="T") m_dataFactor = 1./Units::T;
  else {
    std::ostringstream message;
    message << "BFieldMap: units " << m_dataUnit << " not recognized. Used T instead!";
    logWARNING(message.str());
  }

  double xFactor = 1;
  if      (m_xUnit=="mm") xFactor = Units::mm;
  else if (m_xUnit=="cm") xFactor = Units::cm;
  else if (m_xUnit=="m" ) xFactor = Units::m;
  else {
    std::ostringstream message;
    message << "BFieldMap: X units " << m_xUnit << " not recognized. Used cm instead!";
    logWARNING(message.str());
  }
  m_xMin      *= xFactor;
  m_xMax      *= xFactor;
  m_xBinWidth *= xFactor;

  double yFactor = 1;
  if      (m_yUnit=="mm") yFactor = Units::mm;
  else if (m_yUnit=="cm") yFactor = Units::cm;
  else if (m_yUnit=="m" ) yFactor = Units::m;
  else {
    std::ostringstream message;
    message << "BFieldMap: Y units " << m_yUnit << " not recognized. Used cm instead!";
    logWARNING(message.str());
  }
  m_yMin      *= yFactor;
  m_yMax      *= yFactor;
  m_yBinWidth *= yFactor;

  double zFactor = 1;
  if      (m_zUnit=="mm") zFactor = Units::mm;
  else if (m_zUnit=="cm") zFactor = Units::cm;
  else if (m_zUnit=="m" ) zFactor = Units::m;
  else {
    std::ostringstream message;
    message << "BFieldMap: Z units " << m_zUnit << " not recognized. Used cm instead!";
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
    message << "B field map of type mesh: " << m_fileName << " read in!";
    logINFO(message.str());
  }
  if (m_typeHist) {
    message << "B field map of type histogram: " << m_fileName << " read in!";
    logINFO(message.str());
  }
  if (!m_typeMesh && !m_typeHist) {
    message << "B field map: " << m_fileName << " - binning doesn't correspond to the defined range and the bin size!";
    logERROR(message.str());
  }
}

bool BFieldMap::isInRegionXYZ(double xPos, double yPos, double zPos) const {
  return ((m_xMin <= xPos) && (m_xMax >= xPos) && (m_yMin <= yPos) && (m_yMax >= yPos) && (m_zMin <= zPos) && (m_zMax >= zPos));
}

bool BFieldMap::isInBinRegionXYZ(int& xPos, int& yPos, int& zPos) const {
  return ((xPos>=0) && (xPos<=m_xBinNum) && (yPos>=0) && (yPos<=m_yBinNum) && (zPos>=0) && (zPos<=m_zBinNum));
}

std::vector<double> BFieldMap::calculateBField(double xPos, double yPos, double zPos) const {
  double x     = 0;
  double y     = 0;
  double z     = 0;
  int    x1    = 0;
  int    x2    = 0;
  int    y1    = 0;
  int    y2    = 0;
  int    z1    = 0;
  int    z2    = 0;
  double bAuxField[2][2][2][3];
  std::vector<double> bField;
  for (int i=0; i<3; i++) bField.push_back(0);

  // Check that B field was correctly set
  if (!m_bFieldOK) return bField;

  // Mesh type - refer to mesh points
  if (m_typeMesh) {
    x = (xPos - m_xMin) / m_xBinWidth;
    y = (yPos - m_yMin) / m_yBinWidth;
    z = (zPos - m_zMin) / m_zBinWidth;
  }
  // Histogram type - always refer to the bin center
  if (m_typeHist) {
    x = (xPos - (m_xMin + m_xBinWidth/2.)) / m_xBinWidth;
    y = (yPos - (m_yMin + m_yBinWidth/2.)) / m_yBinWidth;
    z = (zPos - (m_zMin + m_zBinWidth/2.)) / m_zBinWidth;
  }
  // Take the 4 nearest bin centers & interpolate
  x1 = int(floor(x));
  x2 = int(ceil(x));
  y1 = int(floor(y));
  y2 = int(ceil(y));
  z1 = int(floor(z));
  z2 = int(ceil(z));

  // Take care about first and list bins in histogram -> no interpolation
  if (m_typeHist) {

    if ( x>=-1/2. && x<=0 ) x1 = x2;
    if ( y>=-1/2. && y<=0 ) y1 = y2;
    if ( z>=-1/2. && z<=0 ) z1 = z2;

    if ( x>=(m_xBinNum-1) && x<=((m_xBinNum-1)+1/2.) ) x2 = x1;
    if ( y>=(m_yBinNum-1) && y<=((m_yBinNum-1)+1/2.) ) y2 = y1;
    if ( z>=(m_zBinNum-1) && z<=((m_zBinNum-1)+1/2.) ) z2 = z1;
  }

  // Avoid running out of range for max values -> use limit values instead and warn user!
  if (x2>=m_xBinNum) {
    x2 = m_xBinNum-1;
    x1 = m_xBinNum-1;
    std::ostringstream message;
    message << "Value out of range when using: " << m_fileName << "B field calculated for X= " << (x2+1)*(m_xBinWidth)+m_xMin << " instead for X= " << xPos;
    logWARNING(message.str());
    xPos = (x2+1)*(m_xBinWidth)+m_xMin;
  }
  if (y2>=m_yBinNum) {
    y2 = m_yBinNum-1;
    y1 = m_yBinNum-1;
    std::ostringstream message;
    message << "Value out of range when using: " << m_fileName << "B field calculated for Y= " << (y2+1)*(m_yBinWidth)+m_yMin << " instead for Y= " << yPos;
    logWARNING(message.str());
    yPos = (y2+1)*(m_yBinWidth)+m_yMin;
  }
  if (z2>=m_zBinNum) {
    z2 = m_zBinNum-1;
    z1 = m_zBinNum-1;
    std::ostringstream message;
    message << "Value out of range when using: " << m_fileName << "Irradiation calculated for Z= " << (z2+1)*(m_zBinWidth)+m_zMin << " instead for Z= " << zPos;
    logWARNING(message.str());
    zPos = (z2+1)*(m_zBinWidth)+m_zMin;
  }

  if (isInRegionXYZ(xPos, yPos, zPos)) {
    if (isInBinRegionXYZ(x1,y1,z1) && isInBinRegionXYZ(x2,y2,z2)) {

      //if the point is in a intersection of the grid formed by the map bin centers
      if((x1==x2) && (y1==y2) && (z1==z2)) {
        //single value
        for (int i=0; i<3; i++) bField[i] = m_bField[x1][y1][z1][i];
      }
      //interpolate
      else {
        for (int i=0; i<3; i++) bAuxField[1][1][1][i] = m_bField[x1][y1][z1][i];
        for (int i=0; i<3; i++) bAuxField[2][1][1][i] = m_bField[x2][y1][z1][i];
        for (int i=0; i<3; i++) bAuxField[1][2][1][i] = m_bField[x1][y2][z1][i];
        for (int i=0; i<3; i++) bAuxField[1][1][2][i] = m_bField[x1][y1][z2][i];
        for (int i=0; i<3; i++) bAuxField[2][2][1][i] = m_bField[x2][y2][z1][i];
        for (int i=0; i<3; i++) bAuxField[2][1][2][i] = m_bField[x2][y1][z2][i];
        for (int i=0; i<3; i++) bAuxField[1][2][2][i] = m_bField[x1][y2][z2][i];
        for (int i=0; i<3; i++) bAuxField[2][2][2][i] = m_bField[x2][y2][z2][i];
        //bilinear interpolation in x, y and z
        double volume = (x2-x1)*(y2-y1)*(x2-x1);
        for (int i=0; i<3; i++) bField[i] = bAuxField[1][1][1][i]/volume*(x2-x)*(y2-y)*(z2-z) +
                                            bAuxField[2][1][1][i]/volume*(x-x1)*(y2-y)*(z2-z) +
                                            bAuxField[1][2][1][i]/volume*(x2-x)*(y-y1)*(z2-z) +
                                            bAuxField[1][1][2][i]/volume*(x2-x)*(y2-y)*(z-z1) +
                                            bAuxField[2][2][1][i]/volume*(x-x1)*(y-y1)*(z2-z) +
                                            bAuxField[2][1][2][i]/volume*(x-x1)*(y2-y)*(z1-z) +
                                            bAuxField[1][2][2][i]/volume*(x2-x)*(y-y1)*(z-z1) +
                                            bAuxField[2][2][2][i]/volume*(x-x1)*(y-y1)*(z-z1);
      }
    }
    else logERROR("B Field map: (X,Y,Z) bin out of region");
  }
  else logERROR("B Field map: (X,Y,Z) position out of region");

  // Convert output to
  for (int i=0; i<3; i++) bField[i] *= m_dataFactor;
  return bField;
}

bool BFieldMap::drawXZBFieldProj(TCanvas& xzCanvas, std::string name, double minX, double maxX, double minZ, double maxZ)
{

  if (m_bFieldOK) {

    xzCanvas.cd();

    maxX = ceil((maxX - m_xMin)/m_xBinWidth)*m_xBinWidth + m_xMin;
    maxZ = ceil((maxZ - m_zMin)/m_zBinWidth)*m_zBinWidth + m_zMin;

    minX = floor((minX - m_xMin)/m_xBinWidth)*m_xBinWidth + m_xMin;
    minZ = floor((minZ - m_zMin)/m_zBinWidth)*m_zBinWidth + m_zMin;

    double Y = floor((0 - m_yMin)/m_yBinWidth)*m_yBinWidth + m_yMin;

    // Check required range and adapt
    if (minX<m_xMin) {
      minX = m_xMin;
      logWARNING("Drawing B field map in maximum available X range, the required range is too big");
    }
    if (minZ<m_zMin) {
      minZ = m_zMin;
      logWARNING("Drawing B field map in maximum available Z range, the required range is too big");
    }
    if (maxX>m_xMax) {
      maxX = m_xMax;
      logWARNING("Drawing B field map in maximum available X range, the required range is too big");
    }
    if (maxZ>m_zMax) {
      maxZ = m_zMax;
      logWARNING("Drawing B field map in maximum available Z range, the required range is too big");
    }
    if (Y<m_yMin || Y>m_yMax) {
      logERROR("Can't draw B field map (Y=0): Y=0 not provided in B field map range");
      return false;
    }

    int nBinsX = int((maxX - minX)/m_xBinWidth);
    int nBinsZ = int((maxZ - minZ)/m_zBinWidth);
    int nBinsXOffset = int(minX-m_xMin/m_xBinWidth);
    int nBinsZOffset = int(minZ-m_zMin/m_xBinWidth);
    int yBin = (Y - m_yMin)/m_yBinWidth;

    // Draw histogram
    TH2D* his = new TH2D(name.c_str(), std::string("XZ view of B field ["+m_dataUnit+"] (Y=0)").c_str(), nBinsZ, 0, maxZ, nBinsX, 0, maxX);
    his->SetStats(kFALSE);
    his->Draw("COLZ");

    // Get min & max values
    double minBFieldXZ = +std::numeric_limits<double>::max();
    double maxBFieldXZ = -std::numeric_limits<double>::max();
    for (int xBin=0; xBin<nBinsX; xBin++) {
      for (int zBin=0; zBin<nBinsZ; zBin++) {

        double XBField = m_bField[xBin+nBinsXOffset][yBin][zBin+nBinsZOffset][0];
        double ZBField = m_bField[xBin+nBinsXOffset][yBin][zBin+nBinsZOffset][2];
        double vecMag  = sqrt(XBField*XBField+ZBField*ZBField);
        if (vecMag<minBFieldXZ) minBFieldXZ = vecMag;
        if (vecMag>maxBFieldXZ) maxBFieldXZ = vecMag;
      }
    }

    // Fill histogram
    for (int xBin=0; xBin<nBinsX; xBin++) {
      for (int zBin=0; zBin<nBinsZ; zBin++) {

        double bFieldMag   = 0;
        for (int i=0; i<3; i++) bFieldMag   += m_bField[xBin+nBinsXOffset][yBin][zBin+nBinsZOffset][i]*m_bField[xBin+nBinsXOffset][yBin][zBin+nBinsZOffset][i];
        bFieldMag   = sqrt(bFieldMag);

        his->SetBinContent(zBin+1, xBin+1, bFieldMag);

        // Draw B field direction
        double XBField = m_bField[xBin+nBinsXOffset][yBin][zBin+nBinsZOffset][0];
        double ZBField = m_bField[xBin+nBinsXOffset][yBin][zBin+nBinsZOffset][2];
        double vecMag  = sqrt(XBField*XBField+ZBField*ZBField);
        XBField        = XBField/vecMag * ((c_arrowMaxLength-c_arrowMinLength) * m_xBinWidth * (vecMag-minBFieldXZ)/(maxBFieldXZ-minBFieldXZ) + c_arrowMinLength*m_yBinWidth);
        ZBField        = ZBField/vecMag * ((c_arrowMaxLength-c_arrowMinLength) * m_zBinWidth * (vecMag-minBFieldXZ)/(maxBFieldXZ-minBFieldXZ) + c_arrowMinLength*m_yBinWidth);

        double xStart  = (xBin+0.5)*m_xBinWidth + minX - XBField/2.;
        double zStart  = (zBin+0.5)*m_zBinWidth + minZ - ZBField/2.;

        TArrow* direction = new TArrow(zStart, xStart, zStart+ZBField, xStart+XBField);
        direction->Draw("|>");
        direction->SetArrowSize((c_arrowMaxSize-c_arrowMinSize)*(vecMag-minBFieldXZ)/(maxBFieldXZ-minBFieldXZ) + c_arrowMinSize);
        direction->SetAngle(45);

        if ((vecMag-minBFieldXZ)/(maxBFieldXZ-minBFieldXZ)>0.75) {
          direction->SetLineColor(kBlack);
          direction->SetFillColor(kBlack);
        }
        else if ((vecMag-minBFieldXZ)/(maxBFieldXZ-minBFieldXZ)>0.5) {
          direction->SetLineColor(kGray+3);
          direction->SetFillColor(kGray+3);
        }
        else if ((vecMag-minBFieldXZ)/(maxBFieldXZ-minBFieldXZ)>0.25) {
          direction->SetLineColor(kGray+2);
          direction->SetFillColor(kGray+2);
        }
        else {
          direction->SetLineColor(kGray+1);
          direction->SetFillColor(kGray+1);
        }
      }
    }
    his->GetXaxis()->SetTitle(std::string("Z [mm]").c_str());
    his->GetXaxis()->SetTitleOffset(1.2);
    his->GetYaxis()->SetTitle(std::string("X [mm]").c_str());
    his->GetYaxis()->SetTitleOffset(1.2);

    return true;
  }
  else return false;
}

bool BFieldMap::drawYZBFieldProj(TCanvas& yzCanvas, std::string name, double minY, double maxY, double minZ, double maxZ)
{
  if (m_bFieldOK) {

    yzCanvas.cd();

    maxY = ceil((maxY - m_yMin)/m_yBinWidth)*m_yBinWidth + m_yMin;
    maxZ = ceil((maxZ - m_zMin)/m_zBinWidth)*m_zBinWidth + m_zMin;

    minY = floor((minY - m_yMin)/m_yBinWidth)*m_yBinWidth + m_yMin;
    minZ = floor((minZ - m_zMin)/m_zBinWidth)*m_zBinWidth + m_zMin;

    double X = floor((0 - m_xMin)/m_xBinWidth)*m_xBinWidth + m_xMin;

    // Check required range and adapt
    if (minY<m_yMin) {
      minY = m_yMin;
      logWARNING("Drawing B field map in maximum available Y range, the required range is too big");
    }
    if (minZ<m_zMin) {
      minZ = m_zMin;
      logWARNING("Drawing B field map in maximum available Z range, the required range is too big");
    }
    if (maxY>m_yMax) {
      maxY = m_yMax;
      logWARNING("Drawing B field map in maximum available Y range, the required range is too big");
    }
    if (maxZ>m_zMax) {
      maxZ = m_zMax;
      logWARNING("Drawing B field map in maximum available Z range, the required range is too big");
    }

    int nBinsY = int((maxY - minY)/m_yBinWidth);
    int nBinsZ = int((maxZ - minZ)/m_zBinWidth);
    int nBinsYOffset = int(minY-m_yMin/m_yBinWidth);
    int nBinsZOffset = int(minZ-m_zMin/m_xBinWidth);
    int xBin = (X - m_xMin)/m_xBinWidth;

    // Draw histogram
    TH2D* his = new TH2D(name.c_str(), std::string("YZ view of B field ["+m_dataUnit+"] (X=0)").c_str(), nBinsZ, 0, maxZ, nBinsY, 0, maxY);
    his->SetStats(kFALSE);
    his->Draw("COLZ");

    // Get min & max values
    double minBFieldYZ = +std::numeric_limits<double>::max();
    double maxBFieldYZ = -std::numeric_limits<double>::max();
    for (int yBin=0; yBin<nBinsY; yBin++) {
      for (int zBin=0; zBin<nBinsZ; zBin++) {

        double YBField = m_bField[xBin][yBin+nBinsYOffset][zBin+nBinsZOffset][1];
        double ZBField = m_bField[xBin][yBin+nBinsYOffset][zBin+nBinsZOffset][2];
        double vecMag  = sqrt(YBField*YBField+ZBField*ZBField);
        if (vecMag<minBFieldYZ) minBFieldYZ = vecMag;
        if (vecMag>maxBFieldYZ) maxBFieldYZ = vecMag;
      }
    }

    // Fill histogram
    for (int yBin=0; yBin<nBinsY; yBin++) {
      for (int zBin=0; zBin<nBinsZ; zBin++) {

        double bFieldMag   = 0;
        for (int i=0; i<3; i++) bFieldMag   += m_bField[xBin][yBin+nBinsYOffset][zBin+nBinsZOffset][i]*m_bField[xBin][yBin+nBinsYOffset][zBin+nBinsZOffset][i];
        bFieldMag   = sqrt(bFieldMag);

        his->SetBinContent(zBin+1, yBin+1, bFieldMag);

        // Draw B field direction
        double YBField = m_bField[xBin][yBin+nBinsYOffset][zBin+nBinsZOffset][1];
        double ZBField = m_bField[xBin][yBin+nBinsYOffset][zBin+nBinsZOffset][2];
        double vecMag  = sqrt(YBField*YBField+ZBField*ZBField);
        YBField        = YBField/vecMag * ((c_arrowMaxLength-c_arrowMinLength) * m_yBinWidth * (vecMag-minBFieldYZ)/(maxBFieldYZ-minBFieldYZ) + c_arrowMinLength*m_yBinWidth);
        ZBField        = ZBField/vecMag * ((c_arrowMaxLength-c_arrowMinLength) * m_zBinWidth * (vecMag-minBFieldYZ)/(maxBFieldYZ-minBFieldYZ) + c_arrowMinLength*m_zBinWidth);

        double yStart  = (yBin+0.5)*m_yBinWidth + minY - YBField/2.;
        double zStart  = (zBin+0.5)*m_zBinWidth + minZ - ZBField/2.;

        TArrow* direction = new TArrow(zStart, yStart, zStart+ZBField, yStart+YBField);
        direction->Draw("|>");
        direction->SetArrowSize((c_arrowMaxSize-c_arrowMinSize)*(vecMag-minBFieldYZ)/(maxBFieldYZ-minBFieldYZ) + c_arrowMinSize);
        direction->SetAngle(45);

        if ((vecMag-minBFieldYZ)/(maxBFieldYZ-minBFieldYZ)>0.75) {
          direction->SetLineColor(kBlack);
          direction->SetFillColor(kBlack);
        }
        else if ((vecMag-minBFieldYZ)/(maxBFieldYZ-minBFieldYZ)>0.5) {
          direction->SetLineColor(kGray+3);
          direction->SetFillColor(kGray+3);
        }
        else if ((vecMag-minBFieldYZ)/(maxBFieldYZ-minBFieldYZ)>0.25) {
          direction->SetLineColor(kGray+2);
          direction->SetFillColor(kGray+2);
        }
        else {
          direction->SetLineColor(kGray+1);
          direction->SetFillColor(kGray+1);
        }
      }
    }
    his->GetXaxis()->SetTitle(std::string("Z [mm]").c_str());
    his->GetXaxis()->SetTitleOffset(1.2);
    his->GetYaxis()->SetTitle(std::string("Y [mm]").c_str());
    his->GetYaxis()->SetTitleOffset(1.2);
    return true;
  }
  else return false;
}

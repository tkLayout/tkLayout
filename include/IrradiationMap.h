/*
 * IrradiationMap.h
 *
 *  Created on: 18/feb/2014
 *      Author: Stefano Martina
 */

#ifndef IRRADIATIONMAP_H_
#define IRRADIATIONMAP_H_


#include <string>
#include <fstream>
#include <utility>
#include <vector>
#include <cmath>
#include "messageLogger.h"


class IrradiationMap {
public:
  IrradiationMap(std::string irradiationMapFile);
  IrradiationMap();

  /**
   * Populate the map attributes reading the passed file
   * @param irradiationMapFile is the path of the raw file for the irradiation map to be read
   */
  void ingest(std::string irradiationMapFile);

  std::pair<std::pair<double,double>,std::pair<double,double>> region() const;

  std::pair<double,double> binDimension() const;

  /**
   * Get the area of a bin of the map, identifies the resolution of the map
   */
  double binArea() const;

  /**
   * Overriding of the operator < for allowing the sort of a vector of maps by their resolutions
   */
  bool operator < (const IrradiationMap& confrontedMap) const;

  /**
   * Test if a point is inside the area covered by the map
   * @param coordinates is a pair (z,rho) that indicate a point in the plane ZxRho
   */
  bool isInRegion(std::pair<double,double> coordinates) const;

  /**
   * Get the irradiation of the point
   * @param coordinates is a (z,rho) that indicate a point in the plane ZxRho
   * @return the value of the irradiation of the point, or 0 if the point is outside the map
   */
  double calculateIrradiation(std::pair<double,double> coordinates) const;

private:
  const std::string comp_rhoMin = "# R min: ";
  const std::string comp_rhoMax = "# R max: ";
  const std::string comp_rhoBinWidth = "# R bin width: ";
  const std::string comp_rhoBinNum = "# R number of bins: ";
  const std::string comp_zMin = "# Z min: ";
  const std::string comp_zMax = "# Z max: ";
  const std::string comp_zBinWidth = "# Z bin width: ";
  const std::string comp_zBinNum = "# Z number of bins: ";
  const std::string comp_invFemUnit = "# normalization value: ";

  double rhoMin;
  double rhoMax;
  double rhoBinWidth;
  long int rhoBinNum;
  double zMin;
  double zMax;
  double zBinWidth;
  long int zBinNum;
  double invFemUnit;

  //matrix rho X z of the irradiation values
  std::vector< std::vector<double> > irradiation;

  void inizialize();
};


#endif /* IRRADIATIONMAP_H_ */

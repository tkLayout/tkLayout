/**
 * @file IrradiationMap.h
 * @author Stefano Martina
 * @date 18/feb/2014
 */

#ifndef IRRADIATIONMAP_H_
#define IRRADIATIONMAP_H_


#include <string>
#include <fstream>
#include <utility>
#include <vector>
#include <cmath>
#include "MessageLogger.hh"

/**
 * @class IrradiationMap
 * @brief This class represent a single irradiation map.
 * @details It is possible to feed the map with a new file, it can read the values from header.
 * The maps are sortable on resolution with operator <.
 */
class IrradiationMap {
public:
  /**
   * Constructor with feeding
   * @param irradiationMapFile is the path of the new file to feed
   */
  IrradiationMap(std::string irradiationMapFile);

  /**
   * Constructor without feeding
   */
  IrradiationMap();

  /**
   * Populate the map attributes reading the passed file
   * @param irradiationMapFile is the path of the raw file for the irradiation map to be read
   */
  void ingest(std::string irradiationMapFile);

  /**
   * Get the area of a bin of the map, identifies the resolution of the map
   * @return The area of the bin
   */
  double binArea() const;

  /**
   * Overriding of the operator < for allowing the sort of a vector of maps by their resolutions
   */
  bool operator < (const IrradiationMap& confrontedMap) const;

  /**
   * Test if a point is inside the area covered by the map
   * @param coordinates is a pair (z,rho) that indicate a point in the plane ZxRho
   * @return True if the point is inside the map region, false otherwise
   */
  bool isInRegion(const std::pair<double,double>& coordinates) const;

  /**
   * Get the irradiation of the point
   * @param coordinates is a (z,rho) that indicate a point in the plane ZxRho
   * @return the value of the irradiation of the point, or 0 if the point is outside the map
   */
  double calculateIrradiation(const std::pair<double,double>& coordinates) const;

private:
  const std::string comp_rhoMin = "# R min: ";                  /**< Prefix of the line of the header of the feeded file that precedes the value of min rho*/
  const std::string comp_rhoMax = "# R max: ";                  /**< Prefix of the line of the header of the feeded file that precedes the value of max rho*/
  const std::string comp_rhoBinWidth = "# R bin width: ";       /**< Prefix of the line of the header of the feeded file that precedes the value of bin width in rho*/
  const std::string comp_rhoBinNum = "# R number of bins: ";    /**< Prefix of the line of the header of the feeded file that precedes the value of the number of bins in rho*/
  const std::string comp_zMin = "# Z min: ";                    /**< Prefix of the line of the header of the feeded file that precedes the value of min Z*/
  const std::string comp_zMax = "# Z max: ";                    /**< Prefix of the line of the header of the feeded file that precedes the value of max Z*/
  const std::string comp_zBinWidth = "# Z bin width: ";         /**< Prefix of the line of the header of the feeded file that precedes the value of bin width in Z*/
  const std::string comp_zBinNum = "# Z number of bins: ";      /**< Prefix of the line of the header of the feeded file that precedes the value of the number of bins in Z*/
  const std::string comp_invFemUnit = "# normalization value: ";/**< Prefix of the line of the header of the feeded file that precedes the value of the normalization value in fb^-1*/

  double rhoMin;        /**< The value of min rho*/
  double rhoMax;        /**< The value of max rho*/
  double rhoBinWidth;   /**< The value of bin width in rho*/
  long int rhoBinNum;   /**< The value of the number of bins in rho*/
  double zMin;          /**< The value of min Z*/
  double zMax;          /**< The value of max Z*/
  double zBinWidth;     /**< The value of bin width in Z*/
  long int zBinNum;     /**< The value of the number of bins in Z*/
  double invFemUnit;    /**< The value of the normalization value in fb^-1*/

  std::vector< std::vector<double> > irradiation;       /**< The matrix (rho * Z) that contains the irradiation values for each bin of the map*/
};


#endif /* IRRADIATIONMAP_H_ */

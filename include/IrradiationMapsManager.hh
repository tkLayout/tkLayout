/**
 * @file IrradiationMapsManager.h
 * @author Stefano Martina
 * @date 19/feb/2014
 */

#ifndef IRRADIATIONMAPSMANAGER_H_
#define IRRADIATIONMAPSMANAGER_H_

#include <utility>
#include <string>
#include <set>
#include"IrradiationMap.hh"

/**
 * @class IrradiationMapsManager
 * @brief The administrator of the irradiation maps.
 * @details Mantains a set of maps sorted by resolution, when
 * is asked for the irradiation of a point returns the value of the
 * better map that contains this point in his region
 */
class IrradiationMapsManager {
public:
  IrradiationMapsManager();
  ~IrradiationMapsManager();

  /**
   * Add the passed map to internal set
   * @param IrradiationMap is the new map
   */
  void addIrradiationMap(const IrradiationMap& newIrradiationMap);

  /**
   * Create a new map with passed file and ad it to internal set
   * @param newIrradiationMapFile is the path of the file for the new map
   */
  void addIrradiationMap(std::string newIrradiationMapFile);

  /**
   * Get the irradiation of the point in the passed coordinates using the best
   * avaiable map that contains that point
   * @param coordinates represent the point, is a pair (z,r) with coordinates z in Z, and r in Rho
   * @return The value of irradiation in the point
   */
  double calculateIrradiationPower(const std::pair<double,double>& coordinates) const;

private:

  /**
   * The set that contains all the maps ordered by resolution
   */
  std::set<IrradiationMap> irradiationMaps;
};

#endif /* IRRADIATIONMAPSMANAGER_H_ */

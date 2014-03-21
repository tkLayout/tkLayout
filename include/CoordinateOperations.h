/*
 * CordinatesOperations.h
 *
 *  Created on: 20/mar/2014
 *      Author: stefano
 */

#ifndef CORDINATESOPERATIONS_H_
#define CORDINATESOPERATIONS_H_

#include <Math/Vector3Dfwd.h>
#include <vector>
#include "Polygon3d.h"
#include "global_funcs.h"

using ROOT::Math::XYZVector;

namespace CoordinateOperations{

  template<class Polygon> std::vector<XYZVector> computeDistanceVectors(const Polygon& polygon) {
    std::vector<XYZVector> distanceVectors;
    auto vertex0 = polygon.begin();
    for(auto vertex1 = polygon.begin()+1; vertex1 != polygon.end(); vertex0 = vertex1++) {
      XYZVector directionVector = (*vertex0 - *vertex1).Unit();
      if (vertex0->Dot(directionVector) >= 0.0) {
        distanceVectors.push_back(*vertex0);
      }else if (vertex1->Dot(directionVector) <= 0.0) {
        distanceVectors.push_back(*vertex1);
      } else {

        XYZVector appo = *vertex0 - (vertex0->Dot(directionVector) * directionVector);

        distanceVectors.push_back(appo);
      }
    }
    return distanceVectors;
  }

  template<class Polygon> double computeMinZ(const Polygon& polygon) {
    return gmin<double>(polygon.begin(), polygon.end(), [](const XYZVector& v) { return v.Z(); });
  }

  template<class Polygon> double computeMaxZ(const Polygon& polygon) {
    return gmax<double>(polygon.begin(), polygon.end(), [](const XYZVector& v) { return v.Z(); });
  }

  template<class Polygon> double computeMinR(const Polygon& polygon) {
    auto distanceVectors = computeDistanceVectors(polygon);
    return gmin<double>(distanceVectors.begin(), distanceVectors.end(), [](const XYZVector& v) { return v.Rho(); });
  }

  template<class Polygon> double computeMaxR(const Polygon& polygon) {
    auto distanceVectors = computeDistanceVectors(polygon);
    return gmax<double>(distanceVectors.begin(), distanceVectors.end(), [](const XYZVector& v) { return v.Rho(); });
  }

}




#endif /* CORDINATESOPERATIONS_H_ */

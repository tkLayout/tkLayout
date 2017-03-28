/*
 * coordinatesOperations.cpp
 *
 *  Created on: 20/mar/2014
 *      Author: stefano
 */

#include "CoordinateOperations.hh"

namespace CoordinateOperations {

  XYZVector computeDistanceVector(const XYZVector& v0, const XYZVector& v1) {
    XYZVector dirv = (v0 - v1).Unit(); // direction vector of the line passing between v0 and v1
    double projv0 = v0.Dot(dirv); // scalar projection of v0 on the direction vector
    double projv1 = v1.Dot(dirv); // scalar projection of v1 on the direction vector
    if (projv0 <= 0. && projv1 < 0.) {
      return v0;
    } else if (projv0 > 0. && projv1 >= 0.) {
      return v1;
    } else {
      return v0 - (v0.Dot(dirv) * dirv); 
    }
  }

}

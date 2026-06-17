/*
 * CordinatesOperations.h
 *
 *  Created on: 20/mar/2014
 *      Author: stefano
 */

#ifndef CORDINATEOPERATIONS_H
#define CORDINATEOPERATIONS_H

#include <vector>
#include <array>
#include <cstddef>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/VectorUtil.h>
#include <TVector3.h>
#include "Polygon3d.hh"
#include "global_funcs.hh"

using ROOT::Math::XYZVector;


namespace CoordinateOperations {

  /*
   * Given a box center and its three half-axis vectors (already scaled to half-extent),
   * returns the 16 bounding-box candidate points: top & bottom corners plus the top &
   * bottom midpoints of each edge. Used for module-with-hybrids extrema computation.
   * Corner ordering (in the xy plane): 0=(+x,+y), 1=(-x,+y), 2=(-x,-y), 3=(+x,-y).
   */
  inline std::array<XYZVector, 16> boundingBoxCandidates(const XYZVector& center,
                                                         const XYZVector& halfX,
                                                         const XYZVector& halfY,
                                                         const XYZVector& halfZ) {
    std::array<XYZVector, 4> xyCorners;
    for (std::size_t i = 0; i < 4; ++i) {
      double signX = (i == 0 || i == 3) ? 1. : -1.;
      double signY = (i == 0 || i == 1) ? 1. : -1.;
      xyCorners[i] = center + signX * halfX + signY * halfY;
    }
    std::array<XYZVector, 16> pts;
    for (std::size_t i = 0; i < 4; ++i) {
      std::size_t next = (i + 1) % 4;
      XYZVector top = xyCorners[i] + halfZ;
      XYZVector bot = xyCorners[i] - halfZ;
      pts[i*4]     = top;
      pts[i*4 + 1] = bot;
      pts[i*4 + 2] = (top + (xyCorners[next] + halfZ)) * 0.5;
      pts[i*4 + 3] = (bot + (xyCorners[next] - halfZ)) * 0.5;
    }
    return pts;
  }

  /*
   * Convert XYZVector or Polar3DVector into a TVector. 
   * TVector is much more generic, and allows Angle() operation for example, which are not permitted with the old coord classes.
   */
  template<class CoordVector> const TVector3 convertCoordVectorToTVector3(const CoordVector& v) {
    return TVector3(v.X(), v.Y(), v.Z());
  }


  /*
   * Returns projection of v1 on the plane which has v2 as a normal.
   * WARNING!! Here, v2 needs to be unitary: ||v2|| = 1.
   * || || used is the euclidian norm on R3.
   * Dot is the scalar product associated to || ||.
   * NB: Very weird that this is not provided by default, v1.Perp(v2) provides a scalar but not a vector.
   */
  template<class GeoVector> GeoVector projectv1OnPlaneOfNormalUnitv2(const GeoVector& v1, const GeoVector& v2) {
    return (v1 - v1.Dot(v2) * v2);
  }


  XYZVector computeDistanceVector(const XYZVector& v0, const XYZVector& v1); // minimum distance vector (from the origin) of the segment defined by v0 and v1

  template<class Polygon> std::vector<XYZVector> computeDistanceVectors(const Polygon& polygon) {
    std::vector<XYZVector> distanceVectors;
    XYZVector v0 = polygon.getVertex(0);
    for (int i = 1; i < polygon.getNumSides() + 1; i++) {
      XYZVector v1 = polygon.getVertex(i % polygon.getNumSides());
      v0.SetZ(0.0);
      v1.SetZ(0.0);
      distanceVectors.push_back(computeDistanceVector(v0, v1));
      v0 = v1;
    }
    return distanceVectors;
  }

  template<class Polygon> double computeMinZ(const Polygon& polygon) {
    return minget(polygon.begin(), polygon.end(), [](const XYZVector& v) { return v.Z(); });
  }

  template<class Polygon> double computeMaxZ(const Polygon& polygon) {
    return maxget(polygon.begin(), polygon.end(), [](const XYZVector& v) { return v.Z(); });
  }

  template<class Polygon> double computeMinR(const Polygon& polygon) {
    auto distanceVectors = computeDistanceVectors(polygon);
    return minget(distanceVectors.begin(), distanceVectors.end(), [](const XYZVector& v) { return v.Rho(); });
  }

  template<class Polygon> double computeMaxR(const Polygon& polygon) {
    return maxget(polygon.begin(), polygon.end(), [](const XYZVector& v) { return v.Rho(); });
  }

  /**
   * Compute the polygon whose vertices are all the middles of the edges of the polygon sepcified as a parameter.
   */
  template<class Polygon> Polygon* computeMidPolygon(const Polygon& polygon) {
    Polygon* midPoly = new Polygon();

    XYZVector v0 = polygon.getVertex(0);
    for (int i = 1; i < polygon.getNumSides() + 1; i++) {
      XYZVector v1 = polygon.getVertex(i % polygon.getNumSides());
      XYZVector vMid = (v0 + v1) / 2.;
      *midPoly << vMid;
      v0 = v1;   
    }
    return midPoly;
  }


/**
 * Compute the envelope polygon (2n vertices) from a polygon (n vertices) specified as a parameter.
 * param basePolygon
 * param normalOffset
 * return envelopePolygon : poly formed by basePolygon shifted by a - normal offset, and by basePolygon shifted by a + normal offset.
 */
  template<class PolygonA, class PolygonB> PolygonB* computeEnvelopePolygon(const PolygonA& basePolygon, double normalOffset) {
    double innerOffset =  -normalOffset;
    PolygonA* innerPoly = new PolygonA(basePolygon);
    innerPoly->translate(innerPoly->getNormal() * innerOffset);

    double outerOffset =  normalOffset;
    PolygonA* outerPoly = new PolygonA(basePolygon);
    outerPoly->translate(outerPoly->getNormal() * outerOffset);

    PolygonB* envelopePoly = new PolygonB();
    for (int i = 0; i < innerPoly->getNumSides(); i++) *envelopePoly << innerPoly->getVertex(i);
    for (int i = 0; i < outerPoly->getNumSides(); i++) *envelopePoly << outerPoly->getVertex(i);
    // Would have been much nicer to do sth like : *envelopePoly << innerPoly->getVertices() , 
    // unfortunately this apparently does not work and only assign the first two vertices to *envelopePoly.
    // On a more general note, the entire AbstractPolygon and Polygon3d classes should be rewritten.
 
    return envelopePoly;
  }

}




#endif /* CORDINATEOPERATIONS_H */

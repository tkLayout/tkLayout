#ifndef INCLUDE_POLYGON3D_H_
#define INCLUDE_POLYGON3D_H_

#include <vector>
#include <set>
#include <functional>
#include <algorithm>
#include <numeric>

#include <TRandom.h>
#include "math_functions.h"
#include <Math/Vector3Dfwd.h>
#include <Math/GenVector/RotationZ.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/Transform3D.h>
#include "MessageLogger.h"
#include "StringConverter.h"

// Using from namespace
using ROOT::Math::XYZVector;

template<class Polygon> struct PolygonLess {
  bool operator()(const Polygon& t1, const Polygon& t2) const {
    return t1.getDoubleArea() < t2.getDoubleArea();
  }
};

//!
//! Define rotation X
//!
template<class Coords, class FloatType> struct RotateX {};
template<class FloatType> struct RotateX<XYZVector, FloatType> : public std::binary_function<XYZVector, FloatType, XYZVector> {
  XYZVector operator()(const XYZVector& vector, FloatType angle) const {
    return XYZVector(vector.X(),
                     vector.Y()*cos(angle) - vector.Z()*sin(angle),
                     vector.Y()*sin(angle) + vector.Z()*cos(angle));
  }
};

//!
//! Define rotation Y
//!
template<class Coords, class FloatType> struct RotateY {};
template<class FloatType> struct RotateY<XYZVector, FloatType> : public std::binary_function<XYZVector, FloatType, XYZVector> {
  XYZVector operator()(const XYZVector& vector, FloatType angle) const {
    return XYZVector(vector.Z()*sin(angle) + vector.X()*cos(angle),
                     vector.Y(),
                     vector.Z()*cos(angle) - vector.X()*sin(angle));
  }
};

//!
//! Define rotation Z
//!
template<class Coords, class FloatType> struct RotateZ {};
template<class FloatType> struct RotateZ<XYZVector, FloatType> : public std::binary_function<XYZVector, FloatType, XYZVector> {
  XYZVector operator()(const XYZVector& vector, FloatType angle) const {
    return XYZVector(vector.X()*cos(angle) - vector.Y()*sin(angle),
                     vector.X()*sin(angle) + vector.Y()*cos(angle),
                     vector.Z());
  }
};

//!
//! Define cross-product
//!
template<class Coords> Coords crossProduct(const Coords& v1, const Coords& v2);

//!
//! Define unity vector
//!
template<class Coords> Coords unitVector(const Coords& vector);


//!
//! An abstract class defining polygon in 3 dimensions with any number of sides, any number of dimensions, convex or concave
//!
template<int NumSides, class Coords, class Random, class FloatType = double> class AbstractPolygon {

 public:
  AbstractPolygon() : m_numOfAllocVertices(0), m_recalCentre(true), m_recalNormal(true) {}
  AbstractPolygon(const Coords& vertex): m_numOfAllocVertices(0), m_recalCentre(true), m_recalNormal(true) { *this << vertex; }
  virtual ~AbstractPolygon() {};

  //! Fill operators
  virtual AbstractPolygon<NumSides, Coords, Random, FloatType>& operator<<(const Coords& vertex);
  virtual AbstractPolygon<NumSides, Coords, Random, FloatType>& operator<<(const std::vector<Coords>& vertices);
  virtual AbstractPolygon<NumSides, Coords, Random, FloatType>& operator()(const Coords& vertex) { *this << vertex; return *this; }

  const Coords* getVertices() const { return m_vertices;}
  const Coords* begin()       const { return m_vertices; }
  const Coords* end()         const { return (m_vertices+NumSides); }

  //!< Get given vertex
  const Coords& getVertex(int index) const;

  //! Get number of required vertices/sides (3 for triangle, 4 for tetrahedron, etc.)
  int  getNumSides() const { return NumSides; };

  //! Does polygon contain all vertices as requested (3 for triangle, 4 for tetrahedron, etc.)
  inline bool isComplete() const {return m_numOfAllocVertices == NumSides;}

  //! Get polygon area multiplied by 2 (one saves a division by 2 by returning the double area -> for our purposes, i.e. comparison, it's fine)
  double getDoubleArea() const {return m_area; } // we save a division by 2 by returning the double area. for our purposes it's the same

  //! Check if point p (defined by 3D vector) is inside a polygon -> use sum of triangle surfaces with one vertex being the point p. If sum equals to the polygon surface, point p is inside!
  virtual bool isPointInside(const Coords& p) const = 0;

  //  virtual Coords generateRandomPoint(Random* die) const = 0; // Wrongly implemented - TODO: needs check

  //! Get polygon center as 3D vector
  const Coords& getCenter() const;

  //! Get polygon normal as 3D vector
  const Coords& getNormal() const; 

  //! Shift polygon by given 3D vector
  AbstractPolygon<NumSides, Coords, Random, FloatType>& translate(const Coords& vector); 

  //! Scalar multiplication of all vertices vectors by given vector, i.e. mirror would by multiplication by (1,1,-1)
  AbstractPolygon<NumSides, Coords, Random, FloatType>& mirror(const Coords& vector); 

  AbstractPolygon<NumSides, Coords, Random, FloatType>& rotateX(FloatType angle);
  AbstractPolygon<NumSides, Coords, Random, FloatType>& rotateY(FloatType angle);
  AbstractPolygon<NumSides, Coords, Random, FloatType>& rotateZ(FloatType angle);

 protected:

  //!< Compute basic properties of the polygon
  virtual void computeProperties() = 0;

  //!< Recalculate centre position & normal
  void recalGeomProperties() { m_recalCentre = m_recalNormal = true; }

  FloatType m_area;            //!< Polygon area
  Coords m_vertices[NumSides]; //!< Vertices of the polygon
  size_t m_numOfAllocVertices; //!< Number of allocated vertices
  mutable Coords m_center;     //!< Polygon centre position (3D vector)
  mutable Coords m_normal;     //!< Normal vector (3D vector)

  mutable bool m_recalCentre;  //!< Needs to recalculate centre?
  mutable bool m_recalNormal;  //!< Needs to recalculate normal?

}; // Class


template<int NumSides, class Coords, class Random, class FloatType> inline const Coords& AbstractPolygon<NumSides, Coords, Random, FloatType>::getVertex(int index) const {
  try {
    return m_vertices[index];
  } catch (const std::out_of_range& ex) {
    logERROR("AbstractPolygon: number of sides out of index "+ any2str(index) + " required!");
    throw;
  }
}

template<int NumSides, class Coords, class Random, class FloatType> const Coords& AbstractPolygon<NumSides, Coords, Random, FloatType>::getCenter() const {
  if (m_recalCentre) {
    m_center      = std::accumulate(&m_vertices[0], &m_vertices[NumSides], XYZVector())/NumSides;
    m_recalCentre = false;
  }
  return m_center;
}

template<int NumSides, class Coords, class Random, class FloatType> const Coords& AbstractPolygon<NumSides, Coords, Random, FloatType>::getNormal() const {
  if (m_recalNormal) {
    m_normal      = unitVector(crossProduct(m_vertices[1]-m_vertices[0], m_vertices[2]-m_vertices[0]));
    m_recalNormal = false;
  }
  return m_normal;
}

template<int NumSides, class Coords, class Random, class FloatType> AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::translate(const Coords& vector) {
  std::transform(&m_vertices[0], &m_vertices[NumSides], &m_vertices[0], std::bind2nd(std::plus<Coords>(), vector));
  recalGeomProperties();
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::mirror(const Coords& vector) {
  for (int i = 0; i < NumSides; i++) {
    m_vertices[i].SetCoordinates(m_vertices[i].X() * vector.X(), m_vertices[i].Y() * vector.Y(), m_vertices[i].Z() * vector.Z());
  }
  recalGeomProperties();
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::rotateX(FloatType angle) {
  std::transform(&m_vertices[0], &m_vertices[NumSides], &m_vertices[0], std::bind2nd(RotateY<Coords, FloatType>(), angle));
  recalGeomProperties();
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::rotateY(FloatType angle) {
  std::transform(&m_vertices[0], &m_vertices[NumSides], &m_vertices[0], std::bind2nd(RotateY<Coords, FloatType>(), angle));
  recalGeomProperties();
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::rotateZ(FloatType angle) {
  std::transform(&m_vertices[0], &m_vertices[NumSides], &m_vertices[0], std::bind2nd(RotateZ<Coords, FloatType>(), angle));
  recalGeomProperties();
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::operator<<(const Coords& vertex) {
  if (!isComplete()) m_vertices[m_numOfAllocVertices++] = vertex;
  if (isComplete()) computeProperties();
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::operator<<(const std::vector<Coords>& vertices) {
  for(typename std::vector<Coords>::const_iterator it = vertices.begin(); (!isComplete()) && (it != vertices.end()); ++it) {
    *this << *it;
  }
  return *this;
}


//!
//! A general polygon derived from abstract polygon class
//!
template<int NumSides> class Polygon3D : public AbstractPolygon<NumSides, ROOT::Math::XYZVector, TRandom> { // no checks are made on the convexity, but the algorithms in the class only work for convex polygons, so beware!
 public:
  typedef std::multiset<Polygon3D<3>, PolygonLess<Polygon3D<3> > > TriangleSet;

  Polygon3D() : AbstractPolygon<NumSides, ROOT::Math::XYZVector, TRandom>() {}
  Polygon3D(const ROOT::Math::XYZVector& vertex) : AbstractPolygon<NumSides, ROOT::Math::XYZVector, TRandom>(vertex) {}

  const TriangleSet& getTriangulation() const;

  //! Check if point p (defined by 3D vector) is inside a polygon -> use sum of triangle surfaces with one vertex being alwayas the point p.
  //! If sum equals to the polygon surface, point p is inside a polygon!
  bool  isPointInside(const ROOT::Math::XYZVector& point) const;

  //! Check if line at given direction (dir) and from given origin (orig) intersects the polygon -> return true/false
  //! Method to find the intersection with a plane:
  //! * A point p is inside a plane, when: p_vec.n_vec + d = 0 (1) (n_vec = normal & d = c_vec.n_vec, c_vec being e.g the central position of the polygon)
  //! * p_vec can be expressed as: p_vec = o_vec + dir_vec.k (2), where o_vec stands for origin, dir_vec for direction & k is the number we want to find out
  //! -> Hence k = o_vec + (c_vec.n_vec - o_vec.n_vec)/(dir_vec.n_vec) (3)
  bool  isLineIntersecting(const XYZVector& orig, const XYZVector& dir) const;

  //! Check if line at given direction (dir) and from given origin (orig) intersects the polygon -> return true/false + intersection point as 3D vector
  //! Method to find the intersection with a plane:
  //! * A point p is inside a plane, when: p_vec.n_vec + d = 0 (1) (n_vec = normal & d = c_vec.n_vec, c_vec being e.g the central position of the polygon)
  //! * p_vec can be expressed as: p_vec = o_vec + dir_vec.k (2), where o_vec stands for origin, dir_vec for direction & k is the number we want to find out
  //! -> Hence k = o_vec + (c_vec.n_vec - o_vec.n_vec)/(dir_vec.n_vec) (3)
  bool  isLineIntersecting(const XYZVector& orig, const XYZVector& dir, XYZVector& intersection) const;

  // Compute min/max Z/R positions (MinR is tricky, as it lies on the line between 2 vertices)
  double computeMinZ() const {return minget(this->begin(), this->end(), [](const XYZVector& v) { return v.Z();});};
  double computeMaxZ() const {return maxget(this->begin(), this->end(), [](const XYZVector& v) { return v.Z();});};
  double computeMinR() const;
  double computeMaxR() const {return maxget(this->begin(), this->end(), [](const XYZVector& v) { return v.Rho();});}

//  ROOT::Math::XYZVector generateRandomPoint(TRandom* dice) const; // Wrongly implemented - TODO: needs check

protected:

  void                   computeProperties();

  //! Compute
  std::vector<XYZVector> computeDistanceVectors() const;

  TriangleSet m_trianglesByArea; //!< All triangles
}; // Class


//!
//! A triangle: polygon with 3 vertices
//!
class Polygon3D<3>;
typedef Polygon3D<3> Triangle3D;

template<> class Polygon3D<3> : public AbstractPolygon<3, ROOT::Math::XYZVector, TRandom> {  // a triangle can be defined with more than 3 vertices. no error checking is made, simply the additional vertices besides the 3rd one are ignored

public:
  Polygon3D() : AbstractPolygon<3, ROOT::Math::XYZVector, TRandom>() {}
  Polygon3D(const ROOT::Math::XYZVector& vertex) : AbstractPolygon<3, ROOT::Math::XYZVector, TRandom>(vertex) {}
  virtual ~Polygon3D() {};

  bool isPointInside(const ROOT::Math::XYZVector& point) const;
//  ROOT::Math::XYZVector generateRandomPoint(TRandom* dice) const;

private:
  void computeProperties();

}; // Class

//
// Compute minR
//
template<int NumSides> double Polygon3D<NumSides>::computeMinR() const {

  std::vector<XYZVector> interVertexVectors;
  XYZVector v0 = this->getVertex(0);
  for (int i = 1; i<=NumSides; i++) {
     XYZVector v1 = this->getVertex(i % NumSides);
     v0.SetZ(0.0);
     v1.SetZ(0.0);
     XYZVector distVector;

     XYZVector dirv = (v0 - v1).Unit(); // direction vector of the line passing between v0 and v1
     double projv0 = v0.Dot(dirv);      // scalar projection of v0 on the direction vector
     double projv1 = v1.Dot(dirv);      // scalar projection of v1 on the direction vector
     if      (projv0 <= 0. && projv1 < 0.) distVector = v0;
     else if (projv0 > 0. && projv1 >= 0.) distVector = v1;
     else                                  distVector = v0 - (v0.Dot(dirv) * dirv);

     interVertexVectors.push_back(distVector);
     v0 = v1;
   }

  return minget(interVertexVectors.begin(), interVertexVectors.end(), [](const XYZVector& v) { return v.Rho();});
}


//
// Check if point p (defined by 3D vector) is inside a polygon -> use sum of triangle surfaces with one vertex being the point p. If sum equals to the polygon surface, point p is inside!
//
template<int NumSides> bool Polygon3D<NumSides>::isPointInside(const XYZVector& point) const {

  double sum = 0;
  const XYZVector *vbegin = &this->m_vertices[0], *vend = &this->m_vertices[NumSides];

  for (const XYZVector* it = vbegin; it != vend; ++it) {
    Triangle3D triangle;
    triangle << point << *it << (it+1 < vend ? *(it+1) : *vbegin);
    sum += triangle.getDoubleArea();
    if (sum - this->getDoubleArea() > 1e-4) return false;  // If sum area already bigger, quit -> the point can't be within the polygon
  }
  return fabs(this->getDoubleArea() - sum) < 1e-4;
}

//
// Check if line at given direction (dir) and from given origin (orig) intersects the polygon -> return true/false + intersection point as 3D vector
//
template<int NumSides> bool Polygon3D<NumSides>::isLineIntersecting(const XYZVector& orig, const XYZVector& dir, XYZVector& intersection) const {
  double normOrig      = this->getNormal().Dot(orig);
  double normDir       = this->getNormal().Dot(dir);
  double normCentrePos = this->getNormal().Dot(this->getCenter());

  // Check correct hemispere -> remember: shooting a half-line
  if (normDir < 1e-3) return false; // no fabs because if normDir < 0 we want to return false, as we're matching with the module in the opposite direction

  // Calculate intersection vector
  intersection = orig + (((normCentrePos - normOrig)/normDir) * dir); // TODO: slightly inefficient here: we do not need to go 3D and then 2D again...
  return isPointInside(intersection);
}

//
// Compute triangle basic properties
//
template<int NumSides> void Polygon3D<NumSides>::computeProperties() {
  this->m_area = 0;
  for (const XYZVector* it = &this->m_vertices[1]; it != &this->m_vertices[NumSides-1]; ++it) {
    Triangle3D triangle;
    triangle << this->m_vertices[0] << *it << *(it+1);
    double triangleArea = triangle.getDoubleArea();
    m_trianglesByArea.insert(triangle);
    this->m_area += triangleArea;
  }
}

// Wrongly implemented - TODO: needs check
//template<int NumSides> XYZVector Polygon3D<NumSides>::generateRandomPoint(TRandom* dice) const {
//  double binPicker = dice->Uniform(this->getDoubleArea());
//  double prevBin = 0;
//  TriangleSet::const_iterator it;
//  for (it = getTriangulation().begin(); it != getTriangulation().end(); it++) {
//    if (it->getDoubleArea() + prevBin >= binPicker) break;  // found triangle
//    prevBin += it->getDoubleArea();
//  }
//  return it->generateRandomPoint(dice);
//}

template<int NumSides> inline const std::multiset<Triangle3D, PolygonLess<Triangle3D> >& Polygon3D<NumSides>::getTriangulation() const {
  return m_trianglesByArea;
}

//
// Calculate basic properties of triangle: area etc.
//
inline void Triangle3D::computeProperties() {
   this->m_area = sqrt((this->m_vertices[1] - this->m_vertices[0]).Cross(this->m_vertices[2] - this->m_vertices[0]).Mag2());
}

//
// Check if point p (defined by 3D vector) is inside a polygon -> use sum of triangle surfaces with one vertex being the point p. If sum equals to the polygon surface, point p is inside!
//
inline bool Triangle3D::isPointInside(const XYZVector& point) const { // note: this is exactly the same as the non specialized case // inlined to avoid linker issues
  double sum = 0;
  const XYZVector *vbegin = &this->m_vertices[0], *vend = &this->m_vertices[3];
  for (const XYZVector* it = vbegin; it != vend; ++it) {
    Triangle3D triangle;
    triangle << point << *it << (it+1 < vend ? *(it+1) : *vbegin);
    sum += triangle.getDoubleArea();
    if (sum - this->getDoubleArea() > 1e-4) return false;  // If sum area already bigger, quit -> the point can't be within the polygon
  }
  return fabs(this->getDoubleArea() - sum) < 1e-4;
}

//
// Generate random point in a triangle - TODO: WRONG check
//
//inline XYZVector Triangle3D::generateRandomPoint(TRandom* dice) const {
//  double a = dice->Rndm(), b = dice->Rndm();
//  if (a + b > 1) { a = 1 - a; b = 1 - b; }
//  return m_vertices[0] + a*(m_vertices[1]-m_vertices[0]) + b*(m_vertices[2]-m_vertices[0]); // seriously C++??? you've been around for a while now, time to fix this horrid syntax??
//}

#endif // INCLUDE_POLYGON3D_H_




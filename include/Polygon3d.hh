#ifndef POLYGON3D_H
#define POLYGON3D_H

#include <vector>
#include <set>
#include <functional>
#include <algorithm>
#include <numeric>

#include <TRandom.h>
#include <Math/Vector3Dfwd.h>
#include <Math/GenVector/RotationZ.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/Transform3D.h>

using namespace ROOT::Math;

template<class Polygon>
struct PolygonLess {
  bool operator()(const Polygon& t1, const Polygon& t2) const {
    return t1.getDoubleArea() < t2.getDoubleArea();
  }
};

template<class Coords, class FloatType>
struct RotateX {};

template<class FloatType>
struct RotateX<XYZVector, FloatType> : public std::binary_function<XYZVector, FloatType, XYZVector> {
  XYZVector operator()(const XYZVector& vector, FloatType angle) const {
    return XYZVector(vector.X(),
                     vector.Y()*cos(angle) - vector.Z()*sin(angle),
                     vector.Y()*sin(angle) + vector.Z()*cos(angle));
  }
};

template<class Coords, class FloatType>
struct RotateY {};

template<class FloatType>
struct RotateY<XYZVector, FloatType> : public std::binary_function<XYZVector, FloatType, XYZVector> {
  XYZVector operator()(const XYZVector& vector, FloatType angle) const {
    return XYZVector(vector.Z()*sin(angle) + vector.X()*cos(angle),
                     vector.Y(),
                     vector.Z()*cos(angle) - vector.X()*sin(angle));
  }
};

template<class Coords, class FloatType>
struct RotateZ {};

template<class FloatType>
struct RotateZ<XYZVector, FloatType> : public std::binary_function<XYZVector, FloatType, XYZVector> {
  XYZVector operator()(const XYZVector& vector, FloatType angle) const {
    return XYZVector(vector.X()*cos(angle) - vector.Y()*sin(angle),
                     vector.X()*sin(angle) + vector.Y()*cos(angle),
                     vector.Z());
  }
};

template<class Coords> Coords crossProduct(const Coords& v1, const Coords& v2);

template<class Coords> Coords unitVector(const Coords& vector);

template<int NumSides, class Coords, class Random, class FloatType = double>
class AbstractPolygon {  // any number of dimensions, any number of sides, convex or concave. it only has 2 properties
protected:
  FloatType area_;
  Coords v_[NumSides]; // vertices of the polygon
  size_t allocated_;
  virtual void computeProperties() = 0;
  mutable Coords center_, normal_;
  mutable bool centerDirty_, normalDirty_;
  void setGeomDirty(bool dirty) { centerDirty_ = normalDirty_ = dirty; }
public:
  AbstractPolygon() : allocated_(0), centerDirty_(true), normalDirty_(true) {}
  AbstractPolygon(const Coords& vertex): allocated_(0), centerDirty_(true), normalDirty_(true) { *this << vertex; }
  virtual AbstractPolygon<NumSides, Coords, Random, FloatType>& operator<<(const Coords& vertex);
  virtual AbstractPolygon<NumSides, Coords, Random, FloatType>& operator<<(const std::vector<Coords>& vertices);
  virtual AbstractPolygon<NumSides, Coords, Random, FloatType>& operator()(const Coords& vertex) { *this << vertex; return *this; }
  const Coords* getVertices() const;
  const Coords* begin() const { return v_; }
  const Coords* end() const { return (v_+NumSides); }
  const Coords& getVertex(int index) const;
  int getNumSides() const { return NumSides; };
  bool isComplete() const;
  double getDoubleArea() const; // we save a division by 2 by returning the double area. for our purposes it's the same
  virtual bool isPointInside(const Coords& p) const = 0;
  virtual Coords generateRandomPoint(Random* die) const = 0;

  const Coords& getCenter() const;
  const Coords& getNormal() const; 
  AbstractPolygon<NumSides, Coords, Random, FloatType>& translate(const Coords& vector); 
  AbstractPolygon<NumSides, Coords, Random, FloatType>& mirror(const Coords& vector); 
  AbstractPolygon<NumSides, Coords, Random, FloatType>& rotateX(FloatType angle);
  AbstractPolygon<NumSides, Coords, Random, FloatType>& rotateY(FloatType angle);
  AbstractPolygon<NumSides, Coords, Random, FloatType>& rotateZ(FloatType angle);
};

template<int NumSides, class Coords, class Random, class FloatType> 
inline bool AbstractPolygon<NumSides, Coords, Random, FloatType>::isComplete() const { 
  return allocated_ == NumSides; 
}

template<int NumSides, class Coords, class Random, class FloatType> 
inline const Coords* AbstractPolygon<NumSides, Coords, Random, FloatType>::getVertices() const { 
  return v_; 
}

template<int NumSides, class Coords, class Random, class FloatType> 
inline const Coords& AbstractPolygon<NumSides, Coords, Random, FloatType>::getVertex(int index) const { 
  return v_[index]; 
}

template<int NumSides, class Coords, class Random, class FloatType> 
double AbstractPolygon<NumSides, Coords, Random, FloatType>::getDoubleArea() const {
  return area_;
}


template<int NumSides, class Coords, class Random, class FloatType> 
const Coords& AbstractPolygon<NumSides, Coords, Random, FloatType>::getCenter() const {
  if (centerDirty_) {
    center_ = std::accumulate(&v_[0], &v_[NumSides], XYZVector())/NumSides;
    centerDirty_ = false;
  }
  return center_;
}


template<int NumSides, class Coords, class Random, class FloatType> 
const Coords& AbstractPolygon<NumSides, Coords, Random, FloatType>::getNormal() const {
  if (normalDirty_) {
    normal_ = unitVector(crossProduct(v_[1]-v_[0], v_[2]-v_[0]));
    normalDirty_ = false;
  }
  return normal_;
}



template<int NumSides, class Coords, class Random, class FloatType>
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::translate(const Coords& vector) {
  std::transform(&v_[0], &v_[NumSides], &v_[0], std::bind2nd(std::plus<Coords>(), vector));
  setGeomDirty(true);
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType>
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::mirror(const Coords& vector) {
  for (int i = 0; i < NumSides; i++) {
    v_[i].SetCoordinates(v_[i].X() * vector.X(), v_[i].Y() * vector.Y(), v_[i].Z() * vector.Z());
  }
  setGeomDirty(true);
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType>
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::rotateX(FloatType angle) {
  std::transform(&v_[0], &v_[NumSides], &v_[0], std::bind2nd(RotateY<Coords, FloatType>(), angle));
  setGeomDirty(true);
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType>
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::rotateY(FloatType angle) {
  std::transform(&v_[0], &v_[NumSides], &v_[0], std::bind2nd(RotateY<Coords, FloatType>(), angle));
  setGeomDirty(true);
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType>
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::rotateZ(FloatType angle) {
  std::transform(&v_[0], &v_[NumSides], &v_[0], std::bind2nd(RotateZ<Coords, FloatType>(), angle));
  setGeomDirty(true);
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> 
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::operator<<(const Coords& vertex) {
  if (!isComplete()) v_[allocated_++] = vertex;
  if (isComplete()) computeProperties();
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> 
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::operator<<(const std::vector<Coords>& vertices) {
  for(typename std::vector<Coords>::const_iterator it = vertices.begin(); (!isComplete()) && (it != vertices.end()); ++it) {
    *this << *it;
  }
  return *this;
}



template<int NumSides>
class Polygon3d : public AbstractPolygon<NumSides, ROOT::Math::XYZVector, TRandom> { // no checks are made on the convexity, but the algorithms in the class only work for convex polygons, so beware!
public:
  typedef std::multiset<Polygon3d<3>, PolygonLess<Polygon3d<3> > > TriangleSet;
protected:
  TriangleSet trianglesByArea_;
  void computeProperties();
public:
  Polygon3d() : AbstractPolygon<NumSides, ROOT::Math::XYZVector, TRandom>() {}
  Polygon3d(const ROOT::Math::XYZVector& vertex) : AbstractPolygon<NumSides, ROOT::Math::XYZVector, TRandom>(vertex) {}
  const TriangleSet& getTriangulation() const;
  bool isPointInside(const ROOT::Math::XYZVector& p) const;
  bool isLineIntersecting(const XYZVector& orig, const XYZVector& dir) const;
  bool isLineIntersecting(const XYZVector& orig, const XYZVector& dir, XYZVector& intersection) const;
  ROOT::Math::XYZVector generateRandomPoint(TRandom* die) const;
};

template<>
class Polygon3d<3> : public AbstractPolygon<3, ROOT::Math::XYZVector, TRandom> {  // a triangle can be defined with more than 3 vertices. no error checking is made, simply the additional vertices besides the 3rd one are ignored
  void computeProperties();
public:
  Polygon3d() : AbstractPolygon<3, ROOT::Math::XYZVector, TRandom>() {}
  Polygon3d(const ROOT::Math::XYZVector& vertex) : AbstractPolygon<3, ROOT::Math::XYZVector, TRandom>(vertex) {}
  bool isPointInside(const ROOT::Math::XYZVector& p) const;
  ROOT::Math::XYZVector generateRandomPoint(TRandom* die) const;
};

typedef Polygon3d<3> Triangle3d;



template<int NumSides>
bool Polygon3d<NumSides>::isPointInside(const XYZVector& p) const {
  double sum = 0;
  const XYZVector *vbegin = &this->v_[0], *vend = &this->v_[NumSides];
  for (const XYZVector* it = vbegin; it != vend; ++it) {
    Triangle3d t;
    t << p << *it << (it+1 < vend ? *(it+1) : *vbegin);
    sum += t.getDoubleArea();
    if (sum - this->getDoubleArea() > 1e-4) return false;  // early quit if sum area is already bigger
  }
  return fabs(this->getDoubleArea() - sum) < 1e-4;
}

template<int NumSides>
bool Polygon3d<NumSides>::isLineIntersecting(const XYZVector& orig, const XYZVector& dir, XYZVector& intersection) const {
  double normOrig = this->getNormal().Dot(orig);
  double normDir = this->getNormal().Dot(dir);
  double d = this->getCenter().Dot(this->getNormal());
  if (normDir < 1e-3) return false; // no fabs because if normDir < 0 we want to return false, as we're matching with the module in the opposite direction
  intersection = orig + (((d - normOrig)/normDir) * dir); // TODO: slightly inefficient here: we do not need to go 3D and then 2D again...
  return isPointInside(intersection);
}

template<int NumSides>
void Polygon3d<NumSides>::computeProperties() {
  this->area_ = 0;
  for (const XYZVector* it = &this->v_[1]; it != &this->v_[NumSides-1]; ++it) {
    Triangle3d t;
    t << this->v_[0] << *it << *(it+1);
    double triangleArea = t.getDoubleArea();
    trianglesByArea_.insert(t);
    this->area_ += triangleArea;
  }
}

template<int NumSides> 
XYZVector Polygon3d<NumSides>::generateRandomPoint(TRandom* die) const {
  double binPicker = die->Uniform(this->getDoubleArea());
  double prevBin = 0;
  TriangleSet::const_iterator it;
  for (it = getTriangulation().begin(); it != getTriangulation().end(); it++) {
    if (it->getDoubleArea() + prevBin >= binPicker) break;  // found triangle
    prevBin += it->getDoubleArea();
  }
  return it->generateRandomPoint(die);
}


template<int NumSides>
inline const std::multiset<Triangle3d, PolygonLess<Triangle3d> >& Polygon3d<NumSides>::getTriangulation() const {
  return trianglesByArea_;
}


inline void Triangle3d::computeProperties() {
   this->area_ = sqrt((this->v_[1] - this->v_[0]).Cross(this->v_[2] - this->v_[0]).Mag2());
}

inline bool Triangle3d::isPointInside(const XYZVector& p) const { // note: this is exactly the same as the non specialized case // inlined to avoid linker issues
  double sum = 0;
  const XYZVector *vbegin = &this->v_[0], *vend = &this->v_[3];
  for (const XYZVector* it = vbegin; it != vend; ++it) {
    Triangle3d t;
    t << p << *it << (it+1 < vend ? *(it+1) : *vbegin);
    sum += t.getDoubleArea();
    if (sum - this->getDoubleArea() > 1e-4) return false;  // early quit if sum area is already bigger
  }
  return fabs(this->getDoubleArea() - sum) < 1e-4;
}

inline XYZVector Triangle3d::generateRandomPoint(TRandom* die) const {
  double a = die->Rndm(), b = die->Rndm();
  if (a + b > 1) { a = 1 - a; b = 1 - b; }
  return this->v_[0] + a*(this->v_[1]-this->v_[0]) + b*(this->v_[2]-this->v_[0]); // seriously C++??? you've been around for a while now, time to fix this horrid syntax??
}



#endif




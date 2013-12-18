#ifndef GEOMETRIC_MODULE_H
#define GEOMETRIC_MODULE_H

/// ===================================================== GEOMETRIC MODULES =====================================================
//

#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <array>

#include <boost/property_tree/info_parser.hpp>

#include <TVector.h>
#include <TMatrixD.h>

#include "global_funcs.h"
#include "Polygon3d.h"
#include "Property.h"
#include "ModuleBase.h"

using std::vector;
using std::string;
using std::array;

//
namespace ModuleHelpers {
  double polygonAperture(const Polygon3d<4>& poly); 
};


enum ModuleShape { RECTANGULAR, WEDGE };

class GeometricModule : public ModuleDecorable {
protected:
  int numHits_ = 0;
  int tiltAngle_ = 0., skewAngle_ = 0.;
  Polygon3d<4> basePoly_;
  double triangleCross(const XYZVector& P1, const XYZVector& P2, const XYZVector& P3, const XYZVector& PL, const XYZVector& PU);
public:
  Property<double, Default> dsDistance; // a GeometricModule is a purely 2d geometric object represented in 3d space with just a polygon and an additional for thickness value for tracker geometry construction
  GeometricModule() :
      dsDistance("dsDistance", parsedAndChecked(), 0.)
  {}
  void setup() {}

  const Polygon3d<4>& basePoly() const { return basePoly_; }
  const XYZVector& center() const { return basePoly_.getCenter(); }
  const XYZVector& normal() const { return basePoly_.getNormal(); }
  virtual double area() const = 0;
  virtual double length() const = 0;
  virtual double maxWidth() const = 0;
  virtual double minWidth() const = 0;
  virtual double meanWidth() const = 0;
  double thickness() const { return dsDistance() + 0.1; } // for Geometric modules it is assumed they have a 0.1 mm thick generic sensor
 
  double tiltAngle() const { return tiltAngle_; }
  double skewAngle() const { return skewAngle_; }

  void translate(const XYZVector& vector) { basePoly_.translate(vector); }
  void mirror(const XYZVector& vector) { basePoly_.mirror(vector); }
  void rotateX(double angle) { basePoly_.rotateX(angle); tiltAngle_ += angle; }
  void rotateY(double angle) { basePoly_.rotateY(angle); skewAngle_ += angle; }
  void rotateZ(double angle) { basePoly_.rotateZ(angle); }

  virtual GeometricModule* clone() = 0;

  virtual void accept(GeometryVisitor& v) = 0;
  virtual void accept(ConstGeometryVisitor& v) const = 0;
  virtual void build() = 0;

  double trackCross(const XYZVector& PL, const XYZVector& PU);
  int numHits() const { return numHits_; }
  void resetHits() { numHits_ = 0; }

  virtual ModuleShape shape() const = 0;
};




class RectangularModule : public GeometricModule {
  Property<double, NoDefault> length_;
public:
  Property<double, NoDefault> width;
  Property<double, NoDefault> aspectRatio;
  Property<double, Default> waferDiameter;

  RectangularModule() :
      length_("length", parsedOnly()), // not checked because a custom checker is defined for RectangularModules
      width("width", parsedOnly()),  // same here
      aspectRatio("aspectRatio", parsedOnly()), // same here
      waferDiameter("waferDiameter", parsedAndChecked(), 131.)
  {}

  virtual double area() const override { return length()*width(); }
  virtual double length() const override { return length_(); }
  virtual double maxWidth() const override { return width(); }
  virtual double minWidth() const override { return width(); }
  virtual double meanWidth() const override { return width(); }

  void length(double l) { length_(l); }

  virtual void accept(GeometryVisitor& v) { v.visit(*this); v.visit(*(GeometricModule*)this); }
  virtual void accept(ConstGeometryVisitor& v) const { v.visit(*this); v.visit(*(const GeometricModule*)this); }
  virtual RectangularModule* clone() override { return new RectangularModule(*this); }
  virtual void check() override;
  virtual void build() override;

  virtual ModuleShape shape() const { return RECTANGULAR; }
};




class WedgeModule : public GeometricModule {
  double length_, minWidth_, maxWidth_;
  double area_/*dist_*/;
  bool cropped_;
  double amountCropped_;
public:
  Property<float, NoDefault> buildAperture, buildDistance, buildCropDistance;
  Property<float, Default> waferDiameter;

  
  WedgeModule() :
      buildAperture("buildAperture", parsedAndChecked()),
      buildDistance("buildDistance", parsedAndChecked()),
      buildCropDistance("buildCropDistance", parsedAndChecked()),
      waferDiameter("waferDiameter", parsedAndChecked(), 131.)
  {}

  virtual double area() const override { return area_; }
  virtual double length() const override { return length_; }
  virtual double minWidth() const override { return minWidth_; }
  virtual double maxWidth() const override { return maxWidth_; }
  virtual double meanWidth() const override { return (maxWidth_ + minWidth_)/2; }

  virtual void accept(GeometryVisitor& v) { v.visit(*this); v.visit(*(GeometricModule*)this); }
  virtual void accept(ConstGeometryVisitor& v) const { v.visit(*this); v.visit(*(const GeometricModule*)this); }
  virtual WedgeModule* clone() override { return new WedgeModule(*this); }

  virtual void build() override;

  virtual ModuleShape shape() const { return WEDGE; }
};



// ========================================================================================================================================
//
#endif

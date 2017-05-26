#ifndef INCLUDE_GEOMETRIC_MODULE_H_
#define INCLUDE_GEOMETRIC_MODULE_H_

#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <array>

#include <boost/property_tree/info_parser.hpp>

#include <TVector.h>
#include <TMatrixD.h>

#include "Polygon3D.h"
#include "Property.h"
#include "Visitor.h"
#include "Visitable.h"

// Using namespace
using std::vector;
using std::string;
using std::array;

// Definition of global function
namespace ModuleHelpers {
  double polygonAperture(const Polygon3D<4>& poly);
};

// Enum definition related to geometrical shape of modules
enum ModuleShape { RECTANGULAR, WEDGE };

/*
 * @class GeometricModule
 * @brief Base class defining general geometric properties of detector module - shapes
 */
class GeometricModule : public PropertyObject, public Buildable, public Placeable, public Identifiable<int>, public Visitable {

 public:

  //! Constructor
  GeometricModule();

  //! Destructor
  virtual ~GeometricModule() {}

  //! Purely virtual clone function -> needs to be implemented for all derived classes to make cloning of DetectorModule class properly working
  virtual GeometricModule* clone() const = 0;

  //! Translate module by given vector
  void translate(const XYZVector& vector) { m_basePoly.translate(vector); }

  // Rotate module by given angle
  void mirror(const XYZVector& vector) { m_basePoly.mirror(vector); }
  void rotateX(double angle)           { m_basePoly.rotateX(angle); }
  void rotateY(double angle)           { m_basePoly.rotateY(angle); }
  void rotateZ(double angle)           { m_basePoly.rotateZ(angle); }
  void tiltAngle(double angle)         { m_tiltAngle += angle; }
  void skewAngle(double angle)         { m_skewAngle += angle; }
  bool flipped() const                 { return m_flipped; }
  bool flipped(bool newFlip)           { m_flipped=newFlip; return m_flipped; }

  //! GeometryVisitor pattern (purely virtual method) -> module needs to be visitable
  virtual void accept(GeometryVisitor& v) = 0;

  //! GeometryVisitor pattern (purely virtual method - const. option) -> module needs to be visitable
  virtual void accept(ConstGeometryVisitor& v) const = 0;

  //! Build geometrical module
  virtual void build() = 0;

  // Getter methods
  const Polygon3D<4>& basePoly() const { return m_basePoly; }
  const XYZVector& center()      const { return m_basePoly.getCenter(); }
  const XYZVector& normal()      const { return m_basePoly.getNormal(); }
  virtual double area()          const = 0;
  virtual double length()        const = 0;
  virtual double maxWidth()      const = 0;
  virtual double minWidth()      const = 0;
  virtual double meanWidth()     const = 0;
  //double thickness()             const { return dsDistance() + 0.1; } // Defined at the level of detector module (for Geometric modules it is assumed they have a 0.1 mm thick generic sensor)

  double tiltAngle()             const { return m_tiltAngle; }
  double skewAngle()             const { return m_skewAngle; }
  virtual ModuleShape shape()    const = 0;

  // ??? TODO: Still needed?
  //double trackCross(const XYZVector& PL, const XYZVector& PU);
  //void   resetHits() { m_numHits = 0; }
  //int    numHits()               const { return m_numHits; }

  Property<double, Default> dsDistance; // a GeometricModule is a purely 2d geometric object represented in 3d space with just a polygon and an additional for thickness value for tracker geometry construction
  Property<double, Default> physicalLength;

 protected:

  // ??? TODO: Still needed?
  //double triangleCross(const XYZVector& P1, const XYZVector& P2, const XYZVector& P3, const XYZVector& PL, const XYZVector& PU);
  //int          m_numHits   = 0;

  bool         m_flipped   = false;
  double       m_tiltAngle = 0.; //!< Module tilt, i.e. rotation in RZ plane
  double       m_skewAngle = 0.; //!< Module skew, i.e. rotation in XY plane
  Polygon3D<4> m_basePoly;

}; // Class

/*
 * @class RectangularModule
 * @details Class derived from GeometricModule class implementing rectangular shape
 */
class RectangularModule : public GeometricModule, public Clonable<RectangularModule> {

 public:

  //! Constructor
  RectangularModule();

  //! Clone method
  virtual RectangularModule* clone() const { return Clonable<RectangularModule>::clone(); }

  //! GeometryVisitor pattern (purely virtual method) -> module needs to be visitable
  virtual void accept(GeometryVisitor& v)            { v.visit(*this); v.visit(*(GeometricModule*)this); }
  virtual void accept(ConstGeometryVisitor& v) const { v.visit(*this); v.visit(*(const GeometricModule*)this); }

  //! Build rectangular shape
  virtual void build() override;

  //! Cross-check that parameters set correctly
  virtual void check() override;

  // Getter methods
  virtual double area()       const override { return length()*width(); }
  virtual double length()     const override { return m_length(); }
  virtual double maxWidth()   const override { return width(); }
  virtual double minWidth()   const override { return width(); }
  virtual double meanWidth()  const override { return width(); }
  virtual ModuleShape shape() const override { return RECTANGULAR; }

  Property<double, NoDefault> width;         //!< Module width
  Property<double, NoDefault> aspectRatio;   //!< Length to width aspect ratio
  Property<double, Default>   waferDiameter; //!< Build module from wafer of given diameter

private:

  Property<double, NoDefault> m_length;      //!< Module length

}; // Class

/*
 * @class WedgeModule
 * @details Class derived from GeometricModule class implementing wedge shape
 */
class WedgeModule : public GeometricModule, public Clonable<WedgeModule> {

public:


  //! Constructor
  WedgeModule();

  //! Clone method
  virtual WedgeModule* clone() const { return Clonable<WedgeModule>::clone(); }

  //! GeometryVisitor pattern (purely virtual method) -> module needs to be visitable
  virtual void accept(GeometryVisitor& v)            { v.visit(*this); v.visit(*(GeometricModule*)this); }
  virtual void accept(ConstGeometryVisitor& v) const { v.visit(*this); v.visit(*(const GeometricModule*)this); }

  //! Build wedge shape
  virtual void build() override;

  // Getter methods
  virtual double area()       const override { return m_area; }
  virtual double length()     const override { return m_length; }
  virtual double minWidth()   const override { return m_minWidth; }
  virtual double maxWidth()   const override { return m_maxWidth; }
  virtual double meanWidth()  const override { return (m_maxWidth + m_minWidth)/2; }
  virtual ModuleShape shape() const override { return WEDGE; }

  Property<float, NoDefault> buildAperture;     //!< dPhi angle which the wedge-shaped detector has to cover
  Property<float, NoDefault> buildDistance;     //!< Start building the wedge-shaped detector at this radius
  Property<float, NoDefault> buildCropDistance; //!< Stop building the wedge-shaped detector at the given crop distance
  Property<float, Default>   waferDiameter;     //!< Use the wafer of given diameter to build wedge-shaped module

private:

  double m_length;        //!< Wafer length (calculated
  double m_minWidth;      //!< The smaller width of the isosceles trapezoid
  double m_maxWidth;      //!< The bigger width of the isosceles trapezoid
  double m_area/*dist_*/; //!< Trapezoid area
  bool   m_cropped;
  double m_amountCropped;

  const double c_defaultWaferDiam = 131;


}; // Class

#endif /* INCLUDE_GEOMETRIC_MODULE_H_ */

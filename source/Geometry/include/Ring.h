#ifndef INCLUDE_RING_H
#define INCLUDE_RING_H

#include <memory>
#include <vector>
#include <string>
#include <limits.h>

using std::vector;
using std::string;

#include "math_functions.h"
#include "Property.h"
#include "DetectorModule.h"
#include "Visitable.h"

// Typedefs
typedef PtrVector<Ring>            Rings;
typedef PtrVector<TiltedRing>      TiltedRings;
typedef std::map<int, const Ring*> RingIndexMap;

// Enums
enum BuildDirection { TOPDOWN, BOTTOMUP };

/*
 * @class TiltedRing
 * @details Tilted ring keeps information about individual tilted modules positioned in a barrel layer, set in so-called rings.
 * Two rods of modules are set, inner & outer (similar concept to RodPair, with modules positioned at +-smallDelta). Modules are
 * currently built using Tilted Rod, so not here, within the Tilted ring concept. TODO!
 *
 */
class TiltedRing : public PropertyObject, public Buildable, public Identifiable<int>, public Visitable {

 public:

  //! Constructor
  TiltedRing(int id, int numRods, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty);

  //! Build the ring
  void build(double lastThetaEnd);

  //! Cross-check parameters provided from geometry configuration file
  void check() override;

  //! Calculate geometry properties with respect to previous ring or flat part defined by several parameters
  void calculateGeometryProperties(const TiltedRing* prevRing, double flatPartREndInner, double flatPartREndOuter, double flatPartEndModCenterZ, double flatPartZEnd);

  //! GeometryVisitor pattern -> ring visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> ring visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  //! TODO: Return ring modules -> modules defined in tilted rods
  // const BarrelModules& modules() const { return m_modules; }

  // Ring defined by following parameters
  Property<double, NoDefault> innerRadius;          //!< Radius of the inner rod/module (must be specified)
  Property<double, NoDefault> outerRadius;          //!< Radius of the outer rod/module (must be specified
  Property<double, NoDefault> zInner;               //!< Z position of the inner rod/module (define either zInner, zOuter or zOverlap)
  Property<double, NoDefault> zOuter;               //!< Z position of the outer rod/module (define either zInner, zOuter or zOverlap)
  Property<double, NoDefault> ringZOverlap;         //!< Required overlap in Z with respect to the previous ring -> as projected on the current tilted ring module (so effective measure of overlap in Z on the current module)
  Property<double, NoDefault> tiltAngle;            //!< Module tilt angle expected in config file in degrees (0 for BRL, 90deg for ENDCAP) -
  Property<double, NoDefault> theta_g;              //!< Angle between pos. Z axis and line connecting centres of outer & inner modules in degrees (expected in degerees in config file). Ideally 90deg - tiltAngle

  // Calculated quantities
  Property<double, Computable> zErrorInner;         //!< Beam spot coverage in Z by inner rod/module -> only active module parameters taken into account
  Property<double, Computable> zErrorOuter;         //!< Beam spot coverage in Z by outer rod/module -> only active module parameters taken into account
  Property<double, Computable> deltaZInner;         //!< Difference between given ring module center & previous ring module center
  Property<double, Computable> deltaZOuter;         //!< Difference between given ring module center & previous ring module center

  Property<double, Computable> thetaStart;                                            //!< Ring module rods start at given theta, namely outer rod modules starts here
  double thetaEnd() const { return MAX(thetaEndInner_REAL(), thetaEndOuter_REAL()); } //!< To get hermetic coverage, get the worse scenario of where the Ring module rods end-up, i.e. outer rod modules

  Property<double, Computable> tiltAngleIdealInner; //!< Ideal theta angle of the module positioned in inner rod
  Property<double, Computable> deltaTiltIdealInner; //!< Difference between the real angle and ideal angle for inner module;
  Property<double, Computable> tiltAngleIdealOuter; //!< Ideal theta angle of the module positioned in outer rod
  Property<double, Computable> deltaTiltIdealOuter; //!< Difference between the real angle and ideal angle for outer module;
  Property<double, Computable> phiOverlap;          //!< Minimum modules overlap in R-Phi in length units

  Property<double, Computable> rStartOuter_REAL;    //!< Ring modules in the outer rod start at given radius
  Property<double, Computable> rEndOuter_REAL;      //!< Ring modules in the outer rod end-up at given radius
  Property<double, Computable> zStartOuter_REAL;    //!< Ring modules in the outer rod start at given z-pos
  Property<double, Computable> zEndOuter_REAL;      //!< Ring modules in the outer rod end-up at given z-pos
  Property<double, Computable> thetaEndOuter_REAL;  //!< Theta angle for outer rod defined as: tan(radiusEnd/ZEnd)

  Property<double, Computable> rStartInner_REAL;    //!< Ring modules in the inner rod start at given radius
  Property<double, Computable> rEndInner_REAL;      //!< Ring modules in the inner rod end-up at given radius
  Property<double, Computable> zStartInner_REAL;    //!< Ring modules in the inner rod start at given z-pos
  Property<double, Computable> zEndInner_REAL;      //!< Ring modules in the inner rod end-up at given z-pos
  Property<double, Computable> thetaEndInner_REAL;  //!< Theta angle for inner rod defined as: tan(radiusEnd/ZEnd)

  double averageR() const { return (innerRadius() + outerRadius()) / 2.; }
  double averageZ() const { return (zInner() + zOuter()) / 2.; }

  double deltaR()   const { return outerRadius() - innerRadius(); }
  double gapR()     const { return (outerRadius() - innerRadius()) / sin(theta_g() * M_PI / 180.); }
  double deltaZ()   const { return zOuter() - zInner(); }

 private :

  // Helper function building the modules in a ring
  void buildLeftRight(double lastThetaEnd);

  // TODO: Modules defined in tilted rods
  //BarrelModules m_modules; //!< Modules built within a ring
  //MaterialObject materialObject_;

  int    m_numRods;                  //!< Number of rods (modules) in a ring along phi

}; // Class

/*
 * @class Ring
 * @details Ring keeps information about disk modules grouped together at given radius to a ring. It's building procedure
 * is executed automatically via Disk class. Similarly, all its modules are built through build() method. Depending on
 * the configuration, building algorithm uses either the bottom-up or up-down approach. i.e. from the closest ring to the
 * beam-pipe to the furthest one or vice-versa. Bottom-up approach can build ring using wedge-shaped or rectangular modules,
 * the topdown only rectangular modules.
 */
class Ring : public PropertyObject, public Buildable, public Identifiable<int>, public Visitable {

public:

  //! Constructor - specify unique id, build direction & parse geometry config file using boost property tree & read-in module parameters
  Ring(int id, BuildDirection direction, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty);

  //! Build all modules -> use either bottom-up or up-down approach
  void build(double radius, double zOffset);

  //! Setup: link lambda functions to various layer related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Cross-check parameters provided from geometry configuration file
  void check() override;

  //! Limit ring geometry by eta cut
  void cutAtEta(double eta);

  //! GeometryVisitor pattern -> ring visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> ring visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  //! Helper method translating Ring z position by given offset
  void translateZ(double zOffset);

  //! Helper method duplicating the whole Ring from zPos to -zPos or vice versa
  void rotateToNegativeZSide();

  //! Return average disc Z position
  double averageZ()  const { return m_averageZ; }

  //! Return module thickness
  double thickness() const { return smallDelta()*2 + maxModuleThickness(); }

  //! Return ring modules
  const EndcapModules& modules() const { return m_modules; }

  const MaterialObject& materialObject() const {return m_materialObject;}

  ReadonlyProperty<double, Computable>  maxModuleThickness;
  Property<        double, Computable>  minZ;            //!< Minimum rod Z pos
  Property<        double, Computable>  maxZ;            //!< Maximum rod Z pos
  Property<        double, Computable>  minZAllMat;      //!< Minimum rod Z position taking into account all material structures
  Property<        double, Computable>  maxZAllMat;      //!< Maximum rod Z position taking into account all material structures
  Property<        double, Computable>  minR;            //!< Minimum rod radius
  Property<        double, Computable>  maxR;            //!< Maximum rod radius
  Property<        int   , AutoDefault> disk;
  Property<        double, NoDefault>   buildCropRadius; //!< Crop ring at given radius
  Property<        int   , NoDefault>   numModules;      //!< Number of modules -> if set forces the number of modules (in phi) to be exactly numModules

  ReadonlyProperty<double, NoDefault>   smallDelta;      //!< In each rod the modules are positioned in R-Phi zig-zag at average zPos +- smallDelta
  Property<        double, Default>     zRotation;       //!< Rotate ring around Z axis by given angle
  Property<        double, Default>     ringOuterRadius; //!< Define outer ring radius -> if not defined set optimally from disk level
  Property<        double, Default>     ringInnerRadius; //!< Define inner ring radius -> if not defined set optimally from disk level

  Property<        double, Default>     phiOverlap;      //!< Required modules overlap in R-Phi direction (overlap in mm)
  Property<        int   , Default>     phiSegments;     //!< Ring symmetry - number of symmetric segments, e.g 4 -> 90deg symmetry

private:

   EndcapModules  m_modules;        //!< Modules built within a ring
   MaterialObject m_materialObject;

   //! If modules built within a ring using approach bottom-up
   void buildBottomUp(double radius);

   //! If modules built within a ring using approach up-down
   void buildTopDown(double radius);

   //! Helper build method -> build modules based on template module
   void buildModules(EndcapModule* templ, int numMods, double smallDelta);

   //! Helper method calculating for given wedge-shaped modules at given radius their shape (covering phi aperture)
   double computeTentativePhiAperture(double moduleWaferDiameter, double minRadius);

   //! Helper method computing ring optimal parameters when wedge-shaped modules used
   std::pair<double, int> computeOptimalRingParametersWedge(double moduleWaferDiameter, double minRadius);

   //! Helper method computing ring optimal parameters when rectangular-shaped modules used
   std::pair<double, int> computeOptimalRingParametersRectangle(double moduleWidth, double maxRadius);

   template<class T> int roundToOdd(T x) { return round((x-1)/2)*2+1; }
   double solvex(double y);
   double compute_l(double x, double y, double d);
   double compute_d(double x, double y, double l);

   Property<   ModuleShape, NoDefault> m_moduleShape;            //!< Use wedge or rectangular-shaped modules
   Property<BuildDirection, Default>   m_buildDirection;         //!< Use bottom-up or top-down approach -> initialized from disk level if not set by user
   Property<        bool  , Default>   m_requireOddModsPerSlice; //!< Require odd number of modules per segment?
   Property<        int   , Default>   m_additionalModules;      //!< Automatically calculates number of modules and adds this number to it
   Property<        bool  , Default>   m_alignEdges;             //!< Start building modules in the ring with its edge or centre
   Property<        double, Default>   m_ringGap;                //!< Required gap in-between rings
   Property<        int   , Default>   m_smallParity;            //!< In each rod the modules are positioned in R-Phi zig-zag at average zPos +- smallDelta -> use +- or -+smallDelta

   double m_averageZ = 0;      //!< Average Z position of the given ring

   const int c_maxWedgeCalcLoops = 100; //!< Maximum number of loops to be iterated when finding optimal positions for wege-shaped modules

}; // Class

#endif /* INCLUDE_RING_H */

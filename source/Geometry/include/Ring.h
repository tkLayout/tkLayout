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
typedef std::map<int, const Ring*> RingIndexMap;

// Enums
enum BuildDirection { TOPDOWN, BOTTOMUP };





class TiltedRing : public PropertyObject, public Buildable, public Identifiable<int>, public Visitable {

 public:
  typedef PtrVector<BarrelModule> Container;
 private :
  Container modules_;
  //MaterialObject materialObject_;
  
  double thetaOuterUP_, thetaOuterDOWN_, thetaOuter_, tiltAngleIdealOuter_, deltaTiltIdealOuter_;
  double thetaInner_, tiltAngleIdealInner_, deltaTiltIdealInner_;

  double thetaStart_, thetaEnd_;

  double thetaStartInner_, thetaEndInner_;

  double numPhi_, phiOverlap_;

  double rStartOuter_REAL_, zStartOuter_REAL_, rEndOuter_REAL_, zEndOuter_REAL_;
  double rStartInner_REAL_, zStartInner_REAL_, rEndInner_REAL_, zEndInner_REAL_;


 public:
  Property<double, NoDefault> innerRadius;
  Property<double, NoDefault> outerRadius;
  Property<double, NoDefault> zInner;
  Property<double, NoDefault> zOuter;
  Property<double, NoDefault> tiltAngle;
  Property<double, NoDefault> theta_g;
  Property<double, NoDefault> ringZOverlap;

  const Container& modules() const { return modules_; }

 TiltedRing() :
    innerRadius           ("ringInnerRadius"       , parsedAndChecked()),
    outerRadius           ("ringOuterRadius"       , parsedAndChecked()),
    zInner                ("ringInnerZ"            , parsedOnly()),
    zOuter                ("ringOuterZ"            , parsedOnly()),
    tiltAngle             ("tiltAngle"             , parsedAndChecked()),
    theta_g               ("theta_g"               , parsedAndChecked()),
    ringZOverlap          ("ringZOverlap"          , parsedOnly())
      {}

  void build(double lastThetaEnd);
  void buildLeftRight(double lastThetaEnd);
  void check() override;

  void accept(GeometryVisitor& v) {
    v.visit(*this); 
    for (auto& m : modules_) { m.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const {
    v.visit(*this); 
    for (const auto& m : modules_) { m.accept(v); }
    }

  double thetaOuter() const { return thetaOuter_; }
  double thetaInner() const { return thetaInner_; }
  double thetaEnd() const { return thetaEnd_; }


  double tiltAngleIdealOuter() const { return tiltAngleIdealOuter_; }
  double deltaTiltIdealOuter() const { return deltaTiltIdealOuter_; }

  double tiltAngleIdealInner() const { return tiltAngleIdealInner_; }
  double deltaTiltIdealInner() const { return deltaTiltIdealInner_; }

  double averageR() const { return (innerRadius() + outerRadius()) / 2.; }
  double averageZ() const { return (zInner() + zOuter()) / 2.; }

  double deltaR() const { return outerRadius() - innerRadius(); }
  double gapR() const { return (outerRadius() - innerRadius()) / sin(theta_g() * M_PI / 180.); }
  double deltaZ() const { return zOuter() - zInner(); }

  void numPhi(double numPhi) { numPhi_ = numPhi; }
  double numPhi() const { return numPhi_; }
  double phiOverlap() const { return  phiOverlap_; }

  double rStartOuter_REAL() const { return rStartOuter_REAL_; }
  double zStartOuter_REAL() const { return zStartOuter_REAL_; } 
  double rEndOuter_REAL() const { return rEndOuter_REAL_; }
  double zEndOuter_REAL() const { return zEndOuter_REAL_; }
  double thetaEndOuter_REAL () const { return atan(rEndOuter_REAL_ / zEndOuter_REAL_); }

  double rStartInner_REAL() const { return rStartInner_REAL_; }
  double zStartInner_REAL() const { return zStartInner_REAL_; } 
  double rEndInner_REAL() const { return rEndInner_REAL_; }
  double zEndInner_REAL() const { return zEndInner_REAL_; }
  double thetaEndInner_REAL () const { return atan(rEndInner_REAL_ / zEndInner_REAL_); }

};





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

  //! Helper method mirroring the whole Ringc from zPos to -zPos or vice versa
  void mirrorZ();

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

};

#endif /* INCLUDE_RING_H */

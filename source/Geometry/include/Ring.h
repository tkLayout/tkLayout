#ifndef INCLUDE_RING_H
#define INCLUDE_RING_H

#include <vector>
#include <string>
#include <memory>

#include <boost/ptr_container/ptr_vector.hpp>

using std::vector;
using std::string;

#include "global_funcs.h"
#include "Property.h"
#include "DetectorModule.h"
#include "Visitable.h"

// Typedefs
typedef PtrVector<EndcapModule> EndcapModules;

// Enums
enum BuildDirection { TOPDOWN, BOTTOMUP };

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
  Property<        double, Computable>  minR;            //!< Minimum rod radius
  Property<        double, Computable>  maxR;            //!< Maximum rod radius
  Property<        int   , AutoDefault> disk;
  Property<        double, NoDefault>   buildCropRadius; //!< Crop ring at given radius
  Property<        int   , NoDefault>   numModules;      //!< Number of modules -> if set forces the number of modules (in phi) to be exactly numModules

  ReadonlyProperty<double, NoDefault>   smallDelta;      //!< In each rod the modules are positioned in R-Phi zig-zag at average zPos +- smallDelta
  Property<        double, Default>     zRotation;       //!< Rotate ring around Z axis by given angle
  Property<        double, Default>     ringOuterRadius; //!< Define outer ring radius -> if not defined set optimally from disk level
  Property<        double, Default>     ringInnerRadius; //!< Define inner ring radius -> if not defined set optimally from disk level

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
   Property<        double, Default>   m_phiOverlap;             //!< Required modules overlap in R-Phi direction (overlap in mm)
   Property<        bool  , Default>   m_requireOddModsPerSlice; //!< Require odd number of modules per segment?
   Property<        int   , Default>   m_phiSegments;            //!< Ring symmetry - number of symmetric segments, e.g 4 -> 90deg symmetry
   Property<        int   , Default>   m_additionalModules;      //!< Automatically calculates number of modules and adds this number to it
   Property<        bool  , Default>   m_alignEdges;             //!< Start building modules in the ring with its edge or centre
   Property<        double, Default>   m_ringGap;                //!< Required gap in-between rings
   Property<        int   , Default>   m_smallParity;            //!< In each rod the modules are positioned in R-Phi zig-zag at average zPos +- smallDelta -> use +- or -+smallDelta

   double m_averageZ = 0;      //!< Average Z position of the given ring

   const int c_maxWedgeCalcLoops = 100; //!< Maximum number of loops to be iterated when finding optimal positions for wege-shaped modules

};

#endif /* INCLUDE_RING_H */

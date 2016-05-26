#ifndef INCLUDE_RODPAIR_H_
#define INCLUDE_RODPAIR_H_

#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <algorithm>

#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.h"
#include "GeometryFactory.h"
#include "Property.h"
#include "DetectorModule.h"
#include "MessageLogger.h"
#include "Visitable.h"

using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;

typedef PtrVector<BarrelModule>          BarrelModules;

enum class BuildDir   { RIGHT = 1, LEFT = -1 };
enum class StartZMode { MODULECENTER, MODULEEDGE };

/*
 * @class RodPair
 * @brief Base rod class -> use RodStraight or RodTilted
 * @details Rod class holds information about individual layer rods (pairs of pozitive & negative modules, hence RodPair).
 * It's building procedure is executed automatically via Layer class. Similarly, all its components: positive & negative
 * modules in Z are recursively build through private build() method.
 */
class RodPair : public PropertyObject, public Buildable, public Identifiable<int>, public Visitable {

 public:

  //!  Constructor - parse geometry config file using boost property tree & read-in Rod parameters -> use number of modules to build rod
  RodPair(int id, double radius, double rotation, int numModules, const PropertyTree& treeProperty);

  //!  Constructor - parse geometry config file using boost property tree & read-in Rod parameters -> use outerZ to build rod
  RodPair(int id, double radius, double rotation, double outerZ , const PropertyTree& treeProperty);

  //! Position newly individual modules if RodPair cloned from a RodPair (i.e. rotate by respective angle and shift in R) and set rod new id
   void buildClone(int id, double shiftR, double rotation);

  //! Setup: link lambda functions to various rod related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Limit rod geometry by eta cut
  void cutAtEta(double eta);

  // Return rod modules as a pair of modules on positive Z (first) & negative Z (second)
  const std::pair<const BarrelModules&,const BarrelModules&> modules() const { return std::pair<const BarrelModules&,const BarrelModules&>(m_zPlusModules,m_zMinusModules); }

  //! Return rod as a material object
  const MaterialObject& materialObject() const {return m_materialObject; }

  //! GeometryVisitor pattern -> rod visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> rod visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  //! Get total number of modules built
  int numModules() const { return m_zPlusModules.size() + m_zMinusModules.size(); }

  //! Get number of positive (side>0) or negative (side<0) of a rod
  int numModulesSide(int side) const { return side >= 0 ? m_zPlusModules.size() : m_zMinusModules.size(); }

  //! Get rod thickness - purely virtual method
  virtual double thickness() const = 0;

  ReadonlyProperty<double, Computable> maxZ;               //!< Maximum rod Z
  ReadonlyProperty<double, Computable> minZ;               //!< Minimum rod Z
  ReadonlyProperty<double, Computable> minR;               //!< Minimum rod radius
  ReadonlyProperty<double, Computable> maxR;               //!< Maximum rod radius
  ReadonlyProperty<double, Computable> maxModuleThickness; //!< Calculated maximum module thickness in a rod
  Property<        bool  , Default>    beamSpotCover;      //!< Take into account beam spot size when positioning modules

 protected:

  BarrelModules m_zPlusModules, m_zMinusModules; //!< Rod modules in positive/negative Z
  MaterialObject m_materialObject;

  Property<StartZMode, Default> m_startZMode;  //!< Start positioning first Z module at Z=0 centrally or with left edge at Z=0 (first to positive Z, then to negative Z)

  int    m_nModules;      //!< Number of modules to be built in a rod of use outerZ if not defined (passed over from layer)
  double m_outerZ;        //!< If number of modules not defined use rod outerZ value (passed over from layer)
  double m_optimalRadius; //!< Optimal rod radius
  double m_rotation;      //!< Rod rotation in R-Phi

}; // Class

/*
 * @class Straight option of rod class
 * @details Derived class from RodPair implementing straight option of layer rod
 */
class RodPairStraight : public RodPair, public Clonable<RodPairStraight> {

public:

  //! Constructor - parse geometry config file using boost property tree & read-in Rod parameters -> use numModules to build rod
  RodPairStraight(int id, double radius, double rotation, double bigDelta, int bigParity, double smallDelta, int smallParity, int numModules, const PropertyTree& treeProperty);

  //! Constructor - parse geometry config file using boost property tree & read-in Rod parameters -> use outerZ to build rod
  RodPairStraight(int id, double radius, double rotation, double bigDelta, int bigParity, double smallDelta, int smallParity, double outerZ , const PropertyTree& treeProperty);

  //! Build recursively individual modules
  void build();

  //! Get rod thickness
  double thickness() const override { return m_smallDelta*2. + maxModuleThickness(); }

  RangeProperty<std::vector<double> > forbiddenRange;

  Property<double, Default>   zOverlap;            //!< Required overlap of modules in Z (in length units)
  Property<double, NoDefault> zError;              //!< When positioning modules take into account beam spot spread in Z
  Property<bool  , Default>   compressed;          //!< Modules will be compressed in Z if built layer higher than defined outerZ parameter (if number of modules was used to defined  the rod, no compression occurs)
  Property<bool  , Default>   allowCompressionCuts;//!< During compression algorithm cut out modules behind outerZ first

 private:

  //! Helper method calculating module optimal Z position taking into account small/bigDelta, beam spot sizes etc. -> to fully cover eta region
  double computeNextZ(double lastRPos, double newRPos, double lastZPos, double moduleLength, double lastDsDistance, double newDsDistance, BuildDir direction, int parity);

  //! Helper method compressing modules in Z direction (symmetrically from minus/plus Z) such as they fit into -outerZ, +outerZ region
  //! The algorithm first cuts out modules completely exceeding outerZ limit. Later on calculates (+ or -) shift of the last module with
  //! respect to the outerZ position. Finally each module is shifted by a fraction given as shift/moduleCentreZ. To avoid clashes of
  //! of newly positioned modules after their shift, cross-check between neighbouring modules is done.
  void compressToZ(double z);

  PropertyNode<int> m_ringNode;

  static const int   c_nIterations;             //!< Number of iterations allowed when positioning or compressing modules
  static const float c_safetySpaceFactor;       //!< Safety factor used when compressing modules in Z

  double m_smallDelta;     //!< Layer consists of ladders (rods), in which modules are positioned at radius +- smallDelta in Z
  int    m_smallParity;    //!< Algorithm that builds rod modules starts at +smallDelta (positive parity) or -smallDelta (negative parity)
  double m_bigDelta;       //!< Layer consists of ladders (rods), where even/odd rods are positioned at radius +- bigDelta in R-Phi
  int    m_bigParity;      //!< Algorithm that builds rods starts at +bigDelta (positive parity) or -bigDelta (negative parity)

}; // Class
//
//
//struct TiltedModuleSpecs {
//  double r, z, gamma;
//  bool valid() const {
//    return r > 0.0 && fabs(gamma) <= 2*M_PI;
//  }
//};
//
///*
// * @class Tilted option of rod class
// * @details Derived class from RodPair implementing tilted option verison of layer rod
// */
//class TiltedRodPair : public RodPair, public Clonable<TiltedRodPair> {
//  void buildModules(BarrelModules& modules, const RodTemplate& rodTemplate, const vector<TiltedModuleSpecs>& tmspecs, BuildDir direction);
//public:
//  double thickness() const override { std::cerr << "thickness() for tilted rods gives incorrect results as it is calculated as maxR()-minR()\n"; return maxR() - minR(); }
//  void build(const RodTemplate& rodTemplate, const std::vector<TiltedModuleSpecs>& tmspecs);
//
//};

#endif /* INCLUDE_RODPAIR_H_ */

#ifndef INCLUDE_LAYER_H_
#define INCLUDE_LAYER_H_

#include <memory>
#include <vector>
#include <string>
#include <limits.h>

#include "math_functions.h"
#include "ConversionStation.h"
#include "Property.h"
#include "RodPair.h"
#include "Ring.h"
#include "Visitable.h"
#include "MaterialObject.h"

// Used namespace in following classes
using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;
using material::MaterialObject;
using material::ConversionStation;
using material::ConversionStations;

// Typedefs
typedef PtrVector<Layer>      Layers;

enum RadiusMode { SHRINK, ENLARGE, FIXED, AUTO };

/*
 * @class Layer
 * @details Layer class holds information about individual barrel layers. It's building procedure is executed automatically via
 * Barrel class. Similarly, all its components (layer rods (ladders) -> modules) are recursively called through build()
 * method. Depending on the configuration, either straight layers are built using buildStraight() method or tilted using
 * buildTilted() method.
 */
class Layer : public PropertyObject, public Buildable, public Identifiable<int>, public Clonable<Layer>, public Visitable {

 public:

  //! Constructor - parse geometry config file using boost property tree & read-in Layer parameters
  Layer(int id, int barrelNumLayers, bool sameRods, bool barrelMinRFixed, bool barrelMaxRFixed, double barrelRotation,
        const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty);


  //! Build recursively individual subdetector systems: rods -> modules & conversion stations
  void build(int barrelNumLayers, double barrelMinR, double barrelMaxR);

  //! Setup: link lambda functions to various layer related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Limit layer geometry by eta cut
  void cutAtEta(double eta);

  //! Return layer straight rods
  const Rods& flatRods()   const { return m_flatRods; }

  //! Return layer tilted rods
  const Rods& tiltedRods() const { return m_tiltedRods; }

  //! Return tilted ring
  const TiltedRing& tiltedRing(int iRing) const { if (iRing<0 || iRing>=numTiltedRings()) throw std::invalid_argument("Layer::tiltedRing - Ring: "+any2str(iRing)+" out of range: "+any2str(m_tiltedRings.size()-1)+"!!!");
                                                  else                                    return m_tiltedRings.at(iRing); }

  //! Return layer as a material object
  const MaterialObject& materialObject() const { return m_materialObject; }

  //! Return first order conversion station -> TODO: check if needed to be updatable, use PtrVector instead
  ConversionStation* flangeConversionStation() const {return m_flangeConversionStation; }

  //! Return vector of second order conversion stations -> TODO: check if needed to be updatable, use reference instead
  ConversionStations secondConversionStations() const {return m_secondConversionStations; }

  //! GeometryVisitor pattern -> layer visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> layer visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  //! Get number of rods in a layer (by engineering, tilted part must contain the same number of rods as the flat part)
  int numRods() const { return m_flatRods.size(); }

  //! Get number of modules per rod
  int numModulesPerRod() const { return m_flatRods.front().numModules()+(m_tiltedRods.size()>0 ? m_tiltedRods.front().numModules() : 0); }

  //! Get number of modules per rod side (side>0 <=> z+ side, side<0 <=> z- side)
  int numModulesPerRodSide(int side) const { return m_flatRods.front().numModulesSide(side)+(m_tiltedRods.size()>0 ? m_tiltedRods.front().numModulesSide(side) : 0); }

  //! Get total number of modules in a layer
  int totalModules() const { return numModulesPerRod()*numRods(); };

  //! Get number of tilted rings
  int numTiltedRings() const { return m_tiltedRings.size(); };

  //! Get number of flat rings
  //int numFlatRings() const { return buildNumModulesFlat(); };

  //! Get flat rod thickness
  double flatRodThickness() const { return m_flatRods.front().thickness(); }

  //! Get tilted rod thickness
  double tiltedRodThickness() const { return m_tiltedRods.front().thickness(); }

  //! Does flat part has a central module -> build procedure started with module positioned at Z=0 (not by edge)
  bool hasStraightCentralModule() const;

  ReadonlyProperty<double    , Computable>  minFlatZ;           //!< Flat: Minimum layer Z
  ReadonlyProperty<double    , Computable>  maxFlatZ;           //!< Flat: Maximum layer Z
  ReadonlyProperty<double    , Computable>  minFlatR;           //!< Flat: Minimum layer radius (given by different positions of layer modules - big & smallDelta)
  ReadonlyProperty<double    , Computable>  maxFlatR;           //!< Flat: Maximum layer radius (given by different positions of layer modules - big & smallDelta)
  ReadonlyProperty<double    , Computable>  minFlatRAllMat;     //!< Flat: Minimum layer radius taking into account all material structures
  ReadonlyProperty<double    , Computable>  maxFlatRAllMat;     //!< Flat: Maximum layer radius taking into account all material structures
  ReadonlyProperty<double    , Computable>  minTiltedZ;         //!< Tilted: Minimum layer Z
  ReadonlyProperty<double    , Computable>  maxTiltedZ;         //!< Tilted: Maximum layer Z
  ReadonlyProperty<double    , Computable>  minTiltedR;         //!< Tilted: Minimum layer radius (given by different positions of layer modules - big & smallDelta)
  ReadonlyProperty<double    , Computable>  maxTiltedR;         //!< Tilted: Maximum layer radius (given by different positions of layer modules - big & smallDelta)
  ReadonlyProperty<double    , Computable>  minTiltedRAllMat;   //!< Tilted: Minimum layer radius taking into account all material structures
  ReadonlyProperty<double    , Computable>  maxTiltedRAllMat;   //!< Tilted: Maximum layer radius taking into account all material structures
  ReadonlyProperty<double    , Computable>  minZ;               //!< Minimum layer Z
  ReadonlyProperty<double    , Computable>  maxZ;               //!< Maximum layer Z
  ReadonlyProperty<double    , Computable>  minR;               //!< Minimum layer radius (given by different positions of layer modules - big & smallDelta)
  ReadonlyProperty<double    , Computable>  maxR;               //!< Maximum layer radius (given by different positions of layer modules - big & smallDelta)
  ReadonlyProperty<double    , Computable>  minRAllMat;         //!< Minimum layer radius taking into account all material structures
  ReadonlyProperty<double    , Computable>  maxRAllMat;         //!< Maximum layer radius taking into account all material structures
  Property<        double    , NoDefault>   outerZ;             //!< Outer Z position (maximum Z), if buildNumModules not defined, this variable used instead to define layer halfLength
  Property<        int       , AutoDefault> buildNumModules;    //!< Number of modules to be built in each layer, if not defined outerZ variable used instead (default 0)

  //! Defines which algorithm used to optimize layer radius for straight (flat) option: F (fixed), A (auto), E (expand) or S (shrink).
  //! F stands for the positioning at a user-defined radii, i.e. fixed radii.
  //! E stands for an expansion, i.e. optimal radius is bigger then the current one and detectors are positioned to an expanded radius.
  //! S stands for shrinking, i.e. optimal position is lower then the current one and detectors are positioned to a shrunk radius.
  //! A stands for an automatic mode, i.e. the closer configuration is chosen, either E or S. That's the implicit build mode.
  //! In addition, modules are positioned such as all lines (at various eta) going from the primary vertex, defined as (0, 0, +-zError), are always passing through edges of layer modules.
  //! The positioning algorithm starts at Z=0 and takes then into consideration the extreme cases, taking into account all parameters: bigDelta, smallDelta, zError, z/phiOverlap
  Property<RadiusMode, Default>     radiusMode;
  Property<double    , NoDefault>   requestedAvgRadius;   //!< Requested radius at which the layer should be positioned
  Property<double    , Computable>  avgBuildRadius;       //!< Average layer radius (central value) calculated based on position algorithm
  Property<double    , Computable>  avgBuildRadiusFlat;   //!< Average layer radius for flat part only
  Property<double    , Computable>  avgBuildRadiusTilted; //!< Average layer radius for tilted part only
  Property<bool      , Default>     sameParityRods;       //!< When starting to build even/odd rods use the same (not opposite) smallDelta parity (not used for tilted -> smallParity is forced so that central module is at +smallDelta)
  Property<double    , Default>     layerRotation;        //!< Layer rotated by general barrel rotation + this value in R-Phi

  Property<bool      , Default>     isTilted;             //!< Does the layer consists of tilted modules
  Property<bool      , NoDefault>   isTiltedAuto;         //!< Is modules tilt calculated automatically, using internal algorithm? If not, use configuration file -> FILE NOT SUPPORTED ANYMORE
  Property<string    , AutoDefault> tiltedLayerSpecFile;  //!< Configuration file for tilted option -> FILE NOT SUPPORTED ANYMORE
  Property<int       , AutoDefault> buildNumModulesFlat;  //!< A number of modules to be build in the straight part (tilt=0) of each layer rod, i.e. number of rings in flat part (default 0)
  Property<int       , AutoDefault> buildNumModulesTilted;//!< A number of modules to be built in the tilted part of each layer rod, i.e. number of rings in tilted part (default 0)

  Property<double    , NoDefault>   phiOverlap;           //!< Required module overlap in R-Phi (in length units) -> used if flat part built only
  Property<int       , NoDefault>   phiSegments;          //!< Required symmetry in R-Phi - number of symmetric segments (1, 2, 4, ...) -> used if flat part built only
  Property<double    , NoDefault>   bigDelta;             //!< Layer consists of ladders (rods), where even/odd rods are positioned at radius +- bigDelta in R-Phi
  Property<double    , NoDefault>   smallDelta;           //!< Layer consists of ladders (rods), in which modules are positioned at radius +- smallDelta in Z
  Property<int       , Default>     m_bigParity;          //!< Algorithm that builds rods starts at +bigDelta (positive parity) or -bigDelta (negative parity)

  Property<int       , NoDefault>   numberRods;           //!< If tilted part being built, number of rods in phi is not calculated automatically, but must be defined

  Property<double    , Computable>  flatPhiOverlapSmallDeltaMinus; //!< Overlap in R-phi (in length units) between modules in 2 adjacent rods at +smallDelta
  Property<double    , Computable>  flatPhiOverlapSmallDeltaPlus;  //!< Overlap in R-phi (in length units) between modules in 2 adjacent rods at -smallDelta
  Property<double    , Computable>  flatThetaEnd;                  //!< Theta angle at which the flat part ends-up -> important entry parameter for the tilted part
  Property<double    , Computable>  flatREndInner;                 //!< Radial position of last module in inner rod of flat part
  Property<double    , Computable>  flatREndOuter;                 //!< Radial position of last module in outer rod of flat part
  Property<double    , Computable>  flatEndModCentreZ;             //!< Z position of last module of flat part: module center
  Property<double    , Computable>  flatZEnd;                      //!< Z position at which the flat part ends-up

  //!< Beam spot coverage in Z by flat part inner rod -> only active module parameters taken into account
  double zErrorInnerFlat(int iRing)   const { if (iRing<0 || iRing>=buildNumModulesFlat())   throw std::invalid_argument("Layer::zErrorInnerFlat - Ring: "+any2str(iRing)+" out of range!!!"); else return m_zErrorInnerFlat.at(iRing); };

  //!< Beam spot coverage in Z by flat part outer rod -> only active module parameters taken into account
  double zErrorOuterFlat(int iRing)   const { if (iRing<0 || iRing>=buildNumModulesFlat())   throw std::invalid_argument("Layer::zErrorOuterFlat - Ring: "+any2str(iRing)+" out of range!!!"); else return m_zErrorOuterFlat.at(iRing);   };

 private:

  //! Cross-check parameters provided from geometry configuration file
  void check() override;

  //! If straight layer required, build() method internally calls buildStraight()
  void buildStraight(int barrelNumLayers, double barrelMinR, double barrelMaxR);

  //! If tilted layer required, build() method internally calls buildTilted()
  void buildTilted();

  //! Helper function calculating optimal layer radius for straight option
  double calculateOptimalRadius(int numRods, double bigDelta, double smallDelta, double dsDistance, double moduleWidth, double overlap);

  //! Calculate geometry properties of straight part needed to build the tilted part
  void calculateStraightProperties();

  Rods               m_flatRods;                  //!< Layer non-tilted rods
  Rods               m_tiltedRods;                //!< Layer tilted rods (zero if straight layers built only)
  TiltedRings        m_tiltedRings;               //!< Rings, in which the tilted modules will be arranged in a tilted rod
  MaterialObject     m_materialObject;
  ConversionStation* m_flangeConversionStation;   //!< First order layer conversion unit
  ConversionStations m_secondConversionStations;  //!< Vector of second order layer conversion units

  std::map<int, double> m_zErrorInnerFlat;
  std::map<int, double> m_zErrorOuterFlat;

  Property<int   , Default>   m_smallParity;      //!< Algorithm that builds rod modules starts at +smallDelta (positive parity) or -smallDelta (negative parity)
  Property<bool  , Default>   m_useMinMaxRCorrect;//!< Apply smallDelta, bigDelta, detThickness, etc. when calculating minR/maxR for first/last layer? For backwards compatibility of lite version (on) versus older version (off)

  bool   m_sameRods;                              //! Build same geometrical rods across the whole flat barrel

  PropertyNode<int>               m_ringNode;     //!< Property tree node for ring (to grab properties for specific rod modules set in rings)
  PropertyNodeUnique<std::string> m_stationsNode; //!< Property tree nodes for conversion stations (included geometry config file)

}; // Class

#endif /* INCLUDE_LAYER_H_ */

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

typedef std::map<int, TiltedRing*> TiltedRingsTemplate;

class FlatRingsGeometryInfo {
 private:
  std::map<int, double> zErrorInner_;
  std::map<int, double> zErrorOuter_;
 public:
  FlatRingsGeometryInfo();
  void calculateFlatRingsGeometryInfo(std::vector<RodPairStraight*> flatPartRods, double bigParity);
  std::map<int, double> zErrorInner() const { return zErrorInner_; }
  std::map<int, double> zErrorOuter() const { return zErrorOuter_; }
};

// Typedefs
typedef PtrVector<Layer> Layers;

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

  //! Return layer rods
  const Rods& rods() const { return m_rods; }

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

  //! Get number of rods in a layer
  //int numRods() const { return m_rods.size(); }

  //! Get number of modules per rod
  int numModulesPerRod() const { return m_rods.front().numModules(); }

  //! Get number of modules per rod side (side>0 <=> z+ side, side<0 <=> z- side)
  int numModulesPerRodSide(int side) const { return m_rods.front().numModulesSide(side); }

  //! Get total number of modules in a layer
  int totalModules() const { return numModulesPerRod()*numRods(); };

  //! Get rod thickness
  double rodThickness() const { return m_rods.front().thickness(); }

  ReadonlyProperty<double    , Computable>  minZ;               //!< Minimum layer Z
  ReadonlyProperty<double    , Computable>  maxZ;               //!< Maximum layer Z
  ReadonlyProperty<double    , Computable>  minR;               //!< Minimum layer radius (given by different positions of layer modules - big & smallDelta)
  ReadonlyProperty<double    , Computable>  maxR;               //!< Maximum layer radius (given by different positions of layer modules - big & smallDelta)
  ReadonlyProperty<double    , Computable>  minRAllMat;         //!< Minimum layer radius taking into account all material structures
  ReadonlyProperty<double    , Computable>  maxRAllMat;         //!< Maximum layer radius taking into account all material structures
  Property<        double    , NoDefault>   outerZ;             //!< Outer Z position (maximum Z), if buildNumModules not defined, this variable used instead to define layer halfLength
  Property<        int       , AutoDefault> buildNumModules;    //!< Number of modules to be built in each layer, if not defined outerZ variable used instead

  //! Defines which algorithm used to optimize layer radius for straight option: F (fixed), A (auto), E (expand) or S (shrink).
  //! F stands for the positioning at a user-defined radii, i.e. fixed radii.
  //! E stands for an expansion, i.e. optimal radius is bigger then the current one and detectors are positioned to an expanded radius.
  //! S stands for shrinking, i.e. optimal position is lower then the current one and detectors are positioned to a shrunk radius.
  //! A stands for an automatic mode, i.e. the closer configuration is chosen, either E or S. That's the implicit build mode.
  //! In addition, modules are positioned such as all lines (at various eta) going from the primary vertex, defined as (0, 0, +-zError), are always passing through edges of layer modules.
  //! The positioning algorithm starts at Z=0 and takes then into consideration the extreme cases, taking into account all parameters: bigDelta, smallDelta, zError, z/phiOverlap
  Property<RadiusMode, Default>     radiusMode;
  Property<double    , NoDefault>   requestedAvgRadius;  //!< Requested radius at which the layer should be positioned
  Property<double    , Computable>  avgBuildRadius;      //!< Average layer radius (central value) calculated based on position algorithm
  Property<bool      , Default>     sameParityRods;      //!< When starting to build even/odd rods use the same (not opposite) smallDelta parity
  Property<double    , Default>     layerRotation;       //!< Layer rotated by general barrel rotation + this value in R-Phi
  Property<string    , AutoDefault> tiltedLayerSpecFile; //!< Configuration file for tilted option

  Property<double    , Default>     phiOverlap;       //!< Required module overlap in R-Phi (in length units)
  Property<int       , Default>     phiSegments;      //!< Required symmetry in R-Phi - number of symmetric segments (1, 2, 4, ...)
  Property<double    , NoDefault>   bigDelta;         //!< Layer consists of ladders (rods), where even/odd rods are positioned at radius +- bigDelta in R-Phi
  Property<double    , NoDefault>   smallDelta;       //!< Layer consists of ladders (rods), in which modules are positioned at radius +- smallDelta in Z

  Property<int, NoDefault> numRods;

  Property<int, NoDefault> buildNumModulesFlat;
  Property<int, NoDefault> buildNumModulesTilted;
  Property<bool, Default> isTilted;
  Property<bool, NoDefault> isTiltedAuto;


 private:

  class TiltedRingsGeometryInfo {
  private:
    std::map<int, double> deltaZInner_;
    std::map<int, double> deltaZOuter_;
    //std::map<int, double> covInner_;
    std::map<int, double> zErrorInner_;
    std::map<int, double> zErrorOuter_;
  public:
    TiltedRingsGeometryInfo(int numModulesFlat, double, double, double, double, TiltedRingsTemplate tiltedRingsGeometry);
    std::map<int, double> deltaZInner() const { return deltaZInner_; }
    std::map<int, double> deltaZOuter() const { return deltaZOuter_; }
    //std::map<int, double> covInner() const { return covInner_; }
    std::map<int, double> zErrorInner() const { return zErrorInner_; }
    std::map<int, double> zErrorOuter() const { return zErrorOuter_; }
  };

  //! Cross-check parameters provided from geometry configuration file
  void check() override;

  //! If straight layer required, build() method internally calls buildStraight()
  void buildStraight(int barrelNumLayers, double barrelMinR, double barrelMaxR);

  //! If tilted layer required, build() method internally calls buildTilted()
  void buildTilted();

  //! Helper function calculating optimal layer radius for straight option
  double calculateOptimalRadius(int numRods, double bigDelta, double smallDelta, double dsDistance, double moduleWidth, double overlap);

  Rods               m_rods;                      //!< Layer rods
  MaterialObject     m_materialObject;
  ConversionStation* m_flangeConversionStation;   //!< First order layer conversion unit
  ConversionStations m_secondConversionStations;  //!< Vector of second order layer conversion units
  RodTemplate makeRodTemplate();

  Rods m_flat_rods;
  double flatPartPhiOverlapSmallDeltaMinus_;
  double flatPartPhiOverlapSmallDeltaPlus_;
  double flatPartAverageR_;
  FlatRingsGeometryInfo flatRingsGeometryInfo_;
  TiltedRingsTemplate tiltedRingsGeometry_;
  TiltedRingsGeometryInfo tiltedRingsGeometryInfo_ = TiltedRingsGeometryInfo(0,0,0,0,0, tiltedRingsGeometry_);
  TiltedRingsTemplate makeTiltedRingsTemplate(double flatPartThetaEnd);

  Property<int   , Default>   m_smallParity;      //!< Algorithm that builds rod modules starts at +smallDelta (positive parity) or -smallDelta (negative parity)
  Property<int   , Default>   m_bigParity;        //!< Algorithm that builds rods starts at +bigDelta (positive parity) or -bigDelta (negative parity)
  Property<bool  , Default>   m_useMinMaxRCorrect;//!< Apply smallDelta, bigDelta, detThickness, etc. when calculating minR/maxR for first/last layer? For backwards compatibility of lite version (on) versus older version (off)

  bool   m_sameRods;                              //! Build same geometrical rods across the whole barrel

  PropertyNode<int>               m_ringNode;     //!< Property tree node for ring (to grab properties for specific rod modules)
  PropertyNodeUnique<std::string> m_stationsNode; //!< Property tree nodes for conversion stations (included geometry config file)

}; // Class

#endif /* INCLUDE_LAYER_H_ */

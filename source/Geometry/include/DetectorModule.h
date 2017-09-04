#ifndef INCLUDE_DETECTOR_MODULE_H_
#define INCLUDE_DETECTOR_MODULE_H_

#include <boost/ptr_container/ptr_vector.hpp>

#include "GeometricModule.h"
#include "GeometryFactory.h"
#include "MaterialObject.h"
#include "MaterialProperties.h"
#include "math_functions.h"
#include "MessageLogger.h"
#include "Property.h"
#include "Sensor.h"
#include "Visitable.h"
#include "Visitor.h"

using material::MaterialObject;

// Enum definitions related to DetectorModule
enum ModuleSubdetector { BARREL = 1, ENDCAP = 2 };
enum SensorLayout      { NOSENSORS, MONO, PT, STEREO };
enum ZCorrelation      { SAMESEGMENT, MULTISEGMENT };
enum ReadoutType       { READOUT_STRIP, READOUT_PIXEL, READOUT_PT };
enum ReadoutMode       { BINARY, CLUSTER };
enum HitType           { NONE, INNER, OUTER, BOTH = 3, STUB = 7 };

// Typedefs
typedef PtrVector<BarrelModule> BarrelModules;
typedef PtrVector<EndcapModule> EndcapModules;

struct PosRef { int cnt, z, rho, phi; };
struct TableRef { string table; int row, col; };
struct UniRef { string cnt; int layer, ring, phi, side; };

// Forward declaration
class ModuleCap;
class DetectorModule;

using material::ElementsVector;

/*
 * @class BaseModule
 * @brief Module base class to define general cloning functionality for base DetectorModule class and all its derived classes
 */
class BaseModule
{
public:

  //! Copy constructor -> needed for cloning functionalities
  BaseModule(const BaseModule& module) : m_moduleGeom(module.m_moduleGeom->clone()) {}

protected:

  //! Constructor
  BaseModule(GeometricModule* const moduleGeom) : m_moduleGeom(moduleGeom) {}

  //!< Geometrical properties of module, so called geometrical module
  GeometricModule* const m_moduleGeom;
};

/*
 * @class DetectorModule
 * @brief Detector module base class - base to Barrel & Endcap classes
 * @details Implementors of the DetectorModuleInterface must take care of rotating the module based on which part of the subdetector it will be used in (Barrel, EC)
 */
class DetectorModule : public BaseModule, public PropertyObject, public Buildable, public Placeable, public Identifiable<int>, public Visitable
{

public:

  //! Constructor - specify unique id, geometry module defining shape & parse geometry config file using boost property tree & read-in module parameters
  DetectorModule(int id, GeometricModule* moduleGeom, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty);

  //! Constructor - specify unique id, geometry module defining shape & parse geometry config file using boost property tree & read-in module parameters
  DetectorModule(int id, GeometricModule* moduleGeom, const PropertyTree& treeProperty);

  //! Destructor
  ~DetectorModule();

  //! Build detector module with its sensors
  virtual void build();

  //! Setup: link lambda functions to various DetectorModule related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  virtual void setup();

  //! Set materials, so called module-cap to a module
  void setModuleCap(ModuleCap* newCap);

  //! Get materials, i.e. module-cap to be modified
  ModuleCap* getModuleCap() { return m_moduleCap;}

  //! Get materials, i.e. module-cap to be read-only
  const ModuleCap& getModuleCap() const { return *m_moduleCap; }

  //! Geometric module interface -> return reference to a geometrical representation of a module
  const Polygon3D<4>& basePoly() const { return m_moduleGeom->basePoly(); }

  //! Could track at this direction hit the module? (fast method)
  bool couldHit(const XYZVector& direction, double zError) const;

  //! Check if track hit the module -> if yes, return true with passed material, hit position vector & hitType (which module sensor(s) was/were hit)
  bool checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir, Material& hitMaterial, HitType& hitType, XYZVector& hitPos) const;

  //! Calculate min & max module eta coverage taking into account error on beam Z position
  std::pair<double, double> minMaxEtaWithError(double zError) const;

  //! Derive classes use it to identify themselves
  virtual ModuleSubdetector subdet() const = 0;

  //! Return module sensors
  const Sensors& sensors()               const { return m_sensors; }
  const Sensor& innerSensor()            const { return m_sensors.front(); }
  const Sensor& outerSensor()            const { return m_sensors.back(); }

  //! Return module as a material object
  const MaterialObject& materialObject() const { return m_materialObject; }

  //! Get material composition
  ElementsVector& getLocalElements() const {return m_materialObject.getLocalElements(); }

  //! Return module shape (rectangular/wedge)
  ModuleShape shape() const { return m_moduleGeom->shape(); }

  //! Reset number of hits assigned to a module
  void resetHits() { m_numHits = 0; }

  // Get geometry properties
  const XYZVector& center() const { return m_moduleGeom->center(); }
  const XYZVector& normal() const { return m_moduleGeom->normal(); }
  double area()             const { return m_moduleGeom->area(); }
  double dsDistance()       const { return m_moduleGeom->dsDistance(); }

  double thickness()        const { return dsDistance() + sensorThickness(); }
  double thicknessAllMat()  const { double modThick = 0; for (auto& elem : m_materialObject.getLocalElements()) modThick += elem->quantity(); return modThick; }
  double flatMaxZ()         const { return (center().Z() + length()/2.); }
  double length()           const { return m_moduleGeom->length(); }
  double maxWidth()         const { return m_moduleGeom->maxWidth(); }
  double minWidth()         const { return m_moduleGeom->minWidth(); }
  double meanWidth()        const { return m_moduleGeom->meanWidth(); }
  double physicalLength()   const { return m_moduleGeom->physicalLength(); }
  double tiltAngle()        const { return m_moduleGeom->tiltAngle(); }
  double skewAngle()        const { return m_moduleGeom->skewAngle(); }

  double maxEta()           const { return MAX(basePoly().getVertex(0).Eta(), basePoly().getVertex(2).Eta()); }
  double minEta()           const { return MIN(basePoly().getVertex(0).Eta(), basePoly().getVertex(2).Eta()); }

  double maxTheta()         const { return MAX(basePoly().getVertex(0).Theta(), basePoly().getVertex(2).Theta()); }
  double minTheta()         const { return MIN(basePoly().getVertex(0).Theta(), basePoly().getVertex(2).Theta()); }

  double phiAperture()      const { return maxPhi() - minPhi(); }
  double etaAperture()      const { return maxEta() - minEta(); }
  double thetaAperture()    const { return maxTheta() - minTheta(); }

  // Get other properties
  int    numHits()          const { return m_numHits; }
  double totalPowerModule() const { return powerModuleOptical() + powerModuleChip(); }
  double totalPowerStrip()  const { return powerStripOptical() + powerStripChip(); }
  double totalPower()       const { return totalPowerModule() + totalPowerStrip()*outerSensor().numChannels(); }

  int    numStripsAcross()  const { return sensors().front().numStripsAcross(); } // CUIDADO this assumes both sensors have the same number of sensing elements in the transversal direction - typically it is like that
  double sensorThickness()  const { return sensors().front().sensorThickness(); } // CUIDADO this has to be fixed (called in Extractor.cc), sensor thickness can be different for different sensors

  double stripOccupancyPerEvent() const;
  double hitOccupancyPerEvent()   const { return stripOccupancyPerEvent()/2.; }
  double geometricEfficiency()    const;
  double effectiveDsDistance()    const;

  // Set geometry properties
  void translate(const XYZVector& vector) { m_moduleGeom->translate(vector); clearSensorPolys(); }
  void mirror(   const XYZVector& vector) { m_moduleGeom->mirror(vector); clearSensorPolys(); }
  void translateZ(double z)               { m_moduleGeom->translate(XYZVector(0, 0, z)); clearSensorPolys(); }
  void translateR(double radius)          { XYZVector v = m_rAxis.Unit()*radius; m_moduleGeom->translate(v); clearSensorPolys();}
  void rotateToNegativeZSide();
  void rotateX(double angle)              { m_moduleGeom->rotateX(angle); clearSensorPolys(); }
  void rotateY(double angle)              { m_moduleGeom->rotateY(angle); clearSensorPolys(); }
  void rotateZ(double angle)              { m_moduleGeom->rotateZ(angle); clearSensorPolys(); m_rAxis = ROOT::Math::RotationZ(angle)(m_rAxis); }
  void tilt(double angle)                 { rotateX(-angle); m_moduleGeom->tiltAngle(angle); } // Atention: Tilt can only be called BEFORE translating/rotating the module, or they won't work as expected!!
  void skew(double angle)                 { rotateY(-angle); m_moduleGeom->skewAngle(angle); } // Atention: Skew can only be called BEFORE translating/rotating the module, or they won't work as expected!!
  bool flipped() const { return m_moduleGeom->flipped(); } 
  bool flipped(bool newFlip) {
    if (newFlip && numSensors() > 1) {
      m_sensors.front().innerOuter(SensorPosition::UPPER);
      m_sensors.back().innerOuter(SensorPosition::LOWER);
    }
    return m_moduleGeom->flipped(newFlip);
  } 
  void dsDistance(double d)               { m_moduleGeom->dsDistance(d); }

  // Properties
  ReadonlyProperty<double    , UncachedComputable> minZ;       //!< Minimum module Z (use uncached -> calculate its value anytime when needed -> need to be uncached recursively, i.e. also for Sensors)
  ReadonlyProperty<double    , UncachedComputable> maxZ;       //!< Maximum module Z
  ReadonlyProperty<double    , UncachedComputable> minZAllMat; //!< Minimum module Z taking into account all material structures
  ReadonlyProperty<double    , UncachedComputable> maxZAllMat; //!< Maximum module Z taking into account all material structures
  ReadonlyProperty<double    , UncachedComputable> minR;       //!< Minimum module radius
  ReadonlyProperty<double    , UncachedComputable> maxR;       //!< Maximum module radius
  ReadonlyProperty<double    , UncachedComputable> minRAllMat; //!< Minimum module radius taking into account all material structures
  ReadonlyProperty<double    , UncachedComputable> maxRAllMat; //!< Maximum module radius taking into account all material structures
  ReadonlyProperty<double    , UncachedComputable> planarMinZ; //!< Minimum module Z, if module approximated as of zero thickness
  ReadonlyProperty<double    , UncachedComputable> planarMaxZ; //!< Maximum module Z, if module approximated as of zero thickness
  ReadonlyProperty<double    , UncachedComputable> planarMinR; //!< Minimum module R, if module approximated as of zero thickness
  ReadonlyProperty<double    , UncachedComputable> planarMaxR; //!< Maximum module R, if module approximated as of zero thickness
  Property<        double    , Computable        > minPhi;     //!< Module corner at the lowest phi coordinate
  Property<        double    , Computable        > maxPhi;     //!< Module corner at the highest phi coordinate
  ReadonlyProperty<int       , Computable        > maxSegments;
  ReadonlyProperty<int       , Computable        > minSegments;
  ReadonlyProperty<int       , Computable        > totalSegments;
  ReadonlyProperty<int       , Computable        > maxChannels;
  ReadonlyProperty<int       , Computable        > minChannels;
  ReadonlyProperty<int       , Computable        > totalChannels;

  PropertyVector<  string    , ','>           trackingTags;    //!< Tag module -> define in which group(s) of detectors it will be used for tracking

  ReadonlyProperty<int         , AutoDefault> numSensors;      //!< Number of sensors built within detector module
  ReadonlyProperty<double      , Default    > stereoRotation;  //!< Stereo rotation of Z versus R-Phi measurement to calculate resolution (used if non-zero)
  ReadonlyProperty<double      , Computable > resolutionLocalX;//!< TODO: Resolution in RPhi direction (calculated or defined) -> rename to RPhi
  ReadonlyProperty<double      , Computable > resolutionLocalY;//!< TODO: Resolution in Z direction (calculated or defined) -> rename to Z
  ReadonlyProperty<SensorLayout, Default    > sensorLayout;
  ReadonlyProperty<ZCorrelation, NoDefault  > zCorrelation;
  ReadonlyProperty<ReadoutMode , Default    > readoutMode;     //!< Binary or analog (cluster) readout
  ReadonlyProperty<ReadoutType , Default    > readoutType;     //!< Readout strips/pixels/pT (stubs
  ReadonlyProperty<std::string , Default    > moduleType;
  Property<        short       , Default    > plotColor;       //!< Color used to draw module in the final graphical plots

  Property<int16_t, AutoDefault> side;
  
  ReadonlyProperty<int   , Default    > triggerWindow;
  ReadonlyProperty<int   , AutoDefault> numSparsifiedHeaderBits,  numSparsifiedPayloadBits;
  ReadonlyProperty<int   , AutoDefault> numTriggerDataHeaderBits, numTriggerDataPayloadBits;
  Property<        double, AutoDefault> sensorPowerConsumption;  // CUIDADO provide also power per strip (see original module and moduleType methods)
  ReadonlyProperty<double, AutoDefault> powerModuleOptical;
  ReadonlyProperty<double, AutoDefault> powerModuleChip;
  ReadonlyProperty<double, AutoDefault> powerStripOptical;
  ReadonlyProperty<double, AutoDefault> powerStripChip;
  Property<        double, AutoDefault> irradiationPower;
  ReadonlyProperty<double, Default    > triggerErrorX , triggerErrorY;
  ReadonlyProperty<bool  , Default    > reduceCombinatorialBackground;

  Property<double, Default> serviceHybridWidth;
  Property<double, Default> frontEndHybridWidth;
  Property<double, Default> hybridThickness;
  Property<double, Default> supportPlateThickness;

  int16_t cntId() const              { return m_cntId; }
  const std::string& cntName() const { return m_cntName; }
  void cntNameId(const std::string& name, int id) { m_cntName = name; m_cntId = id; }

  virtual PosRef   posRef()     const = 0;
  virtual TableRef tableRef()   const = 0;
  virtual UniRef   uniRef()     const = 0;
  virtual int16_t  moduleRing() const { return -1; }

  // ??? method
  //double trackCross(const XYZVector& PL, const XYZVector& PU) { return m_moduleGeom->trackCross(PL, PU); }

 protected:

  //!< Delete geometrical objects in a module
  void clearSensorPolys() { for (auto& s : m_sensors) s.clearPolys(); }

  MaterialObject         m_materialObject;
  Sensors                m_sensors;             //!< Module sensors
  ModuleCap*             m_moduleCap = nullptr; //!< Module materials assigned to its active part

  XYZVector m_rAxis;
  int       m_numHits   = 0; //!<  Number of hits assigned to a module (module hit Nxtimes)

  std::string m_cntName;
  int16_t     m_cntId;
  mutable double cachedZError_ = -1.;
  mutable std::pair<double,double> cachedMinMaxEtaWithError_;

 private:

  PropertyNode<int> sensorNode;

  double stripOccupancyPerEventBarrel() const;
  double stripOccupancyPerEventEndcap() const;
}; // Class


/*
 * @class BarrelModule
 * @details Derived class from DetectorModule class implementing barrel type detector module
 */
class BarrelModule : public DetectorModule, public Clonable<BarrelModule> {

 public:

  //! Constructor - parse geometry config file using boost property tree & read-in module parameters, specify unique id
  BarrelModule(int id, GeometricModule* moduleGeom, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty);

  //! Build barrel module with its sensors
  void build() override ;

  //! Setup: link lambda functions to various BarrelModule (and DetectorModule - because overriden) related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup() override;

  //! GeometryVisitor pattern -> barrel module visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> barrel module visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  //! Return module type
  virtual ModuleSubdetector subdet() const { return BARREL; }

  Property<int16_t, AutoDefault> layer;
  Property<int16_t, AutoDefault> rod;
  int16_t ring()       const { return (int16_t)myid(); }
  int16_t moduleRing() const { return ring(); }

  PosRef   posRef()   const { return (PosRef){ cntId(), (side() > 0 ? ring() : -ring()), layer(), rod() }; }
  TableRef tableRef() const { return (TableRef){ cntName(), layer(), ring() }; }
  UniRef   uniRef()   const { return  UniRef{ cntName(), layer(), ring(), rod(), side() }; }
}; // Class

/*
 * @class EndcapModule
 * @details Derived class from DetectorModule class implementing endcap type detector module
 */
class EndcapModule : public DetectorModule, public Clonable<EndcapModule> {

 public:

  //! Constructor - parse geometry config file using boost property tree & read-in module parameters, specify unique id
  EndcapModule(int id, GeometricModule* moduleGeom, const PropertyTree& treeProperty);

  //! Setup: link lambda functions to various EndcapModule (and DetectorModule - because overriden) related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup() override;

  //! Build endcap module with its sensors
  void build();

  //! GeometryVisitor pattern -> endcap module visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> endcap module visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  //! Return module type
  virtual ModuleSubdetector subdet() const { return ENDCAP; }

  Property<int16_t, AutoDefault> disk;
  Property<int16_t, AutoDefault> ring;
  int16_t moduleRing() const { return ring(); };
  int16_t blade()      const { return (int16_t)myid(); } // CUIDADO Think of a better name!
  int16_t side()       const { return (int16_t)signum(center().Z()); }

  PosRef posRef()     const { return (PosRef){ cntId(), (side() > 0 ? disk() : -disk()), ring(), blade() }; }
  TableRef tableRef() const { return (TableRef){ cntName(), disk(), ring() }; }
  UniRef uniRef()     const { return UniRef{ cntName(), disk(), ring(), blade(), side() }; }
}; // Class

#endif /* INCLUDE_DETECTOR_MODULE_H_ */

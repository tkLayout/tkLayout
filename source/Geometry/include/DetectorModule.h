#ifndef INCLUDE_DETECTOR_MODULE_H_
#define INCLUDE_DETECTOR_MODULE_H_

#include <boost/ptr_container/ptr_vector.hpp>

#include "Sensor.h"
#include "GeometricModule.h"
#include "GeometryFactory.h"
#include "CoordinateOperations.h"
#include "Property.h"
#include "Visitable.h"
#include "Visitor.h"
#include "MaterialObject.h"
#include "MessageLogger.h"

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
 * @class DetectorModule
 * @brief Detector module base class - base to Barrel & Endcap classes
 * @details Implementors of the DetectorModuleInterface must take care of rotating the module based on which part of the subdetector it will be used in (Barrel, EC)
 */
class DetectorModule : public Decorator<GeometricModule>, public PropertyObject, public Buildable, public Placeable, public Identifiable<int>, public Visitable
{

public:

  //! Constructor - specify unique id, geometry module defining shape & parse geometry config file using boost property tree & read-in module parameters
  DetectorModule(int id, Decorated* decorated, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty);

  //! Constructor - specify unique id, geometry module defining shape & parse geometry config file using boost property tree & read-in module parameters
  DetectorModule(int id, Decorated* decorated, const PropertyTree& treeProperty);

  //! Destructor
  ~DetectorModule();

  //! Build detector module with its sensors
  virtual void build();

  //! Setup: link lambda functions to various DetectorModule related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  virtual void setup();

  //! Set materials, so called module-cap to a module
  void setModuleCap(ModuleCap* newCap);

  //! Get materials, i.e. module-cap to be modified
  ModuleCap* getModuleCap() { return m_myModuleCap;}

  //! Get materials, i.e. module-cap to be read-only
  const ModuleCap& getModuleCap() const { return *m_myModuleCap; }

  //! Geometric module interface -> return reference to a geometrical representation of a module
  const Polygon3d<4>& basePoly() const { return decorated().basePoly(); }

  //! Could track at this direction hit the module?
  bool couldHit(const XYZVector& direction, double zError) const;

  //! Calculate min & max module eta coverage taking into account error on beam Z position
  std::pair<double, double> minMaxEtaWithError(double zError) const;

  //! Module R-Phi-resolution calculated as for a barrel-type module -> transform it to true orientation (rotation by theta angle, skew, tilt)
  double resolutionEquivalentRPhi(double hitRho, double trackR) const;

  //! Module Z-resolution calculated as for a barrel-type module -> transform it to true orientation (rotation by theta angle, skew, tilt)
  double resolutionEquivalentZ   (double hitRho, double trackR, double trackCotgTheta) const;

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
  ModuleShape shape() const { return decorated().shape(); }

  // Get geometry properties
  const XYZVector& center() const { return decorated().center(); }
  const XYZVector& normal() const { return decorated().normal(); }
  double area()             const { const GeometricModule& module = decorated(); return module.area(); }
  double dsDistance()       const { return decorated().dsDistance(); }

  double thickness()        const { return dsDistance() + sensorThickness(); }
  double length()           const { return decorated().length(); }
  double maxWidth()         const { return decorated().maxWidth(); }
  double minWidth()         const { return decorated().minWidth(); }
  double meanWidth()        const { return decorated().meanWidth(); }
  double physicalLength()   const { return decorated().physicalLength(); }
  double tiltAngle()        const { return m_tiltAngle; }
  double skewAngle()        const { return m_skewAngle; }

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
  void translate(const XYZVector& vector) { decorated().translate(vector); clearSensorPolys(); }
  void mirror(   const XYZVector& vector) { decorated().mirror(vector); clearSensorPolys(); }
  void translateZ(double z)               { decorated().translate(XYZVector(0, 0, z)); clearSensorPolys(); }
  void translateR(double radius)          { XYZVector v = m_rAxis.Unit()*radius; decorated().translate(v); clearSensorPolys();}
  void mirrorZ();
  void rotateX(double angle)              { decorated().rotateX(angle); clearSensorPolys(); }
  void rotateY(double angle)              { decorated().rotateY(angle); clearSensorPolys(); }
  void rotateZ(double angle)              { decorated().rotateZ(angle); clearSensorPolys(); m_rAxis = RotationZ(angle)(m_rAxis); }
  void tilt(double angle)                 { rotateX(-angle); m_tiltAngle += angle; } // Atention: Tilt can only be called BEFORE translating/rotating the module, or they won't work as expected!!
  void skew(double angle)                 { rotateY(-angle); m_skewAngle += angle; } // Atention: Skew can only be called BEFORE translating/rotating the module, or they won't work as expected!!
  void dsDistance(double d)               { decorated().dsDistance(d); }

  void resetHits() { m_numHits = 0; }

  // Properties
  ReadonlyProperty<double    , UncachedComputable> minZ;       //!< Minimum module Z (use uncached -> calculate its value anytime when needed -> need to be uncached recursively, i.e. also for Sensors)
  ReadonlyProperty<double    , UncachedComputable> maxZ;       //!< Maximum module Z
  ReadonlyProperty<double    , UncachedComputable> minR;       //!< Minimum module radius
  ReadonlyProperty<double    , UncachedComputable> maxR;       //!< Maximum module radius
  ReadonlyProperty<double    , UncachedComputable> planarMinZ;
  ReadonlyProperty<double    , UncachedComputable> planarMaxZ;
  ReadonlyProperty<double    , UncachedComputable> planarMinR;
  ReadonlyProperty<double    , UncachedComputable> planarMaxR;
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
  Property<        int8_t      , Default    > plotColor;       //!< Color used to draw module in the final graphical plots

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

  int16_t cntId() const { return cntId_; }
  const std::string& cntName() const { return cntName_; }
  void cntNameId(const std::string& name, int id) { cntName_ = name; cntId_ = id; }

  virtual PosRef   posRef()   const = 0;
  virtual TableRef tableRef() const = 0;
  virtual UniRef   uniRef()   const = 0;
  virtual int16_t  moduleRing() const { return -1; }

  double trackCross(const XYZVector& PL, const XYZVector& PU) { return decorated().trackCross(PL, PU); }
  std::pair<XYZVector, HitType> checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir) const;

 protected:

  //!< Delete geometrical objects in a module
  void clearSensorPolys() { for (auto& s : m_sensors) s.clearPolys(); }

  MaterialObject m_materialObject;
  Sensors        m_sensors;           //!< Detector sensors

  ModuleCap* m_myModuleCap = nullptr; //!< Module materials assigned to its active part

  XYZVector m_rAxis;
  double    m_tiltAngle = 0.; //!< Module tilt, i.e. rotation in RZ plane
  double    m_skewAngle = 0.; //!< Module skew, i.e. rotation in XY plane
  int       m_numHits   = 0;

  std::string cntName_;
  int16_t     cntId_;
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
  BarrelModule(int id, Decorated* decorated, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty);

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
  UniRef   uniRef()   const { return UniRef{ cntName(), layer(), ring(), rod(), side() }; }
}; // Class

/*
 * @class EndcapModule
 * @details Derived class from DetectorModule class implementing endcap type detector module
 */
class EndcapModule : public DetectorModule, public Clonable<EndcapModule> {

 public:

  //! Constructor - parse geometry config file using boost property tree & read-in module parameters, specify unique id
  EndcapModule(int id, Decorated* decorated, const PropertyTree& treeProperty);

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

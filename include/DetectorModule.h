#ifndef DETECTOR_MODULE_H
#define DETECTOR_MODULE_H

#include "Sensor.h"
#include "ModuleBase.h"
#include "GeometricModule.h"


//
// ======================================================= DETECTOR MODULES ===============================================================
//


enum ModuleSubdetector { BARREL = 1, ENDCAP = 2 };
enum SensorLayout { NOSENSORS, MONO, PT, STEREO };
enum InefficiencyType { STRIPWISE, EDGEWISE }; 
enum ReadoutType { READOUT_STRIP, READOUT_PIXEL, READOUT_PT };
enum ReadoutMode { BINARY, CLUSTER };



struct PosRef { int cnt, z, rho, phi; };
struct TableRef { string table; int row, col; };
struct UniRef { string cnt; int layer, ring, phi, side; };


class DetectorModule : public Decorator<GeometricModule>, public ModuleBase {// implementors of the DetectorModuleInterface must take care of rotating the module based on which part of the subdetector it will be used in (Barrel, EC)
  PropertyNode<int> sensorNode;
  typedef std::vector<Sensor> Sensors;
protected:
  Sensors sensors_;
  std::string cntName_;
  int16_t cntId_;
  mutable double cachedZError_ = -1.;
  mutable std::pair<double,double> cachedMinMaxEtaWithError_;
  XYZVector rAxis_;
  double tiltAngle_ = 0., skewAngle_ = 0.;
public:
  Property<int16_t, AutoDefault> side;

  Property<double, Computable> minPhi, maxPhi;

  ReadonlyProperty<std::string, Default> moduleType; 
  ReadonlyProperty<int, AutoDefault>     numSensors;
  ReadonlyProperty<SensorLayout, Default> sensorLayout;
  ReadonlyProperty<InefficiencyType, Default> inefficiencyType;
  ReadonlyProperty<ReadoutMode, Default> readoutMode;
  ReadonlyProperty<ReadoutType, Default> readoutType;

  ReadonlyProperty<int, Default> triggerWindow;

  ReadonlyProperty<int, AutoDefault> numSparsifiedHeaderBits,  numSparsifiedPayloadBits;
  ReadonlyProperty<int, AutoDefault> numTriggerDataHeaderBits, numTriggerDataPayloadBits;

  Property<double, AutoDefault> sensorPowerConsumption;  // CUIDADO provide also power per strip (see original module and moduleType methods)
  ReadonlyProperty<double, AutoDefault> powerModuleOptical;
  ReadonlyProperty<double, AutoDefault> powerModuleChip;
  ReadonlyProperty<double, AutoDefault> powerStripOptical;
  ReadonlyProperty<double, AutoDefault> powerStripChip;
  Property<double, AutoDefault> irradiationPower;

  ReadonlyProperty<double, Computable> resolutionLocalX, resolutionLocalY;
  ReadonlyProperty<double, Default>    triggerErrorX , triggerErrorY;

  ReadonlyProperty<double, AutoDefault> stereoRotation;
  
  ReadonlyProperty<bool, Default> reduceCombinatorialBackground;

  PropertyVector<string, ','> trackingTags;

  Property<int8_t, Default> plotColor;

  int16_t cntId() const { return cntId_; }
  const std::string& cntName() const { return cntName_; }
  void cntNameId(const std::string& name, int id) { cntName_ = name; cntId_ = id; }

  DetectorModule(Decorated* decorated) : 
      Decorator<GeometricModule>(decorated),
      sensorNode               ("Sensor"                   , parsedOnly()),
      moduleType               ("moduleType"               , parsedOnly() , string("notype")),
      numSensors               ("numSensors"               , parsedOnly()),
      sensorLayout             ("sensorLayout"             , parsedOnly() , NOSENSORS),
      readoutType              ("readoutType"              , parsedOnly() , READOUT_STRIP), 
      readoutMode              ("readoutMode"              , parsedOnly() , BINARY),
      inefficiencyType         ("inefficiencyType"         , parsedOnly() , EDGEWISE),
      numSparsifiedHeaderBits  ("numSparsifiedHeaderBits"  , parsedOnly()),
      numSparsifiedPayloadBits ("numSparsifiedPayloadBits" , parsedOnly()),
      numTriggerDataHeaderBits ("numTriggerDataHeaderBits" , parsedOnly()),
      numTriggerDataPayloadBits("numTriggerDataPayloadBits", parsedOnly()),
      triggerWindow            ("triggerWindow"            , parsedOnly() , 1),
      powerModuleOptical       ("powerModuleOptical"       , parsedOnly()),
      powerModuleChip          ("powerModuleChip"          , parsedOnly()),
      powerStripOptical        ("powerStripOptical"        , parsedOnly()),
      powerStripChip           ("powerStripChip"           , parsedOnly()),
      triggerErrorX            ("triggerErrorX"            , parsedOnly() , 1.),
      triggerErrorY            ("triggerErrorY"            , parsedOnly() , 1.),
      reduceCombinatorialBackground("reduceCombinatorialBackground", parsedOnly(), false),
      trackingTags             ("trackingTags"             , parsedOnly()),
      resolutionLocalX         ("resolutionLocalX"         , parsedOnly()),
      resolutionLocalY         ("resolutionLocalY"         , parsedOnly()),
      plotColor                ("plotColor"                , parsedOnly(), 0)
  {}

  virtual void setup();

  virtual DetectorModule* clone() = 0;

  virtual void build();
// Geometric module interface
  const Polygon3d<4>& basePoly() const { return decorated().basePoly(); }

  const XYZVector& center() const { return decorated().center(); }
  const XYZVector& normal() const { return decorated().normal(); }
  double area() const { return decorated().area(); }
  double dsDistance() const { return decorated().dsDistance(); }
  void dsDistance(double d) { decorated().dsDistance(d); }
  double thickness() const { return dsDistance() + sensorThickness(); }
  double length()    const { return decorated().length(); }
  double maxWidth()  const { return decorated().maxWidth(); }
  double minWidth()  const { return decorated().minWidth(); }
  double meanWidth() const { return decorated().meanWidth(); }

  double tiltAngle() const { return tiltAngle_; }
  double skewAngle() const { return skewAngle_; }

  double resolutionEquivalentZ   (double hitRho, double trackR, double trackCotgTheta) const;
  double resolutionEquivalentRPhi(double hitRho, double trackR) const;

  void translate(const XYZVector& vector) { decorated().translate(vector); for (auto& s : sensors_) s.poly().translate(vector); }
  void mirror(const XYZVector& vector) { decorated().mirror(vector); for (auto& s : sensors_) s.poly().mirror(vector); }
  void translateZ(double z) { decorated().translate(XYZVector(0, 0, z)); for (auto& s : sensors_) s.poly().translate(XYZVector(0, 0, z)); }
  void translateR(double radius) { 
    XYZVector v = rAxis_.Unit()*radius;
    decorated().translate(v); 
    for (auto& s : sensors_) s.poly().translate(v);
  }
  void mirrorZ() { 
    side(-side());
    for (auto& s : sensors_) s.poly().mirror(XYZVector(1. ,1. , -1.));
    decorated().mirror(XYZVector(1., 1., -1.));
  }

  void rotateX(double angle) { decorated().rotateX(angle); for (auto& s : sensors_) s.poly().rotateX(angle); }
  void rotateY(double angle) { decorated().rotateY(angle); for (auto& s : sensors_) s.poly().rotateY(angle); }
  void rotateZ(double angle) { decorated().rotateZ(angle); for (auto& s : sensors_) s.poly().rotateZ(angle); rAxis_ = RotationZ(angle)(rAxis_); }
  void tilt(double angle) { rotateX(-angle); tiltAngle_ += angle; } // CUIDADO!!! tilt and skew can only be called BEFORE translating/rotating the module, or they won't work as expected!!
  void skew(double angle) { rotateY(-angle); skewAngle_ += angle; }

  ModuleShape shape() const { return decorated().shape(); }
////////

  virtual double maxZ() const = 0;
  virtual double minZ() const = 0;
  virtual double maxR() const = 0;
  virtual double minR() const = 0;

  double phiAperture() const { return maxPhi() - minPhi(); }

  double maxEta() const { return MAX(basePoly().getVertex(0).Eta(), basePoly().getVertex(2).Eta()); }
  double minEta() const { return MIN(basePoly().getVertex(0).Eta(), basePoly().getVertex(2).Eta()); }
  double etaAperture() const { return maxEta() - minEta(); }
  double maxEtaWithError(double zError) const { return minMaxEtaWithError(zError).second; }
  double minEtaWithError(double zError) const { return minMaxEtaWithError(zError).first; }
  std::pair<double, double> minMaxEtaWithError(double zError) const;

  double maxTheta() const { return MAX(basePoly().getVertex(0).Theta(), basePoly().getVertex(2).Theta()); }
  double minTheta() const { return MIN(basePoly().getVertex(0).Theta(), basePoly().getVertex(2).Theta()); }
  double thetaAperture() const { return maxTheta() - minTheta(); }

  const Sensors& sensors() const { return sensors_; }
  const Sensor& innerSensor() const { return sensors_.front(); }
  const Sensor& outerSensor() const { return sensors_.back(); }
  int maxSegments() const { int segm = 0; for (const auto& s : sensors()) { segm = MAX(segm, s.numSegments()); } return segm; }
  int minSegments() const { int segm = 999999; for (const auto& s : sensors()) { segm = MIN(segm, s.numSegments()); } return segm; }
  int totalSegments() const { int cnt = 0; for (const auto& s : sensors()) { cnt += s.numSegments(); } return cnt; }
  int maxChannels() const { int max = 0; for (const auto& s : sensors()) { max = MAX(max, s.numChannels()); } return max; } 
  int minChannels() const { int min = 999999; for (const auto& s : sensors()) { min = MIN(min, s.numChannels()); } return min; } 
  int totalChannels() const { int cnt = 0; for (const auto& s : sensors()) { cnt += s.numChannels(); } return cnt; } 

  double totalPowerModule() const { return powerModuleOptical() + powerModuleChip(); }
  double totalPowerStrip() const { return powerStripOptical() + powerStripChip(); }
  double totalPower() const { return totalPowerModule() + totalPowerStrip()*outerSensor().numChannels(); }

  
  int numStripsAcross() const { return sensors().front().numStripsAcross(); }
  double sensorThickness() const { return sensors().front().sensorThickness(); } // CUIDADO this has to be fixed (called in Extractor.cc), sensor thickness can be different for different sensors

  virtual double stripOccupancyPerEvent() const = 0;
  virtual double hitOccupancyPerEvent() const { return stripOccupancyPerEvent()/2.; }
  virtual double geometricEfficiency() const = 0;
  virtual double effectiveDsDistance() const = 0;

  virtual ModuleSubdetector subdet() const = 0;

  virtual PosRef posRef() const = 0;
  virtual TableRef tableRef() const = 0;
  virtual UniRef uniRef() const = 0;

  bool couldHit(const XYZVector& direction, double zError) const;
  double trackCross(const XYZVector& PL, const XYZVector& PU) { return decorated().trackCross(PL, PU); }
  int numHits() const { return decorated().numHits(); }
  void resetHits() { decorated().resetHits(); }
};






class BarrelModule : public DetectorModule {
public:
  Property<int16_t, AutoDefault> layer;
  int16_t ring() const { return (int16_t)myid(); }
  Property<int16_t, AutoDefault> rod;

  BarrelModule(Decorated* decorated) : DetectorModule(decorated) {}
  BarrelModule* clone() override { return new BarrelModule(*this); }
  void accept(GeometryVisitor& v) { 
    v.visit(*this); 
    v.visit(*(DetectorModule*)this);
    decorated().accept(v); 
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this); 
    v.visit(*(const DetectorModule*)this);
    decorated().accept(v); 
  }

  void setup() override {
    DetectorModule::setup();
    minPhi.setup([&](){return MIN(basePoly().getVertex(0).Phi(), basePoly().getVertex(2).Phi());});
    maxPhi.setup([&](){return MAX(basePoly().getVertex(0).Phi(), basePoly().getVertex(2).Phi());});
  }

  void build();


  double maxZ() const { return MAX(basePoly().getVertex(0).Z(), basePoly().getVertex(2).Z()); } 
  double minZ() const { return MIN(basePoly().getVertex(0).Z(), basePoly().getVertex(2).Z()); } 
  double maxR() const { return MAX(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }
  double minR() const { return center().Rho(); }//MIN(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }

  void mirrorZ() { 
    side(-side());
    for (auto& s : sensors_) s.poly().mirror(XYZVector(1. ,1. , -1.));
    decorated().mirror(XYZVector(1., 1., -1.));
  }

  double stripOccupancyPerEvent() const;
  double geometricEfficiency() const;
  double effectiveDsDistance() const { return dsDistance(); }

  virtual ModuleSubdetector subdet() const { return BARREL; }

  PosRef posRef() const { return (PosRef){ cntId(), (side() > 0 ? ring() : -ring()), layer(), rod() }; }
  TableRef tableRef() const { return (TableRef){ cntName(), layer(), ring() }; }
  UniRef uniRef() const { return UniRef{ cntName(), layer(), ring(), rod(), side() }; }
};



class EndcapModule : public DetectorModule {
public:
  Property<int16_t, AutoDefault> cntId;
  Property<int16_t, AutoDefault> disk;
  Property<int16_t, AutoDefault> ring;
  int16_t blade() const { return (int16_t)myid(); } // CUIDADO Think of a better name!
  int16_t side() const { return (int16_t)signum(center().Z()); }

  EndcapModule(Decorated* decorated) : DetectorModule(decorated) {} 

  void setup() override {
    DetectorModule::setup();
    minPhi.setup([&](){return std::min_element(basePoly().begin(), basePoly().end(), [](const XYZVector& v1, const XYZVector& v2) { return v1.Phi() < v2.Phi(); })->Phi();});
    maxPhi.setup([&](){return std::max_element(basePoly().begin(), basePoly().end(), [](const XYZVector& v1, const XYZVector& v2) { return v1.Phi() < v2.Phi(); })->Phi();});
  }

  EndcapModule* clone() override { EndcapModule* m = new EndcapModule(*this); return m; }
  void build();

  void accept(GeometryVisitor& v) {
    v.visit(*this); 
    v.visit(*(DetectorModule*)this);
    decorated().accept(v); 
  }
  void accept(ConstGeometryVisitor& v) const { 
    v.visit(*this); 
    v.visit(*(const DetectorModule*)this);
    decorated().accept(v); 
  }

  double minZ() const { return center().Z(); }
  double maxZ() const { return center().Z(); }
  double maxR() const { return MAX(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }
  double minR() const { XYZVector side[2];
                        std::partial_sort_copy(basePoly().begin(), basePoly().end(), std::begin(side), std::end(side), [](const XYZVector& v1, const XYZVector& v2) { return v1.Rho() < v2.Rho(); });
                        return ((side[0]+side[1])/2).Rho(); }


  double stripOccupancyPerEvent() const;
  double geometricEfficiency() const;
  double effectiveDsDistance() const { return dsDistance() * center().Rho() / center().Z(); } 

  virtual ModuleSubdetector subdet() const { return ENDCAP; }

  PosRef posRef() const { return (PosRef){ cntId(), (side() > 0 ? disk() : -disk()), ring(), blade() }; }
  TableRef tableRef() const { return (TableRef){ cntName(), disk(), ring() }; }
  UniRef uniRef() const { return UniRef{ cntName(), disk(), ring(), blade(), side() }; }
};



// ===================================================================================================================================
//
#endif

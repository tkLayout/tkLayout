#ifndef DETECTOR_MODULE_H
#define DETECTOR_MODULE_H

#include <limits.h>

#include <boost/ptr_container/ptr_vector.hpp>

#include "Sensor.h"
#include "ModuleBase.h"
#include "GeometricModule.h"
#include "CoordinateOperations.h"
#include "Visitable.h"
#include "MaterialObject.h"

using material::MaterialObject;

//
// ======================================================= DETECTOR MODULES ===============================================================
//


enum ModuleSubdetector { BARREL = 1, ENDCAP = 2 };
enum SensorLayout { NOSENSORS, MONO, PT, STEREO };
enum ZCorrelation { SAMESEGMENT, MULTISEGMENT }; 
enum ReadoutType { READOUT_STRIP, READOUT_PIXEL, READOUT_PT };
enum ReadoutMode { BINARY, CLUSTER };
enum HitType { NONE, INNER, OUTER, BOTH = 3, STUB = 7 };



struct PosRef { int cnt, z, rho, phi; };
struct TableRef { string table; int row, col; };
struct UniRef { string cnt; int layer, ring, phi, side; };

namespace insur {
  class ModuleCap;
}
using insur::ModuleCap;
using material::ElementsVector;

class DetectorModule : public Decorator<GeometricModule>, public ModuleBase {// implementors of the DetectorModuleInterface must take care of rotating the module based on which part of the subdetector it will be used in (Barrel, EC)
  PropertyNode<int> sensorNode;

  typedef PtrVector<Sensor> Sensors;
  double stripOccupancyPerEventBarrel() const;
  double stripOccupancyPerEventEndcap() const;
protected:
  MaterialObject materialObject_;
  Sensors sensors_;
  std::string cntName_;
  int16_t cntId_;
  mutable double cachedZError_ = -1.;
  mutable std::pair<double,double> cachedMinMaxEtaWithError_;
  XYZVector rAxis_;
  double tiltAngle_ = 0., skewAngle_ = 0.;

  int numHits_ = 0;

  void clearSensorPolys() { for (auto& s : sensors_) s.clearPolys(); }
  ModuleCap* myModuleCap_ = NULL;
public:
  void setModuleCap(ModuleCap* newCap) { myModuleCap_ = newCap ; }
  ModuleCap* getModuleCap() { return myModuleCap_ ; }

  //virtual bool hasAnyResolutionLocalXParam();
  //virtual bool hasAnyResolutionLocalYParam();

  Property<int16_t, AutoDefault> side;
  
  Property<double, Computable> minPhi, maxPhi;
  
  ReadonlyProperty<std::string, Default> moduleType; 
  ReadonlyProperty<int, AutoDefault>     numSensors;
  ReadonlyProperty<SensorLayout, Default> sensorLayout;
  ReadonlyProperty<ZCorrelation, NoDefault> zCorrelation;
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

  ReadonlyProperty<double, Default> stereoRotation;
  
  ReadonlyProperty<bool, Default> reduceCombinatorialBackground;

  PropertyVector<string, ','> trackingTags;

  Property<int8_t, Default> plotColor;

  Property<double, Default> serviceHybridWidth;
  Property<double, Default> frontEndHybridWidth;
  Property<double, Default> hybridThickness;
  Property<double, Default> supportPlateThickness;

  int16_t cntId() const { return cntId_; }
  const std::string& cntName() const { return cntName_; }
  void cntNameId(const std::string& name, int id) { cntName_ = name; cntId_ = id; }
  
 DetectorModule(Decorated* decorated) : 
    Decorator<GeometricModule>(decorated),
      materialObject_(MaterialObject::MODULE),
      sensorNode               ("Sensor"                   , parsedOnly()),
      moduleType               ("moduleType"               , parsedOnly() , string("notype")),
      numSensors               ("numSensors"               , parsedOnly()),
      sensorLayout             ("sensorLayout"             , parsedOnly() , NOSENSORS),
      readoutType              ("readoutType"              , parsedOnly() , READOUT_STRIP), 
      readoutMode              ("readoutMode"              , parsedOnly() , BINARY),
      zCorrelation             ("zCorrelation"             , parsedOnly()),
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
      stereoRotation           ("stereoRotation"           , parsedOnly() , 0.),
      reduceCombinatorialBackground("reduceCombinatorialBackground", parsedOnly(), false),
      trackingTags             ("trackingTags"             , parsedOnly()),
      resolutionLocalX         ("resolutionLocalX"         , parsedOnly()),
      resolutionLocalY         ("resolutionLocalY"         , parsedOnly()),
      plotColor                ("plotColor"                , parsedOnly(), 0),
      serviceHybridWidth       ("serviceHybridWidth"       , parsedOnly(), 5),
      frontEndHybridWidth      ("frontEndHybridWidth"      , parsedOnly(), 5),
      hybridThickness          ("hybridThickness"          , parsedOnly(), 1),
      supportPlateThickness    ("supportPlateThickness"    , parsedOnly(), 1)
	{ }

    virtual bool hasAnyResolutionLocalXParam() const = 0;
    virtual bool hasAnyResolutionLocalYParam() const = 0;
    virtual void setup();
    virtual void build();
    // Geometric module interface
    const Polygon3d<4>& basePoly() const { return decorated().basePoly(); }

  const XYZVector& center() const { return decorated().center(); }
  const XYZVector& normal() const { return decorated().normal(); }
  double area() const { 
    //return decorated().area();
    const GeometricModule& module = decorated();
    double area = module.area(); 
    return area;
  }
  double dsDistance() const { return decorated().dsDistance(); }
  void dsDistance(double d) { decorated().dsDistance(d); }
  double thickness() const { return dsDistance() + sensorThickness(); }
  double length()    const { return decorated().length(); }
  double maxWidth()  const { return decorated().maxWidth(); }
  double minWidth()  const { return decorated().minWidth(); }
  double meanWidth() const { return decorated().meanWidth(); }
  double physicalLength() const { return decorated().physicalLength(); }

  double tiltAngle() const { return tiltAngle_; }
  double skewAngle() const { return skewAngle_; }

  double alpha (double trackPhi) const {
    double deltaPhi = center().Phi() + skewAngle() - trackPhi;
    if (fabs(deltaPhi) > M_PI/2.) {
      if (deltaPhi < 0.) deltaPhi = deltaPhi + 2.*M_PI;
      else deltaPhi = deltaPhi - 2.*M_PI;
    }  
    double alpha = deltaPhi + M_PI / 2.;
    return alpha; 
  }
  double beta (double theta) const { return theta + tiltAngle(); }
  virtual double calculateParameterizedResolutionLocalX(double phi) const = 0;
  virtual double calculateParameterizedResolutionLocalY(double theta) const = 0;
  double resolutionEquivalentZ   (double hitRho, double trackR, double trackCotgTheta, double resolutionLocalX, double resolutionLocalY) const;
  double resolutionEquivalentRPhi(double hitRho, double trackR, double resolutionLocalX, double resolutionLocalY) const;

  void translate(const XYZVector& vector) { decorated().translate(vector); clearSensorPolys(); }
  void mirror(const XYZVector& vector) { decorated().mirror(vector); clearSensorPolys(); }
  void translateZ(double z) { decorated().translate(XYZVector(0, 0, z)); clearSensorPolys(); }
  void translateR(double radius) {
    XYZVector v = rAxis_.Unit()*radius;
    decorated().translate(v);
    clearSensorPolys();
  }
  void mirrorZ() { 
    side(-side());
    double zTranslation = -center().Z();
    double zRotation = -center().Phi();
    translateZ(zTranslation);
    rotateZ(zRotation);
    rotateY(M_PI);
    translateZ(zTranslation);
    rotateZ(-zRotation);
    //decorated().mirror(XYZVector(1., 1., -1.));
    clearSensorPolys();
  }

  void rotateX(double angle) { decorated().rotateX(angle); clearSensorPolys(); }
  void rotateY(double angle) { decorated().rotateY(angle); clearSensorPolys(); }
  void rotateZ(double angle) { decorated().rotateZ(angle); clearSensorPolys(); rAxis_ = RotationZ(angle)(rAxis_); }
  void tilt(double angle) { rotateX(-angle); tiltAngle_ += angle; } // CUIDADO!!! tilt and skew can only be called BEFORE translating/rotating the module, or they won't work as expected!!
  void skew(double angle) { rotateY(-angle); skewAngle_ += angle; }

  bool flipped() const { return decorated().flipped(); } 
  bool flipped(bool newFlip) { return decorated().flipped(newFlip); } 
  ModuleShape shape() const { return decorated().shape(); }
////////

  double maxZ() const { return maxget2(sensors_.begin(), sensors_.end(), &Sensor::maxZ); }
  double minZ() const { return minget2(sensors_.begin(), sensors_.end(), &Sensor::minZ); }
  double maxR() const { return maxget2(sensors_.begin(), sensors_.end(), &Sensor::maxR); }
  double minR() const { return minget2(sensors_.begin(), sensors_.end(), &Sensor::minR); }

  double planarMaxZ() const { return CoordinateOperations::computeMaxZ(basePoly()); }
  double planarMinZ() const { return CoordinateOperations::computeMinZ(basePoly()); }
  double planarMaxR() const { return CoordinateOperations::computeMaxR(basePoly()); }
  double planarMinR() const { return CoordinateOperations::computeMinR(basePoly()); }

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
  const MaterialObject& materialObject() const { return materialObject_; }
  const Sensor& innerSensor() const { return sensors_.front(); }
  const Sensor& outerSensor() const { return sensors_.back(); }
  ElementsVector& getLocalElements() const {return materialObject_.getLocalElements(); }
  int maxSegments() const { int segm = 0; for (const auto& s : sensors()) { segm = MAX(segm, s.numSegments()); } return segm; } // CUIDADO NEEDS OPTIMIZATION (i.e. caching or just MAX())
  int minSegments() const { int segm = std::numeric_limits<int>::max(); for (const auto& s : sensors()) { segm = MIN(segm, s.numSegments()); } return segm; }
  int totalSegments() const { int cnt = 0; for (const auto& s : sensors()) { cnt += s.numSegments(); } return cnt; }
  int maxChannels() const { int max = 0; for (const auto& s : sensors()) { max = MAX(max, s.numChannels()); } return max; } 
  int minChannels() const { int min = std::numeric_limits<int>::max(); for (const auto& s : sensors()) { min = MIN(min, s.numChannels()); } return min; } 
  int totalChannels() const { int cnt = 0; for (const auto& s : sensors()) { cnt += s.numChannels(); } return cnt; } 

  double totalPowerModule() const { return powerModuleOptical() + powerModuleChip(); }
  double totalPowerStrip() const { return powerStripOptical() + powerStripChip(); }
  double totalPower() const { return totalPowerModule() + totalPowerStrip()*outerSensor().numChannels(); }

  
  int numStripsAcross() const { return sensors().front().numStripsAcross(); } // CUIDADO this assumes both sensors have the same number of sensing elements in the transversal direction - typically it is like that
  double sensorThickness() const { return sensors().front().sensorThickness(); } // CUIDADO this has to be fixed (called in Extractor.cc), sensor thickness can be different for different sensors

  double stripOccupancyPerEvent() const;
  double hitOccupancyPerEvent() const { return stripOccupancyPerEvent()/2.; }
  double geometricEfficiency() const;
  double effectiveDsDistance() const;

  virtual ModuleSubdetector subdet() const = 0;

  virtual PosRef posRef() const = 0;
  virtual TableRef tableRef() const = 0;
  virtual UniRef uniRef() const = 0;
  virtual int16_t moduleRing() const { return -1; }

  bool couldHit(const XYZVector& direction, double zError) const;
  double trackCross(const XYZVector& PL, const XYZVector& PU) { return decorated().trackCross(PL, PU); }
  std::pair<XYZVector, HitType> checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir);
  int numHits() const { return numHits_; }
  void resetHits() { numHits_ = 0; }

};



class BarrelModule : public DetectorModule, public Clonable<BarrelModule> {
public:
  Property<int16_t, AutoDefault> layer;
  int16_t ring() const { return (int16_t)myid(); }
  int16_t moduleRing() const { return ring(); }
  Property<int16_t, AutoDefault> rod;
  ReadonlyProperty<double, NoDefault> cotalphaLimit;
  ReadonlyProperty<double, NoDefault> resolutionLocalXBarrelParam0Inf;
  ReadonlyProperty<double, NoDefault> resolutionLocalXBarrelParam1Inf;
  ReadonlyProperty<double, NoDefault> resolutionLocalXBarrelParam2Inf;
  ReadonlyProperty<double, NoDefault> resolutionLocalXBarrelParam0Sup;
  ReadonlyProperty<double, NoDefault> resolutionLocalXBarrelParam1Sup;
  ReadonlyProperty<double, NoDefault> resolutionLocalXBarrelParam2Sup;
  ReadonlyProperty<double, NoDefault> resolutionLocalYBarrelParam0;
  ReadonlyProperty<double, NoDefault> resolutionLocalYBarrelParam1;
  ReadonlyProperty<double, NoDefault> resolutionLocalYBarrelParam2;
  ReadonlyProperty<double, NoDefault> resolutionLocalYBarrelParam3;
  ReadonlyProperty<double, NoDefault> resolutionLocalYBarrelParam4;

 BarrelModule(Decorated* decorated) :
  DetectorModule(decorated),
    cotalphaLimit                           ("cotalphaLimit"                           , parsedOnly()),
    resolutionLocalXBarrelParam0Inf         ("resolutionLocalXBarrelParam0Inf"         , parsedOnly()),
    resolutionLocalXBarrelParam1Inf         ("resolutionLocalXBarrelParam1Inf"         , parsedOnly()),
    resolutionLocalXBarrelParam2Inf         ("resolutionLocalXBarrelParam2Inf"         , parsedOnly()),
    resolutionLocalXBarrelParam0Sup         ("resolutionLocalXBarrelParam0Sup"         , parsedOnly()),
    resolutionLocalXBarrelParam1Sup         ("resolutionLocalXBarrelParam1Sup"         , parsedOnly()),
    resolutionLocalXBarrelParam2Sup         ("resolutionLocalXBarrelParam2Sup"         , parsedOnly()),
    resolutionLocalYBarrelParam0            ("resolutionLocalYBarrelParam0"            , parsedOnly()),
    resolutionLocalYBarrelParam1            ("resolutionLocalYBarrelParam1"            , parsedOnly()),
    resolutionLocalYBarrelParam2            ("resolutionLocalYBarrelParam2"            , parsedOnly()),
    resolutionLocalYBarrelParam3            ("resolutionLocalYBarrelParam3"            , parsedOnly()),
    resolutionLocalYBarrelParam4            ("resolutionLocalYBarrelParam4"            , parsedOnly())
      { setup(); }

  bool hasAnyResolutionLocalXParam() const { return (resolutionLocalXBarrelParam0Inf.state() || resolutionLocalXBarrelParam1Inf.state() || resolutionLocalXBarrelParam2Inf.state() || resolutionLocalXBarrelParam0Sup.state() || resolutionLocalXBarrelParam1Sup.state() || resolutionLocalXBarrelParam2Sup.state()); }

  bool hasAnyResolutionLocalYParam() const { return (resolutionLocalYBarrelParam0.state() || resolutionLocalYBarrelParam1.state() || resolutionLocalYBarrelParam2.state() || resolutionLocalYBarrelParam3.state() || resolutionLocalYBarrelParam4.state()); }

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
    minPhi.setup([&](){ return MIN(basePoly().getVertex(0).Phi(), basePoly().getVertex(2).Phi()); });
    maxPhi.setup([&](){ return MAX(basePoly().getVertex(0).Phi(), basePoly().getVertex(2).Phi()); });
    resolutionLocalX.setup([this]() {
	// only set up this if no model parameter specified
	//std::cout <<  "hasAnyResolutionLocalYParam() = " <<  hasAnyResolutionLocalYParam() << std::endl;

	if (!hasAnyResolutionLocalXParam()) {
	  //std::cout << "resolutionLocalX and resolutionLocalXBarrel parameters are all unset. Use of default formulae." << std::endl;
	  double res = 0;
	  for (const Sensor& s : sensors()) res += pow(meanWidth() / s.numStripsAcross() / sqrt(12), 2);
	  return sqrt(res)/numSensors();
	}
	// if model parameters specified, return -1
	else return -1.0;
      });
    resolutionLocalY.setup([this]() {
	// only set up this if no model parameters not specified
	if (!hasAnyResolutionLocalYParam()) {
	  //std::cout << "resolutionLocalY and resolutionLocalYBarrel parameters are all unset. Use of default formulae." << std::endl;
	    if (stereoRotation() != 0.) return resolutionLocalX() / sin(stereoRotation());
	    else {
	      return length() / maxSegments() / sqrt(12); // NOTE: not combining measurements from both sensors. The two sensors are closer than the length of the longer sensing element, making the 2 measurements correlated. considering only the best measurement is then a reasonable approximation (since in case of a PS module the strip measurement increases the precision by only 0.2% and in case of a 2S the sensors are so close that they basically always measure the same thing)
	    }
	  }
	// if model parameters specified, return -1
	else return -1.0;
      });
  }
  
  void check() override;

  void build();


  //double maxZ() const { return MAX(basePoly().getVertex(0).Z(), basePoly().getVertex(2).Z()); } 
  //double minZ() const { return MIN(basePoly().getVertex(0).Z(), basePoly().getVertex(2).Z()); } 
  //double maxR() const { return MAX(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }
  //double minR() const { return center().Rho(); }//MIN(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }

  virtual ModuleSubdetector subdet() const { return BARREL; }

  double calculateParameterizedResolutionLocalX(double trackPhi) const { 
    double resolutionLocalXBarrelParam0, resolutionLocalXBarrelParam1, resolutionLocalXBarrelParam2;
    if ((1./tan(alpha(trackPhi))) < cotalphaLimit()) { resolutionLocalXBarrelParam0 = resolutionLocalXBarrelParam0Inf(); resolutionLocalXBarrelParam1 = resolutionLocalXBarrelParam1Inf(); resolutionLocalXBarrelParam2 = resolutionLocalXBarrelParam2Inf(); }
    else { resolutionLocalXBarrelParam0 = resolutionLocalXBarrelParam0Sup(); resolutionLocalXBarrelParam1 = resolutionLocalXBarrelParam1Sup(); resolutionLocalXBarrelParam2 = resolutionLocalXBarrelParam2Sup(); }
    return resolutionLocalXBarrelParam0 + resolutionLocalXBarrelParam1 * 1./tan(alpha(trackPhi)) + resolutionLocalXBarrelParam2 * pow(1./tan(alpha(trackPhi)), 2); 
}

  double calculateParameterizedResolutionLocalY(double theta) const { return resolutionLocalYBarrelParam0() + resolutionLocalYBarrelParam1() * exp(-resolutionLocalYBarrelParam2() * fabs(1./tan(beta(theta)))) * sin(resolutionLocalYBarrelParam3() * fabs(1./tan(beta(theta))) + resolutionLocalYBarrelParam4()); }

  PosRef posRef() const { return (PosRef){ cntId(), (side() > 0 ? ring() : -ring()), layer(), rod() }; }
  TableRef tableRef() const { return (TableRef){ cntName(), layer(), ring() }; }
  UniRef uniRef() const { return UniRef{ cntName(), layer(), ring(), rod(), side() }; }
};



class EndcapModule : public DetectorModule, public Clonable<EndcapModule> {
public:
  Property<int16_t, AutoDefault> disk;
  Property<int16_t, AutoDefault> ring;
  int16_t moduleRing() const { return ring(); };
  int16_t blade() const { return (int16_t)myid(); } // CUIDADO Think of a better name!
  int16_t side() const { return (int16_t)signum(center().Z()); }
  //bool hasAnyResolutionLocalXParam() override { return (resolutionLocalXEndcapParam0.state() || resolutionLocalXEndcapParam1.state()); }
  //bool hasAnyResolutionLocalYParam() override { return (resolutionLocalYEndcapParam0.state() || resolutionLocalYEndcapParam1.state()); }
  ReadonlyProperty<double, NoDefault> resolutionLocalXEndcapParam0;
  ReadonlyProperty<double, NoDefault> resolutionLocalXEndcapParam1;
  ReadonlyProperty<double, NoDefault> resolutionLocalYEndcapParam0;
  ReadonlyProperty<double, NoDefault> resolutionLocalYEndcapParam1;

 EndcapModule(Decorated* decorated) :
  DetectorModule(decorated),
    resolutionLocalXEndcapParam0            ("resolutionLocalXEndcapParam0"            , parsedOnly()),
    resolutionLocalXEndcapParam1            ("resolutionLocalXEndcapParam1"            , parsedOnly()),
    resolutionLocalYEndcapParam0            ("resolutionLocalYEndcapParam0"            , parsedOnly()),
    resolutionLocalYEndcapParam1            ("resolutionLocalYEndcapParam1"            , parsedOnly())
      { setup(); }

  bool hasAnyResolutionLocalXParam() const { return (resolutionLocalXEndcapParam0.state() || resolutionLocalXEndcapParam1.state()); }
  
  bool hasAnyResolutionLocalYParam() const { return (resolutionLocalYEndcapParam0.state() || resolutionLocalYEndcapParam1.state()); }

  void setup() override {
    DetectorModule::setup();
    minPhi.setup([&](){ return minget2(basePoly().begin(), basePoly().end(), &XYZVector::Phi); });
    maxPhi.setup([&](){ return maxget2(basePoly().begin(), basePoly().end(), &XYZVector::Phi); });
    resolutionLocalX.setup([this]() {
	// only set up this if no model parameter specified
	if (!hasAnyResolutionLocalXParam()) {
	  //std::cout << "resolutionLocalX and resolutionLocalXEndcap parameters are all unset. Use of default formulae." << std::endl;
	    double res = 0;
	    for (const Sensor& s : sensors()) res += pow(meanWidth() / s.numStripsAcross() / sqrt(12), 2);
	    return sqrt(res)/numSensors();
	  }
	// if model parameters specified, return -1
	else return -1.0;
      });
    resolutionLocalY.setup([this]() {
	// only set up this if no model parameters not specified
	if (!hasAnyResolutionLocalYParam()) {
	  //std::cout << "resolutionLocalY and resolutionLocalYEndcap parameters are all unset. Use of default formulae." << std::endl;
	    if (stereoRotation() != 0.) return resolutionLocalX() / sin(stereoRotation());
	    else {
	      return length() / maxSegments() / sqrt(12); // NOTE: not combining measurements from both sensors. The two sensors are closer than the length of the longer sensing element, making the 2 measurements correlated. considering only the best measurement is then a reasonable approximation (since in case of a PS module the strip measurement increases the precision by only 0.2% and in case of a 2S the sensors are so close that they basically always measure the same thing)
	    }
	  }
	// if model parameters specified, return -1
	else return -1.0;
      });
  }

  void check() override;

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

  //double minZ() const { return center().Z(); } // CUIDADO not accounting for sensor placement
  //double maxZ() const { return center().Z(); } // ditto here
  //double maxR() const { return MAX(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }
  //double minR() const { XYZVector side[2];
  //                      std::partial_sort_copy(basePoly().begin(), basePoly().end(), std::begin(side), std::end(side), [](const XYZVector& v1, const XYZVector& v2) { return v1.Rho() < v2.Rho(); });
  //                      return ((side[0]+side[1])/2).Rho(); }


  virtual ModuleSubdetector subdet() const { return ENDCAP; }

  double calculateParameterizedResolutionLocalX(double trackPhi) const { return resolutionLocalXEndcapParam0() + resolutionLocalXEndcapParam1() * 1./tan(alpha(trackPhi)); }

  double calculateParameterizedResolutionLocalY(double theta) const { return resolutionLocalYEndcapParam0() + resolutionLocalYEndcapParam1() * fabs(1./tan(beta(theta))); }

  PosRef posRef() const { return (PosRef){ cntId(), (side() > 0 ? disk() : -disk()), ring(), blade() }; }
  TableRef tableRef() const { return (TableRef){ cntName(), disk(), ring() }; }
  UniRef uniRef() const { return UniRef{ cntName(), disk(), ring(), blade(), side() }; }
};


// ===================================================================================================================================
//
#endif

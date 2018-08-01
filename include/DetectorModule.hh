#ifndef DETECTOR_MODULE_H
#define DETECTOR_MODULE_H

#include <limits.h>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/count.hpp>

#include "global_constants.hh"
#include "SimParms.hh"
#include "MessageLogger.hh"
#include "Sensor.hh"
#include "ModuleBase.hh"
#include "GeometricModule.hh"
#include "CoordinateOperations.hh"
#include "Visitable.hh"
#include "MaterialObject.hh"

namespace insur { class ModuleCap; }
using insur::ModuleCap;
using material::ElementsVector;

using namespace boost::accumulators;
using material::MaterialObject;

// OT CABLING MAP
namespace insur { class OuterBundle; }
using insur::OuterBundle;
namespace insur { class OuterDTC; }
using insur::OuterDTC;
// IT CABLING MAP
namespace insur { class PowerChain; }
using insur::PowerChain;
namespace insur { class HvLine; }
using insur::HvLine;
namespace insur { class GBT; }
using insur::GBT;
namespace insur { class InnerBundle; }
using insur::InnerBundle;
namespace insur { class InnerDTC; }
using insur::InnerDTC;


//
// ======================================================= DETECTOR MODULES ===============================================================
//

enum SensorLayout { NOSENSORS, MONO, PT, STEREO };
enum ZCorrelation { SAMESEGMENT, MULTISEGMENT }; 
enum ReadoutType { READOUT_STRIP, READOUT_PIXEL, READOUT_PT };
enum ReadoutMode { BINARY, CLUSTER };
enum HitType { NONE, INNER, OUTER, BOTH = 3, STUB = 7 };

struct PosRef { int subdetectorId, z, rho, phi; };
struct TableRef { string table; int row, col; };
struct UniRef { string subdetectorName; int layer, ring, phi, side; };


class DetectorModule : public Decorator<GeometricModule>, public ModuleBase, public DetIdentifiable {// implementors of the DetectorModuleInterface must take care of rotating the module based on which part of the subdetector it will be used in (Barrel, EC) 

  PropertyNode<int> sensorNode;
  typedef PtrVector<Sensor> Sensors;
  double stripOccupancyPerEventBarrel() const;
  double stripOccupancyPerEventEndcap() const;

public:
  void setModuleCap(ModuleCap* newCap) { myModuleCap_ = newCap ; }
  ModuleCap* getModuleCap() { return myModuleCap_ ; }
  const ModuleCap* getConstModuleCap() const { return myModuleCap_; }

  Property<int16_t, AutoDefault> side;
  Property<double, AutoDefault> skewAngle;
  
  Property<double, Computable> minPhi, maxPhi;
  
  ReadonlyProperty<std::string, Default> moduleType; 
  ReadonlyProperty<int, AutoDefault>     numSensors;
  ReadonlyProperty<SensorLayout, Default> sensorLayout;
  ReadonlyProperty<ZCorrelation, NoDefault> zCorrelation;
  ReadonlyProperty<ReadoutMode, Default> readoutMode;
  ReadonlyProperty<ReadoutType, Default> readoutType;
  ReadonlyProperty<double, Default> singleHitEfficiency;

  ReadonlyProperty<int, Default> triggerWindow;

  ReadonlyProperty<int, AutoDefault> numSparsifiedHeaderBits,  numSparsifiedPayloadBits;
  ReadonlyProperty<int, AutoDefault> numTriggerDataHeaderBits, numTriggerDataPayloadBits;

  ReadonlyProperty<double, NoDefault> operatingTemp;
  ReadonlyProperty<double, NoDefault> biasVoltage;
  ReadonlyProperty<double, AutoDefault> powerPerModule;
  Property<double, AutoDefault> sensorsIrradiationPowerMean;
  Property<double, AutoDefault> sensorsIrradiationPowerMax;
  Property<double, AutoDefault> sensorsIrradiationMean;
  Property<double, AutoDefault> sensorsIrradiationMax;

  ReadonlyProperty<double, Computable> nominalResolutionLocalX, nominalResolutionLocalY;
  // Local X resolution parameters
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam0;
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam1;
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam2;
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam3;
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam4;
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam5;
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam6;
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam7;
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam8;
  ReadonlyProperty<double, NoDefault> resolutionLocalXParam9;
  // Local Y resolution parameters
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam0;
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam1;
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam2;
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam3;
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam4;
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam5;
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam6;
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam7;
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam8;
  ReadonlyProperty<double, NoDefault> resolutionLocalYParam9;

  ReadonlyProperty<double, Default>    triggerErrorX , triggerErrorY;

  ReadonlyProperty<double, Default> stereoRotation;
  
  ReadonlyProperty<bool, Default> reduceCombinatorialBackground;

  PropertyVector<string, ','> trackingTags;

  Property<int, Default> plotColor;

  Property<double, Default> serviceHybridWidth;
  Property<double, Default> frontEndHybridWidth;
  Property<double, Default> hybridThickness;
  Property<double, Default> supportPlateThickness;
  Property<double, Default> chipThickness;

  Property<bool, Default> removeModule;

  int16_t subdetectorId() const { return subdetectorId_; }
  const std::string& subdetectorName() const { return subdetectorName_; }
  void subdetectorNameId(const std::string& name, const int id) { subdetectorName_ = name; subdetectorId_ = id; }
  
 DetectorModule(Decorated* decorated) : 
    Decorator<GeometricModule>(decorated),
      materialObject_(MaterialObject::MODULE),      sensorNode               ("Sensor"                   , parsedOnly()),
      skewAngle                ("skewAngle"                , parsedOnly()),
      moduleType               ("moduleType"               , parsedOnly() , string("notype")),
      numSensors               ("numSensors"               , parsedOnly()),
      sensorLayout             ("sensorLayout"             , parsedOnly() , NOSENSORS),
      readoutType              ("readoutType"              , parsedOnly() , READOUT_STRIP), 
      singleHitEfficiency      ("singleHitEfficiency"      , parsedOnly() , 1.),
      readoutMode              ("readoutMode"              , parsedOnly() , BINARY),
      zCorrelation             ("zCorrelation"             , parsedOnly()),
      numSparsifiedHeaderBits  ("numSparsifiedHeaderBits"  , parsedOnly()),
      numSparsifiedPayloadBits ("numSparsifiedPayloadBits" , parsedOnly()),
      numTriggerDataHeaderBits ("numTriggerDataHeaderBits" , parsedOnly()),
      numTriggerDataPayloadBits("numTriggerDataPayloadBits", parsedOnly()),
      triggerWindow            ("triggerWindow"            , parsedOnly() , 1),
      operatingTemp            ("operatingTemp"            , parsedAndChecked()),
      biasVoltage              ("biasVoltage"              , parsedAndChecked()),
      powerPerModule           ("powerPerModule"           , parsedOnly()),
      triggerErrorX            ("triggerErrorX"            , parsedOnly() , 1.),
      triggerErrorY            ("triggerErrorY"            , parsedOnly() , 1.),
      stereoRotation           ("stereoRotation"           , parsedOnly() , 0.),
      reduceCombinatorialBackground("reduceCombinatorialBackground", parsedOnly(), false),
      trackingTags             ("trackingTags"             , parsedOnly()),
      nominalResolutionLocalX  ("nominalResolutionLocalX"  , parsedOnly()),
      nominalResolutionLocalY  ("nominalResolutionLocalY"  , parsedOnly()),
      // Local X resolution parameters
      resolutionLocalXParam0            ("resolutionLocalXParam0"            , parsedOnly()),
      resolutionLocalXParam1            ("resolutionLocalXParam1"            , parsedOnly()),
      resolutionLocalXParam2            ("resolutionLocalXParam2"            , parsedOnly()),
      resolutionLocalXParam3            ("resolutionLocalXParam3"            , parsedOnly()),
      resolutionLocalXParam4            ("resolutionLocalXParam4"            , parsedOnly()),
      resolutionLocalXParam5            ("resolutionLocalXParam5"            , parsedOnly()),
      resolutionLocalXParam6            ("resolutionLocalXParam6"            , parsedOnly()),
      resolutionLocalXParam7            ("resolutionLocalXParam7"            , parsedOnly()),
      resolutionLocalXParam8            ("resolutionLocalXParam8"            , parsedOnly()),
      resolutionLocalXParam9            ("resolutionLocalXParam9"            , parsedOnly()),
      // Local Y resolution parameters
      resolutionLocalYParam0            ("resolutionLocalYParam0"            , parsedOnly()),
      resolutionLocalYParam1            ("resolutionLocalYParam1"            , parsedOnly()),
      resolutionLocalYParam2            ("resolutionLocalYParam2"            , parsedOnly()),
      resolutionLocalYParam3            ("resolutionLocalYParam3"            , parsedOnly()),
      resolutionLocalYParam4            ("resolutionLocalYParam4"            , parsedOnly()),
      resolutionLocalYParam5            ("resolutionLocalYParam5"            , parsedOnly()),
      resolutionLocalYParam6            ("resolutionLocalYParam6"            , parsedOnly()),
      resolutionLocalYParam7            ("resolutionLocalYParam7"            , parsedOnly()),
      resolutionLocalYParam8            ("resolutionLocalYParam8"            , parsedOnly()),
      resolutionLocalYParam9            ("resolutionLocalYParam9"            , parsedOnly()),
      plotColor                ("plotColor"                , parsedOnly(), 0),
      serviceHybridWidth       ("serviceHybridWidth"       , parsedOnly(), 0),
      frontEndHybridWidth      ("frontEndHybridWidth"      , parsedOnly(), 0),
      hybridThickness          ("hybridThickness"          , parsedOnly(), 0),
      supportPlateThickness    ("supportPlateThickness"    , parsedOnly(), 0),
      chipThickness            ("chipThickness"            , parsedOnly(), 0),
      removeModule             ("removeModule"             , parsedOnly(), false)
	{ }

  virtual void setup();
  void check() override;
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
  double totalSensorsVolume() const { // Calculate total volume occupied by sensors
    double volume = 0.;
    for (const auto& s : sensors()) volume += area() * s.sensorThickness();
    return volume;
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
  bool isTilted() const { return tiltAngle_ != 0.; }

  // SPATIAL RESOLUTION

  // RETURN GLOBAL SPATIAL RESOLUTION (IN CMS COORDINATES)
  double resolutionEquivalentZ   (double hitRho, double trackR, double trackCotgTheta, double resolutionLocalX, double resolutionLocalY) const;
  double resolutionEquivalentRPhi(double hitRho, double trackR, double resolutionLocalX, double resolutionLocalY) const;

  // RETURN LOCAL SPATIAL RESOLUTION
  // This resolution is either nominal, either from parametrization.
  const double resolutionLocalX(const double phi) const;
  const double resolutionLocalY(const double theta) const;

  // Used to compute local parametrized spatial resolution.
  const bool hasAnyResolutionLocalXParam() const;
  const bool hasAnyResolutionLocalYParam() const;
  const double alpha(const double trackPhi) const;
  const double beta(const double theta) const;
 
  // STATISTICS ON LOCAL SPATIAL RESOLUTION
  accumulator_set<double, features<tag::count, tag::mean, tag::variance, tag::sum, tag::moment<2>>> rollingParametrizedResolutionLocalX;
  accumulator_set<double, features<tag::count, tag::mean, tag::variance, tag::sum, tag::moment<2>>> rollingParametrizedResolutionLocalY;
  

  void translate(const XYZVector& vector) { decorated().translate(vector); clearSensorPolys(); }
  void mirror(const XYZVector& vector) { decorated().mirror(vector); clearSensorPolys(); }
  void translateZ(double z) { decorated().translate(XYZVector(0, 0, z)); clearSensorPolys(); }
  void translateR(double radius) {
    XYZVector v = rAxis_.Unit()*radius;
    decorated().translate(v);
    clearSensorPolys();
  }
  void rotateToNegativeZSide() {
    side(-side());
    rotateY(M_PI);  // Rotation around CMS_Y of angle Pi
    clearSensorPolys();
  }

  void rotateX(double angle) { decorated().rotateX(angle); clearSensorPolys(); }
  void rotateY(double angle) { decorated().rotateY(angle); clearSensorPolys(); }
  void rotateZ(double angle) { decorated().rotateZ(angle); clearSensorPolys(); rAxis_ = RotationZ(angle)(rAxis_); }
  void tilt(double angle) { rotateX(-angle); tiltAngle_ += angle; } // CUIDADO!!! tilt and skew can only be called BEFORE translating/rotating the module, or they won't work as expected!!
  // void skew(double angle) { rotateY(-angle); skewAngle_ += angle; } // This works for endcap modules only !!
  // Skew is now defined at construction time instead, before the module has had a chance to be translated/rotated!

  bool flipped() const { return decorated().flipped(); } 
  bool flipped(bool newFlip) {
    if (newFlip && numSensors() > 1) {
      sensors_.front().innerOuter(SensorPosition::UPPER);
      sensors_.back().innerOuter(SensorPosition::LOWER);
    }
    return decorated().flipped(newFlip);
  } 
  ModuleShape shape() const { return decorated().shape(); }
  ////////

  double maxZ() const { return maxget2(sensors_.begin(), sensors_.end(), &Sensor::maxZ); }
  double minZ() const { return minget2(sensors_.begin(), sensors_.end(), &Sensor::minZ); }
  double maxR() const { return maxget2(sensors_.begin(), sensors_.end(), &Sensor::maxR); }
  double minR() const { return minget2(sensors_.begin(), sensors_.end(), &Sensor::minR); }

  std::map<std::string, double> extremaWithHybrids() const;
  double minZwithHybrids() const { return extremaWithHybrids()["minZ"]; }
  double maxZwithHybrids() const { return extremaWithHybrids()["maxZ"]; }
  double minRwithHybrids() const { return extremaWithHybrids()["minR"]; }
  double maxRwithHybrids() const { return extremaWithHybrids()["maxR"]; }

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
  int maxSegments() const { int segm = 0; for (const auto& s : sensors()) { segm = MAX(segm, s.numSegmentsEstimate()); } return segm; } // CUIDADO NEEDS OPTIMIZATION (i.e. caching or just MAX())
  int minSegments() const { int segm = std::numeric_limits<int>::max(); for (const auto& s : sensors()) { segm = MIN(segm, s.numSegmentsEstimate()); } return segm; }
  int totalSegments() const { int cnt = 0; for (const auto& s : sensors()) { cnt += s.numSegmentsEstimate(); } return cnt; }
  int maxChannels() const { int max = 0; for (const auto& s : sensors()) { max = MAX(max, s.numChannels()); } return max; } 
  int minChannels() const { int min = std::numeric_limits<int>::max(); for (const auto& s : sensors()) { min = MIN(min, s.numChannels()); } return min; } 
  int totalChannels() const { int cnt = 0; for (const auto& s : sensors()) { cnt += s.numChannels(); } return cnt; } 

  double totalPower() const;

  int numStripsAcrossEstimate() const { return sensors().front().numStripsAcrossEstimate(); } // CUIDADO this assumes both sensors have the same number of sensing elements in the transversal direction - typically it is like that
  double pitch() const { return sensors().front().pitch(); }
int numSegmentsEstimate() const { return sensors().front().numSegmentsEstimate(); } // CUIDADO this assumes both sensors have the same number of sensing elements in the transversal direction - typically it is like that
  double stripLength() const { return sensors().front().stripLength(); }
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
  virtual const int diskSurface() const { return -1; }

  inline bool isPixelModule() const { return (moduleType().find(insur::type_pixel) != std::string::npos); }
  inline bool isTimingModule() const { return (moduleType().find(insur::type_timing) != std::string::npos); }

  bool couldHit(const XYZVector& direction, double zError) const;
  double trackCross(const XYZVector& PL, const XYZVector& PU) { return decorated().trackCross(PL, PU); }
  std::pair<XYZVector, HitType> checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir);
  int numHits() const { return numHits_; }
  void resetHits() { numHits_ = 0; }

  std::string summaryType() const;
  std::string summaryFullType() const;

  // OT CABLING
  void setBundle(OuterBundle* bundle) { bundle_ = bundle ; }
  const OuterBundle* getBundle() const { return bundle_; }
  const int isPositiveCablingSide() const;
  const int bundlePlotColor() const; 
  const int opticalChannelSectionPlotColor() const;
  const int powerChannelSectionPlotColor() const; 
  const OuterDTC* getDTC() const;
  const int dtcPlotColor() const;
  const int dtcPhiSectorRef() const;

  // IT CABLING
  void setPowerChain(PowerChain* powerChain) { powerChain_ = powerChain ; }
  const PowerChain* getPowerChain() const { return powerChain_; }
  void setPhiRefInPowerChain(int phiRefInPowerChain) { phiRefInPowerChain_ = phiRefInPowerChain; }
  const int getPhiRefInPowerChain() const { return phiRefInPowerChain_; }
  const int isPositiveZEnd() const;
  const bool isPositiveXSide() const;
  const int powerChainPlotColor() const;
  void setHvLine(HvLine* hvLine) { hvLine_ = hvLine ; }
  const HvLine* getHvLine() const { return hvLine_; }
  void setNumELinks(const int numELinks) { numELinks_ = numELinks; }
  const int numELinks() const { return numELinks_; }
  void setGBT(GBT* myGBT) { GBT_ = myGBT; }
  const GBT* getGBT() const { return GBT_; }
  const int gbtPlotColor() const;
  const InnerBundle* getInnerBundle() const;
  const int innerBundlePlotColor() const;
  const InnerDTC* getInnerDTC() const;
  const int innerDTCPlotColor() const;

protected:
  // Used to compute local parametrized spatial resolution.
  virtual const double calculateParameterizedResolutionLocalX(const double phi) const;
  virtual const double calculateParameterizedResolutionLocalY(const double theta) const;
  const double calculateParameterizedResolutionLocalAxis(const double fabsTanDeepAngle, const bool isLocalXAxis) const;

  MaterialObject materialObject_;
  Sensors sensors_;
  std::string subdetectorName_;
  int16_t subdetectorId_;
  mutable double cachedZError_ = -1.;
  mutable std::pair<double,double> cachedMinMaxEtaWithError_;
  XYZVector rAxis_;
  double tiltAngle_ = 0.;

  int numHits_ = 0;
  
  void clearSensorPolys() { for (auto& s : sensors_) s.clearPolys(); }
  ModuleCap* myModuleCap_ = nullptr;

private:
  // OT CABLING MAP
  OuterBundle* bundle_ = nullptr;
  // IT CABLING MAP
  PowerChain* powerChain_ = nullptr;
  int phiRefInPowerChain_;
  HvLine* hvLine_ = nullptr;
  int numELinks_;
  GBT* GBT_ =  nullptr;
  // The raw pointers are intended. DetectorModule is NOT owning the cabling resources.
  // All cabling ressources are owned by the CablingMap, which is a member variable of Tracker.
  // They get destructed when Tracker is destructed.
};



class BarrelModule : public DetectorModule, public Clonable<BarrelModule> {
public:
  Property<int16_t, AutoDefault> layer;
  int16_t ring() const { return (int16_t)myid(); }
  int16_t moduleRing() const { return ring(); }
  Property<int16_t, AutoDefault> rod;

  BarrelModule(Decorated* decorated) :
    DetectorModule(decorated)
  { setup(); }

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
  void accept(SensorGeometryVisitor& v) {
    v.visit(*this);
    v.visit(*(DetectorModule*)this);
    for (auto& s : sensors_) { s.accept(v); }
  }

  void setup() override {
    DetectorModule::setup();
    minPhi.setup([&](){

      double min = 0;
      // Module corners arranged normally or flipped:
      // 0 |-----|3      3|-----|0
      //   |     |   or   |     |
      // 1 |-----|2      2|-----|1
      //
      //              x (inter. point)
      // --> problem if absolute difference in phi betwwen barrel corners higher than phi
      if (!(fabs(basePoly().getVertex(0).Phi()-basePoly().getVertex(2).Phi())>=M_PI)) {

        min = MIN(basePoly().getVertex(0).Phi(), basePoly().getVertex(2).Phi());
      }
      // Module overlaps the crossline between -pi/2 & +pi/2 -> rotate by 180deg to calculate min
      else {

        Polygon3d<4> polygon = Polygon3d<4>(basePoly());
        polygon.rotateZ(M_PI);

        min = MIN(polygon.getVertex(0).Phi(), polygon.getVertex(2).Phi());

        // Shift by extra 180deg to get back to its original position (i.e. +2*pi with respect to the nominal position)
        min += M_PI;
      }
      // Return value in interval <-pi;+3*pi> instead of <-pi;+pi> to take into account the crossline at pi/2.
      return min;
    });
    maxPhi.setup([&](){

      double max = 0;
      // Module corners arranged normally or flipped:
      // 0 |-----|3      3|-----|0
      //   |     |   or   |     |
      // 1 |-----|2      2|-----|1
      //
      //              x (inter. point)
      // --> problem if absolute difference in phi betwwen barrel corners higher than phi
      if (!(fabs(basePoly().getVertex(0).Phi()-basePoly().getVertex(2).Phi())>=M_PI)) {

        max = MAX(basePoly().getVertex(0).Phi(), basePoly().getVertex(2).Phi());
      }
      // Module overlaps the crossline between -pi/2 & +pi/2 -> rotate by 180deg to calculate min
      else {

        Polygon3d<4> polygon = Polygon3d<4>(basePoly());
        polygon.rotateZ(M_PI);

        max = MAX(polygon.getVertex(0).Phi(), polygon.getVertex(2).Phi());

        // Shift by extra 180deg to get back to its original position (i.e. +2*pi with respect to the nominal position)
        max += M_PI;
      }
      // Return value in interval <-pi;+3*pi> instead of <-pi;+pi> to take into account the crossline at pi/2.
      return max;
    });

    nominalResolutionLocalX.setup([this]() {
	// only set up this if no model parameter specified
	//std::cout <<  "hasAnyResolutionLocalXParam() = " <<  hasAnyResolutionLocalXParam() << std::endl;

	if (!hasAnyResolutionLocalXParam()) {
	  //std::cout << "nominalResolutionLocalX and resolutionLocalXBarrel parameters are all unset. Use of default formulae." << std::endl;
	  double res = 0;
	  for (const Sensor& s : sensors()) res += pow(meanWidth() / s.numStripsAcross() / sqrt(12), 2);
	  return sqrt(res)/numSensors();
	}
	// if model parameters specified, return -1
	else return -1.0;
      });
    nominalResolutionLocalY.setup([this]() {
	// only set up this if no model parameters not specified
	if (!hasAnyResolutionLocalYParam()) {
	  //std::cout << "resolutionLocalY and resolutionLocalYBarrel parameters are all unset. Use of default formulae." << std::endl;
	    if (stereoRotation() != 0.) return nominalResolutionLocalX() / sin(stereoRotation());
	    else {
	      return length() / maxSegments() / sqrt(12); // NOTE: not combining measurements from both sensors. The two sensors are closer than the length of the longer sensing element, making the 2 measurements correlated. considering only the best measurement is then a reasonable approximation (since in case of a PS module the strip measurement increases the precision by only 0.2% and in case of a 2S the sensors are so close that they basically always measure the same thing)
	    }
	  }
	// if model parameters specified, return -1
	else return -1.0;
      });
  }

  void build();


  //double maxZ() const { return MAX(basePoly().getVertex(0).Z(), basePoly().getVertex(2).Z()); } 
  //double minZ() const { return MIN(basePoly().getVertex(0).Z(), basePoly().getVertex(2).Z()); } 
  //double maxR() const { return MAX(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }
  //double minR() const { return center().Rho(); }//MIN(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }

  virtual ModuleSubdetector subdet() const { return BARREL; }

  PosRef posRef() const { return (PosRef){ subdetectorId(), (side() > 0 ? ring() : -ring()), layer(), rod() }; }
  TableRef tableRef() const { return (TableRef){ subdetectorName(), layer(), ring() }; }
  UniRef uniRef() const { return UniRef{ subdetectorName(), layer(), ring(), rod(), side() }; }
};



class EndcapModule : public DetectorModule, public Clonable<EndcapModule> {
public:
  Property<int16_t, AutoDefault> disk;
  Property<int16_t, AutoDefault> ring;
  int16_t moduleRing() const { return ring(); };
  int16_t blade() const { return (int16_t)myid(); } // CUIDADO Think of a better name!
  int16_t side() const { return (int16_t)signum(center().Z()); }
  Property<int, AutoDefault> endcapDiskSurface;
  const int diskSurface() const override { return endcapDiskSurface(); }

  EndcapModule(Decorated* decorated) :
    DetectorModule(decorated)
  { setup(); }


  void setup() override {
    DetectorModule::setup();
    minPhi.setup([&](){

      double min = 0;
      // Module corners arranged normally:
      // 0 |-----|3
      //   |     |
      // 1 |-----|2
      //
      //      x (inter. point)
      if (basePoly().getVertex(1).Phi()<=basePoly().getVertex(0).Phi() &&
          basePoly().getVertex(0).Phi()<=basePoly().getVertex(3).Phi() &&
          basePoly().getVertex(3).Phi()<=basePoly().getVertex(2).Phi()) {

        min=minget2(basePoly().begin(), basePoly().end(), &XYZVector::Phi);
      }
      // Module corners flipped:
      // 3 |-----|0
      //   |     |
      // 2 |-----|1
      //
      //      x (inter. point)
      else if (basePoly().getVertex(2).Phi()<=basePoly().getVertex(3).Phi() &&
               basePoly().getVertex(3).Phi()<=basePoly().getVertex(0).Phi() &&
               basePoly().getVertex(0).Phi()<=basePoly().getVertex(1).Phi()){

        min=minget2(basePoly().begin(), basePoly().end(), &XYZVector::Phi);
      }
      // Module overlaps the crossline between -pi/2 & +pi/2 -> rotate by 180deg to calculate min
      else {

        Polygon3d<4> polygon = Polygon3d<4>(basePoly());
        polygon.rotateZ(M_PI);

        min=minget2(polygon.begin(), polygon.end(), &XYZVector::Phi);

        // Normal arrangement or flipped arrangement
        if (polygon.getVertex(1).Phi()<0 || polygon.getVertex(2).Phi()<0) {

          // Shift by extra 180deg to get back to its original position (i.e. +2*pi with respect to the nominal position)
          min += M_PI;
        }
        else logERROR("Endcap module min calculation failed - algorithm problem. Check algorithm!");

      }
      // Return value in interval <-pi;+3*pi> instead of <-pi;+pi> to take into account the crossline at pi/2.
      return min;
    });
    maxPhi.setup([&](){

      double max = 0;
      // Module corners arranged normally:
      // 0 |-----|3
      //   |     |
      // 1 |-----|2
      //
      //      x (inter. point)
      if (basePoly().getVertex(1).Phi()<=basePoly().getVertex(0).Phi() &&
          basePoly().getVertex(0).Phi()<=basePoly().getVertex(3).Phi() &&
          basePoly().getVertex(3).Phi()<=basePoly().getVertex(2).Phi()) {

        max=maxget2(basePoly().begin(), basePoly().end(), &XYZVector::Phi);
      }
      // Module corners flipped:
      // 3 |-----|0
      //   |     |
      // 2 |-----|1
      //
      //      x (inter. point)
      else if (basePoly().getVertex(2).Phi()<=basePoly().getVertex(3).Phi() &&
               basePoly().getVertex(3).Phi()<=basePoly().getVertex(0).Phi() &&
               basePoly().getVertex(0).Phi()<=basePoly().getVertex(1).Phi()){

        max=maxget2(basePoly().begin(), basePoly().end(), &XYZVector::Phi);
      }
      // Module overlaps the crossline between -pi/2 & +pi/2 -> rotate by 180deg to calculate max.
      else {

        Polygon3d<4> polygon = Polygon3d<4>(basePoly());
        polygon.rotateZ(M_PI);

        max=maxget2(polygon.begin(), polygon.end(), &XYZVector::Phi);

        // Normal arrangement or flipped arrangement
        if (polygon.getVertex(1).Phi()<0 || polygon.getVertex(2).Phi()<0) {

          // Shift by extra 180deg to get back to its original position (i.e. +2*pi with respect to the nominal position)
          max += M_PI;
        }
        else logERROR("Endcap module max calculation failed - algorithm problem. Check algorithm!");
      }
      // Return value in interval <-pi;+3*pi> instead of <-pi;+pi> to take into account the crossline at pi/2.
      return max;
    });
    nominalResolutionLocalX.setup([this]() {
	// only set up this if no model parameter specified
	//std::cout <<  "hasAnyResolutionLocalXParam() = " <<  hasAnyResolutionLocalXParam() << std::endl;
	if (!hasAnyResolutionLocalXParam()) {
	  //std::cout << "nominalResolutionLocalX and resolutionLocalXEndcap parameters are all unset. Use of default formulae." << std::endl;
	    double res = 0;
	    for (const Sensor& s : sensors()) res += pow(meanWidth() / s.numStripsAcross() / sqrt(12), 2);
	    return sqrt(res)/numSensors();
	  }
	// if model parameters specified, return -1
	else return -1.0;
      });
    nominalResolutionLocalY.setup([this]() {
	// only set up this if no model parameters not specified
	if (!hasAnyResolutionLocalYParam()) {
	  //std::cout << "resolutionLocalY and resolutionLocalYEndcap parameters are all unset. Use of default formulae." << std::endl;
	    if (stereoRotation() != 0.) return nominalResolutionLocalX() / sin(stereoRotation());
	    else {
	      return length() / maxSegments() / sqrt(12); // NOTE: not combining measurements from both sensors. The two sensors are closer than the length of the longer sensing element, making the 2 measurements correlated. considering only the best measurement is then a reasonable approximation (since in case of a PS module the strip measurement increases the precision by only 0.2% and in case of a 2S the sensors are so close that they basically always measure the same thing)
	    }
	  }
	// if model parameters specified, return -1
	else return -1.0;
      });
  }

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
  void accept(SensorGeometryVisitor& v) {
    v.visit(*this);
    v.visit(*(DetectorModule*)this);
    for (auto& s : sensors_) { s.accept(v); }
  }

  //double minZ() const { return center().Z(); } // CUIDADO not accounting for sensor placement
  //double maxZ() const { return center().Z(); } // ditto here
  //double maxR() const { return MAX(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }
  //double minR() const { XYZVector side[2];
  //                      std::partial_sort_copy(basePoly().begin(), basePoly().end(), std::begin(side), std::end(side), [](const XYZVector& v1, const XYZVector& v2) { return v1.Rho() < v2.Rho(); });
  //                      return ((side[0]+side[1])/2).Rho(); }


  virtual ModuleSubdetector subdet() const { return ENDCAP; }

  PosRef posRef() const { return (PosRef){ subdetectorId(), (side() > 0 ? disk() : -disk()), ring(), blade() }; }
  TableRef tableRef() const { return (TableRef){ subdetectorName(), disk(), ring() }; }
  UniRef uniRef() const { return UniRef{ subdetectorName(), disk(), ring(), blade(), side() }; }
};


// ===================================================================================================================================
//
#endif

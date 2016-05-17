#ifndef DETECTOR_MODULE_H
#define DETECTOR_MODULE_H

#include <boost/ptr_container/ptr_vector.hpp>

#include "Sensor.h"
#include "ModuleBase.h"
#include "GeometricModule.h"
#include "CoordinateOperations.h"
#include "Visitable.h"
#include "MaterialObject.h"
#include "MessageLogger.h"

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

class DetectorModule;

using insur::ModuleCap;
using material::ElementsVector;

/*
 * Detector module class - base to Barrel & Endcap classes
 */
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

  ModuleCap* m_myModuleCap = nullptr; //!< Module materials assigned to its active part

public:

  //! Set materials, so called module-cap to a module
  void setModuleCap(ModuleCap* newCap);

  //! Get materials, i.e. module-cap
  ModuleCap* getModuleCap() { return m_myModuleCap;}

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
  
  //! Constructor
  DetectorModule(Decorated* decorated);

  //! Destructor
  ~DetectorModule();

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

  double resolutionEquivalentZ   (double hitRho, double trackR, double trackCotgTheta) const;
  double resolutionEquivalentRPhi(double hitRho, double trackR) const;

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
  int minSegments() const { int segm = 999999; for (const auto& s : sensors()) { segm = MIN(segm, s.numSegments()); } return segm; }
  int totalSegments() const { int cnt = 0; for (const auto& s : sensors()) { cnt += s.numSegments(); } return cnt; }
  int maxChannels() const { int max = 0; for (const auto& s : sensors()) { max = MAX(max, s.numChannels()); } return max; } 
  int minChannels() const { int min = 999999; for (const auto& s : sensors()) { min = MIN(min, s.numChannels()); } return min; } 
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


  BarrelModule(Decorated* decorated);

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
  }

  void build();


  //double maxZ() const { return MAX(basePoly().getVertex(0).Z(), basePoly().getVertex(2).Z()); } 
  //double minZ() const { return MIN(basePoly().getVertex(0).Z(), basePoly().getVertex(2).Z()); } 
  //double maxR() const { return MAX(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }
  //double minR() const { return center().Rho(); }//MIN(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }

  virtual ModuleSubdetector subdet() const { return BARREL; }

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

  EndcapModule(Decorated* decorated);

  void setup() override {
    DetectorModule::setup();
    minPhi.setup([&](){

      double min = 0;

      // Rectangular end-cap modules
      if (this->shape()==ModuleShape::RECTANGULAR) {

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
      }

      // Wedge-shaped modules
      //else {
      //}

      // Return value in interval <-pi;+3*pi> instead of <-pi;+pi> to take into account the crossline at pi/2.
      return min;
    });
    maxPhi.setup([&](){

      double max = 0;

      // Rectangular end-cap modules
      if (this->shape()==ModuleShape::RECTANGULAR) {

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
      }

      // Wedge-shaped modules
      //else {
      //}

      // Return value in interval <-pi;+3*pi> instead of <-pi;+pi> to take into account the crossline at pi/2.
      return max;
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

  //double minZ() const { return center().Z(); } // CUIDADO not accounting for sensor placement
  //double maxZ() const { return center().Z(); } // ditto here
  //double maxR() const { return MAX(basePoly().getVertex(0).Rho(), basePoly().getVertex(2).Rho()); }
  //double minR() const { XYZVector side[2];
  //                      std::partial_sort_copy(basePoly().begin(), basePoly().end(), std::begin(side), std::end(side), [](const XYZVector& v1, const XYZVector& v2) { return v1.Rho() < v2.Rho(); });
  //                      return ((side[0]+side[1])/2).Rho(); }


  virtual ModuleSubdetector subdet() const { return ENDCAP; }

  PosRef posRef() const { return (PosRef){ cntId(), (side() > 0 ? disk() : -disk()), ring(), blade() }; }
  TableRef tableRef() const { return (TableRef){ cntName(), disk(), ring() }; }
  UniRef uniRef() const { return UniRef{ cntName(), disk(), ring(), blade(), side() }; }
};


// ===================================================================================================================================
//
#endif

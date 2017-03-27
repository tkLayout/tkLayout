#ifndef RODPAIR_H
#define RODPAIR_H

#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <algorithm>

#include <boost/ptr_container/ptr_vector.hpp>

#include "global_funcs.hh"
#include "Property.hh"
#include "Module.hh"
#include "MessageLogger.hh"
#include "Visitable.hh"

#include <iostream>

using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;


struct TiltedModuleSpecs {
  double r, z, gamma;
  bool valid() const {
    return r > 0.0 && fabs(gamma) <= 2*M_PI;
  }
};


typedef vector<unique_ptr<BarrelModule>> RodTemplate;

class RodPair : public PropertyObject, public Buildable, public Identifiable<int>, public Visitable {
public:
  typedef PtrVector<BarrelModule> Container;
protected:
  Container zPlusModules_, zMinusModules_;
  MaterialObject materialObject_;
public:
  enum class BuildDir { RIGHT = 1, LEFT = -1 };
  enum class StartZMode { MODULECENTER, MODULEEDGE };

private:
  void clearComputables();
public:
  Property<StartZMode, Default> startZMode;
  //Property<double, NoDefault> maxZ;
  Property<double, Computable> maxZ;
  Property<double, Computable> minZ, maxR, minR;
  Property<double, Computable> minZwithHybrids, maxZwithHybrids, minRwithHybrids, maxRwithHybrids;
  ReadonlyProperty<double, Computable> maxModuleThickness;
  Property<bool, Default> beamSpotCover;
  Property<bool, NoDefault> isOuterRadiusRod;

  RodPair() :
      materialObject_(MaterialObject::ROD),
      startZMode("startZMode", parsedAndChecked(), StartZMode::MODULECENTER),
	beamSpotCover("beamSpotCover", parsedAndChecked(), true)
  {}

  void setup() {
    minZ       .setup([&]() { return minget2(zMinusModules_.begin(), zMinusModules_.end(), &Module::minZ); }); // we want the minZ so we don't bother with scanning the zPlus vector
    maxZ .setup([&]() { return maxget2(zPlusModules_.begin(), zPlusModules_.end(), &Module::maxZ); });
    minR       .setup([&]() { return minget2(zPlusModules_.begin(), zPlusModules_.end(), &Module::minR); }); // min and maxR can be found by just scanning the zPlus vector, since the rod pair is symmetrical in R
    maxR       .setup([&]() { return maxget2(zPlusModules_.begin(), zPlusModules_.end(), &Module::maxR); });
    maxModuleThickness.setup([&]() { return maxget2(zPlusModules_.begin(), zPlusModules_.end(), &Module::thickness); });

    maxZwithHybrids       .setup([&]() { return maxget2(zPlusModules_.begin(), zPlusModules_.end(), &Module::maxZwithHybrids); });
    minZwithHybrids       .setup([&]() { return minget2(zMinusModules_.begin(), zMinusModules_.end(), &Module::minZwithHybrids); });
    minRwithHybrids       .setup([&]() { return minget2(zPlusModules_.begin(), zPlusModules_.end(), &Module::minRwithHybrids); });
    maxRwithHybrids       .setup([&]() { return maxget2(zPlusModules_.begin(), zPlusModules_.end(), &Module::maxRwithHybrids); });
  }
  
  virtual double thickness() const = 0;
  virtual bool isTilted() const = 0;

  int numModules() const { return zPlusModules_.size() + zMinusModules_.size(); }
  int numModulesSide(int side) const { return side >= 0 ? zPlusModules_.size() : zMinusModules_.size(); }

  void translate(const XYZVector& translation);
  void translateR(double radius);
  void rotateZ(double angle);

  double Phi() const {
    double phi;
    if (zPlusModules_.size() != 0) phi = zPlusModules_.front().center().Phi();
    else if (zMinusModules_.size() != 0) phi = zMinusModules_.front().center().Phi();
    else phi = std::numeric_limits<double>::quiet_NaN();
    return phi;
  }

  void cutAtEta(double eta);

  const std::pair<const Container&,const Container&> modules() const { return std::pair<const Container&,const Container&>(zPlusModules_,zMinusModules_); }

  void removeModules() { zMinusModules_.erase_if([](DetectorModule& m) { return (m.removeModule()); }); zPlusModules_.erase_if([](DetectorModule& m) { return (m.removeModule()); }); }
  
  void accept(GeometryVisitor& v) {
    v.visit(*this); 
    for (auto& m : zPlusModules_) { m.accept(v); }
    for (auto& m : zMinusModules_) { m.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const {
    v.visit(*this); 
    for (const auto& m : zPlusModules_) { m.accept(v); }
    for (const auto& m : zMinusModules_) { m.accept(v); }
  }
  void accept(SensorGeometryVisitor& v) {
    v.visit(*this); 
    for (auto& m : zPlusModules_) { m.accept(v); }
    for (auto& m : zMinusModules_) { m.accept(v); }
  }

  const MaterialObject& materialObject() const;
};

class StraightRodPair : public RodPair, public Clonable<StraightRodPair> {

  // Templated because they need to work both with forward and reverse iterators (mezzanines are built right to left and the rodTemplate vector is iterated backwards)
  double computeNextZ(double newDsLength, double newDsDistance, double lastDsDistance, double lastZ, BuildDir direction, int parity);
  template<typename Iterator> vector<double> computeZList(Iterator begin, Iterator end, double startZ, BuildDir direction, int smallParity, bool fixedStartZ);
  template<typename Iterator> pair<vector<double>, vector<double>> computeZListPair(Iterator begin, Iterator end, double startZ, int recursionCounter);
  void buildModules(Container& modules, const RodTemplate& rodTemplate, const vector<double>& posList, BuildDir direction, bool isPlusBigDeltaRod, int parity, int side);
  void buildFull(const RodTemplate& rodTemplate, bool isPlusBigDeltaRod); 
  void buildMezzanine(const RodTemplate& rodTemplate, bool isPlusBigDeltaRod); 

public:
 
  RangeProperty<std::vector<double> > forbiddenRange;
  Property<double, NoDefault> smallDelta;
  Property<double, NoDefault> minBuildRadius;
  Property<double, NoDefault> maxBuildRadius;

  Property<double, Default> zOverlap;
  Property<double, NoDefault> zError;
  Property<int, NoDefault> zPlusParity;
  Property<int, NoDefault> buildNumModules;
  Property<bool, Default> mezzanine;
  Property<double, NoDefault> startZ;
  Property<bool, Default> compressed;
  Property<bool, Default> allowCompressionCuts;

  Property<bool, Default> isFlatPart;

  PropertyNode<int> ringNode;
  
  StraightRodPair() :
              forbiddenRange      ("forbiddenRange"      , parsedOnly()),
              zOverlap            ("zOverlap"            , parsedAndChecked() , 1.),
	      zError              ("zError"              , parsedAndChecked()),
	      zPlusParity         ("smallParity"         , parsedOnly()),
              mezzanine           ("mezzanine"           , parsedOnly(), false),
              startZ              ("startZ"              , parsedOnly()),
              compressed          ("compressed"          , parsedOnly(), true),
              allowCompressionCuts("allowCompressionCuts", parsedOnly(), true),
	      ringNode            ("Ring"                , parsedOnly()),
	      isFlatPart          ("isFlatPart"          , parsedOnly(), false)
  {}


  double thickness() const override { return smallDelta()*2. + maxModuleThickness(); }
  bool isTilted() const override { return false; }

  void check() override;
  void build(const RodTemplate& rodTemplate, bool isPlusBigDeltaRod);

  std::set<int> solveCollisionsZPlus();
  std::set<int> solveCollisionsZMinus();
  void compressToZ(double z);

  double thetaEnd_REAL() const {
    double thetaEnd;

    if (zPlusModules_.empty()) { thetaEnd = M_PI/2.; }
    else {
      // findMaxZModule as a function
      auto lastMod = zPlusModules_.back();

      double dsDistance = lastMod.dsDistance();
      double lastR = lastMod.center().Rho();
      
      double rH2ppUP = lastR + 0.5 * dsDistance;  // WARNING !!! FOR THE MOMENT, DOESN T TAKE MODULE WIDTH INTO ACCOUNT, SHOULD BE CHANGED ?

      thetaEnd = atan(rH2ppUP / (lastMod.planarMaxZ()));

      /*std::cout << "lastMod.center().Rho() = " << lastMod.center().Rho() << std::endl;
      std::cout << "lastMod.dsDistance() = " << lastMod.dsDistance() << std::endl;
      std::cout << "lastMod.thickness() = " << lastMod.thickness() << std::endl;
      std::cout << "rH2ppUP = " << rH2ppUP << std::endl;
      std::cout << "thetaEnd = " << thetaEnd << std::endl;
   
      std::cout << "lastMod.planarMaxR() = " << lastMod.planarMaxR() << std::endl;
      std::cout << "lastMod.planarMaxZ() - zOverlap() = " << lastMod.planarMaxZ() - zOverlap() << std::endl;
      double thetaEnd2 = atan(lastMod.planarMaxR() / (lastMod.planarMaxZ() - zOverlap()));
      std::cout << "thetaEnd2 = " << thetaEnd2 << std::endl;*/
    }
    return thetaEnd;
  }

};



class TiltedRodPair : public RodPair, public Clonable<TiltedRodPair> {
 
  void buildModules(Container& modules, const RodTemplate& rodTemplate, const vector<TiltedModuleSpecs>& tmspecs, BuildDir direction, bool flip);

 public :

  double thickness() const override { std::cerr << "thickness() for tilted rods gives incorrect results as it is calculated as maxR()-minR()\n"; return maxR() - minR(); }
  bool isTilted() const override { return true; }
  void check() override;
  void build(const RodTemplate& rodTemplate, const std::vector<TiltedModuleSpecs>& tmspecs, bool flip);

  

}; 



#endif

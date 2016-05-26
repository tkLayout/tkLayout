
#include "DetectorModule.h"
#include "ModuleCap.h"

define_enum_strings(SensorLayout) = { "nosensors", "mono", "pt", "stereo" };
define_enum_strings(ZCorrelation) = { "samesegment", "multisegment" };
define_enum_strings(ReadoutType)  = { "strip", "pixel", "pt" };
define_enum_strings(ReadoutMode)  = { "binary", "cluster" };

//
// Constructor - specify unique id, geometry module defining shape & parse geometry config file using boost property tree & read-in module parameters
//
DetectorModule::DetectorModule(int id, Decorated* decorated, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty) :
 Decorator<GeometricModule>(decorated),
 m_materialObject(MaterialObject::MODULE),
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
{
  // Set the geometry config parameters
  this->myid(id);
  this->store(treeProperty);
  if (nodeProperty.count(id)>0) this->store(nodeProperty.at(id));
}

//
// Constructor - specify unique id, geometry module defining shape & parse geometry config file using boost property tree & read-in module parameters
//
DetectorModule::DetectorModule(int id, Decorated* decorated, const PropertyTree& treeProperty) :
 Decorator<GeometricModule>(decorated),
 m_materialObject(MaterialObject::MODULE),
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
{
  // Set the geometry config parameters
  this->myid(id);
  this->store(treeProperty);
}

//
// Destructor
//
DetectorModule::~DetectorModule()
{
  if (m_myModuleCap!=nullptr) delete m_myModuleCap;
}

//
// Build detector module with its sensors
//
void DetectorModule::build() {

  check();
  if (!decorated().builtok()) {

    decorated().store(propertyTree());
    decorated().build();
  }
  if (numSensors() > 0) {

    for (int iSensor=1; iSensor<=numSensors(); iSensor++) {
      Sensor* s = GeometryFactory::make<Sensor>(iSensor, this, sensorNode, propertyTree());
      s->build();

      m_sensors.push_back(s);
      m_materialObject.sensorChannels[iSensor+1]=s->numChannels();
    }
  } else {

    // Fake sensor to avoid defensive programming when iterating over the sensors and the module is empty
    Sensor* s = GeometryFactory::make<Sensor>(1, this);
    s->build();

    m_sensors.push_back(s);
  }

  m_materialObject.store(propertyTree());
  m_materialObject.build();
}

//
// Setup: link lambda functions to various DetectorModule related properties (use setup functions for ReadOnly Computable properties)
//
void DetectorModule::setup() {

  resolutionLocalX.setup([this]() { 
    double res = 0;
    for (const Sensor& s : sensors()) res += pow(meanWidth() / s.numStripsAcross() / sqrt(12), 2);
    return sqrt(res)/numSensors();
  });

  resolutionLocalY.setup([this]() {
    if (stereoRotation() != 0.) return resolutionLocalX() / sin(stereoRotation());
    else {
      return length() / maxSegments() / sqrt(12); // NOTE: not combining measurements from both sensors. The two sensors are closer than the length of the longer sensing element, making the 2 measurements correlated. considering only the best measurement is then a reasonable approximation (since in case of a PS module the strip measurement increases the precision by only 0.2% and in case of a 2S the sensors are so close that they basically always measure the same thing)
    }
  });

  maxZ.setup([&]() {return maxget2(m_sensors.begin(), m_sensors.end(), &Sensor::maxZ); });
  minZ.setup([&]() {return minget2(m_sensors.begin(), m_sensors.end(), &Sensor::minZ); });
  maxR.setup([&]() {return maxget2(m_sensors.begin(), m_sensors.end(), &Sensor::maxR); });
  minR.setup([&]() {return minget2(m_sensors.begin(), m_sensors.end(), &Sensor::minR); });

  planarMaxZ.setup([&]() { return CoordinateOperations::computeMaxZ(basePoly()); });
  planarMinZ.setup([&]() { return CoordinateOperations::computeMinZ(basePoly()); });
  planarMaxR.setup([&]() { return CoordinateOperations::computeMaxR(basePoly()); });
  planarMinR.setup([&]() { return CoordinateOperations::computeMinR(basePoly()); });

  maxSegments.setup([&]()  { int segm = 0;                               for (const auto& s : m_sensors) { segm = MAX(segm, s.numSegments()); } return segm; });
  minSegments.setup([&]()  { int segm = std::numeric_limits<int>::max(); for (const auto& s : m_sensors) { segm = MIN(segm, s.numSegments()); } return segm; });
  totalSegments.setup([&](){ int cnt  = 0;                               for (const auto& s : m_sensors) { cnt += s.numSegments(); } return cnt; });
  maxChannels.setup([&]()  { int max  = 0;                               for (const auto& s : m_sensors) { max = MAX(max, s.numChannels()); } return max; });
  minChannels.setup([&]()  { int min  = std::numeric_limits<int>::max(); for (const auto& s : m_sensors) { min = MIN(min, s.numChannels()); } return min; });
  totalChannels.setup([&](){ int cnt  = 0;                               for (const auto& s : m_sensors) { cnt += s.numChannels(); } return cnt; });


  // Set the parent for the sensors once again (in case the module's been cloned) -> setup method called after cloning!
  for (Sensor& s : m_sensors) s.parent(this);
}

//
// Could track at this direction hit the module?
//
bool DetectorModule::couldHit(const XYZVector& direction, double zError) const {

  double eta       = direction.Eta();
  double phi       = direction.Phi();
  double shiftPhi  = phi + 2*M_PI;
  bool   withinEta = false;
  bool   withinPhi = false;

  // Eta region covered by module
  if (eta > minMaxEtaWithError(zError).first && eta < minMaxEtaWithError(zError).second) withinEta = true;

  // Phi region is from <-pi;+3*pi> due to crossline at +pi -> need to check phi & phi+2*pi
  if ( (phi     >=minPhi() && phi     <=maxPhi()) ||
       (shiftPhi>=minPhi() && shiftPhi<=maxPhi()) ) withinPhi = true;

  // Checking that hit within a module region works for barrel-type modules only!!!
  if (this->shape()==ModuleShape::RECTANGULAR) return (withinEta && withinPhi);
  // ATTENTION: For wedge shaped modules, min, max procedure will not work correctly -> return true to avoid errors --> will be implemented in the future
  else return true;
}

//
// Module R-Phi-resolution calculated as for a barrel-type module -> transform it to true orientation (rotation by theta angle, skew, tilt)
//
double DetectorModule::resolutionEquivalentRPhi(double hitRho, double trackR) const {

  // Parameters
  double A = hitRho/(2*trackR); // r_i / 2R
  double B = A/sqrt(1-A*A);

  // All modules & its resolution propagated to the resolution of a virtual barrel module (endcap is a tilted module by 90 degrees, barrel is tilted by 0 degrees)
  double resolution = sqrt(pow((B*sin(skewAngle())*cos(tiltAngle()) + cos(skewAngle())) * resolutionLocalX(),2) + pow(B*sin(tiltAngle()) * resolutionLocalY(),2));

  // Return calculated resolution (resolutionLocalX is intrinsic resolution along R-Phi for barrel module)
  return resolution; //resolutionLocalX();
}

//
// Module Z-resolution calculated as for a barrel-type module -> transform it to true orientation (rotation by theta angle, skew, tilt)
//
double DetectorModule::resolutionEquivalentZ(double hitRho, double trackR, double trackCotgTheta) const {

  // Parameters
  double A = hitRho/(2*trackR); 
  double D = trackCotgTheta/sqrt(1-A*A);

  // All modules & its resolution propagated to the resolution of a virtual barrel module (endcap is a tilted module by 90 degrees, barrel is tilted by 0 degrees)
  double resolution = sqrt(pow(((D*cos(tiltAngle()) + sin(tiltAngle()))*sin(skewAngle())) * resolutionLocalX(),2) + pow((D*sin(tiltAngle()) + cos(tiltAngle())) * resolutionLocalY(),2));

  // Return calculated resolution (resolutionLocalY is intrinsic resolution along Z for barrel module)
  return resolution; //resolutionLocalY();
}

void DetectorModule::setModuleCap(ModuleCap* newCap)
{
  if (m_myModuleCap!=nullptr) delete m_myModuleCap;
  m_myModuleCap = newCap ;
}

void DetectorModule::mirrorZ() {
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

std::pair<double, double> DetectorModule::minMaxEtaWithError(double zError) const {
  if (cachedZError_ != zError) {
    cachedZError_ = zError;
    double eta1 = (XYZVector(0., maxR(), maxZ() + zError)).Eta();
    double eta2 = (XYZVector(0., minR(), minZ() - zError)).Eta();
    double eta3 = (XYZVector(0., minR(), maxZ() + zError)).Eta();
    double eta4 = (XYZVector(0., maxR(), minZ() - zError)).Eta();
    cachedMinMaxEtaWithError_ = std::minmax({eta1, eta2, eta3, eta4});
    //cachedMinMaxEtaWithError_ = std::make_pair(MIN(eta1, eta2), MAX(eta1, eta2));
  }
  return cachedMinMaxEtaWithError_;
}

double DetectorModule::stripOccupancyPerEventBarrel() const {
  double rho = center().Rho()/10.;
  double theta = center().Theta();
  double myOccupancyBarrel=(1.63e-4)+(2.56e-4)*rho-(1.92e-6)*rho*rho;
  double factor = fabs(sin(theta))*2; // 2 is a magic adjustment factor
  double dphideta = phiAperture() * etaAperture();
  double minNSegments = minSegments();
  int numStripsAcross = sensors().begin()->numStripsAcross();
  double modWidth = (maxWidth() + minWidth())/2.;

  double occupancy = myOccupancyBarrel / factor / (90/1e3) * (dphideta / minNSegments) * (modWidth / numStripsAcross);

  return occupancy;
}

double DetectorModule::stripOccupancyPerEventEndcap() const {
  double rho = center().Rho()/10.;
  double theta = center().Theta();
  double z = center().Z()/10.;
  double myOccupancyEndcap = (-6.20e-5)+(1.75e-4)*rho-(1.08e-6)*rho*rho+(1.50e-5)*(z);
  double factor=fabs(cos(theta))*2; // 2 is a magic adjustment factor
  double dphideta = phiAperture() * etaAperture();
  double minNSegments = minSegments();
  int numStripsAcross = sensors().begin()->numStripsAcross();
  double modWidth = (maxWidth() + minWidth())/2.;

  double occupancy = myOccupancyEndcap / factor / (90/1e3) * (dphideta / minNSegments) * (modWidth / numStripsAcross);

  return occupancy;
}

double DetectorModule::stripOccupancyPerEvent() const {
  if (fabs(tiltAngle()) < 1e-3) return stripOccupancyPerEventBarrel();
  else if (fabs(tiltAngle()) - M_PI/2. < 1e-3) return stripOccupancyPerEventEndcap();
  else return stripOccupancyPerEventBarrel()*pow(cos(tiltAngle()),2) + stripOccupancyPerEventEndcap()*pow(sin(tiltAngle()),2);
};

double DetectorModule::geometricEfficiency() const {
  double inefficiency = fabs(dsDistance() / (zCorrelation()==SAMESEGMENT ? outerSensor().stripLength() : length()) / tan(center().Theta()+tiltAngle())); // fabs prevents the inefficiency from becoming negative when theta+tilt > 90 deg (meaning the geometrically inefficient area is on the other end of the module)
  return 1-inefficiency;
}

double DetectorModule::effectiveDsDistance() const {
  if (fabs(tiltAngle()) < 1e-3) return dsDistance();
  else return dsDistance()*sin(center().Theta())/sin(center().Theta()+tiltAngle());
}

std::pair<XYZVector, HitType> DetectorModule::checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir) {

  HitType ht = HitType::NONE;
  XYZVector gc; // global coordinates of the hit

  if (numSensors() == 1) {

    auto segm = innerSensor().checkHitSegment(trackOrig, trackDir);
    // <SMe>The following line used to return HitType::BOTH. Changing to INNER in order to avoid double hit counting</SMe>
    if (segm.second > -1) { gc = segm.first; ht = HitType::INNER; } 
  }
  else {

    auto inSegm  = innerSensor().checkHitSegment(trackOrig, trackDir);
    auto outSegm = outerSensor().checkHitSegment(trackOrig, trackDir);
    if (inSegm.second > -1 && outSegm.second > -1) { 
      gc = inSegm.first; // in case of both sensors are hit, the inner sensor hit coordinate is returned
      ht = ((zCorrelation() == SAMESEGMENT && (inSegm.second / (maxSegments()/minSegments()) == outSegm.second)) || zCorrelation() == MULTISEGMENT) ? HitType::STUB : HitType::BOTH;
    }
    else if (inSegm.second > -1)  { gc = inSegm.first;  ht = HitType::INNER; }
    else if (outSegm.second > -1) { gc = outSegm.first; ht = HitType::OUTER; }
  }
  //basePoly().isLineIntersecting(trackOrig, trackDir, gc); // this was just for debug
  if (ht != HitType::NONE) m_numHits++;
  return std::make_pair(gc, ht);
};

BarrelModule::BarrelModule(int id, Decorated* decorated, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty) :
 DetectorModule(id, decorated, nodeProperty, treeProperty)
{}

//
// Build barrel module
//
void BarrelModule::build()
{
  try {

    // Build overriden -> call base class build automatically
    DetectorModule::build();

    // Initialize module: rotations & translations called at higher level for the whole set of modules
    decorated().rotateY(M_PI/2);
    m_rAxis = normal();
    m_tiltAngle = 0.;
    m_skewAngle = 0.;
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }

  cleanup();
  builtok(true);
}

//
// Setup: link lambda functions to various BarrelModule (and DetectorModule - because overriden) related properties (use setup functions for ReadOnly Computable properties)
//
void BarrelModule::setup()
{

  // Setup overriden -> call base class setup automatically
  DetectorModule::setup();

  // Set other lambda functions
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

//
// GeometryVisitor pattern -> barrel module visitable
//
void BarrelModule::accept(GeometryVisitor& v) {
  v.visit(*this);
  v.visit(*(DetectorModule*)this);
  decorated().accept(v);
}

//
// GeometryVisitor pattern -> barrel module visitable (const. option)
//
void BarrelModule::accept(ConstGeometryVisitor& v) const {
  v.visit(*this);
  v.visit(*(const DetectorModule*)this);
  decorated().accept(v);
}

//
// Constructor - parse geometry config file using boost property tree & read-in module parameters, specify unique id
//
EndcapModule::EndcapModule(int id, Decorated* decorated, const PropertyTree& treeProperty) :
 DetectorModule(id, decorated, treeProperty)
{}

//
// Build endcap module with its sensors
//
void EndcapModule::build()
{
  try {

    // Build overriden -> call base class build automatically
    DetectorModule::build();

    // Initialize module: rotations & translations called at higher level for the whole set of modules
    m_rAxis = (basePoly().getVertex(0) + basePoly().getVertex(3)).Unit();
    m_tiltAngle = M_PI/2.;
    m_skewAngle = 0.;
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
  cleanup();
  builtok(true);
}

//
// Setup: link lambda functions to various EndcapModule (and DetectorModule - because overriden) related properties (use setup functions for ReadOnly Computable properties)
//
void EndcapModule::setup()
{

  // Setup overriden -> call base class setup automatically
  DetectorModule::setup();

  // Set other lambda functions
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

//
// GeometryVisitor pattern -> endcap module visitable
//
void EndcapModule::accept(GeometryVisitor& v) {
  v.visit(*this);
  v.visit(*(DetectorModule*)this);
  decorated().accept(v);
}

//
// GeometryVisitor pattern -> endcap module visitable (const. option)
//
void EndcapModule::accept(ConstGeometryVisitor& v) const {
  v.visit(*this);
  v.visit(*(const DetectorModule*)this);
  decorated().accept(v);
}


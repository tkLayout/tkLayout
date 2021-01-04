
#include "DetectorModule.h"
#include "ModuleCap.h"
#include "SimParms.h"

define_enum_strings(SensorLayout) = { "nosensors", "mono", "pt", "stereo" };
define_enum_strings(ZCorrelation) = { "samesegment", "multisegment" };
define_enum_strings(ReadoutType)  = { "strip", "pixel", "pt" };
define_enum_strings(ReadoutMode)  = { "binary", "cluster" };

//
// Constructor - specify unique id, geometry module defining shape & parse geometry config file using boost property tree & read-in module parameters
//
DetectorModule::DetectorModule(int id, GeometricModule* moduleGeom, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty) :
 BaseModule(moduleGeom),
 m_materialObject(MaterialObject::MODULE),
 minZ                     (string("minZ")             ),
 maxZ                     (string("maxZ")             ),
 minR                     (string("minR")             ),
 maxR                     (string("maxR")             ),
 planarMinZ               (string("planarMinZ")       ),
 planarMaxZ               (string("planarMaxZ")       ),
 planarMinR               (string("planarMinR")       ),
 planarMaxR               (string("planarMaxR")       ),
 minPhi                   (string("minPhi")           ),
 maxPhi                   (string("maxPhi")           ),
 maxSegments              (string("maxSegments")      ),
 minSegments              (string("minSegments")      ),
 totalSegments            (string("totalSegments")    ),
 maxChannels              (string("maxChannels")      ),
 minChannels              (string("minChannels")      ),
 totalChannels            (string("totalChannels")    ),
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
 supportPlateThickness    ("supportPlateThickness"    , parsedOnly(), 1),
 m_cntName(""),
 m_cntId(0)
{
  // Set the geometry config parameters
  this->myid(id);
  this->store(treeProperty);
  if (nodeProperty.count(id)>0) this->store(nodeProperty.at(id));
}

//
// Constructor - specify unique id, geometry module defining shape & parse geometry config file using boost property tree & read-in module parameters
//
DetectorModule::DetectorModule(int id, GeometricModule* moduleGeom, const PropertyTree& treeProperty) :
 BaseModule(moduleGeom),
 m_materialObject(MaterialObject::MODULE),
 minZ                     (string("minZ")             ),
 maxZ                     (string("maxZ")             ),
 minR                     (string("minR")             ),
 maxR                     (string("maxR")             ),
 planarMinZ               (string("planarMinZ")       ),
 planarMaxZ               (string("planarMaxZ")       ),
 planarMinR               (string("planarMinR")       ),
 planarMaxR               (string("planarMaxR")       ),
 minPhi                   (string("minPhi")           ),
 maxPhi                   (string("maxPhi")           ),
 maxSegments              (string("maxSegments")      ),
 minSegments              (string("minSegments")      ),
 totalSegments            (string("totalSegments")    ),
 maxChannels              (string("maxChannels")      ),
 minChannels              (string("minChannels")      ),
 totalChannels            (string("totalChannels")    ),
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
 supportPlateThickness    ("supportPlateThickness"    , parsedOnly(), 1),
 m_cntName(""),
 m_cntId(0)
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
  if (m_moduleCap!=nullptr) delete m_moduleCap;
  if (m_moduleGeom!=nullptr) delete m_moduleGeom;
}

//
// Build detector module with its sensors
//
void DetectorModule::build() {

  check();
  if (!m_moduleGeom->builtok()) {

    m_moduleGeom->store(propertyTree());
    m_moduleGeom->build();
  }
  if (numSensors() > 0) {

    for (int iSensor=1; iSensor<=numSensors(); iSensor++) {
      Sensor* s = GeometryFactory::make<Sensor>(iSensor, this, sensorNode, propertyTree());
      s->build();

      m_sensors.push_back(s);
      m_materialObject.sensorChannels[iSensor]=s->numChannels();
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

  maxZAllMat.setup([&]() {return maxget2(m_sensors.begin(), m_sensors.end(), &Sensor::maxZAllMat); });
  minZAllMat.setup([&]() {return minget2(m_sensors.begin(), m_sensors.end(), &Sensor::minZAllMat); });
  maxRAllMat.setup([&]() {return maxget2(m_sensors.begin(), m_sensors.end(), &Sensor::maxRAllMat); });
  minRAllMat.setup([&]() {return minget2(m_sensors.begin(), m_sensors.end(), &Sensor::minRAllMat); });

  planarMaxZ.setup([&]() { return basePoly().computeMaxZ(); });
  planarMinZ.setup([&]() { return basePoly().computeMinZ(); });
  planarMaxR.setup([&]() { return basePoly().computeMaxR(); });
  planarMinR.setup([&]() { return basePoly().computeMinR(); });

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

void DetectorModule::setModuleCap(ModuleCap* newCap)
{
  if (m_moduleCap!=nullptr) delete m_moduleCap;
  m_moduleCap = newCap ;
}

void DetectorModule::rotateToNegativeZSide() {
  side(-side());
  rotateY(M_PI);  // Rotation around FCC_Y of angle Pi
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

//
// Check if track hit the module -> if yes, return true with passed material, hit position vector & hitType (which module sensor(s) was/were hit)
//
bool DetectorModule::checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir, Material& hitMaterial, HitType& hitType, XYZVector& hitPos) const {

  // Initialize: hit was found, material, hitPos & relative hit path length wrt perpendicular passage
  hitMaterial.radiation   = 0.;
  hitMaterial.interaction = 0.;
  hitType                 = HitType::NONE;
  hitPos.SetX(0.);
  hitPos.SetY(0.);
  hitPos.SetZ(0.);

  bool hitFound = false;

  // Detector module consists of 1 sensor
  if (numSensors() == 1) {

    auto segm = innerSensor().checkHitSegment(trackOrig, trackDir);
    if (segm.second > -1) {

      hitPos  = segm.first;
      hitType = HitType::INNER; // The following line used to return HitType::BOTH. Changing to INNER in order to avoid double hit counting
      hitFound= true;
    }
  }
  // 2 sensors detector module
  else {

    auto inSegm  = innerSensor().checkHitSegment(trackOrig, trackDir);
    auto outSegm = outerSensor().checkHitSegment(trackOrig, trackDir);

    // Hit found in both sensors -> use inner sensor coordinates as a hit position
    if (inSegm.second > -1 && outSegm.second > -1) { 

      hitPos  = inSegm.first;
      hitType = ((zCorrelation() == SAMESEGMENT && (inSegm.second / (maxSegments()/minSegments()) == outSegm.second)) || zCorrelation() == MULTISEGMENT) ? HitType::STUB : HitType::BOTH;
      hitFound= true;
    }
    // Hit found in inner sensor only
    else if (inSegm.second > -1)  {

      hitPos  = inSegm.first;
      hitType = HitType::INNER;
      hitFound= true;
    }
    // Hit found in outer sensor only
    else if (outSegm.second > -1) {

      hitPos  = outSegm.first;
      hitType = HitType::OUTER;
      hitFound= true;
    }
  }

  // Update material
  if (hitFound) {

    double theta = trackDir.theta();
    hitMaterial.radiation   = getModuleCap().getRadiationLength();
    hitMaterial.interaction = getModuleCap().getInteractionLength();

    if (subdet() == BARREL) {

      hitMaterial.radiation   /= sin(theta + tiltAngle());
      hitMaterial.interaction /= sin(theta + tiltAngle());

    }
    else if (subdet() == ENDCAP) {

      hitMaterial.radiation   /= cos(theta + tiltAngle() - M_PI/2); // Endcap has tiltAngle = pi/2
      hitMaterial.interaction /= cos(theta + tiltAngle() - M_PI/2); // Endcap has tiltAngle = pi/2

    }
    else {
      logWARNING("DetectorModule::checkTrackHits -> incorrectly scaled material, unknown module type. Neither barrel or endcap");
    }
  }

  return hitFound;
};

BarrelModule::BarrelModule(int id, GeometricModule* moduleGeom, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty) :
 DetectorModule(id, moduleGeom, nodeProperty, treeProperty)
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
    m_moduleGeom->rotateY(M_PI/2);
    m_rAxis = normal();
    m_moduleGeom->tiltAngle(0.);
    m_moduleGeom->skewAngle(0.);
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

      Polygon3D<4> polygon = Polygon3D<4>(basePoly());
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

      Polygon3D<4> polygon = Polygon3D<4>(basePoly());
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
  m_moduleGeom->accept(v);
}

//
// GeometryVisitor pattern -> barrel module visitable (const. option)
//
void BarrelModule::accept(ConstGeometryVisitor& v) const {
  v.visit(*this);
  v.visit(*(const DetectorModule*)this);
  m_moduleGeom->accept(v);
}

//
// Constructor - parse geometry config file using boost property tree & read-in module parameters, specify unique id
//
EndcapModule::EndcapModule(int id, GeometricModule* moduleGeom, const PropertyTree& treeProperty) :
 DetectorModule(id, moduleGeom, treeProperty)
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
    m_moduleGeom->tiltAngle(M_PI/2.);
    m_moduleGeom->skewAngle(0.);
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

        Polygon3D<4> polygon = Polygon3D<4>(basePoly());
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

        Polygon3D<4> polygon = Polygon3D<4>(basePoly());
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
  m_moduleGeom->accept(v);
}

//
// GeometryVisitor pattern -> endcap module visitable (const. option)
//
void EndcapModule::accept(ConstGeometryVisitor& v) const {
  v.visit(*this);
  v.visit(*(const DetectorModule*)this);
  m_moduleGeom->accept(v);
}


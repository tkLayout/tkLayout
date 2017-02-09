#include "Sensor.h"
#include "DetectorModule.h"

define_enum_strings(SensorType) = { "pixel", "largepix", "strip" };

//
// Constructor - specify unique id, const pointer to parent module & parse geometry config file using boost property tree & read-in module parameters
//
Sensor::Sensor(int id, const DetectorModule* parent, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty) :
 minZ           (string("minZ")   ),
 maxZ           (string("maxZ")   ),
 minR           (string("minR")   ),
 maxR           (string("maxR")   ),
 numSegments    ("numSegments"    , parsedOnly()),
 numStripsAcross("numStripsAcross", parsedOnly()),
 numROCX        ("numROCX"        , parsedOnly()),
 numROCY        ("numROCY"        , parsedOnly()),
 sensorThickness("sensorThickness", parsedOnly(), 0.1),
 type           ("sensorType"     , parsedOnly(), SensorType::None),
 m_parent(parent)
{
  // Set unique id & the geometry config parameters
  this->myid(id);
  this->store(treeProperty);
  if (nodeProperty.count(id)>0) this->store(nodeProperty.at(id));
}

//
// Constructor - specify unique id, const pointer to parent module
//
Sensor::Sensor(int id, const DetectorModule* parent) :
 minZ           (string("minZ")   ),
 maxZ           (string("maxZ")   ),
 minR           (string("minR")   ),
 maxR           (string("maxR")   ),
 numSegments    ("numSegments"    , parsedOnly()),
 numStripsAcross("numStripsAcross", parsedOnly()),
 numROCX        ("numROCX"        , parsedOnly()),
 numROCY        ("numROCY"        , parsedOnly()),
 sensorThickness("sensorThickness", parsedOnly(), 0.1),
 type           ("sensorType"     , parsedOnly(), SensorType::None),
 m_parent(parent)
{
  // Set unique id
  this->myid(id);
}

//
// Build -> just check all set values
//
void Sensor::build() {
  try {
    check();
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }

  // Clear polygons for the start -> when asking for polygons a replication based on Detector module polygon will be done.
  clearPolys();

  cleanup();
  builtok(true);
}

//
// Setup: link lambda functions to various DetectorModule related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
//
void Sensor::setup() {

  minR.setup([&]()       { return MIN( CoordinateOperations::computeMinR(lowerEnvelopePoly()), CoordinateOperations::computeMinR(upperEnvelopePoly()) ); });
  maxR.setup([&]()       { return MAX( CoordinateOperations::computeMaxR(lowerEnvelopePoly()), CoordinateOperations::computeMaxR(upperEnvelopePoly()) ); });
  minRAllMat.setup([&]() { return MIN( CoordinateOperations::computeMinR(lowerEnvelopePoly(true)), CoordinateOperations::computeMinR(upperEnvelopePoly(true)) ); });
  maxRAllMat.setup([&]() { return MAX( CoordinateOperations::computeMaxR(lowerEnvelopePoly(true)), CoordinateOperations::computeMaxR(upperEnvelopePoly(true)) ); });
  minZ.setup([&]()       { return MIN( CoordinateOperations::computeMinZ(lowerEnvelopePoly()), CoordinateOperations::computeMinZ(upperEnvelopePoly()) ); });
  maxZ.setup([&]()       { return MAX( CoordinateOperations::computeMaxZ(lowerEnvelopePoly()), CoordinateOperations::computeMaxZ(upperEnvelopePoly()) ); });
  minZAllMat.setup([&]() { return MIN( CoordinateOperations::computeMinZ(lowerEnvelopePoly(true)), CoordinateOperations::computeMinZ(upperEnvelopePoly(true)) ); });
  maxZAllMat.setup([&]() { return MAX( CoordinateOperations::computeMaxZ(lowerEnvelopePoly(true)), CoordinateOperations::computeMaxZ(upperEnvelopePoly(true)) ); });
}

//
// Set reference to Detector module after cloning
//
void Sensor::parent(const DetectorModule* parent) { m_parent = parent; }

//
// Get standard offset wrt module average position -> +-dsDistance if dsDistance defined
//
double Sensor::normalOffset() const {
  return m_parent->numSensors() <= 1 ? 0. : (myid() == 1 ? -m_parent->dsDistance()/2. : m_parent->dsDistance()/2.);
}

//
// Build sensor geometrical representation based on detector module geometrical representation shifted by offset (i.e. by +-thickness/2. to get outer/inner envelope etc.)
//
Polygon3d<4>* Sensor::buildOwnPoly(double polyOffset) const {
  Polygon3d<4>* p = new Polygon3d<4>(m_parent->basePoly());
  p->translate(p->getNormal()*polyOffset);
  return p;
}

//
// Clear sensor polygons -> used when recalculating sensor/module position -> recreating geometrical implementation, i.e. polygon
//
void Sensor::clearPolys() { 

  if (m_hitPoly     !=nullptr) delete m_hitPoly;      m_hitPoly      = nullptr;
  if (m_lowerEnvPoly!=nullptr) delete m_lowerEnvPoly; m_lowerEnvPoly = nullptr;
  if (m_upperEnvPoly!=nullptr) delete m_upperEnvPoly; m_upperEnvPoly = nullptr;
}

//
// Get upper envelope of the sensor (taking into all material if required or just correct sensor Thickness and dsDistance of the module)
//
const Polygon3d<4>& Sensor::upperEnvelopePoly(bool applyAllMaterial) const {

  if (!applyAllMaterial) {
    if (m_upperEnvPoly==nullptr) m_upperEnvPoly = buildOwnPoly(normalOffset() + sensorThickness()/2.);
    return *m_upperEnvPoly;
  }
  else {
    if (m_upperEnvPolyAllMat==nullptr) m_upperEnvPolyAllMat = buildOwnPoly(m_parent->thicknessAllMat()/2.);
    return *m_upperEnvPolyAllMat;
  }
}

//
// Get lower envelope of the sensor (taking into all material if required or just correct sensor Thickness and dsDistance of the module)
//
const Polygon3d<4>& Sensor::lowerEnvelopePoly(bool applyAllMaterial) const {

  if (!applyAllMaterial) {
    if (m_lowerEnvPoly==nullptr) m_lowerEnvPoly = buildOwnPoly(normalOffset() - sensorThickness()/2.);
    return *m_lowerEnvPoly;
  }
  else {
    if (m_lowerEnvPolyAllMat==nullptr) m_lowerEnvPolyAllMat = buildOwnPoly(-m_parent->thicknessAllMat()/2.);
    return *m_lowerEnvPolyAllMat;
  }
}

//
// Get geometry properties
//
const Polygon3d<4>& Sensor::hitPoly() const {
  if (m_hitPoly == nullptr) m_hitPoly = buildOwnPoly(normalOffset());
  return *m_hitPoly;
}

//
// Check whether sensor hitted by a track (coming from trackOrig & shot at given direction)
//
std::pair<XYZVector, int> Sensor::checkHitSegment(const XYZVector& trackOrig, const XYZVector& trackDir) const {
  const Polygon3d<4>& poly = hitPoly();
  XYZVector p;
  if (poly.isLineIntersecting(trackOrig, trackDir, p)) {
    XYZVector v = p - poly.getVertex(0);
    double projL = v.Dot((poly.getVertex(1) - poly.getVertex(0)).Unit());
    return std::make_pair(p, projL / stripLength()); 
  } else return std::make_pair(p, -1);
}

double Sensor::minPitch()    const { return m_parent->minWidth() / (double)numStripsAcross(); }
double Sensor::maxPitch()    const { return m_parent->maxWidth() / (double)numStripsAcross(); }
double Sensor::pitch()       const { return m_parent->meanWidth() / (double)numStripsAcross(); }
double Sensor::stripLength() const { return m_parent->length() / numSegments(); }


#include "Sensor.h"
#include "DetectorModule.h"

//
// Constructor - specify unique id, const pointer to parent module & parse geometry config file using boost property tree & read-in module parameters
//
Sensor::Sensor(int id, const DetectorModule* parent, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty) :
 numSegments(    "numSegments"    , parsedOnly()),
 numStripsAcross("numStripsAcross", parsedOnly()),
 numROCX(        "numROCX"        , parsedOnly()),
 numROCY(        "numROCY"        , parsedOnly()),
 sensorThickness("sensorThickness", parsedOnly(), 0.1),
 type(           "sensorType"     , parsedOnly(), SensorType::None),
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
 numSegments(    "numSegments"    , parsedOnly()),
 numStripsAcross("numStripsAcross", parsedOnly()),
 numROCX(        "numROCX"        , parsedOnly()),
 numROCY(        "numROCY"        , parsedOnly()),
 sensorThickness("sensorThickness", parsedOnly(), 0.1),
 type(           "sensorType"     , parsedOnly(), SensorType::None),
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

  minR.setup([&]() { return CoordinateOperations::computeMinR(envelopePoly()); });
  maxR.setup([&]() { return CoordinateOperations::computeMaxR(envelopePoly()); });
  minZ.setup([&]() { return CoordinateOperations::computeMinZ(envelopePoly()); });
  maxZ.setup([&]() { return CoordinateOperations::computeMaxZ(envelopePoly()); });
}

//
// Set reference to Detector module after cloning
//
void Sensor::parent(const DetectorModule* parent) { m_parent = parent; }

double Sensor::normalOffset() const {
  return m_parent->numSensors() <= 1 ? 0. : (myid() == 1 ? -m_parent->dsDistance()/2. : m_parent->dsDistance()/2.);
}

Polygon3d<4>* Sensor::buildOwnPoly(double polyOffset) const {
  Polygon3d<4>* p = new Polygon3d<4>(m_parent->basePoly());
  p->translate(p->getNormal()*polyOffset);
  return p;
}

void Sensor::clearPolys() { 

  delete m_hitPoly; m_hitPoly = nullptr;
  delete m_envPoly; m_envPoly = nullptr;
}

const Polygon3d<4>& Sensor::hitPoly() const {
  if (m_hitPoly == nullptr) m_hitPoly = buildOwnPoly(normalOffset());
  return *m_hitPoly;
}

const Polygon3d<4>& Sensor::envelopePoly() const {
  if (m_envPoly == nullptr) {
    double envelopeOffset = normalOffset() > 1e-6 ? normalOffset() + sensorThickness()/2. : (normalOffset() < -1e-6 ? normalOffset() - sensorThickness()/2. : 0.);
    m_envPoly = buildOwnPoly(envelopeOffset);
  }
  return *m_envPoly;
}

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

define_enum_strings(SensorType) = { "pixel", "largepix", "strip" };

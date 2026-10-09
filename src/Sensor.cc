#include "Sensor.hh"
#include "DetectorModule.hh"

#include <memory>
#include <Math/Vector3Dfwd.h>

void Sensor::check() {
  PropertyObject::check();
  
  if (!numStripsAcross.state() && !pitchEstimate.state()) throw PathfulException("At least one between numStripsAcross and pitchEstimate must be specified");
  if (numStripsAcross.state() && pitchEstimate.state()) throw PathfulException("Only one between numStripsAcross and pitchEstimate can be specified");
  if (!numSegments.state() && !stripLengthEstimate.state()) throw PathfulException("At least one between numSegments and stripLengthEstimate must be specified");
  if (numSegments.state() && stripLengthEstimate.state()) throw PathfulException("Only one between numSegments and stripLengthEstimate can be specified");
}

const Polygon3d<4>& Sensor::hitPoly() const {
  if (hitPoly_) return *hitPoly_;

  const int numSensors = parent_->numSensors();

  const Polygon3d<4>& basePoly = parent_->basePoly();

  // Single-sensor module, no offset and/or polygon resizing required
  if (numSensors <= 1) {
    if (numSensors <= 0)
      logWARNING("Sensor::hitPoly Warning: DetectorModule::numSensors() <= 0");
    hitPoly_ = new Polygon3d<4>(basePoly);
  }

  // Special case for split-sensor modules
  else if (numSensors == 2 && parent_->sensorLayout() == SensorLayout::MONO) {
    // Account for the dead space between sensors when resizing the module polygon
    const double moduleLength = parent_->length();
    const double deadLength   = parent_->centralDeadAreaLength();
    const double sensorLength = moduleLength / numSensors - 0.5 * deadLength;
    const double scaleFactor  = sensorLength / moduleLength;

    // Module center and unit local Y
    const ROOT::Math::XYZVector moduleCenter = parent_->center();
    const ROOT::Math::XYZVector localY = (basePoly.getVertex(0) - basePoly.getVertex(3)).Unit();

    // Translate the sensor along the local Y
    double offsetY = (moduleLength - sensorLength) * 0.5;
    offsetY = (innerOuter() == SensorPosition::UPPER) ? offsetY : -offsetY;
    ROOT::Math::XYZVector sensorCenter = moduleCenter + offsetY * localY;

    // Resize the Y-length of the module polygon
    std::unique_ptr<Polygon3d<4>> newPoly = std::make_unique<Polygon3d<4>>();
    for (const ROOT::Math::XYZVector* vtx = basePoly.begin(); vtx != basePoly.end(); ++vtx) {
      // Distance vector from the module center to the vertex
      ROOT::Math::XYZVector d = *vtx - moduleCenter;

      // Projection of "d" along the local Y
      ROOT::Math::XYZVector projY = d.Dot(localY) * localY;

      // Start with original "d", subtract the original Y, and add the scaled Y
      ROOT::Math::XYZVector scaledVtx = d + projY * (scaleFactor - 1.0);

      // Update the vertex position and add it to the new polygon
      *newPoly << (sensorCenter + scaledVtx);
    }

    hitPoly_ = newPoly.release();
  }

  // General 3D module case
  else {
    // Calculate the translation along the local Z
    double offsetZ = 0.5 * parent_->dsDistance();
    offsetZ = (innerOuter() == SensorPosition::UPPER) ? offsetZ : -offsetZ;
    ROOT::Math::XYZVector extentZ = offsetZ * parent_->normal();

    // Apply the translation
    std::unique_ptr<Polygon3d<4>> newPoly = std::make_unique<Polygon3d<4>>(basePoly);
    newPoly->translate(extentZ);

    hitPoly_ = newPoly.release();
  }

  return *hitPoly_;
}

const Polygon3d<4>& Sensor::hitMidPoly() const {
  if (hitMidPoly_ == 0) hitMidPoly_ = CoordinateOperations::computeMidPolygon(hitPoly());
  return *hitMidPoly_;
}

const Polygon3d<8>& Sensor::envelopePoly() const {
  double envelopeOffset =  sensorThickness() / 2.;
  if (envelopePoly_ == 0) envelopePoly_ = CoordinateOperations::computeEnvelopePolygon<Polygon3d<4>, Polygon3d<8> >(hitPoly(), envelopeOffset); 
  return *envelopePoly_;
}

const Polygon3d<8>& Sensor::envelopeMidPoly() const {
  double envelopeOffset =  sensorThickness() / 2.;
  if (envelopeMidPoly_ == 0) envelopeMidPoly_ = CoordinateOperations::computeEnvelopePolygon<Polygon3d<4>, Polygon3d<8> >(hitMidPoly(), envelopeOffset); 
  return *envelopeMidPoly_;
}

void Sensor::clearPolys() { 
  delete hitPoly_; 
  hitPoly_ = 0;
  delete hitMidPoly_; 
  hitMidPoly_ = 0;
  delete envelopePoly_;
  envelopePoly_ = 0;
  delete envelopeMidPoly_;
  envelopeMidPoly_ = 0;
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

int Sensor::numStripsAcrossEstimate() const {
  if (numStripsAcross.state()) return numStripsAcross();
  else return floor(parent_->meanWidth() / pitchEstimate() + 0.5);
}
int Sensor::numSegmentsEstimate() const {
  if (numSegments.state()) return numSegments();
  else return floor(parent_->length() / stripLengthEstimate() + 0.5);
}
double Sensor::minPitch() const { return parent_->minWidth() / (double)numStripsAcrossEstimate(); }
double Sensor::maxPitch() const { return parent_->maxWidth() / (double)numStripsAcrossEstimate(); }
double Sensor::pitch() const { return parent_->meanWidth() / (double)numStripsAcrossEstimate(); }
double Sensor::stripLength() const { return parent_->length() / numSegmentsEstimate(); }

double Sensor::alveolaWidth() const { return parent_->meanWidth() / (double)numCrystalsX(); }
double Sensor::alveolaLength() const { return parent_->length() / (double)numCrystalsY(); }

define_enum_strings(SensorType) = { "pixel", "largepix", "strip" };

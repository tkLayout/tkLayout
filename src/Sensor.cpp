#include "Sensor.h"
#include "DetectorModule.h"

void Sensor::check() {
  PropertyObject::check();
  
  if (!numStripsAcross.state() && !pitchEstimate.state()) throw PathfulException("At least one between numStripsAcross and pitchEstimate must be specified");
  if (numStripsAcross.state() && pitchEstimate.state()) throw PathfulException("Only one between numStripsAcross and pitchEstimate can be specified");
  if (!numSegments.state() && !stripLengthEstimate.state()) throw PathfulException("At least one between numSegments and stripLengthEstimate must be specified");
  if (numSegments.state() && stripLengthEstimate.state()) throw PathfulException("Only one between numSegments and stripLengthEstimate can be specified");
}

double Sensor::normalOffset() const {
  return parent_->numSensors() <= 1 ? 0. : (myid() == 1 ? -parent_->dsDistance()/2. : parent_->dsDistance()/2.);
}

Polygon3d<4>* Sensor::buildOwnPoly(double polyOffset) const {
  Polygon3d<4>* p = new Polygon3d<4>(parent_->basePoly());
  p->translate(p->getNormal()*polyOffset);
  return p;
}

void Sensor::clearPolys() { 
  delete hitPoly_; 
  hitPoly_ = 0; 
  delete envPoly_;
  envPoly_ = 0;
}

const Polygon3d<4>& Sensor::hitPoly() const {
  if (hitPoly_ == 0) hitPoly_ = buildOwnPoly(normalOffset());
  return *hitPoly_;
}

const Polygon3d<4>& Sensor::envelopePoly() const {
  if (envPoly_ == 0) {
    double envelopeOffset = normalOffset() > 1e-6 ? normalOffset() + sensorThickness()/2. : (normalOffset() < -1e-6 ? normalOffset() - sensorThickness()/2. : 0.);
    envPoly_ = buildOwnPoly(envelopeOffset); 
  }
  return *envPoly_;
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

define_enum_strings(SensorType) = { "pixel", "largepix", "strip" };

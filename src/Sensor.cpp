#include "Sensor.h"
#include "DetectorModule.h"


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

double Sensor::minPitch() const { return parent_->minWidth() / (double)numStripsAcross(); }
double Sensor::maxPitch() const { return parent_->maxWidth() / (double)numStripsAcross(); }
double Sensor::pitch() const { return parent_->meanWidth() / 2. / (double)numStripsAcross(); }
double Sensor::stripLength() const { return parent_->length() / numSegments(); }

define_enum_strings(SensorType) = { "pixel", "largepix", "strip" };

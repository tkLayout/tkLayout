#include "Sensor.h"
#include "DetectorModule.h"

void Sensor::check() {
  PropertyObject::check();
  
  if (!numStripsAcross.state() && !pitchEstimate.state()) throw PathfulException("At least one between numStripsAcross and pitchEstimate must be specified");
  if (numStripsAcross.state() && pitchEstimate.state()) throw PathfulException("Only one between numStripsAcross and pitchEstimate can be specified");
  if (!numSegments.state() && !stripLengthEstimate.state()) throw PathfulException("At least one between numSegments and stripLengthEstimate must be specified");
  if (numSegments.state() && stripLengthEstimate.state()) throw PathfulException("Only one between numSegments and stripLengthEstimate can be specified");
}

double Sensor::sensorNormalOffset() const {
  double offset;
  if (parent_->numSensors() <= 1) offset = 0.;
  else {
    if (myid() == 1) offset = -parent_->dsDistance()/2.;
    else offset = parent_->dsDistance()/2.;
  }
  return offset;
}

Polygon3d<4>* Sensor::buildHitPoly(double sensorNormalOffset) const {
  Polygon3d<4>* p = new Polygon3d<4>(parent_->basePoly());
  p->translate(p->getNormal() * sensorNormalOffset);
  return p;
}

Polygon3d<8>* Sensor::buildEnvelopePoly(double sensorNormalOffset) const {
  double innerOffset =  sensorNormalOffset - sensorThickness() / 2.;
  Polygon3d<4>* innerPlane = new Polygon3d<4>(parent_->basePoly());
  innerPlane->translate(innerPlane->getNormal() * innerOffset);

  double outerOffset =  sensorNormalOffset + sensorThickness() / 2.;
  Polygon3d<4>* outerPlane = new Polygon3d<4>(parent_->basePoly());
  outerPlane->translate(outerPlane->getNormal() * outerOffset);

  Polygon3d<8>* envelopePoly = new Polygon3d<8>();
  for (int i = 0; i < innerPlane->getNumSides(); i++) *envelopePoly << innerPlane->getVertex(i);
  for (int i = 0; i < outerPlane->getNumSides(); i++) *envelopePoly << outerPlane->getVertex(i);
  // Would have been much nicer to do sth like : *envelopePoly << innerPlane->getVertices() , 
  // unfortunately this apparently does not work and only assign the first two vertices to *envelopePoly.
  // On a more general note, the entire AbstractPolygon class should be rewritten.
 
  return envelopePoly;
}

const Polygon3d<4>& Sensor::hitPoly() const {
  if (hitPoly_ == 0) hitPoly_ = buildHitPoly(sensorNormalOffset());
  return *hitPoly_;
}

const Polygon3d<8>& Sensor::envelopePoly() const {
  if (envPoly_ == 0) envPoly_ = buildEnvelopePoly(sensorNormalOffset()); 
  return *envPoly_;
}

void Sensor::clearPolys() { 
  delete hitPoly_; 
  hitPoly_ = 0; 
  delete envPoly_;
  envPoly_ = 0;
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

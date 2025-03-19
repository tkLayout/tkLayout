#include "Sensor.hh"
#include "DetectorModule.hh"

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

double Sensor::sensorXOffset() const {
  double offset;
  if (parent_->numSensors() <= 1) offset = 0.;
  else {
    if (myid() == 1) offset = -parent_->offsetForSensors();
    else offset = parent_->offsetForSensors();
  }
  return offset;
}


const Polygon3d<4>& Sensor::hitPoly() const {
  if (hitPoly_) return *hitPoly_;
  double offset = sensorNormalOffset();
  hitPoly_ = CoordinateOperations::computeTranslatedPolygon(parent_->basePoly(), offset);
  // Special case for split-sensor modules: sensor's polygon is reduced and shifted along local Y
  if (parent_->numSensors() > 1 && innerOuter() == SensorPosition::NO) {
    const XYZVector unitY(parent_->getLocalY());
    // Resizing the sensor polygon wrt the module polygon, accounting for the dead space between sensors
    double moduleLength = parent_->length();
    double deadLength = parent_->centralDeadAreaLength();
    double sensorLength = moduleLength / parent_->numSensors() - deadLength / 2.;
    Polygon3d<4> *poly = CoordinateOperations::computeResizedPolygon(*hitPoly_, unitY, sensorLength/moduleLength);
    // Shifting the sensor
    double offsetY = -moduleLength/2. + (myid() - 1) * (sensorLength + deadLength) + sensorLength / 2.;
    poly->translate(unitY*offsetY);
    hitPoly_ = poly;
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

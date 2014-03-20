#include "Sensor.h"

Polygon3d<4> Sensor::computeEnvelopePolygon() const {
  if (!hasPoly()) throw std::runtime_error("Tried to perform operation on sensor polygon before creating it");
  Polygon3d<4> temp(*poly_);
  double envelopeTranslation = sensorTranslation_ > 1e-6 ? sensorThickness()/2 : (sensorTranslation_ < -1e-6 ? -sensorThickness()/2 : 0);
  temp.translate(temp.getNormal()*envelopeTranslation);
  return temp; 
}

std::pair<XYZVector, int> Sensor::checkHitSegment(const XYZVector& trackOrig, const XYZVector& trackDir) const {
  XYZVector p;
  if (poly_->isLineIntersecting(trackOrig, trackDir, p)) {
    XYZVector v = p - poly_->getVertex(0);
    double projL = v.Dot((poly_->getVertex(1) - poly_->getVertex(0)).Unit());
    return std::make_pair(p, projL / stripLength()); 
  } else return std::make_pair(p, -1);
}

define_enum_strings(SensorType) = { "pixel", "largepix", "strip" };

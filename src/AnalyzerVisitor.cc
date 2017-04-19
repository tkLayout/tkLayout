#include "AnalyzerVisitor.hh"

void AnalyzerHelpers::drawModuleOnMap(const Module& m, double val, TH2D& map, TH2D& counter) {
  const Polygon3d<4>& poly = m.basePoly();
  XYZVector start = (poly.getVertex(0) + poly.getVertex(3))/2;
  XYZVector end = (poly.getVertex(1) + poly.getVertex(2))/2;
  XYZVector diff = end-start;
  XYZVector point;
  for (double l=0; l<=1; l+=0.1) {
    point = start + l * diff;
    map.Fill(point.Z(), point.Rho(), val);
    counter.Fill(point.Z(), point.Rho(), 1);
  }
}
void AnalyzerHelpers::drawModuleOnMap(const Module& m, double val, TH2D& map) {
  const Polygon3d<4>& poly = m.basePoly();
  XYZVector start = (poly.getVertex(0) + poly.getVertex(3))/2;
  XYZVector end = (poly.getVertex(1) + poly.getVertex(2))/2;
  XYZVector diff = end-start;
  XYZVector point;
  for (double l=0; l<=1; l+=0.1) {
    point = start + l * diff;
    map.Fill(point.Z(), point.Rho(), val);
  }
}

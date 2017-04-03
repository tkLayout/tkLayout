#ifndef VISITOR_H
#define VISITOR_H

class Barrel;
class BarrelModule;
class BeamPipe;
class Detector;
class DetectorModule;
class Disk;
class Endcap;
class EndcapModule;
class GeometricModule;
class Layer;
class RectangularModule;
class Ring;
class TiltedRing;
class RodPair;
class WedgeModule;
class SimParms;
class SupportStructure;
class Tracker;

/*
 * @classGeometryVisitor
 * Visitor design pattern to separate geometry analysis algorithms from geometry implementation - non const. variant
 */
class GeometryVisitor { 
public:
  virtual ~GeometryVisitor() {}
  virtual void visit(Detector&) {}
  virtual void visit(BeamPipe&) {}
  virtual void visit(Tracker&) {}
  virtual void visit(Barrel&) {}
  virtual void visit(Endcap&) {}
  virtual void visit(SupportStructure&) {}
  virtual void visit(Layer&) {}
  virtual void visit(Disk&) {}
  virtual void visit(Ring&) {}
  virtual void visit(TiltedRing&) {}
  virtual void visit(RodPair&) {}
  virtual void visit(BarrelModule&) {}
  virtual void visit(EndcapModule&) {}
  virtual void visit(DetectorModule&) {}
  virtual void visit(RectangularModule&) {}
  virtual void visit(WedgeModule&) {}
  virtual void visit(GeometricModule&) {}
  virtual void visit(SimParms&) {}
};

/*
 * @classConstGeometryVisitor
 * Visitor design pattern to separate geometry analysis algorithms from geometry implementation - const. variant
 */
class ConstGeometryVisitor {
public:
  virtual ~ConstGeometryVisitor() {}
  virtual void visit(const Detector&) {}
  virtual void visit(const BeamPipe&) {}
  virtual void visit(const Tracker&) {}
  virtual void visit(const Barrel&) {}
  virtual void visit(const Endcap&) {}
  virtual void visit(const SupportStructure&) {}
  virtual void visit(const Layer&) {}
  virtual void visit(const Disk&) {}
  virtual void visit(const Ring&) {}
  virtual void visit(const TiltedRing&) {}
  virtual void visit(const RodPair&) {}
  virtual void visit(const BarrelModule&) {}
  virtual void visit(const EndcapModule&) {}
  virtual void visit(const DetectorModule&) {}
  virtual void visit(const RectangularModule&) {}
  virtual void visit(const WedgeModule&) {}
  virtual void visit(const GeometricModule&) {}
  virtual void visit(const SimParms&) {}
};

#endif

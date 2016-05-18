#ifndef VISITOR_H
#define VISITOR_H


class BeamPipe;
class Tracker;
class Barrel;
class Endcap;
class Layer;
class Disk;
class Ring;
class RodPair;
class BarrelModule;
class EndcapModule;
class DetectorModule;
class RectangularModule;
class WedgeModule;
class GeometricModule;
class SimParms;

class GeometryVisitor { 
public:
  virtual ~GeometryVisitor() {}
  virtual void visit(BeamPipe&) {}
  virtual void visit(Tracker&) {}
  virtual void visit(Barrel&) {}
  virtual void visit(Endcap&) {}
  virtual void visit(Layer&) {}
  virtual void visit(Disk&) {}
  virtual void visit(Ring&) {}
  virtual void visit(RodPair&) {}
  virtual void visit(BarrelModule&) {}
  virtual void visit(EndcapModule&) {}
  virtual void visit(DetectorModule&) {}
  virtual void visit(RectangularModule&) {}
  virtual void visit(WedgeModule&) {}
  virtual void visit(GeometricModule&) {}
  virtual void visit(SimParms&) {}
};

class ConstGeometryVisitor {
public:
  virtual ~ConstGeometryVisitor() {}
  virtual void visit(const BeamPipe&) {}
  virtual void visit(const Tracker&) {}
  virtual void visit(const Barrel&) {}
  virtual void visit(const Endcap&) {}
  virtual void visit(const Layer&) {}
  virtual void visit(const Disk&) {}
  virtual void visit(const Ring&) {}
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

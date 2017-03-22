#ifndef VISITOR_H
#define VISITOR_H


class Detector;
class Tracker;
class Barrel;
class Endcap;
class Layer;
class Disk;
class Ring;
class TiltedRing;
class RodPair;
class BarrelModule;
class EndcapModule;
class DetectorModule;
class RectangularModule;
class WedgeModule;
class GeometricModule;
class Sensor;
class SimParms;

class GeometryVisitor { 
public:
  virtual void visit(Detector&) {}
  virtual void visit(Tracker&) {}
  virtual void visit(Barrel&) {}
  virtual void visit(Endcap&) {}
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

class ConstGeometryVisitor {
public:
  virtual void visit(const Detector&) {}
  virtual void visit(const Tracker&) {}
  virtual void visit(const Barrel&) {}
  virtual void visit(const Endcap&) {}
  virtual void visit(const Layer&) {}
  virtual void visit(const Disk&) {}
  virtual void visit(const TiltedRing&) {}
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

class SensorGeometryVisitor { 
public:
  virtual void visit(Detector&) {}
  virtual void visit(Tracker&) {}
  virtual void visit(Barrel&) {}
  virtual void visit(Endcap&) {}
  virtual void visit(Layer&) {}
  virtual void visit(Disk&) {}
  virtual void visit(Ring&) {}
  virtual void visit(TiltedRing&) {}
  virtual void visit(RodPair&) {}
  virtual void visit(BarrelModule&) {}
  virtual void visit(EndcapModule&) {}
  virtual void visit(DetectorModule&) {}
  // virtual void visit(RectangularModule&) {}
  // virtual void visit(WedgeModule&) {}
  //virtual void visit(GeometricModule&) {}
  virtual void visit(Sensor&) {}
  //virtual void visit(SimParms&) {}
};

#endif

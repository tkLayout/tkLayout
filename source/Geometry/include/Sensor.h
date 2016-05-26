#ifndef SENSOR_H
#define SENSOR_H

#include <string>
#include <exception>

#include "global_funcs.h"
#include "Polygon3d.h"
#include "Property.h"
#include "CoordinateOperations.h"

class DetectorModule;

enum class SensorType { Pixel, Largepix, Strip, None };

/*
 * @class Sensor
 * @brief Sensor object is part of detector module, which might consist of several sensors
 * @details Sensor geometrical representation is given by polygon (cloned from Detector module geometrical representation, i.e. basePoly).
 * Anytime the polygon is required, it's deleted first and again replicated based on Detector module polygon.
 */
class Sensor : public PropertyObject, public Buildable, public Identifiable<int> {

 public:

  //! Constructor - specify unique id, const pointer to parent module & parse geometry config file using boost property tree & read-in module parameters
  Sensor(int id, const DetectorModule* parent, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty);

  //! Constructor - specify unique id, const pointer to parent module
  Sensor(int id, const DetectorModule* parent);

  //! Destructor
  virtual ~Sensor() {};

  //! Build -> just check all set values
  void build();

  //! Setup: link lambda functions to various DetectorModule related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Set reference to Detector module after cloning
  void parent(const DetectorModule* parent);

  //! Clear sensor polygons
   void clearPolys();

  // Get geometry properties
  const Polygon3d<4>& hitPoly() const;
  const Polygon3d<4>& envelopePoly() const;

  int    numChannels() const { return numStripsAcross() * numSegments(); }
  double minPitch()    const;
  double maxPitch()    const;
  double pitch()       const;
  double stripLength() const;

  int numROCRows() const { return numStripsAcross() / numROCX(); } 
  int numROCCols() const { return numSegments() / numROCY(); }
  int totalROCs()  const { return numROCX() * numROCY(); }

  double normalOffset() const;

  std::pair<XYZVector, int> checkHitSegment(const XYZVector& trackOrig, const XYZVector& trackDir) const;

  ReadonlyProperty<int   , NoDefault         > numSegments;
  ReadonlyProperty<int   , NoDefault         > numStripsAcross;
  ReadonlyProperty<int   , NoDefault         > numROCX, numROCY;
  ReadonlyProperty<double, Default           > sensorThickness;
  ReadonlyProperty<double, UncachedComputable> minR, maxR; // CUIDADO min/maxR don't take into account the sensor thickness!
  ReadonlyProperty<double, UncachedComputable> minZ, maxZ; // ditto for min/maxZ

  ReadonlyProperty<SensorType, Default> type;

 private:

  const DetectorModule* m_parent; //!< Const pointer to parent detector module

  //! Build sensor geometrical representation based on detector module geometrical representation
  Polygon3d<4>* buildOwnPoly(double polyOffset) const;

  mutable const Polygon3d<4>* m_hitPoly = nullptr;
  mutable const Polygon3d<4>* m_envPoly = nullptr; //! Sensor geometetrical representation

};

#endif

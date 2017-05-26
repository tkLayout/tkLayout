#ifndef SENSOR_H
#define SENSOR_H

#include <exception>
#include <string>

#include "GeometryFactory.h"
#include "Polygon3D.h"
#include "Property.h"

// Forward declaration
class DetectorModule;
class Sensor;

// Typedefs
typedef PtrVector<Sensor> Sensors;

enum SensorPosition { NO, LOWER, UPPER };
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

  //! Check whether sensor hitted by a track (coming from trackOrig & shot at given direction)
  std::pair<XYZVector, int> checkHitSegment(const XYZVector& trackOrig, const XYZVector& trackDir) const;

  //! Set reference to Detector module after cloning
  void parent(const DetectorModule* parent);

  // Get geometry properties
  const Polygon3D<4>& hitPoly() const;

  //! Get upper envelope of the sensor (taking into all material if required or just correct sensor Thickness and dsDistance of the module) -> if taking min/max take min/max(lower, upper)
  const Polygon3D<4>& upperEnvelopePoly(bool applyAllMaterial=false) const;

  //! Get lower envelope of the sensor (taking into all material if required or just correct sensor Thickness and dsDistance of the module) -> if taking min/max take min/max(lower, upper)
  const Polygon3D<4>& lowerEnvelopePoly(bool applyAllMaterial=false) const;

  //! Get standard offset wrt module average position -> +-dsDistance if dsDistance defined
  double normalOffset() const;

  //! Clear sensor polygons -> used when recalculating sensor/module position -> recreating geometrical implementation, i.e. polygon
   void clearPolys();

  int    numChannels() const { return numStripsAcross() * numSegments(); }
  double minPitch()    const;
  double maxPitch()    const;
  double pitch()       const;
  double stripLength() const;

  int numROCRows() const { return numStripsAcross() / numROCX(); } 
  int numROCCols() const { return numSegments() / numROCY(); }
  int totalROCs()  const { return numROCX() * numROCY(); }

  SensorPosition innerOuter(SensorPosition pos) { m_innerOuter = pos; return m_innerOuter; }
  SensorPosition innerOuter() const             { return m_innerOuter; }

  ReadonlyProperty<int   , NoDefault         > numSegments;      //TODO: Number of channels in RPhi -> Rename to number of readout channels in correct direction, is it R-Phi or Z?
  ReadonlyProperty<int   , NoDefault         > numStripsAcross;  //TODO: Number of channels in Z -> Rename to number of readout channels in correct direction, is it R-Phi or Z?
  ReadonlyProperty<int   , NoDefault         > numROCX, numROCY; //!< Number of read-out chips in given direction
  ReadonlyProperty<double, Default           > sensorThickness;  //!< Sensor thickness
  ReadonlyProperty<double, UncachedComputable> minR;             //!< Minimum sensor radius
  ReadonlyProperty<double, UncachedComputable> maxR;             //!< Maximum sensor radius
  ReadonlyProperty<double, UncachedComputable> minRAllMat;       //!< Minimum sensor radius taking into account all module material structures
  ReadonlyProperty<double, UncachedComputable> maxRAllMat;       //!< Maximum sensor radius taking into account all module material structures
  ReadonlyProperty<double, UncachedComputable> minZ;             //!< Minimum sensor Z position
  ReadonlyProperty<double, UncachedComputable> maxZ;             //!< Maximum sensor Z position
  ReadonlyProperty<double, UncachedComputable> minZAllMat;       //!< Minimum sensor Z position taking into account all module material structures
  ReadonlyProperty<double, UncachedComputable> maxZAllMat;       //!< Maximum sensor Z position taking into account all module material structures

  ReadonlyProperty<SensorType, Default> type;                    //!< Default sensor type: pixel, strip, ...

 private:

  const DetectorModule* m_parent; //!< Const pointer to parent detector module

  //! Build sensor geometrical representation based on detector module geometrical representation shifted by offset (i.e. by +-thickness/2. to get outer/inner envelope etc.)
  Polygon3D<4>* buildOwnPoly(double polyOffset) const;
  SensorPosition m_innerOuter = SensorPosition::NO;

  mutable const Polygon3D<4>* m_hitPoly            = nullptr;
  mutable const Polygon3D<4>* m_lowerEnvPoly       = nullptr; //! Lower envelope of sensor geometrical representation -> if taking min/max take min/max(lower, upper)
  mutable const Polygon3D<4>* m_upperEnvPoly       = nullptr; //! Upper envelope of sensor geometrical representation -> if taking min/max take min/max(lower, upper)
  mutable const Polygon3D<4>* m_lowerEnvPolyAllMat = nullptr; //! Lower envelope of sensor full material representation -> if taking min/max take min/max(lower, upper)
  mutable const Polygon3D<4>* m_upperEnvPolyAllMat = nullptr; //! Upper envelope of sensor full material representation -> if taking min/max take min/max(lower, upper)

};

#endif

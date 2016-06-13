/**
 * @file Hit.h
 * @brief This header file defines the hit and track classes used for internal analysis
 */
#ifndef INCLUDE_HIT_H_
#define INCLUDE_HIT_H_

#include <cmath>
#include <vector>

#include "DetectorModule.h"
#include "MaterialProperties.h"

// Forward declaration
class DetectorModule;
class Track;

#undef HIT_DEBUG
#undef HIT_DEBUG_RZ

enum class HitOrientation : short { Undefined, Horizontal, Vertical};  // Hit object orientation
enum class HitKind        : short { Undefined, Active, Inactive };     // Hit object type

/**
 * @class Hit
 * @brief The Hit class is used when analysing a tracker layout to record information about a volume that was hit by a test track.
 *
 * It is used for both active and inactive surface hits. In case of an inactive surface, the pointer to the hit module will be <i>NULL</i>
 * since those volumes only matter with respect to the radiation and interaction lengths they add to the total in the error calculations.
 * All the other information is available to both categories. For convenience, the scaled radiation and interaction lengths are stored in
 * here as well to avoid additional computation and callbacks to the material property objects.
 */
// TODO: Remove pointer to track and find another method without any pointers involved!
class Hit {

public:

  //! Constructor
  Hit();

  //
  Hit(const Hit& h);

  //! Constructor for a hit with no module at a given distance from the origin
  Hit(double myDistance);

  //! Constructor for a hit on a given module at a given distance from the origin
  Hit(double myDistance, const DetectorModule* myModule, HitType activeHitType);

  //! Destructor
  ~Hit();

  //! Given two hits, compare the distance to the z-axis.
  static bool sortSmallerR(Hit* h1, Hit* h2);

  //! Update hit radius (radius in 2D - xy plane), based on new track parameters
  void updateRadius() {m_radius = m_distance * sin(getTrackTheta());};

  // Setter methods
  void setHitModule(const DetectorModule* myModule);
  void setDistance(double newDistance)              { if (newDistance>0) m_distance = newDistance; updateRadius(); };
  void setOrientation(HitOrientation newOrientation){ m_orientation = newOrientation; };
  void setObjectKind(HitKind newObjectKind)         { m_objectKind = newObjectKind;};
  void setTrack(const Track* newTrack)              { m_track = newTrack; updateRadius();};
  void setCorrectedMaterial(RILength newMaterial)   { m_correctedMaterial = newMaterial;};
  void setPixel(bool isPixel)                       { m_isPixel = isPixel;}
  void setTrigger(bool isTrigger)                   { m_isTrigger = isTrigger;}
  void setResolutionRphi(double newRes)             { m_resolutionRphi = newRes; } // Only used for virtual hits on non-modules
  void setResolutionY(double newRes)                { m_resolutionY = newRes; } // Only used for virtual hits on non-modules
  bool setIP(bool newIP)                            { return m_isIP = newIP; }
  void setActiveHitType(HitType activeHitType)      { m_activeHitType = activeHitType; }

  // Getter methods
  const DetectorModule* getHitModule() { return m_hitModule; };

  double         getDistance() const               { return m_distance;};
  double         getRadius() const                 { return m_radius;};
  HitOrientation getOrientation() const            { return m_orientation;};
  HitKind        getObjectKind() const             { return m_objectKind;};
  RILength       getCorrectedMaterial();
  double         getResolutionRphi(double trackR);
  double         getResolutionZ(double trackR);
  double         getD();
  HitType        getActiveHitType() const          { return m_activeHitType; } // NONE, INNER, OUTER, BOTH or STUB -- only meaningful for hits on active elements

  bool isPixel() const   { return m_isPixel; };
  bool isTrigger() const { return m_isTrigger; };
  bool isIP() const      { return m_isIP; };
  bool isSquareEndcap();
  bool isStub() const;

protected:
  
  double         m_distance;      //!< Distance of hit from origin in 3D
  double         m_radius;        //!< Distance of hit from origin in the x/y plane
  HitOrientation m_orientation;   //!< Orientation of the surface
  HitKind        m_objectKind;    //!< Kind of hit object: Active, Inactive, ...
  HitType        m_activeHitType; //!< Hit coming from inner, outer, stub, ... module
  
  const DetectorModule* m_hitModule;  //!< Const pointer to the hit module
  const Track*          m_track;      //!< Const pointer to the track, into which the hit was assigned
  
  RILength m_correctedMaterial;
  
  bool m_isPixel;   //!< Hit coming from the pixel module?
  bool m_isTrigger; //!< Hit comint from the trigger module?
  bool m_isIP;      //!< An artificial hit, simulating IP constraint?

private:
  
  double   getTrackTheta();

  double m_resolutionRphi; // Only used for virtual hits on non-modules
  double m_resolutionY;    // Only used for virtual hits on non-modules

}; // Class

#endif /* INCLUDE_HIT_H_ */

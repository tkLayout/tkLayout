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
class Hit;
class InactiveElement;
class Track;

#undef HIT_DEBUG
#undef HIT_DEBUG_RZ

// Typedefs
typedef std::unique_ptr<Hit> HitPtr;
typedef std::vector<HitPtr>  HitCollection;

enum class HitActivity    : short { Undefined, Active, Inactive };               //!< Hit defined as pure material (inactive) or measurement point (active)
enum class HitPassiveType : short { Undefined, BeamPipe, IP, Service, Support }; //!< Hit defined as IP or pure material hit coming from: support, service, etc.

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

  //! Copy constructor
  Hit(const Hit& h);

  //! Constructor for a hit on an inactive surface at [rPos, zPos] (cylindrical position) from the origin
  Hit(double rPos, double zPos, const InactiveElement* myPassiveElem, HitPassiveType passiveHitType);

  //! Constructor for a hit on a given module at [rPos, zPos] (cylindrical position) from the origin
  Hit(double rPos, double zPos, const DetectorModule* myModule, HitType activeHitType);

  //! Destructor
  ~Hit();

  //! Given two hits, compare the distance to the z-axis based on smaller R
  static bool sortSmallerR(const HitPtr& h1, const HitPtr& h2);

  //! Given two hits, compare the distance to the z-axis based on higher R
  static bool sortHigherR(const HitPtr& h1, const HitPtr& h2);

  // Setter methods
  void setTrack(const Track* newTrack)            { m_track = newTrack;};

  void setAsActive()                              { m_activity = HitActivity::Active;};
  void setAsPassive()                             { m_activity = HitActivity::Inactive;};
  void setAsPixel()                               { m_isPixel    = true;}
  void setActiveHitType(HitType activeHitType)    { m_activeHitType = activeHitType; }
  void setCorrectedMaterial(RILength newMaterial) { m_correctedMaterial = newMaterial;};

  void setTrigger(bool isTrigger)                 { m_isTrigger = isTrigger;}
  void setResolutionRphi(double newRes)           { m_resolutionRPhi = newRes; } // Only used for virtual hits on non-modules
  void setResolutionZ(double newRes)              { m_resolutionZ = newRes; }    // Only used for virtual hits on non-modules
  void setResolutionY(double newRes)              { setResolutionZ(newRes); } // Used for compatibility only -> use setResolutionZ(double newRes) instead

  void setDetName(std::string detName)            { m_detName = detName; }
  void setLayerID(int layerID)                    { m_layerID = layerID; }
  void setDiscID(int discID)                      { m_discID  = discID; }

  // Getter methods
  const DetectorModule*  getHitModule() const         { return m_hitModule; };
  const InactiveElement* getHitPassiveElement() const { return m_hitPassiveElem; }

  double   getDistance() const         { return m_distance;};
  double   getRPos() const             { return m_rPos;};
  double   getZPos() const             { return m_zPos;};
  double   getTilt() const             { if (this->isMeasurable()) return m_hitModule->tiltAngle(); else return 0; };
  bool     isActive() const            { if (m_activity==HitActivity::Active)   return true; else return false;};
  bool     isPassive() const           { if (m_activity==HitActivity::Inactive) return true; else return false;};
  bool     isActivityUndefined() const { if (m_activity==HitActivity::Undefined) return true; else return false;};
  HitType  getActiveHitType() const    { return m_activeHitType; } // NONE, INNER, OUTER, BOTH or STUB -- only meaningful for hits on active elements
  RILength getCorrectedMaterial();
  double   getResolutionRphi(double trackRadius);
  double   getResolutionZ(double trackRadius);
  double   getD();

  std::string getDetName()       const { return m_detName; };
  int         getLayerOrDiscID() const { if(this->isBarrel()) return m_layerID; else if(this->isEndcap()) return m_discID; else return -1;}; //!< Return positive number for layer (barrel hit) or disc (end-cap hit), -1 for beam-pipe or IP

  bool     isBeamPipe() const  { if (m_passiveHitType==HitPassiveType::BeamPipe) return true; else return false; };
  bool     isService() const   { if (m_passiveHitType==HitPassiveType::Service) return true; else return false; };
  bool     isSupport() const   { if (m_passiveHitType==HitPassiveType::Support) return true; else return false; };
  bool     isIP() const        { if (m_passiveHitType==HitPassiveType::IP) return true; else return false; };
  bool     isPixel() const     { return m_isPixel; };
  bool     isBarrel() const    { if (m_hitModule && (m_hitModule->subdet()==BARREL)) return true; else return false;};
  bool     isEndcap() const    { if (m_hitModule && (m_hitModule->subdet()==ENDCAP)) return true; else return false;};
  bool     isMeasurable() const{ if (m_hitModule!=nullptr) return true; else return false;};
  bool     isTrigger() const   { return m_isTrigger; };


  bool     isSquareEndcap();
  bool     isStub() const;

protected:
  
  //! Default constructor
  Hit();

  //! Set pointer to hit module in the constructor
  void setHitModule(const DetectorModule* myModule);

  //! Set pointer to inactive element in the constructor
  void setHitPassiveElement(const InactiveElement* myPassiveElem);

  double         m_distance;      //!< Distance of hit from origin in 3D = sqrt(rPos*rPos + zPos*zPos)
  double         m_rPos;          //!< Distance of hit from origin in the x/y plane (cylindrical coordinates -> r)
  double         m_zPos;          //!< Distance of hit from origin in z (cylindrical coordinates -> z)
  HitActivity    m_activity;      //!< Hit defined as pure material (inactive) or measurement point (active)
  HitType        m_activeHitType; //!< Hit coming from inner, outer, stub, ... module
  HitPassiveType m_passiveHitType;//!< Hit coming from which passive part: beam-pipe, service, support etc.
  
  const DetectorModule*  m_hitModule;     //!< Const pointer to the hit module
  const InactiveElement* m_hitPassiveElem;//!< Const pointer to the hit inactive element
  const Track*           m_track;         //!< Const pointer to the track, into which the hit was assigned
  
  RILength m_correctedMaterial; //!< Material in the way of particle shot at m_track direction, i.e. theta, module tilt angles corrected
  
  bool m_isTrigger; //!< Hit comint from the trigger module?
  bool m_isPixel;   //!< Hit coming from the pixel module?

  std::string m_detName; //!< Detector name, in which the hit has been measured
  int         m_layerID; //!< Corresponding layer ID
  int         m_discID;  //!< Corresponding disc ID

private:
  
  double getTrackTheta();

  double m_resolutionRPhi; // Only used for virtual hits on non-modules
  double m_resolutionZ;    // Only used for virtual hits on non-modules

}; // Class

#endif /* INCLUDE_HIT_H_ */

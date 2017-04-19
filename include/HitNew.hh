/**
 * @file HitNew.h
 * @brief This header file defines the hit and track classes used for internal analysis
 */
#ifndef INCLUDE_HITNEW_H_
#define INCLUDE_HITNEW_H_

#include <cmath>
#include <vector>

#include "DetectorModule.hh"
#include "MaterialProperties.hh"

// Forward declaration
class DetectorModule;
class HitNew;
class TrackNew;

#undef HIT_DEBUG
#undef HIT_DEBUG_RZ

// Typedefs
typedef std::unique_ptr<HitNew> HitNewPtr;
typedef std::vector<HitNewPtr>  HitNewCollection;

enum class HitNewActivity    : short { Undefined, Active, Inactive };     // HitNew defined as pure material (inactive) or measurement point (active)

/**
 * @class HitNew
 * @brief The HitNew class is used when analysing a tracker layout to record information about a volume that was hit by a test track.
 *
 * It is used for both active and inactive surface hits. In case of an inactive surface, the pointer to the hit module will be <i>NULL</i>
 * since those volumes only matter with respect to the radiation and interaction lengths they add to the total in the error calculations.
 * All the other information is available to both categories. For convenience, the scaled radiation and interaction lengths are stored in
 * here as well to avoid additional computation and callbacks to the material property objects.
 */
// TODO: Remove pointer to track and find another method without any pointers involved!
class HitNew {

public:

  //! Copy constructor
  HitNew(const HitNew& h);

  //! Constructor for a hit with no module (for passive materials) at [rPos, zPos] (cylindrical position) from the origin
  HitNew(double rPos, double zPos);

  //! Constructor for a hit on a given module at [rPos, zPos] (cylindrical position) from the origin
  HitNew(double rPos, double zPos, const DetectorModule* myModule, HitType activeHitType);

  //! Destructor
  ~HitNew();

  //! Given two hits, compare the distance to the z-axis based on smaller R
  static bool sortSmallerR(const HitNewPtr& h1, const HitNewPtr& h2);

  //! Given two hits, compare the distance to the z-axis based on higher R
  static bool sortHigherR(const HitNewPtr& h1, const HitNewPtr& h2);

  // Setter methods
  void setHitModule(const DetectorModule* myModule);
  void setTrack(const TrackNew* newTrack)           { m_track = newTrack;};

  void setAsActive()                                { m_activity = HitNewActivity::Active;};
  void setAsPassive()                               { m_activity = HitNewActivity::Inactive;};
  void setAsPixel()                                 { m_isPixel    = true;}
  void setActiveHitType(HitType activeHitType)      { m_activeHitType = activeHitType; }
  void setCorrectedMaterial(RILength newMaterial)   { m_correctedMaterial = newMaterial;};

  void setBeamPipe(bool isBeamPipe)                 { m_detName = "BeamPipe"; m_isBeamPipe = isBeamPipe;}
  void setIP(bool newIP)                            { m_detName = "IP";       m_isIP = newIP; }
  void setTrigger(bool isTrigger)                   { m_isTrigger = isTrigger;}
  void setResolutionRphi(double newRes)             { m_resolutionRPhi = newRes; } // Only used for virtual hits on non-modules
  void setResolutionZ(double newRes)                { m_resolutionZ = newRes; }    // Only used for virtual hits on non-modules
  void setResolutionY(double newRes)                { setResolutionZ(newRes); }    // Used for compatibility only -> use setResolutionZ(double newRes) instead

  // Getter methods
  const DetectorModule* getHitModule() const { return m_hitModule; };

  double   getDistance() const         { return m_distance;};
  double   getRPos() const             { return m_rPos;};
  double   getZPos() const             { return m_zPos;};
  double   getTilt() const             { if (this->isMeasurable()) return m_hitModule->tiltAngle(); else return 0; };
  bool     isActive() const            { if (m_activity==HitNewActivity::Active)   return true; else return false;};
  bool     isPassive() const           { if (m_activity==HitNewActivity::Inactive) return true; else return false;};
  bool     isActivityUndefined() const { if (m_activity==HitNewActivity::Undefined) return true; else return false;};
  HitType  getActiveHitType() const    { return m_activeHitType; } // NONE, INNER, OUTER, BOTH or STUB -- only meaningful for hits on active elements
  RILength getCorrectedMaterial();
  double   getResolutionRphi(double trackRadius);
  double   getResolutionZ(double trackRadius);
  double   getD();

  std::string getDetName()       const { return m_detName; };
  int         getLayerOrDiscID() const;

  bool     isBeamPipe() const  { return m_isBeamPipe; };
  bool     isIP() const        { return m_isIP; };
  bool     isPixel() const     { return m_isPixel; };
  bool     isBarrel() const    { if (m_hitModule && (m_hitModule->subdet()==BARREL)) return true; else return false;};
  bool     isEndcap() const    { if (m_hitModule && (m_hitModule->subdet()==ENDCAP)) return true; else return false;};
  bool     isMeasurable() const{ if (m_hitModule!=nullptr) return true; else return false;};
  bool     isTrigger() const   { return m_isTrigger; };


  bool     isSquareEndcap();
  bool     isStub() const;

protected:
  
  //! Default constructor
  HitNew();

  double         m_distance;      //!< Distance of hit from origin in 3D = sqrt(rPos*rPos + zPos*zPos)
  double         m_rPos;          //!< Distance of hit from origin in the x/y plane (cylindrical coordinates -> r)
  double         m_zPos;          //!< Distance of hit from origin in z (cylindrical coordinates -> z)
  HitNewActivity m_activity;      //!< HitNew defined as pure material (inactive) or measurement point (active)
  HitType        m_activeHitType; //!< HitNew coming from inner, outer, stub, ... module
  
  const DetectorModule* m_hitModule;  //!< Const pointer to the hit module
  const TrackNew*       m_track;      //!< Const pointer to the track, into which the hit was assigned
  
  RILength m_correctedMaterial; //!< Material in the way of particle shot at m_track direction, i.e. theta, module tilt angles corrected
  
  bool m_isBeamPipe;//!< An artificial hit, simulating colision with beam-pipe?
  bool m_isIP;      //!< An artificial hit, simulating IP constraint?
  bool m_isTrigger; //!< Hit coming from the trigger module?
  bool m_isPixel;   //!< Hit coming from the pixel module

  std::string m_detName; //!< Detector name, in which the hit has been measured

private:
  
  double getTrackPhi();
  double getTrackTheta();

  double m_resolutionRPhi; // Only used for virtual hits on non-modules
  double m_resolutionZ;    // Only used for virtual hits on non-modules

}; // Class

#endif /* INCLUDE_HITNEW_H_ */

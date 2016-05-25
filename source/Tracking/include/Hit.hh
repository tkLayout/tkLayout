/**
 * @file hit.hh
 * @brief This header file defines the hit and track classes used for internal analysis
 */


#ifndef _HIT_HH_
#define _HIT_HH_

#include "DetectorModule.h"
#include "PtErrorAdapter.h"
#include <MaterialProperties.h>
#include <cmath>
#include <vector>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <MessageLogger.h>

//using namespace ROOT::Math;
using namespace std;

class Track;

#undef HIT_DEBUG
#undef HIT_DEBUG_RZ

/**
 * @class Hit
 * @brief The Hit class is used when analysing a tracker layout to record information about a volume that was hit by a test track.
 *
 * It is used for both active and inactive surface hits. In case of an inactive surface, the pointer to the hit module will be <i>NULL</i>
 * since those volumes only matter with respect to the radiation and interaction lengths they add to the total in the error calculations.
 * All the other information is available to both categories. For convenience, the scaled radiation and interaction lengths are stored in
 * here as well to avoid additional computation and callbacks to the material property objects.
 */
class Hit {
protected:
  
  double distance_;      // distance of hit from origin in 3D
  double radius_;        // distance of hit from origin in the x/y plane
  int    orientation_;   // orientation of the surface
  int    objectKind_;    // kind of hit object
  
  DetectorModule* hitModule_;    // Pointer to the hit module
  Track*          myTrack_;      // Pointer to the track
  
  RILength correctedMaterial_;
  
  bool isPixel_;
  bool isTrigger_;
  bool isIP_;
  
  HitType activeHitType_;

private:
  
  double myResolutionRphi_; // Only used for virtual hits on non-modules
  double myResolutionY_;    // Only used for virtual hits on non-modules

public:
 
  Hit();
  ~Hit();
  Hit(const Hit& h);
  Hit(double myDistance);
  Hit(double myDistance, DetectorModule* myModule, HitType activeHitType);
  
  DetectorModule* getHitModule() { return hitModule_; };
  void setHitModule(DetectorModule* myModule);
  
  /**
   * @enum An enumeration of the category and orientation constants used within the object
   */
  enum { Undefined, Horizontal, Vertical,  // Hit object orientation 
    Active, Inactive };                    // Hit object type

  double   getDistance() {return distance_;};
  void     setDistance(double newDistance) { if (newDistance>0) distance_ = newDistance; updateRadius(); };
  double   getRadius() {return radius_;};
  void     updateRadius() {radius_ = distance_ * sin(getTrackTheta());};
  int      getOrientation() { return orientation_;};
  void     setOrientation(int newOrientation) { orientation_ = newOrientation; };
  int      getObjectKind() { return objectKind_;};
  void     setObjectKind(int newObjectKind) { objectKind_ = newObjectKind;};
  void     setTrack(Track* newTrack) {myTrack_ = newTrack; updateRadius();};
  double   getTrackTheta();
  RILength getCorrectedMaterial();
  void     setCorrectedMaterial(RILength newMaterial) { correctedMaterial_ = newMaterial;};
  
  bool   isPixel() { return isPixel_; };
  bool   isTrigger() { return isTrigger_; };
  bool   isIP() { return isIP_; };
  void   setPixel(bool isPixel) { isPixel_ = isPixel;}
  void   setTrigger(bool isTrigger) { isTrigger_ = isTrigger;}
  double getResolutionRphi(double trackR);
  double getResolutionZ(double trackR);
  void   setResolutionRphi(double newRes) { myResolutionRphi_ = newRes; } // Only used for virtual hits on non-modules
  void   setResolutionY(double newRes) { myResolutionY_ = newRes; } // Only used for virtual hits on non-modules
  bool   setIP(bool newIP) { return isIP_ = newIP; }

  bool isSquareEndcap();
  double getD();

  void setActiveHitType(HitType activeHitType) { activeHitType_ = activeHitType; }
  HitType getActiveHitType() const { return activeHitType_; } // NONE, INNER, OUTER, BOTH or STUB -- only meaningful for hits on active elements
  bool isStub() const { return activeHitType_ == HitType::STUB; }
};

/**
 * Given two hits, compare the distance to the z-axis.
 */
bool sortSmallerR(Hit* h1, Hit* h2);

#endif

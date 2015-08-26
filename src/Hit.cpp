/**
 * @file Hit.cpp
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "Hit.hh"
#include "Track.hh"
//#include "module.hh"
#include <global_constants.h>
#include <vector>
#include <algorithm>
#include <cstdlib>

using namespace ROOT::Math;
using namespace std;

// bool Track::debugRemoval = false; // debug
//#ifdef HIT_DEBUG_RZ
//bool Track::debugRZCovarianceMatrix = false;  // debug
//bool Track::debugRZCorrelationMatrix = false;  // debug
//bool Track::debugRZErrorPropagation = false;  // debug
//#endif

/**
 * This is a comparator for two Hit objects.
 * @param h1 A pointer to the first hit
 * @param h2 A pointer to the second hit
 * @return The result of the comparison: <i>true</i> if the distance from the z-axis of h1 is smaller than that of h2, false otherwise
 */
bool sortSmallerR(Hit* h1, Hit* h2) {
    return (h1->getDistance() < h2->getDistance());
}

/**
 * Nothing to do for the destructor, as a hit never owns any objects it has pointers to...
 */
Hit::~Hit() {}

/**
 * The default constructor sets the internal parameters to default values.
 */
Hit::Hit() {
    distance_   = 0;
    radius_     = 0;
    objectKind_ = Undefined;
    hitModule_  = NULL;
    orientation_= Undefined;
    myTrack_    = NULL;
    isPixel_    = false;
    isTrigger_  = false;
    isIP_       = false;
    myResolutionRphi_ = 0;
    myResolutionY_    = 0;
    activeHitType_    = HitType::NONE;
}

/**
 * The copy constructor makes sure the new object doesn't point to the old track (the track pointer needs to
 * be set explicitly later). The pointer to the module, on the other hand, stays the same as that of the original.
 */
Hit::Hit(const Hit& h) {
    distance_    = h.distance_;
    radius_      = h.radius_;
    orientation_ = h.orientation_;
    objectKind_  = h.objectKind_;
    hitModule_   = h.hitModule_;
    correctedMaterial_ = h.correctedMaterial_;
    myTrack_     = NULL;
    isPixel_     = h.isPixel_;
    isTrigger_   = h.isTrigger_;
    isIP_        = h.isIP_;
    myResolutionRphi_ = h.myResolutionRphi_;
    myResolutionY_    = h.myResolutionY_;
    activeHitType_    = h.activeHitType_;
}

/**
 * Constructor for a hit with no module at a given distance from the origin
 * @param myDistance distance from the origin
 */
Hit::Hit(double myDistance) {
    distance_    = myDistance;
    objectKind_  = Undefined;
    hitModule_   = NULL;
    orientation_ = Undefined;
    isTrigger_   = false;
    isPixel_     = false;
    isIP_        = false;
    myTrack_     = NULL;
    activeHitType_ = HitType::NONE;
}

/**
 * Constructor for a hit on a given module at a given distance from the origin
 * @param myDistance distance from the origin
 * @param myModule pointer to the module with the hit 
 */
Hit::Hit(double myDistance, Module* myModule, HitType activeHitType) {
    distance_    = myDistance;
    objectKind_  = Active;
    orientation_ = Undefined; 
    isTrigger_   = false;
    isPixel_     = false;
    isIP_        = false;
    setHitModule(myModule);
    myTrack_       = NULL;
    activeHitType_ = activeHitType;
}


/*
 * Setter for the pointer to the active surface that caused the hit.
 * @param myModule A pointer to a barrel or endcap module; may be <i>NULL</i>
 */
void Hit::setHitModule(Module* myModule) {
    if (myModule) {
        hitModule_ = myModule;
        if (myModule->subdet() == BARREL) {
            orientation_ = Horizontal;
        } else {
            orientation_ = Vertical;
        }
    }
}

/**
 * Get the track angle theta.
 * @return The angle from the z-axis of the entire track
 */
double Hit::getTrackTheta() {
    if (myTrack_==NULL) return 0;
    return (myTrack_->getTheta());
};

/**
 * Getter for the final, angle corrected pair of radiation and interaction lengths.
 * @return A copy of the pair containing the requested values; radiation length first, interaction length second
 */
RILength Hit::getCorrectedMaterial() {
    return correctedMaterial_;
}

/**
 * Getter for the rPhi resolution (local x coordinate for a module)
 * If the hit is not active it returns -1
 * If the hit is connected to a module, then the module's resolution
 * is retured (if the hit is trigger-type, then then module's trigger resultion is requested)
 * if there is not any hit module, then the hit's resolution property is read and returned
 * @return the hit's local resolution
 */
double Hit::getResolutionRphi(double trackR) {
  if (objectKind_!=Active) {
    std::cerr << "ERROR: Hit::getResolutionRphi called on a non-active hit" << std::endl;
    return -1;
  } else {
    if (hitModule_) {
      return hitModule_->resolutionEquivalentRPhi(getRadius(), trackR);
     // if (isTrigger_) return hitModule_->resolutionRPhiTrigger();
     // else return hitModule_->resolutionRPhi();
    } else {
      return myResolutionRphi_;
    }
  }
}

/**
 * Getter for the y resolution (local y coordinate for a module)
 * This corresponds to z coord for barrel modules and r coord for end-caps
 * If the hit is not active it returns -1
 * If the hit is connected to a module, then the module's resolution
 * is retured (if the hit is trigger-type, then then module's trigger resultion is requested)
 * if there is not any hit module, then the hit's resolution property is read and returned
 * @return the hit's local resolution
 */
double Hit::getResolutionZ(double trackR) {
  if (objectKind_!=Active) {
    std::cerr << "ERROR: Hit::getResolutionZ called on a non-active hit" << std::endl;
    return -1;
  } else {
    if (hitModule_) {
      return hitModule_->resolutionEquivalentZ(getRadius(), trackR, myTrack_->getCotgTheta());
      //if (isTrigger_) return hitModule_->resolutionYTrigger();
      //else return hitModule_->resolutionY();
    } else {
      return myResolutionY_;
    }
  }
}

/*
 * Checks wether a module belongs to the outer endcap (no pixel allowed)
 * and the hit module is made of a square sensor
 * @return true if the module is in outer endcap and square
 */
bool Hit::isSquareEndcap() {
  bool result = false;
  if (isPixel_) return false;
  //std::cout << "Hit::isSquareEndcap() "; //debug
  if (hitModule_) {
    //std::cout << " hitModule_!= NULL "; //debug
    if (hitModule_->subdet() == ENDCAP && hitModule_->shape() == RECTANGULAR) {
      //std::cout << " getSubdetectorType()==Endcap "; //debug
       //std::cout << " getShape()==Rectangular "; //debug
       result = true;
    }
  }
  //std::cout << std::endl; // debug
  return result;
}

/*
 * Retrieves the module's half width
 * for hit related to endcap modules only
 * @return Modules half width
 */
double Hit::getD() {
  double result = 0;
  //std::cout << "Hit::getD() "; //debug
  if (hitModule_) {
    //std::cout << " hitModule_!= NULL "; //debug
    EndcapModule* myECModule = hitModule_->as<EndcapModule>();
    if (myECModule) {
      //std::cout << " myECModule!= NULL "; //debug
      result = (myECModule->minWidth() + myECModule->maxWidth()) / 2. / 2.;
      //std::cout << " result = " << result; //debug
    }
  }
  //std::cout << std::endl; // debug
  return result;
}



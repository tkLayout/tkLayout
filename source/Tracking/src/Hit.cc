/**
 * @file Hit.cpp
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "Hit.h"

#include <global_constants.h>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "DetectorModule.h"
#include <MessageLogger.h>
#include "Track.h"

//using namespace ROOT::Math;
using namespace std;

// bool Track::debugRemoval = false; // debug
//#ifdef HIT_DEBUG_RZ
//bool Track::debugRZCovarianceMatrix = false;  // debug
//bool Track::debugRZCorrelationMatrix = false;  // debug
//bool Track::debugRZErrorPropagation = false;  // debug
//#endif

/**
 * This is a comparator for two Hit objects based on smaller radius.
 * @param h1 A pointer to the first hit
 * @param h2 A pointer to the second hit
 * @return The result of the comparison: <i>true</i> if the distance from the z-axis of h1 is smaller than that of h2, false otherwise
 */
bool Hit::sortSmallerR(const HitPtr& h1, const HitPtr& h2) { return (h1->getRPos() < h2->getRPos()); }

/**
 * This is a comparator for two Hit objects based on higher radius
 * @param h1 A pointer to the first hit
 * @param h2 A pointer to the second hit
 * @return The result of the comparison: <i>true</i> if the distance from the z-axis of h1 is smaller than that of h2, false otherwise
 */
bool Hit::sortHigherR(const HitPtr& h1, const HitPtr& h2) { return (h1->getRPos() > h2->getRPos()); }

/**
 * Nothing to do for the destructor, as a hit never owns any objects it has pointers to...
 */
Hit::~Hit() {}

/**
 * The default constructor sets the internal parameters to default values.
 */
Hit::Hit() {
    m_detName          = "Undefined";
    m_distance         = 0;
    m_rPos             = 0;
    m_zPos             = 0;
    m_activity         = HitActivity::Undefined;
    m_activeHitType    = HitType::NONE;
    m_hitModule        = nullptr;
    m_track            = nullptr;
    m_isTrigger        = false;
    m_isIP             = false;
    m_isBeamPipe       = false;
    m_layerID          = -1;
    m_discID           = -1;
    m_resolutionRPhi   = 0;
    m_resolutionZ      = 0;
}

/**
 * The copy constructor makes sure the new object doesn't point to the old track (the track pointer needs to
 * be set explicitly later). The pointer to the module, on the other hand, stays the same as that of the original.
 */
Hit::Hit(const Hit& h) {
    m_detName           = h.m_detName;
    m_distance          = h.m_distance;
    m_rPos              = h.m_rPos;
    m_zPos              = h.m_zPos;
    m_activity          = h.m_activity;
    m_activeHitType     = h.m_activeHitType;
    m_hitModule         = h.m_hitModule;
    m_track             = nullptr;
    m_correctedMaterial = h.m_correctedMaterial;
    m_isTrigger         = h.m_isTrigger;
    m_isIP              = h.m_isIP;
    m_isBeamPipe        = h.m_isBeamPipe;
    m_layerID           = h.m_layerID;
    m_discID            = h.m_discID;
    m_resolutionRPhi    = h.m_resolutionRPhi;
    m_resolutionZ       = h.m_resolutionZ;
}

/**
 * Constructor for a hit with no module at a given [rPos, zPos] from the origin
 */
Hit::Hit(double rPos, double zPos) {
    m_detName          = "Undefined";
    m_distance         = sqrt(rPos*rPos + zPos*zPos);
    m_rPos             = rPos;
    m_zPos             = zPos;
    m_activity         = HitActivity::Undefined;
    m_activeHitType    = HitType::NONE;
    m_hitModule        = nullptr;
    m_track            = nullptr;
    m_isTrigger        = false;
    m_isIP             = false;
    m_isBeamPipe       = false;
    m_layerID          = -1;
    m_discID           = -1;
    m_resolutionRPhi   = 0;
    m_resolutionZ      = 0;

}

/**
 * //! Constructor for a hit on a given module at [rPos, zPos] (cylindrical position) from the origin
 * @param myModule pointer to the module with the hit 
 */
Hit::Hit(double rPos, double zPos, const DetectorModule* myModule, HitType activeHitType) {
    m_detName          = "Undefined";
    m_distance         = sqrt(rPos*rPos + zPos*zPos);
    m_rPos             = rPos;
    m_zPos             = zPos;
    m_activity         = HitActivity::Active;
    m_activeHitType    = activeHitType;
    setHitModule(myModule);
    m_track            = nullptr;
    m_isTrigger        = false;
    m_isIP             = false;
    m_isBeamPipe       = false;
    m_layerID          = -1;
    m_discID           = -1;
    m_resolutionRPhi   = 0;
    m_resolutionZ      = 0;

}


/*
 * Setter for the pointer to the active surface that caused the hit.
 * @param myModule A pointer to a barrel or endcap module; may be <i>NULL</i>
 */
void Hit::setHitModule(const DetectorModule* myModule) {

  if (myModule) m_hitModule = myModule;
  else logWARNING("Hit::setHitModule -> can't set module to given hit, pointer null!");
}

/**
 * Get the track angle theta.
 * @return The angle from the z-axis of the entire track
 */
double Hit::getTrackTheta() {

  if (m_track==nullptr) {

    logWARNING("Hit::getTrackTheta -> no track assigned, will return zero!");
    return 0;
  }
  return (m_track->getTheta());
};

/**
 * Getter for the final, angle corrected pair of radiation and interaction lengths.
 * @return A copy of the pair containing the requested values; radiation length first, interaction length second
 */
RILength Hit::getCorrectedMaterial() {
    return m_correctedMaterial;
}

/**
 * Getter for the rPhi resolution (local x coordinate for a module)
 * If the hit is not active it returns -1
 * If the hit is connected to a module, then the module's resolution
 * is retured (if the hit is trigger-type, then then module's trigger resultion is requested)
 * if there is not any hit module, then the hit's resolution property is read and returned
 * @return the hit's local resolution
 */
double Hit::getResolutionRphi(double trackRadius) {

  if (!(this->isActive())) {

    logERROR("Hit::getResolutionRphi called on a non-active hit");
    return -1;
  }
  else {

    if (m_hitModule) {
      return m_hitModule->resolutionEquivalentRPhi(getRPos(), trackRadius);
     // if (isTrigger_) return hitModule_->resolutionRPhiTrigger();
     // else return hitModule_->resolutionRPhi();
    }
    else return m_resolutionRPhi;
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
double Hit::getResolutionZ(double trackRadius) {

  if (!(this->isActive())) {

    logERROR("Hit::getResolutionZ called on a non-active hit");
    return -1;
  }
  else {

    if (m_hitModule) {

      if (m_track==nullptr) {

        logWARNING("Hit::getResolutionZ -> no track assigned, will return zero!");
        return 0;
      }
      else return m_hitModule->resolutionEquivalentZ(getRPos(), trackRadius, m_track->getCotgTheta());
      //if (isTrigger_) return hitModule_->resolutionYTrigger();
      //else return hitModule_->resolutionY();
    }
    else return m_resolutionZ;
  }
}

/*
 * Checks wether a module belongs to the outer endcap (no pixel allowed)
 * and the hit module is made of a square sensor
 * @return true if the module is in outer endcap and square
 */
bool Hit::isSquareEndcap() {

  //std::cout << "Hit::isSquareEndcap() "; //debug

  if (m_hitModule) {
    //std::cout << " hitModule_!= NULL "; //debug
    if (m_hitModule->subdet() == ENDCAP && m_hitModule->shape() == RECTANGULAR) {
      //std::cout << " getSubdetectorType()==Endcap "; //debug
       //std::cout << " getShape()==Rectangular "; //debug
       return true;
    }
  }
  //std::cout << std::endl; // debug
  return false;
}

bool Hit::isStub() const
{
  return m_activeHitType == HitType::STUB;
}

/*
 * Retrieves the module's half width
 * for hit related to endcap modules only
 * @return Modules half width
 */
double Hit::getD() {

  double result = 0;
  //std::cout << "Hit::getD() "; //debug
  if (m_hitModule) {
    //std::cout << " hitModule_!= NULL "; //debug
    try {

      const EndcapModule* myECModule = dynamic_cast<const EndcapModule*>(m_hitModule);//->as<EndcapModule>();
      if (myECModule) {
        //std::cout << " myECModule!= NULL "; //debug
        result = (myECModule->minWidth() + myECModule->maxWidth()) / 2. / 2.;
        //std::cout << " result = " << result; //debug
      }
    }
    catch (exception& e) {}
  }
  //std::cout << std::endl; // debug
  return result;
}



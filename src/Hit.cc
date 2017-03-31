/**
<<<<<<< HEAD:src/Hit.cc:src/Hit.cc
 * @file Hit.cpp
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "Hit.hh"

#include <global_constants.hh>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "Track.hh"
#include "DetectorModule.hh"
#include "MessageLogger.hh"
#include "ModuleCap.hh"
#include "SimParms.hh"

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
    m_detName             = "Undefined";
    m_distance            = 0;
    m_rPos                = 0;
    m_zPos                = 0;
    m_activity            = HitActivity::Undefined;
    m_activeHitType       = HitType::NONE;
    m_hitModule           = nullptr;
    m_track               = nullptr;
    m_isTrigger           = false;
    m_passiveHitType      = HitPassiveType::Undefined;
    m_hitPassiveElem      = nullptr;
    m_isPixel             = false;
    m_resolutionRPhi      = 0;
    m_resolutionZ         = 0;
    m_isPixelIntersticeVol= false;
    m_isPixelTrackingVol  = false;
    m_isIntersticeVol     = false;
    m_isOuterTrackingVol  = false;
    m_isTotalTrackingVol  = false;
}

/**
 * The copy constructor makes sure the new object doesn't point to the old track (the track pointer needs to
 * be set explicitly later). The pointer to the module, on the other hand, stays the same as that of the original.
 */
Hit::Hit(const Hit& h) {
    m_detName             = h.m_detName;
    m_distance            = h.m_distance;
    m_rPos                = h.m_rPos;
    m_zPos                = h.m_zPos;
    m_activity            = h.m_activity;
    m_activeHitType       = h.m_activeHitType;
    m_hitModule           = h.m_hitModule;
    m_track               = nullptr;
    m_correctedMaterial   = h.m_correctedMaterial;
    m_isTrigger           = h.m_isTrigger;
    m_passiveHitType      = h.m_passiveHitType;
    m_hitPassiveElem      = h.m_hitPassiveElem;
    m_isPixel             = h.m_isPixel;
    m_resolutionRPhi      = h.m_resolutionRPhi;
    m_resolutionZ         = h.m_resolutionZ;
    m_isPixelIntersticeVol= false;
    m_isPixelTrackingVol  = false;
    m_isIntersticeVol     = false;
    m_isOuterTrackingVol  = false;
    m_isTotalTrackingVol  = false;
}

/**
 * Constructor for a hit on an inactive surface at a given [rPos, zPos] from the origin
 */
Hit::Hit(double rPos, double zPos, const insur::InactiveElement* myPassiveElem, HitPassiveType passiveHitType) {
    m_detName             = "Undefined";
    m_distance            = sqrt(rPos*rPos + zPos*zPos);
    m_rPos                = rPos;
    m_zPos                = zPos;
    m_activity            = HitActivity::Inactive;
    m_activeHitType       = HitType::NONE;
    m_hitModule           = nullptr;
    m_track               = nullptr;
    m_isTrigger           = false;
    m_passiveHitType      = passiveHitType;
    setHitPassiveElement(myPassiveElem);
    m_isPixel             = false;
    m_resolutionRPhi      = 0;
    m_resolutionZ         = 0;
    m_isPixelIntersticeVol= false;
    m_isPixelTrackingVol  = false;
    m_isIntersticeVol     = false;
    m_isOuterTrackingVol  = false;
    m_isTotalTrackingVol  = false;

    if (m_passiveHitType==HitPassiveType::BeamPipe) m_detName = "BeamPipe";
    if (m_passiveHitType==HitPassiveType::IP)       m_detName = "IP";
    if (m_passiveHitType==HitPassiveType::Support)  m_detName = "Support";
    if (m_passiveHitType==HitPassiveType::Service)  m_detName = "Service";
}

/**
 * Constructor for a hit on a given module at [rPos, zPos] (cylindrical position) from the origin
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
    m_passiveHitType   = HitPassiveType::Undefined;
    m_hitPassiveElem   = nullptr;
    m_isPixel          = false;
    m_resolutionRPhi   = 0;
    m_resolutionZ      = 0;

    if (myModule && myModule->getConstModuleCap()!=nullptr) m_detName = myModule->getConstModuleCap()->getDetName();
}


/*
 * Setter for the pointer to the active surface that caused the hit.
 * @param myModule A pointer to a barrel or endcap module; may be <i>NULL</i>
 */
void Hit::setHitModule(const DetectorModule* myModule) {

  if (myModule) m_hitModule = myModule;
  else logWARNING("Hit::setHitModule -> can't set module to given hit, pointer null!");
}

/*
 * Setter for the pointer to the inactive surface that caused the hit.
 */
void Hit::setHitPassiveElement(const insur::InactiveElement* myPassiveElem) {

  if (myPassiveElem) m_hitPassiveElem = myPassiveElem;
  //else logWARNING("Hit::setHitPassiveElement -> can't set inactive element to given hit, pointer null!");
}

/*
 * Get unique ID of layer or disc to which the hit belongs to (if not link return -1)
 */
int Hit::getLayerOrDiscID() const {

  if (m_hitModule && m_hitModule->getConstModuleCap()!=nullptr) return m_hitModule->getConstModuleCap()->getLayerOrDiscID(); else return -1;
}


/**
 * Get the track angle phi.
 * @return The angle from the z-axis of the entire track
 */
double Hit::getTrackPhi() {

  if (m_track==nullptr) {

    logWARNING("Hit::getTrackPhi -> no track assigned, will return zero!");
    return 0;
  }
  return (m_track->getPhi());
};

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

    // Module hit
    if (m_hitModule) {

      // R-Phi-resolution calculated as for a barrel-type module -> transform local R-Phi res. to a true module orientation (rotation by theta angle, skew, tilt)
      // In detail, take into account a propagation of MS error on virtual barrel plane, on which all measurements are evaluated for consistency (global chi2 fit applied) ->
      // in limit R->inf. propagation along line used, otherwise a very small correction factor coming from the circular shape of particle track is required (similar
      // approach as for local resolutions)
      // TODO: Currently, correction mathematicaly derived only for use case of const magnetic field -> more complex mathematical expression expected in non-const B field
      // (hence correction not applied in such case)
      double A = 0;
      if (SimParms::getInstance().isMagFieldConst()) A = getRPos()/(2*trackRadius); // r_i / 2R
      double B         = A/sqrt(1-A*A);
      double tiltAngle = m_hitModule->tiltAngle();
      double skewAngle = m_hitModule->skewAngle();
      double resLocalX = m_hitModule->resolutionLocalX(getTrackPhi());
      double resLocalY = m_hitModule->resolutionLocalY(getTrackTheta());

      // All modules & its resolution propagated to the resolution of a virtual barrel module (endcap is a tilted module by 90 degrees, barrel is tilted by 0 degrees)
      double resolutionRPhi = sqrt(pow((B*sin(skewAngle)*cos(tiltAngle) + cos(skewAngle)) * resLocalX,2) + pow(B*sin(tiltAngle) * resLocalY,2));

      return resolutionRPhi;
    }
    // IP or beam-constraint etc. hit
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

    // Module hit
    if (m_hitModule) {

      if (m_track==nullptr) {

        logWARNING("Hit::getResolutionZ -> no track assigned, will return zero!");
        return 0;
      }
      else {

        // Z-resolution calculated as for a barrel-type module -> transform local Z res. to a true module orientation (rotation by theta angle, skew, tilt)
        // In detail, take into account a propagation of MS error on virtual barrel plane, on which all measurements are evaluated for consistency (global chi2 fit applied) ->
        // in limit R->inf. propagation along line used, otherwise a very small correction factor coming from the circular shape of particle track is required (similar
        // approach as for local resolutions)
        // TODO: Currently, correction mathematicaly derived only for use case of const magnetic field -> more complex mathematical expression expected in non-const B field
        // (hence correction not applied in such case)
        double A = 0;
        if (SimParms::getInstance().isMagFieldConst()) A = getRPos()/(2*trackRadius);
        double D         = m_track->getCotgTheta()/sqrt(1-A*A);
        double tiltAngle = m_hitModule->tiltAngle();
        double skewAngle = m_hitModule->skewAngle();
        double resLocalX = m_hitModule->resolutionLocalX(getTrackPhi());
        double resLocalY = m_hitModule->resolutionLocalY(getTrackTheta());

        // All modules & its resolution propagated to the resolution of a virtual barrel module (endcap is a tilted module by 90 degrees, barrel is tilted by 0 degrees)
        double resolutionZ = sqrt(pow(((D*cos(tiltAngle) + sin(tiltAngle))*sin(skewAngle)) * resLocalX,2) + pow((D*sin(tiltAngle) + cos(tiltAngle)) * resLocalY,2));

        return resolutionZ;
      }
    }
    // IP or beam-constraint etc. hit
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



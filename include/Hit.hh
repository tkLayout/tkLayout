/**
 * @file hit.hh
 * @brief This header file defines the hit and track classes used for internal analysis
 */


#ifndef _HIT_HH_
#define _HIT_HH_

#include "Module.hh"
#include "PtErrorAdapter.hh"
#include <MaterialProperties.hh>
#include <InactiveElement.hh>
#include <cmath>
#include <vector>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <MessageLogger.hh>


#include <TFile.h>
#include <TProfile.h>
#include <TF1.h>
#include <TAxis.h>
#include <TCanvas.h>

//using namespace ROOT::Math;
using namespace std;

class Track;
namespace insur {
  class InactiveElement;
}
using insur::InactiveElement;

// TODO: why this?
typedef double momentum;  // Track momentum in MeV


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
  double distance_;   // distance of hit from origin in 3D
  double radius_; // distance of hit from origin in the x/y plane
  int orientation_;   // orientation of the surface
  int objectKind_;    // kind of hit object
  int objectCategory_;
  Module* hitModule_; // Pointer to the hit module
  InactiveElement* hitInactiveElement_; // Pointer to the hit inactive element
  //double trackTheta_; // Theta angle of the track
  //Material material_;
  // "Thickness" in terms of radiation_length and interaction_length
  RILength correctedMaterial_;
  Track* myTrack_;
  bool isPixel_;
  bool isPixelIntersticeVolume_;
  bool isPixelTrackingVolume_;
  bool isIntersticeVolume_;
  bool isOuterTrackingVolume_;
  bool isTotalTrackingVolume_;
  bool isTrigger_;
  bool isIP_;
  HitType activeHitType_;

private:
  double resolutionLocalX_; // Only used for hits on active
  double resolutionLocalY_; // Only used for hits on active
  double myResolutionRphi_; // Only used for virtual hits on non-modules
  double myResolutionY_;    // Only used for virtual hits on non-modules

public:
  ~Hit();
  Hit();
  Hit(const Hit& h);
  Hit(double myDistance);
  Hit(double myDistance, Module* myModule, HitType activeHitType);
  Module* getHitModule() { return hitModule_; };
  InactiveElement* getHitInactiveElement() { return hitInactiveElement_; };
  void computeLocalResolution();
  double getResolutionRphi(double trackR);
  double getResolutionZ(double trackR);
  void setHitModule(Module* myModule);
  void setHitInactiveElement(InactiveElement* myInactiveElement);
  /**
   * @enum An enumeration of the category and orientation constants used within the object
   */
  enum { Undefined, Horizontal, Vertical,  // Hit object orientation 
    Active, Inactive };               // Hit object type
  enum {Unknown, Act, BeamPipe, Service, Support };

  double getDistance() {return distance_;};
  void setDistance(double newDistance) { if (newDistance>0) distance_ = newDistance; updateRadius(); };
  double getRadius() {return radius_;};
  void updateRadius() {radius_ = distance_ * sin(getTrackTheta());};
  int getOrientation() { return orientation_;};
  void setOrientation(int newOrientation) { orientation_ = newOrientation; };
  int getObjectKind() { return objectKind_;};
  int getObjectCategory() { return objectCategory_;};
  void setObjectKind(int newObjectKind) { objectKind_ = newObjectKind;};
  void setObjectCategory(int newObjectCategory) { objectCategory_ = newObjectCategory;};
  void setTrack(Track* newTrack) {myTrack_ = newTrack; updateRadius();};
  double getTrackTheta();
  RILength getCorrectedMaterial();
  void setCorrectedMaterial(RILength newMaterial) { correctedMaterial_ = newMaterial;};
  bool isPixel() { return isPixel_; };
  bool isPixelIntersticeVolume() { return isPixelIntersticeVolume_; };
  bool isPixelTrackingVolume() { return isPixelTrackingVolume_; };
  bool isIntersticeVolume() { return isIntersticeVolume_; };
  bool isOuterTrackingVolume() { return isOuterTrackingVolume_; };
  bool isTotalTrackingVolume() { return isTotalTrackingVolume_; };
  bool isTrigger() { return isTrigger_; };
  bool isIP() { return isIP_; };
  void setPixel(bool isPixel) { isPixel_ = isPixel;}
  void setPixelIntersticeVolume(bool isPixelIntersticeVolume) { isPixelIntersticeVolume_ = isPixelIntersticeVolume; }
  void setPixelTrackingVolume(bool isPixelTrackingVolume) { isPixelTrackingVolume_ = isPixelTrackingVolume; }
  void setIntersticeVolume(bool isIntersticeVolume) { isIntersticeVolume_ = isIntersticeVolume; }
  void setOuterTrackingVolume(bool isOuterTrackingVolume) { isOuterTrackingVolume_ = isOuterTrackingVolume; }
  void setTotalTrackingVolume(bool isTotalTrackingVolume) { isTotalTrackingVolume_ = isTotalTrackingVolume; }
  void setTrigger(bool isTrigger) { isTrigger_ = isTrigger;}
  double getResolutionLocalX() { return resolutionLocalX_; }
  double getResolutionLocalY() { return resolutionLocalY_; }
  void setResolutionRphi(double newRes) { myResolutionRphi_ = newRes; } // Only used for virtual hits on non-modules
  void setResolutionY(double newRes) { myResolutionY_ = newRes; } // Only used for virtual hits on non-modules
  bool setIP(bool newIP) { return isIP_ = newIP; }

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

/**
 * @class Track
 * @brief The Track class is essentially a collection of consecutive hits.
 *
 * Once those hits have been stored, though, it also provides a series of error analysis functions that use the information about
 * radiation and interaction length from the hits as a basis for the calculations.
 */
class Track {
protected:
  double theta_;
  double phi_;
  double cotgTheta_, eta_; // calculated from theta and then cached
  std::vector<Hit*> hitV_;
  // Track resolution as a function of momentum
// WARNING: Use new TrackNew class for tracking  (class currently kept for compatibility)
//  TMatrixTSym<double> correlations_;
//  TMatrixT<double> covariances_;
//  TMatrixTSym<double> correlationsRZ_;
//  TMatrixT<double> covariancesRZ_;
//  double deltarho_;
//  double deltaphi_;
//  double deltad_;
//  double deltaCtgTheta_;
//  double deltaZ0_;
//  double deltaP_;
  void computeLocalResolution();
// WARNING: Use new TrackNew class for tracking  (class currently kept for compatibility)
//  void computeCorrelationMatrixRZ();
//  void computeCovarianceMatrixRZ();
//  void computeCorrelationMatrix();
//  void computeCovarianceMatrix();
  
  std::set<std::string> tags_;
  double transverseMomentum_;
public:
  Track();
  Track(const Track& t);
  ~Track();
  Track& operator=(const Track &t);
  std::vector<Hit*> getHitV() const { return hitV_; }
  bool noHits() { return hitV_.empty(); }
  int nHits() { return hitV_.size(); }
  double setTheta(double& newTheta);
  double getTheta() const {return theta_;}
  double getEta() const { return eta_; } // calculated when theta is set, then cached
  double getCotgTheta() const { return cotgTheta_; } // ditto here
  double setPhi(double& newPhi);
  double getPhi() const {return phi_;}
// WARNING: Use new TrackNew class for tracking  (class currently kept for compatibility)
//  TMatrixTSym<double>& getCorrelations() { return correlations_; }
//  TMatrixT<double>& getCovariances() { return covariances_; }
//  const double& getDeltaRho() const { return deltarho_; }
//  const double& getDeltaPhi() const { return deltaphi_; }
//  const double& getDeltaD() const { return deltad_; }
//  const double& getDeltaCtgTheta() const { return deltaCtgTheta_; }
//  const double& getDeltaZ0() const { return deltaZ0_; }
//  const double& getDeltaP() const { return deltaP_; }

  Hit* addHit(Hit* newHit);
  const std::set<std::string>& tags() const { return tags_; }
  void sort();
  void assignTrackingVolumesToHits();
// WARNING: Use new TrackNew class for tracking  (class currently kept for compatibility)
//  void computeErrors();
//  void printErrors();
  void print();
  void removeMaterial();
  int nActiveHits(bool usePixels = false, bool useIP = true) const;
  std::vector<double> hadronActiveHitsProbability(bool usePixels = false);
  double hadronActiveHitsProbability(int nHits, bool usePixels = false);
  void addEfficiency();
  void keepTriggerOnly();
  void keepTaggedOnly(const string& tag);
  void setTriggerResolution(bool isTrigger);
  // static bool debugRemoval; // debug
  double expectedTriggerPoints(const double& triggerMomentum) const;
#ifdef HIT_DEBUG_RZ
  static bool debugRZCovarianceMatrix;  // debug
  static bool debugRZCorrelationMatrix;  // debug
  static bool debugRZErrorPropagation;  // debug
#endif
  void addIPConstraint(double dr, double dz);
  RILength getCorrectedMaterial();
  std::vector<std::pair<Module*, HitType>> getHitModules() const;

// WARNING: Use new TrackNew class for tracking  (class currently kept for compatibility)
//  void setTransverseMomentum(const double newPt);
//  double getTransverseMomentum() const { return transverseMomentum_; }

//  void pruneHits();
};
#endif

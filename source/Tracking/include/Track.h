/**
 * @file Track.hh
 * @brief This header file defines track class used for internal analysis
 */

#ifndef INCLUDE_TRACK_H_
#define INCLUDE_TRACK_H_

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

class Hit;

#undef TRACK_DEBUG_RZ

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
  double cotgTheta_; 
  double eta_;                // calculated from theta and then cached
  double transverseMomentum_;
  double magField_;

  std::vector<Hit*> hitV_;

  // Track resolution as a function of momentum
  TMatrixTSym<double> correlations_;
  TMatrixT<double>    covariances_;
  TMatrixTSym<double> correlationsRZ_;
  TMatrixT<double>    covariancesRZ_;
  double deltarho_;
  double deltaphi_;
  double deltad_;
  double deltaCtgTheta_;
  double deltaZ0_;
  double deltaP_;
  void computeCorrelationMatrixRZ();
  void computeCovarianceMatrixRZ();
  void computeCorrelationMatrix();
  void computeCovarianceMatrix();
  
  std::set<std::string> tags_;

public:

  Track();
  Track(const Track& t);
  ~Track();
  Track& operator=(const Track &t);

  bool   noHits() { return hitV_.empty(); }
  int    nHits() { return hitV_.size(); }
  
  double setTheta(double& newTheta);
  double getTheta() const {return theta_;}
  double getEta() const { return eta_; } // calculated when theta is set, then cached
  double getCotgTheta() const { return cotgTheta_; } // ditto here
  double setPhi(double& newPhi);
  double getPhi() const {return phi_;}
  void   setTransverseMomentum(const double newPt) { transverseMomentum_ = newPt; }
  double getTransverseMomentum() const { return transverseMomentum_; }
  bool   pruneHits(); // Remove hits that don't follow the parabolic approximation used in tracking - still needs to be updated (not all approximations taken into account)
  void   setMagField(const double & magField) { magField_ = magField; }
  double getMagField() const;
  double getRho() const;
  double getRadius() const;
  
  TMatrixTSym<double>& getCorrelations() { return correlations_; }
  TMatrixT<double>&    getCovariances()  { return covariances_; }
  
  const double& getDeltaRho() const { return deltarho_; }
  const double& getDeltaPhi() const { return deltaphi_; }
  const double& getDeltaD() const { return deltad_; }
  const double& getDeltaCtgTheta() const { return deltaCtgTheta_; }
  const double& getDeltaZ0() const { return deltaZ0_; }
  const double& getDeltaP() const { return deltaP_; }

  Hit* addHit(Hit* newHit);
  const std::set<std::string>& tags() const { return tags_; }
  void sort();
  void computeErrors();
  void printErrors();
  void print();
  void removeMaterial();
  int nActiveHits(bool usePixels = false, bool useIP = true) const;
  std::vector<double> hadronActiveHitsProbability(bool usePixels = false);
  double hadronActiveHitsProbability(int nHits, bool usePixels = false);
  void addEfficiency(double efficiency, bool alsoPixel = false);
  void keepTriggerOnly();
  void keepTaggedOnly(const string& tag);
  void setTriggerResolution(bool isTrigger);
  // static bool debugRemoval; // debug
  double expectedTriggerPoints(const double& triggerMomentum) const;
#ifdef TRACK_DEBUG_RZ
  static bool debugRZCovarianceMatrix;  // debug
  static bool debugRZCorrelationMatrix;  // debug
  static bool debugRZErrorPropagation;  // debug
#endif
  void addIPConstraint(double dr, double dz);
  RILength getCorrectedMaterial();
  std::vector<std::pair<DetectorModule*, HitType>> getHitModules() const;

};
#endif /* INCLUDE_TRACK_H_ */

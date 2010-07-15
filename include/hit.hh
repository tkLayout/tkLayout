#ifndef _HIT_HH_
#define _HIT_HH_

#include "module.hh"
#include <cmath>
#include <vector>
#include <map>
#include<TMatrixT.h>
#include <TMatrixTSym.h>

//using namespace ROOT::Math;
using namespace std;

class Track;

typedef double momentum;  // Track momentum in MeV

class Hit {
protected:
  double distance_;   // distance of hit from origin in 3D
  double radius_; // distance of hit from origin in the x/y plane
  int orientation_;   // orientation of the surface
  int objectKind_;    // kind of hit object
  Module* hitModule_; // Pointer to the hit module
  //double trackTheta_; // Theta angle of the track
  //pair<double, double> material_;
  // "Thickness" in terms of radiation_length and interaction_length
  pair<double, double> correctedMaterial_;
  Track* myTrack_;
  // Error of the measurement due to scattering in mm
  // as a function of the track momentum
  std::map<momentum ,double> rphiError_ ;

private:
  
public:
  ~Hit();
  Hit();
  Hit(const Hit& h);
  Hit(double myDistance);
  Hit(double myDistance, Module* myModule);
  Module* getHitModule() { return hitModule_; };
  void setHitModule(Module* myModule);
  enum { Undefined, Horizontal, Vertical,  // Hit object orientation 
	 Active, Inactive };               // Hit object type

  double getDistance() {return distance_;};
  void setDistance(double newDistance) { if (newDistance>0) distance_ = newDistance; updateRadius(); };
  double getRadius() {return radius_;};
  void updateRadius() {radius_ = distance_ * sin(getTrackTheta());};
  int getOrientation() { return orientation_;};
  void setOrientation(int newOrientation) { orientation_ = newOrientation; };
  int getObjectKind() { return objectKind_;};
  void setObjectKind(int newObjectKind) { objectKind_ = newObjectKind;};
  void setTrack(Track* newTrack) {myTrack_ = newTrack; updateRadius();};
  double getTrackTheta();
  pair<double, double> getCorrectedMaterial();
  void setCorrectedMaterial(pair<double, double> newMaterial) { correctedMaterial_ = newMaterial;};
  void setError(double rphiError, momentum p) { rphiError_[p] = rphiError; };
};

bool sortSmallerR(Hit* h1, Hit* h2);

class Track {
protected:
  double theta_;
  std::vector<Hit*> hitV_;
  // Track resolution as a function of momentum
  map<momentum, TMatrixTSym<double> > correlations_;
  map<momentum, TMatrixT<double> > covariances_;
  map<momentum, double> deltarho_;
  map<momentum, double> deltaphi_;
  map<momentum, double> deltad_;
public:
  Track();
  Track(const Track& t);
  ~Track();
  double setTheta(double& newTheta);
  double getTheta() {return theta_;};
  map<momentum, TMatrixTSym<double> >& getCorrelations() { return correlations_; };
  map<momentum, TMatrixT<double> >& getCovariances() { return covariances_; };
  map<momentum, double>& getDeltaRho() { return deltarho_; };
  map<momentum, double>& getDeltaPhi() { return deltaphi_; };
  map<momentum, double>& getDeltaD() { return deltad_; };
  Hit* addHit(Hit* newHit) {hitV_.push_back(newHit); newHit->setTrack(this); return newHit;};
  void sort();
  void computeCorrelationMatrix(const vector<double>& momenta); //TODO: move to private after testing
  void computeCovarianceMatrix(const map<double, TMatrixTSym<double> >& correlations); //TODO: move to private after testing
  void computeErrors(const std::vector<momentum>& momentaList); //TODO: should call computeCorrelationMatrix() and computeCovarianceMatrix() internally
};
#endif

#ifndef _HIT_HH_
#define _HIT_HH_

#include "module.hh"
#include <vector>

//using namespace ROOT::Math;
using namespace std;

class Track;

class Hit {
protected:
  double distance_;   // distance of hit from origin
  int orientation_;   // orientation of the surface
  int objectKind_;    // kind of hit object
  Module* hitModule_; // Pointer to the hit module
  //double trackTheta_; // Theta angle of the track
  //pair<double, double> material_;
  pair<double, double> correctedMaterial_;
  Track* myTrack_;

private:
  
public:
  ~Hit();
  Hit();
  Hit(double myDistance);
  Hit(double myDistance, Module* myModule);
  void setHitModule(Module* myModule);
  enum { Undefined, Horizontal, Vertical,  // Hit object orientation 
	 Active, Inactive };      // Hit object type

  double getDistance() {return distance_;};
  void setDistance(double newDistance) { if (newDistance>0) distance_ = newDistance;};
  int getOrientation() { return orientation_;};
  void setOrientation(int newOrientation) { orientation_ = newOrientation; };
  int getObjectKind() { return objectKind_;};
  void setObjectKind(int newObjectKind) { objectKind_ = newObjectKind;};
  void setTrack(Track* newTrack) {myTrack_ = newTrack;};
  double getTrackTheta();
  pair<double, double> getCorrectedMaterial();
  void setCorrectedMaterial(pair<double, double> newMaterial) { correctedMaterial_ = newMaterial;};
 
};

bool sortSmallerR(Hit* h1, Hit* h2) {
  return (h1->getDistance() < h2->getDistance());
}

class Track {
protected:
  double theta_;
  std::vector<Hit*> hitV_;
public:
  Track() {};
  ~Track() {};
  double setTheta(double& newTheta) {theta_ = newTheta; return theta_;};
  double getTheta() {return theta_;};
  Hit* addHit(Hit* newHit) {hitV_.push_back(newHit); newHit->setTrack(this); return newHit;};
  void sort();
};

#endif

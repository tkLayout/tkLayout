#ifndef _HIT_HH_
#define _HIT_HH_

#include "module.hh"

//using namespace ROOT::Math;
using namespace std;

class Hit {
protected:
  double distance_;   // distance of hit from origin
  int orientation_;   // orientation of the surface
  int objectKind_;    // kind of hit object
  Module* hitModule_; // Pointer to the hit module
  double trackTheta_; // Theta angle of the track
  //pair<double, double> material_;
  pair<double, double> correctedMaterial_;
private:
  
public:
  ~Hit();
  Hit();
  Hit(double myDistance);
  Hit(double myDistance, Module* myModule);
  Hit(double myDistance, Module* myModule, double trackTheta_);
  void setHitModule(Module* myModule);
  enum { Undefined, Horizontal, Vertical,  // Hit object orientation 
	 Active, Inactive };      // Hit object type

  double getDistance() {return distance_;};
  void setDistance(double newDistance) { if (newDistance>0) distance_ = newDistance;};
  int getOrientation() { return orientation_;};
  void setOrientation(int newOrientation) { orientation_ = newOrientation; };
  int getObjectKind() { return objectKind_;};
  void setObjectKind(int newObjectKind) { objectKind_ = newObjectKind;};
  double getTrackTheta() { return trackTheta_;};
  void setTrackTheta(double newTrackTheta) { trackTheta_ = newTrackTheta;};
  //pair<double, double> getBareMaterial();
  //void setBareMaterial(pair<double, double> newMaterial) { material_ = newMaterial;};
  pair<double, double> getCorrectedMaterial();
  void setCorrectedMaterial(pair<double, double> newMaterial) { correctedMaterial_ = newMaterial;};
 
};

#endif

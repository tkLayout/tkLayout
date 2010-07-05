#include "hit.hh"
#include "module.hh"

using namespace ROOT::Math;
using namespace std;

Hit::~Hit() {
}

Hit::Hit() {
  distance_ = 0;
  objectKind_ = Undefined;
  hitModule_ = NULL;
  orientation_ = Undefined;
  trackTheta_ = 0;
}

Hit::Hit(double myDistance) {
  distance_ = myDistance;
  objectKind_ = Undefined;
  hitModule_ = NULL;
  orientation_ = Undefined;
  trackTheta_ = 0;
}

Hit::Hit(double myDistance, Module* myModule) {
  distance_ = myDistance;
  objectKind_ = Active;
  hitModule_ = myModule;
  orientation_ = Undefined;
  trackTheta_ = 0;
  setHitModule(myModule);
}

Hit::Hit(double myDistance, Module* myModule, double newTrackTheta) {
  distance_ = myDistance;
  objectKind_ = Active;
  hitModule_ = myModule;
  orientation_ = Undefined;
  trackTheta_ = newTrackTheta;
  setHitModule(myModule);
}

void Hit::setHitModule(Module* myModule) {
  if (myModule) {
    int subDetectorType = myModule->getSubdetectorType();
    if (subDetectorType==Module::Barrel) {
      orientation_ = Horizontal;
    } else if (subDetectorType==Module::Endcap) {
      orientation_ = Vertical;
    } else {
      cerr << "ERROR: a generic module was assigned to a hit. This should not happen!" << endl;
    }
  }
}

//pair<double, double> Hit::getBareMaterial() {
//  return material_;
//}

pair<double, double> Hit::getCorrectedMaterial() {
 return correctedMaterial_;
 /* pair<double, double> correctedMaterial;
  double factor=0;

  if (orientation_ == Horizontal) {
    if (sin(trackTheta_)==0) {
      cerr << "ERROR: impossible to have a hit at theta=0 or theta=M_PI" << endl;
      factor = 0;
    } else factor = 1 / sin(trackTheta_);
  } else if (orientation_ == Vertical) {
    if (cos(trackTheta_)==0) {
      cerr << "ERROR: impossible to have a hit on a vertical surface at theta=M_PI/2" << endl;
      factor = 0;
    } else factor = 1 / cos(trackTheta_);
  }
  
  correctedMaterial.first = material_.first * factor;
  correctedMaterial.second = material_.second * factor;

  return correctedMaterial; */
}

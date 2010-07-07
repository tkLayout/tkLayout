#include "hit.hh"
#include "module.hh"
#include <vector>
#include <algorithm>

using namespace ROOT::Math;
using namespace std;

bool sortSmallerR(Hit* h1, Hit* h2) {
  return (h1->getDistance() < h2->getDistance());
}

Hit::~Hit() {
}

Hit::Hit() {
  distance_ = 0;
  objectKind_ = Undefined;
  hitModule_ = NULL;
  orientation_ = Undefined;
  //trackTheta_ = 0;
  myTrack_ = NULL;
}

Hit::Hit(const Hit& h) {
    distance_ = h.distance_;
    orientation_ = h.orientation_;
    objectKind_ = h.objectKind_;
    hitModule_ = h.hitModule_;
    correctedMaterial_ = h.correctedMaterial_;
    myTrack_ = NULL;
    //TODO: rphiError_ ???
}

Hit::Hit(double myDistance) {
  distance_ = myDistance;
  objectKind_ = Undefined;
  hitModule_ = NULL;
  orientation_ = Undefined;
  //trackTheta_ = 0;
  myTrack_ = NULL;
}

Hit::Hit(double myDistance, Module* myModule) {
  distance_ = myDistance;
  objectKind_ = Active;
  orientation_ = Undefined;
  //trackTheta_ = 0;
  setHitModule(myModule);
  myTrack_ = NULL;
}



void Hit::setHitModule(Module* myModule) {
  if (myModule) {
    hitModule_ = myModule;
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

double Hit::getTrackTheta() {
  if (myTrack_==NULL)
    return 0;
  return (myTrack_->getTheta());
};

pair<double, double> Hit::getCorrectedMaterial() {
 return correctedMaterial_;
}

Track::Track() {
  theta_ = 0;
}

Track::Track(const Track& t) {
    theta_ = t.theta_;
    vector<Hit*>::const_iterator iter, guard = t.hitV_.end();
    for (iter = t.hitV_.begin(); iter != guard; iter++) {
        Hit* h = new Hit(*(*iter));
        addHit(h);
    }
}

Track::~Track() {
  std::vector<Hit*>::iterator hitIt;
  for (hitIt=hitV_.begin(); hitIt!=hitV_.end(); hitIt++) {
    if ((*hitIt)!=NULL) {
      delete (*hitIt);
    }
  }
  hitV_.clear();
}

void Track::sort() {
  std::sort(hitV_.begin(), hitV_.end(), sortSmallerR);
}

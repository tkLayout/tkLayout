/**
 * @file Hit.cc
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "Hit.hh"
//#include "module.hh"
#include <global_constants.hh>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "SimParms.hh"

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
    distance_ = 0;
    radius_ = 0;
    objectKind_ = Undefined;
    objectCategory_ = Unknown;
    hitModule_ = NULL;
    hitInactiveElement_ = NULL;
    orientation_ = Undefined;
    myTrack_ = NULL;
    isPixel_ = false;
    isPixelIntersticeVolume_ = false;
    isPixelTrackingVolume_ = false;
    isIntersticeVolume_ = false;
    isOuterTrackingVolume_ = false;
    isTotalTrackingVolume_ = false;
    isTrigger_ = false;
    isIP_ = false;
    resolutionLocalX_ = 0;
    resolutionLocalY_ = 0;
    myResolutionRphi_ = 0;
    myResolutionY_ = 0;
    activeHitType_ = HitType::NONE;
}

/**
 * The copy constructor makes sure the new object doesn't point to the old track (the track pointer needs to
 * be set explicitly later). The pointer to the module, on the other hand, stays the same as that of the original.
 */
Hit::Hit(const Hit& h) {
    distance_ = h.distance_;
    radius_ = h.radius_;
    orientation_ = h.orientation_;
    objectKind_ = h.objectKind_;
    objectCategory_ = h.objectCategory_;
    hitModule_ = h.hitModule_;
    hitInactiveElement_ = h.hitInactiveElement_;
    correctedMaterial_ = h.correctedMaterial_;
    myTrack_ = NULL;
    isPixel_ = h.isPixel_;
    isPixelIntersticeVolume_ = h.isPixelIntersticeVolume_;
    isPixelTrackingVolume_ = h.isPixelTrackingVolume_;
    isIntersticeVolume_ = h.isIntersticeVolume_;
    isOuterTrackingVolume_ = h.isOuterTrackingVolume_;
    isTotalTrackingVolume_ = h.isTotalTrackingVolume_;
    isTrigger_ = h.isTrigger_;
    isIP_ = h.isIP_;
    resolutionLocalX_ = h.resolutionLocalX_;
    resolutionLocalY_ = h.resolutionLocalY_;
    myResolutionRphi_ = h.myResolutionRphi_;
    myResolutionY_ = h.myResolutionY_;
    activeHitType_ = h.activeHitType_;
}

/**
 * Constructor for a hit with no module at a given distance from the origin
 * @param myDistance distance from the origin
 */
Hit::Hit(double myDistance) {
    distance_ = myDistance;
    objectKind_ = Undefined;
    objectCategory_ = Unknown;
    hitModule_ = NULL;
    hitInactiveElement_ = NULL;
    orientation_ = Undefined;
    isTrigger_ = false;
    isPixel_ = false;
    isPixelIntersticeVolume_ = false;
    isPixelTrackingVolume_ = false;
    isIntersticeVolume_ = false;
    isOuterTrackingVolume_ = false;
    isTotalTrackingVolume_ = false;
    isIP_ = false;
    myTrack_ = NULL;
    activeHitType_ = HitType::NONE;
}

/**
 * Constructor for a hit on a given module at a given distance from the origin
 * @param myDistance distance from the origin
 * @param myModule pointer to the module with the hit 
 */
Hit::Hit(double myDistance, Module* myModule, HitType activeHitType) {
    distance_ = myDistance;
    objectKind_ = Active;
    objectCategory_ = Act;
    hitInactiveElement_ = NULL;
    orientation_ = Undefined; 
    isTrigger_ = false;
    isPixel_ = false;
    isPixelIntersticeVolume_ = false;
    isPixelTrackingVolume_ = false;
    isIntersticeVolume_ = false;
    isOuterTrackingVolume_ = false;
    isTotalTrackingVolume_ = false;
    isIP_ = false;
    setHitModule(myModule);
    myTrack_ = NULL;
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

/*
 * Setter for the pointer to the inactive surface that caused the hit.
 * @param myInactiveElement A pointer to an inactive element (service or support structure); may be <i>NULL</i>
 */
void Hit::setHitInactiveElement(InactiveElement* myInactiveElement) {
    if (myInactiveElement) {
        hitInactiveElement_ = myInactiveElement;
        if (myInactiveElement->isVertical()) {
            orientation_ = Vertical;
        } else {
            orientation_ = Horizontal;
        }
    }
}

/**
 * Get the track angle theta.
 * @return The angle from the z-axis of the entire track
 */
double Hit::getTrackTheta() {
    if (myTrack_==NULL)
        return 0;
    return (myTrack_->getTheta());
};

/**
 * Getter for the final, angle corrected pair of radiation and interaction lengths.
 * @return A copy of the pair containing the requested values; radiation length first, interaction length second
 */
RILength Hit::getCorrectedMaterial() {
    return correctedMaterial_;
}


void Hit::computeLocalResolution() {
  if (objectKind_!= Active) {
    std::cerr << "ERROR: Hit::computeLocalResolution called on a non-active hit" << std::endl;
  } else {
    if (hitModule_) {
      resolutionLocalX_ = hitModule_->resolutionLocalX(myTrack_->getPhi());
      resolutionLocalY_ = hitModule_->resolutionLocalY(myTrack_->getTheta());
      
      if (hitModule_->hasAnyResolutionLocalXParam()) hitModule_->rollingParametrizedResolutionLocalX(resolutionLocalX_);
      if (hitModule_->hasAnyResolutionLocalYParam()) hitModule_->rollingParametrizedResolutionLocalY(resolutionLocalY_);
    }
  }
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
      return hitModule_->resolutionEquivalentRPhi(getRadius(), trackR, resolutionLocalX_, resolutionLocalY_);
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
      return hitModule_->resolutionEquivalentZ(getRadius(), trackR, myTrack_->getCotgTheta(), resolutionLocalX_, resolutionLocalY_);
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


/**
 * The default constructor sets the parameter for the track angle to zero.
 */
Track::Track() {
    theta_ = 0;
}

/**
 * The copy constructor creates a deep copy of the vector of hits.
 */
Track::Track(const Track& t) {
  theta_ = t.theta_;
  phi_ = t.phi_;
  cotgTheta_ = t.cotgTheta_;
  eta_ = t.eta_;
//  correlations_.ResizeTo(t.correlations_);
//  correlations_ = t.correlations_;
//  covariances_.ResizeTo(t.covariances_);
//  covariances_ = t.covariances_;
//  correlationsRZ_.ResizeTo(t.correlationsRZ_);
//  correlationsRZ_ = t.correlationsRZ_;
//  covariancesRZ_.ResizeTo(t.covariancesRZ_);
//  covariancesRZ_ = t.covariancesRZ_;
//  deltarho_ = t.deltarho_;
//  deltaphi_ = t.deltaphi_;
//  deltad_ = t.deltad_;
//  deltaCtgTheta_ = t.deltaCtgTheta_;
//  deltaZ0_ = t.deltaZ0_;
//  deltaP_ = t.deltaP_;
  vector<Hit*>::const_iterator iter, guard = t.hitV_.end();
  for (iter = t.hitV_.begin(); iter != guard; iter++) {
    Hit* h = new Hit(*(*iter));
    addHit(h);
  }
  transverseMomentum_ = t.transverseMomentum_;
  tags_ = t.tags_;
}

Track& Track::operator= (const Track &t) {
  // check for self-assignment by comparing the address of the
  // implicit object and the parameter
  if (this == &t)
    return *this;
  
  // do the copy
  theta_ = t.theta_;
  phi_ = t.phi_;
  cotgTheta_ = t.cotgTheta_;
  eta_ = t.eta_;
//  correlations_.ResizeTo(t.correlations_);
//  correlations_ = t.correlations_;
//  covariances_.ResizeTo(t.covariances_);
//  covariances_ = t.covariances_;
//  correlationsRZ_.ResizeTo(t.correlationsRZ_);
//  correlationsRZ_ = t.correlationsRZ_;
//  covariancesRZ_.ResizeTo(t.covariancesRZ_);
//  covariancesRZ_ = t.covariancesRZ_;
//  deltarho_ = t.deltarho_;
//  deltaphi_ = t.deltaphi_;
//  deltad_ = t.deltad_;
//  deltaCtgTheta_ = t.deltaCtgTheta_;
//  deltaZ0_ = t.deltaZ0_;
//  deltaP_ = t.deltaP_;
  vector<Hit*>::const_iterator iter, guard = t.hitV_.end();
  for (iter = t.hitV_.begin(); iter != guard; iter++) {
    Hit* h = new Hit(*(*iter));
    addHit(h);
  }
  transverseMomentum_ = t.transverseMomentum_;
  tags_ = t.tags_;
 
  // return the existing object
  return *this;
}

/**
 * Gives the number of active hits
 * @param usePixels take into account also pixel hits
 * @return how many active hits there are in a track
 */
int Track::nActiveHits (bool usePixels /* = false */, bool useIP /* = true */ ) const {
  std::vector<Hit*>::const_iterator hitIt;
  Hit* myHit;
  int result=0;
  for (hitIt=hitV_.begin();
       hitIt!=hitV_.end();
       ++hitIt) {
    myHit=(*hitIt);
    if (myHit) {
      if ((useIP) || (!myHit->isIP())) {
	if ( (usePixels) || (!myHit->isPixel()) ) {
	  if (myHit->getObjectKind()==Hit::Active)
	    result++;
	}
      }
    }
  }
  return result;
}


/**
 * Gives the probabilty of having "clean" hits
 * for nuclear-interacting particles
 * @param usePixels take into account also pixel hits
 * @return a vector with the probabilities of hits
 */
std::vector<double> Track::hadronActiveHitsProbability(bool usePixels /*= false */) {
  std::vector<Hit*>::iterator hitIt;
  std::vector<double> result;
  double probability=1;
  Hit* myHit;
  RILength myMaterial;
  sort();
  // int debugCount = 0; // debug
  for (hitIt=hitV_.begin();
       hitIt!=hitV_.end();
       ++hitIt) {
    myHit=(*hitIt);
    if (myHit) {
      if ( (usePixels) || (!myHit->isPixel()) ) {
	if (myHit->getObjectKind()==Hit::Active) {
	  result.push_back(probability);
	}
      }
      // DEBUG:
      // std::cerr << "Hit " << debugCount++ 
      // << ((myHit->getObjectKind()==Hit::Active) ? "Active" : "Inactive")
      // << " probability = " << probability << endl;

      // Decrease the probability that the
      // next hit is a clean one
      myMaterial = myHit->getCorrectedMaterial();
      probability /= exp(myMaterial.interaction);
    }
  }
  return result;
}

/**
 * Gives the probability of having a given number of "clean" hits
 * for nuclear-interacting particles
 * @param nHits the required number of clean hits
 * @param usePixels take into account also pixel hits
 * @return a vector with the probabilities of hits
 */
double Track::hadronActiveHitsProbability(int nHits, bool usePixels /* = false */ ) {
  std::vector<Hit*>::iterator hitIt;
  double probability=1;
  Hit* myHit;
  RILength myMaterial;
  int goodHits=0;
  sort();
  for (hitIt=hitV_.begin();
       hitIt!=hitV_.end();
       ++hitIt) {
    myHit=(*hitIt);
    if (myHit) {
      if ( (usePixels) || (!myHit->isPixel()) ) {
	if (myHit->getObjectKind()==Hit::Active)
	  goodHits++;
      }
      // If I reached the requested number of hits
      if (goodHits==nHits) 
	return probability;
      // Decrease the probability that the
      // next hit is a clean one
      myMaterial = myHit->getCorrectedMaterial();
      probability /= exp(myMaterial.interaction);
    }
  }
  // If I did not reach the requestd number of active hits
  // The probability is zero
  return 0;
}



/**
 * Modifies the hits to remove the material
 */
void Track::removeMaterial() {
  std::vector<Hit*>::iterator it;
  RILength nullMaterial;
  for (it = hitV_.begin(); it!=hitV_.end(); ++it) {
    (*it)->setCorrectedMaterial(nullMaterial);
  }
}

/**
 * The destructor makes sure that the hit vector is cleaned up properly.
 */
Track::~Track() {
    std::vector<Hit*>::iterator hitIt;
    for (hitIt=hitV_.begin(); hitIt!=hitV_.end(); hitIt++) {
        if ((*hitIt)!=NULL) {
            delete (*hitIt);
        }
    }
    hitV_.clear();
}

/**
 * Setter for the track azimuthal angle.
 * @param newTheta A reference to the value of the angle from the z-axis of the track
 */
double Track::setTheta(double& newTheta) {
    theta_ = newTheta;
    cotgTheta_ = 1/tan(newTheta);
    eta_ = -log(tan(theta_/2));
    std::vector<Hit*>::iterator iter, guard = hitV_.end();
    for (iter = hitV_.begin(); iter != guard; iter++) (*iter)->updateRadius();
    return theta_;
};

/**
 * Setter for the track polar angle.
 * @param newTheta A reference to the value of the angle from the z-axis of the track
 */
double Track::setPhi(double& newPhi) {
    phi_ = newPhi;
    //std::vector<Hit*>::iterator iter, guard = hitV_.end();
    //for (iter = hitV_.begin(); iter != guard; iter++) (*iter)->updateRadius();
    return phi_;
};


/**
 * Adds a new hit to the track
 * @param newHit a pointer to the new hit to be added
 */
// TODO: maybe updateradius is not necessary here. To be checked
Hit* Track::addHit(Hit* newHit) {
  hitV_.push_back(newHit);
  if (newHit->getHitModule() != NULL) {
    tags_.insert(newHit->getHitModule()->trackingTags.begin(), newHit->getHitModule()->trackingTags.end()); 
  }
  newHit->setTrack(this); 
  newHit->updateRadius(); 
  return newHit;
}

/**
 * This function sorts the hits in the internal vector by their distance to the z-axis.
 */
void Track::sort() {
    std::stable_sort(hitV_.begin(), hitV_.end(), sortSmallerR);
}

void Track::assignTrackingVolumesToHits() {
  sort();

  double firstActiveHitPixelDistance, firstActiveHitOuterDistance, lastActiveHitPixelDistance, lastActiveHitOuterDistance;
  auto firstActiveHitPixel = std::find_if(hitV_.begin(), hitV_.end(), [&](Hit* hit) { return (hit->isPixel() && hit->getObjectKind()==Hit::Active); } );
  if (firstActiveHitPixel != hitV_.end()) firstActiveHitPixelDistance = (*firstActiveHitPixel)->getDistance();
  auto firstActiveHitOuter = std::find_if(hitV_.begin(), hitV_.end(), [&](Hit* hit) { return (!hit->isPixel() && hit->getObjectKind()==Hit::Active); } );
  if (firstActiveHitOuter != hitV_.end()) firstActiveHitOuterDistance = (*firstActiveHitOuter)->getDistance();

  auto lastActiveHitPixel = std::find_if(hitV_.rbegin(), hitV_.rend(), [&](Hit* hit) { return (hit->isPixel() && hit->getObjectKind()==Hit::Active); } );
  if (lastActiveHitPixel != hitV_.rend()) lastActiveHitPixelDistance = (*lastActiveHitPixel)->getDistance();
  auto lastActiveHitOuter = std::find_if(hitV_.rbegin(), hitV_.rend(), [&](Hit* hit) { return (!hit->isPixel() && hit->getObjectKind()==Hit::Active); } );
  if (lastActiveHitOuter != hitV_.rend()) lastActiveHitOuterDistance = (*lastActiveHitOuter)->getDistance();


  //std::vector<Hit*>::reverse_iterator lastActiveHitTotal = std::find_if(hitV_.rbegin(), hitV_.rend(), 
  //[&](Hit* hit) { return (hit->getObjectKind()==Hit::Active); }
  //);
  //for (auto it = lastActiveHitTotal; it != hitV_.rend(); ++it) {
  // (*it)->setTotalTrackingVolume(true);
  //}


  for (auto& it : hitV_) {
    double distance = it->getDistance();
    if (distance < firstActiveHitPixelDistance) it->setPixelIntersticeVolume(true);
    if (distance >= firstActiveHitPixelDistance && distance <= lastActiveHitPixelDistance) it->setPixelTrackingVolume(true);
    if (distance > lastActiveHitPixelDistance && distance < firstActiveHitOuterDistance) it->setIntersticeVolume(true);
    if (distance >= firstActiveHitOuterDistance && distance <= lastActiveHitOuterDistance) it->setOuterTrackingVolume(true);
    if (distance <= lastActiveHitPixelDistance || distance <= lastActiveHitOuterDistance) it->setTotalTrackingVolume(true);
  }

  
 

}



/*bool Track::isHitInTrackingVolume(double distance) {
  sort();

std::vector<Hit*>::reverse_iterator lastActiveHitTotal = std::find_if(hitV_.rbegin(), hitV_.rend(), 
									[&](Hit* hit) { return (hit->getObjectKind()==Hit::Active); }
									);

									}*/



/**
 * Compute the correlation matrices of the track hits for a series of different energies.
 * @param momenta A reference of the list of energies that the correlation matrices should be calculated for
 */
//void Track::computeCorrelationMatrix() {
//
//  // matrix size
//  int n = hitV_.size();
//  correlations_.ResizeTo(n,n);
//
//  // pre-compute the squares of the scattering angles
//  std::vector<double> thetasq;
//  // pre-fetch the error on ctg(theta)
//  // will be zero, if not known
//  double deltaCtgT = deltaCtgTheta_;
//
//  // precompute the curvature in mm^-1
//  double rho = 1E-3 * SimParms::getInstance().magField() * 0.3 / transverseMomentum_;
//  for (int i = 0; i < n - 1; i++) {
//    double th = hitV_.at(i)->getCorrectedMaterial().radiation;
//    //#ifdef HIT_DEBUG
//    //	    std::cerr << "material (" << i << ") = " << th << "\t at r=" << hitV_.at(i)->getRadius() << std::endl;
//    //#endif
//    if (th>0)
//      th = (13.6 * 13.6) / (1000 * 1000 * transverseMomentum_ * transverseMomentum_) * th * (1 + 0.038 * log(th)) * (1 + 0.038 * log(th));
//    else
//      th = 0;
//    thetasq.push_back(th);
//  }
//  // correlations: c is column, r is row
//  for (int c = 0; c < n; c++) {
//    // dummy value for correlations involving inactive surfaces
//    if (hitV_.at(c)->getObjectKind() == Hit::Inactive) {
//      for (int r = 0; r <= c; r++) correlations_(r, c) = 0.0;
//    }
//    // one of the correlation factors refers to an active surface
//    else {
//      for (int r = 0; r <= c; r++) {
//        // dummy value for correlation involving an inactive surface
//        if (hitV_.at(r)->getObjectKind() == Hit::Inactive) correlations_(r, c) = 0.0;
//        // correlations between two active surfaces
//        else {
//          double sum = 0.0;
//          for (int i = 0; i < r; i++)
//            sum = sum + (hitV_.at(c)->getRadius() - hitV_.at(i)->getRadius()) * (hitV_.at(r)->getRadius() - hitV_.at(i)->getRadius()) * thetasq.at(i);
//          if (r == c) {
//            double prec = hitV_.at(r)->getResolutionRphi(pt2radius(transverseMomentum_, SimParms::getInstance().magField())); // if Bmod = getResoX natural
//            sum = sum + prec * prec;
//          }
//          correlations_(r, c) = sum;
//          if (r != c) correlations_(c, r) = sum;
//        }
//      }
//    }
//  }
//
//  // remove zero rows and columns
//  int ia = -1;
//  bool look_for_active = false;
//  for (int i = 0; i < n; i++) {
//    if ((hitV_.at(i)->getObjectKind() == Hit::Inactive) && (!look_for_active)) {
//      ia = i;
//      look_for_active = true;
//    }
//    else if ((hitV_.at(i)->getObjectKind() == Hit::Active) && (look_for_active)) {
//      for (int j = 0; j < n; j++) {
//        correlations_(ia, j) = correlations_(i, j);
//        correlations_(j, ia) = correlations_(j, i);
//      }
//      correlations_(ia, ia) = correlations_(i, i);
//      ia++;
//    }
//  }
//  // resize matrix if necessary
//  if (ia != -1) correlations_.ResizeTo(ia, ia);
//
//  // check if matrix is sane and worth keeping
//  if (!((correlations_.GetNoElements() > 0) && (correlations_.Determinant() != 0.0))) {
//    logERROR(Form("A singular matrix was found (this is unexpected: all analyzed tracks should have >= 3 hits). nElements=%d, determinant = %f", correlations_.GetNoElements(), correlations_.Determinant()));
//  }
//}

/**
 * Compute the covariance matrices of the track hits from a series of previously calculated correlation matrices.
 * @param A reference to the map of correlation matrices per energy  that serves as the value source for the computation
 */
//void Track::computeCovarianceMatrix() {
//  unsigned int offset = 0;
//  unsigned int nhits = hitV_.size();
//  int n = correlations_.GetNrows();
//  TMatrixT<double> C(correlations_); // Local copy to be inverted
//  TMatrixT<double> diffsT(3, n);
//  TMatrixT<double> diffs(n, 3);
//  covariances_.ResizeTo(3, 3);
//
//  // set up partial derivative matrices diffs and diffsT
//  for (unsigned int i = 0; i < nhits; i++) {
//    if (hitV_.at(i)->getObjectKind()  == Hit::Active) {
//      diffs(i - offset, 0) = 0.5 * hitV_.at(i)->getRadius() * hitV_.at(i)->getRadius();
//      diffs(i - offset, 1) = - hitV_.at(i)->getRadius();
//      diffs(i - offset, 2) = 1;
//    }
//    else offset++;
//  }
//  diffsT.Transpose(diffs);
//  covariances_ = diffsT * C.Invert() * diffs;
//}

void Track::computeLocalResolution() {
  int n = hitV_.size();
  for (int i = 0; i < n; i++) {
    if (hitV_.at(i)->getObjectKind() != Hit::Inactive) {
      hitV_.at(i)->computeLocalResolution();
      //std::cout << hitV_.at(i)->getResolutionLocalX() << std::endl;
    }
  }
}

/**
 * Compute the correlation matrices of the track hits for a series of different energies.
 * @param momenta A reference of the list of energies that the correlation matrices should be calculated for
 */
//void Track::computeCorrelationMatrixRZ() {
//
//  // matrix size
//  int n = hitV_.size();
//  double ctgTheta = 1/tan(theta_);
//  correlationsRZ_.ResizeTo(n,n);
//
//  // set up correlation matrix
//  double curvatureR = pt2radius(transverseMomentum_, SimParms::getInstance().magField());
//  // pre-compute the squares of the scattering angles
//  // already divided by sin^2 (that is : we should use p instead of p_T here
//  // but the result for theta^2 differ by a factor 1/sin^2, which is exactly the
//  // needed factor to project the scattering angle on an horizontal surface
//  std::vector<double> thetaOverSin_sq;
//  for (int i = 0; i < n - 1; i++) {
//    double th = hitV_.at(i)->getCorrectedMaterial().radiation;
//    if (th>0)
//      // equivalent to p=transverseMomentum_/sin(theta_); and then computing th/sin(theta)/sin(theta) using p in place of p_T
//      th = (13.6 * 13.6) / (1000 * 1000 * transverseMomentum_ * transverseMomentum_ ) * th * (1 + 0.038 * log(th)) * (1 + 0.038 * log(th));
//    else
//      th = 0;
//    thetaOverSin_sq.push_back(th);
//  }
//  // correlations: c is column, r is row
//  for (int c = 0; c < n; c++) {
//      // dummy value for correlations involving inactive surfaces
//    if (hitV_.at(c)->getObjectKind() == Hit::Inactive) {
//      for (int r = 0; r <= c; r++) correlationsRZ_(r, c) = 0.0;
//    }
//    // one of the correlation factors refers to an active surface
//    else {
//      for (int r = 0; r <= c; r++) {
//        // dummy value for correlation involving an inactive surface
//        if (hitV_.at(r)->getObjectKind() == Hit::Inactive) correlationsRZ_(r, c) = 0.0;
//        // correlations between two active surfaces
//        else {
//          double sum = 0.0;
//          for (int i = 0; i < r; i++)
//            sum += thetaOverSin_sq.at(i)
//              * (hitV_.at(c)->getDistance() - hitV_.at(i)->getDistance())
//              * (hitV_.at(r)->getDistance() - hitV_.at(i)->getDistance());
//          if (r == c) {
//            double prec = hitV_.at(r)->getResolutionZ(curvatureR);
//            sum = sum + prec * prec;
//          }
//          correlationsRZ_(r, c) = sum;
//          if (r != c) correlationsRZ_(c, r) = sum;
//#undef CORRELATIONS_OFF_DEBUG
//#ifdef CORRELATIONS_OFF_DEBUG
//          if (r!=c) {
//            correlationsRZ_(c, r)=0;
//            correlationsRZ_(r, c)=0;
//          }
//#endif
//        }
//      }
//    }
//  }
//  // remove zero rows and columns
//  int ia = -1;
//  bool look_for_active = false;
//  for (int i = 0; i < n; i++) {
//    if ((hitV_.at(i)->getObjectKind() == Hit::Inactive) && (!look_for_active)) {
//      ia = i;
//      look_for_active = true;
//    }
//    else if ((hitV_.at(i)->getObjectKind() == Hit::Active) && (look_for_active)) {
//      for (int j = 0; j < n; j++) {
//        correlationsRZ_(ia, j) = correlationsRZ_(i, j);
//        correlationsRZ_(j, ia) = correlationsRZ_(j, i);
//      }
//      correlationsRZ_(ia, ia) = correlationsRZ_(i, i);
//      ia++;
//    }
//  }
//  // resize matrix if necessary
//  if (ia != -1) correlationsRZ_.ResizeTo(ia, ia);
//
//  // check if matrix is sane and worth keeping
//  if (!((correlationsRZ_.GetNoElements() > 0) && (correlationsRZ_.Determinant() != 0.0))) {
//    std::cerr << "WARNING: this should be handled properly" << std::endl;
//  }
//}

/**
 * Compute the covariance matrices of the track hits from a series of previously calculated correlation matrices.
 * @param A reference to the map of correlation matrices per energy  that serves as the value source for the computation
 */
//void Track::computeCovarianceMatrixRZ() {
//  unsigned int offset = 0;
//  unsigned int nhits = hitV_.size();
//  int n = correlationsRZ_.GetNrows();
//  TMatrixT<double> C(correlationsRZ_); // Local copy to be inverted
//  TMatrixT<double> diffsT(2, n);
//  TMatrixT<double> diffs(n, 2);
//  covariancesRZ_.ResizeTo(2,2);
//
//  // set up partial derivative matrices diffs and diffsT
//  for (unsigned int i = 0; i < nhits; i++) {
//    if (hitV_.at(i)->getObjectKind()  == Hit::Active) {
//      // partial derivatives for x = p[0] * y + p[1]
//      diffs(i - offset, 0) = hitV_.at(i)->getRadius();
//      diffs(i - offset, 1) = 1;
//    }
//    else offset++;
//    }
//  diffsT.Transpose(diffs);
//  // Invert the C matrix
//  // TODO: check if this matrix can be inverted
//  C.Invert();
//  // compute covariancesRZ_ from diffsT, the correlation matrix and diffs
//  covariancesRZ_ = diffsT * C * diffs;
//}


/**
 * Calculate the errors of the track curvature radius, the propagation direction at the point of closest approach and the
 * distance of closest approach to the origin, all of them for each momentum of the test particle.
 * @param momentaList A reference of the list of energies that the errors should be calculated for
 */
//void Track::computeErrors() {
//  deltarho_ = 0 ;
//  deltaphi_ = 0 ;
//  deltad_ = 0 ;
//  deltaCtgTheta_ = 0 ;
//  deltaZ0_ = 0 ;
//  deltaP_ = 0 ;
//
//
//  // Compute spatial resolution for all active hits
//  computeLocalResolution();
//
//  // Compute the relevant matrices (RZ plane)
//  computeCorrelationMatrixRZ();
//  computeCovarianceMatrixRZ();
//  TMatrixT<double> dataRz(covariancesRZ_); // Local copy to be inverted
//  double err;
//  dataRz = dataRz.Invert();
//
//  if (dataRz(0, 0) >= 0) err = sqrt(dataRz(0, 0));
//  else err = -1;
//  deltaCtgTheta_ = err;
//
//  if (dataRz(1, 1) >= 0) err = sqrt(dataRz(1, 1));
//  else err = -1;
//  deltaZ0_ = err;
//
//  // rPhi plane
//  computeCorrelationMatrix();
//  computeCovarianceMatrix();
//
//  // calculate delta rho, delta phi and delta d maps from covariances_ matrix
//  TMatrixT<double> data(covariances_);
//  data = data.Invert();
//  if (data(0, 0) >= 0) err = sqrt(data(0, 0));
//  else err = -1;
//  deltarho_ = err;
//  if (data(1, 1) >= 0) err = sqrt(data(1, 1));
//  else err = -1;
//  deltaphi_ = err;
//  if (data(2, 2)) err = sqrt(data(2, 2));
//  else err = -1;
//  deltad_ = err;
//
//  // Combining into p measurement
//  double ptErr = deltarho_;
//  double R = transverseMomentum_ / SimParms::getInstance().magField() / 0.3 * 1E3; // curvature radius in mm
//  ptErr *= R; // fractional dpT/pT = dRho / Rho = dRho * R
//  // dp/p = dp_t/p_t + A / (1+A^2) * dA // with A = ctg(theta)
//  // dp/p = dp_t/p_t + sin(theta)*cos(theta) //
//  // double A = 1 / tan(theta_);
//  // double pErr = ptErr + A / (1+A*A) * ctgThetaErr;
//  deltaP_ = ptErr + sin(theta_) * cos(theta_) * deltaCtgTheta_;
//}

/**
 * Print the values in the correlation and covariance matrices and the drho, dphi and dd vectors per momentum.
 */
//void Track::printErrors() {
//    std::cout << "Overview of track errors:" << std::endl;
//    std::cout << "Hit correlation matrix: " << std::endl;
//    correlations_.Print();
//    std::cout << "Covariance matrix: " << std::endl;
//    covariances_.Print();
//    std::cout << "Rho errors by momentum: " << deltarho_ << std::endl;
//    std::cout << "Phi errors by momentum: " << deltaphi_ << std::endl;
//    std::cout << "D errors by momentum: " << deltad_ << std::endl;
//}

void Track::print() {
  std::cout << "******************" << std::endl;
  std::cout << "Track eta=" << eta_ << std::endl;
  for (const auto& it:hitV_) {
    std::cout << "    Hit"
              << " r=" << it->getRadius()
              << " d=" << it->getDistance()
              << " rl=" << it->getCorrectedMaterial().radiation
              << " il=" << it->getCorrectedMaterial().interaction
              << " getObjectKind()=" << it->getObjectKind();
    if (it->getObjectKind()==Hit::Active) {
      std::cout << " activeHitType_=" << it->getActiveHitType();
    }
    std::cout << std::endl;
  }
}

/**
 * Changes some active hits into inactive
 * according to the efficiency 
 * @param efficiency the modules active fraction
 * @param alsoPixel true if the efficiency removal applies to the pixel hits also
 */
void Track::addEfficiency() {
  for (std::vector<Hit*>::iterator it = hitV_.begin(); it!=hitV_.end(); ++it) {
    if ((*it)->getObjectKind() == Hit::Active) {
      double efficiency = (*it)->getHitModule()->singleHitEfficiency();
      if (efficiency!=1) {
        if ((double(random())/RAND_MAX) > efficiency) { // This hit is LOST
          (*it)->setObjectKind(Hit::Inactive);
        }
      }
    }
  }
}

/**
 * Makes all non-trigger hits inactive
 */
void Track::keepTriggerOnly() {
  // int iRemove=0;
  for (std::vector<Hit*>::iterator it = hitV_.begin(); it!=hitV_.end(); ++it) {
    // if (debugRemoval) std::cerr << "Hit number "
    //	                           << iRemove++ << ": ";
    // if (debugRemoval) std::cerr << "r = " << (*it)->getRadius() << ", ";
    // if (debugRemoval) std::cerr << "d = " << (*it)->getDistance() << ", ";
    if ((*it)->getObjectKind() == Hit::Active) {
      // if (debugRemoval) std::cerr << "active ";
      if ((*it)->isPixel()) {
	// if (debugRemoval) std::cerr << "pixel: removed";
	(*it)->setObjectKind(Hit::Inactive);
      } else {
	Module* myModule = (*it)->getHitModule();
	if (myModule) {
	  // if (debugRemoval) std::cerr << "module ";
	  if (myModule->sensorLayout() != PT) {
	    // if (debugRemoval) std::cerr << "non-pt: removed";
	    (*it)->setObjectKind(Hit::Inactive);
	  } else {
	    // if (debugRemoval) std::cerr << "pt: kept";
	  }
	} else {
	  // if (debugRemoval) std::cerr << "active without module: kept";
	}
      }
    } else {
      // if (debugRemoval) std::cerr << "inactive";
    }
    // if (debugRemoval) std::cerr << std::endl;
  }

  // debugRemoval=false;
}


void Track::keepTaggedOnly(const string& tag) {
  for (auto h : hitV_) {
    Module* m = h->getHitModule();
    if (!m) continue;
    if (std::count_if(m->trackingTags.begin(), m->trackingTags.end(), [&tag](const string& s){ return s == tag; })) h->setObjectKind(Hit::Active);
    else h->setObjectKind(Hit::Inactive);
  }
}

/**
 * Sets all the hits to their trigger resolution
 */
void Track::setTriggerResolution(bool isTrigger) {
  Hit* myHit;
  for (std::vector<Hit*>::iterator it = hitV_.begin(); it!=hitV_.end(); ++it) {
    myHit = (*it);
    if (myHit->getObjectKind() == Hit::Active) {
      myHit->setTrigger(isTrigger);
    }
  }
}


/**
 * Adds the constraint of the IP in the form of a virtual module
 */
void Track::addIPConstraint(double dr, double dz) {
  // This modeling of the IP constraint waas validated:
  // By placing dr = 0.5 mm and dz = 1 mm one obtains
  // sigma(d0) = 0.5 mm and sigma(z0) = 1 mm
  Hit* newHit = new Hit(dr);
  newHit->setIP(true);
  RILength emptyMaterial;
  emptyMaterial.radiation = 0;
  emptyMaterial.interaction = 0;
  newHit->setPixel(false);
  newHit->setCorrectedMaterial(emptyMaterial);
  newHit->setOrientation(Hit::Horizontal);
  newHit->setObjectKind(Hit::Active);
  newHit->setResolutionRphi(dr);
  newHit->setResolutionY(dz);
  this->addHit(newHit);
}

RILength Track::getCorrectedMaterial() {
  std::vector<Hit*>::const_iterator hitIt;
  Hit* myHit;
  RILength result;
  result.radiation = 0;
  result.interaction = 0;
  for (hitIt=hitV_.begin();
       hitIt!=hitV_.end();
       ++hitIt) {
    myHit=(*hitIt);
    result += myHit->getCorrectedMaterial();
  }

  return result;
}

double Track::expectedTriggerPoints(const double& triggerMomentum) const {
  std::vector<Hit*>::const_iterator hitIt;
  Hit* myHit;
  double result=0;

  for (hitIt=hitV_.begin();
       hitIt!=hitV_.end();
       ++hitIt) {
    myHit=(*hitIt);
    if ((myHit) &&
	(myHit->isTrigger()) &&
	(!myHit->isIP()) &&
	(myHit->getObjectKind()==Hit::Active)) {
      // We've got a possible trigger here
      // Let's find the corresponding module
      Module* myModule = myHit->getHitModule();
      if (myModule) {
	result += PtErrorAdapter(*myModule).getTriggerProbability(triggerMomentum);
      } else {
	// Whoops: problem here: an active hit is not linked to any module
	std::cerr << "ERROR: this SHOULD NOT happen. in expectedTriggerPoints() an active hit does not correspond to any module!" << std::endl;
      }
    }
  }
  return result;
}


std::vector<std::pair<Module*, HitType>> Track::getHitModules() const {
  std::vector<Hit*>::const_iterator hitIt;
  Hit* myHit;
  std::vector<std::pair<Module*, HitType>> result;

  for (hitIt=hitV_.begin(); hitIt!=hitV_.end(); ++hitIt) {
    myHit=(*hitIt);
    if ((myHit) &&
        (myHit->isTrigger()) &&
        (!myHit->isIP()) &&
        (myHit->getObjectKind()==Hit::Active)) {
      // We've got a possible trigger here
      // Let's find the corresponding module
      Module* myModule = myHit->getHitModule();
      if (myModule) {
        result.push_back(std::make_pair(myModule, myHit->getActiveHitType()));
      } else {
        // Whoops: problem here: an active hit is not linked to any module
        std::cerr << "ERROR: this SHOULD NOT happen. in expectedTriggerPoints() an active hit does not correspond to any module!" << std::endl;
      }
    }
  }
  return result;
}

//void Track::setTransverseMomentum(const double newPt) {
//  transverseMomentum_ = newPt;
//}

//void Track::pruneHits() {
//  double R = transverseMomentum_ / SimParms::getInstance().magField() / 0.3 * 1E3; // curvature radius in mm
//  std::vector<Hit*> hitN;
//  for (auto hitIt=hitV_.begin(); hitIt!=hitV_.end(); ++hitIt) {
//    if (((*hitIt)->getRadius()) < 2*R) {  hitN.push_back(*hitIt); }
//  }
//  hitV_ = hitN;
//}

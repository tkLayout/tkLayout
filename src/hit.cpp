/**
 * @file hit.cpp
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "hit.hh"
#include "module.hh"
#include <global_constants.h>
#include <vector>
#include <algorithm>
#include <cstdlib>

using namespace ROOT::Math;
using namespace std;

// bool Track::debugRemoval = false; // debug
#ifdef HIT_DEBUG_RZ
bool Track::debugRZCovarianceMatrix = false;  // debug
bool Track::debugRZCorrelationMatrix = false;  // debug
bool Track::debugRZErrorPropagation = false;  // debug
#endif

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
    hitModule_ = NULL;
    orientation_ = Undefined;
    myTrack_ = NULL;
    isPixel_ = false;
    isTrigger_ = false;
    isIP_ = false;
    myResolutionRphi_ = 0;
    myResolutionY_ = 0;
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
    hitModule_ = h.hitModule_;
    correctedMaterial_ = h.correctedMaterial_;
    myTrack_ = NULL;
    isPixel_ = h.isPixel_;
    isTrigger_ = h.isTrigger_;
    isIP_ = h.isIP_;
    myResolutionRphi_ = h.myResolutionRphi_;
    myResolutionY_ = h.myResolutionY_;
}

/**
 * Constructor for a hit with no module at a given distance from the origin
 * @param myDistance distance from the origin
 */
Hit::Hit(double myDistance) {
    distance_ = myDistance;
    objectKind_ = Undefined;
    hitModule_ = NULL;
    orientation_ = Undefined;
    isTrigger_ = false;
    isPixel_ = false;
    isIP_ = false;
    myTrack_ = NULL;
}

/**
 * Constructor for a hit on a given module at a given distance from the origin
 * @param myDistance distance from the origin
 * @param myModule pointer to the module with the hit 
 */
Hit::Hit(double myDistance, Module* myModule) {
    distance_ = myDistance;
    objectKind_ = Active;
    orientation_ = Undefined; 
    isTrigger_ = false;
    isPixel_ = false;
    isIP_ = false;
    setHitModule(myModule);
    myTrack_ = NULL;
}


/*
 * Setter for the pointer to the active surface that caused the hit.
 * @param myModule A pointer to a barrel or endcap module; may be <i>NULL</i>
 */
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
Material Hit::getCorrectedMaterial() {
    return correctedMaterial_;
}

/**
 * Getter for the rPhi resolution (local x coordinate for a module)
 * If the hit is not active it returns -1
 * If the hit is connected to a module, then the module's resolution
 * is retured (if the hit is trigger-type, then then module's trigger resultion is requested)
 * if there is not any hit module, then the hit's resolution property is read and returned
 * @return the hit's local resolution
 */
double Hit::getResolutionRphi() {
  if (objectKind_!=Active) {
    std::cerr << "ERROR: Hit::getResolutionRphi called on a non-active hit" << std::endl;
    return -1;
  } else {
    if (hitModule_) {
      if (isTrigger_) return hitModule_->getResolutionRphiTrigger();
      else return hitModule_->getResolutionRphi();
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
double Hit::getResolutionY() {
  if (objectKind_!=Active) {
    std::cerr << "ERROR: Hit::getResolutionY called on a non-active hit" << std::endl;
    return -1;
  } else {
    if (hitModule_) {
      if (isTrigger_) return hitModule_->getResolutionYTrigger();
      else return hitModule_->getResolutionY();
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
    if (hitModule_->getSubdetectorType()==Module::Endcap) {
      //std::cout << " getSubdetectorType()==Endcap "; //debug
      if (hitModule_->getShape()==Module::Rectangular) {
       //std::cout << " getShape()==Rectangular "; //debug
       result = true;
      }
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
    EndcapModule* myECModule = NULL;
    myECModule = dynamic_cast<EndcapModule*>(hitModule_);
    if (myECModule) {
      //std::cout << " myECModule!= NULL "; //debug
      result = (myECModule->getWidthLo() + myECModule->getWidthHi()) / 2. / 2.;
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
    correlations_ = t.correlations_;
    covariances_ = t.covariances_;
    deltarho_ = t.deltarho_;
    deltaphi_ = t.deltaphi_;
    deltad_ = t.deltad_;
    deltaCtgTheta_ = t.deltaCtgTheta_;
    deltaZ0_ = t.deltaZ0_;
    deltaP_ = t.deltaP_;
    vector<Hit*>::const_iterator iter, guard = t.hitV_.end();
    for (iter = t.hitV_.begin(); iter != guard; iter++) {
        Hit* h = new Hit(*(*iter));
        addHit(h);
    }
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
  Material myMaterial;
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
  Material myMaterial;
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
  // If I did not reach the requested number of active hits
  // The probability is zero
  return 0;
}



/**
 * Modifies the hits to remove the material
 */
void Track::removeMaterial() {
  std::vector<Hit*>::iterator it;
  Material nullMaterial;
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
 * This function sorts the hits in the internal vector by their distance to the z-axis.
 */
void Track::sort() {
    std::stable_sort(hitV_.begin(), hitV_.end(), sortSmallerR);
}

/**
 * Compute the correlation matrices of the track hits for a series of different energies.
 * @param momenta A reference of the list of energies that the correlation matrices should be calculated for
 */
void Track::computeCorrelationMatrix(const vector<double>& momenta) {
    // reset map
    correlations_.clear();
    // matrix size
    int n = hitV_.size();

    
#ifdef HIT_DEBUG
    std::cerr << std::endl
	      << std::endl
	      << "=== Track::computeCorrelationMatrix() == " << std::endl
	      << " theta = " << theta_ << std::endl;
#endif
    


    // set up correlation matrix
    for (unsigned int p = 0; p < momenta.size(); p++) {

#ifdef HIT_DEBUG
      std::cerr << " p = " << momenta.at(p) << std::endl;
#endif

        TMatrixTSym<double> corr(n);
        // pre-compute the squares of the scattering angles
        std::vector<double> thetasq;
	// pre-fetch the error on ctg(theta)
	// will be zero, if not known
	double deltaCtgT = deltaCtgTheta_[momenta.at(p)];

	// precompute the curvature in mm^-1
	double rho = 1E-3 * insur::magnetic_field * 0.3 / momenta.at(p);
        for (int i = 0; i < n - 1; i++) {
            double th = hitV_.at(i)->getCorrectedMaterial().radiation;
#ifdef HIT_DEBUG
	    std::cerr << "material (" << i << ") = " << th << "\t at r=" << hitV_.at(i)->getRadius() << std::endl;
#endif
	    if (th>0)
	      th = (13.6 * 13.6) / (1000 * 1000 * momenta.at(p) * momenta.at(p)) * th * (1 + 0.038 * log(th)) * (1 + 0.038 * log(th));
	    else
	      th = 0;
            thetasq.push_back(th);
        }
        // correlations: c is column, r is row
        for (int c = 0; c < n; c++) {
            // dummy value for correlations involving inactive surfaces
            if (hitV_.at(c)->getObjectKind() == Hit::Inactive) {
                for (int r = 0; r <= c; r++) corr(r, c) = 0.0;
            }
            // one of the correlation factors refers to an active surface
            else {
                for (int r = 0; r <= c; r++) {
                    // dummy value for correlation involving an inactive surface
                    if (hitV_.at(r)->getObjectKind() == Hit::Inactive) corr(r, c) = 0.0;
                    // correlations between two active surfaces
                    else {
                        double sum = 0.0;
                        for (int i = 0; i < r; i++)
                            sum = sum + (hitV_.at(c)->getRadius() - hitV_.at(i)->getRadius()) * (hitV_.at(r)->getRadius() - hitV_.at(i)->getRadius()) * thetasq.at(i);
                        if (r == c) {
                            double prec = hitV_.at(r)->getResolutionRphi();
#ifdef HIT_DEBUG
			    std::cerr << "Hit precision: " << prec << std::endl;
			    std::cerr << "Radius: " << hitV_.at(r)->getRadius() << std::endl;
#endif
			    if (hitV_.at(r)->getOrientation()==Hit::Vertical) {
			      // I have to introduce an additional error in the position
			      // to account for the uncertainty on r

                              // TODO: IMPORTANT
                              // error on ctgTheta is correlated!

			      // The component due to ctgTheta is
			      double deltar_ctg = hitV_.at(c)->getRadius() * tan(theta_) * deltaCtgT;
			      // The intrinsic r measurement resolution is
			      double deltar_y =  hitV_.at(c)->getResolutionY();
                              // We must combine this information: we get a bit more precise
                              double deltar_tot_sq = 1 / (
							  (1/deltar_ctg/deltar_ctg)
							  + (1/deltar_y/deltar_y));
			      // This is equivalent to a (squared) rPhi error of
			      double delta_rPhi_sq = pow(rho * hitV_.at(c)->getRadius(),2) * deltar_tot_sq;

                              // If the module is square things get slightly worse:
                              // + [ 1/63 (d/r)^3 ] * Delta_r^2
                              if (hitV_.at(c)->isSquareEndcap()) {
                                  //std::cout << "p= " << momenta.at(p); 
                                  //std::cout << " r= " << hitV_.at(c)->getRadius();
                                  //std::cout << " z= " << hitV_.at(c)->getRadius()/tan(theta_);
                                  //std::cout << " before= " << sqrt(delta_rPhi_sq)*1000; 
                                  double d_over_r = hitV_.at(c)->getD() / hitV_.at(c)->getRadius(); // TODO: make this nicer IMPORTANT!
                                  delta_rPhi_sq += ( 1/63.*pow(d_over_r,3) ) * deltar_tot_sq;
                                  //std::cout << " add= " << sqrt(( 1/63.*pow(d_over_r,3) ) * deltar_tot_sq)*1000; 
                                  //std::cout << std::endl; 
                              }
			      
			      // Which finally composes to the actual r-Phi error as:
			      prec = sqrt(prec*prec+delta_rPhi_sq);
			      
			    }
                            sum = sum + prec * prec;
                        }
                        corr(r, c) = sum;
                        if (r != c) corr(c, r) = sum;
#undef CORRELATIONS_OFF_DEBUG
#ifdef CORRELATIONS_OFF_DEBUG
if (r!=c) {
corr(c, r)=0;
corr(r, c)=0;
}
#endif
                    }
                }
            }
        }

        // remove zero rows and columns
        int ia = -1;
        bool look_for_active = false;
        for (int i = 0; i < n; i++) {
            if ((hitV_.at(i)->getObjectKind() == Hit::Inactive) && (!look_for_active)) {
                ia = i;
                look_for_active = true;
            }
            else if ((hitV_.at(i)->getObjectKind() == Hit::Active) && (look_for_active)) {
                for (int j = 0; j < n; j++) {
                    corr(ia, j) = corr(i, j);
                    corr(j, ia) = corr(j, i);
                }
                corr(ia, ia) = corr(i, i);
                ia++;
            }
        }
        // resize matrix if necessary
        if (ia != -1) corr.ResizeTo(ia, ia);

#ifdef HIT_DEBUG
	std::cerr << "Correlation matrix: " << std::endl;
	corr.Print();
#endif

        // check if matrix is sane and worth keeping
        if ((corr.GetNoElements() > 0) && (corr.Determinant() != 0.0)) {
            pair<double, TMatrixTSym<double> > par(momenta.at(p), corr);
            correlations_.insert(par);
        }
    }
}

/**
 * Compute the covariance matrices of the track hits from a series of previously calculated correlation matrices.
 * @param A reference to the map of correlation matrices per energy  that serves as the value source for the computation
 */
void Track::computeCovarianceMatrix() {
    map<momentum, TMatrixTSym<double> >::const_iterator iter, guard = correlations_.end();
    covariances_.clear();

#ifdef HIT_DEBUG
    std::cerr << std::endl
	      << std::endl
	      << "=== Track::computeCovarianceMatrix() == " << std::endl
	      << " theta = " << theta_ << std::endl;
#endif

    for (iter = correlations_.begin(); iter != guard; iter++) {
        unsigned int offset = 0;
        unsigned int nhits = hitV_.size();
        int n = iter->second.GetNrows();
        TMatrixT<double> C(correlations_[iter->first]);
        TMatrixT<double> diffsT(3, n);
        TMatrixT<double> diffs(n, 3);
        TMatrixT<double> cov(3, 3);

        // set up partial derivative matrices diffs and diffsT
        for (unsigned int i = 0; i < nhits; i++) {
            if (hitV_.at(i)->getObjectKind()  == Hit::Active) {
                diffs(i - offset, 0) = 0.5 * hitV_.at(i)->getRadius() * hitV_.at(i)->getRadius();
                diffs(i - offset, 1) = - hitV_.at(i)->getRadius();
                diffs(i - offset, 2) = 1;
            }
            else offset++;
        }
        diffsT.Transpose(diffs);
#ifdef HIT_DEBUG
	diffs.Print();
#endif
        // compute cov from diffsT, the correlation matrix and diffs
        cov = diffsT * C.Invert() * diffs;
#ifdef HIT_DEBUG
	cov.Print();
#endif
        pair<momentum, TMatrixT<double> > par(iter->first, cov);
        covariances_.insert(par);
    }
}

/**
 * Compute the correlation matrices of the track hits for a series of different energies.
 * @param momenta A reference of the list of energies that the correlation matrices should be calculated for
 */
void Track::computeCorrelationMatrixRZ(const vector<double>& momenta) {
    // reset map
    correlationsRZ_.clear();
    // matrix size
    int n = hitV_.size();
    
#ifdef HIT_DEBUG_RZ
    if (debugRZCorrelationMatrix)
      std::cerr << std::endl
		<< std::endl
		<< "=== Track::computeCorrelationMatrixRZ() == " << std::endl
		<< " theta = " << theta_ << std::endl;
#endif

    double ctgTheta = 1/tan(theta_);
    
    // set up correlation matrix
    for (unsigned int p = 0; p < momenta.size(); p++) {

#ifdef HIT_DEBUG_RZ
    if (debugRZCorrelationMatrix)
      std::cerr << " p = " << momenta.at(p) << std::endl;
#endif

        TMatrixTSym<double> corr(n);
        // pre-compute the squares of the scattering angles
	// already divided by sin^2 (that is : we should use p instead of p_T here
	// but the result for theta^2 differ by a factor 1/sin^2, which is exactly the
	// needed factor to project the scattering angle on an horizontal surface
        std::vector<double> thetaOverSin_sq;
        for (int i = 0; i < n - 1; i++) {
            double th = hitV_.at(i)->getCorrectedMaterial().radiation;
#ifdef HIT_DEBUG_RZ
    if (debugRZCorrelationMatrix)
	    std::cerr << "material (" << i << ") = " << th
		      << "\t at r=" << hitV_.at(i)->getRadius()
		      << "\t " << ((hitV_.at(i)->getOrientation()==Hit::Horizontal) ? "Horizontal" : "Vertical")
		      << "\t " << ((hitV_.at(i)->getObjectKind()==Hit::Active) ? "Active" : "Inactive")
		      << std::endl;
#endif
	    if (th>0)
	      // equivalent to p=momenta.at(p)/sin(theta_); and then computing th/sin(theta)/sin(theta) using p in place of p_T
	      th = (13.6 * 13.6) / (1000 * 1000 * momenta.at(p) * momenta.at(p) ) * th * (1 + 0.038 * log(th)) * (1 + 0.038 * log(th));
	    else
	      th = 0;
            thetaOverSin_sq.push_back(th);
        }
        // correlations: c is column, r is row
        for (int c = 0; c < n; c++) {
            // dummy value for correlations involving inactive surfaces
            if (hitV_.at(c)->getObjectKind() == Hit::Inactive) {
                for (int r = 0; r <= c; r++) corr(r, c) = 0.0;
            }
            // one of the correlation factors refers to an active surface
            else {
                for (int r = 0; r <= c; r++) {
                    // dummy value for correlation involving an inactive surface
                    if (hitV_.at(r)->getObjectKind() == Hit::Inactive) corr(r, c) = 0.0;
                    // correlations between two active surfaces
                    else {
                        double sum = 0.0;
                        for (int i = 0; i < r; i++)
			  sum += thetaOverSin_sq.at(i)
			    * (hitV_.at(c)->getDistance() - hitV_.at(i)->getDistance())
			    * (hitV_.at(r)->getDistance() - hitV_.at(i)->getDistance());
                        if (r == c) {
			  double prec = hitV_.at(r)->getResolutionY();
			  if (hitV_.at(r)->getOrientation()==Hit::Vertical) prec *= ctgTheta;
#ifdef HIT_DEBUG_RZ
			    if (debugRZCorrelationMatrix) {
			      std::cerr << "Hit precision: " << prec << "\t";
			      std::cerr << "Distance: " << hitV_.at(r)->getDistance() << std::endl;
			    }
#endif
                            sum = sum + prec * prec;
                        }
                        corr(r, c) = sum;
                        if (r != c) corr(c, r) = sum;
#undef CORRELATIONS_OFF_DEBUG
#ifdef CORRELATIONS_OFF_DEBUG
if (r!=c) {
corr(c, r)=0;
corr(r, c)=0;
}
#endif
                    }
                }
            }
        }
        // remove zero rows and columns
        int ia = -1;
        bool look_for_active = false;
        for (int i = 0; i < n; i++) {
            if ((hitV_.at(i)->getObjectKind() == Hit::Inactive) && (!look_for_active)) {
                ia = i;
                look_for_active = true;
            }
            else if ((hitV_.at(i)->getObjectKind() == Hit::Active) && (look_for_active)) {
                for (int j = 0; j < n; j++) {
                    corr(ia, j) = corr(i, j);
                    corr(j, ia) = corr(j, i);
                }
                corr(ia, ia) = corr(i, i);
                ia++;
            }
        }
        // resize matrix if necessary
        if (ia != -1) corr.ResizeTo(ia, ia);

#ifdef HIT_DEBUG_RZ
	if (debugRZCorrelationMatrix) {
	  std::cerr << "Correlation matrix: " << std::endl;
	  corr.Print();
	  debugRZCorrelationMatrix = false;
	}
#endif

        // check if matrix is sane and worth keeping
        if ((corr.GetNoElements() > 0) && (corr.Determinant() != 0.0)) {
            pair<double, TMatrixTSym<double> > par(momenta.at(p), corr);
            correlationsRZ_.insert(par);
        }
    }

}

/**
 * Compute the covariance matrices of the track hits from a series of previously calculated correlation matrices.
 * @param A reference to the map of correlation matrices per energy  that serves as the value source for the computation
 */
void Track::computeCovarianceMatrixRZ() {
    map<momentum, TMatrixTSym<double> >::const_iterator iter, guard = correlationsRZ_.end();
    covariancesRZ_.clear();

#ifdef HIT_DEBUG_RZ
    if (debugRZCovarianceMatrix) {
      std::cerr << std::endl
		<< std::endl
		<< "=== Track::computeCovarianceMatrixRZ() == " << std::endl
		<< " theta = " << theta_ << std::endl;
    }
#endif

    for (iter = correlationsRZ_.begin(); iter != guard; iter++) {
        unsigned int offset = 0;
        unsigned int nhits = hitV_.size();
        int n = iter->second.GetNrows();
        TMatrixT<double> C(correlationsRZ_[iter->first]);
        TMatrixT<double> diffsT(2, n);
        TMatrixT<double> diffs(n, 2);
        TMatrixT<double> cov(2, 2);

        // set up partial derivative matrices diffs and diffsT
        for (unsigned int i = 0; i < nhits; i++) {
	  if (hitV_.at(i)->getObjectKind()  == Hit::Active) {
	    // partial derivatives for x = p[0] * y + p[1]
	    diffs(i - offset, 0) = hitV_.at(i)->getRadius();
	    diffs(i - offset, 1) = 1;
	  }
	  else offset++;
        }
        diffsT.Transpose(diffs);
#ifdef HIT_DEBUG_RZ
	if (debugRZCovarianceMatrix) {
	  std::cerr << "Partial derivatives matrix" << std::endl;
	  diffs.Print();
	  diffsT.Print();
	  std::cerr << "Error correlation matrix:" << std::endl;
	  C.Print();
	}
#endif
	// Invert the C matrix
	C.Invert();
#ifdef HIT_DEBUG_RZ
	if (debugRZCovarianceMatrix) {
	  std::cerr << "Error correlation matrix (inverted):" << std::endl;	  
	  C.Print();
	}
#endif
        // compute cov from diffsT, the correlation matrix and diffs
        cov = diffsT * C * diffs;
#ifdef HIT_DEBUG_RZ
	if (debugRZCovarianceMatrix) {
	  cov.Print();
	  debugRZCovarianceMatrix = false;
	}
#endif
        pair<momentum, TMatrixT<double> > par(iter->first, cov);
        covariancesRZ_.insert(par);
    }
}

/**
 * Calculate the errors of the track curvature radius, the propagation direction at the point of closest approach and the
 * distance of closest approach to the origin, all of them for each momentum of the test particle.
 * @param momentaList A reference of the list of energies that the errors should be calculated for
 */
void Track::computeErrors(const std::vector<momentum>& momentaList) {
    map<momentum, TMatrixT<double> >::const_iterator iter, guard;
    deltarho_.clear();
    deltaphi_.clear();
    deltad_.clear();
    deltaCtgTheta_.clear();
    deltaZ0_.clear();
    deltaP_.clear();

    // Compute the relevant matrices (RZ plane)
    computeCorrelationMatrixRZ(momentaList);
    computeCovarianceMatrixRZ();
    guard = covariancesRZ_.end();
    for (iter = covariancesRZ_.begin(); iter != guard; iter++) {
        TMatrixT<double> data(iter->second);
        pair<momentum, double> err;
        err.first = iter->first;
        data = data.Invert();
#ifdef HIT_DEBUG_RZ
	if (debugRZErrorPropagation) {
	  std::cerr << "Matrix S" << std::endl;
	  data.Print();
	}
#endif
        if (data(0, 0) >= 0) err.second = sqrt(data(0, 0));
        else err.second = -1;
        deltaCtgTheta_.insert(err);
#ifdef HIT_DEBUG_RZ
	if (debugRZErrorPropagation) {
	  std::cerr << err.second << std::endl;
	}
#endif
        if (data(1, 1) >= 0) err.second = sqrt(data(1, 1));
        else err.second = -1;
        deltaZ0_.insert(err);

#ifdef HIT_DEBUG_RZ
	if (debugRZErrorPropagation) {
	  std::cerr << err.second << std::endl;
	  debugRZErrorPropagation = false;
	}
#endif
    }

    // rPhi plane
    computeCorrelationMatrix(momentaList);
    computeCovarianceMatrix();
    // calculate delta rho, delta phi and delta d maps from covariances_ matrix
    guard = covariances_.end();
    for (iter = covariances_.begin(); iter != guard; iter++) {
        TMatrixT<double> data(iter->second);
        pair<momentum, double> err;
        err.first = iter->first;
        data = data.Invert();
#ifdef HIT_DEBUG
	std::cerr << "Matrix S" << std::endl;
	data.Print();
#endif
        if (data(0, 0) >= 0) err.second = sqrt(data(0, 0));
        else err.second = -1;
        deltarho_.insert(err);
#ifdef HIT_DEBUG
	std::cerr << err.second << std::endl;
#endif
        if (data(1, 1) >= 0) err.second = sqrt(data(1, 1));
        else err.second = -1;
        deltaphi_.insert(err);
        if (data(2, 2)) err.second = sqrt(data(2, 2));
        else err.second = -1;
        deltad_.insert(err);
    }

    // Combining into p measurement
    for (unsigned int i=0; i<momentaList.size(); ++i) {
      double ptErr = deltarho_[momentaList[i]];
      double R = momentaList[i] / insur::magnetic_field / 0.3 * 1E3; // radius in mm
      ptErr *= R; // fractional dpT/pT = dRho / Rho = dRho * R
      double ctgThetaErr = deltaCtgTheta_[momentaList[i]];
      // dp/p = dp_t/p_t + A / (1+A^2) * dA // with A = ctg(theta)
      // dp/p = dp_t/p_t + sin(theta)*cos(theta) //
      //double A = 1 / tan(theta_);
      //double pErr = ptErr + A / (1+A*A) * ctgThetaErr;
      double pErr = ptErr + sin(theta_) * cos(theta_) * ctgThetaErr;
      deltaP_[momentaList[i]]=pErr;
    }

}

/**
 * Print the values in the correlation and covariance matrices and the drho, dphi and dd vectors per momentum.
 */
void Track::printErrors() {
    std::map<momentum, double>::const_iterator iter, guard;
    std::map<momentum, TMatrixT<double> >::const_iterator miter, mguard;
    std::map<momentum, TMatrixTSym<double> >::const_iterator siter, sguard;
    std::cout << "Overview of track errors:" << std::endl;
    std::cout << "Hit correlation matrices: " << correlations_.size() << (correlations_.size() == 1 ? " entry." : " entries.") << std::endl;
    sguard = correlations_.end();
    for (siter = correlations_.begin(); siter != sguard; siter++) {
      std::cout << "momentum = " << siter->first << ":";
      siter->second.Print();
    }
    std::cout << "Covariance matrices: " << covariances_.size() << (covariances_.size() == 1 ? " entry." : " entries.") << std::endl; // debug
    mguard = covariances_.end();
    for (miter = covariances_.begin(); miter != mguard; miter++) {
      std::cout << "momentum = " << miter->first << ":";
      miter->second.Print();
    }
    std::cout << "Rho errors by momentum: " << deltarho_.size() << (deltarho_.size() == 1 ? " entry." : " entries.") << std::endl;
    guard = deltarho_.end();
    for (iter = deltarho_.begin(); iter != guard; iter++)
        std::cout << "momentum = " << iter->first << ": deltaRho = " << iter->second << std::endl;
    std::cout << "Phi errors by momentum: " << deltaphi_.size() << (deltaphi_.size() == 1 ? " entry." : " entries.") << std::endl;
    guard = deltaphi_.end();
    for (iter = deltaphi_.begin(); iter != guard; iter++)
        std::cout << "momentum = " << iter->first << ": deltaPhi = " << iter->second << std::endl;
    std::cout << "D errors by momentum: " << deltad_.size() << (deltad_.size() == 1 ? " entry." : " entries.") << std::endl;
    guard = deltad_.end();
    for (iter = deltad_.begin(); iter != guard; iter++)
        std::cout << "momentum = " << iter->first << ": deltaD = " << iter->second << std::endl;
}

/**
 * Changes some active hits into inactive
 * according to the efficiency 
 * @param efficiency the modules active fraction
 * @param alsoPixel true if the efficiency removal applies to the pixel hits also
 */
void Track::addEfficiency(double efficiency, bool pixel /* = false */ ) {
  for (std::vector<Hit*>::iterator it = hitV_.begin(); it!=hitV_.end(); ++it) {
    if ((*it)->getObjectKind() == Hit::Active) {
      if ((pixel)&&(*it)->isPixel()) {
	if ((double(random())/RAND_MAX) > efficiency) { // This hit is LOST
	  (*it)->setObjectKind(Hit::Inactive);
	}
      }
      if ((!pixel)&&(!(*it)->isPixel())) {
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
	  if (myModule->getReadoutType() != Module::Pt) {
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
  Material emptyMaterial;
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

Material Track::getCorrectedMaterial() {
  std::vector<Hit*>::const_iterator hitIt;
  Hit* myHit;
  Material result;
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
	result += myModule->getTriggerProbability(triggerMomentum);
      } else {
	// Whoops: problem here: an active hit is not linked to any module
	std::cerr << "ERROR: this SHOULD NOT happen. in expectedTriggerPoints() an active hit does not correspond to any module!" << std::endl;
      }
    }
  }
  return result;
}

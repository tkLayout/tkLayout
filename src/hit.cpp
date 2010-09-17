/**
 * @file hit.cpp
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "hit.hh"
#include "module.hh"
#include <vector>
#include <algorithm>

using namespace ROOT::Math;
using namespace std;

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
    //trackTheta_ = 0;
    myTrack_ = NULL;
    isPixel_ = false;
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
}

/**
 * //TODO
 */
Hit::Hit(double myDistance) {
    distance_ = myDistance;
    objectKind_ = Undefined;
    hitModule_ = NULL;
    orientation_ = Undefined;
    //trackTheta_ = 0;
    myTrack_ = NULL;
}

/**
 * //TODO
 */
Hit::Hit(double myDistance, Module* myModule) {
    distance_ = myDistance;
    objectKind_ = Active;
    orientation_ = Undefined;
    //trackTheta_ = 0;
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
pair<double, double> Hit::getCorrectedMaterial() {
    return correctedMaterial_;
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
    vector<Hit*>::const_iterator iter, guard = t.hitV_.end();
    for (iter = t.hitV_.begin(); iter != guard; iter++) {
        Hit* h = new Hit(*(*iter));
        addHit(h);
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
 * Setter for the track angle.
 * @param newTheta A reference to the value of the angle from the z-axis of the track
 */
double Track::setTheta(double& newTheta) {
    theta_ = newTheta;
    std::vector<Hit*>::iterator iter, guard = hitV_.end();
    for (iter = hitV_.begin(); iter != guard; iter++) (*iter)->updateRadius();
    return theta_;
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

    std::cerr << std::endl
	      << std::endl
	      << "=== Track::computeCorrelationMatrix() == " << std::endl
	      << " theta = " << theta_ << std::endl; // debug


    // set up correlation matrix
    for (unsigned int p = 0; p < momenta.size(); p++) {

        std::cerr << " p = " << momenta.at(p) << std::endl; // debug

        TMatrixTSym<double> corr(n);
        // pre-compute the squares of the scattering angles
        std::vector<double> thetasq;
        for (int i = 0; i < n - 1; i++) {
            double th = hitV_.at(i)->getCorrectedMaterial().first;
            th = (13.6 * 13.6) / (1000 * 1000 * momenta.at(p) * momenta.at(p)) * th * (1 + 0.038 * log(th)) * (1 + 0.038 * log(th));
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
                            double prec = hitV_.at(r)->getHitModule()->getPrecisionRho();
			    std::cerr << "Hit precision: " << prec << std::endl; // debug
                            sum = sum + prec * prec;
                        }
                        corr(r, c) = sum;
                        if (r != c) corr(c, r) = sum;
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

	std::cerr << "Correlation matrix: " << std::endl; // debug
	corr.Print(); // debug

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
void Track::computeCovarianceMatrix(const map<double, TMatrixTSym<double> >& correlations) {
    map<momentum, TMatrixTSym<double> >::const_iterator iter, guard = correlations.end();
    covariances_.clear();

    std::cerr << std::endl
	      << std::endl
	      << "=== Track::computeCovarianceMatrix() == " << std::endl
	      << " theta = " << theta_ << std::endl; // debug

    for (iter = correlations.begin(); iter != guard; iter++) {
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
	diffs.Print(); // debug
        // compute cov from diffsT, the correlation matrix and diffs
        cov = diffsT * C.Invert() * diffs;
	cov.Print(); // debug
        pair<momentum, TMatrixT<double> > par(iter->first, cov);
        covariances_.insert(par);
    }
}

/**
 * Calculate the errors of the track curvature radius, the propagation direction at the point of closest approach and the
 * distance of closest approach to the origin, all of them for each momentum of the test particle.
 * @param momentaList A reference of the list of energies that the errors should be calculated for
 */
void Track::computeErrors(const std::vector<momentum>& momentaList) {
    // preliminary work
    computeCorrelationMatrix(momentaList);
    computeCovarianceMatrix(correlations_);
    // calculate delta rho, delta phi and delta d maps from covariances_ matrix
    map<momentum, TMatrixT<double> >::const_iterator iter, guard = covariances_.end();
    deltarho_.clear();
    deltaphi_.clear();
    deltad_.clear();
    for (iter = covariances_.begin(); iter != guard; iter++) {
        TMatrixT<double> data(iter->second);
        pair<momentum, double> err;
        err.first = iter->first;
        data = data.Invert();
	std::cerr << "Matrix S" << std::endl; // debug
	data.Print(); // debug
        if (data(0, 0) >= 0) err.second = sqrt(data(0, 0));
        else err.second = -1;
        deltarho_.insert(err);
	std::cerr << err.second << std::endl;
        if (data(1, 1) >= 0) err.second = sqrt(data(1, 1));
        else err.second = -1;
        deltaphi_.insert(err);
        if (data(2, 2)) err.second = sqrt(data(2, 2));
        else err.second = -1;
        deltad_.insert(err);
    }
    
    //for (deltarhoIt_=deltarho_.begin(); deltarhoIt_!=deltarho_.end(); ++deltarhoIt_)
    //  std::cerr << "deltarho ( " << (*deltarhoIt_).first << ") = " << (*deltarhoIt_).second << ", ";
    //std::cerr << std::endl; //debug
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
    std::cout << "Covariance matrices: " << covariances_.size() << (covariances_.size() == 1 ? " entry." : " entries.") << std::endl;
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


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
    radius_ = 0;
    objectKind_ = Undefined;
    hitModule_ = NULL;
    orientation_ = Undefined;
    //trackTheta_ = 0;
    myTrack_ = NULL;
}

Hit::Hit(const Hit& h) {
    distance_ = h.distance_;
    radius_ = h.radius_;
    orientation_ = h.orientation_;
    objectKind_ = h.objectKind_;
    hitModule_ = h.hitModule_;
    correctedMaterial_ = h.correctedMaterial_;
    myTrack_ = NULL;
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

Track::~Track() {
    std::vector<Hit*>::iterator hitIt;
    for (hitIt=hitV_.begin(); hitIt!=hitV_.end(); hitIt++) {
        if ((*hitIt)!=NULL) {
            delete (*hitIt);
        }
    }
    hitV_.clear();
}

double Track::setTheta(double& newTheta) {
    theta_ = newTheta;
    std::vector<Hit*>::iterator iter, guard = hitV_.end();
    for (iter = hitV_.begin(); iter != guard; iter++) (*iter)->updateRadius();
    return theta_;
};

void Track::sort() {
    std::stable_sort(hitV_.begin(), hitV_.end(), sortSmallerR);
}

void Track::computeCorrelationMatrix(const vector<double>& momenta) {
    // reset map
    correlations_.clear();
    // matrix size
    int n = hitV_.size();
    // set up correlation matrix
    for (unsigned int p = 0; p < momenta.size(); p++) {
        TMatrixTSym<double> corr(n);
        // pre-compute the squares of the scattering angles
        std::vector<double> thetasq;
        for (int i = 0; i < n - 1; i++) {
            double th = hitV_.at(i)->getCorrectedMaterial().second;
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
                            double pitch = (hitV_.at(r)->getHitModule()->getLowPitch() + hitV_.at(r)->getHitModule()->getHighPitch()) / 2.0;
                            sum = sum + pitch * pitch / 12.0;
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
        // check if matrix is sane and worth keeping
        if ((corr.GetNoElements() > 0) && (corr.Determinant() != 0.0)) {
            pair<double, TMatrixTSym<double> > par(momenta.at(p), corr);
            correlations_.insert(par);
        }
    }
}

void Track::computeCovarianceMatrix(const map<double, TMatrixTSym<double> >& correlations) {
    map<momentum, TMatrixTSym<double> >::const_iterator iter, guard = correlations.end();
    covariances_.clear();
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
        // compute cov from diffsT, the correlation matrix and diffs
        cov = diffsT * C.Invert() * diffs;
        pair<momentum, TMatrixT<double> > par(iter->first, cov);
        covariances_.insert(par);
    }
}

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
        if (data(0, 0) >= 0) err.second = sqrt(data(0,0));
        else err.second = -1;
        deltarho_.insert(err);
        if (data(1, 1) >= 0) err.second = sqrt(data(1,1));
        else err.second = -1;
        deltaphi_.insert(err);
        if (data(2, 2)) err.second = sqrt(data(2,2));
        else err.second = -1;
        deltad_.insert(err);
    }
}

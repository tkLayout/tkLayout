/**
 * @file Track.cpp
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "Track.h"

#include <global_constants.h>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "Hit.h"
#include "MessageLogger.h"
#include "MaterialProperties.h"
#include "PtErrorAdapter.h"
#include "SimParms.h"
#include "Units.h"


using namespace ROOT::Math;
using namespace std;

//
// Track constructor -> need to use setter methods to set: 2 of these [theta, phi, eta, cot(theta)] & 2 of these [mag. field, transv. momentum, radius]
//
Track::Track() :
  m_theta(0),
  m_phi(0),
  m_cotgTheta(0),
  m_eta(0),
  m_pt(0),
  m_magField(0),
  m_radius(0),
  m_deltaRho(0),
  m_deltaPhi(0),
  m_deltaD0(0),
  m_deltaCtgTheta(0),
  m_deltaZ0(0),
  m_deltaPt(0),
  m_deltaP(0)
{}

//
// Track copy-constructor -> creates deep copy of hit vector
//
Track::Track(const Track& track) {

  m_theta        = track.m_theta;
  m_phi          = track.m_phi;
  m_cotgTheta    = track.m_cotgTheta;
  m_eta          = track.m_eta;
  m_pt           = track.m_pt;
  m_magField     = track.m_magField;
  m_radius       = track.m_radius;
  m_deltaRho     = track.m_deltaRho;
  m_deltaPhi     = track.m_deltaPhi;
  m_deltaD0      = track.m_deltaD0;
  m_deltaCtgTheta= track.m_deltaCtgTheta;
  m_deltaZ0      = track.m_deltaZ0;
  m_deltaPt      = track.m_deltaPt;
  m_deltaP       = track.m_deltaP;

  m_origin       = track.m_origin;
  m_direction    = track.m_direction;

  m_correlationsRPhi.ResizeTo(track.m_correlationsRPhi);
  m_correlationsRPhi = track.m_correlationsRPhi;
  m_covariancesRPhi.ResizeTo(track.m_covariancesRPhi);
  m_covariancesRPhi  = track.m_covariancesRPhi;
  m_correlationsRZ.ResizeTo(track.m_correlationsRZ);
  m_correlationsRZ   = track.m_correlationsRZ;
  m_covariancesRZ.ResizeTo(track.m_covariancesRZ);
  m_covariancesRZ    = track.m_covariancesRZ;

  for (auto iHit : track.m_hits) {
    Hit* hit = new Hit(*iHit);
    addHit(hit);
  }
  m_tags = track.m_tags;
}

//
// Assign operator with deep copy of hit vector
//
Track& Track::operator= (const Track& track) {

  // check for self-assignment by comparing the address of the
  // implicit object and the parameter
  if (this == &track) return *this;
  
  // Do the copy
  m_theta        = track.m_theta;
  m_phi          = track.m_phi;
  m_cotgTheta    = track.m_cotgTheta;
  m_eta          = track.m_eta;
  m_pt           = track.m_pt;
  m_radius       = track.m_radius;
  m_magField     = track.m_magField;
  m_deltaRho     = track.m_deltaRho;
  m_deltaPhi     = track.m_deltaPhi;
  m_deltaD0      = track.m_deltaD0;
  m_deltaCtgTheta= track.m_deltaCtgTheta;
  m_deltaZ0      = track.m_deltaZ0;
  m_deltaPt      = track.m_deltaPt;
  m_deltaP       = track.m_deltaP;

  m_origin       = track.m_origin;
  m_direction    = track.m_direction;

  m_correlationsRPhi.ResizeTo(track.m_correlationsRPhi);
  m_correlationsRPhi = track.m_correlationsRPhi;
  m_covariancesRPhi.ResizeTo(track.m_covariancesRPhi);
  m_covariancesRPhi  = track.m_covariancesRPhi;
  m_correlationsRZ.ResizeTo(track.m_correlationsRZ);
  m_correlationsRZ   = track.m_correlationsRZ;
  m_covariancesRZ.ResizeTo(track.m_covariancesRZ);
  m_covariancesRZ    = track.m_covariancesRZ;

  for (auto iHit : track.m_hits) {
    Hit* hit = new Hit(*iHit);
    addHit(hit);
  }
  m_tags = track.m_tags;

  // Return the existing object
  return *this;
}

//
// Destructor
//
Track::~Track() {

  for (auto iHit : m_hits) {

    if (iHit!=nullptr) delete iHit;
  }
  m_hits.clear();
}

//
// Main method calculating track parameters using Karimaki approach & parabolic approximation in R-Phi plane: 1/R, D0, phi parameters
// and using linear fit in s-Z plane: cotg(theta), Z0 parameters
//
void Track::computeErrors() {

  // Initialize first
  m_deltaRho      = 0;
  m_deltaPhi      = 0;
  m_deltaD0       = 0;
  m_deltaCtgTheta = 0;
  m_deltaZ0       = 0;
  m_deltaPt       = 0;
  m_deltaP        = 0;

  // Compute the relevant matrices in RZ plane first
  computeCorrelationMatrixRZ();
  computeCovarianceMatrixRZ();
  double err;

  // delta(cotgTheta)
  if (m_covariancesRZ(0, 0)>=0) err = sqrt(m_covariancesRZ(0, 0));
  else err = -1;
  m_deltaCtgTheta = err;

  // delta(Z0)
  if (m_covariancesRZ(1, 1)>=0) err = sqrt(m_covariancesRZ(1, 1));
  else err = -1;
  m_deltaZ0 = err;

  // Compute the relevant matrices in R-Phi plane
  computeCorrelationMatrixRPhi();
  computeCovarianceMatrixRPhi();

  // delta(1/R) & delta(pT)
  if (m_covariancesRPhi(0, 0)>=0) err = sqrt(m_covariancesRPhi(0, 0));
  else err = -1;
  m_deltaRho = err;
  if (m_deltaRho!=-1) m_deltaPt  = m_deltaRho * m_radius; // dpT/pT = dRho / Rho = dRho * R
  else                m_deltaPt  = -1;

  // delta(phi)
  if (m_covariancesRPhi(1, 1) >= 0) err = sqrt(m_covariancesRPhi(1, 1));
  else err = -1;
  m_deltaPhi = err;

  // delta(D0)
  if (m_covariancesRPhi(2, 2)) err = sqrt(m_covariancesRPhi(2, 2));
  else err = -1;
  m_deltaD0 = err;

  // Combining into p measurement
  // dp/p = dp_t/p_t + A / (1+A^2) * dA // with A = ctg(theta)
  // dp/p = dp_t/p_t + sin(theta)*cos(theta)
  if (m_deltaPt!=-1 && m_deltaCtgTheta!=-1) m_deltaP = sqrt(m_deltaPt*m_deltaPt + sin(m_theta)*sin(m_theta) * cos(m_theta)*cos(m_theta) * m_deltaCtgTheta*m_deltaCtgTheta);
}

//
// Adds a new hit to the track (hit radius automatically updated)
//
Hit* Track::addHit(Hit* newHit) {

  // Add new hit
  m_hits.push_back(newHit);

  // Add tracking tags
  if (newHit->getHitModule() != nullptr) {
    m_tags.insert(newHit->getHitModule()->trackingTags.begin(), newHit->getHitModule()->trackingTags.end());
  }
  newHit->setTrack(this);

  return newHit;
}

//
// Add IP constraint to the track, technically new hit is assigned: with no material and hit resolution in R-Phi as dr, in s-Z as dz
//
Hit* Track::addIPConstraint(double dr, double dz) {

  // This modeling of the IP constraint was validated:
  // By placing dr = 0.5 mm and dz = 1 mm one obtains
  // sigma(d0) = 0.5 mm and sigma(z0) = 1 mm
  Hit* newHit = new Hit(dr);
  newHit->setIP(true);

  RILength emptyMaterial;
  emptyMaterial.radiation   = 0;
  emptyMaterial.interaction = 0;

  newHit->setPixel(false);
  newHit->setCorrectedMaterial(emptyMaterial);
  newHit->setOrientation(HitOrientation::Horizontal);
  newHit->setObjectKind(HitKind::Active);
  newHit->setResolutionRphi(dr);
  newHit->setResolutionY(dz);
  this->addHit(newHit);

  return newHit;
}

//
// Sort internally all hits assigned to this track -> sorting algorithm based on hit radius (the smaller, the sooner)
//
void Track::sortHits() {std::stable_sort(m_hits.begin(), m_hits.end(), Hit::sortSmallerR); }

//
// Remove hits that don't follow the parabolic approximation used in tracking - TODO: still needs to be updated (not all approximations taken into account)
//
bool Track::pruneHits() {

  double helixRadius = getRadius();
  bool   isPruned    = false;

  std::vector<Hit*> newHits;
  for (auto iHit : m_hits) {

    if (iHit->getRadius()<2*helixRadius) newHits.push_back(iHit);
    else {

      // Clear memory
      delete iHit;
      iHit = nullptr;

      isPruned = true;
    }
  }

  m_hits = newHits;
  return isPruned;
}

//
// Set active only hits with the given tag
//
void Track::keepTaggedOnly(const string& tag) {

  for (auto iHit : m_hits) {

    const DetectorModule* module = iHit->getHitModule();
    if (!module) continue;

    if (std::count_if(module->trackingTags.begin(), module->trackingTags.end(), [&tag](const string& s){ return s == tag; })) iHit->setObjectKind(HitKind::Active);
    else iHit->setObjectKind(HitKind::Inactive);
  }
}

//
// Remove material from all assigned hits -> modify all hits such as they are without any material
//
void Track::removeMaterial() {

  // Material object with no material assigned
  RILength nullMaterial;

  // Reset all material assigned to hits
  for (auto iHit : m_hits) iHit->setCorrectedMaterial(nullMaterial);
}

//
// Helper method printing track covariance matrices in R-Phi
//
void Track::printErrors() {

  std::cout << "Overview of track errors:" << std::endl;
  std::cout << "Hit correlation matrix: "  << std::endl;
  m_correlationsRPhi.Print();

  std::cout << "Covariance matrix: " << std::endl;
  m_covariancesRPhi.Print();

  std::cout << "Rho errors by momentum: " << m_deltaRho << std::endl;
  std::cout << "Phi errors by momentum: " << m_deltaPhi << std::endl;
  std::cout << "D errors by momentum: "   << m_deltaD0  << std::endl;
}

//
// Helper method printing track hits
//
void Track::printHits() {

  std::cout << "******************" << std::endl;
  std::cout << "Track eta=" << m_eta << std::endl;

  for (const auto& it : m_hits) {
    std::cout << "    Hit"
              << " r="  << it->getRadius()
              << " d="  << it->getDistance()
              << " rl=" << it->getCorrectedMaterial().radiation
              << " il=" << it->getCorrectedMaterial().interaction
              << " getObjectKind()=" << static_cast<short>(it->getObjectKind());
    if (it->getObjectKind()==HitKind::Active) {
      std::cout << " activeHitType_=" << static_cast<short>(it->getActiveHitType());
    }
    std::cout << std::endl;
  }
}

//
// Set track polar angle - theta, azimuthal angle - phi, particle transverse momentum - pt
// (magnetic field obtained automatically from SimParms singleton class)Setter for the track azimuthal angle.
//
const Polar3DVector& Track::setThetaPhiPt(const double& newTheta, const double& newPhi, const double& newPt) {

  m_theta     = newTheta;
  m_cotgTheta = 1/tan(newTheta);
  m_eta       = -log(tan(m_theta/2));
  m_phi       = newPhi;
  m_pt        = newPt;
  m_magField  = SimParms::getInstance()->magneticField();
  m_radius    = m_pt / (0.3 * m_magField);

  m_direction.SetCoordinates(m_radius/sin(m_theta), m_theta, m_phi);

  if (m_magField<=0) {
    logERROR("Track::setThetaPhiPt -> Magnetic field not defined!!!");
    EXIT_FAILURE;
  }

  for (auto iHit : m_hits) iHit->updateRadius();

  return m_direction;
};

//
// Get number of active hits assigned to track for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
//
int Track::getNActiveHits (std::string tag, bool useIP /* = true */ ) const {

  // Result variable
  int nHits=0;

  for (auto iHit : m_hits) {
    if (iHit) {
      if ((useIP) || (!iHit->isIP())) {
        for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
          if ((tag==*it) && (iHit->getObjectKind()==HitKind::Active)) nHits++;
        }
      }
    }
  } // For

  return nHits;
}

//
// Get the probabilty of having "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
//
std::vector<double> Track::getHadronActiveHitsProbability(std::string tag) {

  // Result variable
  std::vector<double> probabilities;
  double probability = 1;

  // Sort hits first
  sortHits();

  for (auto iHit : m_hits) {
    if (iHit) {
      for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
         if ((tag==*it) && (iHit->getObjectKind()==HitKind::Active)) probabilities.push_back(probability);
      }

      // Decrease the probability that the next hit is a clean one
      RILength myMaterial = iHit->getCorrectedMaterial();
      probability /= exp(myMaterial.interaction);
    }
  } // For

  return probabilities;
}

//
// Get the probabilty of having a given number of "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
//
double Track::getHadronActiveHitsProbability(std::string tag, int nHits) {

  // Probability
  double probability = 1;

  // Number of clean hits
  int goodHits = 0;

  // Sort hits first
  sortHits();

  for (auto iHit : m_hits) {

    if (iHit) {
      for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
        if ((tag==*it) && (iHit->getObjectKind()==HitKind::Active)) goodHits++;
      }


      // If I reached the requested number of hits
      if (goodHits==nHits) return probability;

      // Decrease the probability that the
      // next hit is a clean one
      RILength myMaterial = iHit->getCorrectedMaterial();
      probability /= exp(myMaterial.interaction);
    }
  }

  // If I did not reach the requested number of active hits
  // The probability is zero
  return 0;
}

//
// Get track material
//
RILength Track::getMaterial() {

  RILength totalMaterial;
  totalMaterial.radiation   = 0;
  totalMaterial.interaction = 0;

  for (auto iHit : m_hits) totalMaterial += iHit->getCorrectedMaterial();

  return totalMaterial;
}

//
// Get a vector of pairs: Detector module & hit type
//
std::vector<std::pair<const DetectorModule*, HitType>> Track::getHitModules() const {

  std::vector<std::pair<const DetectorModule*, HitType>> result;

  for (auto iHit : m_hits) {

    if ((iHit) && (iHit->isTrigger()) && (!iHit->isIP()) && (iHit->getObjectKind()==HitKind::Active)) {

      // We've got a possible trigger here
      // Let's find the corresponding module
      const DetectorModule* myModule = iHit->getHitModule();
      if (myModule) result.push_back(std::make_pair(myModule, iHit->getActiveHitType()));
      else {
        // Whoops: problem here: an active hit is not linked to any module
        std::cerr << "ERROR: this SHOULD NOT happen. in expectedTriggerPoints() an active hit does not correspond to any module!" << std::endl;
      }
    }
  }
  return result;
}


//
// Compute the correlation matrix of the track parameters in R-Phi projection
//
void Track::computeCorrelationMatrixRPhi() {

  // Matrix size
  int n = m_hits.size();
  m_correlationsRPhi.ResizeTo(n,n);

  // Pre-fetch the error on ctg(theta) -> will use zero value, if not known
  double deltaCtgT = m_deltaCtgTheta;

  // Precompute the curvature
  double rho = getRho();

  // Get contributions from Multiple Couloumb scattering
  std::vector<double> msThetaSq;

  for (int i = 0; i < n - 1; i++) {

    // MS theta
    double msTheta = 0.0;

    // Material in terms of rad. lengths
    double XtoX0 = m_hits.at(i)->getCorrectedMaterial().radiation;
    //std::cout << std::fixed << std::setprecision(4) << "Material (" << i << ") = " << XtoX0 << "\t at r=" << m_hits.at(i)->getRadius() << "\t of type=" << m_hits.at(i)->getObjectKind() << std::endl;

    if (XtoX0>0) {
      msTheta = (13.6*Units::MeV * 13.6*Units::MeV) / (m_pt/Units::MeV * m_pt/Units::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));
    }
    else {
      msTheta = 0;
    }
    msThetaSq.push_back(msTheta);
  }

  // Correlations: c is column, r is row
  for (int c = 0; c < n; c++) {

    // Dummy value for correlations involving inactive surfaces
    if (m_hits.at(c)->getObjectKind() == HitKind::Inactive) {
      for (int r = 0; r <= c; r++) m_correlationsRPhi(r, c) = 0.0;
    }
    // One of the correlation factors refers to an active surface
    else {

      for (int r = 0; r <= c; r++) {
        // Dummy value for correlation involving an inactive surface
        if (m_hits.at(r)->getObjectKind() == HitKind::Inactive) m_correlationsRPhi(r, c) = 0.0;

        // Correlations between two active surfaces
        else {

          double sum = 0.0;

          for (int i = 0; i < r; i++) {

            sum = sum + (m_hits.at(c)->getRadius() - m_hits.at(i)->getRadius()) * (m_hits.at(r)->getRadius() - m_hits.at(i)->getRadius()) * msThetaSq.at(i);
          //std::cout << ">> " << std::fixed << std::setprecision(4) << i << " " << c << " : " << m_hits.at(c)->getRadius()
          //            << " " << m_hist.at(i)->getRadius() << " " << msThetaSq.at(i) << " " << sum << std::endl;
          }
          if (r == c) {

            double prec = m_hits.at(r)->getResolutionRphi(m_radius); // if Bmod = getResoX natural
            sum = sum + prec * prec;
            //std::cout << ">>> " << sum << std::endl;
          }
          m_correlationsRPhi(r, c) = sum;
          if (r != c) m_correlationsRPhi(c, r) = sum;
        }
      }
    }
  } // Correlations: c is column, r is row

  // Remove zero rows and columns
  int  nResized        = -1;
  bool look_for_active = false;

  for (int i = 0; i < n; i++) {

    if ((m_hits.at(i)->getObjectKind() == HitKind::Inactive) && (!look_for_active)) {
      nResized = i;
      look_for_active = true;
    }
    else if ((m_hits.at(i)->getObjectKind() == HitKind::Active) && (look_for_active)) {

      for (int j = 0; j < n; j++) {
        m_correlationsRPhi(nResized, j) = m_correlationsRPhi(i, j);
        m_correlationsRPhi(j, nResized) = m_correlationsRPhi(j, i);
      }
      m_correlationsRPhi(nResized, nResized) = m_correlationsRPhi(i, i);
      nResized++;
    }
  }

  // Resize matrix if necessary
  if (nResized != -1) m_correlationsRPhi.ResizeTo(nResized, nResized);
//  std::cout << std::endl;
//  for (int i = 0; i<nResized; i++) {
//    std::cout << "(";
//    for (int j=0; j<nResized;j++) {
//
//      std::cout << " " << std::fixed << std::setprecision(4) << m_correlationsRPhi(i,j);
//    }
//    std::cout << ")" << std::endl;
//  }
//  std::cout << std::endl;

  // Check if matrix is sane and worth keeping
  if (!((m_correlationsRPhi.GetNoElements() > 0) && (m_correlationsRPhi.Determinant() != 0.0))) {
    logWARNING("Correlation matrix in R-Phi -> zero determinat or zero number of elements");
  }
}

//
// Compute the covariance matrix of the track parameters in R-Phi projection
//
void Track::computeCovarianceMatrixRPhi() {

  unsigned int offset = 0;
  unsigned int nHits  = m_hits.size();

  int n = m_correlationsRPhi.GetNrows();

  TMatrixT<double> C(m_correlationsRPhi); // Local copy to be inverted
  TMatrixT<double> diffsT(3, n);          // Derivatives of track parameters transposed (in R-Phi -> 3 track parameters)
  TMatrixT<double> diffs(n, 3);           // Derivatives of track parameters (in R-Phi -> 3 track parameters)

  m_covariancesRPhi.ResizeTo(3, 3);

  // Set up partial derivative matrices diffs and diffsT -> using Karimaki approach & parabolic aproximations to define these matrices
  for (auto i = 0; i < nHits; i++) {

    if (m_hits.at(i)->getObjectKind()  == HitKind::Active) {
      diffs(i - offset, 0) = 0.5 * m_hits.at(i)->getRadius() * m_hits.at(i)->getRadius();
      diffs(i - offset, 1) = - m_hits.at(i)->getRadius();
      diffs(i - offset, 2) = 1;
    }
    else offset++;
  }

  // Transpose
  diffsT.Transpose(diffs);

  // Get covariance matrix using global chi2 fit: cov(i,j) = (D^T * C^-1 * D)^-1
  m_covariancesRPhi = diffsT * C.Invert() * diffs;
  m_covariancesRPhi.Invert();
}

//
// Compute the correlation matrix of the track parameters in R-Z projection
//
void Track::computeCorrelationMatrixRZ() {

  // Matrix size
  int n = m_hits.size();
  m_correlationsRZ.ResizeTo(n,n);

  // Pre-compute the squares of the scattering angles
  // already divided by sin^2 (that is : we should use p instead of p_T here
  // but the result for theta^2 differ by a factor 1/sin^2, which is exactly the
  // needed factor to project the scattering angle on an horizontal surface
  std::vector<double> msThetaOverSinSq;

  for (int i = 0; i < n - 1; i++) {

    // MS theta
    double msTheta = 0.0;

    // Material in terms of rad. lengths
    double XtoX0 = m_hits.at(i)->getCorrectedMaterial().radiation;

    if (XtoX0>0) {
      // Equivalent to p=transverseMomentum_/sin(theta_); and then computing th/sin(theta)/sin(theta) using p in place of p_T
      msTheta = (13.6*Units::MeV * 13.6*Units::MeV) / (m_pt/Units::MeV * m_pt/Units::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));
    }
    else {
      msTheta = 0;
    }
    msThetaOverSinSq.push_back(msTheta);
  }

  // Correlations: c is column, r is row
  for (int c = 0; c < n; c++) {

    // Dummy value for correlations involving inactive surfaces
    if (m_hits.at(c)->getObjectKind() == HitKind::Inactive) {
      for (int r = 0; r <= c; r++) m_correlationsRZ(r, c) = 0.0;
    }
    // One of the correlation factors refers to an active surface
    else {

      for (int r = 0; r <= c; r++) {
        // Dummy value for correlation involving an inactive surface
        if (m_hits.at(r)->getObjectKind() == HitKind::Inactive) m_correlationsRZ(r, c) = 0.0;

        // Correlations between two active surfaces
        else {

          double sum = 0.0;

          for (int i = 0; i < r; i++) sum += msThetaOverSinSq.at(i)
                                            * (m_hits.at(c)->getDistance() - m_hits.at(i)->getDistance())
                                            * (m_hits.at(r)->getDistance() - m_hits.at(i)->getDistance());

          if (r == c) {
            double prec = m_hits.at(r)->getResolutionZ(m_radius);
            sum = sum + prec * prec;
          }

          m_correlationsRZ(r, c) = sum;
          if (r != c) m_correlationsRZ(c, r) = sum;
#undef CORRELATIONS_OFF_DEBUG
#ifdef CORRELATIONS_OFF_DEBUG
          if (r!=c) {
            m_correlationsRZ(c, r)=0;
            m_correlationsRZ(r, c)=0;
          }
#endif
        }
      }
    }
  } // Correlations: c is column, r is row

  // Remove zero rows and columns
  int  nResized        = -1;
  bool look_for_active = false;

  for (int i = 0; i < n; i++) {

    if ((m_hits.at(i)->getObjectKind() == HitKind::Inactive) && (!look_for_active)) {
      nResized = i;
      look_for_active = true;
    }
    else if ((m_hits.at(i)->getObjectKind() == HitKind::Active) && (look_for_active)) {

      for (int j = 0; j < n; j++) {
        m_correlationsRZ(nResized, j) = m_correlationsRZ(i, j);
        m_correlationsRZ(j, nResized) = m_correlationsRZ(j, i);
      }

      m_correlationsRZ(nResized, nResized) = m_correlationsRZ(i, i);
      nResized++;
    }
  }
  
  // Resize matrix if necessary
  if (nResized!=-1) m_correlationsRZ.ResizeTo(nResized, nResized);

  // Check if matrix is sane and worth keeping
  if (!((m_correlationsRZ.GetNoElements() > 0) && (m_correlationsRZ.Determinant() != 0.0))) {
    logWARNING("Correlation matrix in s-Z -> zero determinat or zero number of elements");
  }
}

//
// Compute the covariance matrix of the track parameters in R-Z projection
//
void Track::computeCovarianceMatrixRZ() {

  unsigned int offset = 0;
  unsigned int nHits  = m_hits.size();

  int n = m_correlationsRZ.GetNrows();

  TMatrixT<double> C(m_correlationsRZ); // Local copy to be inverted
  TMatrixT<double> diffsT(2, n);        // Derivatives of track parameters transposed (in R-Z -> 2 track parameters)
  TMatrixT<double> diffs(n, 2);         // Derivatives of track parameters (in R-Phi -> 3 track parameters)

  m_covariancesRZ.ResizeTo(2,2);
  
  // Set up partial derivative matrices diffs and diffsT -> line fit in s-Z to define these matrices
  for (auto i = 0; i < nHits; i++) {

    if (m_hits.at(i)->getObjectKind()  == HitKind::Active) {

      // Partial derivatives for x = p[0] * y + p[1]
      diffs(i - offset, 0) = m_hits.at(i)->getRadius();
      diffs(i - offset, 1) = 1;
    }
    else offset++;
  }

  // Transpose
  diffsT.Transpose(diffs);

  // Get covariance matrix using global chi2 fit: cov(i,j) = (D^T * C^-1 * D)^-1
  // TODO: check if this matrix can be inverted
  m_covariancesRZ = diffsT * C.Invert() * diffs;
  m_covariancesRZ.Invert();
}


//
//
///**
// * Changes some active hits into inactive
// * according to the efficiency
// * @param efficiency the modules active fraction
// * @param alsoPixel true if the efficiency removal applies to the pixel hits also
// */
//void Track::addEfficiency(double efficiency, bool pixel /* = false */ ) {
//  for (std::vector<Hit*>::iterator it = hitV_.begin(); it!=hitV_.end(); ++it) {
//    if ((*it)->getObjectKind() == Hit::Active) {
//      if ((pixel)&&(*it)->isPixel()) {
//	if ((double(random())/RAND_MAX) > efficiency) { // This hit is LOST
//	  (*it)->setObjectKind(Hit::Inactive);
//	}
//      }
//      if ((!pixel)&&(!(*it)->isPixel())) {
//	if ((double(random())/RAND_MAX) > efficiency) { // This hit is LOST
//	  (*it)->setObjectKind(Hit::Inactive);
//	}
//      }
//    }
//  }
//}
//
///**
// * Makes all non-trigger hits inactive
// */
//void Track::keepTriggerOnly() {
//  // int iRemove=0;
//  for (std::vector<Hit*>::iterator it = hitV_.begin(); it!=hitV_.end(); ++it) {
//    // if (debugRemoval) std::cerr << "Hit number "
//    //	                           << iRemove++ << ": ";
//    // if (debugRemoval) std::cerr << "r = " << (*it)->getRadius() << ", ";
//    // if (debugRemoval) std::cerr << "d = " << (*it)->getDistance() << ", ";
//    if ((*it)->getObjectKind() == Hit::Active) {
//      // if (debugRemoval) std::cerr << "active ";
//      if ((*it)->isPixel()) {
//	// if (debugRemoval) std::cerr << "pixel: removed";
//	(*it)->setObjectKind(Hit::Inactive);
//      } else {
//	DetectorModule* myModule = (*it)->getHitModule();
//	if (myModule) {
//	  // if (debugRemoval) std::cerr << "module ";
//	  if (myModule->sensorLayout() != PT) {
//	    // if (debugRemoval) std::cerr << "non-pt: removed";
//	    (*it)->setObjectKind(Hit::Inactive);
//	  } else {
//	    // if (debugRemoval) std::cerr << "pt: kept";
//	  }
//	} else {
//	  // if (debugRemoval) std::cerr << "active without module: kept";
//	}
//      }
//    } else {
//      // if (debugRemoval) std::cerr << "inactive";
//    }
//    // if (debugRemoval) std::cerr << std::endl;
//  }
//
//  // debugRemoval=false;
//}
//
//
//
//
///**
// * Sets all the hits to their trigger resolution
// */
//void Track::setTriggerResolution(bool isTrigger) {
//  Hit* myHit;
//  for (std::vector<Hit*>::iterator it = hitV_.begin(); it!=hitV_.end(); ++it) {
//    myHit = (*it);
//    if (myHit->getObjectKind() == Hit::Active) {
//       myHit->setTrigger(isTrigger);
//    }
//  }
//}
//
//
//
//
//
//
//double Track::expectedTriggerPoints(const double& triggerMomentum) const {
//  std::vector<Hit*>::const_iterator hitIt;
//  Hit* myHit;
//  double result=0;
//
//  for (hitIt=hitV_.begin();
//       hitIt!=hitV_.end();
//       ++hitIt) {
//    myHit=(*hitIt);
//    if ((myHit) &&
//	(myHit->isTrigger()) &&
//	(!myHit->isIP()) &&
//	(myHit->getObjectKind()==Hit::Active)) {
//      // We've got a possible trigger here
//      // Let's find the corresponding module
//      DetectorModule* myModule = myHit->getHitModule();
//      if (myModule) {
//	result += PtErrorAdapter(*myModule).getTriggerProbability(triggerMomentum);
//      } else {
//	// Whoops: problem here: an active hit is not linked to any module
//	std::cerr << "ERROR: this SHOULD NOT happen. in expectedTriggerPoints() an active hit does not correspond to any module!" << std::endl;
//      }
//    }
//  }
//  return result;
//}

/**
 * @file TrackNew.cc
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "TrackNew.hh"

#include <algorithm>
#include <cstdlib>

#include <global_constants.hh>
#include "HitNew.hh"
#include "MessageLogger.hh"
#include "MaterialProperties.hh"
#include "SimParms.hh"
#include "Units.hh"


using namespace ROOT::Math;
using namespace std;

//
// TrackNew constructor -> need to use setter methods to set: 2 of these [theta, phi, eta, cot(theta)] & 2 of these [mag. field, transv. momentum, radius]
//
TrackNew::TrackNew() :
  m_theta(0),
  m_phi(0),
  m_cotgTheta(0),
  m_eta(0),
  m_pt(0),
  m_reSortHits(true),
  m_covRPhiDone(false),
  m_covRZDone(false),
  m_refPointRPosCache(0),
  m_propagOutInCache(true)
{}

//
// TrackNew copy-constructor -> creates deep copy of hit vector
//
TrackNew::TrackNew(const TrackNew& track) {

  m_theta        = track.m_theta;
  m_phi          = track.m_phi;
  m_cotgTheta    = track.m_cotgTheta;
  m_eta          = track.m_eta;
  m_pt           = track.m_pt;

  m_origin       = track.m_origin;
  m_direction    = track.m_direction;

  m_reSortHits   = track.m_reSortHits;
  m_covRPhiDone  = track.m_covRPhiDone;
  m_covRZDone    = track.m_covRZDone;
  m_refPointRPosCache = track.m_refPointRPosCache;
  m_propagOutInCache  = track.m_propagOutInCache;

  m_varMatrixRPhi.ResizeTo(track.m_varMatrixRPhi);
  m_varMatrixRPhi = track.m_varMatrixRPhi;
  m_covMatrixRPhi.ResizeTo(track.m_covMatrixRPhi);
  m_covMatrixRPhi = track.m_covMatrixRPhi;

  m_varMatrixRZ.ResizeTo(track.m_varMatrixRZ);
  m_varMatrixRZ = track.m_varMatrixRZ;
  m_covMatrixRZ.ResizeTo(track.m_covMatrixRZ);
  m_covMatrixRZ = track.m_covMatrixRZ;

  for (auto& iHit : track.m_hits) {
    HitNewPtr hit(new HitNew(*iHit));
    addHit(std::move(hit));
  }
  m_tags = track.m_tags;
}

//
// Assign operator with deep copy of hit vector
//
TrackNew& TrackNew::operator= (const TrackNew& track) {

  // check for self-assignment by comparing the address of the
  // implicit object and the parameter
  if (this == &track) return *this;
  
  // Do the copy
  m_theta        = track.m_theta;
  m_phi          = track.m_phi;
  m_cotgTheta    = track.m_cotgTheta;
  m_eta          = track.m_eta;
  m_pt           = track.m_pt;

  m_origin       = track.m_origin;
  m_direction    = track.m_direction;

  m_reSortHits   = track.m_reSortHits;
  m_covRPhiDone  = track.m_covRPhiDone;
  m_covRZDone    = track.m_covRZDone;
  m_refPointRPosCache = track.m_refPointRPosCache;
  m_propagOutInCache  = track.m_propagOutInCache;

  m_varMatrixRPhi.ResizeTo(track.m_varMatrixRPhi);
  m_varMatrixRPhi = track.m_varMatrixRPhi;
  m_covMatrixRPhi.ResizeTo(track.m_covMatrixRPhi);
  m_covMatrixRPhi = track.m_covMatrixRPhi;

  m_varMatrixRZ.ResizeTo(track.m_varMatrixRZ);
  m_varMatrixRZ = track.m_varMatrixRZ;
  m_covMatrixRZ.ResizeTo(track.m_covMatrixRZ);
  m_covMatrixRZ = track.m_covMatrixRZ;

  for (auto& iHit : track.m_hits) {
    HitNewPtr hit(new HitNew(*iHit));
    addHit(std::move(hit));
  }
  m_tags = track.m_tags;

  // Return the existing object
  return *this;
}

//
// Destructor
//
TrackNew::~TrackNew() {

  // Clear memory
  m_hits.clear();
}

//
// Main method calculating track parameters in s-z plane only, using linear fit with parameters: cotg(theta), z0 -> internally calling computation of covMatrixRZ
// As the Multiple scattering effects must be set in a way to have then correct propagation of errors up-to ref. point [rPos,zPos] (including all dead materials between the
// last measurement and the ref. point). E.g. for standard estimation of D0,Z0 parameters one calculates MS inside->out from rPos=0 (zPos can be calculated from rPos using theta).
// PropagOutIn variable defines whether detectors at higher R than the ref. point (true) should be used for error calculation/propagation (e.g. d0,z0) or whether detectors at
// lower R (false).
// Return true if errors correctly calculated
//
bool TrackNew::computeErrorsRZ(double refPointRPos/*=0*/, bool propagOutIn/*=true*/) {

  // Sort hits based on particle direction: in-out or out-in (if needed)
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  // Compute the relevant 2x2 covariance matrix in RZ plane
  m_covRZDone         = computeCovarianceMatrixRZ(fabs(refPointRPos),propagOutIn);
  m_refPointRPosCache = refPointRPos;
  m_propagOutInCache  = propagOutIn;

  return m_covRZDone;
}

//
// Main method calculating track parameters in r-phi plane only, using parabolic track approximation in R-Phi plane: 1/R, d0, phi0 parameters -> internally calling computation
// of covMatrixRPhi. As the Multiple scattering effects must be set in a way to have then correct propagation of errors up-to the ref. point [rPos,zPos] (including all dead
// materials between the last measurement and the ref. point). E.g. for standard estimation of D0,Z0 parameters one calculates MS inside->out from rPos=0 (zPos can be
// calculated from rPos using theta). PropagOutIn variable defines whether detectors at higher R than the ref. point (true) should be used for error calculation/propagation
// (e.g. d0,z0) or whether detectors at lower R (false).
// Return true if errors correctly calculated
//
bool TrackNew::computeErrorsRPhi(double refPointRPos/*=0*/, bool propagOutIn/*=true*/) {

  // Sort hits based on particle direction: in-out or out-in (if needed)
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  // Compute the relevant 3x3 covariance matrix in RPhi plane
  m_covRPhiDone       = computeCovarianceMatrixRPhi(fabs(refPointRPos),propagOutIn);
  m_refPointRPosCache = refPointRPos;
  m_propagOutInCache  = propagOutIn;

  return m_covRPhiDone;
}

//
// Calculate magnetic field at given z, assuming B = B(z).e_z + 0.e_x + 0 e_y
//
double TrackNew::getMagField(double z) const {

  double magField = 0;

//  // Option 1: Const mag. field across the detector: Bz = const (assumed true for now)
  if (SimParms::getInstance().isMagFieldConst()) {

    magField = SimParms::getInstance().magField();
  }
// TODO:
//  // Option 2: Mag. field is a function in Z: B = B(z)
//  else {
//
//    for (unsigned int i=0; i<SimParms::getInstance().getNMagFieldRegions(); i++) {
//
//      // Magnetic field regions are considered to be defined in metres
//      if (z<SimParms::getInstance().magFieldZRegions[i]) {
//
//        // Magnetic field is considered to be in Tesla
//        magField = SimParms::getInstance().magField[i];
//        break;
//      }
//    }
//  }

  return magField;
}

//
// Get DeltaRho (error on 1/R) at refPoint [rPos, zPos].
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
// Using 3x3 covariance propagator in case [r,z]!=[0,0]
//
double TrackNew::getDeltaRho(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in R-Phi if something changed
  if (!m_covRPhiDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeErrorsRPhi(refPointRPos, propagOutIn);

  double deltaRho = -1.;

  if (m_covRPhiDone && m_covMatrixRPhi(0,0)>=0) deltaRho = sqrt(m_covMatrixRPhi(0,0));

  // TODO: Not working correctly for B = B(z) & [r,z]!=[0,0], so print warning ...
  if (refPointRPos!=0. && !SimParms::getInstance().isMagFieldConst()) {

    logWARNING("TrackNew::getDeltaRho(): Mathematical method to get deltaRho at [r,z]!=[0,0] in non const. B field not implemented, hence returned value at [r,z]=[0,0].");
  }

  return deltaRho;
}

//
// Get DeltaPtOvePt at refPoint [rPos, zPos] (utilize the calculated deltaRho quantity)
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
//
double TrackNew::getDeltaPtOverPt(double refPointRPos, bool propagOutIn/*=true*/) {

  double deltaPtOverPt = -1.;

  // delta(1/R) & delta(pT) -> estimated at point [r,z] = [0,0] (important for use case, when B != const -> B = B(z))
  double deltaRho = getDeltaRho(refPointRPos,propagOutIn);
  double radius   = getRadius(refPointRPos*m_cotgTheta);       // Approximative transformation from rPos to zPos using tan(theta)
  if (deltaRho!=-1) deltaPtOverPt = deltaRho * radius; // dpT(z)/pT(z) = dRho(z) / Rho(z) = dRho(z) * R(z)

  return deltaPtOverPt;
}

//
// Get DeltaPOverP at refPoint [rPos, zPos] (utilize deltaRho & deltaCotgTheta quantities)
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
//
double TrackNew::getDeltaPOverP(double refPointRPos, bool propagOutIn/*=true*/) {

  double deltaPOverP = -1.;

  // Combining into p measurement
  // dp/p = dp_t/p_t + A / (1+A^2) * dA // with A = ctg(theta)
  // dp/p = dp_t/p_t + sin(theta)*cos(theta)*dcotg(theta)
  double deltaPtOverPt = getDeltaPtOverPt(refPointRPos,propagOutIn);
  double deltaCtgTheta = getDeltaCtgTheta(refPointRPos,propagOutIn);
  if (deltaPtOverPt!=-1 && deltaCtgTheta!=-1) deltaPOverP = sqrt(deltaPtOverPt*deltaPtOverPt + sin(m_theta)*sin(m_theta) * cos(m_theta)*cos(m_theta) * deltaCtgTheta*deltaCtgTheta);

  return deltaPOverP;
}

//
// Get DeltaPhi(0) at refPoint [rPos, zPos] ([0,0])
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
// Using 3x3 covariance propagator in case [r,z]!=[0,0]
//
double TrackNew::getDeltaPhi(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in R-Phi if something changed
  if (!m_covRPhiDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeErrorsRPhi(refPointRPos, propagOutIn);

  double deltaPhi0 = -1.;

  // No covariance propagation necessary at [0,0,0] point
  if (refPointRPos==0.) {

    if (m_covRPhiDone && m_covMatrixRPhi(1,1)>=0) deltaPhi0 = sqrt(m_covMatrixRPhi(1, 1));
  }
  else {

    if (m_covRPhiDone) {

      double covRhoRho    = m_covMatrixRPhi(0,0);
      double covRhoPhi0   = m_covMatrixRPhi(0,1);
      double covRhoD0     = m_covMatrixRPhi(0,2);
      double covPhi0Phi0  = m_covMatrixRPhi(1,1);
      double covPhi0D0    = m_covMatrixRPhi(1,2);
      double covD0D0      = m_covMatrixRPhi(2,2);

      double rho          = getRho(refPointRPos*m_cotgTheta);
      double rho2         = rho*rho;
      double rho4         = rho2*rho2;
      double refPointRPos2= refPointRPos*refPointRPos;

      double deltaPhi0Sq  = refPointRPos2*covRhoRho   + covPhi0Phi0                   + rho4*refPointRPos2*covD0D0;
             deltaPhi0Sq += 2*refPointRPos*covRhoPhi0 - 2*rho2*refPointRPos2*covRhoD0 - 2*rho2*refPointRPos*covPhi0D0;

      if (deltaPhi0Sq>=0) deltaPhi0 = sqrt(deltaPhi0Sq);
    }
  }

  // TODO: Not working correctly for B = B(z) & [r,z]!=[0,0], so print warning ...
  if (refPointRPos!=0. && !SimParms::getInstance().isMagFieldConst()) {

    logWARNING("TrackNew::getDeltaPhi(): Mathematical method to get deltaPhi at [r,z]!=[0,0] in non const. B field not implemented, hence returned value at [r,z]=[0,0].");
  }

  return deltaPhi0;
}

//
// Get DeltaD(0) at refPoint [rPos, zPos] ([0,0])
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
// Using 3x3 covariance propagator in case [r,z]!=[0,0]
//
double TrackNew::getDeltaD(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in R-Phi if something changed
  if (!m_covRPhiDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeErrorsRPhi(refPointRPos, propagOutIn);

  double deltaD0 = -1.;

  // No covariance propagation necessary at [0,0,0] point
  if (refPointRPos==0.) {

    if (m_covRPhiDone && m_covMatrixRPhi(2,2)>=0) deltaD0 = sqrt(m_covMatrixRPhi(2,2));
  }
  else {

    if (m_covRPhiDone) {

      double covRhoRho     = m_covMatrixRPhi(0,0);
      double covRhoPhi0    = m_covMatrixRPhi(0,1);
      double covRhoD0      = m_covMatrixRPhi(0,2);
      double covPhi0Phi0   = m_covMatrixRPhi(1,1);
      double covPhi0D0     = m_covMatrixRPhi(1,2);
      double covD0D0       = m_covMatrixRPhi(2,2);

      double refPointRPos2 = refPointRPos*refPointRPos;
      double refPointRPos3 = refPointRPos2*refPointRPos;
      double refPointRPos4 = refPointRPos3*refPointRPos;

      double deltaD0Sq  = refPointRPos4/4.*covRhoRho + refPointRPos3*covRhoPhi0 + refPointRPos2*covRhoD0;
             deltaD0Sq += refPointRPos2*covPhi0Phi0  + 2*refPointRPos*covPhi0D0 + covD0D0;

      if (deltaD0Sq>=0) deltaD0 = sqrt(deltaD0Sq);
    }
  }

  // TODO: Not working correctly for B = B(z) & [r,z]!=[0,0], so print warning ...
  if (refPointRPos!=0. && !SimParms::getInstance().isMagFieldConst()) {

    logWARNING("TrackNew::getDeltaD(): Mathematical method to get deltaD0 at [r,z]!=[0,0] in non const. B field not implemented, hence returned value at [r,z]=[0,0].");
  }

  return deltaD0;
}

//
// Get DeltaCtgTheta at refPoint [rPos, zPos]
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
//
double TrackNew::getDeltaCtgTheta(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in s-Z if something changed
  if (!m_covRZDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeErrorsRZ(refPointRPos, propagOutIn);

  double deltaCtgTheta = -1.;
  if (m_covRZDone && m_covMatrixRZ(0, 0)>=0) deltaCtgTheta = sqrt(m_covMatrixRZ(0, 0));

  return deltaCtgTheta;
}

//
// Get DeltaZ(0) at refPoint [rPos, zPos] ([0,0])
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
// Using 2x2 covariance propagator in case [r,z]!=[0,0]
//
double TrackNew::getDeltaZ(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in s-Z if something changed
  if (!m_covRZDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeErrorsRZ(refPointRPos, propagOutIn);

  double deltaZ0 = -1.;

  // No covariance propagation necessary at [0,0,0] point
  if (refPointRPos==0.) {

    if (m_covRZDone && m_covMatrixRZ(1, 1)>=0) deltaZ0 = sqrt(m_covMatrixRZ(1, 1));
  }
  else {

    if (m_covRZDone) {

      double covZ0Z0       = m_covMatrixRZ(1,1);
      double covZ0CtgTheta = m_covMatrixRZ(0,1);
      double covCtgThCtgTh = m_covMatrixRZ(0,0);

      //std::cout << "<< " << rPos << " " << covZ0Z0 << " " << covZ0CtgTheta << "  " << covCtgThCtgTh << std::endl;

      double deltaZ0Sq = covZ0Z0 + 2*refPointRPos*covZ0CtgTheta + refPointRPos*refPointRPos*covCtgThCtgTh;
      if (deltaZ0Sq>=0) deltaZ0 = sqrt(deltaZ0Sq);
    }
  }

  return deltaZ0;
}

//
// Get DeltaCTau for secondary particles coming from the primary vertex at ~ [0,0] -> an important quantity to estimate the
// resolution of secondary vertices
//
double TrackNew::getDeltaCTau() {

  double deltaCTau = -1;
  double deltaD0   = getDeltaD0();
  double deltaZ0   = getDeltaZ0();

  if (deltaD0>0 && deltaZ0>0) {
    deltaCTau = (cos(m_theta)*cos(m_theta)+1)/deltaD0/deltaD0 + sin(m_theta)*sin(m_theta)/deltaZ0/deltaZ0;
    deltaCTau = sqrt(2/deltaCTau);
  }

  return deltaCTau;
}

//
// Adds a new hit to the track (hit radius automatically updated)
//
void TrackNew::addHit(HitNewPtr newHit) {

  // Add tracking tags
  if (newHit->getHitModule() != nullptr) {
    m_tags.insert(newHit->getHitModule()->trackingTags.begin(), newHit->getHitModule()->trackingTags.end());
  }
  newHit->setTrack(this);

  // Add new hit if it follows the parabolic approximation -> hits practically found at high pT limit, so in reality don't have to lie on the track
  if (followsParabolicApprox(newHit->getRPos(), newHit->getZPos())) m_hits.push_back(std::move(newHit));
  else newHit.reset(nullptr);

  // Hits need to be re-sorted & cov. matrices recalculated
  m_reSortHits  = true;
  m_covRPhiDone = false;
  m_covRZDone   = false;
}

//
// Add IP constraint to the track, technically new hit is assigned: with no material and hit resolution in R-Phi as dr, in s-Z as dz
//
void TrackNew::addIPConstraint(double dr, double dz) {

  // This modeling of the IP constraint was validated:
  // By placing dr = 0.5 mm and dz = 1 mm one obtains
  // sigma(d0) = 0.5 mm and sigma(z0) = 1 mm
  HitNewPtr newHit(new HitNew(0,0)); //(dr,dz)); // TODO: Cross-check, should be Hit(0,0) ???
  newHit->setIP(true);

  RILength emptyMaterial;
  emptyMaterial.radiation   = 0;
  emptyMaterial.interaction = 0;

  newHit->setCorrectedMaterial(emptyMaterial);
  newHit->setAsActive();
  newHit->setResolutionRphi(dr);
  newHit->setResolutionZ(dz);

  // Add new hit if it follows the parabolic approximation -> hits practically found at high pT limit, so in reality don't have to lie on the track
  if (followsParabolicApprox(newHit->getRPos(), newHit->getZPos())) m_hits.push_back(std::move(newHit));
  else newHit.reset(nullptr);

  // Hits need to be re-sorted & cov. matrices recalculated
  m_reSortHits  = true;
  m_covRPhiDone = false;
  m_covRZDone   = false;
}

//
// Simulate efficiency by changing some active hits to non-active hits (passive)
//
void TrackNew::addEfficiency() {
  for (auto& iHit : m_hits) {
    if (iHit->isActive()) {
      double efficiency = iHit->getHitModule()->singleHitEfficiency();
      if (efficiency!=1) {
        if ((double(random())/RAND_MAX)>efficiency) iHit->setAsPassive(); // This hit is LOST
      }
    }
  }
}

//
// Set track polar angle - theta, azimuthal angle - phi, particle transverse momentum - pt
// (magnetic field obtained automatically from SimParms singleton class)Setter for the track azimuthal angle.
//
const Polar3DVector& TrackNew::setThetaPhiPt(const double& newTheta, const double& newPhi, const double& newPt) {

  m_theta     = newTheta;
  m_cotgTheta = 1/tan(newTheta);
  m_eta       = -log(tan(m_theta/2));
  m_phi       = newPhi;
  m_pt        = newPt;

  if (m_pt>=0) m_direction.SetCoordinates(+1, m_theta, m_phi); // Particle inside-out
  else         m_direction.SetCoordinates(-1, m_theta, m_phi); // Particle outside-in

  // Clear all previously assigned hits -> hits need to be recalculated
  m_hits.clear();

  // Hits need to be re-sorted & cov. matrices recalculated
  m_reSortHits  = true;
  m_covRPhiDone = false;
  m_covRZDone   = false;

  return m_direction;
}

//
// Re-set transverse momentum + resort hits (if changing direction) + initiate recalc of cov matrices + prune hits (otherwise they may not lie on the new track, originally found at high pT limit)
void TrackNew::resetPt(double newPt) {

  if (newPt*m_pt<0) m_reSortHits = true;
  m_covRPhiDone = false;
  m_covRZDone   = false;

  m_pt = newPt;

  pruneHits();
}

//
// Sort internally all hits assigned to this track -> sorting algorithm based on hit radius - by smaller radius sooner or vice-versa (inner-2-outer approach or vice-versa)
//
void TrackNew::sortHits(bool bySmallerR) { bySmallerR ? std::stable_sort(m_hits.begin(), m_hits.end(), HitNew::sortSmallerR) : std::stable_sort(m_hits.begin(), m_hits.end(), HitNew::sortHigherR); }

//
// Remove hits that don't follow the parabolic approximation used in tracking - TODO: still needs to be updated (not all approximations taken into account here)
//
bool TrackNew::pruneHits() {

  bool isPruned = false;

  HitNewCollection newHits;
  for (auto& iHit : m_hits) {

    if (followsParabolicApprox(iHit->getRPos(),iHit->getZPos())) newHits.push_back(std::move(iHit));
    else {

      // Clear memory
      iHit.reset();

      isPruned = true;
    }
  }

  m_hits.clear();
  for (auto& iHit : newHits) m_hits.push_back(std::move(iHit));

  return isPruned;
}

//
// Set active only hits with the given tag
//
void TrackNew::keepTaggedHitsOnly(const string& tag, bool useIP /*=true*/) {

  for (auto& iHit : m_hits) {

    // IP constraint hit
    if (tag=="all" && iHit->isIP() && useIP) iHit->setAsActive();

    // Measurement hit
    if (iHit->isMeasurable()) {
      if (tag=="all") iHit->setAsActive();
      else {
        if (std::count_if(iHit->getHitModule()->trackingTags.begin(), iHit->getHitModule()->trackingTags.end(), [&tag](const string& s){ return s == tag; })) iHit->setAsActive();
        else iHit->setAsPassive();
      }
    }
  }

  // Cov. matrices need to be recalculated
  m_covRPhiDone = false;
  m_covRZDone   = false;
}

//
// Remove material from all assigned hits -> modify all hits such as they are without any material
//
void TrackNew::removeMaterial() {

  // Material object with no material assigned
  RILength nullMaterial;

  // Reset all material assigned to hits
  for (auto& iHit : m_hits) iHit->setCorrectedMaterial(nullMaterial);

  // Cov. matrices need to be recalculated
  m_covRPhiDone = false;
  m_covRZDone   = false;
}

//
// Helper method printing track covariance matrices in R-Phi
//
void TrackNew::printErrors() {

  std::cout << "Overview of track errors:" << std::endl;
  std::cout << "Hit variance matrix: "  << std::endl;
  m_varMatrixRPhi.Print();

  std::cout << "Covariance matrix: " << std::endl;
  m_covMatrixRPhi.Print();

  // Print errors @ [r,z]=[0,0]
  double rPos = 0.0;

  std::cout << "Rho errors by momentum: " << getDeltaRho(rPos) << std::endl;
  std::cout << "Phi0 errors by momentum: "<< getDeltaPhi0()    << std::endl;
  std::cout << "D0 errors by momentum: "  << getDeltaD0()      << std::endl;
}

//
// Helper method printing symmetric matrix
//
void TrackNew::printSymMatrix(const TMatrixTSym<double>& matrix) const {

  std::cout << std::endl;

  int nCols = matrix.GetNcols();
  int nRows = matrix.GetNrows();

  for (int i = 0; i<nRows; i++) {
    std::cout << "(";
    for (int j=0; j<nCols;j++) {

      std::cout << " " << std::scientific << std::setprecision(4) << matrix(i,j);
    }
    std::cout << ")" << std::endl;
  }
  std::cout << std::endl;
}

//
// Helper method printing matrix
//
void TrackNew::printMatrix(const TMatrixT<double>& matrix) const {

  std::cout << std::endl;

  int nCols = matrix.GetNcols();
  int nRows = matrix.GetNrows();

  for (int i = 0; i<nRows; i++) {
    std::cout << "(";
    for (int j=0; j<nCols;j++) {

      std::cout << " " << std::fixed << std::setprecision(5) << matrix(i,j);
    }
    std::cout << ")" << std::endl;
  }
  std::cout << std::endl;
}

//
// Helper method printing track hits
//
void TrackNew::printHits() const {

  std::cout << "******************" << std::endl;
  std::cout << "Track eta=" << m_eta << std::endl;

  for (const auto& it : m_hits) {
    std::cout << "    Hit";
    if (it->isActive())   std::cout << " r="  << it->getRPos() << " +- " << it->getResolutionRphi(getRadius(it->getZPos()));
    else                  std::cout << " r="  << it->getRPos();
    if (it->isActive())   std::cout << " z="  << it->getZPos() << " +- " << it->getResolutionZ(getRadius(it->getZPos()));
    else                  std::cout << " z="  << it->getZPos();
    std::cout << " d="  << it->getDistance()
              << " rl=" << it->getCorrectedMaterial().radiation
              << " il=" << it->getCorrectedMaterial().interaction;
    if (it->isActive())   std::cout << " active";
    else                  std::cout << " inactive";
    if (it->isBarrel())   std::cout << " barrel";
    if (it->isEndcap())   std::cout << " endcap";
    if (it->isBeamPipe()) std::cout << " beam-pipe";
    if (it->isIP())       std::cout << " ip";
    if (it->getLayerOrDiscID()!=-1) std::cout << " " << it->getDetName() << " L/D_id= " << it->getLayerOrDiscID();

    if (it->isActive()) {
      std::cout << " activeHitType_=" << static_cast<short>(it->getActiveHitType());
    }
    std::cout << std::endl;
  }
}

//
// Helper method printing track hits
//
void TrackNew::printActiveHits() const {

  std::cout << "******************" << std::endl;
  std::cout << "Track eta=" << m_eta << std::endl;

  for (const auto& it : m_hits) {
    if (it->isActive()) {

      std::cout << "    Hit";
      std::cout << " r="  << it->getRPos() << " +- " << it->getResolutionRphi(getRadius(it->getZPos()));
      std::cout << " z="  << it->getZPos() << " +- " << it->getResolutionZ(getRadius(it->getZPos()));
      std::cout << " d="  << it->getDistance()
                << " rl=" << it->getCorrectedMaterial().radiation
                << " il=" << it->getCorrectedMaterial().interaction;
      if (it->isActive())   std::cout << " active";
      else                  std::cout << " inactive";
      if (it->isBarrel())   std::cout << " barrel";
      if (it->isEndcap())   std::cout << " endcap";
      if (it->isBeamPipe()) std::cout << " beam-pipe";
      if (it->isIP())       std::cout << " ip";
      if (it->getLayerOrDiscID()!=-1) std::cout << " " << it->getDetName() << " L/D_id= " << it->getLayerOrDiscID();
      std::cout << " activeHitType_=" << static_cast<short>(it->getActiveHitType());
      std::cout << std::endl;
    }
  }
}

//
// Get number of active hits assigned to track for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file). If tag specified as "all" no extra tag required
//
int TrackNew::getNActiveHits (std::string tag, bool useIP /* = true */ ) const {

  // Result variable
  int nHits=0;

  for (auto& iHit : m_hits) {
    if (iHit && iHit->isActive()){
      if (iHit->isIP() && useIP) {
        nHits++;
      }
      else if (!iHit->isIP()) {

        // Check tag for non-IP assigned hits
        bool tagOK = false;
        for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
          if (tag==*it || tag=="all") tagOK = true;
        }

        if (tagOK) nHits++;
      }
    }
  } // For

  return nHits;
}

//
// Get number of active hits coming from measurement planes or IP constraint assigned to track for given tag. If tag specified as "all", all module & IP hits assigned.
//
int TrackNew::getNMeasuredHits(std::string tag, bool useIP /*=true*/) const {

  // Result variable
  int nHits=0;

  for (auto& iHit : m_hits) {
    if (iHit && iHit->isActive()) {
      if (iHit->isIP() && useIP) {
        nHits++;
      }
      else if (iHit->isMeasurable()) {

        // Check tag for non-IP assigned hits
        bool tagOK = false;
        for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
          if (tag==*it || tag=="all") tagOK = true;
        }
        if (tagOK) nHits++;
      }
    }
  } // For

  return nHits;

}

//
// Get reference to a hit, which can be measured, i.e. coming from measurement plane (active or inactive) or IP constraint
//
const HitNew* TrackNew::getMeasurableOrIPHit(int iHit) {

  int   hitCounter = 0;
  const HitNew* pHit  = nullptr;

  // Sort hits based on particle direction: in-out or out-in
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  for (auto& hit : m_hits) {
    if (hit && (hit->isIP() || hit->isMeasurable())) {

      // Hit we're looking for!
      if (hitCounter==iHit) {

        pHit = hit.get();
        break;
      }
      hitCounter++;
    }
  }

  return pHit;
}

//
// Reverse search - Get reference to a hit, which can be measured, i.e. coming from measurement plane (active or inactive) or IP constraint
//
const HitNew* TrackNew::getRMeasurableOrIPHit(int iHit) {

  int   hitCounter = 0;
  const HitNew* pHit  = nullptr;

  // Sort hits based on particle direction: in-out or out-in
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  for (auto hit = m_hits.rbegin(); hit != m_hits.rend(); ++hit) {
    if (*hit && ((*hit)->isIP() || (*hit)->isMeasurable())) {

      // Hit we're looking for!
      if (hitCounter==iHit) {

        pHit = (*hit).get();
        break;
      }
      hitCounter++;
    }
  }

  return pHit;
}

//
// Get the probabilty of having "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
// If tag specified as "all" no extra tag required
//
std::vector<double> TrackNew::getHadronActiveHitsProbability(std::string tag) {

  // Result variable
  std::vector<double> probabilities;
  double probability = 1;

  // Sort hits based on particle direction: in-out or out-in
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  for (auto& iHit : m_hits) {
    if (iHit) {
      if (iHit->isActive()){

        // Check tag
        bool tagOK = false;
        for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
          if (tag==*it || tag=="all") tagOK = true;
        }

        if (tagOK) probabilities.push_back(probability);
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
// If tag specified as "all" no extra tag required
//
double TrackNew::getHadronActiveHitsProbability(std::string tag, int nHits) {

  // Probability
  double probability = 1;

  // Number of clean hits
  int goodHits = 0;

  // Sort hits based on particle direction: in-out or out-in
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  for (auto& iHit : m_hits) {

    if (iHit) {
      if (iHit->isActive()) {

        // Check tag
        bool tagOK = false;
        for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
          if (tag==*it || tag=="all") tagOK = true;
        }

        if (tagOK) goodHits++;
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
RILength TrackNew::getMaterial() const {

  RILength totalMaterial;
  totalMaterial.radiation   = 0;
  totalMaterial.interaction = 0;

  for (auto& iHit : m_hits) totalMaterial += iHit->getCorrectedMaterial();

  return totalMaterial;
}

//
// Get a vector of pairs: Detector module & hit type for Trigger hits
//
std::vector<std::pair<const DetectorModule*, HitType>> TrackNew::getHitModules() const {

  std::vector<std::pair<const DetectorModule*, HitType>> result;

  for (auto& iHit : m_hits) {

    if ((iHit) && (iHit->isTrigger()) && (!iHit->isIP()) && (iHit->isActive())) {

      // We've got a possible trigger here
      // Let's find the corresponding module
      const DetectorModule* myModule = iHit->getHitModule();
      if (myModule) result.push_back(std::make_pair(myModule, iHit->getActiveHitType()));
      else {
        // Whoops: problem here: an active hit is not linked to any module
        logERROR("Track::getHitModules: This SHOULD NOT happen. In expectedTriggerPoints() an active hit does not correspond to any module!");
      }
    }
  }
  return result;
}

//
// Compute 3x3 covariance matrix of the track parameters in R-Phi projection, using NxN (N hits = K+L: K active hits on detectors + L passive (artificial) hits due to material)
// Ref. point dictates, whether Multiple scattering effects need to be calculated inside-out or outside-in. MS effect is symmetric as regards track fitting. PropagOutIn variable
// defines whether detectors at higher R than the ref. point (true) should be used for error calculation/propagation (e.g. d0,z0) or whether detectors at lower R (false).
// Return true if V invertable
//
bool TrackNew::computeCovarianceMatrixRPhi(double refPointRPos, bool propagOutIn) {

  // Matrix size
  int nHits = m_hits.size();
  m_varMatrixRPhi.ResizeTo(nHits,nHits);
  m_varMatrixRPhi.Zero();

  //
  // Find hit index ranges relevant for MS effects calculation based on requirements on refPoint & propagation direction
  int iStart          = nHits;
  int iEnd            = nHits-1;
  int nHitsUsed       = 0;
  int nActiveHitsUsed = 0;

  bool bySmallerR = true;
  bool useIP      = false;

  // Particle traverses inside-out with error propagation outside-in -> hits already sorted correctlly
  if (m_pt>=0 && propagOutIn) {

    for (auto it=m_hits.rbegin(); it!=m_hits.rend(); it++) {

      if (refPointRPos<(*it)->getRPos()) {
        if ((*it)->isMeasurable()) nActiveHitsUsed++;
        iStart--;
      }
      else break;
    }
  }
  // Particle traverses inside-out with error propagation inside-out -> sort hits by higher radius & resort backwards after matrix calculated
  else if (m_pt>=0 && !propagOutIn) {

    sortHits(!bySmallerR);
    for (auto it=m_hits.rbegin(); it!=m_hits.rend(); it++) {

      if (refPointRPos>(*it)->getRPos()) {
        if ((*it)->isMeasurable()) nActiveHitsUsed++;
        iStart--;
      }
      else break;
    }
  }
  // TODO: Implement
  else {

    logWARNING("Variance matrix V(NxN) in R-Phi -> calculations of particle traversing outside-in not yet implemented.");
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); // Important -> resort back
    return false;
  }

  // Check that enough hits for RPhi calculation, i.e. >=3
  nHitsUsed = nHits - iStart;
  if ((nActiveHitsUsed)<3) {

    std::string message = "Variance matrix V(NxN) in R-Phi -> refPointRPos[mm]="+any2str(refPointRPos/Units::mm,1);
                message+= " in combination with propagator direction doesn't provide sufficient number of hits for tracking!";
    logWARNING(message);
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); // Important -> resort back
    return false;
  }

  // Get contributions from Multiple Couloumb scattering
  std::vector<double> msThetaOverSinSq;

  for (int i=iStart; i<iEnd; i++) {

    // MS theta
    double msTheta = 0.0;

    // Material in terms of rad. lengths
    double XtoX0 = m_hits.at(i)->getCorrectedMaterial().radiation;
    //std::cout << std::fixed << std::setprecision(4) << "Material (" << i << ") = " << XtoX0 << "\t at r=" << m_hits.at(i)->getRadius(m_hits.at(i)->getZPos()) << "\t of type=" << m_hits.at(i)->getObjectKind() << std::endl;

    if (XtoX0>0) {

      // MS error depends on path length = deltaR/sin(theta), so one can precalculate msTheta_real as msTheta/sin^2(theta), which practically means using pT
      // instead of p & then one has to multiply the msTheta by deltaR to get MS error
      msTheta = (13.6*Units::MeV * 13.6*Units::MeV) / (m_pt/Units::MeV * m_pt/Units::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));

      // Take into account a propagation of MS error on virtual barrel plane, on which all measurements are evaluated for consistency (global chi2 fit applied) ->
      // in limit R->inf. propagation along line used, otherwise a very small correction factor coming from the circular shape of particle track is required (similar
      // approach as for local resolutions)
      // TODO: Currently, correction mathematicaly derived only for use case of const magnetic field -> more complex mathematical expression expected in non-const B field
      // (hence correction not applied in such case)
      double A = 0;
      if (SimParms::getInstance().isMagFieldConst()) A = m_hits.at(i)->getRPos()/2./getRadius(m_hits.at(i)->getZPos());     // r_i/2R
      double corrFactor = 1 + A*A*cos(m_theta)*cos(m_theta)/(1-A*A);

      msTheta *= corrFactor;
    }
    else {
      msTheta = 0;
    }
    msThetaOverSinSq.push_back(msTheta);
  }

  //
  // Calculate the variance matrix with all correlation terms: c is column, r is row (hits are assumed to be sorted)
  for (int c=iStart ; c<=iEnd; c++) {

    // Dummy value for correlations involving inactive surfaces
    if (m_hits.at(c)->isPassive()) {
      for (int r = 0; r <= c; r++) m_varMatrixRPhi(r, c) = 0.0;
    }
    // One of the correlation factors refers to an active surface
    else {

      for (int r=iStart; r <= c; r++) {
        // Dummy value for correlation involving an inactive surface
        if (m_hits.at(r)->isPassive()) m_varMatrixRPhi(r, c) = 0.0;

        // Correlations between two active surfaces
        else {

          double sum = 0.0;

          for (int i=iStart; i<r; i++) {

            sum += msThetaOverSinSq.at(i-(nHits-nHitsUsed))
                   * (m_hits.at(c)->getRPos() - m_hits.at(i)->getRPos())
                   * (m_hits.at(r)->getRPos() - m_hits.at(i)->getRPos());

          }
          if (r == c) {

            double prec = m_hits.at(r)->getResolutionRphi(getRadius(m_hits.at(r)->getZPos()));
            //std::cout << ">>> " << sqrt(sum) << " " << prec << std::endl;
            sum = sum + prec * prec;

          }
          m_varMatrixRPhi(r, c) = sum;
          if (r != c) m_varMatrixRPhi(c, r) = sum;
        }
      }
    }
  } // Correlation terms: c is column, r is row

  // Print variance matrix
  //std::cout << "Variance matrix in R-Phi (with zero cols/rows): " << std::endl;
  //printSymMatrix(m_varMatrixRPhi);

  //
  // Remove zero rows and columns in covariance matrix
  int  rActual = -1;          // Row, at which to move the active row due to a sequence of zero rows or inactive hits inbetween
  bool lookForActive = false; // Start looking for shift of active rows, after first passive row found

  for (int r=0; r<nHits; r++) {

    // Keep actual row @ zero for first N passive (zero) layers
    if ((m_hits.at(r)->isPassive() || r<iStart) && (!lookForActive)) {

      // Next hit has to be active (set as active and considered in track fitting (see iStart))
      if ((r+1)<nHits && m_hits.at(r+1)->isActive() && (r+1)>=iStart) lookForActive = true;

      // Previous hit has to be passive or not being considered in track fitting (see iStart))
      if (!((r-1)>=0 && (m_hits.at(r-1)->isPassive() || (r-1)<iStart)) ) rActual = r;
    }
    // Shift active layer to zero-th row + i active layers, which have already been shifted by number of zero layers
    else if ((m_hits.at(r)->isActive()) && (lookForActive)) {
      for (int c=0; c<nHits; c++) {
        m_varMatrixRPhi(rActual, c) = m_varMatrixRPhi(r, c);
        m_varMatrixRPhi(c, rActual) = m_varMatrixRPhi(c, r);
      }

      m_varMatrixRPhi(rActual, rActual) = m_varMatrixRPhi(r, r);
      rActual++;
    }
  }
  // If some rows/colums were zero -> matrix rank needs to be adjusted
  int nResized = rActual;
  if (nResized!=0) m_varMatrixRPhi.ResizeTo(nResized, nResized);

  // Print variance matrix
  //std::cout << "Variance matrix in R-Phi: " << std::endl;
  //printSymMatrix(m_varMatrixRPhi);

  // Check if matrix is sane and worth keeping
  if (!((m_varMatrixRPhi.GetNoElements() > 0) && (m_varMatrixRPhi.Determinant() != 0.0))) {
    logWARNING("Variance matrix V(NxN) in R-Phi -> zero determinat or zero number of elements");
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); //  Important -> resort back
    return false;
  }

  //
  // Compute 3x3 covariance matrix of the track parameters in R-Phi projection
  unsigned int offset = iStart;

  int n = m_varMatrixRPhi.GetNrows();

  TMatrixT<double> V(m_varMatrixRPhi); // Local copy to be inverted
  TMatrixT<double> diffsT(3, n);       // Derivatives of track parameters transposed (in R-Phi -> 3 track parameters)
  TMatrixT<double> diffs(n, 3);        // Derivatives of track parameters (in R-Phi -> 3 track parameters)

  m_covMatrixRPhi.ResizeTo(3, 3);

  // Set up partial derivative matrices diffs and diffsT -> using Karimaki approach & parabolic aproximations to define these matrices
  for (auto i=iStart; i<=iEnd; i++) {

    if (m_hits.at(i)->isActive()) {
      diffs(i - offset, 0) = computeDfOverDRho(m_hits.at(i)->getRPos(),m_hits.at(i)->getZPos());
      diffs(i - offset, 1) = +m_hits.at(i)->getRPos(); // No impact of sign on results, but from analytical derivation point of view correct with a plus sign!!! Was minus sign here!!!
      diffs(i - offset, 2) = 1;
    }
    else offset++;
  }

  // Transpose
  diffsT.Transpose(diffs);

  // Print A matrix
  //std::cout << "A matrix in R-Phi: " << std::endl;
  //printMatrix(diffs);

  // Get covariance matrix using global chi2 fit: C = cov(i,j) = (D^T * V^-1 * D)^-1
  m_covMatrixRPhi = diffsT * V.Invert() * diffs;
  m_covMatrixRPhi.Invert();

  // Sort-back hits based on particle direction if they were resorted
  if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR);

  return true;
}

//
// Compute 2x2 covariance matrix of track parameters in R-Z (s-z), using the NxN variance matrix (N hits = K+L: K active hits on detectors + L passive (artificial) hits due to material).
// Ref. point dictates, whether Multiple scattering effects need to be calculated inside-out or outside-in. MS effect is symmetric as regards track fitting. PropagOutIn variable
// defines whether detectors at higher R than the ref. point (true) should be used for error calculation/propagation (e.g. d0,z0) or whether detectors at lower R (false).
// Return true if V invertable
//
bool TrackNew::computeCovarianceMatrixRZ(double refPointRPos, bool propagOutIn) {

  // Matrix size
  int nHits = m_hits.size();

  m_varMatrixRZ.ResizeTo(nHits,nHits);
  m_varMatrixRZ.Zero();

  //
  // Find hit index ranges relevant for MS effects calculation based on requirements on refPoint & propagation direction
  int iStart          = nHits;
  int iEnd            = nHits-1;
  int nHitsUsed       = 0;
  int nActiveHitsUsed = 0;

  bool bySmallerR = true;
  bool useIP      = false;

  // Particle traverses inside-out with error propagation outside-in -> hits already sorted correctlly
  if (m_pt>=0 && propagOutIn) {

    for (auto it=m_hits.rbegin(); it!=m_hits.rend(); it++) {

      if (refPointRPos<(*it)->getRPos()) {
        if ((*it)->isMeasurable()) nActiveHitsUsed++;
        iStart--;
      }
      else break;
    }
  }
  // Particle traverses inside-out with error propagation inside-out -> sort hits by higher radius & resort backwards after matrix calculated
  else if (m_pt>=0 && !propagOutIn) {

    sortHits(!bySmallerR);
    for (auto it=m_hits.rbegin(); it!=m_hits.rend(); it++) {

      if (refPointRPos>(*it)->getRPos()) {
        if ((*it)->isMeasurable()) nActiveHitsUsed++;
        iStart--;
      }
      else break;
    }
  }
  // TODO: Implement
  else {

    logWARNING("Variance matrix V(NxN) in R-Z (s-Z) -> calculations of particle traversing outside-in not yet implemented.");
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); // Important -> resort back!
    return false;
  }

  // Check that enough hits for RZ calculation, i.e. >=2
  nHitsUsed = nHits - iStart;
  if ((nActiveHitsUsed)<2) {

    std::string message = "Variance matrix V(NxN) in R-Z (s-Z) -> refPointRPos[mm]="+any2str(refPointRPos/Units::mm,1);
                message+= " in combination with propagator direction doesn't provide sufficient number of hits for tracking!";
    logWARNING(message);
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); // Important -> resort back!
    return false;
  }

  // Pre-compute the squares of the scattering angles
  // already divided by sin^2 (that is : we should use p instead of p_T here
  // but the result for theta^2 differ by a factor 1/sin^2, which is exactly the
  // needed factor to project the scattering angle on an horizontal surface
  std::vector<double> msThetaOverSinSq;

  for (int i=iStart; i<iEnd; i++) {

    // MS theta
    double msTheta = 0.0;

    // Material in terms of rad. lengths
    double XtoX0 = m_hits.at(i)->getCorrectedMaterial().radiation;

    if (XtoX0>0) {
      // MS error depends on path length = deltaR/sin(theta), so one can precalculate msTheta_real as msTheta/sin^2(theta), which practically means using pT
      // instead of p & then one has to multiply the msTheta by deltaR to get MS error
      msTheta = (13.6*Units::MeV * 13.6*Units::MeV) / (m_pt/Units::MeV * m_pt/Units::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));

      // Take into account a propagation of MS error on virtual barrel plane, on which all measurements are evaluated for consistency (global chi2 fit applied) ->
      // in limit R->inf. propagation along line used, otherwise a very small correction factor coming from the circular shape of particle track is required (similar
      // approach as for local resolutions)
      // TODO: Currently, correction mathematicaly derived only for use case of const magnetic field -> more complex mathematical expression expected in non-const B field
      // (hence correction not applied in such case)
      double A = 0;
      if (SimParms::getInstance().isMagFieldConst()) A = m_hits.at(i)->getRPos()/2./getRadius(m_hits.at(i)->getZPos());  // r_i/2R
      double corrFactor = pow( cos(m_theta)*cos(m_theta)/sin(m_theta)/sqrt(1-A*A) + sin(m_theta) ,2); // Without correction it would be 1/sin(theta)^2

      msTheta *=corrFactor;
    }
    else {
      msTheta = 0;
    }
    msThetaOverSinSq.push_back(msTheta);
  }

  //
  // Calculate the variance matrix with all correlation terms: c is column, r is row (hits are assumed to be sorted)
  for (int c=iStart ; c<=iEnd; c++) {

    // Dummy value for correlations involving inactive surfaces
    if (m_hits.at(c)->isPassive()) {
      for (int r=0; r<=c; r++) m_varMatrixRZ(r, c) = 0.0;
    }
    // One of the correlation factors refers to an active surface
    else {

      for (int r=iStart; r<=c; r++) {
        // Dummy value for correlation involving an inactive surface
        if (m_hits.at(r)->isPassive()) m_varMatrixRZ(r, c) = 0.0;

        // Correlations between two active surfaces
        else {

          double sum = 0.0;

          for (int i=iStart; i<r; i++) sum += msThetaOverSinSq.at(i-(nHits-nHitsUsed))
                                            * (m_hits.at(c)->getRPos() - m_hits.at(i)->getRPos())
                                            * (m_hits.at(r)->getRPos() - m_hits.at(i)->getRPos());

          if (r == c) {
            double prec = m_hits.at(r)->getResolutionZ(getRadius(m_hits.at(r)->getZPos()));
            sum = sum + prec * prec;
          }

          m_varMatrixRZ(r, c) = sum;
          if (r != c) m_varMatrixRZ(c, r) = sum;
#undef CORRELATIONS_OFF_DEBUG
#ifdef CORRELATIONS_OFF_DEBUG
          if (r!=c) {
            m_varMatrixRZ(c, r)=0;
            m_varMatrixRZ(r, c)=0;
          }
#endif
        }
      }
    }
  } // Calculate correlation terms: c is column, r is row

  // Print variance matrix
  //std::cout << "Variance matrix in R-Z (with zero cols/rows): " << std::endl;
  //printSymMatrix(m_varMatrixRZ);

  //
  // Remove zero rows and columns in covariance matrix
  int  rActual = -1;          // Row, at which to move the active row due to a sequence of zero rows or inactive hits inbetween
  bool lookForActive = false; // Start looking for shift of active rows, after first passive row found

  for (int r=0; r<nHits; r++) {

    // Keep actual row @ zero for first N passive (zero) layers
    if ((m_hits.at(r)->isPassive() || r<iStart) && (!lookForActive)) {

      // Next hit has to be active (set as active and considered in track fitting (see iStart))
      if ((r+1)<nHits && m_hits.at(r+1)->isActive() && (r+1)>=iStart) lookForActive = true;

      // Previous hit has to be passive or not being considered in track fitting (see iStart))
      if (!((r-1)>=0 && (m_hits.at(r-1)->isPassive() || (r-1)<iStart)) ) rActual = r;
    }
    // Shift active layer to zero-th row + i active layers, which have already been shifted by number of zero layers
    else if ((m_hits.at(r)->isActive()) && (lookForActive)) {
      for (int c=0; c<nHits; c++) {
        m_varMatrixRZ(rActual, c) = m_varMatrixRZ(r, c);
        m_varMatrixRZ(c, rActual) = m_varMatrixRZ(c, r);
      }

      m_varMatrixRZ(rActual, rActual) = m_varMatrixRZ(r, r);
      rActual++;
    }
  }
  // If some rows/colums were zero -> matrix rank needs to be adjusted
  int nResized = rActual;
  if (nResized!=0) m_varMatrixRZ.ResizeTo(nResized, nResized);

  // Print variance matrix
  //std::cout << "Variance matrix in R-Z: " << std::endl;
  //printSymMatrix(m_varMatrixRZ);

  // Check if matrix is sane and worth keeping
  if (!((m_varMatrixRZ.GetNoElements() > 0) && (m_varMatrixRZ.Determinant() != 0.0))) {
    logWARNING("Variance matrix V(NxN) in R-Z (s-Z) -> zero determinat or zero number of elements");
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); // Important - resort back!
    return false;
  }

  //
  // Compute 2x2 covariance matrix of the track parameters in R-Z (s-Z) projection
  unsigned int offset       = iStart;
  unsigned int varMatrixDim = m_varMatrixRZ.GetNrows();

  TMatrixT<double> V(m_varMatrixRZ);        // Local copy to be inverted
  TMatrixT<double> diffsT(2, varMatrixDim); // Derivatives of track parameters transposed (in R-Z -> 2 track parameters)
  TMatrixT<double> diffs(varMatrixDim, 2);  // Derivatives of track parameters (in R-Z -> 2 track parameters)

  m_covMatrixRZ.ResizeTo(2,2);

  // Set up partial derivative matrices diffs and diffsT -> line fit in s-Z to define these matrices
  for (auto i=iStart; i<=iEnd; i++) {

    if (m_hits.at(i)->isActive()) {

      // Partial derivatives for x = p[0] * y + p[1]
      diffs(i - offset, 0) = m_hits.at(i)->getRPos();
      diffs(i - offset, 1) = 1;
    }
    else offset++;
  }

  // Transpose
  diffsT.Transpose(diffs);

  // Print
  //std::cout << "Diff matrix in R-Z: " << std::endl;
  //printMatrix(diffsT);

  // Get covariance matrix using global chi2 fit: C = cov(i,j) = (D^T * V^-1 * D)^-1
  m_covMatrixRZ = diffsT * V.Invert() * diffs;
  m_covMatrixRZ.Invert();

  // Sort-back hits based on particle direction if they were resorted
  if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR);

  return true;
}

//
// Helper fce returning derivative: df(rho, d0, phi0)/drho, where f approximates
// a helix by set of parabolas. In general, N connected parabolas used, for const B
// field only one parabola applied.
//
double TrackNew::computeDfOverDRho(double rPos, double zPos) {

  double DfOverDRho = 0;

  // Option 1: Const mag. field across the detector: Bz = const -> for now assumed to be const.
  if (SimParms::getInstance().isMagFieldConst()) {

    DfOverDRho = 0.5 * rPos*rPos;
  }
//  // Option 2: Mag. field is a function in Z: B = B(z)
//  else {
//
//    int nRegions = SimParms::getInstance().getNMagFieldRegions();
//
//    // Find i-th region corresponding to the current zPos
//    int iRegion = 0;
//    for (iRegion=0; iRegion < nRegions; iRegion++) {
//
//      if (zPos<(SimParms::getInstance().magFieldZRegions[iRegion])) break;
//    }
//
//    // Check that zPos not beyond Z-range, in which B field has been defined
//    if (iRegion==nRegions) {
//
//      std::ostringstream message;
//      message << "Track::computeDfOverDRho(): Hit z-position: " << zPos/Units::mm << " beyond defined B field Z-range: [0," << SimParms::getInstance().magFieldZRegions[nRegions-1]/Units::mm << "]!";
//      logERROR(message.str());
//      exit(1);
//    }
//
//    // Z pos. in the first region or only 1 region defined (const mag. field)
//    if (iRegion==0) {
//      DfOverDRho = 0.5 * rPos*rPos;
//    }
//    // Z pos in i-th region (generally N regions defined)
//    else {
//
//      // Get reference magnetic field B0 (at [r,z] = [0,0]
//      double B0   = SimParms::getInstance().magField[0];
//
//      double Bi   = 0.; // B-field in ith z-region
//      double Bi_1 = 0.; // B-field in (i-1)th z-region
//      double xi   = 0.; // x-position corresponding to the ith z-region
//
//      // Sum-up all contributions across the regions: 0th - ith
//      for (int i=1; i<=iRegion; i++) {
//
//        // Get current value of magnetic field B_i & B_i-1
//        Bi   = SimParms::getInstance().magField[i];
//        Bi_1 = SimParms::getInstance().magField[i-1];
//        xi   = SimParms::getInstance().magFieldZRegions[i-1] * tan(m_theta); // (z0,z1,z2...) -> intervals defined as 0-z0, z0-z1, z1-z2
//
//        // Add dB/dz terms
//        DfOverDRho+= -1.0*(Bi-Bi_1)/B0 * xi*rPos;
//        DfOverDRho+= +0.5*(Bi-Bi_1)/B0 * xi*xi;
//
//        // Add Bi/B0 term
//        if (i==iRegion) DfOverDRho += +0.5*Bi/B0 * rPos*rPos;
//      }
//    }
//  }

  return DfOverDRho;
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
//    // if (debugRemoval) std::cerr << "r = " << (*it)->getRadius((*it)->getZPos()) << ", ";
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

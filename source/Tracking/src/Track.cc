/**
 * @file Track.cpp
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "Track.h"

#include <algorithm>
#include <cstdlib>

#include <global_constants.h>
#include "Hit.h"
#include "MessageLogger.h"
#include "MaterialProperties.h"
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
  m_pt(0)
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

  m_origin       = track.m_origin;
  m_direction    = track.m_direction;

  m_varMatrixRPhi.ResizeTo(track.m_varMatrixRPhi);
  m_varMatrixRPhi = track.m_varMatrixRPhi;
  m_covMatrixRPhi.ResizeTo(track.m_covMatrixRPhi);
  m_covMatrixRPhi = track.m_covMatrixRPhi;

  m_varMatrixRZ.ResizeTo(track.m_varMatrixRZ);
  m_varMatrixRZ = track.m_varMatrixRZ;
  m_covMatrixRZ.ResizeTo(track.m_covMatrixRZ);
  m_covMatrixRZ = track.m_covMatrixRZ;

  for (auto& iHit : track.m_hits) {
    HitPtr hit(new Hit(*iHit));
    addHit(std::move(hit));
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

  m_origin       = track.m_origin;
  m_direction    = track.m_direction;

  m_varMatrixRPhi.ResizeTo(track.m_varMatrixRPhi);
  m_varMatrixRPhi = track.m_varMatrixRPhi;
  m_covMatrixRPhi.ResizeTo(track.m_covMatrixRPhi);
  m_covMatrixRPhi = track.m_covMatrixRPhi;

  m_varMatrixRZ.ResizeTo(track.m_varMatrixRZ);
  m_varMatrixRZ = track.m_varMatrixRZ;
  m_covMatrixRZ.ResizeTo(track.m_covMatrixRZ);
  m_covMatrixRZ = track.m_covMatrixRZ;

  for (auto& iHit : track.m_hits) {
    HitPtr hit(new Hit(*iHit));
    addHit(std::move(hit));
  }
  m_tags = track.m_tags;

  // Return the existing object
  return *this;
}

//
// Destructor
//
Track::~Track() {

  // Clear memory
  m_hits.clear();
}

//
// Calculate magnetic field at given z, assuming B = B(z).e_z + 0.e_x + 0 e_y
//
double Track::getMagField(double z) const {

  double magField = 0;

  // Option 1: Const mag. field across the detector: Bz = const
  if (SimParms::getInstance().isMagFieldConst()) {

    magField = SimParms::getInstance().magField[0];
  }
  // Option 2: Mag. field is a function in Z: B = B(z)
  else {

    for (unsigned int i=0; i<SimParms::getInstance().getNMagFieldRegions(); i++) {

      // Magnetic field regions are considered to be defined in metres
      if (z<SimParms::getInstance().magFieldZRegions[i]) {

        // Magnetic field is considered to be in Tesla
        magField = SimParms::getInstance().magField[i];
        break;
      }
    }
  }

  return magField;
}

//
// Main method calculating track parameters using Karimaki approach & parabolic approximation in R-Phi plane: 1/R, D0, phi parameters
// and using linear fit in s-Z plane: cotg(theta), Z0 parameters
//
void Track::computeErrors() {

  // Compute the relevant 2x2 covariance matrix in RZ plane first (check that V matrix can be inverted)
  if (computeVarianceMatrixRZ()) computeCovarianceMatrixRZ();

  // Compute the relevant 3x3 covariance matrix in R-Phi plane (check that V matrix can be inverted)
  if (computeVarianceMatrixRPhi()) computeCovarianceMatrixRPhi();
}

//
// Get DeltaRho (error on 1/R) at path length s projected to XY plane, i.e. at [r,z] sXY ~ r
// Using 3x3 covariance propagator in case [r,z]!=[0,0]
//
double Track::getDeltaRho(double rPos) const {

  double deltaRho = -1.;
  if (m_covMatrixRPhi(0, 0)>=0) deltaRho = sqrt(m_covMatrixRPhi(0, 0));

  // TODO: Not working correctly for B = B(z) & [r,z]!=[0,0], so print warning ...
  if (rPos!=0. && !SimParms::getInstance().isMagFieldConst()) {

    logWARNING("Track::getDeltaRho(): Mathematical method to get deltaRho at [r,z]!=[0,0] in non const. B field not implemented, hence returned value at [r,z]=[0,0].");
  }

  return deltaRho;
}

//
// Get DeltaPtOvePt at path length s projected to XY plane, i.e. at [r,z] sXY ~ r (utilize the calculated deltaRho quantity)
//
double Track::getDeltaPtOverPt(double rPos) const {

  double deltaPtOverPt = -1.;

  // delta(1/R) & delta(pT) -> estimated at point [r,z] = [0,0] (important for use case, when B != const -> B = B(z))
  double deltaRho = getDeltaRho(rPos);
  double radius   = getRadius(rPos*m_cotgTheta);       // Approximative transformation from rPos to zPos using tan(theta)
  if (deltaRho!=-1) deltaPtOverPt = deltaRho * radius; // dpT(z)/pT(z) = dRho(z) / Rho(z) = dRho(z) * R(z)

  return deltaPtOverPt;
}

//
// Get DeltaPOverP at path length s projected to XY plane, i.e. at [r,z] sXY ~ r (utilize deltaRho & deltaCotgTheta quantities)
//
double Track::getDeltaPOverP(double rPos) const {

  double deltaPOverP = -1.;

  // Combining into p measurement
  // dp/p = dp_t/p_t + A / (1+A^2) * dA // with A = ctg(theta)
  // dp/p = dp_t/p_t + sin(theta)*cos(theta)
  double deltaPtOverPt = getDeltaPtOverPt(rPos);
  double deltaCtgTheta = getDeltaCtgTheta();
  if (deltaPtOverPt!=-1 && deltaCtgTheta!=-1) deltaPOverP = sqrt(deltaPtOverPt*deltaPtOverPt + sin(m_theta)*sin(m_theta) * cos(m_theta)*cos(m_theta) * deltaCtgTheta*deltaCtgTheta);

  return deltaPOverP;
}

//
// Get DeltaPhi0 at point (r,z) at path length s projected to XY plane, i.e. at [r,z] sXY ~ r
// Using 3x3 covariance propagator in case [r,z]!=[0,0]
//
double Track::getDeltaPhi0(double rPos) const {

  double deltaPhi0 = -1.;
  if (m_covMatrixRPhi(1, 1) >= 0) deltaPhi0 = sqrt(m_covMatrixRPhi(1, 1));

  // No covariance propagation necessary at [0,0,0] point
  if (rPos==0.) return deltaPhi0;
  else {

    double covRhoRho    = m_covMatrixRPhi(0,0);
    double covRhoPhi0   = m_covMatrixRPhi(0,1);
    double covRhoD0     = m_covMatrixRPhi(0,2);
    double covPhi0Phi0  = m_covMatrixRPhi(1,1);
    double covPhi0D0    = m_covMatrixRPhi(1,2);
    double covD0D0      = m_covMatrixRPhi(2,2);
    double rho          = getRho(rPos*m_cotgTheta);

    double deltaPhi0Sq  = rPos*rPos*covRhoRho + covPhi0Phi0                   + rho*rho*rho*rho*rPos*rPos*covD0D0;
           deltaPhi0Sq += 2*rPos*covRhoPhi0   - 2*rho*rho*rPos*rPos*covRhoD0  - 2*rho*rho*rPos*covPhi0D0;

    if (deltaPhi0Sq>=0) return sqrt(deltaPhi0Sq);
    else                return -1.;
  }

  // TODO: Not working correctly for B = B(z) & [r,z]!=[0,0], so print warning ...
  if (rPos!=0. && !SimParms::getInstance().isMagFieldConst()) {

    logWARNING("Track::getDeltaPhi(): Mathematical method to get deltaPhi at [r,z]!=[0,0] in non const. B field not implemented, hence returned value at [r,z]=[0,0].");
  }

  return deltaPhi0;
}

//
// Get DeltaD0 at path length s projected to XY plane, i.e. at [r,z] sXY ~ r
// Using 3x3 covariance propagator in case [r,z]!=[0,0]
//
double Track::getDeltaD0(double rPos) const {

  double deltaD0 = -1.;
  if (m_covMatrixRPhi(2, 2)) deltaD0 = sqrt(m_covMatrixRPhi(2, 2));

  // No covariance propagation necessary at [0,0,0] point
  if (rPos==0.) return deltaD0;
  else {

    double covPhi0Phi0  = m_covMatrixRPhi(1,1);
    double covPhi0D0    = m_covMatrixRPhi(1,2);
    double covD0D0      = m_covMatrixRPhi(2,2);

    double deltaD0Sq  = rPos*rPos*covPhi0Phi0 + 2*rPos*covPhi0D0 + covD0D0;

    if (deltaD0Sq>=0) return sqrt(deltaD0Sq);
    else              return -1.;
  }

  // TODO: Not working correctly for B = B(z) & [r,z]!=[0,0], so print warning ...
  if (rPos!=0. && !SimParms::getInstance().isMagFieldConst()) {

    logWARNING("Track::getDeltaD0(): Mathematical method to get deltaD0 at [r,z]!=[0,0] in non const. B field not implemented, hence returned value at [r,z]=[0,0].");
  }

  return deltaD0;
}

//
// Get DeltaCtgTheta at path length s projected to XY plane, i.e. at [r,z] sXY ~ r (independent on sXY)
//
double Track::getDeltaCtgTheta() const {

  double deltaCtgTheta = -1.;
  if (m_covMatrixRZ(0, 0)>=0) deltaCtgTheta = sqrt(m_covMatrixRZ(0, 0));

  return deltaCtgTheta;
}

//
// Get DeltaZ0 at path length s projected to XY plane, i.e. at [r,z] sXY ~ r
// Using 2x2 covariance propagator in case [r,z]!=[0,0]
//
double Track::getDeltaZ0(double rPos) const {

  double deltaZ0 = -1.;
  if (m_covMatrixRZ(1, 1)>=0) deltaZ0 = sqrt(m_covMatrixRZ(1, 1));

  // No covariance propagation necessary at [0,0,0] point
  if (rPos==0.) return deltaZ0;
  else {

    double covZ0Z0       = m_covMatrixRZ(1,1);
    double covZ0CtgTheta = m_covMatrixRZ(0,1);
    double covCtgThCtgTh = m_covMatrixRZ(0,0);

    double deltaZ0Sq = covZ0Z0 + 2*rPos*covZ0CtgTheta + rPos*rPos*covCtgThCtgTh;
    if (deltaZ0Sq>=0) return sqrt(deltaZ0Sq);
    else              return -1.;
  }
}

//
// Adds a new hit to the track (hit radius automatically updated)
//
void Track::addHit(HitPtr newHit) {

  // Add tracking tags
  if (newHit->getHitModule() != nullptr) {
    m_tags.insert(newHit->getHitModule()->trackingTags.begin(), newHit->getHitModule()->trackingTags.end());
  }
  newHit->setTrack(this);

  // Add new hit
  m_hits.push_back(std::move(newHit));
}

//
// Add IP constraint to the track, technically new hit is assigned: with no material and hit resolution in R-Phi as dr, in s-Z as dz
//
void Track::addIPConstraint(double dr, double dz) {

  // This modeling of the IP constraint was validated:
  // By placing dr = 0.5 mm and dz = 1 mm one obtains
  // sigma(d0) = 0.5 mm and sigma(z0) = 1 mm
  HitPtr newHit(new Hit(dr,dz)); // TODO: Cross-check, should be Hit(0,0) ???
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
  this->addHit(std::move(newHit));

}

//
// Sort internally all hits assigned to this track -> sorting algorithm based on hit radius - by smaller radius sooner or vice-versa (inner-2-outer approach or vice-versa)
//
void Track::sortHits(bool bySmallerR) { bySmallerR ? std::stable_sort(m_hits.begin(), m_hits.end(), Hit::sortSmallerR) : std::stable_sort(m_hits.begin(), m_hits.end(), Hit::sortHigherR); }

//
// Remove hits that don't follow the parabolic approximation used in tracking - TODO: still needs to be updated (not all approximations taken into account)
//
bool Track::pruneHits() {

  bool isPruned = false;

  HitCollection newHits;
  for (auto& iHit : m_hits) {

    if ( iHit->getRPos()<2*getRadius(iHit->getZPos()) ) newHits.push_back(std::move(iHit));
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
void Track::keepTaggedOnly(const string& tag) {

  for (auto& iHit : m_hits) {

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
  for (auto& iHit : m_hits) iHit->setCorrectedMaterial(nullMaterial);
}

//
// Helper method printing track covariance matrices in R-Phi
//
void Track::printErrors() {

  std::cout << "Overview of track errors:" << std::endl;
  std::cout << "Hit variance matrix: "  << std::endl;
  m_varMatrixRPhi.Print();

  std::cout << "Covariance matrix: " << std::endl;
  m_covMatrixRPhi.Print();

  // Print errors @ [r,z]=[0,0]
  double rPos = 0.0;

  std::cout << "Rho errors by momentum: " << getDeltaRho(rPos) << std::endl;
  std::cout << "Phi0 errors by momentum: "<< getDeltaPhi0(rPos)<< std::endl;
  std::cout << "D0 errors by momentum: "  << getDeltaD0(rPos)  << std::endl;
}

//
// Helper method printing symmetric matrix
//
void Track::printSymMatrix(const TMatrixTSym<double>& matrix) {

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
void Track::printMatrix(const TMatrixT<double>& matrix) {

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
void Track::printHits() {

  std::cout << "******************" << std::endl;
  std::cout << "Track eta=" << m_eta << std::endl;

  for (const auto& it : m_hits) {
    std::cout << "    Hit"
              << " r="  << it->getRPos()
              << " z="  << it->getZPos()
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

  m_direction.SetCoordinates(1, m_theta, m_phi);

  return m_direction;
};

//
// Get number of active hits assigned to track for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file). If tag specified as "all" no extra tag required
//
int Track::getNActiveHits (std::string tag, bool useIP /* = true */ ) const {

  // Result variable
  int nHits=0;

  for (auto& iHit : m_hits) {
    if (iHit) {
      if ((useIP) || (!iHit->isIP())) {
        if (iHit->getObjectKind()==HitKind::Active){

          // Check tag
          bool tagOK = false;
          for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
            if (tag==*it || tag=="all") tagOK = true;
          }

          if (tagOK) nHits++;
        }
      }
    }
  } // For

  return nHits;
}

//
// Get the probabilty of having "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
// If tag specified as "all" no extra tag required
//
std::vector<double> Track::getHadronActiveHitsProbability(std::string tag) {

  // Result variable
  std::vector<double> probabilities;
  double probability = 1;

  // Sort hits first
  bool bySmallerR = true;
  sortHits(bySmallerR);

  for (auto& iHit : m_hits) {
    if (iHit) {
      if (iHit->getObjectKind()==HitKind::Active){

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
double Track::getHadronActiveHitsProbability(std::string tag, int nHits) {

  // Probability
  double probability = 1;

  // Number of clean hits
  int goodHits = 0;

  // Sort hits first
  bool bySmallerR = true;
  sortHits(bySmallerR);

  for (auto& iHit : m_hits) {

    if (iHit) {
      if (iHit->getObjectKind()==HitKind::Active) {

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
RILength Track::getMaterial() {

  RILength totalMaterial;
  totalMaterial.radiation   = 0;
  totalMaterial.interaction = 0;

  for (auto& iHit : m_hits) totalMaterial += iHit->getCorrectedMaterial();

  return totalMaterial;
}

//
// Get a vector of pairs: Detector module & hit type for Trigger hits
//
std::vector<std::pair<const DetectorModule*, HitType>> Track::getHitModules() const {

  std::vector<std::pair<const DetectorModule*, HitType>> result;

  for (auto& iHit : m_hits) {

    if ((iHit) && (iHit->isTrigger()) && (!iHit->isIP()) && (iHit->getObjectKind()==HitKind::Active)) {

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
// Compute the variance matrix in R-Phi: NxN (N hits = K+L: K active hits on detectors + L passive (artificial) hits due to material)
//
bool Track::computeVarianceMatrixRPhi() {

  // Variance matrix size
  int n = m_hits.size();
  m_varMatrixRPhi.ResizeTo(n,n);

  // Get contributions from Multiple Couloumb scattering
  std::vector<double> msThetaOverSinSq;

  for (int i = 0; i < n - 1; i++) {

    // MS theta
    double msTheta = 0.0;

    // Material in terms of rad. lengths
    double XtoX0 = m_hits.at(i)->getCorrectedMaterial().radiation;
    //std::cout << std::fixed << std::setprecision(4) << "Material (" << i << ") = " << XtoX0 << "\t at r=" << m_hits.at(i)->getRadius(m_hits.at(i)->getZPos()) << "\t of type=" << m_hits.at(i)->getObjectKind() << std::endl;

    if (XtoX0>0) {

      // MS error depends on path length = deltaR/sin(theta), so one can precalculate msTheta_real as msTheta/sin^2(theta), which practically means using pT
      // instead of p & then one has to multiply the msTheta by deltaR to get MS error
      msTheta = (13.6*Units::MeV * 13.6*Units::MeV) / (m_pt/Units::MeV * m_pt/Units::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));

      // Take into account a very small correction factor coming from the circular shape of particle track (similar approach as for local resolutions)
      double A = m_hits.at(i)->getRPos()/2./getRadius(m_hits.at(i)->getZPos());     // r_i/2R
           A = 0; // For testing purposes only -> don't apply helix correction
      double corrFactor = 1 + A*A*cos(m_theta)*cos(m_theta)/(1-A*A);

      msTheta *= corrFactor;
    }
    else {
      msTheta = 0;
    }
    msThetaOverSinSq.push_back(msTheta);
  }

  // Calculate correlation terms: c is column, r is row (hits are assumed to be sorted)
  for (int c = 0; c < n; c++) {

    // Dummy value for correlations involving inactive surfaces
    if (m_hits.at(c)->getObjectKind() == HitKind::Inactive) {
      for (int r = 0; r <= c; r++) m_varMatrixRPhi(r, c) = 0.0;
    }
    // One of the correlation factors refers to an active surface
    else {

      for (int r = 0; r <= c; r++) {
        // Dummy value for correlation involving an inactive surface
        if (m_hits.at(r)->getObjectKind() == HitKind::Inactive) m_varMatrixRPhi(r, c) = 0.0;

        // Correlations between two active surfaces
        else {

          double sum = 0.0;

          for (int i = 0; i < r; i++) {

            sum += msThetaOverSinSq.at(i)
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
        m_varMatrixRPhi(nResized, j) = m_varMatrixRPhi(i, j);
        m_varMatrixRPhi(j, nResized) = m_varMatrixRPhi(j, i);
      }
      m_varMatrixRPhi(nResized, nResized) = m_varMatrixRPhi(i, i);
      nResized++;
    }
  }

  // Resize matrix if necessary
  if (nResized != -1) m_varMatrixRPhi.ResizeTo(nResized, nResized);

  // Print
  //std::cout << "Variance matrix in R-Phi: " << std::endl;
  //printSymMatrix(m_varMatrixRPhi);

  // Check if matrix is sane and worth keeping
  if (!((m_varMatrixRPhi.GetNoElements() > 0) && (m_varMatrixRPhi.Determinant() != 0.0))) {
    logWARNING("Variance matrix V(NxN) in R-Phi -> zero determinat or zero number of elements");
    return false;
  }
  else return true;
}

//
// Helper fce returning derivative: df(rho, d0, phi0)/drho, where f approximates
// a helix by set of parabolas. In general, N connected parabolas used, for const B
// field only one parabola applied.
//
double Track::computeDfOverDRho(double rPos, double zPos) {

  double DfOverDRho = 0;

  // Option 1: Const mag. field across the detector: Bz = const
  if (SimParms::getInstance().isMagFieldConst()) {

    DfOverDRho = 0.5 * rPos*rPos;
  }
  // Option 2: Mag. field is a function in Z: B = B(z)
  else {

    int nRegions = SimParms::getInstance().getNMagFieldRegions();

    // Find i-th region corresponding to the current zPos
    int iRegion = 0;
    for (iRegion=0; iRegion < nRegions; iRegion++) {

      if (zPos<(SimParms::getInstance().magFieldZRegions[iRegion])) break;
    }

    // Check that zPos not beyond Z-range, in which B field has been defined
    if (iRegion==nRegions) {

      std::ostringstream message;
      message << "Track::computeDfOverDRho(): Hit z-position: " << zPos/Units::mm << " beyond defined B field Z-range: [0," << SimParms::getInstance().magFieldZRegions[nRegions-1]/Units::mm << "]!";
      logERROR(message.str());
      exit(1);
    }

    // Z pos. in the first region or only 1 region defined (const mag. field)
    if (iRegion==0) {
      DfOverDRho = 0.5 * rPos*rPos;
    }
    // Z pos in i-th region (generally N regions defined)
    else {

      // Get reference magnetic field B0 (at [r,z] = [0,0]
      double B0   = SimParms::getInstance().magField[0];

      double Bi   = 0.; // B-field in ith z-region
      double Bi_1 = 0.; // B-field in (i-1)th z-region
      double xi   = 0.; // x-position corresponding to the ith z-region

      // Sum-up all contributions across the regions: 0th - ith
      for (int i=1; i<=iRegion; i++) {

        // Get current value of magnetic field B_i & B_i-1
        Bi   = SimParms::getInstance().magField[i];
        Bi_1 = SimParms::getInstance().magField[i-1];
        xi   = SimParms::getInstance().magFieldZRegions[i-1] * tan(m_theta); // (z0,z1,z2...) -> intervals defined as 0-z0, z0-z1, z1-z2

        // Add dB/dz terms
        DfOverDRho+= -1.0*(Bi-Bi_1)/B0 * xi*rPos;
        DfOverDRho+= +0.5*(Bi-Bi_1)/B0 * xi*xi;

        // Add Bi/B0 term
        if (i==iRegion) DfOverDRho += +0.5*Bi/B0 * rPos*rPos;
      }
    }
  }

  return DfOverDRho;
}

//
// Compute 3x3 covariance matrix of the track parameters in R-Phi projection
//
void Track::computeCovarianceMatrixRPhi() {

  unsigned int offset = 0;
  unsigned int nHits  = m_hits.size();

  int n = m_varMatrixRPhi.GetNrows();

  TMatrixT<double> V(m_varMatrixRPhi); // Local copy to be inverted
  TMatrixT<double> diffsT(3, n);       // Derivatives of track parameters transposed (in R-Phi -> 3 track parameters)
  TMatrixT<double> diffs(n, 3);        // Derivatives of track parameters (in R-Phi -> 3 track parameters)

  m_covMatrixRPhi.ResizeTo(3, 3);

  // Set up partial derivative matrices diffs and diffsT -> using Karimaki approach & parabolic aproximations to define these matrices
  for (auto i = 0; i < nHits; i++) {

    if (m_hits.at(i)->getObjectKind()  == HitKind::Active) {
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
}

//
// Compute the variance matrix in R-Z (s-Z): NxN (N hits = K+L: K active hits on detectors + L passive (artificial) hits due to material)
//
bool Track::computeVarianceMatrixRZ() {

  // Matrix size
  int n = m_hits.size();
  m_varMatrixRZ.ResizeTo(n,n);

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
      // MS error depends on path length = deltaR/sin(theta), so one can precalculate msTheta_real as msTheta/sin^2(theta), which practically means using pT
      // instead of p & then one has to multiply the msTheta by deltaR to get MS error
      msTheta = (13.6*Units::MeV * 13.6*Units::MeV) / (m_pt/Units::MeV * m_pt/Units::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));

      // Take into account a propagation of MS error on virtual barrel plane, on which all measurements are evaluated for consistency (global chi2 fit) ->
      // in limit R->inf. propagation along line, otherwise very small correction factor coming from the circular shape of particle track is used (similar
      // approach as for local resolutions)
      double A = m_hits.at(i)->getRPos()/2./getRadius(m_hits.at(i)->getZPos());  // r_i/2R
           A = 0; // For testing purposes only -> don't apply helix correction
      double corrFactor = pow( cos(m_theta)*cos(m_theta)/sin(m_theta)/sqrt(1-A*A) + sin(m_theta) ,2); // Without correction it would be 1/sin(theta)^2

      msTheta *=corrFactor;
    }
    else {
      msTheta = 0;
    }
    msThetaOverSinSq.push_back(msTheta);
  }

  // Calculate correlation terms: c is column, r is row (hits are assumed to be sorted)
  for (int c = 0; c < n; c++) {

    // Dummy value for correlations involving inactive surfaces
    if (m_hits.at(c)->getObjectKind() == HitKind::Inactive) {
      for (int r = 0; r <= c; r++) m_varMatrixRZ(r, c) = 0.0;
    }
    // One of the correlation factors refers to an active surface
    else {

      for (int r = 0; r <= c; r++) {
        // Dummy value for correlation involving an inactive surface
        if (m_hits.at(r)->getObjectKind() == HitKind::Inactive) m_varMatrixRZ(r, c) = 0.0;

        // Correlations between two active surfaces
        else {

          double sum = 0.0;

          for (int i = 0; i < r; i++) sum += msThetaOverSinSq.at(i)
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
        m_varMatrixRZ(nResized, j) = m_varMatrixRZ(i, j);
        m_varMatrixRZ(j, nResized) = m_varMatrixRZ(j, i);
      }

      m_varMatrixRZ(nResized, nResized) = m_varMatrixRZ(i, i);
      nResized++;
    }
  }
  
  // Resize matrix if necessary
  if (nResized!=-1) m_varMatrixRZ.ResizeTo(nResized, nResized);

  // Check if matrix is sane and worth keeping
  if (!((m_varMatrixRZ.GetNoElements() > 0) && (m_varMatrixRZ.Determinant() != 0.0))) {
    logWARNING("Variance matrix V(NxN) in R-Z (s-Z) -> zero determinat or zero number of elements");
    return false;
  }
  else return true;
}

//
// Compute 2x2 covariance matrix of the track parameters in R-Z (s-Z) projection
//
void Track::computeCovarianceMatrixRZ() {

  unsigned int offset = 0;
  unsigned int nHits  = m_hits.size();

  int n = m_varMatrixRZ.GetNrows();

  TMatrixT<double> V(m_varMatrixRZ); // Local copy to be inverted
  TMatrixT<double> diffsT(2, n);     // Derivatives of track parameters transposed (in R-Z -> 2 track parameters)
  TMatrixT<double> diffs(n, 2);      // Derivatives of track parameters (in R-Z -> 2 track parameters)

  m_covMatrixRZ.ResizeTo(2,2);
  
  // Set up partial derivative matrices diffs and diffsT -> line fit in s-Z to define these matrices
  for (auto i = 0; i < nHits; i++) {

    if (m_hits.at(i)->getObjectKind()  == HitKind::Active) {

      // Partial derivatives for x = p[0] * y + p[1]
      diffs(i - offset, 0) = m_hits.at(i)->getRPos();
      diffs(i - offset, 1) = 1;
    }
    else offset++;
  }

  // Transpose
  diffsT.Transpose(diffs);

  // Get covariance matrix using global chi2 fit: C = cov(i,j) = (D^T * V^-1 * D)^-1
  m_covMatrixRZ = diffsT * V.Invert() * diffs;
  m_covMatrixRZ.Invert();
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

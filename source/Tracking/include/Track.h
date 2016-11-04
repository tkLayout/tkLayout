/**
 * @file Track.h
 * @brief This header file defines track class used for internal analysis
 */

#ifndef INCLUDE_TRACK_H_
#define INCLUDE_TRACK_H_

#include <cmath>
#include <memory>
#include <set>
#include <vector>

#include "DetectorModule.h"
#include "Hit.h"
#include <Math/Vector3D.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>

// Forward declaration
class Track;
class RILength;

#undef TRACK_DEBUG_RZ

//Typedef
typedef std::unique_ptr<Track> TrackPtr;
typedef std::vector<TrackPtr>  TrackCollection;

/**
 * @class Track
 * @details The Track class is essentially a collection of consecutive hits used to calculate track parameters
 * with correlations among them. These hits can bet set as active (used in tracking) or inactive (unused in
 * tracking). The track model assumes mag. field to be variable only in z: B = B(z).e_z + 0.e_x + 0.e_y, further
 * the parabolic approximation for track is considered (a la Karimaki + Multiple Coulomb scattering is taken into
 * account), disentangling the problem to fit 5 track parameters simultaneously into a problem of fitting track
 * independenly in R-Phi and s-Z plane. The covariance matrix is then 3x3 + 2x2. Parametrization used is the
 * following, one fits: phi, 1/R, impact parameter in R-Phi plane D0 & cotg(theta), impact parameter in s-Z
 * plane Z0. To provide tracking independently on tracker geometry, each track is assigned a tag, all sub-trackers
 * with the same tag are used then to find resolution of such a geometry.
 */
class Track {

public:

  //! Track constructor -> need to use setter methods to set: 2 of these [theta, phi, eta, cot(theta)] & 2 of these [mag. field, transv. momentum, radius]
  Track();

  //! Track copy-constructor -> creates deep copy of hit vector
  Track(const Track& track);

  //! Assign operator with deep copy of hit vector
  Track& operator=(const Track &track);

  //! Destructor
  ~Track();

  //! Main method calculating track parameters using Karimaki approach & parabolic approximation in R-Phi plane: 1/R, D0, phi parameters
  //! and using linear fit in s-Z plane: cotg(theta), Z0 parameters
  void computeErrors();

  //! Add new hit to track and return a pointer to that hit
  void addHit(HitPtr newHit);

  //! Add IP constraint to the track, technically new hit is assigned: with no material and hit resolution in R-Phi as dr, in s-Z as dz
  void addIPConstraint(double dr, double dz);

  //! Sort internally all hits assigned to this track -> sorting algorithm based on hit radius - by smaller radius sooner or vice-versa (inner-2-outer approach or vice-versa)
  void sortHits(bool bySmallerR);

  //! Does track contain no hits?
  bool hasNoHits() const {return m_hits.empty(); }

  //! Remove hits that don't follow the parabolic approximation used in tracking - TODO: still needs to be updated (not all approximations taken into account)
  bool pruneHits();
  
  //! Set active only hits with the given tag
  void keepTaggedOnly(const string& tag);

  //! Remove material from all assigned hits -> modify all hits such as they are without any material
  void removeMaterial();

  //! Helper method printing track covariance matrices in R-Phi
  void printErrors();

  //! Helper method printing track hits
  void printHits();

  //! Set track polar angle - theta, azimuthal angle - phi, particle transverse momentum - pt
  const Polar3DVector& setThetaPhiPt(const double& newTheta, const double& newPhi, const double& newPt);

  //! Set track origin
  const XYZVector& setOrigin(const double& X, const double& Y, const double& Z) {m_origin.SetCoordinates(X,Y,Z); return m_origin;}

  //! Re-set transverse momentum, pT
  void resetPt(double newPt) {m_pt = newPt;}

  // Getter methods
  double getTheta() const              { return m_theta;}
  double getEta() const                { return m_eta; }
  double getCotgTheta() const          { return m_cotgTheta; }
  double getPhi() const                { return m_phi;}
  double getTransverseMomentum() const { return m_pt; }
  double getPt() const                 { return m_pt; }
  
  // Calculate magnetic field at given z, assuming B = B(z).e_z + 0.e_x + 0 e_y
  double getMagField(double z) const;

  // Calculate radius or 1/R at given z, assuming B = B(z).e_z + 0.e_x + 0 e_y
  double getRho(double z) const        { return (getRadius(z)>0 ? 1/getRadius(z) : 0);}
  double getRadius(double z) const     { return m_pt / (0.3 * getMagField(z)); }

  const Polar3DVector& getDirection() const { return m_direction; }
  const XYZVector& getOrigin() const        { return m_origin; }

  //! Get number of active hits assigned to track for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file). If tag specified as "all" no extra tag required
  int getNActiveHits(std::string tag, bool useIP = true) const;

  //! Get the probabilty of having "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
  //! If tag specified as "all" no extra tag required
  std::vector<double> getHadronActiveHitsProbability(std::string tag);

  //! Get the probabilty of having a given number of "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
  double getHadronActiveHitsProbability(std::string tag, int nHits);

  //! Get number of hits assigned
  int getNHits() const { return m_hits.size(); }

  //! Get track material
  RILength getMaterial();

  //! Get all tags assigned to the track (see m_tags variable below
  const std::set<std::string>& getTags() const { return m_tags; }

  //! Get a vector of pairs: Detector module & hit type for Trigger hits
  std::vector<std::pair<const DetectorModule*, HitType>> getHitModules() const;

  const double& getDeltaRho() const      { return m_deltaRho; }
  const double& getDeltaPhi() const      { return m_deltaPhi; }
  const double& getDeltaD0() const       { return m_deltaD0; }
  const double& getDeltaCtgTheta() const { return m_deltaCtgTheta; }
  const double& getDeltaZ0() const       { return m_deltaZ0; }
  const double& getDeltaPtOverPt() const { return m_deltaPt; }
  const double& getDeltaPOverP()  const  { return m_deltaP; }

//  void addEfficiency(double efficiency, bool alsoPixel = false);
//  void keepTriggerOnly();
//  void setTriggerResolution(bool isTrigger);
//  // static bool debugRemoval; // debug
//  double expectedTriggerPoints(const double& triggerMomentum) const;
//#ifdef TRACK_DEBUG_RZ
//  static bool debugRZCovarianceMatrix;  // debug
//  static bool debugRZCorrelationMatrix;  // debug
//  static bool debugRZErrorPropagation;  // debug
//#endif

protected:

  //! Compute the correlation matrix of the track parameters in R-Z projection
  void computeCorrelationMatrixRZ();
  //! Compute the covariance matrix of the track parameters in R-Z projection
  void computeCovarianceMatrixRZ();

  //! Compute the correlation matrix of the track parameters in R-Phi projection
  void computeCorrelationMatrixRPhi();
  //! Compute the covariance matrix of the track parameters in R-Phi projection
  void computeCovarianceMatrixRPhi();

  double m_theta;              //!< Track shot at given theta & phi, i.e. theta at primary vertex
  double m_phi;                //!< Track shot at given theta & phi, i.e. phi at primary vertex
  double m_cotgTheta;          //!< Automatically calculated from theta at [0,0]
  double m_eta;                //!< Automatically calculated from eta at [0,0]
  double m_pt;                 //!< Particle transverse momentum (assuming B = fce of z only -> pT doesn't change along the path, only radius changes)

  Polar3DVector  m_direction;  //!< Track parameters as a 3-vector: R, theta, phi
  XYZVector      m_origin;     //!< Track origin as a 3-vector: X, Y, Z

  HitCollection         m_hits;//!< Track assigned hits
  std::set<std::string> m_tags;//!< Which subdetectors to be used in the track (each subdetector is tagged by a set of tags, e.g. pixel, fwd, tracker -> used in tracking of pixels, fwd tracking & full tracker)

  // Track resolution as a function of momentum
  TMatrixTSym<double>  m_correlationsRPhi; //!< Correlation matrix in R-Phi: pT (1/R), phi, D0
  TMatrixT<double>     m_covariancesRPhi;  //!< Covariance matrix in R-Phi: pT (1/R), phi, D0
  TMatrixTSym<double>  m_correlationsRZ;   //!< Correlation matrix in R-Phi: Z0, cotg(theta)
  TMatrixT<double>     m_covariancesRZ;    //!< Covariance matrix in R-Phi: Z0, cotg(theta)

  double m_deltaRho;       //!< Tracking error on 1/R
  double m_deltaPhi;       //!< Tracking error on Phi angle
  double m_deltaD0;        //!< Tracking error on impact parameter in R-Phi plane: d0
  double m_deltaCtgTheta;  //!< Tracking error on cotg(theta)
  double m_deltaZ0;        //!< Tracking error on impact parameter in s-z plane: Z0
  double m_deltaPt;        //!< Tracking error on transverse momentum
  double m_deltaP;         //!< Tracking error on total momentum

}; // Class

#endif /* INCLUDE_TRACK_H_ */

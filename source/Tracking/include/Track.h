/**
 * Track.h
 *
 * Created on: 20. 4. 2016
 *     Author: Drasal (CERN)
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
 * together with their correlation matrix. These hits can bet set as active (used in tracking) or inactive (unused in
 * tracking). The track model assumes mag. field to be variable only in z: B = B(z).e_z + 0.e_x + 0.e_y, further
 * the parabolic approximation for track is considered (a la Karimaki parametrization + Multiple Coulomb scattering
 * is taken into account), disentangling the problem to fit 5 track parameters simultaneously into a problem of
 * fitting track independenly in R-Phi and s-Z plane. The direction, in which Multiple scattering is calculated is
 * defined by particle direction or track propagation requirements. Either inside-out or outside-in can be set. The final
 * covariance matrix is 3x3 (rho, phi0, d0) + 2x2 (cotg(theta), z0). The parametrization used is the following, one fits:
 * phi0, 1/R, impact parameter in R-Phi plane d0 & cotg(theta), impact parameter in s-Z plane z0. To provide tracking independently
 * on complexity of defined tracker geometry, each track is assigned a tag, all sub-trackers being assigned the same tag
 * are then used in the track for tracking. In order to extract individual track parameters at arbitrary reference
 * point (not only at [0,0,0]), the standard technique of covariance propagator is applied (for details see getDelta***()
 * methods). The propagators are currently implemented for const Bz, not for case of B = B(z).
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

  //
  // Define track properties

  //! Add new hit to track and return a pointer to that hit (+ prune hits + resort hits)
  void addHit(HitPtr newHit);

  //! Add IP constraint (+ prune hits + resort hits) to the track, technically new hit is assigned: with no material and hit resolution in R-Phi as dr, in s-Z as dz
  void addIPConstraint(double dr, double dz);

  //! Set track polar angle - theta, azimuthal angle - phi, particle transverse momentum - pt (signed: + -> particle in-out, - -> particle out-in)
  const Polar3DVector& setThetaPhiPt(const double& newTheta, const double& newPhi, const double& newPt);

  //! Set track origin
  const XYZVector& setOrigin(const double& X, const double& Y, const double& Z) {m_origin.SetCoordinates(X,Y,Z); return m_origin;}

  //! Re-set transverse momentum + resort hits (if changing direction) + initiate recalc of cov matrices + prune hits (otherwise they may not lie on the new track, originally found at high pT limit)
  void resetPt(double newPt);

  //! Set active only hits with the given tag. If tag="all" all hits coming from measurement planes or IP will be set as active
  void keepTaggedHitsOnly(const string& tag, bool useIP = true);

  //! Keep only first N hits coming from measurement planes or IP as active -> return true if possible N>= size of hits vector
  bool keepFirstNHitsActive(signed int N, bool useIP = true);

  //! Keep only last N hits coming from measurement planes or IP as active -> return true if possible N>= size of hits vector
  bool keepLastNHitsActive(signed int N, bool useIP = true);

  //! Remove material from all assigned hits -> modify all hits such as they are without any material
  void removeMaterial();

  //
  // Print methods

  //! Helper method printing track covariance matrices in R-Phi
  void printErrors();

  //! Helper method printing symmetric matrix
  void printSymMatrix(const TMatrixTSym<double>&);

  //! Helper method printing general matrix
  void printMatrix(const TMatrixT<double>&);

  //! Helper method printing track hits
  void printHits();

  //! Helper method printing active track hits
  void printActiveHits();

  //! Does track contain no hits?
  bool hasNoHits() const {return m_hits.empty(); }

  //
  // Getter methods
  double getTheta() const              { return m_theta;}
  double getEta() const                { return m_eta; }
  double getCotgTheta() const          { return m_cotgTheta; }
  double getPhi() const                { return m_phi;}
  double getTransverseMomentum() const { return m_pt; }
  double getPt() const                 { return m_pt; }
  
  //! Get DeltaRho (error on 1/R) at path length s projected to XY plane, i.e. at [r,z] sXY ~ r
  //! Using 3x3 covariance propagator in case [r,z]!=[0,0]
  double getDeltaRho(double rPos);

  //! Get DeltaPtOvePt at path length s projected to XY plane, i.e. at [r,z] sXY ~ r (utilize the calculated deltaRho quantity)
  double getDeltaPtOverPt(double rPos);

  //! Get DeltaPOverP at path length s projected to XY plane, i.e. at [r,z] sXY ~ r (utilize deltaRho & deltaCotgTheta quantities)
  double getDeltaPOverP(double rPos);

  //! Get DeltaPhi0 at point (r,z) at path length s projected to XY plane, i.e. at [r,z] sXY ~ r
  //! Using 3x3 covariance propagator in case [r,z]!=[0,0]
  double getDeltaPhi(double rPos);
  double getDeltaPhi0() { return getDeltaPhi(0.0); };

  //! Get DeltaD0 at path length s projected to XY plane, i.e. at [r,z] sXY ~ r
  //! Using 3x3 covariance propagator in case [r,z]!=[0,0]
  double getDeltaD(double rPos);
  double getDeltaD0() { return getDeltaD(0.0); };

  //! Get DeltaCtgTheta at path length s projected to XY plane, i.e. at [r,z] sXY ~ r (independent on sXY)
  double getDeltaCtgTheta();

  //! Get DeltaZ0 at path length s projected to XY plane, i.e. at [r,z] sXY ~ r
  //! Using 2x2 covariance propagator in case [r,z]!=[0,0]
  double getDeltaZ(double rPos);
  double getDeltaZ0() { return getDeltaZ(0.0); }

  // Calculate magnetic field at given z, assuming B = B(z).e_z + 0.e_x + 0 e_y
  double getMagField(double zPos) const;

  // Calculate radius or 1/R at given z, assuming B = B(z).e_z + 0.e_x + 0 e_y
  double getRho(double zPos) const    { return (getRadius(zPos)!=0 ? 1/getRadius(zPos) : 0);}
  double getRadius(double zPos) const { return fabs(m_pt / (0.3 * getMagField(zPos))); }

  const Polar3DVector& getDirection() const { return m_direction; }
  const XYZVector&     getOrigin() const    { return m_origin; }

  //! Get number of active hits assigned to track for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file). If tag specified as "all" no extra tag required
  int getNActiveHits(std::string tag, bool useIP = true) const;

  //! Get number of hits, set as active & coming from measurement planes or IP constraint. All hits must be assigned to tracker with given tag. If tag specified as "all", all module & IP hits assigned.
  int getNMeasuredHits(std::string tag, bool useIP = true) const;

  //! Get reference to a hit, which can be measured, i.e. coming from measurement plane (active or inactive) or IP constraint
  const Hit* getMeasurableOrIPHit(int iHit) const;

  //! Reverse search -> Get reversely reference to a hit, which can be measured, i.e. coming from measurement plane (active or inactive) or IP constraint
  const Hit* getRMeasurableOrIPHit(int iHit) const;

  //! Get the probabilty of having "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
  //! If tag specified as "all" no extra tag required
  std::vector<double> getHadronActiveHitsProbability(std::string tag);

  //! Get the probabilty of having a given number of "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
  double getHadronActiveHitsProbability(std::string tag, int nHits);

  //! Get number of hits assigned
  int getNHits() const { return m_hits.size(); }

  //! Hits collection iterators
  HitCollection::const_iterator getBeginHits() { return m_hits.cbegin(); }
  HitCollection::const_iterator getEndHits()   { return m_hits.cend(); }
  //HitCollection::const_iterator getRBeginHits(){ return m_hits.crbegin(); }
  //HitCollection::const_iterator getREndHits()  { return m_hits.crend(); }

  //! Get track material
  RILength getMaterial() const;

  //! Get all tags assigned to the track (see m_tags variable below
  const std::set<std::string>& getTags() const { return m_tags; }

  //! Get a vector of pairs: Detector module & hit type for Trigger hits
  std::vector<std::pair<const DetectorModule*, HitType>> getHitModules() const;

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

  //! Main method calculating track parameters in s-z plane only, using linear fit: cotg(theta), z0 parameters -> internally calling computation of varMatrixRZ & covMatrixRZ
  void computeErrorsRZ();
  //! Compute the variance matrix in R-Z (s-Z): NxN (N hits = K+L: K active hits on detectors + L passive (artificial) hits due to material)
  //! Return true if V invertable
  bool computeVarianceMatrixRZ();
  //! Compute 2x2 covariance matrix of the track parameters in R-Z (s-Z) projection
  void computeCovarianceMatrixRZ();

  //! Main method calculating track parameters in r-phi plane only, using Karimaki parametrization & parabolic track
  //! approximation in R-Phi plane: 1/R, d0, phi parameters -> call this method before calling getDelta***() methods.
  //! -> internally calling computation of varMatrixRPhi & covMatrixRPhi
  void computeErrorsRPhi();
  //! Compute the variance matrix in R-Phi: NxN (N hits = K+L: K active hits on detectors + L passive (artificial) hits due to material)
  //! Return true if V invertable
  bool computeVarianceMatrixRPhi();
  //! Compute 3x3 covariance matrix of the track parameters in R-Phi projection
  void computeCovarianceMatrixRPhi();

  //! Helper fce returning derivative: d[f(rho, d0, phi0)]/d[rho], where f approximates
  //! a helix by set of parabolas. In general, N connected parabolas used, for const B
  //! field only one parabola assumed.
  double computeDfOverDRho(double rPos, double zPos);

  //! Sort internally all hits assigned to this track -> sorting algorithm based on hit radius - by smaller radius sooner or vice-versa (inner-2-outer approach or vice-versa)
  //! Sorting naturally affects whether in->out approach is being used or out->in is being used. Approach set directly by dedicated setter methods.
  void sortHits(bool bySmallerR);

  //! Remove hits that don't follow the parabolic approximation used in tracking - TODO: still needs to be updated (not all approximations taken into account)
  bool followsParabolicApprox(double rPos, double zPos) { return rPos<2*getRadius(zPos); }
  bool pruneHits();

  double m_theta;               //!< Track shot at given theta & phi, i.e. theta at primary vertex
  double m_phi;                 //!< Track shot at given theta & phi, i.e. phi at primary vertex
  double m_pt;                  //!< Particle transverse momentum (assuming B = fce of z only -> pT doesn't change along the path, only radius changes), pT sign: + -> particle traverses inside-out, - -> particle traverses outside-in
  double m_cotgTheta;           //!< Automatically calculated from theta at [0,0]
  double m_eta;                 //!< Automatically calculated from eta at [0,0]

  Polar3DVector  m_direction;   //!< Track parameters as a 3-vector: R, theta, phi
  XYZVector      m_origin;      //!< Track origin as a 3-vector: X, Y, Z TODO: For tracking model origin assumed to be at [0,0,0]

  bool           m_reSortHits;  //!< Caching whether necessary to resort hits (sorting will be done again if a new hit added or direction changed)

  bool           m_covRPhiDone; //!< Caching whether errors in R-Phi already calculated (will be recalculated, if direction of propagation changed, or added new hit etc.)
  bool           m_covRZDone;   //!< Caching whether errors in R-Z already calculated (will be recalculated, if direction of propagation changed, or added new hit etc.)
//  double         m_ZPointPropag;//!< Caching Z coordinate, with respect to which the errors calculation have been done last time
//  double         m_DPointPropag;//!< Caching R-Phi coordinate (D), with respect to which the errors calculation have been done last time
//  bool           m_InOutPropag; //!< Caching whether MS effects propagated in-out (true) or out-in (false),

  HitCollection         m_hits;         //!< Hits assigned to track
  std::set<std::string> m_tags;         //!< Which subdetectors to be used in tracking (each subdetector is tagged by a set of tags, e.g. pixel, fwd, tracker -> used in tracking of pixels, fwd tracking & full tracker)

  // Track parameters covariance matrices
  TMatrixTSym<double>   m_varMatrixRPhi; //!< NxN (hits) Variance matrix in R-Phi: V(NxN) (N hits = K+L: K active hits on detectors + L passive (artificial) hits due to material)
  TMatrixT<double>      m_covMatrixRPhi; //!< 3x3 Covariance matrix in R-Phi: pT (1/R), phi, d0 at ref. point [r,z] = [0,0]: C(3x3)
  TMatrixTSym<double>   m_varMatrixRZ;   //!< NxN (hits) Variance matrix in R-Z V(NxN) (N hits = K+L: K active hits on detectors + L passive (artificial) hits due to material)
  TMatrixT<double>      m_covMatrixRZ;   //!< 2x2 Covariance matrix in R-Phi: z0, cotg(theta) at ref. point [r,z] = [0,0]: C(2x2)

}; // Class

#endif /* INCLUDE_TRACK_H_ */

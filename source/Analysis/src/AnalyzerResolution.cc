/*
 * AnalyzerResolution.cc
 *
 *  Created on: 20. 4. 2016
 *      Author: drasal
 */
#include <AnalyzerResolution.h>

#include "BeamPipe.h"
#include "DetectorModule.h"
#include "MainConfigHandler.h"
#include "MaterialProperties.h"
#include "MessageLogger.h"
#include "ModuleCap.h"
#include "rootweb.h"
#include "SimParms.h"
#include <TRandom3.h>
#include "Tracker.h"
#include "Track.h"
#include "Units.h"


AnalyzerResolution::AnalyzerResolution(std::vector<const Tracker*> trackers, const BeamPipe* beamPipe) : AnalyzerUnit("AnalyzerResolution", trackers, beamPipe),
 m_nTracks(0),
 m_etaMin(-1*geom_max_eta_coverage),
 m_etaMax(+1*geom_max_eta_coverage)
{};

bool AnalyzerResolution::init(int nTracks)
{
  // Set nTracks
  m_nTracks = nTracks;

  if (m_nTracks<=0 || m_trackers.size()==0) {

    std::ostringstream message;
    message << "AnalyzerResolution::init(): Number of simulation tracks zero or negative: " << m_nTracks << " or zero number of trackers to be initialized";
    logERROR(message.str());
    return false;
  }
  else {

    m_isInitOK = true;
    return m_isInitOK;
  }
}

bool AnalyzerResolution::analyze()
{
  // Check that initialization OK
  if (!m_isInitOK) return false;

  // Random generator
  TRandom3 myDice;
  myDice.SetSeed(random_seed);

  // Initialize
  double efficiency  = SimParms::getInstance()->efficiency();

  // Tracks pruned
  bool isPruned = false;

  for (int iTrack = 0; iTrack < m_nTracks; iTrack++) {

    // Define track
    Track matTrack;

    double eta   = 0 + m_etaMax/m_nTracks*(iTrack+0.5);
    double theta = 2 * atan(exp(-eta));
    double phi   = myDice.Rndm() * M_PI * 2.0;
    double pT    = 100*Units::TeV; // Arbitrarily high number

    matTrack.setThetaPhiPt(theta, phi, pT);
    matTrack.setOrigin(0, 0, 0); // TODO: Not assuming z-error when calculating Material budget

    // Assign material to the track
    MatTrackVisitor matVisitor(matTrack);
    for (auto iTracker : m_trackers) iTracker->accept(matVisitor);  // Assign to material track hits corresponding to modules
    m_beamPipe->accept(matVisitor);                                 // Assign to material track hit corresponding to beam-pipe

    // Go through tracks with non-zero number of hits
    if (!matTrack.hasNoHits()) {

      for (string tag : matTrack.getTags()) {

        // Analyze only hits within track coming from the defined bunch of detectors (the same tag)
        matTrack.keepTaggedOnly(tag);

        // Add IP constraint
        if (SimParms::getInstance()->useIPConstraint()) matTrack.addIPConstraint(SimParms::getInstance()->rError(), SimParms::getInstance()->zErrorCollider());

        // Sort hits
        bool smallerRadius = true;
        matTrack.sortHits(true);
        //matTrack.printHits();

        // Remove some hits randomly based on inefficiency parameter
        //if (efficiency!=1) matTrack.addEfficiency(pxlEfficiency, stripEfficiency);

        // For each momentum/transverse momentum compute the tracks error
        for (const auto& pIter : MainConfigHandler::getInstance().getMomenta()) {

          int    parameter = pIter/Units::MeV; // Store p or pT in MeV as int (key to the map)
          double momentum  = pIter;

          // Case I) Initial momentum is equal to pT
          double pT = momentum;

          // Active+passive material
          TrackPtr trackPt(new Track(matTrack));
          trackPt->setThetaPhiPt(theta, phi, pT);

          // Remove tracks with less than 3 hits
          bool pruned = trackPt->pruneHits();
          if (pruned) isPruned = true;

          if (trackPt->getNActiveHits(tag, true)>2) {

            trackPt->computeErrors();
            std::map<int, TrackCollection>& myMap        = taggedTrackPtCollectionMap[tag];
            TrackCollection&                myCollection = myMap[parameter];
            myCollection.push_back(std::move(trackPt));
          }

          // Ideal (no material)
          TrackPtr idealTrackPt(new Track(matTrack));
          idealTrackPt->setThetaPhiPt(theta, phi, pT);

          // Remove tracks with less than 3 hits
          pruned = idealTrackPt->pruneHits();
          if (pruned) isPruned = true;

          // Remove material
          idealTrackPt->removeMaterial();
          if (idealTrackPt->getNActiveHits(tag, true)>2) {

            idealTrackPt->computeErrors();
            std::map<int, TrackCollection>& myMapIdeal        = taggedTrackPtCollectionMapIdeal[tag];
            TrackCollection&                myCollectionIdeal = myMapIdeal[parameter];
            myCollectionIdeal.push_back(std::move(idealTrackPt));
          }

          // Case II) Initial momentum is equal to p
          pT = momentum*sin(theta);

          // Active+passive material
          TrackPtr trackP(new Track(matTrack));
          trackP->setThetaPhiPt(theta, phi, pT);

          // Remove tracks with less than 3 hits
          pruned = trackP->pruneHits();
          if (pruned) isPruned = true;

          if (trackP->getNActiveHits(tag, true)>2) {

            trackP->computeErrors();
            std::map<int, TrackCollection>& myMapII        = taggedTrackPCollectionMap[tag];
            TrackCollection&                myCollectionII = myMapII[parameter];
            myCollectionII.push_back(std::move(trackP));
          }

          // Ideal (no material)
          TrackPtr idealTrackP(new Track(matTrack));
          idealTrackP->setThetaPhiPt(theta, phi, pT);

          // Remove tracks with less than 3 hits
          pruned = idealTrackP->pruneHits();
          if (pruned) isPruned = true;

          // Remove material
          idealTrackP->removeMaterial();

          if (idealTrackP->getNActiveHits(tag, true)>2) {

            idealTrackP->computeErrors();
            std::map<int, TrackCollection>& myMapIdealII        = taggedTrackPCollectionMapIdeal[tag];
            TrackCollection&                myCollectionIdealII = myMapIdealII[parameter];
            myCollectionIdealII.push_back(std::move(idealTrackP));
          }
        } // For momenta
      } // For tags
    }
  } // For tracks

  // Log pruning procedure -> no result mode might occur if starting momenta wrongly set
  if (isPruned) {

    std::string message = std::string("Resolution - some tracks pruned! Hits that don't follow the parabolic approximation removed. Check momenta if no result appears!");
    logWARNING(message);
  }

  return true;
}

bool AnalyzerResolution::visualize(RootWSite& webSite)
{
  // Check that initialization & analysis OK
  if (!m_isInitOK && !m_isAnalysisOK) return false;

  return true;
}

//
// Material budget visitor - constructor
//
MatTrackVisitor::MatTrackVisitor(Track& matTrack) :
    m_matTrack(matTrack)
{
}

//
// Destructor
//
MatTrackVisitor::~MatTrackVisitor()
{
}

//
// Visit BeamPipe -> update track with beam pipe hit
//
void MatTrackVisitor::visit(const BeamPipe& bp)
{
  // Add hit corresponding with beam-pipe
  double theta    = m_matTrack.getTheta();
  double distance = (bp.radius()+bp.thickness())/2./sin(theta);
  HitPtr hit(new Hit(distance));
  hit->setOrientation(HitOrientation::Horizontal);
  hit->setObjectKind(HitKind::Inactive);

  Material material;
  material.radiation   = bp.radLength()/sin(theta);
  material.interaction = bp.intLength()/sin(theta);
  hit->setCorrectedMaterial(material);
  hit->setBeamPipe(true);
  m_matTrack.addHit(std::move(hit));
}

//
// Visit BarrelModule (no limits on Rods, Layers or Barrels)
//
void MatTrackVisitor::visit(const BarrelModule& m)
{
  analyzeModuleMB(m);
}

//
// Visit EndcapModule (no limits on Rings or Endcaps)
//
void MatTrackVisitor::visit(const EndcapModule& m)
{
  analyzeModuleMB(m);
}

//
// Analyze if module crossed by given track & how much material is in the way
//
void MatTrackVisitor::analyzeModuleMB(const DetectorModule& m)
{
  // Collision detection: material tracks being shot in z+ only, so consider only modules that lie on +Z side
  if (m.maxZ() > 0) {

    XYZVector direction(m_matTrack.getDirection());

    auto pair    = m.checkTrackHits(m_matTrack.getOrigin(), direction);
    auto hitRho  = pair.first.rho();
    auto hitType = pair.second;

    if (hitType!=HitType::NONE) {

      Material material;
      material.radiation   = m.getModuleCap().getRadiationLength();
      material.interaction = m.getModuleCap().getInteractionLength();

      // Fill material map
      double theta = m_matTrack.getTheta();
      double rho   = hitRho;
      double z     = rho/tan(theta);

      // Treat barrel & endcap modules separately
      double tiltAngle = m.tiltAngle();

      if (m.subdet() == BARREL) {

        material.radiation   /= sin(theta + tiltAngle);
        material.interaction /= sin(theta + tiltAngle);

      }
      else if (m.subdet() == ENDCAP) {

        material.radiation   /= cos(theta + tiltAngle - M_PI/2); // Endcap has tiltAngle = pi/2
        material.interaction /= cos(theta + tiltAngle - M_PI/2); // Endcap has tiltAngle = pi/2

      }
      else {
        logWARNING("MatTrackVisitor::analyzeModuleMB -> incorrectly scaled material, unknown module type. Neither barrel or endcap");
      }

      // Create Hit object with appropriate parameters, add to Track t
      HitPtr hit(new Hit(pair.first.R(), &m, hitType));
      hit->setCorrectedMaterial(material);
      m_matTrack.addHit(std::move(hit));
    }
  } // Z>0
}


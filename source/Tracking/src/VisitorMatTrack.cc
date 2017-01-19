/*
 * VisitorMatTrack.cc
 *
 *  Created on: 17. 1. 2017
 *      Author: drasal
 */
#include "../include/VisitorMatTrack.h"

#include "Barrel.h"
#include "BeamPipe.h"
#include "DetectorModule.h"
#include "InactiveElement.h"
#include "MaterialProperties.h"
#include "ModuleCap.h"
#include "SupportStructure.h"
#include "Track.h"

//
// Material track visitor - constructor
//
VisitorMatTrack::VisitorMatTrack(Track& matTrack) :
    m_matTrack(matTrack)
{
}

//
// Destructor
//
VisitorMatTrack::~VisitorMatTrack()
{
}

//
// Visit BeamPipe -> update track with beam pipe hit
//
void VisitorMatTrack::visit(const BeamPipe& bp)
{
  // Add hit corresponding with beam-pipe
  double theta = m_matTrack.getTheta();
  double eta   = m_matTrack.getEta();
  double rPos  = (bp.radius()+bp.thickness()/2.);
  double zPos  = rPos/tan(theta);

  HitPtr hit(new Hit(rPos, zPos));
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
// Visit Barrel
//
void VisitorMatTrack::visit(const Barrel& b)
{
  for (auto& element : b.services()) {
    analyzeInactiveElement(element);
  }
}

//
// Visit BarrelModule (no limits on Rods, Layers or Barrels)
//
void VisitorMatTrack::visit(const BarrelModule& m)
{
  analyzeModuleMB(m);
}

//
// Visit EndcapModule (no limits on Rings or Endcaps)
//
void VisitorMatTrack::visit(const EndcapModule& m)
{
  analyzeModuleMB(m);
}

//
// Visit Support strucutre
//
void VisitorMatTrack::visit(const SupportStructure& s)
{
  for (auto& elem : s.inactiveElements()) {
    analyzeInactiveElement(elem);
  }
}

//
// Analyze if module crossed by given track & how much material is in the way
//
void VisitorMatTrack::analyzeModuleMB(const DetectorModule& m)
{
  // Collision detection: material tracks being shot in z+ only, so consider only modules that lie on +Z side
  if (m.maxZ() > 0) {

    XYZVector direction(m_matTrack.getDirection());

    auto pair    = m.checkTrackHits(m_matTrack.getOrigin(), direction);
    //auto hitDistance = pair.first.R();
    auto hitRPos = pair.first.rho();
    auto hitZPos = pair.first.z();
    auto hitType = pair.second;

    if (hitType!=HitType::NONE) {

      Material material;
      material.radiation   = m.getModuleCap().getRadiationLength();
      material.interaction = m.getModuleCap().getInteractionLength();

      // Fill material map
      double theta = m_matTrack.getTheta();

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
        logWARNING("VisitorMatTrack::analyzeModuleMB -> incorrectly scaled material, unknown module type. Neither barrel or endcap");
      }

      // Create Hit object with appropriate parameters, add to Track t
      HitPtr hit(new Hit(hitRPos, hitZPos, &m, hitType));
      hit->setCorrectedMaterial(material);
      m_matTrack.addHit(std::move(hit));
    }
  } // Z>0
}

//
// Helper method - analyse inactive element & estimate how much material is in the way
//
void VisitorMatTrack::analyzeInactiveElement(const insur::InactiveElement& e)
{
  // Work-out only if inactive element has non-zero material assigned, otherwise no effect of inactive material
  if (e.getRadiationLength()!=0 || e.getInteractionLength()!=0) {

    // Collision detection: rays are shot in z+ only, so only volumes in z+ need to be considered
    // only volumes of the requested category, or those without one (which should not exist) are examined
    if ((e.getZOffset() + e.getZLength()) > 0) {

      // collision detection: check eta range
      auto etaMinMax = e.getEtaMinMax();

      // Volume hit
      double eta   = m_matTrack.getEta();
      double theta = m_matTrack.getTheta();

      if ((etaMinMax.first < eta) && (etaMinMax.second > eta)) {

        /*
        if (eta<0.01) {
          std::cout << "Hitting an inactive surface at z=("
                    << iter->getZOffset() << " to " << iter->getZOffset()+iter->getZLength()
                    << ") r=(" << iter->getInnerRadius() << " to " << iter->getInnerRadius()+iter->getRWidth() << ")" << std::endl;
          const std::map<std::string, double>& localMasses = iter->getLocalMasses();
          const std::map<std::string, double>& exitingMasses = iter->getExitingMasses();
          for (auto massIt : localMasses) std::cerr   << "       localMass" <<  massIt.first << " = " << any2str(massIt.second) << " g" << std::endl;
          for (auto massIt : exitingMasses) std::cerr << "     exitingMass" <<  massIt.first << " = " << any2str(massIt.second) << " g" << std::endl;
        }
        */

        // Initialize
        Material material;
        material.radiation   = 0.0;
        material.interaction = 0.0;

        double rPos = 0.0;
        double zPos = 0.0;

        // Radiation and interaction lenth scaling for vertical volumes
        if (e.isVertical()) {

          zPos = e.getZOffset() + e.getZLength() / 2.0;
          rPos = zPos * tan(theta);

          material.radiation   = e.getRadiationLength();
          material.interaction = e.getInteractionLength();

          // Special treatment for user-defined supports as they can be very close to z=0
          if (e.getCategory() == MaterialProperties::u_sup) {

            double s = e.getZLength() / cos(theta);
            if (s > (e.getRWidth() / sin(theta))) s = e.getRWidth() / sin(theta);

            // add the hit if it's declared as inside the tracking volume, add it to 'others' if not
            //if (e.track()) {}

            material.radiation  *= s / e.getZLength();
            material.interaction*= s / e.getZLength();
          }
          else {

            material.radiation   /= cos(theta);
            material.interaction /= cos(theta);
          }
        }
        // Radiation and interaction length scaling for horizontal volumes
        else {

          rPos = e.getInnerRadius() + e.getRWidth() / 2.0;
          zPos = rPos/tan(theta);

          material.radiation   = e.getRadiationLength();
          material.interaction = e.getInteractionLength();

          // Special treatment for user-defined supports; should not be necessary for now
          // as all user-defined supports are vertical, but just in case...
          if (e.getCategory() == MaterialProperties::u_sup) {

            double s = e.getZLength() / sin(theta);

            if (s > (e.getRWidth()/cos(theta))) s = e.getRWidth() / cos(theta);

            // add the hit if it's declared as inside the tracking volume, add it to 'others' if not
            //if (e.track()) {}

            material.radiation  *= s / e.getZLength();
            material.interaction*= s / e.getZLength();
          }
          else {

            material.radiation   /= sin(theta);
            material.interaction /= sin(theta);
          }
        }

        // Create Hit object with appropriate parameters, add to Track t
        HitPtr hit(new Hit(rPos, zPos));
        if (e.isVertical()) hit->setOrientation(HitOrientation::Vertical);
        else                hit->setOrientation(HitOrientation::Horizontal);
        hit->setObjectKind(HitKind::Inactive);
        hit->setCorrectedMaterial(material);
        m_matTrack.addHit(std::move(hit));

      } // Eta min max
    } // +Z check
  } // Non-zero material assigned
}




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
#include "Disk.h"
#include "Endcap.h"
#include "InactiveElement.h"
#include "Layer.h"
#include "MaterialProperties.h"
#include "ModuleCap.h"
#include "SupportStructure.h"
#include "Track.h"
#include "Tracker.h"

//
// Material track visitor - constructor
//
VisitorMatTrack::VisitorMatTrack(Track& matTrack) :
    m_matTrack(matTrack)
{
  m_detName  = "Undefined";
  m_brlName  = "Undefined";
  m_ecapName = "Undefined";
  m_layerID  = -1;
  m_discID   = -1;
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

  HitPtr hit(new Hit(rPos, zPos, nullptr, HitPassiveType::BeamPipe));

  Material material;
  material.radiation   = bp.radLength()/sin(theta);
  material.interaction = bp.intLength()/sin(theta);
  hit->setCorrectedMaterial(material);
  m_matTrack.addHit(std::move(hit));
}

//
// Visit Tracker
//
void VisitorMatTrack::visit(const Tracker& t)
{
  m_detName = t.myid();
  m_layerID    = -1;
  m_discID     = -1;
}

//
// Visit Barrel
//
void VisitorMatTrack::visit(const Barrel& b)
{
  m_brlName  = b.myid();
  m_layerID  = -1;
  m_discID   = -1;

  for (auto& element : b.services()) {
    analyzeInactiveElement(element);
  }
}

//
// Visit Endcap
//
void VisitorMatTrack::visit(const Endcap& e)
{
  m_ecapName = e.myid();
  m_layerID  = -1;
  m_discID   = -1;
}

//
// Visit Layer
//
void VisitorMatTrack::visit(const Layer& l)
{
  m_layerID = l.myid();
  m_discID  = -1;
}

//
// Visit Disc
//
void VisitorMatTrack::visit(const Disk& d)
{
  m_layerID = -1;
  m_discID  = d.myid();
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
  // Collision detection: material tracks being shot in z+ only, so consider only modules that lie on +Z side after correction on primary vertex origin
  if ((m.maxZ()-m_matTrack.getOrigin().Z())>0) {

    XYZVector direction(m_matTrack.getDirection());
    Material  hitMaterial;
    XYZVector hitPosition;
    HitType   hitType;

    if (m.checkTrackHits(m_matTrack.getOrigin(), direction, hitMaterial, hitType, hitPosition)) {

      auto hitRPos = hitPosition.rho();
      auto hitZPos = hitPosition.z();

      // Create Hit object with appropriate parameters, add to Track t
      HitPtr hit(new Hit(hitRPos, hitZPos, &m, hitType));
      hit->setCorrectedMaterial(hitMaterial);
      if (m.subdet() == BARREL) {
        hit->setDetName(m_detName+"_"+m_brlName);
        hit->setLayerID(m_layerID);
      }
      else {
        hit->setDetName(m_detName+"_"+m_ecapName);
        hit->setDiscID(m_discID);
      }
      m_matTrack.addHit(std::move(hit));
    }
  } // (module_Z-orig_Z)>0
}

//
// Helper method - analyse inactive element & estimate how much material is in the way
//
void VisitorMatTrack::analyzeInactiveElement(const InactiveElement& e)
{

  XYZVector direction(m_matTrack.getDirection());
  Material  hitMaterial;
  XYZVector hitPosition;

  // Hit found -> create
  if (e.checkTrackHits(m_matTrack.getOrigin(),direction, hitMaterial, hitPosition)) {

    auto hitRPos = hitPosition.rho();
    auto hitZPos = hitPosition.z();

    // Create Hit object with appropriate parameters, add to Track t
    if (e.getCategory() == MaterialProperties::b_sup ||
        e.getCategory() == MaterialProperties::e_sup ||
        e.getCategory() == MaterialProperties::u_sup ||
        e.getCategory() == MaterialProperties::t_sup ) {

      HitPtr hit(new Hit(hitRPos, hitZPos, &e, HitPassiveType::Support));
      hit->setCorrectedMaterial(hitMaterial);
      m_matTrack.addHit(std::move(hit));
    }
    else if (e.getCategory() == MaterialProperties::b_ser ||
             e.getCategory() == MaterialProperties::e_ser ) {

      HitPtr hit(new Hit(hitRPos, hitZPos, &e, HitPassiveType::Service));
      hit->setCorrectedMaterial(hitMaterial);
      m_matTrack.addHit(std::move(hit));
    }
  } // Hit found
}




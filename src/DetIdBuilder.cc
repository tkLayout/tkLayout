#include <cmath>
#include <cstdint>
#include <map>

#include "DetIdBuilder.hh"
#include "Layer.hh"
#include "RodPair.hh"
#include "DetectorModule.hh"
#include "Sensor.hh"
#include "Endcap.hh"
#include "Disk.hh"
#include "Ring.hh"
#include "global_funcs.hh"

namespace {
  constexpr size_t DETECTOR_LEVEL = 0;
  constexpr size_t SUBDETECTOR_LEVEL = 1;

  // Barrel geometry hierarchy

  namespace Barrels { // Pixel and Outer common levels
    constexpr size_t UNUSED_LEVEL = 2;
    constexpr size_t LAYER_LEVEL = 3;
  }

  namespace PixelBarrel {
    constexpr size_t LADDER_LEVEL = 4;
    constexpr size_t MODULE_LEVEL = 5;
    constexpr size_t SENSOR_LEVEL = 6;
  }

  namespace TOB {
    constexpr size_t RING_LEVEL = 4;
    constexpr size_t LADDER_LEVEL = 5;
    constexpr size_t MODULE_LEVEL = 6;
    constexpr size_t SENSOR_LEVEL = 7;
  }

  // Endcap geometry hierarchy

  namespace Endcaps { // Pixel and Outer common levels
    constexpr size_t SIDE_LEVEL = 2;
    constexpr size_t RING_LEVEL = 5;
    constexpr size_t PANEL_LEVEL = 6;
    constexpr size_t SENSOR_LEVEL = 8;
  }

  namespace PixelEndcap {
    constexpr size_t DISK_LEVEL = 3;
    constexpr size_t SUBDISK_LEVEL = 4;
    constexpr size_t MODULE_LEVEL = 7;
  }

  namespace TID {
    constexpr size_t UNUSED_LEVEL = 3;
    constexpr size_t DISK_LEVEL = 4;
    constexpr size_t MODULE_LEVEL = 7;
  }

  // CMSSW DetId scheme bit values for the Tracker and subdetectors DetIds
  namespace DetId {
    constexpr uint32_t Fixed_0 = 0;
    constexpr uint32_t Fixed_1 = 1;
    constexpr uint32_t Tracker = 1;
    constexpr uint32_t PixelBarrel = 1;
    constexpr uint32_t PixelEndcap = 2;
    constexpr uint32_t OuterEndcap = 4;
    constexpr uint32_t OuterTracker = 5;

    // Outer Tracker ring level bits based on tiltedness
    namespace OTRing {
      constexpr uint32_t NegativeZ = 1;
      constexpr uint32_t PositiveZ = 2;
      constexpr uint32_t Flat = 3;
    }

    // Endcap disk level bits based on the Z-coordinate of the disk
    namespace EndcapDisk {
      constexpr uint32_t NegativeZ = 1;
      constexpr uint32_t PositiveZ = 2;
    }
  }
}

    //************************************//
    //*               Visitor             //
    //*          BarrelDetIdBuilder       //
    //*                                   //
    //************************************//

BarrelDetIdBuilder::BarrelDetIdBuilder(bool isPixelTracker, std::vector<int> geometryHierarchySizes) : 
  isPixelTracker_(isPixelTracker),
  geometryHierarchySizes_(geometryHierarchySizes)
{}

void BarrelDetIdBuilder::visit(Barrel& b) {
  geometryHierarchyIds_[DETECTOR_LEVEL] = DetId::Tracker;
  geometryHierarchyIds_[SUBDETECTOR_LEVEL] = isPixelTracker_ ? DetId::PixelBarrel : DetId::OuterTracker;
  geometryHierarchyIds_[Barrels::UNUSED_LEVEL] = DetId::Fixed_0;
}

void BarrelDetIdBuilder::visit(Layer& l) {
  // Increasing in radius
  geometryHierarchyIds_[Barrels::LAYER_LEVEL] = l.layerNumber();

  // Store information for visit(RodPair&) and visit(BarrelModule&)
  numRods_ = l.numRods();
  numFlatRings_ = l.buildNumModulesFlat();
  numRings_ = l.isTilted() ? l.buildNumModules() : numFlatRings_;
}

void BarrelDetIdBuilder::visit(RodPair& r) {
  double phiSegment = 2 * M_PI / numRods_;          // Phi interval between 2 consecutive rods.
  double startAngle = femod( r.Phi(), phiSegment);  // Calculates the smallest positive rod's Phi in the Layer.

  // Calculates Phi position identifier.
  phiRef_ = 1 + round(femod(r.Phi() - startAngle, 2 * M_PI) / phiSegment);  
  // First rod with Phi > 0 has phiRef 1. 
  // Next rod (with bigger Phi) has phiRef 2. 
  // phiRef increments with increasing Phi.

  isCentered_ = (r.startZMode() == RodPair::StartZMode::MODULECENTER);
}

void BarrelDetIdBuilder::visit(BarrelModule& m) {
  const bool isPositiveZ = m.side() > 0;

  // Z+ modules are copied to the Z- side, reversing their order in |Z|
  uint32_t flatRingRef = isPositiveZ ? m.ring() + numFlatRings_ : 1 + numFlatRings_ - m.ring();
  // Remove the offset
  flatRingRef = isCentered_ && isPositiveZ ? flatRingRef - 1 : flatRingRef;

  if (isPixelTracker_) { // Pixel Barrel
    if (m.isTilted()) {
      logWARNING("No CMSSW DetId scheme exists for Pixel Barrel with tilted modules.");
    }
    else {
      // Increasing in phi
      geometryHierarchyIds_[PixelBarrel::LADDER_LEVEL] = phiRef_;
      // Increasing in |z|
      geometryHierarchyIds_[PixelBarrel::MODULE_LEVEL] = flatRingRef;
    }
  }
  else { // Outer Tracker Barrel
    if (m.isTilted()) {
      uint32_t tiltedRingRef = isPositiveZ ? m.ring() - numFlatRings_ : 1 + numRings_ - m.ring();

      geometryHierarchyIds_[TOB::RING_LEVEL] = isPositiveZ ? DetId::OTRing::PositiveZ : DetId::OTRing::NegativeZ;
      // Increasing in |z(rings)|
      geometryHierarchyIds_[TOB::LADDER_LEVEL] = tiltedRingRef;
      // Increasing in phi(barrel)
      geometryHierarchyIds_[TOB::MODULE_LEVEL] = phiRef_;
    }
    else {
      geometryHierarchyIds_[TOB::RING_LEVEL] = DetId::OTRing::Flat;
      // Increasing in phi(barrel)
      geometryHierarchyIds_[TOB::LADDER_LEVEL] = phiRef_;
      // Increasing in |z(rings)|
      geometryHierarchyIds_[TOB::MODULE_LEVEL] = flatRingRef;
    }
  }

  // Sensor level is irrelevant at module level, assign a fixed value
  geometryHierarchyIds_[isPixelTracker_ ? PixelBarrel::SENSOR_LEVEL : TOB::SENSOR_LEVEL] = DetId::Fixed_0;

  m.buildDetId(geometryHierarchyIds_, geometryHierarchySizes_);
}

void BarrelDetIdBuilder::visit(Sensor& s) {
  if (s.subdet() != ModuleSubdetector::BARREL) return;

  geometryHierarchyIds_[isPixelTracker_ ? PixelBarrel::SENSOR_LEVEL : TOB::SENSOR_LEVEL] = s.innerOuter();

  s.buildDetId(geometryHierarchyIds_, geometryHierarchySizes_);
}



    //************************************//
    //*               Visitor             //
    //*          EndcapDetIdBuilder       //
    //*                                   //
    //************************************//

EndcapDetIdBuilder::EndcapDetIdBuilder(bool isPixelTracker, bool hasSubDisks, std::vector<int> geometryHierarchySizes) : 
  isPixelTracker_(isPixelTracker),
  hasSubDisks_(hasSubDisks),
  geometryHierarchySizes_(geometryHierarchySizes)
{}

void EndcapDetIdBuilder::visit(Endcap& e) {
  geometryHierarchyIds_[DETECTOR_LEVEL] = DetId::Tracker;
  geometryHierarchyIds_[SUBDETECTOR_LEVEL] = isPixelTracker_ ? DetId::PixelEndcap : DetId::OuterEndcap;
}

void EndcapDetIdBuilder::visit(Disk& d) {
  geometryHierarchyIds_[Endcaps::SIDE_LEVEL] = d.side() ? DetId::EndcapDisk::PositiveZ : DetId::EndcapDisk::NegativeZ;
  // Increasing in |z|
  geometryHierarchyIds_[isPixelTracker_ ? PixelEndcap::DISK_LEVEL : TID::DISK_LEVEL] = d.diskNumber();

  // Store information for visit(Ring&)
  numEmptyRings_ = d.numEmptyRings();
}
   
void EndcapDetIdBuilder::visit(Ring& r) {
  geometryHierarchyIds_[Endcaps::PANEL_LEVEL] = DetId::Fixed_1;
  // Increasing in radius
  geometryHierarchyIds_[Endcaps::RING_LEVEL] = r.myid() - numEmptyRings_;

  // Store information for visit(EndcapModule&)
  numModules_ = r.numModules();
}
 
void EndcapDetIdBuilder::visit(EndcapModule& m) {
  double phiSegment = 2 * M_PI / numModules_;                // Phi interval between 2 consecutive modules.
  double startAngle = femod( m.center().Phi(), phiSegment);  // Calculates the smallest positive Phi in the Ring.

  // Calculates 1-based Phi identifier (range: 1 to numModules_)
  uint32_t phiRef_ = 1 + round(femod(m.center().Phi() - startAngle, 2*M_PI) / phiSegment);

  if (isPixelTracker_) { // Pixel Endcap
    if (hasSubDisks_) {
      const bool isOddModule = (phiRef_ % 2 == 1);
      // Increasing in |z|
      geometryHierarchyIds_[PixelEndcap::SUBDISK_LEVEL] = isOddModule ? DetId::Fixed_0 : DetId::Fixed_1;
      // Increasing in phi
      geometryHierarchyIds_[PixelEndcap::MODULE_LEVEL] = isOddModule ? ((phiRef_ + 1) / 2) : (phiRef_ / 2);
    }
    else {
      logWARNING("No CMSSW DetId scheme exists for Pixel Endcap with Single Disk.");
    }
  }
  else { // Outer Tracker Endcap
    if (hasSubDisks_) {
      logWARNING("No CMSSW DetId scheme exists for Outer Endcap with SubDisks.");
    }
    else {
      geometryHierarchyIds_[TID::UNUSED_LEVEL] = DetId::Fixed_0;
      // Increasing in phi
      geometryHierarchyIds_[TID::MODULE_LEVEL] = phiRef_;
    }
  }

  // First module with Phi > 0 has phiRef 1. 
  // Next module (with bigger Phi) has phiRef 2. 
  // phiRef increments with increasing Phi.
  // If we have subdisks, apply the same logic but to 
  // modules in individual subdisks

  // Sensor level is irrelevant at module level, assign a fixed value
  geometryHierarchyIds_[Endcaps::SENSOR_LEVEL] = DetId::Fixed_0;

  m.buildDetId(geometryHierarchyIds_, geometryHierarchySizes_);
}

void EndcapDetIdBuilder::visit(Sensor& s) {  
  if (s.subdet() != ModuleSubdetector::ENDCAP) return;

  geometryHierarchyIds_[Endcaps::SENSOR_LEVEL] = isPixelTracker_ ? DetId::Fixed_0 : s.innerOuter();

  s.buildDetId(geometryHierarchyIds_, geometryHierarchySizes_);
}

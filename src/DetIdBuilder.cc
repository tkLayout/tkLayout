#include <DetIdBuilder.hh>


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

  // !!! Tracker level
  geometryHierarchyIds_[0] = 1;   // always 1.

  // !!! Barrel level
  if (!isPixelTracker_) { geometryHierarchyIds_[1] = 205 % 100; }  // Outer Tracker : always 5. 
                                                        // In CMSSW, Phase 2 Outer Tracker Barrel is identified by 205, and its id is 205 % 100.
  else { geometryHierarchyIds_[1] = 201 % 100; }                   // Inner Tracker : always 1.
                                                        // In CMSSW, Phase 2 Inner Tracker Barrel is identified by 201, and its id is 201 % 100.

  // !!! Unused hierarchy level			   
  geometryHierarchyIds_[2] = 0;  // always 0.
}

void BarrelDetIdBuilder::visit(Layer& l) {
  // !!! Layer level
  geometryHierarchyIds_[3] = l.layerNumber();

  isTiltedLayer_ = l.isTilted();
  numRods_ = l.numRods();
  numFlatRings_ = l.buildNumModulesFlat();
  numRings_ = l.buildNumModules();
  if (!isTiltedLayer_) numRings_ = numFlatRings_;
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
  int side = m.uniRef().side;
  uint32_t ringRef;

  // Flat layer, or flat part of a tilted layer
  if (!m.isTilted()) {
    if (!isPixelTracker_) {
      // !!! Category level
      geometryHierarchyIds_[4] = 3;  // flat part : always 3.

      // !!! Rod level
      geometryHierarchyIds_[5] = phiRef_;
    }
    // !!! Inner Tracker : Rod level
    else { geometryHierarchyIds_[4] = phiRef_; }

    // Calculated ring identifier. 
    // Ring numbering starts from 1 from lowest Z, and increments with increasing Z.
    if (isCentered_) ringRef = (side > 0 ? (m.uniRef().ring + numFlatRings_ - 1) : (1 + numFlatRings_ - m.uniRef().ring));
    else ringRef = (side > 0 ? (m.uniRef().ring + numFlatRings_) : (1 + numFlatRings_ - m.uniRef().ring));

    // !!! Module level
    if (!isPixelTracker_) { geometryHierarchyIds_[6] = ringRef; }

    // !!! Inner Tracker : Ring level
    else { geometryHierarchyIds_[5] = ringRef; }
  }

  // Tilted part
  else {
    if (!isPixelTracker_) {
      // !!! Category level
      uint32_t category = (side > 0 ? 2 : 1); // tilted part : 2 on +Z side, 1 on -Z side.
      geometryHierarchyIds_[4] = category;

      // !!! Ring level
      ringRef = (side > 0 ? (m.uniRef().ring - numFlatRings_) : (1 + numRings_ - m.uniRef().ring));
      geometryHierarchyIds_[5] = ringRef;

      // !!! Module level
      geometryHierarchyIds_[6] = phiRef_;
    }
    // No scheme for tilted Inner Tracker exists in CMSSW.
    else { logWARNING("Tilted Pixel DetIds not supported yet."); }
  }

  // !!! Assign a null ref to sensor level.
  uint32_t sensorRef = 0;
  if (!isPixelTracker_) { geometryHierarchyIds_[7] = sensorRef; }
  else { geometryHierarchyIds_[6] = sensorRef; };

  // NOW THAT ALL NECESSARY REFS ARE COMPUTED, BUILD MODULE DETID !!
  m.buildDetId(geometryHierarchyIds_, geometryHierarchySizes_);
}

void BarrelDetIdBuilder::visit(Sensor& s) {
  if (s.subdet() == ModuleSubdetector::BARREL) {

    // !!! Sensor level
    if (!isPixelTracker_) {
      uint32_t sensorRef = (s.innerOuter() == SensorPosition::LOWER ? 1 : 2); // Lower sensor 1, upper sensor 2.
      geometryHierarchyIds_[7] = sensorRef;
    }
    // !!! Inner tracker : Sensor level
    else { geometryHierarchyIds_[6] = 0; }

    // NOW THAT ALL NECESSARY REFS ARE COMPUTED, BUILD SENSOR DETID !!
    s.buildDetId(geometryHierarchyIds_, geometryHierarchySizes_);
  }  
}



    //************************************//
    //*               Visitor             //
    //*          EndcapDetIdBuilder       //
    //*                                   //
    //************************************//

EndcapDetIdBuilder::EndcapDetIdBuilder(bool isPixelTracker, std::vector<int> geometryHierarchySizes) : 
  isPixelTracker_(isPixelTracker),
  geometryHierarchySizes_(geometryHierarchySizes)
{}

void EndcapDetIdBuilder::visit(Endcap& e) {
  // !!! Tracker level
  geometryHierarchyIds_[0] = 1;   // always 1.

  // !!! Endcap level
  if (!isPixelTracker_) { geometryHierarchyIds_[1] = 204 % 100; } // Outer Tracker : always 4. 
                                                       // In CMSSW, Phase 2 Outer Tracker Endcap is identified by 204, and its id is 204 % 100.
  else { geometryHierarchyIds_[1] = 202 % 100; }                  // Inner Tracker : always 2. 
                                                       // In CMSSW, Phase 2 Inner Tracker Endcap is identified by 202, and its id is 202 % 100.
}

void EndcapDetIdBuilder::visit(Disk& d) {
  bool side = d.side();

  // !!! Z side level
  uint32_t sideRef = (side ? 2 : 1);  // 2 for Z+ side, 1 for -Z side
  geometryHierarchyIds_[2] = sideRef;

  // !!! Unused hierarchy level	
  geometryHierarchyIds_[3] = 0;   // always 0

  // !!! Disk level
  uint32_t diskRef = d.diskNumber();
  geometryHierarchyIds_[4] = diskRef;

  numEmptyRings_ = d.numEmptyRings();
}
   
void EndcapDetIdBuilder::visit(Ring& r) {
  // !!! Ring level
  uint32_t ringRef = r.myid() - numEmptyRings_;
  geometryHierarchyIds_[5] = ringRef;

  // !!! Unused hierarchy level	
  geometryHierarchyIds_[6] = 1;   // always 1

  numModules_ = r.numModules();
}
 
void EndcapDetIdBuilder::visit(EndcapModule& m) {
  double phiSegment = 2 * M_PI / numModules_;                // Phi interval between 2 consecutive modules.
  double startAngle = femod( m.center().Phi(), phiSegment);  // Calculates the smallest positive Phi in the Ring.

  // Calculates Phi identifier.
  uint32_t phiRef_ = 1 + round(femod(m.center().Phi() - startAngle, 2*M_PI) / phiSegment);
  geometryHierarchyIds_[7] = phiRef_;
  // First module with Phi > 0 has phiRef 1. 
  // Next module (with bigger Phi) has phiRef 2. 
  // phiRef increments with increasing Phi.

  // !!! Assign a null ref to sensor level.
  uint32_t sensorRef = 0;
  geometryHierarchyIds_[8] = sensorRef;

  // NOW THAT ALL NECESSARY REFS ARE COMPUTED, BUILD MODULE DETID !!
  m.buildDetId(geometryHierarchyIds_, geometryHierarchySizes_);
}

void EndcapDetIdBuilder::visit(Sensor& s) {  
  if (s.subdet() == ModuleSubdetector::ENDCAP) {

    // !!! Sensor level
    if (!isPixelTracker_) {
      uint32_t sensorRef = (s.innerOuter() == SensorPosition::LOWER ? 1 : 2);  // Lower sensor 1, upper sensor 2.
      geometryHierarchyIds_[8] = sensorRef;
    }
    // !!! Inner tracker : Sensor level
    else { geometryHierarchyIds_[8] = 0; }

    // NOW THAT ALL NECESSARY REFS ARE COMPUTED, BUILD SENSOR DETID !!
    s.buildDetId(geometryHierarchyIds_, geometryHierarchySizes_);

    /*for (int a = 0; a < geometryHierarchyIds_.size(); a++) {
      std::cout << "values = " << std::endl;
      std::cout << geometryHierarchyIds_.at(a) << std::endl;
      std::cout << "scheme = " << std::endl;
      std::cout << geometryHierarchySizes_.at(a) << std::endl;
      }*/
    //std::bitset<32> test(s.myDetId());
    //std::cout << s.myDetId() << " " << test << " " << "rho = " <<  s.hitPoly().getCenter().Rho() << " z = " <<  s.hitPoly().getCenter().Z() << " phi = " <<  (s.hitPoly().getCenter().Phi() * 180. / M_PI) << std::endl;
  }
}

#include "Disk.h"
#include "messageLogger.h"
#include "ConversionStation.h"

const std::vector<double> Disk::scanSmallDeltas() const {
  std::vector<double> ringsSmallDeltas;

  for (int i = 1; i <= numRings(); i++) {
    Ring* ringTemplate = GeometryFactory::make<Ring>();
    ringTemplate->store(propertyTree());
    if (ringNode.count(i) > 0) ringTemplate->store(ringNode.at(i));
    ringsSmallDeltas.push_back(ringTemplate->smallDelta());
  }
  return ringsSmallDeltas;
}


const std::vector<double> Disk::scanDsDistances() const {
  std::vector<double> ringsDsDistances;

  for (int i = 1; i <= numRings(); i++) {
    RectangularModule* modTemplate = GeometryFactory::make<RectangularModule>();
    modTemplate->store(propertyTree());
    if (ringNode.count(i) > 0) modTemplate->store(ringNode.at(i));
    ringsDsDistances.push_back(modTemplate->dsDistance());
  }
  return ringsDsDistances;
}


const ScanDiskInfo Disk::scanPropertyTree() const {

  const std::vector<double>& ringsSmallDeltas = scanSmallDeltas();
  const std::vector<double>& ringsDsDistances = scanDsDistances();
  const ScanDiskInfo& diskInfo = std::make_pair(ringsSmallDeltas, ringsDsDistances);
  return diskInfo; 
}


void Disk::check() {
  PropertyObject::check();
  
  if (numRings.state() && innerRadius.state()) throw PathfulException("Only one between numRings and innerRadius can be specified");
  else if (!numRings.state() && !innerRadius.state()) throw PathfulException("At least one between numRings and innerRadius must be specified"); 
}

double Disk::getSmallDelta(const vector<double>& diskSmallDeltas, int ringNumber) const {
  if (ringNumber > diskSmallDeltas.size()) throw PathfulException("Tries to access information from ring which is not in disk.");
  return diskSmallDeltas.at(ringNumber - 1);
}

double Disk::getDsDistance(const vector<double>& diskDsDistances, int ringNumber) const {
  if (ringNumber > diskDsDistances.size()) throw PathfulException("Tries to access information from ring which is not in disk.");
  return diskDsDistances.at(ringNumber - 1);
}

void Disk::cutAtEta(double eta) { 
  for (auto& r : rings_) r.cutAtEta(eta); 
  //  rings_.erase_if([](const Ring& r) { return r.numModules() == 0; }); // get rid of rods which have been completely pruned
  rings_.erase_if([](const Ring& r) { return r.modules().size() == 0; }); // get rid of rods which have been completely pruned
  numRings(rings_.size());
  minZ.clear();
  minR.clear();
  maxZ.clear();
  maxR.clear();
}



std::pair<double> Disk::computeStringentZ(int i, int parity, const ScanEndcapInfo& extremaDisksInfo) {
  double lastZ;     // Z of the most stringent point in Ring i+1
  double newZ;      // Z of the most stringent point in Ring i

  const ScanDiskInfo& innermostDiskInfo = extremaDisksInfo.first;
  const vector<double>& innermostDiskSmallDeltas = innermostDiskInfo.first;
  const vector<double>& innermostDiskDsDistances = innermostDiskInfo.second;

  const ScanDiskInfo& outermostDiskInfo = extremaDisksInfo.second;
  const vector<double>& outermostDiskSmallDeltas = outermostDiskInfo.first;
  const vector<double>& outermostDiskDsDistances = outermostDiskInfo.second;

  double lastSmallDelta; // smallDelta (Ring i+1)
  double newSmallDelta;   // smallDelta (Ring i)
  double lastDsDistance; // dsDistance (Ring i+1)
  double newDsDistance;  // dsDistance (Ring i)
     
  // Case where Ring (i+1) is the innermost ring, and Ring (i) is the outermost ring.
  if (parity > 0) {
    // In this case, the innermost disk of the Endcaps block is the most stringent.
    // As a result, geometry information is taken from that disk.
    lastSmallDelta = getSmallDelta(innermostDiskSmallDeltas, i+1);
    newSmallDelta = getSmallDelta(innermostDiskSmallDeltas, i);
    lastDsDistance = getDsDistance(innermostDiskDsDistances, i+1);
    newDsDistance = getDsDistance(innermostDiskDsDistances, i);

    // Calculates Z of the most stringent points
    lastZ = buildZ() - zHalfLength() - bigDelta() - lastSmallDelta - lastDsDistance / 2.;
    newZ  = buildZ() - zHalfLength() + bigDelta() + newSmallDelta + newDsDistance / 2.;
  }

  // Case where Ring (i+1) is the outermost ring, and Ring (i) is the innermost ring.
  else {
    // In this case, the outermost disk of the Endcaps block is the most stringent.
    // As a result, geometry information is taken from that disk.
    lastSmallDelta = getSmallDelta(outermostDiskSmallDeltas, i+1);
    newSmallDelta = getSmallDelta(outermostDiskSmallDeltas, i);
    lastDsDistance = getDsDistance(outermostDiskDsDistances, i+1);
    newDsDistance = getDsDistance(outermostDiskDsDistances, i);

    // Calculates Z of the most stringent points
    lastZ = buildZ() + zHalfLength() + bigDelta() - lastSmallDelta - lastDsDistance / 2.;
    newZ  = buildZ() + zHalfLength() - bigDelta() + newSmallDelta + newDsDistance / 2.;
  }

  return std::make_pair(lastZ, newZ);
}


void Disk::buildTopDown(const ScanEndcapInfo& extremaDisksInfo) {

  double lastRho;
  
  for (int i = numRings(), parity = -bigParity(); i > 0; i--, parity *= -1) {

    // CREATES RING
    Ring* ring = GeometryFactory::make<Ring>();
    ring->buildDirection(Ring::TOPDOWN);
    ring->store(propertyTree());
    if (ringNode.count(i) > 0) ring->store(ringNode.at(i));

    // INITIALIZATION
    if (i == numRings()) {
      ring->buildStartRadius(outerRadius());
    }
    
    // ALL THIS IS TO AUTOMATICALLY CALCULATE RING RADIUSHIGH FOR RING (i) FROM RIMG (i+1)
    // In what follows, 'innermost' corresponds to 'lower Z', and 'outermost' corresponds to 'bigger Z'.
    else {

      // 1) FIND THE Z OF THE MOST STRINGENT POINTS IN RING (i+1) AND RING (i)
      std::pair<double> stringentZ = computeStringentZ(i, parity, extremaDisksInfo);
      double lastZ = stringentZ.first;           // Z of the most stringent point in Ring i+1
      double newZ = stringentZ.second;           // Z of the most stringent point in Ring i
      

      // 2) CALCULATES RING RADIUSHIGH FOR RING (i) FROM RING (i+1)
      
      // Case A : Consider zError
      double zError   = parity > 0 ? +zError() : -zError();
      double nextRhoWithZError   = lastRho / (lastZ - zError) * (newZ - zError);

      // Case B : Consider rOverlap
      double nextRhoWithROverlap  = (lastRho + rOverlap()) / lastZ * newZ;
      
      // Takes the most stringent of cases A and B
      double nextRho = MAX(nextRhoWithZError, nextRhoWithROverlap);

      ring->buildStartRadius(nextRho);
    }

    // NOW THAT RADIUS HAS BEEN CALCULATED, FINISHES BUILDING THE RING AND STORE IT
    ring->myid(i);
    ring->build();
    ring->translateZ(parity > 0 ? bigDelta() : -bigDelta());
    rings_.insert(rings_.begin(), ring);
    ringIndexMap_[i] = ring;

    // Keep for next calculation
    lastRho        = ring->minR();
  }
}


void Disk::build(const ScanEndcapInfo& extremaDisksInfo) {
  ConversionStation* conversionStation;
  materialObject_.store(propertyTree());
  materialObject_.build();

  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    if (numRings.state()) buildTopDown(extremaDisksInfo);
    translateZ(placeZ());
  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  for (auto& currentStationNode : stationsNode) {
    conversionStation = new ConversionStation();
    conversionStation->store(currentStationNode.second);
    conversionStation->check();
    conversionStation->build();
      
    if(conversionStation->stationType() == ConversionStation::Type::FLANGE) {
      if(flangeConversionStation_ == nullptr) { //take only first defined flange station
        flangeConversionStation_ = conversionStation;
      }
    } else if(conversionStation->stationType() == ConversionStation::Type::SECOND) {
      secondConversionStations_.push_back(conversionStation);
    }
  }

  cleanup();
  builtok(true);
}


void Disk::translateZ(double z) { averageZ_ += z; for (auto& r : rings_) r.translateZ(z); }

void Disk::mirrorZ() {
  averageZ_ = -averageZ_;
  for (auto& r : rings_) r.mirrorZ();
}

const MaterialObject& Disk::materialObject() const {
  return materialObject_;
}

ConversionStation* Disk::flangeConversionStation() const {
  return flangeConversionStation_;
}

const std::vector<ConversionStation*>& Disk::secondConversionStations() const {
  return secondConversionStations_;
}

std::vector<std::set<const Module*>> Disk::getModuleSurfaces() const {
  class SurfaceSplitterVisitor : public ConstGeometryVisitor {
  private:
    int sideIndex;
    double ringAvgZ;
    int iRing=0;
  public:
    double diskAverageZ;
    std::set<const Module*> mod[4];
    void visit(const Ring& r) {
      sideIndex = 0;
      ringAvgZ = r.averageZ();
      if (ringAvgZ>diskAverageZ) sideIndex+=2;
    }
    void visit(const Module& m) {
      int deltaSide=0;
      if (m.center().Z()>ringAvgZ) deltaSide+=1;
      mod[sideIndex+deltaSide].insert(&m);
    }
  };
  SurfaceSplitterVisitor v;
  v.diskAverageZ = averageZ_;
  this->accept(v);
  std::vector<std::set<const Module*>> result;
  for (int i=0; i<4; ++i) result.push_back(v.mod[i]);
  return result;
}

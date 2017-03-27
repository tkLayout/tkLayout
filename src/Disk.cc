#include "Disk.hh"
#include "MessageLogger.hh"
#include "ConversionStation.hh"

/** Scan Property Tree and returns the vector of rings' smallDeltas.
 */
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

/** Scan Property Tree and returns the vector of rings' dsDistances.
 */
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

/** Scan Property Tree and gathers relevant info which is needed for building a disk.
 */
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


const double Disk::getRingInfo(const vector<double>& ringsInfo, int ringNumber) const {
  if (ringNumber > ringsInfo.size()) {
    throw PathfulException(Form("When building disk, tried to access information from ring %i", ringNumber));
  }
  return ringsInfo.at(ringNumber - 1);
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


/** This is used for building ring (i+1) from ring (i).
    It returns the Z of the most stringent points in Ring (i+1) and Ring (i).
    These Z belong either to the innermost, either to the outermost Disk in the Endcap.
    In what follows, 'innermost' corresponds to 'lower Z', and 'outermost' corresponds to 'bigger Z' (Z+ side).
    * @param i : ringNumber
    * @param parity : ring's parity
    * @param extremaDisksInfo : info from the innermost and outermost disks in the Endcap
    * return : Z in Ring (i+1), and Z in Ring (i)
*/
std::pair<double, double> Disk::computeStringentZ(int i, int parity, const ScanEndcapInfo& extremaDisksInfo) {
  double lastZ;     // Z of the most stringent point in Ring (i+1)
  double newZ;      // Z of the most stringent point in Ring (i)

  double lastSmallDelta; // smallDelta (Ring i+1)
  double newSmallDelta;   // smallDelta (Ring i)
  double lastDsDistance; // dsDistance (Ring i+1)
  double newDsDistance;  // dsDistance (Ring i)
     
  // Case where Ring (i+1) is the innermost ring, and Ring (i) is the outermost ring.
  if (parity > 0) {
    // IN THIS CASE, THE INNERMOST DISK OF THE ENDCAPS BLOCK IS THE MOST STRINGENT.
    // AS A RESULT, GEOMETRY INFORMATION IS TAKEN FROM THAT DISK.
    const ScanDiskInfo& innermostDiskInfo = extremaDisksInfo.first;
    const vector<double>& innermostDiskSmallDeltas = innermostDiskInfo.first;
    const vector<double>& innermostDiskDsDistances = innermostDiskInfo.second;

    lastSmallDelta = getRingInfo(innermostDiskSmallDeltas, i+1);
    newSmallDelta = getRingInfo(innermostDiskSmallDeltas, i);
    lastDsDistance = getRingInfo(innermostDiskDsDistances, i+1);
    newDsDistance = getRingInfo(innermostDiskDsDistances, i);

    // Calculates Z of the most stringent points
    lastZ = buildZ() - zHalfLength() - bigDelta() - lastSmallDelta - lastDsDistance / 2.;
    newZ  = buildZ() - zHalfLength() + bigDelta() + newSmallDelta + newDsDistance / 2.;
  }

  // Case where Ring (i+1) is the outermost ring, and Ring (i) is the innermost ring.
  else {
    // IN THIS CASE, THE OUTERMOST DISK OF THE ENDCAPS BLOCK IS THE MOST STRINGENT.
    // AS A RESULT, GEOMETRY INFORMATION IS TAKEN FROM THAT DISK.
    const ScanDiskInfo& outermostDiskInfo = extremaDisksInfo.second;
    const vector<double>& outermostDiskSmallDeltas = outermostDiskInfo.first;
    const vector<double>& outermostDiskDsDistances = outermostDiskInfo.second;

    lastSmallDelta = getRingInfo(outermostDiskSmallDeltas, i+1);
    newSmallDelta = getRingInfo(outermostDiskSmallDeltas, i);
    lastDsDistance = getRingInfo(outermostDiskDsDistances, i+1);
    newDsDistance = getRingInfo(outermostDiskDsDistances, i);

    // Calculates Z of the most stringent points
    lastZ = buildZ() + zHalfLength() + bigDelta() - lastSmallDelta - lastDsDistance / 2.;
    newZ  = buildZ() + zHalfLength() - bigDelta() + newSmallDelta + newDsDistance / 2.;
  }

  return std::make_pair(lastZ, newZ);
}


/** Calculates ring (i) radiusHigh, using ring (i+1).
    The most stringent of zError and rOverlap is used.
 */
double Disk::computeNextRho(int parity, double lastZ, double newZ, double lastRho) {

  // Case A : Consider zError
  double zErrorShift   = (parity > 0 ? zError() : - zError());
  double nextRhoWithZError   = lastRho / (lastZ - zErrorShift) * (newZ - zErrorShift);

  // Case B : Consider rOverlap
  double nextRhoWithROverlap  = (lastRho + rOverlap()) / lastZ * newZ;
      
  // Takes the most stringent of cases A and B
  double nextRho = MAX(nextRhoWithZError, nextRhoWithROverlap);
  return nextRho;
}



/** In disk, build ring (i) from ring (i+1).
    i is ringNumber.
 */
void Disk::buildTopDown(const ScanEndcapInfo& extremaDisksInfo) {

  double lastRho;
  
  for (int i = numRings(), parity = -bigParity(); i > 0; i--, parity *= -1) {

    // CREATES RING (NOT PLACED AT ANY RADIUS YET)
    Ring* ring = GeometryFactory::make<Ring>();
    ring->buildDirection(Ring::TOPDOWN);
    ring->store(propertyTree());
    if (ringNode.count(i) > 0) ring->store(ringNode.at(i));

    // INITIALIZATION
    if (i == numRings()) {
      ring->buildStartRadius(outerRadius());
    }
    
    // ALL THIS IS TO AUTOMATICALLY CALCULATE RING RADIUSHIGH FOR RING (i) FROM RIMG (i+1)
    else {

      // 1) FIND THE Z OF THE MOST STRINGENT POINTS IN RING (i+1) AND RING (i)
      std::pair<double, double> stringentZ = computeStringentZ(i, parity, extremaDisksInfo);
      double lastZ = stringentZ.first;           // Z of the most stringent point in Ring (i+1)
      double newZ = stringentZ.second;           // Z of the most stringent point in Ring (i)
      
      // 2) CALCULATES RING (i) RADIUSHIGH USING RING (i+1)
      double nextRho = computeNextRho(parity, lastZ, newZ, lastRho);

      // 3) NOW, CAN ASSIGN THE CALCULATED RADIUS TO RING (i) ! 
      ring->buildStartRadius(nextRho);
    }

    // NOW THAT RADIUS HAS BEEN CALCULATED, FINISH BUILDING THE RING AND STORE IT
    ring->myid(i);
    ring->build();
    ring->translateZ(parity > 0 ? bigDelta() : -bigDelta());
    rings_.insert(rings_.begin(), ring);
    ringIndexMap_[i] = ring;

    // Keep for next calculation
    lastRho = ring->minR();
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

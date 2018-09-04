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

/** Scan Property Tree and returns the rings' sensorThickness.
 */
const double Disk::scanSensorThickness() const {
  double ringsSensorThickness;

  int i = numRings(); // Look at only 1 Ring, assumes sensorThickness is identical everywhere in the Disk!
  RectangularModule* modTemplate = GeometryFactory::make<RectangularModule>();
  modTemplate->store(propertyTree());
  if (ringNode.count(i) > 0) modTemplate->store(ringNode.at(i));
  ringsSensorThickness = modTemplate->sensorThickness();

  return ringsSensorThickness;
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
  const int numRings = ringsInfo.size();
  if (ringNumber > numRings) {
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
    rSafetyMargin is also taken into account.
 */
double Disk::computeNextRho(const int parity, const double zError, const double rSafetyMargin, const double lastZ, const double newZ, const double lastRho, const double oneBeforeLastRho) {

  // Case A: Consider zError.
  double zErrorShift   = (parity > 0 ? zError : - zError);
  double nextRho = lastRho / (lastZ - zErrorShift) * (newZ - zErrorShift);

  // Case B : Consider rOverlap
  if (rOverlap.state()) {
    double nextRhoWithROverlap  = (lastRho + rOverlap()) / lastZ * newZ;
    // Takes the most stringent of cases A and B
    nextRho = MAX(nextRho, nextRhoWithROverlap);
  }

  // If relevant, consider rSafetyMargin.
  if (oneBeforeLastRho > 1.) {
    double nextRhoSafe = oneBeforeLastRho - rSafetyMargin;
    nextRho = MIN(nextRho, nextRhoSafe);
  }    

  return nextRho;
}



/** In disk, build ring (i) from ring (i+1).
    i is ringNumber.
 */
void Disk::buildTopDown(const ScanEndcapInfo& extremaDisksInfo) {

  double oneBeforeLastRho = 0.;
  double lastRho = 0.;

  for (int i = numRings(), parity = -bigParity(); i > 0; i--, parity *= -1) {

    // CREATES RING (NOT PLACED AT ANY RADIUS YET)
    Ring* ring = GeometryFactory::make<Ring>();
    ring->buildDirection(Ring::TOPDOWN);
    ring->store(propertyTree());
    if (ringNode.count(i) > 0) ring->store(ringNode.at(i));
    ring->subdetectorName(subdetectorName());

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
      double zError = ring->zError();
      double rSafetyMargin = ring->rSafetyMargin();
      double nextRho = computeNextRho(parity, zError, rSafetyMargin, lastZ, newZ, lastRho, oneBeforeLastRho);

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
    oneBeforeLastRho = lastRho;
    lastRho = ring->minR();
  }
}


void Disk::build(const ScanEndcapInfo& extremaDisksInfo) {
  ConversionStation* conversionStation;
  materialObject_.store(propertyTree());
  materialObject_.matSubdetectorName(subdetectorName());
  std::cout << "Disk::build()  : materialObject_.matSubdetectorName() = " << materialObject_.matSubdetectorName() << std::endl;
  materialObject_.build();

  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    if (numRings.state()) buildTopDown(extremaDisksInfo);
    translateZ(placeZ());
  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  for (auto& currentStationNode : stationsNode) {
    conversionStation = new ConversionStation();
    conversionStation->store(propertyTree());
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


const std::map<int, std::vector<const Module*> > Disk::getSurfaceModules() const {

  class SurfaceSplitterVisitor : public ConstGeometryVisitor {
  public:
    const std::map<int, std::vector<const Module*> > surfaceModules() const { return surfaceModules_; }
    void visit(const EndcapModule& m) {
      surfaceModules_[m.diskSurface()].push_back(&m);    
    }
  private:
    std::map<int, std::vector<const Module*> > surfaceModules_;
  };

  SurfaceSplitterVisitor v;
  accept(v);
  const std::map<int, std::vector<const Module*> >& surfaceModules = v.surfaceModules();
  return surfaceModules;
}


/** Binds the 2 points which are provided as arguments, and returns corresponding Z of intersection with (Z) axis.
 */
const std::pair<double, bool> Disk::computeIntersectionWithZAxis(double lastZ, double lastRho, double newZ, double newRho) const {
  // Slope of the line that binds the most stringent points of ring (i) and ring (i+1).
  double slope = (newRho - lastRho) / (newZ - lastZ);
  bool isPositiveSlope = (slope > 0.);

  // Used to calculate the coverage in Z of ring (i) with respect to ring (i+1).
  double zIntersection;
  zIntersection = newZ - newRho / slope;  // Intersection of the line with (Z) axis.

  return std::make_pair(zIntersection, isPositiveSlope);
}


/** This computes the actual coverage in Z of a disk (after it is built).
    It calculates the actual zError, using the relevant coordinates of the disk.
*/
void Disk::computeActualZCoverage() {
  
  double lastMinRho;
  double lastMinZ;
  const double ringsSensorThickness = scanSensorThickness();

  for (int i = numRings(), parity = -bigParity(); i > 0; i--, parity *= -1) {

    if (i != numRings()) {
      // CALCULATE THE COVERAGE IN Z OF RING (i) WITH RESPECT TO RING (i+1).
      // RING (i+1) HAS BIGGER RADIUS, ring (i) HAS SMALLER RADIUS.
     
      // Find the radii and Z of the most stringent points in ring (i).
      double newMaxRho = rings_.at(i-1).buildStartRadius();
      double newMaxZ = rings_.at(i-1).maxZ() - ringsSensorThickness / 2.;

      // Calculation : Min coordinates of ring (i+1) with max coordinates of ring (i)
      std::pair<double, bool> intersectionWithZAxis = computeIntersectionWithZAxis(lastMinZ, lastMinRho, newMaxZ, newMaxRho);
      double zIntersection = intersectionWithZAxis.first;
      bool isPositiveSlope = intersectionWithZAxis.second;
      
      double zErrorCoverage = 0.;
      // CASE WHERE RING (i+1) HAS SMALLER Z, AND RING (i) HAS BIGGER Z.
      if (parity > 0.) {
	if (isPositiveSlope) zErrorCoverage = zIntersection;
	else zErrorCoverage = -std::numeric_limits<double>::infinity();
      }

      // CASE WHERE RING (i+1) HAS BIGGER Z, AND RING (i) HAS SMALLER Z.
      else {
	if (isPositiveSlope) zErrorCoverage = -zIntersection;
	else zErrorCoverage = std::numeric_limits<double>::infinity();
      }
      
      // STORE THE RESULT
      rings_.at(i-1).actualZError(zErrorCoverage);
      ringIndexMap_[i]->actualZError(zErrorCoverage);
    }

    // Keep for next calculation : radii and Z of the most stringent point in ring (i+1).
    lastMinRho = rings_.at(i-1).minR();
    lastMinZ = rings_.at(i-1).minZ() + ringsSensorThickness / 2.;
  }
}


/** This computes the actual coverage of a disk (after it is built).
 */
void Disk::computeActualCoverage() {
  // Actual coverage in Z
  computeActualZCoverage();

  // Actual coverage in Phi
  for (auto& r : rings_) r.computeActualPhiCoverage();
}

void Disk::translateZ(double z) { averageZ_ += z; for (auto& r : rings_) r.translateZ(z); }

void Disk::rotateToNegativeZSide() {
  averageZ_ = -averageZ_;
  for (auto& r : rings_) r.rotateToNegativeZSide();
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

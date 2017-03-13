#include "Disk.h"
#include "messageLogger.h"
#include "ConversionStation.h"

const std::vector<double> Disk::getSmallDeltasFromTree() const {
  std::vector<double> smallDeltas;

  for (int i = 1; i <= numRings(); i++) {
    Ring* ringTemplate = GeometryFactory::make<Ring>();
    ringTemplate->store(propertyTree());
    if (ringNode.count(i) > 0) ringTemplate->store(ringNode.at(i));
    smallDeltas.push_back(ringTemplate->smallDelta());
  }
  return smallDeltas;
}

const std::vector<double> Disk::getDsDistancesFromTree() const {
  std::vector<double> dsDistances;

  for (int i = 1; i <= numRings(); i++) {
    Ring* ringTemplate = GeometryFactory::make<Ring>();
    ringTemplate->store(propertyTree());
    if (ringNode.count(i) > 0) ringTemplate->store(ringNode.at(i));

    RectangularModule* modTemplate = GeometryFactory::make<RectangularModule>();
    modTemplate->store(propertyTree());
    if (ringNode.count(i) > 0) modTemplate->store(ringNode.at(i));
    dsDistances.push_back(modTemplate->dsDistance());
  }
  return dsDistances;
}

void Disk::check() {
  PropertyObject::check();
  
  if (numRings.state() && innerRadius.state()) throw PathfulException("Only one between numRings and innerRadius can be specified");
  else if (!numRings.state() && !innerRadius.state()) throw PathfulException("At least one between numRings and innerRadius must be specified"); 
}

double Disk::getSmallDelta(const vector<double>& diskSmallDeltas, int ringIndex) const {
  //std::cout << " ringIndex = " << ringIndex << std::endl;
  //std::cout << "diskSmallDeltas.size() = " << diskSmallDeltas.size() << std::endl;
  if (ringIndex > diskSmallDeltas.size()) throw PathfulException("Tries to access information from ring which is not in disk.");
  return diskSmallDeltas.at(ringIndex - 1);
}

double Disk::getDsDistance(const vector<double>& diskDsDistances, int ringIndex) const {
  //std::cout << " ringIndex = " << ringIndex << std::endl;
  //std::cout << "diskDsDistances.size() = " << diskDsDistances.size() << std::endl;
  if (ringIndex > diskDsDistances.size()) throw PathfulException("Tries to access information from ring which is not in disk.");
  return diskDsDistances.at(ringIndex - 1);
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

/*void Disk::buildBottomUp(const vector<double>& buildDsDistances) {
//  auto childmap = propertyTree().getChildMap<int>("Ring");
  double lastRho = 0;
  double lastSmallDelta;

  for (int i = 1, parity = -bigParity(); lastRho < outerRadius(); i++, parity *= -1) {
    Ring* ring = GeometryFactory::make<Ring>();
    ring->buildDirection(Ring::BOTTOMUP);
    ring->buildCropRadius(outerRadius());
    ring->store(propertyTree());
    if (ringNode.count(i) > 0) ring->store(ringNode.at(i)); //childmap.count(i) > 0 ? childmap[i] : propertyTree()
    if (i == 1) {
      ring->buildStartRadius(innerRadius());
    } else {

      // Calcaluate new and last position in Z
      double newZ  = buildZ() + (parity > 0 ? + bigDelta() : - bigDelta()) - ring->smallDelta() - getDsDistance(buildDsDistances, i)/2;  
      double lastZ = buildZ() + (parity > 0 ? - bigDelta() : + bigDelta()) + lastSmallDelta     + getDsDistance(buildDsDistances, i)/2;
      //double originZ = parity > 0 ? -zError() : +zError();



      // Calculate shift in Z position of extreme cases (either innemost or outermost disc)
      // Remember that disc put always in to the centre of endcap
      double shiftZ        = parity > 0 ? -zHalfLength() : +zHalfLength();

      // Calculate shift in Z due to beam spot
      double errorShiftZ   = parity > 0 ? -zError() : +zError();

      // Calculate next rho taking into account overlap in extreme cases of innermost or outermost disc
      double nextRhoWOverlap  = (lastRho - rOverlap())/(lastZ - shiftZ)*(newZ - shiftZ);
      // Calculate next rho taking into account overlap zError in extreme cases of innermost or outermost disc
      double nextRhoWZError   = (lastRho)/(lastZ - shiftZ - errorShiftZ)*(newZ - shiftZ - errorShiftZ);

      // New next rho as a minimum
      double nextRho = MIN(nextRhoWOverlap, nextRhoWZError);
      ring->buildStartRadius(nextRho);

       // For debug test only
       //std::cout << ">>> noOverlap:        " << (lastRho)/lastZ * newZ              << " New [ ; ]: " << (lastRho)/(lastZ - zHalfLength()) * (newZ - zHalfLength()) << " " << (lastRho)/(lastZ + zHalfLength()) * (newZ + zHalfLength()) << std::endl;
       //std::cout << ">>> yesOverlap:       " << (lastRho - rOverlap())/lastZ * newZ << " New [ ; ]: " << (lastRho - rOverlap())/(lastZ - zHalfLength()) * (newZ - zHalfLength()) << " " << (lastRho - rOverlap())/(lastZ + zHalfLength()) * (newZ + zHalfLength()) << std::endl;
       //std::cout << ">>> yesZErr:          " << (lastRho)/lastZ * newZ              << " New [ ; ]: " << (lastRho)/(lastZ - zError() - zHalfLength()) * (newZ - zError() - zHalfLength()) << " " << (lastRho)/(lastZ + zError() + zHalfLength()) * (newZ + zError() + zHalfLength()) << std::endl;
       //
       //std::cout << ">>> ShiftZ: " << shiftZ << " ErrorZ: " << errorShiftZ << " nR: "<< nextRhoWOverlap << " nRShifted: " << nextRhoWZError << std::endl;
    }
    ring->myid(i);
    ring->build();
    ring->translateZ(parity > 0 ? bigDelta() : -bigDelta());
    rings_.push_back(ring);
    //ringIndexMap_.insert(i, ring);
    ringIndexMap_[i] = ring;

    // Keep for next calculation
    lastRho        = ring->maxR();
    lastSmallDelta = ring->smallDelta();
  }
  numRings(rings_.size());
}*/

void Disk::buildTopDown(const vector<double>& firstDiskSmallDeltas, const vector<double>& lastDiskSmallDeltas, const vector<double>& firstDiskDsDistances, const vector<double>& lastDiskDsDistances) {
  double lastRho;
  //double lastSmallDelta;
  for (int i = numRings(), parity = -bigParity(); i > 0; i--, parity *= -1) {
    Ring* ring = GeometryFactory::make<Ring>();
    ring->buildDirection(Ring::TOPDOWN);
    ring->store(propertyTree());
    if (ringNode.count(i) > 0) ring->store(ringNode.at(i));

    if (i == numRings()) {
      ring->buildStartRadius(outerRadius());
    } else {

      double lastZ, newZ;
      // Calcalate new and last position in Z
      if (parity > 0) {
	lastZ = buildZ() - zHalfLength() - bigDelta() - getSmallDelta(firstDiskSmallDeltas, i+1) - getDsDistance(firstDiskDsDistances, i+1) / 2.;
	newZ  = buildZ() - zHalfLength() + bigDelta() + getSmallDelta(firstDiskSmallDeltas, i) + getDsDistance(firstDiskDsDistances, i) / 2.;

	std::cout << "OOOOOOO firstDiskSmallDeltas.size() = " << firstDiskSmallDeltas.size() << std::endl;

	std::cout << "!!!!!! DISK i = " << myid() << std::endl;
	std::cout << "!!!!!! RING i = " << i << std::endl;
	std::cout << "parity = " << parity << std::endl;
	std::cout << "lastZ = " << lastZ << std::endl;
	std::cout << "last smallDelta = " << getSmallDelta(firstDiskSmallDeltas, i+1) << std::endl;
	std::cout << "last dsDistance = " << getDsDistance(firstDiskDsDistances, i+1) << std::endl;
	std::cout << "newZ = " << newZ << std::endl;
	std::cout << "new small delta = " << getSmallDelta(firstDiskSmallDeltas, i) << std::endl;
	std::cout << "new dsDistance = " << getDsDistance(firstDiskDsDistances, i) << std::endl;


      }
      else {
	lastZ = buildZ() + zHalfLength() + bigDelta() - getSmallDelta(lastDiskSmallDeltas, i+1) - getDsDistance(lastDiskDsDistances, i+1) / 2.;
        newZ  = buildZ() + zHalfLength() - bigDelta() + getSmallDelta(lastDiskSmallDeltas, i) + getDsDistance(lastDiskDsDistances, i) / 2.;

	std::cout << "OOOOOOO lastDiskSmallDeltas.size() = " << lastDiskSmallDeltas.size() << std::endl;

	std::cout << "!!!!!! DISK i = " << myid() << std::endl;
	std::cout << "!!!!!! RING i = " << i << std::endl;
	std::cout << "parity = " << parity << std::endl;
	std::cout << "lastZ = " << lastZ << std::endl;
	std::cout << "last smallDelta = " << getSmallDelta(lastDiskSmallDeltas, i+1) << std::endl;
	std::cout << "last dsDistance = " << getDsDistance(lastDiskDsDistances, i+1) << std::endl;
	std::cout << "newZ = " << newZ << std::endl;
	std::cout << "new small delta = " << getSmallDelta(lastDiskSmallDeltas, i) << std::endl;
	std::cout << "new dsDistance = " << getDsDistance(lastDiskDsDistances, i) << std::endl;



      }


      // Calculate shift in Z position of extreme cases (either innemost or outermost disc)
      // Remember that disc put always in to the centre of endcap
      //double shiftZ = parity > 0 ? +zHalfLength() : -zHalfLength();

      // Calculate shift in Z due to beam spot
      double errorShiftZ   = parity > 0 ? +zError() : -zError();

      // Calculate next rho taking into account overlap in extreme cases of innermost or outermost disc
      double nextRhoWOverlap  = (lastRho + rOverlap()) / lastZ * newZ;
      // Calculate next rho taking into account overlap zError in extreme cases of innermost or outermost disc
      double nextRhoWZError   = lastRho / (lastZ - errorShiftZ) * (newZ - errorShiftZ);

      double nextRho = MAX(nextRhoWOverlap, nextRhoWZError);
      std::cout << "FINALLLLLLLL nextRho " << nextRho << std::endl;
      ring->buildStartRadius(nextRho);
    }

    ring->myid(i);
    ring->build();
    ring->translateZ(parity > 0 ? bigDelta() : -bigDelta());
    //rings_.push_back(ring);
    rings_.insert(rings_.begin(), ring);
    //ringIndexMap_.insert(i, ring);
    ringIndexMap_[i] = ring;

    // Keep for next calculation
    lastRho        = ring->minR();
    //lastSmallDelta = ring->smallDelta();
  }
}

void Disk::build(const vector<double>& firstDiskSmallDeltas, const vector<double>& lastDiskSmallDeltas, const vector<double>& firstDiskDsDistances, const vector<double>& lastDiskDsDistances) {
  ConversionStation* conversionStation;
  materialObject_.store(propertyTree());
  materialObject_.build();

  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    if (numRings.state()) buildTopDown(firstDiskSmallDeltas, lastDiskSmallDeltas, firstDiskDsDistances, lastDiskDsDistances);
    //else buildBottomUp(buildDsDistances);
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


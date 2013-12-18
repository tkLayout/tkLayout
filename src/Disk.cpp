#include "Disk.h"

void Disk::check() {
  PropertyObject::check();
  
  if (numRings.state() && innerRadius.state()) throw PathfulException("Only one between numRings and innerRadius can be specified");
  else if (!numRings.state() && !innerRadius.state()) throw PathfulException("At least one between numRings and innerRadius must be specified"); 
}

double Disk::getDsDistance(const vector<double>& buildDsDistances, int rindex) const {
  return rindex-1 >= buildDsDistances.size() ? buildDsDistances.back() : buildDsDistances.at(rindex-1);
}

void Disk::cutAtEta(double eta) { 
  for (auto& r : rings_) r.cutAtEta(eta); 
  rings_.erase_if([](const Ring& r) { return r.numModules() == 0; }); // get rid of rods which have been completely pruned
  numRings(rings_.size());
}

void Disk::buildBottomUp(const vector<double>& buildDsDistances) {
//  auto childmap = propertyTree().getChildMap<int>("Ring");
  double lastRho = 0;
  double lastSmallDelta;

  for (int i = 1, parity = -bigParity(); lastRho < outerRadius(); i++, parity *= -1) {
    Ring* ring = new Ring();
    ring->setup();
    ring->buildDirection(Ring::BOTTOMUP);
    ring->buildCropRadius(outerRadius());
    ring->store(propertyTree());
    if (ringNode.count(i) > 0) ring->store(ringNode.at(i)); /*childmap.count(i) > 0 ? childmap[i] : propertyTree()*/
    if (i == 1) {
      ring->buildStartRadius(innerRadius());
    } else {
      double newZ  = buildZ() + (parity > 0 ? + bigDelta() : - bigDelta()) - ring->smallDelta() - getDsDistance(buildDsDistances, i)/2;  
      double lastZ = buildZ() + (parity > 0 ? - bigDelta() : + bigDelta()) + lastSmallDelta + getDsDistance(buildDsDistances, i)/2;
      double originZ = parity > 0 ? -zError() : +zError();
      double nextRhoOrigin = (lastRho - minRingOverlap())/lastZ * newZ;
      double nextRhoShifted = lastRho/(lastZ - originZ) * (newZ - originZ);
      double nextRho = MIN(nextRhoOrigin, nextRhoShifted);
      ring->buildStartRadius(nextRho);
    }
    ring->myid(i);
    ring->build();
    ring->translateZ(parity > 0 ? bigDelta() : -bigDelta());
    rings_.push_back(ring);
    lastRho = ring->maxR();
    lastSmallDelta = ring->smallDelta();
  }
  numRings(rings_.size());
}

void Disk::buildTopDown(const vector<double>& buildDsDistances) {
  double lastRho;
  double lastSmallDelta;
  for (int i = numRings(), parity = -bigParity(); i > 0; i--, parity *= -1) {
    Ring* ring = new Ring();
    ring->setup();
    ring->buildDirection(Ring::TOPDOWN);
    ring->store(propertyTree());
    if (ringNode.count(i) > 0) ring->store(ringNode.at(i));
    if (i == numRings()) {
      ring->buildStartRadius(outerRadius());
    } else {
      double newZ  = buildZ() + (parity > 0 ? + bigDelta() : - bigDelta()) + ring->smallDelta() + getDsDistance(buildDsDistances, i)/2; // CUIDADO was + smallDelta + dsDistances[nRing-1]/2;
      double lastZ = buildZ() + (parity > 0 ? - bigDelta() : + bigDelta()) - lastSmallDelta - getDsDistance(buildDsDistances, i+1)/2; // CUIDADO was - smallDelta - dsDistances[nRing-1]/2; // try with prevRing here
      double originZ = parity > 0 ? zError() : -zError();
      double nextRhoOrigin = (lastRho + minRingOverlap())/lastZ * newZ;
      double nextRhoShifted = lastRho/(lastZ - originZ) * (newZ - originZ);
      double nextRho = MAX(nextRhoOrigin, nextRhoShifted);
      ring->buildStartRadius(nextRho);
    }
    ring->myid(i);
    ring->build();
    ring->translateZ(parity > 0 ? bigDelta() : -bigDelta());
    rings_.push_back(ring);
    lastRho = ring->minR();
    lastSmallDelta = ring->smallDelta();
  }
}

void Disk::build(const vector<double>& buildDsDistances) {
  try {
    std::cout << ">>> Building " << fullid(*this) << " <<<" << std::endl;
    if (numRings.state()) buildTopDown(buildDsDistances);
    else buildBottomUp(buildDsDistances);
    translateZ(placeZ());
  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }
  cleanup();
  builtok(true);
}


void Disk::translateZ(double z) { averageZ_ += z; for (auto& r : rings_) r.translateZ(z); }

void Disk::mirrorZ() {
  averageZ_ = -averageZ_;
  for (auto& r : rings_) r.mirrorZ();
}


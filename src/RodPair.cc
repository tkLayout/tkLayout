#include "RodPair.hh"
#include "MessageLogger.hh"

void RodPair::clearComputables() { 
}

void RodPair::translate(const XYZVector& translation) {
  for (auto& m : zPlusModules_) { m.translate(translation); }
  for (auto& m : zMinusModules_) { m.translate(translation); }
  clearComputables();
}

void RodPair::translateR(double radius) {
  for (auto& m : zPlusModules_) { m.translateR(radius); }
  for (auto& m : zMinusModules_) { m.translateR(radius); }
  clearComputables();
}

void RodPair::rotateZ(double angle) {
  for (auto& m : zPlusModules_) { m.rotateZ(angle); }
  for (auto& m : zMinusModules_) { m.rotateZ(angle); }
  clearComputables();
}

void RodPair::cutAtEta(double eta) { 
  zPlusModules_.erase_if([eta](BarrelModule& m) { return fabs(m.center().Eta()) > eta; }); 
  zMinusModules_.erase_if([eta](BarrelModule& m) { return fabs(m.center().Eta()) > eta; }); 
}

const MaterialObject& RodPair::materialObject() const{
  return materialObject_;
}

void StraightRodPair::compressToZ(double z) {

  if (zPlusModules_.empty()) return;

  logINFO("Layer slated for compression");

  auto findMaxZModule = [&]() { return *std::max_element(zPlusModules_.begin(), zPlusModules_.end(), [](const BarrelModule& m1, const BarrelModule& m2) { return m1.planarMaxZ() < m2.planarMaxZ(); }); };
  auto findMinZModule = [&]() { return *std::min_element(zMinusModules_.begin(), zMinusModules_.end(), [](const BarrelModule& m1, const BarrelModule& m2) { return m1.planarMinZ() < m2.planarMinZ(); }); };

  double Deltap =  fabs(z) - findMaxZModule().planarMaxZ();
  double Deltam = -fabs(z) - findMinZModule().planarMinZ();

  logINFO("Rod length cannot exceed " + any2str(z));
  logINFO("  Z+ rod will be compressed by " + any2str(fabs(Deltap)));
  logINFO("  Z- rod will be compressed by " + any2str(fabs(Deltam)));

  if (allowCompressionCuts()) {
    logINFO("Compression algorithm is allowed to cut modules which fall out entirely from the maximum z line");
    int zPlusOrigSize = zPlusModules_.size();
    for (auto it = zPlusModules_.begin(); it != zPlusModules_.end();) {
      if (it->planarMinZ() > z) it = zPlusModules_.erase(it); 
      else ++it;
    }
    logINFO("  " + any2str(zPlusOrigSize - zPlusModules_.size()) + " modules were cut from the Z+ rod");
    if (zPlusOrigSize - zPlusModules_.size()) { 
      Deltap = fabs(z) - findMaxZModule().planarMaxZ();  // we have to use findMaxZModule instead of checking only the last module as there might be inversions at high Z
      logINFO("  Z+ rod now exceeding by " + any2str(fabs(Deltap)));
    }
    int zMinusOrigSize = zMinusModules_.size();
    for (auto it = zMinusModules_.begin(); it != zMinusModules_.end();) {
      if (it->planarMaxZ() < -z) it = zMinusModules_.erase(it); 
      else ++it;
    }
    logINFO("  " + any2str(zMinusOrigSize - zMinusModules_.size()) + " modules were cut from the Z- rod");
    if (zMinusOrigSize - zMinusModules_.size()) { 
      Deltam = -fabs(z) - findMinZModule().planarMinZ(); 
      logINFO("  Z- rod now exceeding by " + any2str(fabs(Deltam)));
    }
  }

  logINFO("Iterative compression of Z+ rod");
  int i;
  static const int maxIterations = 50;
  for (i = 0; fabs(Deltap) > 0.1 && i < maxIterations; i++) {
    double deltap=Deltap/((findMaxZModule().center()).Z());
    std::map<int, double> zGuards;
    zGuards[1] = zGuards[-1] = std::numeric_limits<double>::lowest();
    int parity = zPlusParity();
    for (auto it = zPlusModules_.begin(); it < zPlusModules_.end(); ++it, parity = -parity) {
      double translation = deltap*it->center().Z();
      double minPhysZ = MIN(it->planarMinZ(), it->center().Z() - it->physicalLength()/2);
      if (minPhysZ + translation < zGuards[parity]) 
        translation = zGuards[parity] - minPhysZ;
      it->translateZ(translation); 
      double maxPhysZ = MAX(it->planarMaxZ(), it->center().Z() + it->physicalLength()/2);
      zGuards[parity] = maxPhysZ;
    }
    Deltap =  fabs(z) - findMaxZModule().planarMaxZ();
  }
  if (i == maxIterations) {
    logINFO("Iterative compression didn't terminate after " + any2str(maxIterations) + " iterations. Z+ rod still exceeds by " + any2str(fabs(Deltap))); 
    logWARNING("Failed to compress Z+ rod. Still exceeding by " + any2str(fabs(Deltap)) + ". Check info tab.");
  } else logINFO("Z+ rod successfully compressed after " + any2str(i) + " iterations. Rod now only exceeds by " + any2str(fabs(Deltap)) + " mm.");

  logINFO("Iterative compression of Z- rod");
  for (i = 0; fabs(Deltam) > 0.1 && i < maxIterations; i++) {
    double deltam=Deltam/((findMinZModule().center()).Z()); // this can be optimized
    std::map<int, double> zGuards;
    zGuards[1] = zGuards[-1] = std::numeric_limits<double>::max();
    int parity = -zPlusParity();
    for (auto it = zMinusModules_.begin(); it < zMinusModules_.end(); ++it, parity = -parity) {
      double translation = deltam*it->center().Z();
      double maxPhysZ = MAX(it->planarMaxZ(), it->center().Z() + it->physicalLength()/2);
      if (maxPhysZ + translation > zGuards[parity]) 
        translation = zGuards[parity] - maxPhysZ;
      it->translateZ(translation); 
      double minPhysZ = MIN(it->planarMinZ(), it->center().Z() - it->physicalLength()/2);
      zGuards[parity] = minPhysZ;
    }
    Deltam = -fabs(z) - findMinZModule().planarMinZ(); 
  }  
  if (i == 20) {
    logINFO("Iterative compression didn't terminate after " + any2str(maxIterations) + " iterations. Z- rod still exceeds by " + any2str(fabs(Deltam))); 
    logWARNING("Failed to compress Z- rod. Still exceeding by " + any2str(fabs(Deltam)) + ". Check info tab.");
  } else logINFO("Z- rod successfully compressed after " + any2str(i) + " iterations. Rod now only exceeds by " + any2str(fabs(Deltam)) + " mm.");

}

std::set<int> StraightRodPair::solveCollisionsZPlus() {
  std::map<int, double> zGuards;
  std::set<int> collisions;
  zGuards[1] = zGuards[-1] = std::numeric_limits<double>::lowest();
  int parity = zPlusParity();
  int n = 1;
  bool found = false;
  logINFO("Checking for collisions in Z+ rod");
  for (auto& m : zPlusModules_) {
    double minPhysZ = m.center().Z() - MAX(m.physicalLength(), m.length())/2;
    if (zGuards[parity] > minPhysZ) {
      double offset = zGuards[parity] - minPhysZ;
      m.translateZ(offset);
      collisions.insert(n);
      logINFO("  Module " + any2str(n) + " collides with previous " + (parity > 0 ? "outer" : "inner") + " module. Translated by " + any2str(offset) + " mm");
      found = true;
    }
    zGuards[parity] = m.center().Z() + MAX(m.physicalLength(), m.length())/2;
    parity = -parity;
    n++;
  }
  return collisions;
}

std::set<int> StraightRodPair::solveCollisionsZMinus() {
  std::map<int, double> zGuards;
  std::set<int> collisions;
  zGuards[1] = zGuards[-1] = std::numeric_limits<double>::max();
  int parity = -zPlusParity();
  int n = -1;
  logINFO("Checking for collisions in Z- rod");
  for (auto& m : zMinusModules_) {
    double maxPhysZ = m.center().Z() + MAX(m.physicalLength(), m.length())/2;
    if (zGuards[parity] < maxPhysZ) {
      double offset = zGuards[parity] - maxPhysZ;
      m.translateZ(offset);
      collisions.insert(n);
      logINFO("  Module " + any2str(n) + " collides with previous " + (parity > 0 ? "outer" : "inner") + " module. Translated by " + any2str(offset) + " mm");
    }
    zGuards[parity] = m.center().Z() - MAX(m.physicalLength(), m.length())/2;
    parity = -parity;
    n--;
  }
  return collisions;
}


double StraightRodPair::computeNextZ(double newDsLength, double newDsDistance, double lastDsDistance, double lastZ, BuildDir direction, int parity) {
  double d = smallDelta();
  double dz = zError();
  double ov = zOverlap();
  double maxr = maxBuildRadius();
  double minr = minBuildRadius();
 
  // Case A : zOverlap is considered for computing next Z
  double lastRA;
  double newRA;
  if (parity > 0) { // Going upwards
    lastRA = maxr - d + lastDsDistance/2;
    newRA = maxr + d - newDsDistance/2; // use of maxr because outer rod wins, from -d to +d because going upwards
  } else { // Going downwards
    lastRA = minr + d + lastDsDistance/2;
    newRA = minr - d - newDsDistance/2; // use of minr because inner rod wins, from +d to -d because going downwards
  }

  // Case B : zError is considered for computing next Z
  double lastRB;
  double newRB; 
  if (parity > 0) { // Going upwards
    if ((direction == BuildDir::RIGHT && lastZ > dz) || (direction == BuildDir::LEFT && lastZ < -dz)) { // lastZ is not in [-zError; zError]
      lastRB = maxr - d + lastDsDistance/2;
      newRB = maxr + d - newDsDistance/2;  // use of maxr because outer rod wins, from -d to +d because going upwards
    }
    else { // lastZ is in [-zError; zError] : special case which is often forgotten
      lastRB = minr - d - lastDsDistance/2;
      newRB = minr + d + newDsDistance/2; // use of minr because inner rod wins, from -d to +d because going upwards
    }
  }
  else { // Going downwards
    lastRB = minr + d + lastDsDistance/2;
    newRB = minr - d - newDsDistance/2; // use of minr because inner rod wins, from +d to -d because going downwards
  }

  double newZ = lastZ;
  if (!beamSpotCover()) dz = 0;
  if (direction == BuildDir::RIGHT) {
    double newZA = (newZ - ov) * newRA/lastRA; // newZ associated to case A
    double originZ = parity > 0 ? dz : -dz;
    double newZB = (newZ - originZ) * newRB/lastRB + originZ; // newZ associated to case B
    if (beamSpotCover()) newZ = MIN(newZA, newZB); // Take the most stringent of cases A and B
    else newZ = newZA;
    if (forbiddenRange.state()) {
      double forbiddenRange_begin,forbiddenRange_end; 
      forbiddenRange_begin=(forbiddenRange[0]+forbiddenRange[1])/2;
      forbiddenRange_end=forbiddenRange[1];
      if (newZ-lastZ >= (forbiddenRange_begin - newDsLength) && newZ - lastZ <= (forbiddenRange_end - newDsLength)) {
        newZ = lastZ + forbiddenRange_begin - newDsLength;
      } else {
         forbiddenRange_begin=forbiddenRange[0];
         forbiddenRange_end=(forbiddenRange[0]+forbiddenRange[1])/2;
	 if (newZ-lastZ >= (forbiddenRange_begin - newDsLength) && newZ - lastZ <= (forbiddenRange_end - newDsLength)) {
            newZ = lastZ + forbiddenRange_begin - newDsLength;
	 }
      }
    }
  } 
  else {
    double newZA = (newZ + ov) * newRA/lastRA; // newZ associated to case A
    double originZ = parity > 0 ? -dz : dz;
    double newZB = (newZ - originZ) * newRB/lastRB + originZ; // newZ associated to case B
    if (beamSpotCover()) newZ = MAX(newZA, newZB); // Take the most stringent of cases A and B
    else newZ = newZA;
    if (forbiddenRange.state()) {
      double forbiddenRange_begin,forbiddenRange_end;              
      forbiddenRange_begin=(forbiddenRange[0]+forbiddenRange[1])/2;
      forbiddenRange_end=forbiddenRange[1];
      if (lastZ - newZ >= (forbiddenRange_begin - newDsLength) && lastZ - newZ <= (forbiddenRange_end - newDsLength)){
        newZ = lastZ - forbiddenRange_begin + newDsLength;
      } else {
        forbiddenRange_begin=forbiddenRange[0];
        forbiddenRange_end=(forbiddenRange[0]+forbiddenRange[1])/2;
        if (lastZ - newZ >= (forbiddenRange_begin - newDsLength) && lastZ - newZ <= (forbiddenRange_end - newDsLength)){
          newZ = lastZ - forbiddenRange_begin + newDsLength;
        }
      }
    }
  }
  return newZ;
}



template<typename Iterator> vector<double> StraightRodPair::computeZList(Iterator begin, Iterator end, double startZ, BuildDir direction, int smallParity, bool fixedStartZ) {

  vector<double> zList;

  double targetZ = maxZ.state() ? maxZ() : std::numeric_limits<double>::max(); // unreachable target in case maxZ not set
  int targetMods = buildNumModules.state() ? buildNumModules() : std::numeric_limits<int>::max(); // unreachable target in case numModules not set
  // If we rely on numModules and we start from the center, then the number of modules on the left side of the rod must be lower
  if (!mezzanine() && (startZMode() == StartZMode::MODULECENTER) && (direction == BuildDir::LEFT)) targetMods--; 

  double newZ = startZ; // + lengthOffset/2;

  int parity = smallParity;
  BarrelModule* lastm = begin->get();

  int n = 0;

  if (fixedStartZ) {
    zList.push_back(newZ);
    newZ += (direction == BuildDir::RIGHT ? (*begin)->length() : (*begin)->length());
    parity = -parity;
    ++begin;
    n = 1;
  }

  for (auto& curm : pair2range(make_pair(begin, end))) {
    if (abs(newZ) >= targetZ || n++ >= targetMods) break;

    newZ = computeNextZ(curm->length(), curm->dsDistance(), lastm->dsDistance(), newZ, direction, parity);
    zList.push_back(newZ);
    newZ += (direction == BuildDir::RIGHT ? curm->length() : -curm->length());
    lastm = curm.get();
    parity = -parity;
  } 

  double lengthOffset = lastm->length();
  double physicalLengthOffset = lastm->physicalLength();

  for (; abs(newZ) < targetZ && n < targetMods; n++) {  // in case the rodtemplate vector is finished but we haven't hit the targetZ or targetMods yet, we keep on using the last module for dsDistance and length
    newZ = computeNextZ(lastm->length(), lastm->dsDistance(), lastm->dsDistance(), newZ, direction, parity);
    zList.push_back(newZ);
    newZ += (direction == BuildDir::RIGHT ? lastm->length() : -lastm->length());
    parity = -parity;
  }

  return zList;

}


template<typename Iterator> pair<vector<double>, vector<double>> StraightRodPair::computeZListPair(Iterator begin, Iterator end, double startZ, int recursionCounter) {
  bool fixedStartZ = true;
  vector<double> zPlusList = computeZList(begin, end, startZ, BuildDir::RIGHT, zPlusParity(), fixedStartZ);
  vector<double> zMinusList = computeZList(begin, end, startZ, BuildDir::LEFT, -zPlusParity(), !fixedStartZ);

  double zUnbalance = 0.;
  if (!zPlusList.empty() && !zMinusList.empty()) { zUnbalance = (zPlusList.back()+(*(end-1))->length()/2) + (zMinusList.back()-(*(end-1))->length()/2); } // balancing uneven pos/neg strings

  if (++recursionCounter == 100) { // this stops infinite recursion if the balancing doesn't converge
    std::ostringstream tempSS;
    tempSS << "Balanced module placement in rod pair at avg build radius " << (maxBuildRadius()+minBuildRadius())/2. << " didn't converge!! Layer is skewed";
    tempSS << "Unbalance is " << zUnbalance << " mm";
    logWARNING(tempSS);

    return std::make_pair(zPlusList, zMinusList);
  }

  if (fabs(zUnbalance) > 0.1) { // 0.1 mm unbalance is tolerated
    return computeZListPair(begin, end,
                            startZ-zUnbalance/2, // countering the unbalance by displacing the startZ (by half the inverse unbalance, to improve convergence)
                            recursionCounter);
  } else {
    std::ostringstream tempSS;
    tempSS << "Balanced module placement in rod pair at avg build radius " << (maxBuildRadius()+minBuildRadius())/2. << " converged after " << recursionCounter << " step(s).\n" 
           << "   Residual Z unbalance is " << zUnbalance << ".\n"
           << "   Positive string has " << zPlusList.size() << " modules, negative string has " << zMinusList.size() << " modules.\n";
    if (!zPlusList.empty()) { tempSS << "Z+ rod starts at " << zPlusList.front() << ".\n"; }
    if (!zMinusList.empty()) { tempSS << "Z- rod starts at " << zMinusList.front() << ".\n"; }
    logINFO(tempSS);
    return std::make_pair(zPlusList, zMinusList);
  }
}

void StraightRodPair::buildModules(Container& modules, const RodTemplate& rodTemplate, const vector<double>& posList, BuildDir direction, bool isPlusBigDeltaRod, int parity, int side) {
  for (int i=0; i<(int)posList.size(); i++, parity = -parity) {
    BarrelModule* mod;
    if (!mezzanine() && (startZMode() == StartZMode::MODULECENTER) && (direction == BuildDir::LEFT)) { // skips the central module information.
      mod = GeometryFactory::make<BarrelModule>(i < (rodTemplate.size() - 1) ? *rodTemplate[i+1].get() : *rodTemplate.rbegin()->get());
      mod->myid(i+2);
    }
    else {
      mod = GeometryFactory::make<BarrelModule>(i < rodTemplate.size() ? *rodTemplate[i].get() : *rodTemplate.rbegin()->get());
      mod->myid(i+1);
    }
    mod->side(side);
    //mod->store(propertyTree());
    //if (ringNode.count(i+1) > 0) mod->store(ringNode.at(i+1)); 
    //mod->build();
    mod->translateR(parity > 0 ? smallDelta() : -smallDelta());
    if (smallDelta() != 0) { mod->flipped(parity != 1); } // When smallDelta() != 0, the flip is alternated.
    else { mod->flipped(!isPlusBigDeltaRod); } // When smallDelta() == 0, the flip only depends whether the rod is located at + BigDelta or at - BigDelta.
    mod->translateZ(posList[i] + (direction == BuildDir::RIGHT ? mod->length()/2 : -mod->length()/2));
    // mod->translate(XYZVector(parity > 0 ? smallDelta() : -smallDelta(), 0, posList[i])); // CUIDADO: we are now translating the center instead of an edge as before
    modules.push_back(mod);
  }
}

void StraightRodPair::buildFull(const RodTemplate& rodTemplate, bool isPlusBigDeltaRod) {
  double startZ = startZMode() == StartZMode::MODULECENTER ? -(*rodTemplate.begin())->length()/2. : 0.;
  /*std::cout << "startZ = " << startZ << std::endl;
  std::cout << "rodTemplate.size() = " << rodTemplate.size() << std::endl;*/
  auto zListPair = computeZListPair(rodTemplate.begin(), rodTemplate.end(), startZ, 0);

    // actual module creation
    // CUIDADO log rod balancing effort
  buildModules(zPlusModules_, rodTemplate, zListPair.first, BuildDir::RIGHT, isPlusBigDeltaRod, zPlusParity(), 1);
  double currMaxZ = zPlusModules_.size() > 1 ? MAX(zPlusModules_.rbegin()->planarMaxZ(), (zPlusModules_.rbegin()+1)->planarMaxZ()) : (!zPlusModules_.empty() ? zPlusModules_.rbegin()->planarMaxZ() : 0.); 
  // CUIDADO this only checks the positive side... the negative side might actually have a higher fabs(maxZ) if the barrel is long enough and there's an inversion
  buildModules(zMinusModules_, rodTemplate, zListPair.second, BuildDir::LEFT, isPlusBigDeltaRod, -zPlusParity(), -1);

  auto collisionsZPlus = solveCollisionsZPlus();
  auto collisionsZMinus = solveCollisionsZMinus();
  if (!collisionsZPlus.empty() || !collisionsZMinus.empty()) logWARNING("Some modules have been translated to avoid collisions. Check info tab");

  if (compressed() && maxZ.state() && currMaxZ > maxZ()) compressToZ(maxZ());
  currMaxZ = zPlusModules_.size() > 1 ? MAX(zPlusModules_.rbegin()->planarMaxZ(), (zPlusModules_.rbegin()+1)->planarMaxZ()) : (!zPlusModules_.empty() ? zPlusModules_.rbegin()->planarMaxZ() : 0.);
  maxZ(currMaxZ);
}

void StraightRodPair::buildMezzanine(const RodTemplate& rodTemplate, bool isPlusBigDeltaRod) {
  // compute Z list (only once since the second mezzanine has just inverted signs for z) 
  vector<double> zList = computeZList(rodTemplate.rbegin(), rodTemplate.rend(), startZ(), BuildDir::LEFT, zPlusParity(), false);
  vector<double> zListNeg;
  std::transform(zList.begin(), zList.end(), std::back_inserter(zListNeg), std::negate<double>());

  buildModules(zPlusModules_, rodTemplate, zList, BuildDir::LEFT, isPlusBigDeltaRod, zPlusParity(), 1); // CUIDADO mezzanine layer rings in reverse order????
  maxZ(startZ());

  buildModules(zMinusModules_, rodTemplate, zListNeg, BuildDir::RIGHT, isPlusBigDeltaRod, zPlusParity(), -1);
}

void StraightRodPair::check() {
  PropertyObject::check();

  if (mezzanine()) {
    if (startZMode.state()) logWARNING("Ignoring startZMode (set to modulecenter by default) for a mezzanine.");
  }
  if (!isFlatPart() && !zPlusParity.state()) throw PathfulException("Straight rod : smallParity must be specified.");
}

void StraightRodPair::build(const RodTemplate& rodTemplate, bool isPlusBigDeltaRod) {
  materialObject_.store(propertyTree());
  materialObject_.build();

  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();
    if (!mezzanine()) buildFull(rodTemplate, isPlusBigDeltaRod);
    else buildMezzanine(rodTemplate, isPlusBigDeltaRod);
  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  cleanup();
  builtok(true);
}

void TiltedRodPair::buildModules(Container& modules, const RodTemplate& rodTemplate, const vector<TiltedModuleSpecs>& tmspecs, BuildDir direction, bool flip) {
  auto it = rodTemplate.begin();
  int side = (direction == BuildDir::RIGHT ? 1 : -1);
  if (tmspecs.empty()) return;
  int i = 0;
  if (direction == BuildDir::LEFT && fabs(tmspecs[0].z) < 0.5) { i = 1; it++; } // this skips the first module if we're going left (i.e. neg rod) and z=0 because it means the pos rod has already got a module there
  for (; i < tmspecs.size(); i++, ++it) {
    //std::cout << "i = " << i << std::endl;
    //std::cout << "tmspecs[i].r = " << tmspecs[i].r << std::endl;
    //std::cout << "tmspecs[i].z = " << tmspecs[i].z << std::endl;
    BarrelModule* mod = GeometryFactory::make<BarrelModule>(**it);
    mod->myid(i+1);
    mod->side(side);
    mod->tilt(side * tmspecs[i].gamma);
    mod->translateR(tmspecs[i].r);
    if (tmspecs[i].gamma == 0) { mod->flipped(i%2); } // flat part of the tilted rod, i is the ring number
    else { mod->flipped(flip); } // tilted part of the tilted rod
    mod->translateZ(side * tmspecs[i].z);
    modules.push_back(mod);
  }
}

void TiltedRodPair::check() {
  PropertyObject::check();

  if (startZMode() != StartZMode::MODULECENTER) throw PathfulException("Tilted layer : only startZMode = modulecenter can be specified.");
}


void TiltedRodPair::build(const RodTemplate& rodTemplate, const std::vector<TiltedModuleSpecs>& tmspecs, bool flip) {
  materialObject_.store(propertyTree());
  materialObject_.build();

  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();
    buildModules(zPlusModules_, rodTemplate, tmspecs, BuildDir::RIGHT, flip);
    buildModules(zMinusModules_, rodTemplate, tmspecs, BuildDir::LEFT, flip);
    //maxZ(zPlusModules_.back().maxZ());
    //if (!zPlusModules_.empty()) {
    //std::cout << "zPlusModules_.back().maxZ() = " << zPlusModules_.back().maxZ() << std::endl;
    //}

  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }
  cleanup();
  builtok(true);
}


define_enum_strings(RodPair::StartZMode) = { "modulecenter", "moduleedge" };

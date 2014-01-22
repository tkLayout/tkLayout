#include "RodPair.h"

void RodPair::clearComputables() { 
  minAperture.clear(); 
  maxAperture.clear(); 
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


void StraightRodPair::compressToZ(double z) {

  if (zPlusModules_.empty()) return;

  const BarrelModule& maxBarrelModule = *std::max_element(zPlusModules_.begin(), zPlusModules_.end(), [](const BarrelModule& m1, const BarrelModule& m2) { return m1.maxZ() < m2.maxZ(); });
  const BarrelModule& minBarrelModule = *std::min_element(zMinusModules_.begin(), zMinusModules_.end(), [](const BarrelModule& m1, const BarrelModule& m2) { return m1.minZ() < m2.minZ(); });

  double maxZ = fabs(maxBarrelModule.maxZ());
  double minZ = fabs(minBarrelModule.minZ());

  double Deltap =  fabs(z) - maxZ;
  double Deltam = -fabs(z) + minZ;

  double deltam=0;
  double deltap=0;
  XYZVector modShift;
  deltam=Deltam/((minBarrelModule.center()).Z());
  deltap=Deltap/((maxBarrelModule.center()).Z());
  for (auto& m : zPlusModules_) {
    double myMeanZ = m.center().Z();
    m.translateZ(myMeanZ > 0 ? deltap*myMeanZ : deltam*myMeanZ);
  }
  for (auto& m : zMinusModules_) {
    double myMeanZ = m.center().Z();
    m.translateZ(myMeanZ > 0 ? deltap*myMeanZ : deltam*myMeanZ);
  }

}


double StraightRodPair::computeNextZ(double newDsDistance, double lastDsDistance, double lastZ, BuildDirection direction, int parity) {
  double d = smallDelta();
  double dz = zError();
  double ov = minModuleOverlap();
  double maxr = maxBuildRadius();
  double minr = minBuildRadius();

  double newR = (parity > 0 ? maxr + d : minr - d) - newDsDistance/2;
  double lastR = (parity > 0 ? maxr - d : minr + d) + lastDsDistance/2;

  double newZ = lastZ;

  if (direction == BuildDirection::RIGHT) {
    double originZ = parity > 0 ? dz : -dz;
    double newZorigin = (newZ - ov) * newR/lastR;
    double newZshifted = (newZ - originZ) * newR/lastR + originZ;
    newZ = MIN(newZorigin, newZshifted);
  } else {
    double originZ = parity > 0 ? -dz : dz;
    double newZorigin = (newZ + ov) * newR/lastR;
    double newZshifted = (newZ - originZ) * newR/lastR + originZ;
    newZ = MAX(newZorigin, newZshifted);
  }

  return newZ;
}



template<typename Iterator> vector<double> StraightRodPair::computeZList(Iterator begin, Iterator end, double startZ, BuildDirection direction, int smallParity, bool looseStartZ) {

  vector<double> zList;

  double targetZ = maxZ.state() ? maxZ() : std::numeric_limits<double>::max(); // unreachable target in case maxZ not set
  int targetMods = buildNumModules.state() ? buildNumModules() : std::numeric_limits<int>::max(); // unreachable target in case numModules not set

  double lengthOffset = direction == BuildDirection::RIGHT ? (*begin)->length() : -(*begin)->length();

  double newZ = startZ; // + lengthOffset/2;

  int parity = smallParity;
  BarrelModule* lastm = begin->get();

  int n = 0;

  if (!looseStartZ) {
    zList.push_back(newZ);
    newZ = newZ + lengthOffset; 
    parity = -parity;
    ++begin;
    n = 1;
  }

  for (auto& curm : pair2range(make_pair(begin, end))) {
    if (abs(newZ) >= targetZ || n++ >= targetMods) break;
    newZ = computeNextZ(curm->dsDistance(), lastm->dsDistance(), newZ, direction, parity);
    zList.push_back(newZ);
    newZ += lengthOffset; 
    lastm = curm.get();
    parity = -parity;
  } 

  for (; abs(newZ) < targetZ && n < targetMods; n++) {  // in case the rodtemplate vector is finished but we haven't hit the targetZ or targetMods yet, we keep on using the last module for dsDistance and length
    newZ = computeNextZ(lastm->dsDistance(), lastm->dsDistance(), newZ, direction, parity);
    zList.push_back(newZ);
    newZ += lengthOffset; 
    parity = -parity;
  }
  //return listZ[listZ.size()-1] + (direction > 0 ? modLengthZ : -modLengthZ);
  // CUIDADO ???

  return zList;

}


template<typename Iterator> pair<vector<double>, vector<double>> StraightRodPair::computeZListPair(Iterator begin, Iterator end, double startZ, int recursionCounter) {

  bool looseStartZ = false;
  vector<double> zPlusList = computeZList(begin, end, startZ, BuildDirection::RIGHT, zPlusParity(), looseStartZ);
  vector<double> zMinusList = computeZList(begin, end, startZ, BuildDirection::LEFT, -zPlusParity(), !looseStartZ);

  double zUnbalance = (zPlusList.back()+(*(end-1))->length()/2) + (zMinusList.back()-(*(end-1))->length()/2); // balancing uneven pos/neg strings
  if (fabs(zUnbalance) > 0.1 && ++recursionCounter < 100) { // 0.1 mm unbalance is tolerated
    return computeZListPair(begin, end,
                            startZ-zUnbalance/2, // countering the unbalance by displacing the startZ (by half the inverse unbalance, to improve convergence)
                            recursionCounter);
  } else {
    // CUIDADO HANDLE RECURSION COUNTER HITTING 100 ERROR
    return std::make_pair(zPlusList, zMinusList);
  }
}

void StraightRodPair::buildModules(Container& modules, const RodTemplate& rodTemplate, const vector<double>& posList, BuildDirection direction, int parity, int side) {
  for (int i=0; i<(int)posList.size(); i++, parity = -parity) {
    BarrelModule* mod = new BarrelModule(i < rodTemplate.size() ? *rodTemplate[i].get() : *rodTemplate.rbegin()->get());
    mod->setup();
    mod->myid(i+1);
    mod->side(side);
    //mod->store(propertyTree());
    //if (ringNode.count(i+1) > 0) mod->store(ringNode.at(i+1)); 
    //mod->build();
    mod->translateR(parity > 0 ? smallDelta() : -smallDelta());
    mod->translateZ(posList[i] + (direction == BuildDirection::RIGHT ? mod->length()/2 : -mod->length()/2));
   // mod->translate(XYZVector(parity > 0 ? smallDelta() : -smallDelta(), 0, posList[i])); // CUIDADO: we are now translating the center instead of an edge as before
    modules.push_back(mod);
  }
}

void StraightRodPair::buildFull(const RodTemplate& rodTemplate) {

  auto zListPair = computeZListPair(rodTemplate.begin(), rodTemplate.end(), 0., 0);

    // actual module creation
    // CUIDADO log rod balancing effort
  buildModules(zPlusModules_, rodTemplate, zListPair.first, BuildDirection::RIGHT, zPlusParity(), 1);
  double currMaxZ = zPlusModules_.size() > 1 ? MAX(zPlusModules_.rbegin()->maxZ(), (zPlusModules_.rbegin()+1)->maxZ()) : (!zPlusModules_.empty() ? zPlusModules_.rbegin()->maxZ() : 0.); 
  // CUIDADO this only checks the positive side... the negative side might actually have a higher fabs(maxZ) if the barrel is long enough and there's an inversion
  buildModules(zMinusModules_, rodTemplate, zListPair.second, BuildDirection::LEFT, -zPlusParity(), -1);

  if (maxZ.state() && currMaxZ > maxZ()) compressToZ(maxZ());
  else maxZ(currMaxZ);
}

void StraightRodPair::buildMezzanine(const RodTemplate& rodTemplate) {
  // compute Z list (only once since the second mezzanine has just inverted signs for z) 
  vector<double> zList = computeZList(rodTemplate.rbegin(), rodTemplate.rend(), startZ(), BuildDirection::LEFT, zPlusParity(), false);
  vector<double> zListNeg;
  std::transform(zList.begin(), zList.end(), std::back_inserter(zListNeg), std::negate<double>());

  buildModules(zPlusModules_, rodTemplate, zList, BuildDirection::LEFT, zPlusParity(), 1); // CUIDADO mezzanine layer rings in reverse order????
  maxZ(startZ());

  buildModules(zMinusModules_, rodTemplate, zListNeg, BuildDirection::RIGHT, zPlusParity(), -1);
}


void StraightRodPair::build(const RodTemplate& rodTemplate) {
  try {
    std::cout << ">>> Building " << fullid(*this) << " <<<" << std::endl;
    check();
    if (!mezzanine()) buildFull(rodTemplate);
    else buildMezzanine(rodTemplate);
  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  cleanup();
  builtok(true);
}

void TiltedRodPair::buildModules(Container& modules, const RodTemplate& rodTemplate, const vector<TiltedModuleSpecs>& tmspecs, BuildDirection direction) {
  auto it = rodTemplate.begin();
  int side = (direction == BuildDirection::RIGHT ? 1 : -1);
  for (int i = 0; i < tmspecs.size(); i++, ++it) {
    BarrelModule* mod = new BarrelModule(**it);
    mod->setup();
    mod->myid(i+1);
    mod->side(side);
    mod->tilt(side * tmspecs[i].gamma);
    mod->translateR(tmspecs[i].r);
    mod->translateZ(side * tmspecs[i].z);
    if (i == 0 
        && !modules.empty() 
        && fabs(mod->center().Z() - modules.front().center().Z()) < 1e-3 
        && fabs(mod->center().Rho() - modules.front().center().Rho()) < 1e-3) {
      delete mod; // we skip the creation of the first module if it's at the same coordinates as the first module of the opposite Z rod in the pair 
    } else {
      modules.push_back(mod);
    }
  }
}


void TiltedRodPair::build(const RodTemplate& rodTemplate, const std::vector<TiltedModuleSpecs>& tmspecs) {
  try {
    std::cout << ">>> Building " << fullid(*this) << " <<<" << std::endl;
    check();
    buildModules(zPlusModules_, rodTemplate, tmspecs, BuildDirection::RIGHT);
    buildModules(zMinusModules_, rodTemplate, tmspecs, BuildDirection::LEFT);

  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }
  cleanup();
  builtok(true);
}

#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "global_funcs.h"
#include "layer.hh"
#include "Math/RotationZ.h"
#include "Math/Vector3D.h"

using namespace ROOT::Math;

/******************/
/*                */
/* Generic module */
/*                */
/******************/

Layer::~Layer() {
  ModuleVector::iterator modIt;
  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    if ((*modIt)!=NULL) {
      delete (*modIt);
    }
  }
  moduleSet_.clear();
}

Layer::Layer() {
  setDefaultParameters();
}

void Layer::setDefaultParameters() {
  layerName_ = "Layer";
}

void Layer::translate(XYZVector Delta) {
  ModuleVector::iterator modIt;

  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    (*modIt)->translate(Delta);
  }
}


// TODO: tidy up the next two functions: the Endcap-specific onw should
// just call the generic Layer::rotateY_PI(), and multiply the averageZ_ by -1
void BarrelLayer::rotateY_PI() {
  ModuleVector::iterator modIt;

  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    (*modIt)->rotateY_PI();
  }
}

void BarrelLayer::reflectZ() {
  ModuleVector::iterator modIt;

  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    (*modIt)->reflectZ();
  }
}

void BarrelLayer::shiftRho(double Delta) {
  ModuleVector::iterator modIt;

  averageRadius_ = InvalidRadius;
  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    (*modIt)->shiftRho(Delta);
  }
}

void EndcapLayer::rotateY_PI() {
  ModuleVector::iterator modIt;

  averageZ_ *= -1;

  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    (*modIt)->rotateY_PI();
  }
}


void Layer::shapeVolume(TGeoVolume* container,
                        TGeoMedium* medium,
                        TGeoManager* geom) {
  // TODO:
  // Should build only the containing volume
}

void Layer::shapeModuleVolumes(TGeoVolume* container,
                               TGeoMedium* medium,
                               TGeoManager* geom) {

  // TODO:
  // Should build only the containing volume

}

double Layer::getMaxZ() {
  double maxZ;

  ModuleVector::iterator modIt;
  maxZ=(*moduleSet_.begin())->getMaxZ();
  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    if ((*modIt)->getMaxZ()>maxZ) maxZ=(*modIt)->getMaxZ();
  }

  return maxZ;
};

double Layer::getMinZ() {
  double minZ;

  ModuleVector::iterator modIt;
  minZ=(*moduleSet_.begin())->getMinZ();
  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    if ((*modIt)->getMinZ()<minZ) minZ=(*modIt)->getMinZ();
  }

  return minZ;
};

double Layer::getMaxRho() {
  double maxRho;

  ModuleVector::iterator modIt;
  maxRho=(*moduleSet_.begin())->getMaxRho();
  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    if ((*modIt)->getMaxRho()>maxRho) maxRho=(*modIt)->getMaxRho();
  }

  return maxRho;
};

double Layer::getMinRho() {
  double minRho;

  ModuleVector::iterator modIt;
  minRho=(*moduleSet_.begin())->getMinRho();
  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    if ((*modIt)->getMinRho()<minRho) minRho=(*modIt)->getMinRho();
  }

  return minRho;
};

int Layer::getNModules() const {
  return moduleSet_.size();
}

/******************/
/*                */
/* Barrel layer   */
/*                */
/******************/

BarrelLayer::~BarrelLayer() {
  if (sampleModule_) delete sampleModule_;
}

BarrelLayer::BarrelLayer() {
  sampleModule_ = NULL;
}

BarrelLayer::BarrelLayer(BarrelLayer& inputLayer) {
  ModuleVector::iterator modIt;
  ModuleVector inputModuleV = inputLayer.moduleSet_;
  BarrelModule* aModule;
  BarrelModule* stdModule;

  setName(inputLayer.layerName_, inputLayer.layerIndex_);
  setContainerName(inputLayer.containerName_);
  averageRadius_ = inputLayer.averageRadius_;
  nOfRods_ = inputLayer.getRods();
  nModsOnString_ = inputLayer.getModulesOnRod();

  if ( (sampleModule_=dynamic_cast<BarrelModule*>(inputLayer.getSampleModule())) ) {
    sampleModule_ = new BarrelModule(*(inputLayer.getSampleModule()));
  }

  for (modIt=inputModuleV.begin(); modIt!=inputModuleV.end(); modIt++) {
    if ( (stdModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
      aModule = new BarrelModule(*stdModule);
      moduleSet_.push_back(aModule);
    } else {
      logERROR("In BarrelLayer::BarrelLayer(BarrelLayer& inputLayer) I found a non-barrel module in the source");
    }
  }
}

BarrelLayer::BarrelLayer(double waferDiameter, double heightOverWidth) {
  sampleModule_ = new BarrelModule(waferDiameter, heightOverWidth);
}

BarrelLayer::BarrelLayer(double heightOverWidth ) {
  sampleModule_ = new BarrelModule(heightOverWidth);
}

BarrelLayer::BarrelLayer(const BarrelModule& mySample) {
  sampleModule_ = new BarrelModule(mySample);
}

BarrelLayer::BarrelLayer(BarrelModule* mySample) {
  sampleModule_ = new BarrelModule(*mySample);
}


double BarrelLayer::getMaxModuleThickness() {
  ModuleVector::iterator iter, guard = moduleSet_.end();
  double res = 0.0;
  for (iter = moduleSet_.begin(); iter != guard; iter++) {
    if ((*iter)->getModuleThickness() > res) res = (*iter)->getModuleThickness();
  }
  return res;
}



double BarrelLayer::PlaceWithMaxZ::operator()( // MAX Z VERSION
  std::vector<double>& listZ,
  double startZ,
  std::pair<double, double> worstCaseRadii,
  double bigDelta,
  double smallDelta,
  const vector<double>& dsDistances, // sensor spacings
  double modLengthZ,
  double originDeltaZ,
  double baseOverlapZ,
  int parity,
  int direction,
  bool looseStartZ) const {

  double newZ = startZ;

  int i = 0;
  if (!looseStartZ) {
    listZ.push_back(startZ);
    newZ = startZ + (direction > 0 ? modLengthZ : -modLengthZ);
    i++;
    parity = -parity;
  }

  while (abs(newZ) < maxZ_) { // maxZ is always > 0 even for Neg Z rods
    if (i >= (int)dsDistances.size()) { logERROR("computeListZ(): trying to place more modules than the dsDistances vector specifies. Rod will be shorter"); break; }
    double newR = (parity > 0 ? worstCaseRadii.second + bigDelta + smallDelta : worstCaseRadii.first - bigDelta - smallDelta) - dsDistances[i]/2; 
    double lastR = (parity > 0 ? worstCaseRadii.second + bigDelta - smallDelta : worstCaseRadii.first - bigDelta + smallDelta) + dsDistances[i > 0 ? i-1 : 0]/2;// (if looseStart) the previous module (opposite rod) has the same dsDistance
    if (direction > 0) {
      double originZ = parity > 0 ? originDeltaZ : -originDeltaZ;
      double newZorigin = (newZ - baseOverlapZ)*newR/lastR;
      double newZshifted = (newZ - originZ)*newR/lastR + originZ;
      newZ = newZorigin < newZshifted ? newZorigin : newZshifted;
      listZ.push_back(newZ);
      newZ += modLengthZ;
    } else {
      double originZ = parity > 0 ? -originDeltaZ : originDeltaZ;
      double newZorigin = (newZ + baseOverlapZ)*newR/lastR;
      double newZshifted = (newZ - originZ)*newR/lastR + originZ;
      newZ = newZorigin > newZshifted ? newZorigin : newZshifted;
      listZ.push_back(newZ);
      newZ -= modLengthZ;
    }
    parity = -parity;
    i++;
  }

  /*double excess = abs(newZ) - maxZ_;
  if (excess > 0) { // rod is too long, we have to compress it
    for (i = 0; i < (int)listZ.size(); i++) {
      listZ[i] -= (direction > 0 ? excess : -excess) / listZ.size() * i;
    }
  }*/

  //return listZ[listZ.size()-1] + (direction > 0 ? modLengthZ : -modLengthZ);
  //
  return newZ;
}



double BarrelLayer::PlaceWithNumModules::operator()( //NUM MODULES VERSION
  std::vector<double>& listZ,
  double startZ,
  std::pair<double, double> worstCaseRadii,
  double bigDelta,
  double smallDelta,
  const vector<double>& dsDistances, // sensor spacings
  double modLengthZ,
  double originDeltaZ,
  double baseOverlapZ,
  int parity,
  int direction,
  bool looseStartZ) const {

  double newZ = startZ;

  int i = 0;
  if (!looseStartZ) {
    listZ.push_back(startZ);
    newZ = startZ + (direction > 0 ? modLengthZ : -modLengthZ);
    i++;
    parity = -parity;
  }

  for (; i<numModules_; i++) {
    double newR = (parity > 0 ? worstCaseRadii.second + bigDelta + smallDelta : worstCaseRadii.first - bigDelta - smallDelta) - dsDistances[i]/2; 
    double lastR = (parity > 0 ? worstCaseRadii.second + bigDelta - smallDelta : worstCaseRadii.first - bigDelta + smallDelta) + dsDistances[i > 0 ? i-1 : 0]/2;// (if looseStart) the previous module (opposite rod) has the same dsDistance
    if (direction > 0) {
      double originZ = parity > 0 ? originDeltaZ : -originDeltaZ;
      double newZorigin = (newZ - baseOverlapZ)*newR/lastR;
      double newZshifted = (newZ - originZ)*newR/lastR + originZ;
      newZ = newZorigin < newZshifted ? newZorigin : newZshifted;
      listZ.push_back(newZ);
      newZ += modLengthZ;
    } else {
      double originZ = parity > 0 ? -originDeltaZ : originDeltaZ;
      double newZorigin = (newZ + baseOverlapZ)*newR/lastR;
      double newZshifted = (newZ - originZ)*newR/lastR + originZ;
      newZ = newZorigin > newZshifted ? newZorigin : newZshifted;
      listZ.push_back(newZ);
      newZ -= modLengthZ;
    }
    parity = -parity;
  }

  return newZ;
}

void BarrelLayer::buildStringPair(
  ModuleVector& thisModuleSet,
  double averageRadius,
  std::pair<double, double> worstCaseRadii,
  double bigDelta,
  double smallDelta,
  const vector<double>& dsDistances,
  double baseOverlap,
  double zDelta,
  const ModulePlacementStrategy& computeListZ,
  int smallParity,
  BarrelModule* sampleModule) {

  buildStringPairRecursion(thisModuleSet,
                           averageRadius,
                           worstCaseRadii,
                           bigDelta, 
                           smallDelta, 
                           dsDistances, 
                           baseOverlap,
                           zDelta, 
                           0, // startZ
                           computeListZ,
                           smallParity,
                           0, // recursion counter
                           sampleModule);
} 

void BarrelLayer::buildStringPairRecursion(
  ModuleVector& thisModuleSet,
  double averageRadius,
  std::pair<double, double> worstCaseRadii,
  double bigDelta,
  double smallDelta,
  const vector<double>& dsDistances,
  double baseOverlap,
  double zDelta,
  double startZ,
  const ModulePlacementStrategy& computeListZ,
  int smallParity,
  int recursionCounter,
  BarrelModule* sampleModule) {

  if (recursionCounter++ == 100) { // this stops infinite recursion if the balancing doesn't converge
    tempString.str("");
    tempString << "Balanced module placement in string at avg radius " << averageRadius << " didn't converge!! String badly skewed";
    tempString << "Unbalance is " << startZ << " mm";
    logWARNING(tempString.str());
    return;
  }  
  // create Z lists and balance them

  // pre-compute parameters
  double modLengthZ = sampleModule->getMaxZ() - sampleModule->getMinZ();
  //double radiiRatio = (averageRadius+smallDelta)/(averageRadius-smallDelta);

  // compute Z list for positive strings
  std::vector<double> listZPos;
  //double farthestPosZ = computeListZ(listZPos, startZ, averageRadius-smallDelta, averageRadius+smallDelta, modLengthZ, DSDISTANCE, zDelta, baseOverlap, numModules, smallParity, 1, false);
  double farthestPosZ = computeListZ(listZPos, startZ, worstCaseRadii, bigDelta, smallDelta, dsDistances, modLengthZ, zDelta, baseOverlap, smallParity, 1, false);
  
  // compute Z list for negative strings
  std::vector<double> listZNeg;
  double farthestNegZ = computeListZ(listZNeg, startZ, worstCaseRadii, bigDelta, smallDelta, dsDistances, modLengthZ, zDelta, baseOverlap, -smallParity, -1, true);

  double zUnbalance = farthestPosZ + farthestNegZ; // balancing uneven pos/neg strings
  if (abs(zUnbalance) > 0.1) { // 0.1 mm unbalance is tolerated
    buildStringPairRecursion(thisModuleSet,
                             averageRadius,
                             worstCaseRadii,
                             bigDelta,
                             smallDelta,
                             dsDistances,
                             baseOverlap,
                             zDelta,
                             startZ-zUnbalance/2, // countering the unbalance by displacing the startZ (by half the inverse unbalance, to improve convergence)
                             computeListZ,
                             smallParity,
                             recursionCounter,
                             sampleModule);
  } else {

    ostringstream tempSS;
    // actual module creation
    tempSS << "Balanced module placement in string at avg radius "
      << averageRadius << " converged after "
      << recursionCounter << " step(s)."
      << " Residual Z unbalance is " << zUnbalance << "."
      << " Positive string has " << listZPos.size() << " modules, negative string has " << listZNeg.size() << " modules.";
    logINFO(tempSS);
    recursionCounter = 0; // we made it so we reset the recursion counter for the next string

    // positive Z string
    for (int i=0, parity = smallParity; i<(int)listZPos.size(); i++, parity = -parity) {
      BarrelModule* m = new BarrelModule(*sampleModule);
      m->translate(XYZVector(0, (parity > 0 ? smallDelta : -smallDelta) + averageRadius, 0));
      m->setEdgeZ(listZPos[i], 1);
      m->setRing(i+1);
      m->setZSide(1);
      thisModuleSet.push_back(m);        
    }

    // negative Z string
    for (int i=0, parity = -smallParity; i<(int)listZNeg.size(); i++, parity = -parity) {
      BarrelModule* m = new BarrelModule(*sampleModule);
      m->translate(XYZVector(0, (parity > 0 ? smallDelta : -smallDelta) + averageRadius, 0));
      m->setEdgeZ(listZNeg[i], -1);
      m->setRing(i+1);
      m->setZSide(-1);
      thisModuleSet.push_back(m);        
    }
  }
}


void BarrelLayer::buildMezzanineStringPair(
  ModuleVector& thisModuleSet,
  double averageRadius,
  std::pair<double, double> worstCaseRadii,
  double bigDelta,
  double smallDelta,
  const vector<double>& dsDistances,
  double baseOverlap,
  double zDelta,
  double startZ, // the first mezzanine string of the pair will start at startZ with direction, the second at -startZ, both with a direction of -sign(startZ)
  const ModulePlacementStrategy& computeListZ,
  int smallParity,
  BarrelModule* sampleModule) {

  // pre-compute parameters
  double modLengthZ = sampleModule->getMaxZ() - sampleModule->getMinZ();
  //double radiiRatio = (averageRadius+smallDelta)/(averageRadius-smallDelta);

  vector<double> reversedDsDistances(dsDistances.size());
  std::reverse_copy(dsDistances.begin(), dsDistances.end(), reversedDsDistances.begin()); // Mezzanine modules are placed from farthest to closest, therefore the dsDistances must be reversed as they are ordered from closest to farthest

  int direction = -((startZ>0) || (startZ<0));

  // compute Z list (only once since the second mezzanine has just inverted signs for z) 
  std::vector<double> listZ;
  computeListZ(listZ, startZ, worstCaseRadii, bigDelta, smallDelta, reversedDsDistances, modLengthZ, zDelta, baseOverlap, smallParity, direction, false);

  for (int i=0, parity = smallParity; i<(int)listZ.size(); i++, parity = -parity) {
    BarrelModule* m = new BarrelModule(*sampleModule);
    m->translate(XYZVector(0, (parity > 0 ? smallDelta : -smallDelta) + averageRadius, 0));
    m->setEdgeZ(listZ[i], -1);
    m->setRing(i+1);
    m->setZSide(1);
    thisModuleSet.push_back(m);        
  }

  for (int i=0, parity = smallParity; i<(int)listZ.size(); i++, parity = -parity) {
    BarrelModule* m = new BarrelModule(*sampleModule);
    m->translate(XYZVector(0, (parity > 0 ? smallDelta : -smallDelta) + averageRadius, 0));
    m->setEdgeZ(-listZ[i], 1);
    m->setRing(i+1);
    m->setZSide(-1);
    thisModuleSet.push_back(m);        
  }
}



// Computes the optimal radius of the inner of two
// layers, which should cover each other
double BarrelLayer::layerRadius(const double& nMod,  // number of strings in a layer
                                const double& g,     // Gap between inner and outer
                                const double& o,     // Needed overlap in mm
                                const double& b      // Module's half width
                               ) {
  double result;

  double f = b - o;

  double t = tan(2*M_PI/nMod);
  double B = t*g - f - b;
  double C = -1 * (t*f*b+f*g);

  // The solution is always the one with plus
  result = -1*B+sqrt(pow(B, 2)-4*t*C);
  result /= 2*t;

  return result;
}



// The following formula is used:
//     alpha + beta = gamma
//     
//     where: alpha is the outer rod top sensor coverage angle (unknown, depends on the radius)
//            beta  is the inner rod top sensor coverage angle (unknown, depends on the radius)
//            gamma is the target coverage angle (known given the number of modules)
//
//              f                 f         2 * PI
//     atan( ------- ) + atan( ------- ) =  ------
//            r + R             r + S       nmods
//
//     where: r is the average radius (unknown)
//            f is the half length of a module minus the overlap (known) with the next one.
//            r + R is the radius of the top sensor of the outer rod, expressed as distance (known) from the average radius (unknown)
//            r + S is the radius of the top sensor of the inner rod, expressed as distance (known) from the average radius (unknown)
//
//      the above formula results in a quadratic equation to solve for r. Only the positive solution is taken (negative radius impossible)
//
double BarrelLayer::layerRadius(int numMods,
                                double bigD,
                                double smallD,
                                double dsDist,
                                double modW,
                                double overlap) {

  double d = dsDist/2;

  double f = (modW/2) - (overlap/2);

  double R = bigD + smallD + d;
  double S = -bigD + smallD + d;

  double T = tan(2*M_PI/numMods);

  double a = T;
  double b = R*T + S*T - 2*f;
  double c = R*S*T - R*f - S*f - T*f*f;

  double r = (-b + sqrt(b*b - 4*a*c))/(2*a);

  return r;
}


std::pair<double, int> BarrelLayer::computeRadius(double avgR,
                                                  double bigD,
                                                  double smallD,
                                                  double dsDist,
                                                  double modW,
                                                  double overlap,
                                                  int optimal,
                                                  int sliceMods) { // number of mods forming a slice. the ring will be composed of a certain number of slices each with sliceMods modules, so that the total number of mods will be multiple of this number

  std::pair<double, int> result;

  double d = dsDist/2;
  double f = (modW/2) - (overlap/2);

  double gamma = atan(f/(avgR + bigD + smallD + d)) + atan(f/(avgR - bigD + smallD + d));

  double slices = 2*M_PI/(gamma * sliceMods); // number of slices

  switch (optimal) {
  case SHRINK:
    slices = floor(slices);
    result.first = layerRadius(slices*sliceMods, bigD, smallD, dsDist, modW, overlap);
    break;
  case ENLARGE:
    slices = ceil(slices);
    result.first = layerRadius(slices*sliceMods, bigD, smallD, dsDist, modW, overlap);
    break;
  case FIXED:
    result.first= avgR;
    slices = ceil(slices);
    break;
  case AUTO:
    double nLo, nHi;
    double rLo, rHi;

    nLo = floor(slices);
    nHi = ceil(slices);
    rLo = layerRadius(nLo*sliceMods, bigD, smallD, dsDist, modW, overlap);
    rHi = layerRadius(nHi*sliceMods, bigD, smallD, dsDist, modW, overlap);

    if (fabs(rHi-avgR)<fabs(rLo-avgR)) {
      result.first = rHi;
      slices = nHi;
    } else {
      result.first = rLo;
      slices = nLo;
    }
    break;
  }
  result.second = slices * sliceMods;

  return result;
}



// Computes the optimal radius of the inner of two
// layers, which should cover each other, using the desired
// optimization
std::pair<double, int> BarrelLayer::computeRadius(const double& x,    // radius of the inner module
                                                  const double& g,    // Gap between inner and outer
                                                  const double& o,    // Needed overlap in mm
                                                  const double& b,    // Module's half width
                                                  const int& optimal, // wether to shrink or enlarge, fix or auto
                                                  const int& base )   // Base number of modules
{

  std::pair<double, int> result;

  // Module's half width minus the desired overlap
  double f = b - o;

  // (modules=2*nHalf)
  //  double nHalf = 2*M_PI/(2*atan((f*(x+g)+b*x)/(((x+g)*x)-f*b)));

  // Tentative number of modules on a base
  // (2-> half shell)
  // (4-> quarter shell)
  double nBase = 2*M_PI/(base*atan((f*(x+g)+b*x)/(((x+g)*x)-f*b)));

  // Here the behaviour depends on the pushDirection
  switch (optimal) {
  case SHRINK:
    nBase = floor(nBase);
    result.first = layerRadius(base*nBase, g, o, b);
    break;
  case ENLARGE:
    nBase = ceil(nBase);
    result.first = layerRadius(base*nBase, g, o, b);
    break;
  case FIXED:
    result.first=x;
    nBase = ceil(nBase);
    break;
  case AUTO:
    double nLo, nHi;
    double rLo, rHi;

    nLo = floor(nBase);
    nHi = ceil(nBase);
    rLo = layerRadius(base*nLo, g, o, b);
    rHi = layerRadius(base*nHi, g, o, b);

    if (fabs(rHi-x)<fabs(rLo-x)) {
      result.first = rHi;
      nBase = nHi;
    } else {
      result.first = rLo;
      nBase = nLo;
    }
    break;
  }

  result.second=int(nBase)*base;

  return result;

}


// Returns the optimal average radius
// along with the total number of strings in a layer
// It only requires the number of strings in a layer to be
// multiple of 2
std::pair<double, int> BarrelLayer::layerPhi(double avRad,          // Desired layer average radius
                                             double bigDelta,       // String half gap
                                             double smallDelta,     // Half distance between modules in the same string
                                             double maxDsDistance,  // Maximum sensor spacing in a rod
                                             double o,              // overlap
                                             double b,              // Module width
                                             int optim,             // wether to shrink or enlarge, fix or auto
                                             int base,              // base number of modules (nMod=base*m)
                                             bool stringSameParity  // Wether the strings should have
                                             // the same or opposite module parity (in/out)
                                            ) {

  // I chose OppdInner for opposite-parity strings because it's the worse condition
  // (of the smartest configuration).

  std::pair<double, int> result;

  if (stringSameParity) {
    // SameOuter; (it's always like this)
    //result = computeRadius(avRad-bigDelta+smallDelta, 2*bigDelta, o, b, optim, base);
    //result.first -= -bigDelta+smallDelta;

    result = computeRadius(avRad, bigDelta, smallDelta, maxDsDistance, b, o, optim, base);  


    // SameInner;
    // result = computeRadius(avRad-bigDelta-smallDelta, 2*bigDelta, o, b, optim, base);
    // result.first -= bigDelta+smallDelta;
  } else {
    // Here the worst condition depends on the parameters
    // OppdInner;
    //result = computeRadius(avRad-bigDelta+smallDelta, 2*(bigDelta-smallDelta), o, b, optim, base);
    //result.first -= -bigDelta+smallDelta;

    result = computeRadius(avRad, bigDelta, smallDelta, maxDsDistance, b, o, optim, base);  

    // OppdOuter;
    // result = computeRadius(avRad-bigDelta-smallDelta, 2*(bigDelta-smallDelta), o, b, optim, base);
    // result.first -= bigDelta+smallDelta;
  }

  return result;
}


void BarrelLayer::buildLayer(double averageRadius,
                             std::pair<double, double> worstCaseRadii,
                             double bigDelta,
                             double smallDelta,
                             const vector<double>& dsDistances,
                             double overlap,
                             double safetyOrigin,
                             const ModulePlacementStrategy& moduleStrategy,
                             int pushDirection,      // Direction where the average
                             // radius is pushed to optimize module usage
                             // +1: towards outside, 0 no optimize, -1 towards inside
                             int base,               // Modules base number
                             bool stringSameParity,  // Wether the strings should have
                             // the same or opposite module parity (in/out)
                             BarrelModule* sampleModule,
                             int sectioned, // = NoSection  (deprecated)
                             double farthestZ // = 0 not zero in case you build a mezzanine
                            ) {

  //  int stringParity;      // Parity of string in the shell (inner, outer)
  //  int stringSmallParity; // Parity of module in the string
  std::pair <double, int> optimalBarrel;
  double goodRadius;
  int nStrings;

  nOfRods_ = 0;
  nModsOnString_ = 0;
  optimalBarrel = layerPhi(averageRadius   , // tentativeX    
                           bigDelta,
                           smallDelta,
                           *std::max_element(dsDistances.begin(), dsDistances.end()), // maximum sensor spacing in the rod
                           overlap,
                           sampleModule->getWidth() , // Module half width
                           pushDirection, // wether to shrink or enlarge (>0 means enlarge)
                           base,
                           stringSameParity
                          );

  if (worstCaseRadii <= std::make_pair(0.0, 0.0)) { // if we pass 0.0 as worstCaseRadii, the layer is built the traditional way, i.e. around the radius found by layerPhi (which is close to averageRadius)
    worstCaseRadii = std::make_pair(optimalBarrel.first, optimalBarrel.first);
  }                                                 // else the layer is built placing the modules taking the supplied worstCaseRadii into account. This is to create multiple (over-hermetic) layers with identical rods
    
  goodRadius = optimalBarrel.first;
  nStrings = optimalBarrel.second;
  averageRadius_ = optimalBarrel.first;

  tempString.str(""); tempString << "GoodRadius: " << goodRadius
    << ", nStrings:   " << nStrings;
  logINFO(tempString);

  if (nStrings%2!=0) {
    logWARNING("You just asked for a layer with a number of strings not multiple of 2");
  }

  double stringPhiShift = 2*M_PI/double(nStrings);
  double rightAngle=0;
  std::vector<Module*>::iterator itMod;

  // build archetypical strings
  // build two strings at average radius (one with inverted parity), traslate one inwards and one outwards by big delta

  int smallParity = (smallDelta > 0) - (smallDelta < 0);  // extract sign

  std::vector<Module*> archetypeIn, archetypeOut;
  if (farthestZ==0) { // Build a standard (double) string
    buildStringPair(archetypeIn,
                    goodRadius,
                    worstCaseRadii, // wanna-be in
                    bigDelta,
                    smallDelta,
                    dsDistances,
                    overlap,
                    safetyOrigin,
                    moduleStrategy,
                    smallParity,
                    sampleModule);
    buildStringPair(archetypeOut, // wanna-be out
                    goodRadius,
                    worstCaseRadii,
                    bigDelta,
                    smallDelta,
                    dsDistances,
                    overlap,
                    safetyOrigin,
                    moduleStrategy,
                    stringSameParity ? smallParity : -smallParity,
                    sampleModule);
  } else { // Build a mezzanine (double) string
    buildMezzanineStringPair(archetypeIn, // wanna-be in
                             goodRadius,
                             worstCaseRadii,
                             bigDelta,
                             smallDelta,
                             dsDistances,
                             overlap,
                             safetyOrigin,
                             farthestZ,
                             moduleStrategy,
                             smallParity,
                             sampleModule);
    buildMezzanineStringPair(archetypeOut, // wanna-be out
                             goodRadius,
                             worstCaseRadii,
                             bigDelta,
                             smallDelta,
                             dsDistances,
                             overlap,
                             safetyOrigin,
                             farthestZ,
                             moduleStrategy,
                             stringSameParity ? smallParity : -smallParity,
                             sampleModule);
  }

  for (std::vector<Module*>::iterator it = archetypeIn.begin(); it < archetypeIn.end(); ++it) { // move archetype inward
    (*it)->translate(XYZVector(0, -bigDelta, 0)); 
  }

  for (std::vector<Module*>::iterator it = archetypeOut.begin(); it < archetypeOut.end(); ++it) { // move archetype outward
    (*it)->translate(XYZVector(0, bigDelta, 0)); 
  }


  for (int i=0; i<nStrings; i++) {
    /* stringParity = (i%2)*2-1; 
       if (stringSameParity) {
       stringSmallParity = 1;
       } else {
       stringSmallParity = -1 * stringParity;
       }
    // std::cout << "Building a string" << std::endl; // debug
    */
    std::vector<Module*> aString = (i%2) ? archetypeOut : archetypeIn;

    /*     if (farthestZ==0) { // Build a standard (double) string
           buildStringPair(aString,
           goodRadius+stringParity*bigDelta,
           stringSmallParity*smallDelta,
           overlap,
           safetyOrigin,
           nModules,
           sampleModule);
           } else { // Build a mezzanine (single) string
           buildMezzanineString(aString,
           goodRadius+stringParity*bigDelta,
           stringSmallParity*smallDelta,
           overlap,
           safetyOrigin,
           nModules,
           sampleModule, farthestZ);
           }
           */   


    if (i==0) {
      //std::cerr << "Computing min edge: "; // debug
      edge phiEdge;
      BarrelModule* stdBarrelMod;
      for (itMod=aString.begin(); itMod!=aString.end(); itMod++) {
        if ( (stdBarrelMod=dynamic_cast<BarrelModule*>(*itMod)) ) {
          phiEdge=((BarrelModule*)(*itMod))->getEdgePhiSide(-1);
          if (itMod==aString.begin()) {
            rightAngle=phiEdge.first;
            //std::cerr << "first angle is " << rightAngle << " and then: "; // debug
          } else {
            //std::cerr << phiEdge.first << " ";// debug
            if (phiEdge.first<rightAngle) rightAngle=phiEdge.first;
          }
        } else {
          // This should never happen
          logERROR("In building a string: found a non-barrel module");
        }
      }
      //std::cerr << "done. The best is: " << rightAngle << std::endl; // debug
    }

    bool firstOnes;
    int aSection;
    std::vector<Module* >::iterator secondMod;
    secondMod=aString.begin()++;

    int zIndex = 0;

    for (itMod=aString.begin(); itMod!=aString.end(); itMod++) {
      aSection=NoSection;
      if ((itMod==aString.begin())||(itMod==secondMod)) {
        firstOnes=true;
      } else {
        firstOnes=false;
      }

      // TODO: find a way to conjugate the need for a simple geometry
      // for a fast simulation (first string up) and a real geometry with modules
      // not crossing the XZ plane
      //(*itMod)->rotatePhi(i*stringPhiShift-rightAngle);

      BarrelModule* tempMod = new BarrelModule((BarrelModule&)**itMod);
      tempMod->setContainerId(getContainerId());
      tempMod->rotatePhi(i*stringPhiShift);
      tempMod->setPhiIndex(i);
      if (i==0) {
        aSection |= YZSection;
      }
      if (firstOnes) {
        aSection |= XYSection;
        //  std::cout << "it's one of the first: it belongs to section" << XYSection << std::endl; // debug
      }

      switch (sectioned) {
      case NoSection:
        moduleSet_.push_back(tempMod);
        break;
      case XYSection:
        if (firstOnes) {
          moduleSet_.push_back(tempMod);
        }
        break;
      case YZSection:
        if (i==0) {
          //(*itMod)->rotatePhi(rightAngle);
          moduleSet_.push_back(tempMod);
        }
        break;
      }

      tempMod->setSection(aSection);

      zIndex++;

      //  if (firstOnes)
      //  std::cout << "its section is " << (*itMod)->getSection() << std::endl; //debug
    }
  }
  nOfRods_ = nStrings;
  nModsOnString_ = archetypeIn.size(); // archetypeIn and Out will have the same number of modules?
  if (archetypeIn.size() != archetypeOut.size()) { logWARNING("buildLayer(): Archetypes for inner and outer rods of layer at radius " + any2str(averageRadius) + " have different lengths"); }

  while (!archetypeIn.empty()) {
    delete archetypeIn.back();
    archetypeIn.pop_back();
  }
  while (!archetypeOut.empty()) {
    delete archetypeOut.back();
    archetypeOut.pop_back();
  }
}


// Always look for this plot when changing the geometry!
// Execute the program with 2> aaa
// and then "sort -u aaa"
void BarrelLayer::neededModulesPlot(double smallDelta, // Half distance between modules in the same string
                                    double bigDelta,   // String half gap
                                    double o,          // overlap
                                    double b,          // Module half width
                                    int base
                                   ) {
  double minR=400;
  double maxR=1100;
  int steps=50000;


  //   TProfile* nmod[2];

  //   char name[100];
  //   sprintf(name, "NeededModulesParitySame");
  //   nmod[0] = new TProfile(name, name, steps+1, minR, maxR);
  //   nmod[0]->SetLineColor(1);

  //   sprintf(name, "NeededModulesParityOppd");
  //   nmod[1] = new TProfile(name, name, steps+1, minR, maxR);
  //   nmod[1]->SetLineColor(2);

  std::pair<double, int> SameOuter, SameInner, OppdOuter, OppdInner;


  int maxnSame, maxnOppd;

  for (double x=minR; x<maxR; x+=(maxR-minR)/double(steps)) {
    SameOuter = computeRadius(x-bigDelta+smallDelta, 2*bigDelta, o, b, FIXED, base);
    SameInner = computeRadius(x-bigDelta-smallDelta, 2*bigDelta, o, b, FIXED, base);
    OppdInner = computeRadius(x-bigDelta+smallDelta, 2*(bigDelta-smallDelta), o, b, FIXED, base);
    OppdOuter = computeRadius(x-bigDelta-smallDelta, 2*(bigDelta-smallDelta), o, b, FIXED, base);

    maxnSame = (SameOuter.second > SameInner.second) ? SameOuter.second : SameInner.second ;
    maxnOppd = (OppdOuter.second > OppdInner.second) ? OppdOuter.second : OppdInner.second ;

    if (SameOuter.second > SameInner.second) {
      logDEBUG("SameOuter is bigger than SameInner");
    }
    if (SameOuter.second < SameInner.second) {
      logDEBUG("SameInner is bigger than SameOuter");
    }

    if (OppdOuter.second > OppdInner.second) {
      logDEBUG("OppdOuter is bigger than OppdInner");
    }
    if (OppdOuter.second < OppdInner.second) {
      logDEBUG("OppdInner is bigger than OppdOuter");
    }

    if (maxnSame>maxnOppd) logDEBUG("Same is bigger than Oppd");
    if (maxnSame<maxnOppd) logDEBUG("Oppd is bigger than Same");
  }

  return;
}

int BarrelLayer::cutOverEta(double etaCut) {
  int nCut = 0;
  ModuleVector::iterator modIt;

  double theta;
  double eta;

  //double zmax = Layer::getMaxZ();
  //double zmin = getMinZ();

  int minRingDelete=0;
  bool firstCutFound=false;
  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    theta=(*modIt)->getMeanTheta();
    eta = -1*log(tan(theta/2.));
    if (fabs(eta)>etaCut) {
      if (!firstCutFound) {
        minRingDelete=(*modIt)->getRing();
        firstCutFound=true;
      } else {
        if ((*modIt)->getRing() < minRingDelete) minRingDelete = (*modIt)->getRing();
      }
    }
  }

  if (firstCutFound) {
    // Bookkeeping
    nModsOnString_ = minRingDelete -1;
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); ) {
      if ((*modIt)->getRing()>=minRingDelete) {
        // Erase the useless module!
        delete (*modIt);
        moduleSet_.erase(modIt);
        nCut++;
      } else {
        modIt++;
      }
    }
  }
  return nCut;
}

double BarrelLayer::getMaxZ(int direction) {
  double maxZ;
  //double minZ;
  double aZ;
  ModuleVector::iterator modIt;
  BarrelModule* aBarrelModule;
  bool firstModule=true;

  if (direction==0) {
    logERROR("BarrelLayer::getMaxZ was called with direction == 0");
    return 0;
  }
  direction/=int(fabs(direction));

  maxZ=0;

  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
      aZ = aBarrelModule->getEdgeZSide(direction).first;
      if (firstModule) {
        firstModule=false;
        maxZ=aZ;
      }
      if (direction*aZ>direction*maxZ) {
        maxZ=aZ;
      }
    }
  }

  return maxZ;
}



void BarrelLayer::compressToZ(double newMaxZ) {
  double maxZ;
  double minZ;
  double aZ;
  ModuleVector::iterator modIt;
  BarrelModule* aBarrelModule;
  BarrelModule* maxBarrelModule=NULL;
  BarrelModule* minBarrelModule=NULL;
  bool firstModule=true;

  maxZ=0;
  minZ=0;

  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
      aZ = aBarrelModule->getEdgeZSide(+1).first;
      if ((aZ>maxZ)||firstModule) {
        maxBarrelModule = aBarrelModule;
        maxZ=aZ;
      }
      aZ = aBarrelModule->getEdgeZSide(-1).first;
      if ((aZ<minZ)||firstModule) {
        minBarrelModule = aBarrelModule;
        minZ=aZ;
      }
      if (firstModule) firstModule=false;
    }
  }

  maxZ=fabs(maxZ);
  minZ=fabs(minZ);

  //   std::cout << "The layer's highest z+ is " << maxZ << std::endl; //debug
  //   std::cout << "The layer's highest z- is " << minZ << std::endl; //debug

  // Measure the Delta of the farthest module
  double Deltap;
  double Deltam;
  Deltap = fabs(newMaxZ) - maxZ;
  Deltam = -1*fabs(newMaxZ) + minZ;

  // std::cerr << "Delta-: " << Deltam << std::endl; // debug
  // std::cerr << "Delta+: " << Deltap << std::endl; // debug
  // TODO: raise an alarm if Delta > 0

  // delta is the ratio between Delta and the module's position
  double deltam=0;
  double deltap=0;
  double myMeanZ;
  XYZVector modShift;
  if (maxBarrelModule!=NULL) {
    deltam=Deltam/((minBarrelModule->getMeanPoint()).Z());
    deltap=Deltap/((maxBarrelModule->getMeanPoint()).Z());
    // std::cerr << "delta-: " << deltam << std::endl; // debug
    // std::cerr << "delta+: " << deltap << std::endl; // debug

    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
      if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
        myMeanZ = ((aBarrelModule->getMeanPoint()).Z());
        if (myMeanZ>0) {
          modShift=XYZVector(0, 0, deltap*myMeanZ);
        } else {
          modShift=XYZVector(0, 0, deltam*myMeanZ);
        }
        aBarrelModule->translate(modShift);
      }
    }
  } else {
    logERROR("int BarrelLayer::compactToZ couldn't find the maxZ module");
  }


  // TODO: add exit code
}

// This procedure assumes there is no module with a border < newMinZ
void BarrelLayer::compressExceeding(double newMaxZ, double newMinZ) {
  double maxZ;
  double aZ;
  ModuleVector::iterator modIt;
  BarrelModule* aBarrelModule;
  BarrelModule* maxBarrelModule=NULL;
  bool firstModule=true;

  int direction = int(newMaxZ/newMaxZ);

  maxZ=0;

  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
      aZ = aBarrelModule->getEdgeZSide(direction).first;
      if ((aZ>maxZ)||firstModule) {
        maxBarrelModule = aBarrelModule;
        maxZ=aZ;
      }
      if (firstModule) firstModule=false;
    }
  }

  // Measure the Delta of the farthest module
  double Delta;
  Delta = newMaxZ - maxZ;

  // TODO: if Deltam < 0 raise an error

  // delta is the ratio between Delta and the module's loewst Z w.r.t. newMinZ
  double delta=0;
  double myMinZ;
  XYZVector modShift;
  if (maxBarrelModule!=NULL) {
    delta=Delta/((maxBarrelModule->getEdgeZSide(direction*-1)).first-newMinZ);

    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
      if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
        myMinZ = ((aBarrelModule->getEdgeZSide(direction*-1)).first-newMinZ);
        if (myMinZ>0) {
          modShift=XYZVector(0, 0, delta*myMinZ);
        }
        aBarrelModule->translate(modShift);
      }
    }
  } else {
    logERROR("int BarrelLayer::compressExceeding couldn't find the maxZ module");
  }


  // TODO: add exit code
}

double BarrelLayer::computeAverageRadius() {
  ModuleVector::iterator modIt;
  BarrelModule* aBarrelModule;
  double myRadius=0;
  int nModules=0;

  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
      myRadius+=aBarrelModule->getMeanPoint().Rho();
      nModules++;
    }
  }

  myRadius /= nModules;
  averageRadius_=myRadius;
  return myRadius;
}

/******************/
/*                */
/* Endcap layer   */
/*                */
/******************/

EndcapLayer::~EndcapLayer() {
  if (sampleModule_) delete sampleModule_;
}

EndcapLayer::EndcapLayer() {
  sampleModule_ = NULL;
}

EndcapLayer::EndcapLayer(EndcapLayer& inputLayer) {
  ModuleVector::iterator modIt;
  ModuleVector inputModuleV = inputLayer.moduleSet_;
  EndcapModule* aModule;
  EndcapModule* stdModule;

  setName(inputLayer.layerName_, inputLayer.layerIndex_);
  setContainerName(inputLayer.containerName_);
  averageZ_  = inputLayer.averageZ_;
  nOfRings_ = inputLayer.getRings();
  nModsOnRing_.clear();
  for (unsigned int i = 0; i < inputLayer.getModulesOnRing().size(); i++) {
    nModsOnRing_.push_back(inputLayer.getModulesOnRing().at(i));
  }

  if ( (sampleModule_=dynamic_cast<EndcapModule*>(inputLayer.getSampleModule())) ) {
    sampleModule_ = new EndcapModule(*(inputLayer.getSampleModule()));
  }

  for (modIt=inputModuleV.begin(); modIt!=inputModuleV.end(); modIt++) {
    if ( (stdModule=dynamic_cast<EndcapModule*>(*modIt)) ) {
      aModule = new EndcapModule(*stdModule);
      moduleSet_.push_back(aModule);
    } else {
      logERROR("in EndcapLayer::EndcapLayer(EndcapLayer& inputLayer) I found a non-endcap module in the source");
    }
  }
}


EndcapLayer::EndcapLayer(const EndcapModule& sample, double alpha, double d) {
  sampleModule_ = new EndcapModule(sample, alpha, d, -1);
}

EndcapLayer::EndcapLayer(double alpha, double d) {
  sampleModule_ = new EndcapModule(alpha, d, -1);
}

EndcapLayer::EndcapLayer(const EndcapModule& mySample) {
  sampleModule_ = new EndcapModule(mySample);
}

//   void rotateY(double alpha);
void EndcapLayer::translateZ(const double& zShift) {
  ModuleVector::iterator modIt;

  averageZ_+=zShift;

  XYZVector myShift(0, 0, zShift);
  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
    (*modIt)->translate(myShift);
  }
}


#define MAX_LOOPS 10000

// y=r^2/l^2
// R: circle radius
// L: line origin (negative)
// x=sin^2(alpha)
// alpha is the line angle wrt x axis
double EndcapLayer::solvex(double y) {
  return(1/8.*(3*y+2-sqrt(9*pow(y, 2)-4*y+4)));
}

double EndcapLayer::gamma1(double x, double y, double r) {
  double result;
  result = sqrt(1-x) - sqrt(y-x);
  result *= r / sqrt(y);
  return result;
}

double EndcapLayer::gamma2(double x, double y, double r) {
  double result;
  result = sqrt(1-x) + sqrt(y-x);
  result *= r / sqrt(y);
  return result;
}

double EndcapLayer::Area(double x, double y, double r) {
  double result;
  result = sqrt((1-(x/y))*x/y);
  result *= 4*pow(r, 2)*(1-x);
  return result;
}

double EndcapLayer::compute_l(double x, double y, double d) {
  double result = d;
  result /= (1-x)-sqrt((1-x)*(y-x));
  return result;
}

double EndcapLayer::compute_d(double x, double y, double l) {
  double result = l;
  result *= (1-x)-sqrt((1-x)*(y-x));
  return result;
}


void EndcapLayer::buildSingleDisk(int nRings,     // top-to-bottom fixed nRings version
                                  double maxRadius,
                                  double smallDelta,
                                  double bigDelta,
                                  double diskZ,
                                  double overlap,
                                  double zError,
                                  const std::vector<double>& dsDistances,
                                  int phiSegments,
                                  bool oddSegments, bool alignEdges,
                                  std::map<int, EndcapModule*> sampleModules,
                                  std::map<int, int> ringDirectives,
                                  int diskParity,
                                  int sectioned /*=NoSection*/) {

  averageZ_=diskZ;

  std::map<int, int>::iterator aDirective;

  int ringParity;
  int nearDirection = int(diskZ/fabs(diskZ))*-1;
  int nRing;
  int addModules;

  double lastRho = 0;
  double nextRho = maxRadius;   
  //double nextRho = minRadius;
  //double destZ;
  XYZVector shiftThis;
  edge aSide;
  edge minSafetyEdge;
  std::ostringstream tag;

//  if ((nRings % 2) == 0) diskParity *= -1; // if rings are even, since we're going top to bottom, we have to invert the parity of the topmost to end up with the correct parity for the bottom one

  ringParity = diskParity;

  for (nRing=nRings; nRing>0; nRing--, ringParity*=-1) {
    EndcapModule* sampleModule = sampleModules[nRing];
    if (sampleModule==NULL) sampleModule = sampleModules[0];
    if (sampleModule==NULL) {
      std::cerr << "ERROR: a module prototype is not available for ring " << nRing
        << " i will not populate this ring" << std::endl;
      // TODO (important!): put a proper error handling here
      // For the moment it will just segfault (hi hi hi)
    }
    tempString.str(""); tempString  << "Ring number " << nRing;
    logDEBUG(tempString);
    // Ring parity is 1, -1, 1, -1 for
    //      nRing=    1,  2, 3,  4
    // if diskParity==1 and the opposite, otherwise

    sampleModule->setRing(nRing);

    // Debug
    tempString.str(""); tempString  << "Looking for directives of ring " << nRing;
    logINFO(tempString);

    aDirective = ringDirectives.find(nRing);
    if (aDirective!=ringDirectives.end()) {
      addModules = ringDirectives[nRing];
      logINFO("Found a directive");
    } else {
      addModules = 0;
      logINFO("Found no directive");
    }

    nextRho -= sampleModule->getHeight(); // we do this because the buildRing builds modules from bottom to top, so we tell it to start building the ring from the min radius
    // ringParity = 1 means the ring is nearer to the interaction point
    //
    int numPlacedModules; // set by buildring
    lastRho = buildRing(nextRho,
                        smallDelta,
                        ringParity*bigDelta,
                        diskZ,
                        overlap,
                        phiSegments,
                        oddSegments, alignEdges,
                        nearDirection,
                        sampleModule,
                        numPlacedModules,
                        maxRadius,
                        addModules,
                        sectioned);
    if (nModsOnRing_.size() < (unsigned)nRing) nModsOnRing_.resize(nRing);
    nModsOnRing_[nRing-1] = numPlacedModules;

    tempString.str(""); tempString << "smallDelta: " << smallDelta;
    logDEBUG(tempString);
    tempString.str(""); tempString << "bigDelta: " << bigDelta;
    logDEBUG(tempString);
    tempString.str(""); tempString << "lastrho: " << lastRho;
    logDEBUG(tempString);
    tempString.str(""); tempString << "aRingModule = new EndcapModule(*sampleModule, "
      << 100/lastRho << ", " << lastRho << ", -1);";
    logDEBUG(tempString);

    double newZ  = diskZ + (ringParity > 0 ? + bigDelta : - bigDelta) + smallDelta + dsDistances[nRing-1]/2;
    double lastZ = diskZ + (ringParity > 0 ? - bigDelta: + bigDelta) - smallDelta - dsDistances[nRing-1]/2;
    double originZ = ringParity > 0 ? zError : -zError;
    double nextRhoOrigin = (nextRho + overlap)/lastZ * newZ;
    double nextRhoShifted = nextRho/(lastZ - originZ) * (newZ - originZ);

    nextRho = nextRhoOrigin > nextRhoShifted ? nextRhoOrigin : nextRhoShifted;

  }

  nOfRings_ = nRings;

}

#define NEWENDCAP

#ifdef NEWENDCAP

void EndcapLayer::buildSingleDisk(double minRadius,
                                  double maxRadius,
                                  double smallDelta,
                                  double bigDelta,
                                  double diskZ,
                                  double overlap,
                                  double zError,
                                  const std::vector<double>& dsDistances,
                                  int phiSegments,
                                  bool oddSegments, bool alignEdges,
                                  std::map<int, EndcapModule*> sampleModules,
                                  std::map<int, int> ringDirectives,
                                  int diskParity,
                                  int sectioned /*=NoSection*/) {

  averageZ_=diskZ;

  //EndcapModule* aRingModule;

  //EndcapModule* trialModule[3];

  std::map<int, int>::iterator aDirective;


  int ringParity;
  int nearDirection = int(diskZ/fabs(diskZ))*-1;
  int nRing;
  int addModules;

  double lastRho = 0;
  double nextRho = minRadius;

  for (nRing=1; lastRho<maxRadius; nRing++) {
    EndcapModule* sampleModule = sampleModules[nRing];
    if (sampleModule==NULL) sampleModule = sampleModules[0];
    if (sampleModule==NULL) {
      std::cerr << "ERROR: a module prototype is not available for ring " << nRing
        << " i will not populate this ring" << std::endl;
      // TODO (important!): put a proper error handling here
      // For the moment it will just segfault (hi hi hi)
    }
    tempString.str(""); tempString  << "Ring number " << nRing;
    logDEBUG(tempString);
    // Ring parity is 1, -1, 1, -1 for
    //      nRing=    1,  2, 3,  4
    // if diskParity==1 and the opposite, otherwise
    ringParity = diskParity*(((nRing%2)*2)-1);

    sampleModule->setRing(nRing);

    // Debug
    tempString.str(""); tempString  << "Looking for directives of ring " << nRing;
    logINFO(tempString);

    aDirective = ringDirectives.find(nRing);
    if (aDirective!=ringDirectives.end()) {
      addModules = ringDirectives[nRing];
      logINFO("Found a directive");
    } else {
      addModules = 0;
      logINFO("Found no directive");
    }

    int numPlacedModules;
    // ringParity = 1 means the ring is nearer to the interaction point
    lastRho = buildRing(nextRho,
                        smallDelta,
                        ringParity*bigDelta,
                        diskZ,
                        overlap,
                        phiSegments,
                        oddSegments, alignEdges,
                        nearDirection,
                        sampleModule,
                        numPlacedModules,
                        maxRadius,
                        addModules,
                        sectioned);

    if (nModsOnRing_.size() < (unsigned)nRing) nModsOnRing_.resize(nRing);
    nModsOnRing_[nRing-1] = numPlacedModules;

    tempString.str(""); tempString << "smallDelta: " << smallDelta;
    logDEBUG(tempString);
    tempString.str(""); tempString << "bigDelta: " << bigDelta;
    logDEBUG(tempString);
    tempString.str(""); tempString << "lastrho: " << lastRho;
    logDEBUG(tempString);
    tempString.str(""); tempString << "aRingModule = new EndcapModule(*sampleModule, "
      << 100/lastRho << ", " << lastRho << ", -1);";
    logDEBUG(tempString);

    double newZ  = diskZ + (ringParity > 0 ? + bigDelta : - bigDelta) - smallDelta - dsDistances[nRing-1]/2;  
    double lastZ = diskZ + (ringParity > 0 ? - bigDelta : + bigDelta) + smallDelta + dsDistances[nRing-1]/2;
    double originZ = ringParity > 0 ? -zError : +zError;
    double nextRhoOrigin = (lastRho - overlap)/lastZ * newZ;
    double nextRhoOrigin2 = (lastRho)/lastZ * newZ - overlap;
    double nextRhoShifted = lastRho/(lastZ - originZ) * (newZ - originZ);

    nextRho = nextRhoOrigin < nextRhoShifted ? nextRhoOrigin : nextRhoShifted;
  
    nextRhoOrigin2 = nextRhoOrigin2;
  }

  nOfRings_ = nRing - 1;

}

#else

void EndcapLayer::buildSingleDisk(double minRadius,
                                  double maxRadius,
                                  double smallDelta,
                                  double bigDelta,
                                  double diskZ,
                                  double overlap,
                                  double zError,
                                  const std::vector<double>& dsDistances,
                                  int phiSegments,
                                  bool oddSegments, bool alignEdges,
                                  std::map<int, EndcapModule*> sampleModules,
                                  std::map<int, int> ringDirectives,
                                  int diskParity,
                                  int sectioned /*=NoSection*/) {

  averageZ_=diskZ;

  EndcapModule* aRingModule;

  EndcapModule* trialModule[3];

  std::map<int, int>::iterator aDirective;


  int ringParity;
  int nearDirection = int(diskZ/fabs(diskZ))*-1;
  int nRing;
  int addModules;

  double lastRho = 0;
  double nextRho = minRadius;
  double destZ;
  XYZVector shiftThis;
  edge aSide;
  edge minSafetyEdge;
  std::ostringstream tag;

  for (nRing=1; lastRho<maxRadius; nRing++) {
    EndcapModule* sampleModule = sampleModules[nRing];
    if (sampleModule==NULL) sampleModule = sampleModules[0];
    if (sampleModule==NULL) {
      std::cerr << "ERROR: a module prototype is not available for ring " << nRing
        << " i will not populate this ring" << std::endl;
      // TODO (important!): put a proper error handling here
      // For the moment it will just segfault (hi hi hi)
    }
    tempString.str(""); tempString  << "Ring number " << nRing;
    logDEBUG(tempString);
    // Ring parity is 1, -1, 1, -1 for
    //      nRing=    1,  2, 3,  4
    // if diskParity==1 and the opposite, otherwise
    ringParity = diskParity*(((nRing%2)*2)-1);

    sampleModule->setRing(nRing);

    // Debug
    tempString.str(""); tempString  << "Looking for directives of ring " << nRing;
    logINFO(tempString);

    aDirective = ringDirectives.find(nRing);
    if (aDirective!=ringDirectives.end()) {
      addModules = ringDirectives[nRing];
      logINFO("Found a directive");
    } else {
      addModules = 0;
      logINFO("Found no directive");
    }

    int numPlacedModules;
    // ringParity = 1 means the ring is nearer to the interaction point
    lastRho = buildRing(nextRho,
                        smallDelta,
                        ringParity*bigDelta,
                        diskZ,
                        overlap,
                        phiSegments,
                        oddSegments, alignEdges,
                        nearDirection,
                        sampleModule,
                        numPlacedModules,
                        maxRadius,
                        addModules,
                        sectioned);

    if (nModsOnRing_.size() < (unsigned)nRing) nModsOnRing_.resize(nRing);
    nModsOnRing_[nRing-1] = numPlacedModules;

    tempString.str(""); tempString << "smallDelta: " << smallDelta;
    logDEBUG(tempString);
    tempString.str(""); tempString << "bigDelta: " << bigDelta;
    logDEBUG(tempString);
    tempString.str(""); tempString << "lastrho: " << lastRho;
    logDEBUG(tempString);
    tempString.str(""); tempString << "aRingModule = new EndcapModule(*sampleModule, "
      << 100/lastRho << ", " << lastRho << ", -1);";
    logDEBUG(tempString);

    aRingModule = new EndcapModule(*sampleModule, 100/lastRho, lastRho, -1);

    aRingModule->setRing(nRing);
    tempString.str("");
    tempString << "aRingModule = " << std::endl;
    tempString << (*aRingModule);
    tempString << std::endl;
    tempString << "sampleModule = " << std::endl;
    tempString << (*sampleModule);
    tempString << std::endl;
    tempString << "aRingModule->getEdgeRhoSide(-1).first = "
      << aRingModule->getEdgeRhoSide(-1).first;
    logDEBUG(tempString);


    // Place the mock-up module so as to simulate the worst position for this ring
    shiftThis = XYZVector(0,
                          0,
                          diskZ + -1*nearDirection*smallDelta + ringParity*nearDirection*bigDelta);

    // Compute the destination z
    destZ = diskZ + nearDirection*smallDelta + -1*ringParity*nearDirection*bigDelta;

    aRingModule->translate(shiftThis);
    tempString.str(""); tempString << "Pushing from z="
      << diskZ + -1*nearDirection*smallDelta + ringParity*nearDirection*bigDelta
      << " to z=" << destZ;
    logDEBUG(tempString);
    tempString.str(""); tempString << "aRingModule->getEdgeRhoSide(-1).first = " << aRingModule->getEdgeRhoSide(-1).first;
    logDEBUG(tempString);

    for (int i=0; i<3; i++) {
      trialModule[i] = new EndcapModule(*aRingModule);
    }
    delete aRingModule;

    // Three tests: with a border and projecting form +- zError

    int direct;

    for (int i=0; i<3; i++) {
      direct=i-1;
      aSide = trialModule[i]->getEdgeRhoSide(-1);
      trialModule[i]->projectSideZ(aSide.second, destZ, direct*zError);
      tempString.str(""); tempString << "Fake module base: " << trialModule[i]->getEdgeRhoSide(-1).first;
      logDEBUG(tempString);

      if (i==1) {
        // trialModule[1] is the one with no zError, but overlap
        // TODO: here I will use the right setedgeside and getedgeside of the
        // endcapmodule, after they're debugged
        trialModule[i]->translate(XYZVector(0., -1*overlap, 0.));
      }

    }

    minSafetyEdge = trialModule[0]->getEdgeRhoSide(-1);
    minSafetyEdge.second = 0;

    // Debug:
    for (int i=0; i<3; i++) {
      tempString.str("");
      tempString << "safetyEdge["
        << i << "]="
        << trialModule[i]->getEdgeRhoSide(-1).first;
      logDEBUG(tempString);
    }
    for (int i=1; i<3; i++) {
      aSide = trialModule[i]->getEdgeRhoSide(-1);
      if (aSide.first<minSafetyEdge.first) {
        minSafetyEdge.first = aSide.first;
        minSafetyEdge.second = i;
      }
    }

    // Now we just keep the safest, while we delete the others
    for (int i=0; i<3; i++) {
      delete trialModule[i];
    }

    tempString.str(""); tempString << "Ring computation: using safety rule #" << minSafetyEdge.second;
    nextRho = minSafetyEdge.first;
    tempString << " next radius at rho: " << nextRho;
    logDEBUG(tempString); 
  }

  nOfRings_ = nRing - 1;

}

#endif

// Returns a 'standard' module
double EndcapLayer::buildRing(double minRadius,
                              double smallDelta,
                              double bigDelta,
                              double diskZ,
                              double overlap,
                              int phiSegments,
                              bool oddSegments, bool alignEdges,
                              int nearDirection,
                              EndcapModule* sampleModule,
                              int& numPlacedModules, 
                              double maxRadius /* = -1 */,
                              int addModules, /* = 0 */
                              int sectioned /* = NoSection */) {


  double alpha; // Module angular aperture
  bool wedges=false;

  // Used only for barrel modules
  double modMaxRadius=1;

  if ( sampleModule->getShape()==Module::Wedge ) {
    wedges=true;
    double r = sampleModule->getDiameter()/2.;

    tempString.str("");
    tempString << "Desired distance: " << minRadius << std::endl;
    tempString << "Desired overlap: " << overlap << std::endl;
    double l=(minRadius-r);
    tempString << "l: " << l << std::endl;
    tempString << "r: " << r << std::endl;
    double y=pow(r/l, 2);
    tempString << "y: " << y << std::endl;
    double x=solvex(y);
    tempString << "x: " << x << std::endl;
    tempString << "gamma1: " << gamma1(x, y, r) << std::endl;
    tempString << "gamma2: " << gamma2(x, y, r) << std::endl;

    // Max area
    tempString << "Area : " << Area(x, y, r) << std::endl;
    tempString << "Optimization: " << std::endl;


    double tempd;
    bool computeFinish=false;

    for (int i=0; i < MAX_LOOPS; i++) {
      l=compute_l(x, y, minRadius);
      y=pow(r/l, 2);
      x=solvex(y);

      tempd=compute_d(x, y, l);

      if (fabs(minRadius-tempd)<1e-15) {
        computeFinish=true;
        break;
      }

    }

    if (!computeFinish) {
      tempString << "Computation not finished" << std::endl;
      tempString << "Delta_d is now: " << minRadius-tempd << std::endl;
    }
    tempString << "x: " << x << std::endl;
    tempString << "Area : " << Area(x, y, r) << std::endl;

    alpha=asin(sqrt(x))*2;

    tempString << "mod aperture: " << alpha*180./M_PI;
    logDEBUG(tempString);
  } else if ( sampleModule->getShape()==Module::Rectangular ) { // it's a barrel module
    modMaxRadius=minRadius + sampleModule->getHeight();
    // widthHi or widthLo would give the same number: it's a square module!
    alpha=2*asin(sampleModule->getWidthHi()/2. / modMaxRadius);
  } else { // It's not a barrel or an endcap module
    tempString.str("");
    tempString << "The sample EndcapModule was not Barrel nor Endcap... this should never happen!";
    logERROR(tempString);
    alpha = 0;
  }

  // The needed overlap becomes an angle delta by
  // checking the unsafest point (r=r_min)
  double delta, effectiveAlpha;

  if (wedges) {
    delta = overlap / minRadius;
  } else {
    delta = overlap / modMaxRadius;
  }
  effectiveAlpha = alpha - delta;


  tempString.str(""); tempString << "overlap delta: " << delta*180./M_PI;
  logINFO(tempString);
  tempString.str(""); tempString << "mod aperture (effective): " << effectiveAlpha*180./M_PI;
  logINFO(tempString);


  double n;
  n = 2*M_PI / effectiveAlpha;

  tempString.str(""); tempString << "Number of modules: " << n;
  logINFO(tempString);

  // Optimal number of modules 
  int nOpt;

  if (wedges) {
    // We round to the nearest multiple of phiSegments
    //nOpt = int(floor((n/double(phiSegments))+.5)) * phiSegments;
    nOpt = round(n/double(phiSegments), oddSegments) * phiSegments;
    tempString.str(""); tempString << "I would use " << nOpt << " modules (optimized), ";
    nOpt += (addModules*phiSegments);
    tempString << "i will use " << nOpt << " modules (user request)";
  } else {
    // For square modules we can only increase the number of sensors to
    // cover the whole area
    nOpt = int(ceil(n/double(phiSegments))) * phiSegments;
    tempString.str("");
    if (oddSegments) {
      int nInSegment = nOpt/phiSegments;
      if ((nInSegment%2)==0) { // Number of modules in a segment is not zero (as it was required)
        nInSegment++; // The nearest possible even number
        tempString << "Number of modules in this ring was increased from " << nOpt;
        nOpt = nInSegment * phiSegments;
        tempString << " to " << nOpt << " to comply with user request of oddSegments.";
        logINFO(tempString);
      }
    }
    tempString.str("");
    tempString << "I will use " << nOpt << " modules";
  }
  logINFO(tempString);

  // Additional rotation to align the module edges
  double alignmentRotation;
  if (alignEdges) alignmentRotation = 0.5;
  else alignmentRotation = 0;
  //if (alignEdges) {
  //  if ((nOpt/phiSegments)%2==0) {
  //alignmentRotation = 0.5;
  // }
  //}

  double goodAlpha;
  goodAlpha = 2*M_PI/double(nOpt);
  goodAlpha += delta;

  //   XYZVector diskShift   = XYZVector(0, 0, diskZ);
  //   XYZVector petalShift  = XYZVector(0, 0, smallDelta);
  //   XYZVector ringShift   = XYZVector(0, 0, bigDelta);
  int ringParity;
  EndcapModule *myModule;
  //BarrelModule *myBarrelModule;
  //Module* myModule;

  //   std::cout << "diskShift:  " << diskZ << std::endl
  //      << "petalShift: " << smallDelta << std::endl
  //      << "ringShift:  " << bigDelta << std::endl;

  bool YZSectionTaken = false;
  for (int i=0; i<nOpt; i++) {
    ringParity = ((i%2)*2)-1;
    if (wedges) {
      myModule = new EndcapModule(*sampleModule, goodAlpha, minRadius, maxRadius);
    } else {
      myModule = new EndcapModule(*sampleModule, minRadius);
    }
    myModule->setContainerId(getContainerId());
    myModule->setPhiIndex(i);
    myModule->setZSide(signum(diskZ));
    myModule->rotatePhi(M_PI/2+2.*M_PI*((i+alignmentRotation)/double(nOpt)));
    XYZVector shift = XYZVector(0, 0, diskZ + nearDirection*ringParity*smallDelta + nearDirection*bigDelta);
    myModule->translate(shift);

    // Find the correct module for the YZ section
    double phase = double(i)/nOpt;
    bool inYZ = ((phase>=0.75)&&(!YZSectionTaken));
    //std::cout << "Phase = " << phase << " , inYZ = "  << inYZ << std::endl;
    if (inYZ) YZSectionTaken = true;

    if (inYZ) myModule->setSection(YZSection);
    if (sectioned == NoSection) {
      moduleSet_.push_back(myModule);
    } else {
      if (sectioned == YZSection) {
        if (inYZ) moduleSet_.push_back(myModule);
      }
    }
  }
  numPlacedModules = nOpt;

  double lastRho;
  if (wedges) {
    EndcapModule* aRingModule = new EndcapModule(*sampleModule, goodAlpha, minRadius, maxRadius);
    lastRho = minRadius + aRingModule->getHeight();
    tempString.str(""); tempString << "Actual module area: " << aRingModule->getArea();
    logINFO(tempString);
    // Just to be sure...!
    if (aRingModule->wasCut()) {
      tempString.str(""); tempString << "The ring Module was cut, losing: " << aRingModule->getLost();
      logINFO(tempString);
      lastRho=maxRadius;
    }
    delete aRingModule;
  } else {
    lastRho = modMaxRadius;
  }

  return lastRho;
}

int EndcapLayer::cutOverEta(double etaCut) {
  int nCut = 0;
  ModuleVector::iterator modIt;

  double theta;
  double eta;

  for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); ) {
    theta=(*modIt)->getMeanTheta();
    eta = -1*log(tan(theta/2.));
    if (fabs(eta)>etaCut) {
      // Bookkeeping
      decreaseModCount((*modIt)->getRing() - 1);
      // Erase the useless module!
      delete (*modIt);
      moduleSet_.erase(modIt);
      nCut++;
    } else {
      modIt++;
    }
  }

  return nCut;
}

double EndcapLayer::getMaxModuleThickness() {
  ModuleVector::iterator iter, guard = moduleSet_.end();
  double res = 0.0;
  for (iter = moduleSet_.begin(); iter != guard; iter++) {
    if ((*iter)->getModuleThickness() > res) res = (*iter)->getModuleThickness();
  }
  return res;
}

// Rounds to the nearest (odd) integer
// @param x number to round
// @param odd if true rounds to the nearest odd natural number
// (integer) instead of just the nearest integer
// @return the nearest (odd) integer to x
int Layer::round(const double& x, const bool& odd) {
  if (!odd) return int(floor(x+.5));
  else return round((x-1)/2, false)*2+1;
}

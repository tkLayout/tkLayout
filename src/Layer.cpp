#include "Layer.h"
#include "RodPair.h"
#include "messageLogger.h"
#include "ConversionStation.h"

void Layer::check() {
  PropertyObject::check();

  if (buildNumModules() > 0 && maxZ.state()) throw PathfulException("Only one between numModules and maxZ can be specified");
  else if (buildNumModules() == 0 && !maxZ.state()) throw PathfulException("At least one between numModules and maxZ must be specified");

}

void Layer::cutAtEta(double eta) { 
  for (auto& r : rods_) r.cutAtEta(eta); 
  rods_.erase_if([](const RodPair& r) { return r.numModules() == 0; }); // get rid of rods which have been completely pruned
}

double Layer::calculatePlaceRadius(int numRods,
                                   double bigDelta,
                                   double smallDelta,
                                   double dsDistance,
                                   double moduleWidth,
                                   double overlap) {

  double d = dsDistance/2;

  double f = (moduleWidth/2) - (overlap/2);

  double R = bigDelta + smallDelta + d;
  double S = -bigDelta + smallDelta + d;

  double T = tan(2*M_PI/numRods);

  double a = T;
  double b = R*T + S*T - 2*f;
  double c = R*S*T - R*f - S*f - T*f*f;

  double r = (-b + sqrt(b*b - 4*a*c))/(2*a);

  return r;
}

std::pair<float, int> Layer::calculateOptimalLayerParms(const RodTemplate& rodTemplate) {
                                                              
  // CUIDADO fix placeRadiusHint!!!!!
  double maxDsDistance = (*std::max_element(rodTemplate.begin(), 
                                            rodTemplate.end(), 
                                            [](const unique_ptr<BarrelModule>& m1, const unique_ptr<BarrelModule>& m2) { return m1->dsDistance() > m2->dsDistance(); } ))->dsDistance();
  float moduleWidth = (*rodTemplate.rbegin())->minWidth();
  float f = moduleWidth/2 - phiOverlap()/2;
  float gamma = atan(f/(placeRadiusHint() + bigDelta() + smallDelta() + maxDsDistance/2)) + atan(f/(placeRadiusHint() - bigDelta() + smallDelta() + maxDsDistance/2));
  float tentativeModsPerSegment = 2*M_PI/(gamma * phiSegments());

  float optimalRadius;
  int optimalModsPerSegment;

  switch (radiusMode()) {
  case SHRINK:
    optimalModsPerSegment = floor(tentativeModsPerSegment);
    optimalRadius = calculatePlaceRadius(optimalModsPerSegment*phiSegments(), bigDelta(), smallDelta(), maxDsDistance, moduleWidth, phiOverlap());
    break;
  case ENLARGE:
    optimalModsPerSegment = ceil(tentativeModsPerSegment);
    optimalRadius = calculatePlaceRadius(optimalModsPerSegment*phiSegments(), bigDelta(), smallDelta(), maxDsDistance, moduleWidth, phiOverlap());
    break;
  case FIXED:
    optimalModsPerSegment = ceil(tentativeModsPerSegment);
    optimalRadius = placeRadiusHint();
    break;
  case AUTO: {
    int modsPerSegLo = floor(tentativeModsPerSegment);
    int modsPerSegHi = ceil(tentativeModsPerSegment);
    float radiusLo = calculatePlaceRadius(modsPerSegLo*phiSegments(), bigDelta(), smallDelta(), maxDsDistance, moduleWidth, phiOverlap());
    float radiusHi = calculatePlaceRadius(modsPerSegHi*phiSegments(), bigDelta(), smallDelta(), maxDsDistance, moduleWidth, phiOverlap());

    if (fabs(radiusHi - placeRadiusHint()) < fabs(radiusLo - placeRadiusHint())) {
      optimalRadius = radiusHi;
      optimalModsPerSegment = modsPerSegHi;
    } else {
      optimalRadius = radiusLo;
      optimalModsPerSegment = modsPerSegLo;
    }
    break;
             }
  default:
    throw PathfulException("Invalid value for enum radiusMode");
  }

  return std::make_pair(optimalRadius, optimalModsPerSegment*phiSegments());
}



RodTemplate Layer::makeRodTemplate() {
  RodTemplate rodTemplate(buildNumModules() > 0 ? buildNumModules() : (!ringNode.empty() ? ringNode.rbegin()->first + 1 : 1)); // + 1 to make room for a default constructed module to use when building rods in case the rodTemplate vector doesn't have enough elements
  for (int i = 0; i < rodTemplate.size(); i++) {
    rodTemplate[i] = std::move(unique_ptr<BarrelModule>(GeometryFactory::make<BarrelModule>(GeometryFactory::make<RectangularModule>())));
    rodTemplate[i]->store(propertyTree());
    if (ringNode.count(i+1) > 0) rodTemplate[i]->store(ringNode.at(i+1));
    rodTemplate[i]->build();
  }
  return rodTemplate;
}


void Layer::buildStraight() {

  RodTemplate rodTemplate = makeRodTemplate();

  std::pair<double, int> optimalLayerParms = calculateOptimalLayerParms(rodTemplate);
  placeRadius_ = optimalLayerParms.first; 
  numRods_ = optimalLayerParms.second;
  if (!minBuildRadius.state() || !maxBuildRadius.state()) {
    minBuildRadius(placeRadius_);
    maxBuildRadius(placeRadius_);
  }

  float rodPhiRotation = 2*M_PI/numRods_;

  // FIRST ROD : assign common properties
  StraightRodPair* first = GeometryFactory::make<StraightRodPair>();
  first->myid(1);
  first->minBuildRadius(minBuildRadius()-bigDelta());
  first->maxBuildRadius(maxBuildRadius()+bigDelta());
  if (buildNumModules() > 0) first->buildNumModules(buildNumModules());
  else if (maxZ.state()) first->maxZ(maxZ());
  first->smallDelta(smallDelta());
  //first->ringNode = ringNode; // we need to pass on the contents of the ringNode to allow the RodPair to build the module decorators
  first->store(propertyTree());
  // SECOND ROD : copy first rod
  logINFO(Form("Copying rod %s", fullid(*this).c_str()));
  StraightRodPair* second = GeometryFactory::clone(*first);
  second->myid(2);

  // FIRST ROD : build and store
  bool isPlusBigDeltaRod = (bigParity() > 0);
  first->build(rodTemplate, isPlusBigDeltaRod);
  first->translateR(placeRadius_ + (isPlusBigDeltaRod ? bigDelta() : -bigDelta()));
  rods_.push_back(first);

  // SECOND ROD : assign other properties, build and store 
  if (!sameParityRods()) second->zPlusParity(first->zPlusParity()*-1);
  isPlusBigDeltaRod = (bigParity() < 0);
  second->build(rodTemplate, isPlusBigDeltaRod);
  second->translateR(placeRadius_ + (isPlusBigDeltaRod ? bigDelta() : -bigDelta()));
  second->rotateZ(rodPhiRotation);
  rods_.push_back(second);

  // All other Rods
  for (int i = 2; i < numRods_; i++) {
    RodPair* rod = i%2 ? GeometryFactory::clone(*second) : GeometryFactory::clone(*first); // clone rods
    rod->myid(i+1);
    rod->rotateZ(rodPhiRotation*(i%2 ? i-1 : i));
    rods_.push_back(rod);
  }
}

void Layer::buildTilted() {
  std::ifstream ifs(tiltedLayerSpecFile());
  if (ifs.fail()) throw PathfulException("Cannot open tilted modules spec file \"" + tiltedLayerSpecFile() + "\"");

  string line;
  vector<TiltedModuleSpecs> tmspecs1, tmspecs2;
  while(getline(ifs, line).good()) {
    if (line.empty()) continue;
    auto tokens = split<double>(line, " ", false);
    if (tokens.size() < 7) { logERROR("Failed parsing tilted barrel line: " + line); continue; };
    TiltedModuleSpecs t1{ tokens[0], tokens[1], tokens[2]*M_PI/180. };
    TiltedModuleSpecs t2{ tokens[3], tokens[4], tokens[5]*M_PI/180. };
    if (t1.valid()) tmspecs1.push_back(t1);
    if (t2.valid()) tmspecs2.push_back(t2);
    numRods_ = tokens[6]; // this assumes every row of the spec file has the same value for the last column (num rods in phi) 
  }
  ifs.close();

  buildNumModules(tmspecs1.size());

  RodTemplate rodTemplate = makeRodTemplate();

  float rodPhiRotation = 2*M_PI/numRods_;

  TiltedRodPair* first = GeometryFactory::make<TiltedRodPair>();
  first->myid(1);
  first->store(propertyTree());
  first->build(rodTemplate, tmspecs1, 1);
  rods_.push_back(first);

  TiltedRodPair* second = GeometryFactory::make<TiltedRodPair>();
  second->myid(2);
  second->store(propertyTree());
  second->build(rodTemplate, tmspecs2, 0);
  second->rotateZ(rodPhiRotation);
  rods_.push_back(second);

  for (int i = 2; i < numRods_; i++) {
    RodPair* rod = i%2 ? GeometryFactory::clone(*second) : GeometryFactory::clone(*first); // clone rods
    rod->myid(i+1);
    rod->rotateZ(rodPhiRotation*(i%2 ? i-1 : i));
    rods_.push_back(rod);
  }

  // computing the layer's place radius as the average of all the modules' radii
  placeRadius_  = std::accumulate(tmspecs1.begin(), tmspecs1.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
  placeRadius_ += std::accumulate(tmspecs2.begin(), tmspecs2.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
  placeRadius_ /= tmspecs1.size() + tmspecs2.size();

}

void Layer::build() {
  ConversionStation* conversionStation;

  try { 
    materialObject_.store(propertyTree());
    materialObject_.build();

    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    if (tiltedLayerSpecFile().empty()) buildStraight();
    else buildTilted();

    for (auto& currentStationNode : stationsNode) {
      conversionStation = new ConversionStation();
      conversionStation->store(currentStationNode.second);
      conversionStation->check();
      conversionStation->build();
      
      if(conversionStation->stationType() == ConversionStation::Type::FLANGE) {
        if(flangeConversionStation_ == nullptr) { //take only first defined flange station
          flangeConversionStation_ = conversionStation;
        }
      }else if(conversionStation->stationType() == ConversionStation::Type::SECOND) {
        secondConversionStations_.push_back(conversionStation);
      }
    }

        
    cleanup();
    builtok(true);

  } catch (PathfulException& pe) { 
    pe.pushPath(fullid(*this)); 
    throw; 
  }
}

const MaterialObject& Layer::materialObject() const{
  return materialObject_;
}

ConversionStation* Layer::flangeConversionStation() const {
  return flangeConversionStation_;
}

const std::vector<ConversionStation*>& Layer::secondConversionStations() const {
  return secondConversionStations_;
}


define_enum_strings(Layer::RadiusMode) = { "shrink", "enlarge", "fixed", "auto" };

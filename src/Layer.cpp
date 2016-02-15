#include "Layer.h"
#include "RodPair.h"
#include "messageLogger.h"
#include "ConversionStation.h"

void Layer::check() {
  PropertyObject::check();

  if (buildNumModules() > 0 && maxZ.state()) throw PathfulException("Only one between numModules and maxZ can be specified");
  else if (buildNumModules() == 0 && !maxZ.state()) throw PathfulException("At least one between numModules and maxZ must be specified");

  //if ((isTilted() && isTiltedAuto()) && !buildNumModulesFlat()) throw PathfulException("Automatic tilted layer : numModulesFlat must be specified");
  //if ((isTilted() && isTiltedAuto()) && !buildNumModulesTilted()) throw PathfulException("Automatic tilted layer : numModulesTilted must be specified");
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
  std::cout << "rodTemplate.size() = " << rodTemplate.size() << std::endl;
  for (int i = 0; i < rodTemplate.size(); i++) {
    rodTemplate[i] = std::move(unique_ptr<BarrelModule>(GeometryFactory::make<BarrelModule>(GeometryFactory::make<RectangularModule>())));
    rodTemplate[i]->store(propertyTree());
    if (ringNode.count(i+1) > 0) rodTemplate[i]->store(ringNode.at(i+1));
    rodTemplate[i]->build();
  }
  return rodTemplate;
}



TiltedRodTemplate Layer::makeTiltedRodTemplate(double flatPartThetaEnd) {
  TiltedRodTemplate tiltedRodTemplate;
  for (int i = 0; i < buildNumModulesTilted(); i++) {
    //std::cout << "i = " << i << std::endl;

    TiltedRing* tiltedRing = GeometryFactory::make<TiltedRing>();
    tiltedRing->myid(buildNumModulesFlat()+i+1);
    tiltedRing->store(propertyTree());
    if (ringNode.count(buildNumModulesFlat()+i+1) > 0) tiltedRing->store(ringNode.at(buildNumModulesFlat()+i+1));
    tiltedRing->numPhi(numModulesPhi());


    if ( i == 0 ) {
      std::cout << "flatPartThetaEnd = " << flatPartThetaEnd * 180. / M_PI << std::endl;
      tiltedRing->build(flatPartThetaEnd); 
    }
    else {
      //std::cout << "(tiltedRodTemplate.at(i-1))->thetaEnd() = " << (tiltedRodTemplate.at(i-1))->thetaEnd()<< std::endl;
      tiltedRing->build((tiltedRodTemplate.at(i-1))->thetaEnd()); 
    }


    
    //std::cout << "tiltedRing->outerRadius() = " << tiltedRing->outerRadius() << std::endl;
    //std::cout << "tiltedRing->zOuter() = " << tiltedRing->zOuter() << std::endl;
    //std::cout << "tiltedRing->tiltAngle() = " << tiltedRing->tiltAngle() << std::endl;
    //std::cout << "tiltedRodTemplate.at(i-1))->thetaEnd() = " << (tiltedRodTemplate.at(i-1))->thetaEnd() << std::endl;

    tiltedRodTemplate.push_back(tiltedRing);
  }
  return tiltedRodTemplate;
}



void Layer::buildStraight(bool isFlatPart) {

  RodTemplate rodTemplate = makeRodTemplate();

  if (!isFlatPart) {
    std::pair<double, int> optimalLayerParms = calculateOptimalLayerParms(rodTemplate);
    placeRadius_ = optimalLayerParms.first; 
    numRods_ = optimalLayerParms.second;
    numModulesPhi(numRods_);
  }
  else {
    placeRadius_ = placeRadiusHint();
    numRods_ = numModulesPhi();
  }
  if (!minBuildRadius.state() || !maxBuildRadius.state()) {
    minBuildRadius(placeRadius_);
    maxBuildRadius(placeRadius_);
  }

  float rodPhiRotation = 2*M_PI/numRods_;

  StraightRodPair* first = GeometryFactory::make<StraightRodPair>();
  first->myid(1);
  first->minBuildRadius(minBuildRadius()-bigDelta());
  first->maxBuildRadius(maxBuildRadius()+bigDelta());
  if (buildNumModules() > 0) first->buildNumModules(buildNumModules());
  else if (maxZ.state()) first->maxZ(maxZ());
  first->smallDelta(smallDelta());
  //first->ringNode = ringNode; // we need to pass on the contents of the ringNode to allow the RodPair to build the module decorators
  if (isFlatPart && buildNumModulesFlat() > 0) { first->zPlusParity( pow(-1, buildNumModulesFlat()) ); std::cout << pow(-1, buildNumModulesFlat()) << std::endl; }
  first->store(propertyTree());
  first->build(rodTemplate);

  logINFO(Form("Copying rod %s", fullid(*this).c_str()));
  StraightRodPair* second = GeometryFactory::clone(*first);
  second->myid(2);
  if (!sameParityRods()) second->zPlusParity(first->zPlusParity()*-1);

  first->translateR(placeRadius_ + (bigParity() > 0 ? bigDelta() : -bigDelta()));
  //first->translate(XYZVector(placeRadius_+bigDelta(), 0, 0));
  /*vector<TiltedModuleSpecs> coords = first->giveZPlusModulesCoords(1);
  for (int i = 0; i < coords.size(); i++) {
    std::cout << "tmspecs1[i].r = " << coords[i].r << "tmspecs1[i].z = " << coords[i].z << std::endl;
    }*/
  rods_.push_back(first);
  //if (isFlatPart) { flatRods_.push_back(first); }

  second->translateR(placeRadius_ + (bigParity() > 0 ? -bigDelta() : bigDelta()));
  //second->translate(XYZVector(placeRadius_-bigDelta(), 0, 0));
  second->rotateZ(rodPhiRotation);
  /*vector<TiltedModuleSpecs> coords2 = second->giveZPlusModulesCoords(1);
  for (int i = 0; i < coords2.size(); i++) {
  std::cout << "tmspecs2[i].r = " << coords2[i].r << "tmspecs2[i].z = " << coords2[i].z << std::endl;
  }*/
  flatPartThetaEnd_ = second->thetaEnd();
  rods_.push_back(second);
  //if (isFlatPart) { flatRods_.push_back(second); }

  if (!isFlatPart) {
    for (int i = 2; i < numRods_; i++) {
      RodPair* rod = i%2 ? GeometryFactory::clone(*second) : GeometryFactory::clone(*first); // clone rods
      rod->myid(i+1);
      rod->rotateZ(rodPhiRotation*(i%2 ? i-1 : i));
      rods_.push_back(rod);
      //if (isFlatPart) { flatRods_.push_back(rod); }
    }
  }
}

void Layer::buildTilted() {

  vector<TiltedModuleSpecs> tmspecs1, tmspecs2;


  if (!isTiltedAuto()) {
    std::ifstream ifs(tiltedLayerSpecFile());
    if (ifs.fail()) throw PathfulException("Cannot open tilted modules spec file \"" + tiltedLayerSpecFile() + "\"");

    string line;
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
  }

  else {

    flatPartThetaEnd_ = M_PI / 2.;

    if (buildNumModulesFlat() != 0) {
      buildNumModules(buildNumModulesFlat());
      buildStraight(true);

      //StraightRodPair* flatRod1 = rods_.begin();
      vector<TiltedModuleSpecs> info1 = (rods_.begin())->giveZPlusModulesCoords(buildNumModulesFlat());
      tmspecs1.insert(tmspecs1.end(), info1.begin(), info1.end());
      //StraightRodPair* flatRod2 = rods_.begin() + 1;
      vector<TiltedModuleSpecs> info2 = (rods_.begin() + 1)->giveZPlusModulesCoords(buildNumModulesFlat());
      tmspecs2.insert(tmspecs2.end(), info2.begin(), info2.end());

      rods_.clear();
    }

    numRods_ = numModulesPhi();

    std::cout << "myid() = " << myid() << std::endl;
    std::cout << "numRods_ = " << numRods_ << std::endl;
    //std::cout << "tmspecs1.end().r = " << tmspecs1[tmspecs1.size()-1].r << "tmspecs1.end().z = " << tmspecs1[tmspecs1.size()-1].z << std::endl;
    //std::cout << "tmspecs2.end().r = " << tmspecs2[tmspecs2.size()-1].r << "tmspecs2.end().z = " << tmspecs2[tmspecs2.size()-1].z << std::endl;
    

    TiltedRodTemplate tiltedRodTemplate = makeTiltedRodTemplate(flatPartThetaEnd_);

    for (int i = 0; i < tiltedRodTemplate.size(); i++) {
      TiltedModuleSpecs ti{ tiltedRodTemplate[i]->innerRadius(), tiltedRodTemplate[i]->zInner(), tiltedRodTemplate[i]->tiltAngle()*M_PI/180. };
      TiltedModuleSpecs to{ tiltedRodTemplate[i]->outerRadius(), tiltedRodTemplate[i]->zOuter(), tiltedRodTemplate[i]->tiltAngle()*M_PI/180. };
      /*std::cout << "i = " << i << std::endl;
      std::cout << "tiltAngle = " << tiltedRodTemplate[i]->tiltAngle()*M_PI/180. << std::endl;
      std::cout << "innerRadius = " << tiltedRodTemplate[i]->innerRadius() << std::endl;
      std::cout << "zInner = " << tiltedRodTemplate[i]->zInner() << std::endl;
      std::cout << "outerRadius = " << tiltedRodTemplate[i]->outerRadius() << std::endl;
      std::cout << "zOuter = " << tiltedRodTemplate[i]->zOuter() << std::endl;*/




      /*std::cout << "theta2 = " << tiltedRodTemplate[i]->thetaOuter() * 180. / M_PI << std::endl;
      std::cout << "idealTilt2 = " << tiltedRodTemplate[i]->tiltAngleIdealOuter() << std::endl;
      std::cout << "gap = " << tiltedRodTemplate[i]->gapR() << std::endl;
      std::cout << "avR = " << tiltedRodTemplate[i]->averageR() << std::endl;
      if (i >= 1) { std::cout << "cov1 = " << (tiltedRodTemplate[i]->thetaStartInner() - tiltedRodTemplate[i-1]->thetaEndInner()) * 180. / M_PI << std::endl; }
      if (i >= 1) { std::cout << "deltaz2 = " << tiltedRodTemplate[i]->zOuter() - tiltedRodTemplate[i-1]->zOuter() << std::endl; }

      if (i >= 1) {
	double zErrorAngle = atan( (tiltedRodTemplate[i]->rStartOuter_REAL() - tiltedRodTemplate[i-1]->rEndOuter_REAL()) / (tiltedRodTemplate[i-1]->zEndOuter_REAL() - tiltedRodTemplate[i]->zStartOuter_REAL()) );
	std::cout << "zError = " << tiltedRodTemplate[i]->zStartOuter_REAL() + tiltedRodTemplate[i]->rStartOuter_REAL() / tan(zErrorAngle) << std::endl; 
      }

      std::cout << atan(tiltedRodTemplate[i]->rEndOuter_REAL() / tiltedRodTemplate[i]->zEndOuter_REAL()) * 180. / M_PI << std::endl;*/

      if (ti.valid()) { tmspecs1.push_back(ti); }
      if (to.valid()) { tmspecs2.push_back(to); }
    }
    

  }




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

    //if (tiltedLayerSpecFile().empty()) buildStraight();
    if (!isTilted()) buildStraight(false);
    else buildTilted();

    for (auto& currentStationNode : stationsNode) {
      conversionStation = new ConversionStation();
      conversionStation->store(currentStationNode.second);
      conversionStation->check();
      conversionStation->build();
      
      if (conversionStation->stationType() == ConversionStation::Type::FLANGE) {
        if (flangeConversionStation_ == nullptr) { //take only first defined flange station
          flangeConversionStation_ = conversionStation;
        }
      } else if(conversionStation->stationType() == ConversionStation::Type::SECOND) {
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

#include "Layer.h"
#include "RodPair.h"
#include "messageLogger.h"
#include "ConversionStation.h"

Layer::TiltedRingsGeometryInfo::TiltedRingsGeometryInfo(int numModulesFlat, TiltedRingsTemplate tiltedRingsGeometry) {
  for (int i = (numModulesFlat + 2); i < (numModulesFlat + tiltedRingsGeometry.size() + 1); i++) {
    deltaZOuter_[i] = tiltedRingsGeometry[i]->zOuter() - tiltedRingsGeometry[i-1]->zOuter();

    //covInner_[i] = (tiltedRingsGeometry[i]->thetaStartInner() - tiltedRingsGeometry[i-1]->thetaEndInner());
    double zErrorInnerAngle = atan( (tiltedRingsGeometry[i]->rStartInner_REAL() - tiltedRingsGeometry[i-1]->rEndInner_REAL()) / (tiltedRingsGeometry[i]->zStartInner_REAL() - tiltedRingsGeometry[i-1]->zEndInner_REAL()) );
    zErrorInner_[i] = tiltedRingsGeometry[i]->zStartInner_REAL() - tiltedRingsGeometry[i]->rStartInner_REAL() / tan(zErrorInnerAngle);

    double zErrorOuterAngle = atan( (tiltedRingsGeometry[i]->rStartOuter_REAL() - tiltedRingsGeometry[i-1]->rEndOuter_REAL()) / (tiltedRingsGeometry[i]->zStartOuter_REAL() - tiltedRingsGeometry[i-1]->zEndOuter_REAL()) );
    zErrorOuter_[i] = tiltedRingsGeometry[i]->zStartOuter_REAL() - tiltedRingsGeometry[i]->rStartOuter_REAL() / tan(zErrorOuterAngle);
  }
}

void Layer::check() {
  PropertyObject::check();

  if (buildNumModules() > 0 && maxZ.state()) throw PathfulException("Only one between numModules and maxZ can be specified");
  else if (buildNumModules() == 0 && !maxZ.state()) throw PathfulException("At least one between numModules and maxZ must be specified");

  if (!isTilted()) {
    if (!phiOverlap.state()) throw PathfulException("Straight layer : phiOverlap must be specified.");
    if (numRods.state()) throw PathfulException("Straight layer : numRods should not be specified.");
    if (isTiltedAuto.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify isTiltedAuto. Not used.");
  }

  if (!isTilted() || (isTilted() && !isTiltedAuto())) {
    if (buildNumModulesFlat.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify numModulesFlat. Not used.");
    if (buildNumModulesTilted.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify numModulesTilted. Not used.");
  }

  if (isTilted()) {
    if (!isTiltedAuto.state()) throw PathfulException("Tilted layer : isTiltedAuto must be specified.");
    if (phiOverlap.state()) throw PathfulException("Tilted layer : phiOverlap should not be specified.");
    if (isTiltedAuto()) {    
      if (!buildNumModulesFlat.state()) throw PathfulException("Tilted layer with automatic placement : numModulesFlat must be specified.");
      if (!buildNumModulesTilted.state()) throw PathfulException("Tilted layer with automatic placement : numModulesTilted must be specified.");
      if (!numRods.state()) throw PathfulException("Tilted layer with automatic placement : numRods must be specified.");
    }
  }
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



TiltedRingsTemplate Layer::makeTiltedRingsTemplate(double flatPartThetaEnd) {
  TiltedRingsTemplate tiltedRingsGeometry;

  for (int i = (buildNumModulesFlat() + 1); i < (buildNumModulesFlat() + buildNumModulesTilted() + 1); i++) {

    TiltedRing* tiltedRing = GeometryFactory::make<TiltedRing>();
    tiltedRing->myid(i);
    tiltedRing->store(propertyTree());
    if (ringNode.count(i) > 0) tiltedRing->store(ringNode.at(i));
    tiltedRing->numPhi(numRods());

    double lastThetaEnd;
    if (i == (buildNumModulesFlat() + 1)) lastThetaEnd = flatPartThetaEnd; 
    else lastThetaEnd = tiltedRingsGeometry[i-1]->thetaEnd();
 
    tiltedRing->build(lastThetaEnd); 
    tiltedRingsGeometry[i] = tiltedRing;
  }

  return tiltedRingsGeometry;
}



void Layer::buildStraight(bool isFlatPart) {

  RodTemplate rodTemplate = makeRodTemplate();

  if (!isFlatPart) {
    std::pair<double, int> optimalLayerParms = calculateOptimalLayerParms(rodTemplate);
    placeRadius_ = optimalLayerParms.first; 
    numRods(optimalLayerParms.second);
  }
  else {
    placeRadius_ = placeRadiusHint();
  }

  if (!minBuildRadius.state() || !maxBuildRadius.state()) {
    minBuildRadius(placeRadius_);
    maxBuildRadius(placeRadius_);
  }

  float rodPhiRotation = 2*M_PI/numRods();

  StraightRodPair* first = GeometryFactory::make<StraightRodPair>();
  first->myid(1);
  first->minBuildRadius(minBuildRadius()-bigDelta());
  first->maxBuildRadius(maxBuildRadius()+bigDelta());
  if (buildNumModules() > 0) first->buildNumModules(buildNumModules());
  else if (maxZ.state()) first->maxZ(maxZ());
  first->smallDelta(smallDelta());
  //first->ringNode = ringNode; // we need to pass on the contents of the ringNode to allow the RodPair to build the module decorators
  if (isFlatPart) { first->zPlusParity( pow(-1, buildNumModulesFlat()) ); }
  first->store(propertyTree());
  first->build(rodTemplate);

  logINFO(Form("Copying rod %s", fullid(*this).c_str()));
  StraightRodPair* second = GeometryFactory::clone(*first);
  second->myid(2);
  if (!sameParityRods()) second->zPlusParity(first->zPlusParity()*-1);

  first->translateR(placeRadius_ + (bigParity() > 0 ? bigDelta() : -bigDelta()));
  if (!isFlatPart) { rods_.push_back(first); }
  else { flatPartRods_.push_back(first); }

  second->translateR(placeRadius_ + (bigParity() > 0 ? -bigDelta() : bigDelta()));
  second->rotateZ(rodPhiRotation);
  if (!isFlatPart) { rods_.push_back(second); }
  else { flatPartRods_.push_back(second); }

  


  for (int i = 2; i < numRods(); i++) {
    StraightRodPair* rod = i%2 ? GeometryFactory::clone(*second) : GeometryFactory::clone(*first); // clone rods
    rod->myid(i+1);
    rod->rotateZ(rodPhiRotation*(i%2 ? i-1 : i));
    if (!isFlatPart) { rods_.push_back(rod); }
    else { flatPartRods_.push_back(rod); }
  }

}

void Layer::buildTilted() {

  vector<TiltedModuleSpecs> tmspecsi, tmspecso;


  if (!isTiltedAuto()) {
    std::ifstream ifs(tiltedLayerSpecFile());
    if (ifs.fail()) throw PathfulException("Cannot open tilted modules spec file \"" + tiltedLayerSpecFile() + "\"");

    string line;
    while(getline(ifs, line).good()) {
      if (line.empty()) continue;
      auto tokens = split<double>(line, " ", false);
      if (tokens.size() < 7) { logERROR("Failed parsing tilted barrel line: " + line); continue; };
      TiltedModuleSpecs ti{ tokens[0], tokens[1], tokens[2]*M_PI/180. };
      TiltedModuleSpecs to{ tokens[3], tokens[4], tokens[5]*M_PI/180. };
      if (ti.valid()) tmspecsi.push_back(ti);
      if (to.valid()) tmspecso.push_back(to);
      numRods(tokens[6]); // this assumes every row of the spec file has the same value for the last column (num rods in phi) 
    }
    ifs.close();
  }

  else {

    double flatPartThetaEnd = M_PI / 2.;

    if (buildNumModulesFlat() != 0) {

      buildNumModules(buildNumModulesFlat());
      buildStraight(true);


      if (flatPartRods_.size() >= 2) {
	StraightRodPair* flatPartRod1 = flatPartRods_.front();
	const auto& zPlusModules1 = flatPartRod1->modules().first;
	for (const auto& m : zPlusModules1) {
	  TiltedModuleSpecs t1{m.center().Rho(), m.center().Z(), 0.0};
	  if (t1.valid()) (bigParity() > 0 ? tmspecso.push_back(t1) : tmspecsi.push_back(t1));
	}

	StraightRodPair* flatPartRod2 = flatPartRods_.at(1);
	const auto& zPlusModules2 = flatPartRod2->modules().first;
	for (const auto& m : zPlusModules2) {
	  TiltedModuleSpecs t2{m.center().Rho(), m.center().Z(), 0.0};
	  if (t2.valid()) (bigParity() > 0 ? tmspecsi.push_back(t2) : tmspecso.push_back(t2));
	}

	flatPartThetaEnd = (bigParity() > 0 ? flatPartRod1->thetaEnd() : flatPartRod2->thetaEnd());
      }
      else { logERROR(to_string(flatPartRods_.size()) + " straight rod was built for the whole flat part."); }     
    }


    if (tmspecsi.size() != buildNumModulesFlat()) {
      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " but flat part inner rod has " + to_string(tmspecsi.size()) + " module(s).");
    }
    if (tmspecso.size() != buildNumModulesFlat()) {
      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " but flat part outer rod has " + to_string(tmspecso.size()) + " module(s).");
    }
    



    tiltedRingsGeometry_ = makeTiltedRingsTemplate(flatPartThetaEnd);

    if (tiltedRingsGeometry_.size() == buildNumModulesTilted()) {
      for (int i = 0; i < buildNumModulesTilted(); i++) {
	int ringNumber = buildNumModulesFlat() + 1 + i;
	TiltedModuleSpecs ti{ tiltedRingsGeometry_[ringNumber]->innerRadius(), tiltedRingsGeometry_[ringNumber]->zInner(), tiltedRingsGeometry_[ringNumber]->tiltAngle()*M_PI/180. };
	TiltedModuleSpecs to{ tiltedRingsGeometry_[ringNumber]->outerRadius(), tiltedRingsGeometry_[ringNumber]->zOuter(), tiltedRingsGeometry_[ringNumber]->tiltAngle()*M_PI/180. };
	/*std::cout << "i = " << i << std::endl;
	  std::cout << "tiltAngle = " << tiltedRingsGeometry_[i]->tiltAngle()*M_PI/180. << std::endl;
	  std::cout << "innerRadius = " << tiltedRingsGeometry_[i]->innerRadius() << std::endl;
	  std::cout << "zInner = " << tiltedRingsGeometry_[i]->zInner() << std::endl;
	  std::cout << "outerRadius = " << tiltedRingsGeometry_[i]->outerRadius() << std::endl;
	  std::cout << "zOuter = " << tiltedRingsGeometry_[i]->zOuter() << std::endl;*/




	/*std::cout << "theta2 = " << tiltedRingsGeometry_[i]->thetaOuter() * 180. / M_PI << std::endl;
	  std::cout << "idealTilt2 = " << tiltedRingsGeometry_[i]->tiltAngleIdealOuter() << std::endl;
	  std::cout << "gap = " << tiltedRingsGeometry_[i]->gapR() << std::endl;
	  std::cout << "avR = " << tiltedRingsGeometry_[i]->averageR() << std::endl;
	  if (i >= 1) { std::cout << "cov1 = " << (tiltedRingsGeometry_[ringNumber]->thetaStartInner() - tiltedRingsGeometry_[ringNumber-1]->thetaEndInner()) * 180. / M_PI << std::endl; }
	  if (i >= 1) { std::cout << "deltaz2 = " << tiltedRingsGeometry_[i]->zOuter() - tiltedRingsGeometry_[i-1]->zOuter() << std::endl; }

	  if (i >= 1) {
	  double zErrorAngle = atan( (tiltedRingsGeometry_[i]->rStartOuter_REAL() - tiltedRingsGeometry_[i-1]->rEndOuter_REAL()) / (tiltedRingsGeometry_[i-1]->zEndOuter_REAL() - tiltedRingsGeometry_[i]->zStartOuter_REAL()) );
	  std::cout << "zError = " << tiltedRingsGeometry_[i]->zStartOuter_REAL() + tiltedRingsGeometry_[i]->rStartOuter_REAL() / tan(zErrorAngle) << std::endl; 
	  }

	  std::cout << "cov2 = " << atan(tiltedRingsGeometry_[i]->rEndOuter_REAL() / tiltedRingsGeometry_[i]->zEndOuter_REAL()) * 180. / M_PI << std::endl;*/

	if (ti.valid()) tmspecsi.push_back(ti);
	if (to.valid()) tmspecso.push_back(to);
      }
      tiltedRingsGeometryInfo_ = TiltedRingsGeometryInfo(buildNumModulesFlat(), tiltedRingsGeometry_);
    } 
    else {
      logERROR("numModulesTilted = " + to_string(buildNumModulesTilted()) + " but tilted rings geometry template has " + to_string(tiltedRingsGeometry_.size()) + " elements.");
    }


    buildNumModules(buildNumModulesFlat() + buildNumModulesTilted());

    if (tmspecsi.size() != buildNumModules()) {
      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " and numModulesTilted = " + to_string(buildNumModulesTilted()) + " but tilted rod 1 has " + to_string(tmspecsi.size()) + " module(s) in total.");
    }
    if (tmspecso.size() != buildNumModules()) {
      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " and numModulesTilted = " + to_string(buildNumModulesTilted()) + " but tilted rod 2 has " + to_string(tmspecso.size()) + " module(s) in total.");
    }

  }

  RodTemplate rodTemplate = makeRodTemplate();

  float rodPhiRotation = 2*M_PI/numRods();

  TiltedRodPair* first = GeometryFactory::make<TiltedRodPair>();
  first->myid(1);
  first->store(propertyTree());
  first->build(rodTemplate, tmspecsi, 1);
  rods_.push_back(first);

  TiltedRodPair* second = GeometryFactory::make<TiltedRodPair>();
  second->myid(2);
  second->store(propertyTree());
  second->build(rodTemplate, tmspecso, 0);
  second->rotateZ(rodPhiRotation);
  rods_.push_back(second);

  for (int i = 2; i < numRods(); i++) {
    RodPair* rod = i%2 ? GeometryFactory::clone(*second) : GeometryFactory::clone(*first); // clone rods
    rod->myid(i+1);
    rod->rotateZ(rodPhiRotation*(i%2 ? i-1 : i));
    rods_.push_back(rod);
    }

  // computing the layer's place radius as the average of all the modules' radii
  placeRadius_  = std::accumulate(tmspecsi.begin(), tmspecsi.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
  placeRadius_ += std::accumulate(tmspecso.begin(), tmspecso.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
  placeRadius_ /= tmspecsi.size() + tmspecso.size();

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

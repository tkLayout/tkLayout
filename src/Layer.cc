#include "Layer.hh"
#include "RodPair.hh"
#include "MessageLogger.hh"
#include "ConversionStation.hh"

FlatRingsGeometryInfo::FlatRingsGeometryInfo() {}

void FlatRingsGeometryInfo::calculateFlatRingsGeometryInfo(std::vector<StraightRodPair*> flatPartRods, double bigParity) {
 
  StraightRodPair* minusBigDeltaRod = (bigParity > 0 ? flatPartRods.at(1) : flatPartRods.front());
  const auto& minusBigDeltaModules = minusBigDeltaRod->modules().first;
  StraightRodPair* plusBigDeltaRod = (bigParity > 0 ? flatPartRods.front() : flatPartRods.at(1));
  const auto& plusBigDeltaModules = plusBigDeltaRod->modules().first;  

  int i = 0;
  double rStartInner;
  double zStartInner_REAL;
  double rEndInner;
  double zEndInner_REAL;
  int smallParity = minusBigDeltaRod->zPlusParity();
  for (const auto& m : minusBigDeltaModules) {
    if (i > 0) {
      zStartInner_REAL = m.planarMinZ();
      // Special case where ring has been built going upwards, and with zEndInner_REAL < zError
      if ((zStartInner_REAL < zEndInner_REAL) && (smallParity > 0)) {
	rEndInner -= m.dsDistance(); 
	rStartInner = m.center().Rho() + 0.5 * m.dsDistance();
      }
	else { // Standard case
	rStartInner = m.center().Rho() - 0.5 * m.dsDistance();
      }

      if (rStartInner != rEndInner) {
	double fact = (((rStartInner - rEndInner) > 0) ? 1. : -1.);
	
	if (zStartInner_REAL != zEndInner_REAL) {
	  // Let's call (H1pUP, H1ppDOWN) the line that binds H1pUP of the previous module, with H1ppDOWN of the next module.
	  // Standard Case : (H1pUP, H1ppDOWN) is secant with (Z).
	  // Let's call P the intersection point.
	  // zError is defined as the Z of point P.
	  double zErrorInnerAngle = atan( (rStartInner - rEndInner) / (zStartInner_REAL - zEndInner_REAL) );
	  zErrorInner_[i] = fact * (zStartInner_REAL - rStartInner / tan(zErrorInnerAngle));
	}	
	else { // Case where (H1pUP, H1ppDOWN) is orthogonal to (Z).
	  zErrorInner_[i] = fact * zStartInner_REAL;
	}
      }
      else { // Case where consecutive modules are placed at the same r within a rod, for the Pixel for example.
	if (zStartInner_REAL != zEndInner_REAL) { // If consecutive modules (along a rod) are not touching, zError is undefined.
	  zErrorInner_[i] = std::numeric_limits<double>::quiet_NaN();
	}
	else { // If consecutive modules (along a rod) are touching, zError is infinity.
	  zErrorInner_[i] = std::numeric_limits<double>::infinity();
	}
      }

    }

    zEndInner_REAL = m.planarMaxZ();
    rEndInner = m.center().Rho() + 0.5 * m.dsDistance();
   
    i++;
    smallParity = -smallParity;
  }

  i = 0;
  double rStartOuter;
  double zStartOuter_REAL;
  double rEndOuter;
  double zEndOuter_REAL;
  smallParity = plusBigDeltaRod->zPlusParity();
  for (const auto& m : plusBigDeltaModules) {
    if (i > 0) {
      zStartOuter_REAL = m.planarMinZ();
      if ((zStartOuter_REAL < zEndOuter_REAL) && (smallParity > 0)) {
	  rEndOuter -= m.dsDistance(); 
	  rStartOuter = m.center().Rho() + 0.5 * m.dsDistance();
	}
	else {
	  rStartOuter = m.center().Rho() - 0.5 * m.dsDistance();
	}

      if (rStartOuter != rEndOuter) {
	double fact = (((rStartOuter - rEndOuter) > 0) ? 1. : -1.);
	if (zStartOuter_REAL != zEndOuter_REAL) {
	  double zErrorOuterAngle = atan( (rStartOuter - rEndOuter) / (zStartOuter_REAL - zEndOuter_REAL) );
	  zErrorOuter_[i] = fact * (zStartOuter_REAL - rStartOuter / tan(zErrorOuterAngle));
	}
	else {
	  zErrorOuter_[i] = fact * zStartOuter_REAL;
	}
      }
      else {
	if (zStartOuter_REAL != zEndOuter_REAL) {
	  zErrorOuter_[i] = std::numeric_limits<double>::quiet_NaN();
	}
	else {
	  zErrorOuter_[i] = std::numeric_limits<double>::infinity();
	}
      }

    }
 
    zEndOuter_REAL = m.planarMaxZ();
    rEndOuter = m.center().Rho() + 0.5 * m.dsDistance();

    i++;
    smallParity = -smallParity;
  }  
}

Layer::TiltedRingsGeometryInfo::TiltedRingsGeometryInfo(int numModulesFlat, double flatPartrEndInner, double flatPartrEndOuter, double flatPartzEnd,  double flatPartzEnd_REAL, TiltedRingsTemplate tiltedRingsGeometry) {
  for (int i = (numModulesFlat + 1); i < (numModulesFlat + tiltedRingsGeometry.size() + 1); i++) {

    if (i == (numModulesFlat + 1)) {
      deltaZInner_[i] = tiltedRingsGeometry[i]->zInner() - flatPartzEnd;
      deltaZOuter_[i] = tiltedRingsGeometry[i]->zOuter() - flatPartzEnd;

      double zErrorInnerAngle = atan( (tiltedRingsGeometry[i]->rStartInner_REAL() - flatPartrEndInner) / (tiltedRingsGeometry[i]->zStartInner_REAL() - flatPartzEnd_REAL) );
      zErrorInner_[i] = tiltedRingsGeometry[i]->zStartInner_REAL() - tiltedRingsGeometry[i]->rStartInner_REAL() / tan(zErrorInnerAngle);

      double zErrorOuterAngle = atan( (tiltedRingsGeometry[i]->rStartOuter_REAL() - flatPartrEndOuter) / (tiltedRingsGeometry[i]->zStartOuter_REAL() - flatPartzEnd_REAL) );
      zErrorOuter_[i] = tiltedRingsGeometry[i]->zStartOuter_REAL() - tiltedRingsGeometry[i]->rStartOuter_REAL() / tan(zErrorOuterAngle);

    }

    else {
      deltaZInner_[i] = tiltedRingsGeometry[i]->zInner() - tiltedRingsGeometry[i-1]->zInner();
      deltaZOuter_[i] = tiltedRingsGeometry[i]->zOuter() - tiltedRingsGeometry[i-1]->zOuter();

      //covInner_[i] = (tiltedRingsGeometry[i]->thetaStartInner() - tiltedRingsGeometry[i-1]->thetaEndInner());

      double zErrorInnerAngle = atan( (tiltedRingsGeometry[i]->rStartInner_REAL() - tiltedRingsGeometry[i-1]->rEndInner_REAL()) / (tiltedRingsGeometry[i]->zStartInner_REAL() - tiltedRingsGeometry[i-1]->zEndInner_REAL()) );
      zErrorInner_[i] = tiltedRingsGeometry[i]->zStartInner_REAL() - tiltedRingsGeometry[i]->rStartInner_REAL() / tan(zErrorInnerAngle);

      double zErrorOuterAngle = atan( (tiltedRingsGeometry[i]->rStartOuter_REAL() - tiltedRingsGeometry[i-1]->rEndOuter_REAL()) / (tiltedRingsGeometry[i]->zStartOuter_REAL() - tiltedRingsGeometry[i-1]->zEndOuter_REAL()) );
      zErrorOuter_[i] = tiltedRingsGeometry[i]->zStartOuter_REAL() - tiltedRingsGeometry[i]->rStartOuter_REAL() / tan(zErrorOuterAngle);
    }
  }
}

void Layer::check() {
  PropertyObject::check();

  //if (buildNumModules() > 0 && maxZ.state()) throw PathfulException("Only one between numModules and maxZ can be specified");
  //else if (buildNumModules() == 0 && !maxZ.state()) throw PathfulException("At least one between numModules and maxZ must be specified");


  if (!isTilted()) {
    if (buildNumModules() > 0 && maxZ.state()) throw PathfulException("Only one between numModules and maxZ can be specified");
    if (buildNumModules() == 0 && !maxZ.state()) throw PathfulException("At least one between numModules and maxZ must be specified");
    if (!phiOverlap.state()) throw PathfulException("Flat layer : phiOverlap must be specified.");
    if (!phiSegments.state()) throw PathfulException("Flat layer : phiSegments must be specified.");
    if (numRods.state()) throw PathfulException("Flat layer : numRods should not be specified.");
    if (isTiltedAuto.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify isTiltedAuto. Not used.");
  }

  if (!isTilted() || (isTilted() && !isTiltedAuto())) {
    if (buildNumModulesFlat.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify numModulesFlat. Not used.");
    if (buildNumModulesTilted.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify numModulesTilted. Not used.");
  }

  if (isTilted()) {
    if (maxZ.state()) logERROR("Tilted layer : maxZ was specified. Routing of services will be forced to be at Z = maxZ.");
    if (!isTiltedAuto.state()) throw PathfulException("Tilted layer : isTiltedAuto must be specified.");
    if (phiOverlap.state()) throw PathfulException("Tilted layer : phiOverlap should not be specified.");
    if (phiSegments.state()) throw PathfulException("Tilted layer : phiSegments should not be specified.");
    if (isTiltedAuto()) {    
      if (!buildNumModulesFlat.state()) throw PathfulException("Tilted layer with automatic placement : numModulesFlat must be specified.");
      if (!buildNumModulesTilted.state()) throw PathfulException("Tilted layer with automatic placement : numModulesTilted must be specified.");
      if (buildNumModules() > 0 && buildNumModulesFlat.state() && buildNumModulesTilted.state()) {
	if (buildNumModules() !=  (buildNumModulesFlat() + buildNumModulesTilted())) {
	  throw PathfulException("Tilted layer : numModules != (numModulesFlat + numModulesTilted). Anyway, for automatic placement, it is not needed to specify numModules. Please specify numModulesFlat and numModulesTilted only.");
	}
      }
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
  //std::cout << "rodTemplate.size() = " << rodTemplate.size() << std::endl;
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
    else {
      lastThetaEnd = tiltedRingsGeometry[i-1]->thetaEndOuter_REAL();
    }
 
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

  // FIRST ROD : assign common properties
  StraightRodPair* first = GeometryFactory::make<StraightRodPair>();
  first->myid(1);
  first->minBuildRadius(minBuildRadius()-bigDelta());
  first->maxBuildRadius(maxBuildRadius()+bigDelta());
  if (buildNumModules() > 0) first->buildNumModules(buildNumModules());
  else if (maxZ.state()) first->maxZ(maxZ());
  first->smallDelta(smallDelta());
  //first->ringNode = ringNode; // we need to pass on the contents of the ringNode to allow the RodPair to build the module decorators
  if (isFlatPart) { first->isFlatPart(true); first->zPlusParity( pow(-1, buildNumModulesFlat()) ); }
  first->store(propertyTree());
  // SECOND ROD : copy first rod
  logINFO(Form("Copying rod %s", fullid(*this).c_str()));
  StraightRodPair* second = GeometryFactory::clone(*first);
  second->myid(2);


  // FIRST ROD : build and store
  bool isPlusBigDeltaRod = (bigParity() > 0);
  first->isOuterRadiusRod(isPlusBigDeltaRod);
  first->build(rodTemplate, isPlusBigDeltaRod);
  first->translateR(placeRadius_ + (isPlusBigDeltaRod ? bigDelta() : -bigDelta()));
  if (!isFlatPart) { rods_.push_back(first); buildNumModulesFlat(first->numModulesSide(1)); }
  else { flatPartRods_.push_back(first); }

  // SECOND ROD : assign other properties, build and store 
  if (!sameParityRods()) second->zPlusParity(first->zPlusParity()*-1);
  isPlusBigDeltaRod = (bigParity() < 0);
  second->isOuterRadiusRod(isPlusBigDeltaRod);
  second->build(rodTemplate, isPlusBigDeltaRod);
  second->translateR(placeRadius_ + (isPlusBigDeltaRod ? bigDelta() : -bigDelta()));
  second->rotateZ(rodPhiRotation);
  if (!isFlatPart) { rods_.push_back(second); }
  else { flatPartRods_.push_back(second); }

  // All other Rods
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
    int numModulesFlat = 0;
    int numModulesTilted = 0;
    while(getline(ifs, line).good()) {
      if (line.empty()) continue;
      auto tokens = split<double>(line, " ", false);
      if (tokens.size() < 7) { logERROR("Failed parsing tilted barrel line: " + line); continue; };
      TiltedModuleSpecs ti{ tokens[0], tokens[1], tokens[2]*M_PI/180. };
      TiltedModuleSpecs to{ tokens[3], tokens[4], tokens[5]*M_PI/180. };
      if (ti.valid()) tmspecsi.push_back(ti);
      if (to.valid()) tmspecso.push_back(to);
      if (ti.valid() || to.valid()) {
	if (tokens[2] == 0. && tokens[5] == 0.) numModulesFlat++;
	else numModulesTilted++;
      }
      numRods(tokens[6]); // this assumes every row of the spec file has the same value for the last column (num rods in phi) 
    }
    buildNumModulesFlat(numModulesFlat);
    buildNumModulesTilted(numModulesTilted);
    buildNumModules(numModulesFlat + numModulesTilted);
    ifs.close();
  }

  else {

    double flatPartThetaEnd = M_PI / 2.;
    double flatPartrEndInner = 0;
    double flatPartrEndOuter = 0;
    double flatPartzEnd = 0;
    double flatPartzEnd_REAL = 0;

    double flatPartrInnerSmall = std::numeric_limits<double>::max();
    double flatPartrOuterSmall = std::numeric_limits<double>::max();
    double flatPartrInnerBig = 0;
    double flatPartrOuterBig = 0;

    if (buildNumModulesFlat() != 0) {

      buildNumModules(buildNumModulesFlat());
      buildStraight(true);


      if (flatPartRods_.size() >= 2) {
	StraightRodPair* flatPartRod1 = flatPartRods_.front();
	const auto& zPlusModules1 = flatPartRod1->modules().first;
	for (const auto& m : zPlusModules1) {
	  TiltedModuleSpecs t1{m.center().Rho(), m.center().Z(), 0.0};
	  if (t1.valid()) (bigParity() > 0 ? tmspecso.push_back(t1) : tmspecsi.push_back(t1));
	  (bigParity() > 0 ? flatPartrOuterSmall = MIN(flatPartrOuterSmall, m.center().Rho() + 0.5*m.dsDistance()) : flatPartrInnerSmall = MIN(flatPartrInnerSmall, m.center().Rho() + 0.5*m.dsDistance()));
	  (bigParity() > 0 ? flatPartrOuterBig = MAX(flatPartrOuterBig, m.center().Rho() + 0.5*m.dsDistance()) : flatPartrInnerBig = MAX(flatPartrInnerBig, m.center().Rho() + 0.5*m.dsDistance()));
	}

	StraightRodPair* flatPartRod2 = flatPartRods_.at(1);
	const auto& zPlusModules2 = flatPartRod2->modules().first;
	for (const auto& m : zPlusModules2) {
	  TiltedModuleSpecs t2{m.center().Rho(), m.center().Z(), 0.0};
	  if (t2.valid()) (bigParity() > 0 ? tmspecsi.push_back(t2) : tmspecso.push_back(t2));
	  (bigParity() > 0 ? flatPartrInnerSmall = MIN(flatPartrInnerSmall, m.center().Rho() + 0.5*m.dsDistance()) : flatPartrOuterSmall = MIN(flatPartrOuterSmall, m.center().Rho() + 0.5*m.dsDistance()));
	  (bigParity() > 0 ? flatPartrInnerBig = MAX(flatPartrInnerBig, m.center().Rho() + 0.5*m.dsDistance()) : flatPartrOuterBig = MAX(flatPartrOuterBig, m.center().Rho() + 0.5*m.dsDistance()));
	}

	flatPartThetaEnd = (bigParity() > 0 ? flatPartRod1->thetaEnd_REAL() : flatPartRod2->thetaEnd_REAL());
	auto lastMod1 = zPlusModules1.back();
	auto lastMod2 = zPlusModules2.back();	
	flatPartrEndInner = (bigParity() > 0 ? lastMod2.center().Rho() + 0.5* lastMod2.dsDistance() : lastMod1.center().Rho() + 0.5* lastMod1.dsDistance());
	flatPartrEndOuter = (bigParity() > 0 ? lastMod1.center().Rho() + 0.5* lastMod1.dsDistance() : lastMod2.center().Rho() + 0.5* lastMod2.dsDistance());
	flatPartzEnd = (bigParity() > 0 ? lastMod2.center().Z() : lastMod1.center().Z());	
	flatPartzEnd_REAL = (bigParity() > 0 ? lastMod1.planarMaxZ() : lastMod2.planarMaxZ());	//TAKE CAREEEEEE : REMOVE FLAT PART OVERLAP ?


	RectangularModule* flatPartrmod = GeometryFactory::make<RectangularModule>();
	flatPartrmod->store(propertyTree());
	if (ringNode.count(1) > 0) flatPartrmod->store(ringNode.at(1));
	flatPartrmod->build();
	double width = flatPartrmod->width();
	double dsDistance = flatPartrmod->dsDistance();
	double T = tan(2.*M_PI / numRods());
	double A = 1. / (2. * flatPartrInnerSmall);
	double B = 1. / (2. * flatPartrOuterSmall);
	double a = T * A * B;
	double b = - (A + B);
	double c = - T;
	double s = (-b - sqrt(b*b - 4*a*c))/(2*a);
	flatPartPhiOverlapSmallDeltaMinus_ = width + s;
	A = 1. / (2. * flatPartrInnerBig);
	B = 1. / (2. * flatPartrOuterBig);
	a = T * A * B;
	b = - (A + B);
	c = - T;
	s = (-b - sqrt(b*b - 4*a*c))/(2*a);
	flatPartPhiOverlapSmallDeltaPlus_ = width + s;

	flatPartAverageR_ = (flatPartrInnerSmall + flatPartrInnerBig + flatPartrOuterSmall + flatPartrOuterBig) / 4. - 0.5 * dsDistance;

	flatRingsGeometryInfo_.calculateFlatRingsGeometryInfo(flatPartRods_, bigParity());


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
	/*std::cout << "ringNumber = " << ringNumber << std::endl;
	  std::cout << "tiltAngle = " << tiltedRingsGeometry_[ringNumber]->tiltAngle()*M_PI/180. << std::endl;
	  std::cout << "innerRadius = " << tiltedRingsGeometry_[ringNumber]->innerRadius() << std::endl;
	  std::cout << "zInner = " << tiltedRingsGeometry_[ringNumber]->zInner() << std::endl;
	  std::cout << "outerRadius = " << tiltedRingsGeometry_[ringNumber]->outerRadius() << std::endl;
	  std::cout << "zOuter = " << tiltedRingsGeometry_[ringNumber]->zOuter() << std::endl;*/




	/*std::cout << "theta2 = " << tiltedRingsGeometry_[ringNumber]->thetaOuter() * 180. / M_PI << std::endl;
	  std::cout << "idealTilt2 = " << tiltedRingsGeometry_[ringNumber]->tiltAngleIdealOuter() << std::endl;
	  std::cout << "gap = " << tiltedRingsGeometry_[ringNumber]->gapR() << std::endl;
	  std::cout << "avR = " << tiltedRingsGeometry_[ringNumber]->averageR() << std::endl;
	  if (i >= 1) { std::cout << "cov1 = " << (tiltedRingsGeometry_[ringNumber]->thetaStartInner() - tiltedRingsGeometry_[ringNumber-1]->thetaEndInner()) * 180. / M_PI << std::endl; }
	  if (i >= 1) { std::cout << "deltaz2 = " << tiltedRingsGeometry_[ringNumber]->zOuter() - tiltedRingsGeometry_[i-1]->zOuter() << std::endl; }

	  if (i >= 1) {
	  double zErrorAngle = atan( (tiltedRingsGeometry_[ringNumber]->rStartOuter_REAL() - tiltedRingsGeometry_[i-1]->rEndOuter_REAL()) / (tiltedRingsGeometry_[i-1]->zEndOuter_REAL() - tiltedRingsGeometry_[ringNumber]->zStartOuter_REAL()) );
	  std::cout << "zError = " << tiltedRingsGeometry_[ringNumber]->zStartOuter_REAL() + tiltedRingsGeometry_[ringNumber]->rStartOuter_REAL() / tan(zErrorAngle) << std::endl; 
	  }

	  std::cout << "cov2 = " << atan(tiltedRingsGeometry_[ringNumber]->rEndOuter_REAL() / tiltedRingsGeometry_[ringNumber]->zEndOuter_REAL()) * 180. / M_PI << std::endl;*/

	if (ti.valid()) tmspecsi.push_back(ti);
	if (to.valid()) tmspecso.push_back(to);
      }
      tiltedRingsGeometryInfo_ = TiltedRingsGeometryInfo(buildNumModulesFlat(), flatPartrEndInner, flatPartrEndOuter, flatPartzEnd,  flatPartzEnd_REAL, tiltedRingsGeometry_);
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
  first->isOuterRadiusRod(false);
  first->store(propertyTree());
  first->build(rodTemplate, tmspecsi, 1);
  rods_.push_back(first);

  TiltedRodPair* second = GeometryFactory::make<TiltedRodPair>();
  second->myid(2);
  second->isOuterRadiusRod(true);
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

    if (!isTilted()) {
      buildStraight(false);
      if (buildNumModules() > 0 ) buildNumModulesFlat(buildNumModules());
      buildNumModulesTilted(0);
    }
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

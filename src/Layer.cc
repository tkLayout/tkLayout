#include "Layer.hh"
#include "RodPair.hh"
#include "MessageLogger.hh"
#include "ConversionStation.hh"


define_enum_strings(Layer::RadiusMode) = { "shrink", "enlarge", "fixed", "auto" };


void FlatRingsGeometryInfo::calculateFlatRingsGeometryInfo(std::vector<StraightRodPair*> flatPartRods, double bigParity) {
 
  StraightRodPair* minusBigDeltaRod = (bigParity > 0 ? flatPartRods.at(1) : flatPartRods.front());
  const auto& minusBigDeltaModules = minusBigDeltaRod->modules().first;
  StraightRodPair* plusBigDeltaRod = (bigParity > 0 ? flatPartRods.front() : flatPartRods.at(1));
  const auto& plusBigDeltaModules = plusBigDeltaRod->modules().first;  

  int i = 0;
  double rStartInner;
  double zStartInner_REAL;
  double rEndInner = 0.; // It gets set at the first iteration
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
  double rEndOuter = 0.; // It gets set at the first iteration
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


TiltedRingsGeometryInfo::TiltedRingsGeometryInfo(int numModulesFlat, double flatPartrEndInner, double flatPartrEndOuter, double flatPartzEnd,  double flatPartzEnd_REAL, TiltedRingsTemplate tiltedRingsGeometry) {
  const int numTiltedRings = tiltedRingsGeometry.size();
  for (int i = (numModulesFlat + 1); i < (numModulesFlat + numTiltedRings + 1); i++) {

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


/* PUBLIC */

/*
 * Check that properties assignment (from cfg files) make some sense.
 */
void Layer::check() {
  PropertyObject::check();

  // UNTILTED LAYER
  if (!isTilted() && !isSkewedForInstallation()) {
    if (buildNumModules() > 0 && maxZ.state()) throw PathfulException("Only one between numModules and maxZ can be specified");
    if (buildNumModules() == 0 && !maxZ.state()) throw PathfulException("At least one between numModules and maxZ must be specified");
    if (numRods.state() && (phiOverlap.state() || phiSegments.state())) throw PathfulException("Flat layer : Only one between numRods  and (phiOverlap + phiSegments) can be specified.");
    if (!numRods.state() && !phiOverlap.state()) throw PathfulException("Flat layer : phiOverlap must be specified.");
    if (!numRods.state() && !phiSegments.state()) throw PathfulException("Flat layer : phiSegments must be specified.");
    if (phiForbiddenRanges.state() && rotateLayerByRodsDeltaPhiHalf()) throw PathfulException("Flat layer : Only one between phiForbiddenRanges and rotateLayerByRodsDeltaPhiHalf must be specified.");
    if (isTiltedAuto.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify isTiltedAuto. Not used.");
  }

  if (!isTilted() || (isTilted() && !isTiltedAuto())) {
    if (buildNumModulesFlat.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify numModulesFlat. Not used.");
    if (buildNumModulesTilted.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify numModulesTilted. Not used.");
  }

  // TILTED LAYER
  if (isTilted()) {
    if (maxZ.state()) logERROR("Tilted layer : maxZ was specified. Routing of services will be forced to be at Z = maxZ.");
    if (!isTiltedAuto.state()) throw PathfulException("Tilted layer : isTiltedAuto must be specified.");
    if (phiOverlap.state()) throw PathfulException("Tilted layer : phiOverlap should not be specified.");
    if (phiSegments.state()) throw PathfulException("Tilted layer : phiSegments should not be specified.");
    // TILTED AUTO PLACEMENT MODE
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

  // SKEWED LAYER
  if (isSkewedForInstallation()) {
    if (!skewedModuleEdgeShift.state()) throw PathfulException("Skewed layer mode: skewedModuleEdgeShift must be specified.");
    if (!numRods.state()) throw PathfulException("Skewed layer mode: numRods must be specified.");
    if (radiusMode() != RadiusMode::FIXED) throw PathfulException("Skewed layer mode: the (average) radii of layers must be specified.");
    if (isTilted()) throw PathfulException("A layer was set to both skewed and tilted: this is not presently supported.");
    if (phiForbiddenRanges.state()) throw PathfulException("Skewed layer mode: phiForbiddenRange is not supported.");
    if (rotateLayerByRodsDeltaPhiHalf()) throw PathfulException("Skewed layer mode: rotateLayerByRodsDeltaPhiHalf is not supported.");
    if (phiOverlap.state()) throw PathfulException("Skewed layer mode: phiOverlap should not be specified.");
    if (phiSegments.state()) throw PathfulException("Skewed layer mode: phiSegments should not be specified.");
  }
}


/*
 * Build a layer.
 */
void Layer::build() {
  ConversionStation* conversionStation;

  try { 
    materialObject_.store(propertyTree());
    materialObject_.build();

    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    // UNTILTED LAYER
    if (!isTilted()) {
      buildStraight();
      if (buildNumModules() > 0 ) buildNumModulesFlat(buildNumModules());
      buildNumModulesTilted(0);
    }
    // TILTED LAYER
    else buildTilted();

    for (auto& currentStationNode : stationsNode) {
      conversionStation = new ConversionStation(subdetectorName());
      conversionStation->store(propertyTree());
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


void Layer::cutAtEta(double eta) {
  for (auto& r : rods_) r.cutAtEta(eta); 
  rods_.erase_if([](const RodPair& r) { return r.numModules() == 0; }); // get rid of rods which have been completely pruned
}


/* PRIVATE */

/* STRAIGHT LAYER */

/*
 * Build a straight layer (with opposition to tilted layer).
 * The ladders can be skewed around the (CMS_Z) axis if requested.
 */
void Layer::buildStraight() {

  // COMPUTES A ROD TEMPLATE
  RodTemplate rodTemplate = makeRodTemplate();

  // MID-RADIUS PLACEMENT AND NUMRODS
  computePlaceRadiiAndNumRods(rodTemplate);  

  // NON-SKEWED LAYER: compute shift in Phi between consecutive rods centers.
  const double rodCenterPhiShift = computeRodCenterPhiShift();

  // SKEWED LAYER: compute parameters of interest.
  const SkewedLayerPhiShifts& phiShifts = (isSkewedForInstallation() ? buildSkewed() : SkewedLayerPhiShifts{ 0.,0.,0.} );
  const double installationMinusBigDeltaRodCenterPhiShift = phiShifts.installationMinusBigDeltaRodCenterPhiShift;
  const double commonRodCenterPhiShift = phiShifts.commonRodCenterPhiShift;
  const double skewedRodCenterPhiShift = phiShifts.skewedRodCenterPhiShift;


  // FIRST ROD : assign common properties
  StraightRodPair* firstRod = GeometryFactory::make<StraightRodPair>(subdetectorName());
  assignRodCommonProperties(firstRod);
  
  // CLONES FROM FIRST ROD
  logINFO(Form("Copying rod %s", fullid(*this).c_str()));
  // second rod case
  StraightRodPair* secondRod = GeometryFactory::clone(*firstRod);
  // skewed rod case
  StraightRodPair* skewedRod = GeometryFactory::clone(*firstRod);

  // FIRST ROD : place and store
  const bool isFirstRodAtPlusBigDelta = placeAndStoreFirstRod(firstRod, rodTemplate,
							      rodCenterPhiShift, installationMinusBigDeltaRodCenterPhiShift);

  // SECOND ROD : assign other properties, place and store
  const double firstRodCenterPhi = firstRod->Phi();
  placeAndStoreSecondRod(secondRod, rodTemplate, 
			 isFirstRodAtPlusBigDelta, firstRod->zPlusParity(), firstRodCenterPhi, 
			 rodCenterPhiShift, commonRodCenterPhiShift);

  // NON-SKEWED INSTALLATION MODE
  if (!isSkewedForInstallation()) { 
    buildAndStoreClonedRodsInNonSkewedMode(firstRod, secondRod,
					   rodCenterPhiShift);
  }

  // SKEWED INSTALLATION MODE
  else {
    buildAndStoreClonedRodsInSkewedMode(firstRod, secondRod, skewedRod,
					commonRodCenterPhiShift, skewedRodCenterPhiShift);
  }

}


/* generic optimizations methods */

/*
 * Create a template untilted rod (can be skewed).
 */
RodTemplate Layer::makeRodTemplate(const double skewAngle) {
  RodTemplate rodTemplate(buildNumModules() > 0 ? buildNumModules() : (!ringNode.empty() ? ringNode.rbegin()->first + 1 : 1)); // + 1 to make room for a default constructed module to use when building rods in case the rodTemplate vector doesn't have enough elements
  const int numModules = rodTemplate.size();
  for (int i = 0; i < numModules; i++) {
    rodTemplate[i] = std::move(unique_ptr<BarrelModule>(GeometryFactory::make<BarrelModule>(GeometryFactory::make<RectangularModule>(), subdetectorName())));
    rodTemplate[i]->store(propertyTree());
    if (ringNode.count(i+1) > 0) rodTemplate[i]->store(ringNode.at(i+1));
    if (isSkewedForInstallation()) rodTemplate[i]->skewAngle(skewAngle);
    rodTemplate[i]->build();
  }
  return rodTemplate;
}


/*
 * Compute the layer mid-radii and numRods if needed.
 */
void Layer::computePlaceRadiiAndNumRods(const RodTemplate& rodTemplate) {

  // Is it a non-tilted layer, or the flat part of a tilted layer?
  // Indeed, Layer::buildStraight() is also used to build the flat part of a tilted layer.
  const bool isFlatPart = isTilted();  

  // FLAT LAYER
  if (!isFlatPart) {
    // if numRods is not specified, compute it from coverage info.
    if (!numRods.state()) {
      std::pair<double, int> optimalLayerParms = calculateOptimalLayerParms(rodTemplate);
      placeRadius_ = optimalLayerParms.first; 
      numRods(optimalLayerParms.second);
    }
    // otherwise, just place the layer at the specified radius.
    else { placeRadius_ = placeRadiusHint(); }
  }

  // FLAT PART OF A TILTED LAYER
  else { placeRadius_ = placeRadiusHint(); } // place the layer at the specified radius.

  if (!minBuildRadius.state() || !maxBuildRadius.state()) {
    minBuildRadius(placeRadius_);
    maxBuildRadius(placeRadius_);
  }
}


/*
 * Straight layer: compute parameters of interest.
 */
std::pair<float, int> Layer::calculateOptimalLayerParms(const RodTemplate& rodTemplate) {
                                                              
  // CUIDADO fix placeRadiusHint!!!!!
  double maxDsDistance = (*std::max_element(rodTemplate.begin(), 
                                            rodTemplate.end(), 
                                            [](const unique_ptr<BarrelModule>& m1, const unique_ptr<BarrelModule>& m2) { return m1->dsDistance() > m2->dsDistance(); } ))->dsDistance();
  float moduleWidth = (*rodTemplate.rbegin())->minWidth();
  float f = moduleWidth/2 - phiOverlap()/2;
  float gamma = atan(f/(placeRadiusHint() + bigDelta() + smallDelta() + maxDsDistance/2)) + atan(f/(placeRadiusHint() - bigDelta() + smallDelta() + maxDsDistance/2));

  float tentativeModsPerSegment = 2. * M_PI / (gamma * phiSegments());

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


/*
 * Straight layer: compute mean radius from other available parameters.
 */
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


/*
 * Non-Skewed mode:
 * Compute shift in Phi between consecutive rods centers.
 */
const double Layer::computeRodCenterPhiShift() const {
  double rodCenterPhiShift = 0;

  if (!phiForbiddenRanges.state()) { rodCenterPhiShift = 2. * M_PI / numRods(); }

  // PHI FORBIDDEN RANGE: forbid rod placement in phi, in a given phi range.
  else {
    const double forbiddenPhiLowerA = phiForbiddenRanges.at(0) * M_PI / 180.;
    const double forbiddenPhiUpperA = phiForbiddenRanges.at(1) * M_PI / 180.;
    const double forbiddenPhiLowerB = phiForbiddenRanges.at(2) * M_PI / 180.;
    const double forbiddenPhiUpperB = phiForbiddenRanges.at(3) * M_PI / 180.;
    const double forbiddenPhiCircumference = forbiddenPhiUpperA - forbiddenPhiLowerA + forbiddenPhiUpperB - forbiddenPhiLowerB;
    rodCenterPhiShift = (2*M_PI - forbiddenPhiCircumference) / numRods();
  }

  return rodCenterPhiShift;
}


/*
 * Assign rod common properties & property tree (use before cloning).
 */
void Layer::assignRodCommonProperties(StraightRodPair* rod) const {
  rod->myid(1);
  rod->minBuildRadius(minBuildRadius() - bigDelta());
  rod->maxBuildRadius(maxBuildRadius() + bigDelta());
  if (buildNumModules() > 0) rod->buildNumModules(buildNumModules());
  else if (maxZ.state()) rod->maxZ(maxZ());
  rod->smallDelta(smallDelta());
  const bool isFlatPart = isTilted();
  if (isFlatPart) { rod->isFlatPart(true); rod->zPlusParity( pow(-1, buildNumModulesFlat()) ); }
  //rod->ringNode = ringNode; // we need to pass on the contents of the ringNode to allow the RodPair to build the module decorators
  rod->store(propertyTree());
}


/* not specific to skewed mode */

/*
 * Place and store the first rod.
 */
const bool Layer::placeAndStoreFirstRod(StraightRodPair* firstRod, const RodTemplate& rodTemplate,
					const double rodCenterPhiShift, const double installationMinusBigDeltaRodCenterPhiShift) {

  // COMPUTE PHI POSITION
  // Starts at phi = 0 (Warning: full barrel can have a rotation, which will add up to these).
  double firstRodCenterPhi = 0.;

  // One (and only one) of the following special cases can be present.
  // Special case A: 
  // If requested in cfg, shift the layer by rodsDeltaPhi / 2
  if (rotateLayerByRodsDeltaPhiHalf()) { firstRodCenterPhi = rodCenterPhiShift / 2.; }

  // Special case B: phi forbidden range.
  // Forbid rod placement in phi, in a given phi range. Starts placement at specified value.
  else if (phiForbiddenRanges.state()) {
    const double forbiddenPhiUpperA = phiForbiddenRanges.at(1) * M_PI / 180.;
    firstRodCenterPhi = (forbiddenPhiUpperA + rodCenterPhiShift / 2.);
    // NB : should not be rodCenterPhiShift / 2. (deltaPhi between consecutive rods) 
    // but rod phi aperture (phi covered by a rod).
  }
  
  // Special case C: skewed mode.
  // The first rod has a specific shift in Phi.
  else if (isSkewedForInstallation()) { firstRodCenterPhi = installationMinusBigDeltaRodCenterPhiShift; }


  // PLACE AND STORE
  const bool isFirstRodAtPlusBigDelta = (!isSkewedForInstallation() ? (bigParity() > 0) : false);
  placeAndStoreRod(firstRod, rodTemplate, isFirstRodAtPlusBigDelta, firstRodCenterPhi);
  const bool isFlatPart = isTilted();
  if (!isFlatPart) { buildNumModulesFlat(firstRod->numModulesSide(1)); }

  return isFirstRodAtPlusBigDelta;
}


/*
 * Place and store the second rod.
 */
void Layer::placeAndStoreSecondRod(StraightRodPair* secondRod, const RodTemplate& rodTemplate, 
				   const bool isFirstRodAtPlusBigDelta, const int firstRodZPlusParity, const double firstRodCenterPhi, 
				   const double rodCenterPhiShift, const double commonRodCenterPhiShift) {

  secondRod->myid(2);
  if (!sameParityRods()) secondRod->zPlusParity(-1 * firstRodZPlusParity);
  
  // COMPUTE PHI POSITION
  double secondRodCenterPhi = firstRodCenterPhi;
  if (!isSkewedForInstallation()) { secondRodCenterPhi += rodCenterPhiShift; }
  // skewed mode: the second rod has a specific shift in Phi.
  else { secondRodCenterPhi += commonRodCenterPhiShift; }

  // PLACE AND STORE
  const bool isSecondRodAtPlusBigDelta = !isFirstRodAtPlusBigDelta;
  placeAndStoreRod(secondRod, rodTemplate, isSecondRodAtPlusBigDelta, secondRodCenterPhi);

}


/*
 * Place a rod radially and in phi + Store it.
 */
void Layer::placeAndStoreRod(StraightRodPair* rod, const RodTemplate& rodTemplate, const bool isPlusBigDeltaRod, const double rodCenterPhi) {
  rod->isOuterRadiusRod(isPlusBigDeltaRod);
  rod->build(rodTemplate, isPlusBigDeltaRod);
  rod->translateR(placeRadius_ + (isPlusBigDeltaRod ? bigDelta() : -bigDelta()));
  
  rod->rotateZ(rodCenterPhi);
  if (!isSkewedForInstallation()) {
    const bool isFlatPart = isTilted();
    if (!isFlatPart) { rods_.push_back(rod); buildNumModulesFlat(rod->numModulesSide(1)); }
    else { flatPartRods_.push_back(rod); }
  }
}


/*
 * Non-Skewed mode:
 * General method to build & store all other rods (apart from first 2 rods).
 */
void Layer::buildAndStoreClonedRodsInNonSkewedMode(const StraightRodPair* firstRod, const StraightRodPair* secondRod,
						   const double rodCenterPhiShift) {
  // ALL OTHER RODS
  for (int i = 2; i < numRods(); i++) {

    // CLONE RODS
    StraightRodPair* rod = i%2 ? GeometryFactory::clone(*secondRod) : GeometryFactory::clone(*firstRod);
    rod->myid(i+1);

    // PHI ROTATION
    double deltaPhi = rodCenterPhiShift * (i%2 ? i-1 : i);  // Extra Phi, added to the copy of secondRod or firstRod respectively.
    // phiForbidden Range
    if (phiForbiddenRanges.state()) {
      const double firstRodCenterPhi = firstRod->Phi();
      const double secondRodCenterPhi = secondRod->Phi();
      double iRodPhi = (i%2 ? secondRodCenterPhi : firstRodCenterPhi) + deltaPhi;  // Total Phi of rod number i.
      const double forbiddenPhiLowerB = phiForbiddenRanges.at(2) * M_PI / 180.;
      const double forbiddenPhiUpperB = phiForbiddenRanges.at(3) * M_PI / 180.;
      if (moduloComp(forbiddenPhiLowerB, iRodPhi, 2. * M_PI)) {
	deltaPhi += (forbiddenPhiUpperB - forbiddenPhiLowerB);  // Shift all Phi : + Phi forbidden range B circumference when relevant.
      }
    }
    rod->rotateZ(deltaPhi);

    // STORE
    const bool isFlatPart = isTilted();
    if (!isFlatPart) { rods_.push_back(rod); }
    else { flatPartRods_.push_back(rod); }
  }
}


/* dedicated to skewed mode */

/*
 * Skewed layer: Compute all parameters of interest.
 */
const SkewedLayerPhiShifts Layer::buildSkewed() {

  // any module
  RectangularModule* mod = GeometryFactory::make<RectangularModule>();
  mod->store(propertyTree());
  if (ringNode.count(1) > 0) mod->store(ringNode.at(1));  // no different module types allowed along a rod for the moment!!!
  mod->build();
  const double moduleWidth = mod->width();

  // COMPUTE ALL SKEWED INFO 
  const SkewedLayerInfo& info = computeSkewedLayerInfo(placeRadiusHint(), bigDelta(), numRods(), moduleWidth, 
						       skewedModuleEdgeShift(), installationOverlapRatio());

  // STORE SKEWED INFO OF INTEREST IN THE LAYER
  skewedModuleCenterRho(info.skewedModuleCenterRho);
  skewAngle(info.skewAngle);
  unitPhiOverlapLength(info.unitPhiOverlapLength);
  installationHorizontalOverlapLength(info.installationHorizontalOverlapLength);

  // RETURN PHI SHIFTS INFO
  // This info is not stored in the Layer, because it does not have interest in itself.
  // It will just be used once, to compute rods phi shifts.
  const SkewedLayerPhiShifts& phiShifts = info.phiShifts;
  return phiShifts; 
}


/*
 * Skewed layer: geometry helper.
 * To understand all formulae, look at draft images on tkLayout's wiki.
 * The present code relies on 2 assumptions: 
 - there are only 2 skewed rods.
 - a skewed rod is always placed at +bigDelta.
 * NB: Sensor thickness is not taken into account here, as it is not needed!
 */
const SkewedLayerInfo Layer::computeSkewedLayerInfo(const double layerCenterRho, const double bigDelta, const int numRods, const double moduleWidth, 
						    const double skewedModuleEdgeShift, const double installationOverlapRatio) {

  // skewAperture is the angle between (XZ) plane and the sensor plane.
  const double skewAperture = 2. * asin( 0.5 * skewedModuleEdgeShift / moduleWidth);

  // NON-SKEWED RODS: COMPUTE RADII AND PHI COVERAGE
  const double minusBigDeltaRodCenterRho = layerCenterRho - bigDelta;
  const double plusBigDeltaRodCenterRho = layerCenterRho + bigDelta;
  const double minusBigDeltaRodMaxRho = sqrt( pow(minusBigDeltaRodCenterRho, 2.) + pow(moduleWidth/2. , 2.));
  const double minusBigDeltaRodPhiCov = 2. * atan( 0.5*moduleWidth / minusBigDeltaRodCenterRho);
  const double plusBigDeltaRodPhiCov = 2. * atan( 0.5*moduleWidth / plusBigDeltaRodCenterRho);

  // SKEWED RODS: COMPUTE RADII AND PHI COVERAGE
  const double skewedModuleEdgeMinRho = sqrt( pow(plusBigDeltaRodCenterRho, 2.) + pow(moduleWidth/2. , 2.)); // the skewed rod is at +bigDelta
  const double skewedModuleCenterRho = sqrt( pow(skewedModuleEdgeMinRho, 2) + pow(0.5 * moduleWidth, 2.) 
					     + skewedModuleEdgeMinRho * moduleWidth * sin(skewAperture - 0.5 * plusBigDeltaRodPhiCov)
					     ); // Al-Kashi
  const double skewedRodLowerHalfPhiCov = acos( (pow(skewedModuleEdgeMinRho, 2.) + pow(skewedModuleCenterRho, 2.) - pow((0.5 * moduleWidth), 2.)) 
					       / ( 2. * skewedModuleEdgeMinRho * skewedModuleCenterRho)
						); // Al-Kashi
  const double skewAngle = skewAperture - (plusBigDeltaRodPhiCov / 2. - skewedRodLowerHalfPhiCov);

  const double beta = atan( (moduleWidth/2. - skewedModuleEdgeShift * sin(skewAperture/2.)) 
			    / (plusBigDeltaRodCenterRho + skewedModuleEdgeShift * cos(skewAperture/2.)) 
			    );   // the skewed rod is at +bigDelta  
  const double skewedRodTotalPhiCov = plusBigDeltaRodPhiCov / 2. + beta;

  const double skewedModuleMaxRho = (plusBigDeltaRodCenterRho + skewedModuleEdgeShift * cos(skewAperture/2.)) / cos(beta); // the skewed rod is at +bigDelta
  /*const double skewedRodUpperHalfPhiCov = acos( 
					       (pow(skewedModuleMaxRho, 2.) + pow(skewedModuleCenterRho, 2.) - pow((0.5 * moduleWidth), 2.)) 
					       / ( 2. * skewedModuleMaxRho * skewedModuleCenterRho)
					       );*/

  // NUMBER OF RODS
  const int numSkewedRods = 2; // Harcoded number of skewed rods!!
  const int numMinusBigDeltaRods = numRods / 2;
  const int numPlusBigDeltaRods = numRods / 2 - numSkewedRods; // the skewed rods are at +bigDelta

  // COMPUTE PHI OVERLAPS!!
  const double totalPhiCov = numMinusBigDeltaRods * minusBigDeltaRodPhiCov 
    + numPlusBigDeltaRods * plusBigDeltaRodPhiCov 
    + numSkewedRods * skewedRodTotalPhiCov;
  const double totalPhiOverlap = totalPhiCov - 2. * M_PI;
  const double unitPhiOverlapAngle = totalPhiOverlap / (numRods - numSkewedRods + numSkewedRods * installationOverlapRatio);
  const double unitPhiOverlapLength = unitPhiOverlapAngle * layerCenterRho;

  
  const double installationPhiOverlapAngle = installationOverlapRatio * unitPhiOverlapAngle;
  const double skewedInstallationPhiOverlapAngle = atan( (minusBigDeltaRodMaxRho * sin(installationPhiOverlapAngle)) 
							 / (skewedModuleMaxRho + minusBigDeltaRodMaxRho * cos(installationPhiOverlapAngle)) 
							 );
  const double unskewedInstallationPhiOverlapAngle = installationPhiOverlapAngle - skewedInstallationPhiOverlapAngle;
  const double installationHorizontalOverlapLength = skewedModuleMaxRho * sin(skewedInstallationPhiOverlapAngle);
    


  // COMPUTE PHI SHIFTS (TO BE USED FOR ROD PLACEMENT IN PHI)
  //installationSkewedRodCenterPhiShift = skewedRodTotalPhiCov / 2. - skewedInstallationPhiOverlapAngle;
  //nextRodCenterPhiShift = skewedRodTotalPhiCov / 2. + minusBigDeltaRodPhiCov / 2. - unitPhiOverlapAngle;
  const double installationMinusBigDeltaRodCenterPhiShift = minusBigDeltaRodPhiCov / 2. - unskewedInstallationPhiOverlapAngle;
  const double commonRodCenterPhiShift = minusBigDeltaRodPhiCov / 2. + plusBigDeltaRodPhiCov / 2. - unitPhiOverlapAngle;
  const double skewedRodCenterPhiShift = skewedRodLowerHalfPhiCov + minusBigDeltaRodPhiCov / 2. - unitPhiOverlapAngle; // the skewed rod is at +bigDelta 


  // GATHER ALL RESULTS OF INTEREST
  const SkewedLayerPhiShifts& phiShifts = SkewedLayerPhiShifts{ installationMinusBigDeltaRodCenterPhiShift, 
								commonRodCenterPhiShift, 
								skewedRodCenterPhiShift };
  const SkewedLayerInfo& info = SkewedLayerInfo{ skewedModuleCenterRho, skewAngle, 
						 unitPhiOverlapLength, installationHorizontalOverlapLength, phiShifts };

  return info;
}


/*
 * Skewed mode:
 * General method to build & store all rods.
 */
void Layer::buildAndStoreClonedRodsInSkewedMode(const StraightRodPair* firstRod, const StraightRodPair* secondRod, StraightRodPair* skewedRod,
						const double commonRodCenterPhiShift, const double skewedRodCenterPhiShift) {
  
  const int numRodsPerXSide = numRods() / 2;
  double lastNonSkewedRodCenterPhi = 0.;  // Phi of the center of the last non-skewed rod which has been computed.
                                          // This will in turn be used to compute the phi of the skewed rod.

  // LOOP ON ALL RODS
  // NB: Works on one (X) side.
  // The other (X) side is automatically computed on the fly, as a rotation in Phi by Pi.
  for (int rodId = 1; rodId <= (numRodsPerXSide); rodId++) {
    
    // NON-SKEWED RODS
    if (rodId != numRodsPerXSide) {
      lastNonSkewedRodCenterPhi = buildAndStoreNonSkewedRodsInSkewedMode(rodId, numRodsPerXSide,
									 firstRod, secondRod,
									 commonRodCenterPhiShift);
    }

    // SKEWED ROD
    else {
      buildAndStoreSkewedRods(numRodsPerXSide,
			      skewedRod,
			      lastNonSkewedRodCenterPhi, skewedRodCenterPhiShift);
    }

  }
}


/*
 * Skewed mode:
 * Method to build & store all non-skewed rods.
 */
double Layer::buildAndStoreNonSkewedRodsInSkewedMode(const int rodId, const int numRodsPerXSide,
						     const StraightRodPair* firstRod, const StraightRodPair* secondRod,
						     const double commonRodCenterPhiShift) {
  double lastNonSkewedRodCenterPhi = 0.;

  // CLONE RODS
  // if rodId is odd: clone firstRod, otherwise clone secondRod.
  StraightRodPair* rod = (rodId-1)%2 ? GeometryFactory::clone(*secondRod) : GeometryFactory::clone(*firstRod);
  const double firstRodCenterPhi = firstRod->Phi();
  const double secondRodCenterPhi = secondRod->Phi();
  const double rodClonedPhi = (rodId-1)%2 ? secondRodCenterPhi : firstRodCenterPhi;  
  rod->myid(rodId);
  lastNonSkewedRodCenterPhi = rodClonedPhi;

  // FIRST 2 RODS
  // if bigParity() > 0, mirror all rods phi placements through the (XZ) plane.
  // In practice, need to mirror through (YZ) plane (ie, consider M_PI - angle instead of angle).
  // Then everytheing end up being rotated by a global Pi/2 rotation around CMS_Z.
  if (rodId <= 2 && bigParity() > 0) {
    const double rodPhiPosition = M_PI - rodClonedPhi; // mirror through (XZ) plane.
    rod->rotateZ(-rodClonedPhi + rodPhiPosition);      // -rodClonedPhi is to remove whatever phi position the rod already has.
    lastNonSkewedRodCenterPhi = rodPhiPosition;
  }

  // CLONE THE FIRST 2 RODS, PLACE AND STORE.
  if (rodId >= 3) {
    // Phi rotation
    const double rodPhiShift = secondRodCenterPhi + (rodId - 2) * commonRodCenterPhiShift;
    const double rodPhiPosition = (bigParity() < 0 ? rodPhiShift : M_PI - rodPhiShift);
    rod->rotateZ(-rodClonedPhi + rodPhiPosition);      // -rodClonedPhi is to remove whatever phi position the rod already has.
    lastNonSkewedRodCenterPhi = rodPhiPosition;
  }

  // STORE
  rods_.push_back(rod);
      
  // THE OTHER (X) SIDE HALF IS BUILT ON THE FLY HERE.
  // Build Rod at + PI.
  StraightRodPair* rotatedByPiInPhiRod = buildRotatedByPiInPhiRod(rod, numRodsPerXSide);
  // Store Rod at + PI.
  rods_.push_back(rotatedByPiInPhiRod);

  return lastNonSkewedRodCenterPhi;
}


/*
 * Skewed mode:
 * Method to build & store the skewed rods.
 */
void Layer::buildAndStoreSkewedRods(const int numRodsPerXSide,
				    StraightRodPair* skewedRod,
				    const double lastNonSkewedRodCenterPhi, const double skewedRodCenterPhiShift) {

  // BUILD A ROD WITH SKEWED MODULES.
  const double orientedSkewAngle = bigParity() * skewAngle();  // orientation of the skew angle in the (XY) plane
  RodTemplate skewedRodTemplate = makeRodTemplate(orientedSkewAngle);

  // PLACE SKEWED ROD RADIALLY
  const bool isPlusBigDeltaRod = true;                         // the skewed rod is at +bigDelta
  skewedRod->isOuterRadiusRod(isPlusBigDeltaRod);
  skewedRod->build(skewedRodTemplate, isPlusBigDeltaRod);
  skewedRod->translateR(skewedModuleCenterRho());
  skewedRod->myid(numRodsPerXSide);

  // PLACE SKEWED ROD IN PHI
  const double skewedRodPhi = lastNonSkewedRodCenterPhi - bigParity() * skewedRodCenterPhiShift;
  skewedRod->rotateZ(skewedRodPhi);

  // STORE
  rods_.push_back(skewedRod);

  // THE OTHER (X) SIDE HALF IS BUILT ON THE FLY HERE.
  // Build Rod at + PI.
  StraightRodPair* rotatedByPiInPhiRod = buildRotatedByPiInPhiRod(skewedRod, numRodsPerXSide);
  // Store Rod at + PI.
  rods_.push_back(rotatedByPiInPhiRod);

  skewedModuleMinRho(skewedRod->minR());  // takes sensor thickness into account. WARNING: min Rho is not compulsory reached at the skewed sensor edge!!
  skewedModuleMaxRho(skewedRod->maxR());  // takes sensor thickness into account
}


/*
 * Build & return a clone of initialRod, but located at +Pi in Phi.
 */
StraightRodPair* Layer::buildRotatedByPiInPhiRod(const StraightRodPair* initialRod, const int numRodsPerXSide) const {
  if (!initialRod) {
    throw PathfulException("Tried to clone a rod which does not exist.");
  }
  else {
    StraightRodPair* rotatedByPiInPhiRod = GeometryFactory::clone(*initialRod);
    const int initialRodId = initialRod->myid();
    rotatedByPiInPhiRod->myid(initialRodId + numRodsPerXSide);
    rotatedByPiInPhiRod->rotateZ(M_PI);

    return rotatedByPiInPhiRod;
  }
}


/* TILTED LAYER */

/*
 * Build a tilted layer.
 */
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
      buildStraight();


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


    const int numInnerUntiltedModules = tmspecsi.size();
    if (numInnerUntiltedModules != buildNumModulesFlat()) {
      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " but flat part inner rod has " + to_string(tmspecsi.size()) + " module(s).");
    }
    const int numOuterUntiltedModules = tmspecso.size();
    if (numOuterUntiltedModules != buildNumModulesFlat()) {
      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " but flat part outer rod has " + to_string(tmspecso.size()) + " module(s).");
    }
    



    tiltedRingsGeometry_ = makeTiltedRingsTemplate(flatPartThetaEnd);

    const int numTiltedModules = tiltedRingsGeometry_.size();
    if (numTiltedModules == buildNumModulesTilted()) {
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

    const int numInnerModules = tmspecsi.size();
    if (numInnerModules != buildNumModules()) {
      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " and numModulesTilted = " + to_string(buildNumModulesTilted()) + " but tilted rod 1 has " + to_string(tmspecsi.size()) + " module(s) in total.");
    }
    const int numOuterModules = tmspecso.size();
    if (numOuterModules != buildNumModules()) {
      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " and numModulesTilted = " + to_string(buildNumModulesTilted()) + " but tilted rod 2 has " + to_string(tmspecso.size()) + " module(s) in total.");
    }

  }

  RodTemplate rodTemplate = makeRodTemplate();

  const double rodCenterPhiShift = 2 * M_PI / numRods();
  const double firstRodCenterPhi = (rotateLayerByRodsDeltaPhiHalf()) ? rodCenterPhiShift / 2. : 0;

  for (int i = 0; i < numRods(); i++) {
    TiltedRodPair* rod = GeometryFactory::make<TiltedRodPair>(subdetectorName());
    rod->store(propertyTree());
    rod->myid(i+1);
    bool isRodAtOuterRadius;
    // For positive bigParity, even rods are at outer radius
    if (i%2==0) isRodAtOuterRadius = (bigParity() > 0 ? true : false);
    // For positive bigParity, odd rods are are lower radius
    else isRodAtOuterRadius = (bigParity() > 0 ? false : true);
    const std::vector<TiltedModuleSpecs> myRodSensorsCenters = (isRodAtOuterRadius ? tmspecso : tmspecsi);
    rod->isOuterRadiusRod(isRodAtOuterRadius); 
    rod->build(rodTemplate, myRodSensorsCenters, !isRodAtOuterRadius);
    rod->rotateZ(firstRodCenterPhi + rodCenterPhiShift * i);
    rods_.push_back(rod);
  }

  // computing the layer's place radius as the average of all the modules' radii
  placeRadius_  = std::accumulate(tmspecsi.begin(), tmspecsi.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
  placeRadius_ += std::accumulate(tmspecso.begin(), tmspecso.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
  placeRadius_ /= tmspecsi.size() + tmspecso.size();

}


/*
 * Create a template tilted rod.
 */
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

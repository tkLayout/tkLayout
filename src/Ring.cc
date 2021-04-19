#include "Ring.hh"
#include "MessageLogger.hh"

inline void Ring::check() {
  PropertyObject::check();
  if (moduleShape() == ModuleShape::WEDGE && buildDirection() == TOPDOWN) throw PathfulException("Wedge-module rings cannot be built top down");
}

inline double Ring::solvex(double y) {
  return(1/8.*(3*y+2-sqrt(9*pow(y, 2)-4*y+4)));
}

inline double Ring::compute_l(double x, double y, double d) {
  double result = d;
  result /= (1-x)-sqrt((1-x)*(y-x));
  return result;
}

inline double Ring::compute_d(double x, double y, double l) {
  double result = l;
  result *= (1-x)-sqrt((1-x)*(y-x));
  return result;
}


void Ring::cutAtEta(double eta) { modules_.erase_if([eta](EndcapModule& m) { return fabs(m.center().Eta()) > eta; }); }

double Ring::computeTentativePhiAperture(double moduleWaferDiameter, double minRadius) {
  double r = moduleWaferDiameter/2;

  double l = (minRadius-r);
  double y = pow(r/l, 2);
  double x = solvex(y);


  double tempd;

  int i = 0;
  for (; i < MAX_WEDGE_CALC_LOOPS; i++) {
    l = compute_l(x, y, minRadius);
    y = pow(r/l, 2);
    x = solvex(y);

    tempd = compute_d(x, y, l);

    if (fabs(minRadius - tempd)<1e-15) break;
  }

  if (i >= MAX_WEDGE_CALC_LOOPS) {
    //logWarning("Maximum number of iterations hit while computing wedge geometry");
  }

  double alpha = asin(sqrt(x)) * 2;

  return alpha;
}

std::pair<double, int> Ring::computeOptimalRingParametersWedge(double moduleWaferDiameter, double minRadius) {
  double delta = phiOverlap()/minRadius;// SM: The needed overlap becomes an angle delta by
  //     checking the unsafest point (r=r_min)

  double tentativeAlpha = computeTentativePhiAperture(moduleWaferDiameter, minRadius) - delta;
  float tentativeNumMods = 2*M_PI / tentativeAlpha; 
  int optimalNumMods = (!requireOddModsPerSlice()? round(tentativeNumMods/phiSegments()) : roundToOdd(tentativeNumMods/phiSegments()) + additionalModules()) * phiSegments();
  float optimalAlpha = 2*M_PI/optimalNumMods + delta;

  return std::make_pair(optimalAlpha, optimalNumMods);
}

std::pair<double, int> Ring::computeOptimalRingParametersRectangle(double moduleWidth, double highRadius) {
  double effectiveWidth = moduleWidth - phiOverlap();
  double optimalAlpha = 2 * atan( effectiveWidth / (2. * highRadius));
  double tentativeNumMods = 2 * M_PI/ optimalAlpha;
  int modsPerSlice = ceil(tentativeNumMods/phiSegments());
  if ((modsPerSlice % 2) == 0 && requireOddModsPerSlice()) modsPerSlice++;
  int optimalNumMods = modsPerSlice * phiSegments();

  return std::make_pair(optimalAlpha, optimalNumMods);
}


void Ring::buildModules(EndcapModule* templ, int numMods, double smallDelta, double phiShift, double modTranslateX) {
  double alignmentRotation = alignEdges() ? 0.5 : 0.;
  double deltaPhiNom = 2.*M_PI/numMods;
  double deltaPhiLarge = deltaPhiNom+2*phiShift;
  double nominalZRot = 0.;
  if (deltaPhiLarge!=deltaPhiNom){ deltaPhiNom = deltaPhiNom-(4*phiShift/(numMods-2));}
  templ->store(propertyTree());

  for (int i = 0, parity = smallParity(); i < numMods; i++, parity *= -1) {
    if (moduleNode.count(i) > 0) {
     templ->store(moduleNode.at(i));
    }
    EndcapModule* mod = GeometryFactory::clone(*templ);
    mod->build();
    mod->translate(XYZVector(modTranslateX, 0, 0));
    mod->myid(i+1);
    if(mod->yawAngleFromConfig() > -999){
      mod->notInRegularRing();
      double tmp_r = mod->center().Rho();
      mod->translateR(-tmp_r); //perform yaw angle rotation with the module at the centre
      mod->yaw(mod->yawAngleFromConfig());
      if(mod->manualRhoCentre() > 0 ){// For irregular rings (with yawed modules) that need to stay tangent to the same inner/outer radius, the module centre is shifted
        mod->translateR(mod->manualRhoCentre()-mod->center().Rho());
      } else {
        mod->translateR(tmp_r-mod->center().Rho());
      }
    }
    if (deltaPhiLarge == deltaPhiNom) {
      nominalZRot = (i + alignmentRotation)*deltaPhiNom;
    } else {
      mod->notInRegularRing();
      if (((i + alignmentRotation)*deltaPhiNom + alignmentRotation*(deltaPhiLarge-deltaPhiNom)) < M_PI ){
        nominalZRot = (i + alignmentRotation)*deltaPhiNom + alignmentRotation*(deltaPhiLarge-deltaPhiNom);
      } else {
        nominalZRot = (i + alignmentRotation)*deltaPhiNom + (1+alignmentRotation)*(deltaPhiLarge-deltaPhiNom);
      } 
    }
    mod->rotateZ(nominalZRot);
    mod->rotateZ(zRotation());
    mod->translateZ(parity*smallDelta);
    mod->setIsSmallerAbsZModuleInRing(parity < 0);
    const bool isFlipped = (!isRingOn4Dees() ? (parity < 0) // Case where ring modules are placed on both sides of a same dee.
			                                    // Half of the ring modules are flipped: 
			                                    // the ones placed with parity < 0, hence on the small |Z| side.
			    : isSmallerAbsZRingInDisk()); // Case where ring modules are placed on 4 dees.
                                                          // If the ring is on the small |Z| side with respect to the disk:
                                                          // all ring modules are flipped.
    mod->flipped(isFlipped);
    modules_.push_back(mod);  
  }
}


void Ring::buildBottomUp() {
  double startRadius = buildStartRadius()+ringGap();
  if (ringOuterRadius()>0){
    logWARNING("outer radius was set for a bottom-up endcap building. Ignoring ringOuterRadius.");
  }
  if (ringInnerRadius()>0)  {
    if (ringGap()!=0) logWARNING("innerRadius and ringGap were both specified. Ignoring ringGap.");
    startRadius = ringInnerRadius();
  }
  buildStartRadius(startRadius);

  int numMods;
  double modLength;

  EndcapModule* emod = nullptr;
  if (moduleShape() == ModuleShape::WEDGE) {
    WedgeModule* wmod = GeometryFactory::make<WedgeModule>();
    wmod->store(propertyTree());

    auto optimalRingParms = computeOptimalRingParametersWedge(wmod->waferDiameter(), buildStartRadius());
    double alpha = optimalRingParms.first;
    numMods = optimalRingParms.second;

    wmod->buildAperture(alpha);
    wmod->buildDistance(buildStartRadius());
    wmod->buildCropDistance(buildCropRadius());

    wmod->build();

    modLength = wmod->length();
    emod = GeometryFactory::make<EndcapModule>(wmod, subdetectorName());

  } else {

    RectangularModule* rmod = GeometryFactory::make<RectangularModule>();
    rmod->store(propertyTree());
    rmod->build();

    auto optimalRingParms = computeOptimalRingParametersRectangle(rmod->width(), buildStartRadius() + rmod->length());
    numMods = optimalRingParms.second;

    modLength = rmod->length();

    emod = GeometryFactory::make<EndcapModule>(rmod, subdetectorName()); 

  }

  minRadius_ = buildStartRadius();
  maxRadius_ = buildStartRadius() + modLength;

  if (numModules.state()) numMods = numModules();
  else numModules(numMods);
  buildModules(emod, numMods, smallDelta(),phiShift(),buildStartRadius()+modLength/2);

  delete emod;
}


void Ring::buildTopDown() {
  double startRadius = buildStartRadius()-ringGap();
  if (ringInnerRadius()>0){
    logWARNING("inner radius was set for a top-down endcap building. Ignoring ringInnerRadius.");
  }
  if (ringOuterRadius()>0)  {
    if (ringGap()!=0) logWARNING("outerRadius and ringGap were both specified. Ignoring ringGap.");
    startRadius = ringOuterRadius();
  }
  buildStartRadius(startRadius);

  RectangularModule* rmod = GeometryFactory::make<RectangularModule>();
  rmod->store(propertyTree());

  auto optimalRingParms = computeOptimalRingParametersRectangle(rmod->width(), buildStartRadius());
  int numMods = optimalRingParms.second;

  EndcapModule* emod = GeometryFactory::make<EndcapModule>(rmod, subdetectorName());

  minRadius_ = buildStartRadius() - rmod->length();
  maxRadius_ = buildStartRadius();

  if (numModules.state()) numMods = numModules();
  else numModules(numMods);
  buildModules(emod, numMods, smallDelta(),phiShift(),buildStartRadius() - rmod->length()/2);

  delete emod;
}


void Ring::build() {
  materialObject_.store(propertyTree());
  //std::cout << "Ring::build()  : materialObject_.subdetectorName() = " << materialObject_.subdetectorName() << std::endl;
  materialObject_.build();

  if(materialObject_.isPopulated()) {
    logUniqueWARNING("Is not possible to define rod material for disks, material ignored.");
  }

  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();
    if (buildDirection() == BOTTOMUP) buildBottomUp();
    else buildTopDown();
  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  cleanup();
  builtok(true);
}

void Ring::translateZ(double z) {
  for (auto& m : modules_) {
    m.translateZ(z);
  }
}

void Ring::rotateToNegativeZSide() {
  for (auto& m : modules_) {
    m.rotateToNegativeZSide();
  }
}

/** This computes the actual coverage in Phi of a disk (after it is built).
    It calcualtes the actual phiOverlap, using the relevant coordinates of the disk.
 */
void Ring::computeActualPhiCoverage() {
  double width = modules_.at(0).minWidth();
  double radiusHigh = modules_.at(0).minR() + modules_.at(0).length();
  int numMods = numModules();

  double phiOverlap = width - 2. * radiusHigh * tan(M_PI / numMods);

  // STORE THE RESULT
  actualPhiOverlap(phiOverlap);
}

const MaterialObject& Ring::materialObject() const{
  return materialObject_;
}

define_enum_strings(Ring::BuildDirection) = { "topdown", "bottomup" };















void TiltedRing::check() {
  PropertyObject::check();
  if (ringZOverlap.state()) {
    if (zInner.state() || zOuter.state()) throw PathfulException("Only one parameter among ringZOverlap, ringInnerZ, and ringOuterZ can be specified.");
  }
  else {
    if (!zInner.state() && !zOuter.state()) throw PathfulException("At least one parameter among ringZOverlap, ringInnerZ, and ringOuterZ must be specified.");
    if (zInner.state() && zOuter.state()) throw PathfulException("Only one parameter among ringZOverlap, ringInnerZ, and ringOuterZ can be specified.");
  }
}




void TiltedRing::buildLeftRight(double lastThetaEnd) {

  thetaStart_ = lastThetaEnd;
  double tilt = tiltAngle() * M_PI / 180.;
  double theta_gRad = theta_g() * M_PI / 180.;

  RectangularModule* rmod = GeometryFactory::make<RectangularModule>();
  rmod->store(propertyTree());
  rmod->build();
  double dsDistance = rmod->dsDistance();
  double length = rmod->length();
  double lengthEff;
  double width = rmod->width();
  

  

  if (thetaStart_ == (M_PI / 2.)) {
    throw PathfulException("Start building tilted rings at thetaStart = M_PI/2.");
    //thetaOuterUP_ = M_PI / 2.;
    //thetaOuterDOWN_ = M_PI / 2.;
    //thetaOuter_ = M_PI / 2.;
    //zOuter(0.0);
    //zInner(0.0);
  }

  else {

    // CASE A : ZOVERLAP IS SPECIFIED IN INPUT
    if (ringZOverlap.state()) {

      // Calculate lengthEff
      lengthEff = length - 2.*ringZOverlap();

      // Calculate thetaOuter
      thetaOuterUP_ = atan( outerRadius() / ( outerRadius()/tan(thetaStart_) + dsDistance*cos(tilt)/(2.*tan(thetaStart_)) + lengthEff*sin(tilt)/(2.*tan(thetaStart_)) - dsDistance/2.*sin(tilt) + lengthEff/2.*cos(tilt) ));

      thetaOuterDOWN_ = atan( outerRadius() / ( outerRadius()/tan(thetaStart_) - dsDistance*cos(tilt)/(2.*tan(thetaStart_)) + lengthEff*sin(tilt)/(2.*tan(thetaStart_)) + dsDistance/2.*sin(tilt) + lengthEff/2.*cos(tilt) ));

      thetaOuter_ = MAX(thetaOuterUP_, thetaOuterDOWN_);

      // Calculate zOuter
      zOuter(outerRadius() / tan(thetaOuter_));
    
      /*std::cout << "thetaStart_ * 180. / M_PI = " << thetaStart_ * 180. / M_PI << std::endl;
	std::cout << "outerRadius() = " << outerRadius() << std::endl;
	std::cout << "lengthEff = " << lengthEff << std::endl;
	std::cout << "ringZOverlap() = " << ringZOverlap() << std::endl;
	std::cout << "thetaOuterUP_  * 180. / M_PI = " << thetaOuterUP_  * 180. / M_PI << std::endl;
	std::cout << "thetaOuterDOWN_  * 180. / M_PI = " << thetaOuterDOWN_  * 180. / M_PI << std::endl;
	std::cout << "thetaOuter_  * 180. / M_PI = " << thetaOuter_  * 180. / M_PI << std::endl;
	std::cout << "zOuter() = " << zOuter() << std::endl;*/

      // Calculate zInner
      zInner( zOuter() - (outerRadius() - innerRadius()) / tan(theta_gRad));
    }


    // CASE B : ZINNER OR ZOUTER IS SPECIFIED IN INPUT
    else {

      // If zOuter is set, calculate zInner
      if (zOuter.state()) {
	if (zOuter() == 0.) throw PathfulException("Start building tilted rings at zOuter = 0.");
	zInner( zOuter() - (outerRadius() - innerRadius()) / tan(theta_gRad));
      }

      // If zInner is set, calculate zOuter
      if (zInner.state()) zOuter( zInner() + (outerRadius() - innerRadius()) / tan(theta_gRad));

      // Calculate thetaOuter
      thetaOuter_ = atan( outerRadius() / zOuter());

      // Calculate ringZOverlap
      double ringZOverlapUP = 0.5 * ( length - (outerRadius()/tan(thetaOuter_) - outerRadius()/tan(thetaStart_) - dsDistance*cos(tilt)/(2.*tan(thetaStart_)) + dsDistance*sin(tilt)/2. ) / ( sin(tilt)/(2.*tan(thetaStart_)) + cos(tilt)/2. ) );

      double ringZOverlapDOWN = 0.5 * ( length - (outerRadius()/tan(thetaOuter_) - outerRadius()/tan(thetaStart_) + dsDistance*cos(tilt)/(2.*tan(thetaStart_)) - dsDistance*sin(tilt)/2. ) / ( sin(tilt)/(2.*tan(thetaStart_)) + cos(tilt)/2. ) );

      //std::cout << " ringZOverlapUP = " <<  ringZOverlapUP <<  "ringZOverlapDOWN = " << ringZOverlapDOWN << std::endl;
      
      ringZOverlap( MIN(ringZOverlapUP, ringZOverlapDOWN) );

      // Calculate lengthEff
      lengthEff = length - 2.*ringZOverlap();
    }
     
  }




  // MODULE 2 (OUTER MODULE)

  tiltAngleIdealOuter_ = 90. - thetaOuter_ * 180. / M_PI;
  deltaTiltIdealOuter_ = tiltAngle() - tiltAngleIdealOuter_;

  //std::cout << "zOuter() = " << zOuter() << std::endl;
  //std::cout << "zInner() = " << zInner() << std::endl;


  //double zH2p = zOuter() - 0.5 * lengthEff * cos(tilt);
  //double rH2p = outerRadius() + 0.5 * lengthEff * sin(tilt);
  //double zH2pp = zOuter() + 0.5 * lengthEff * cos(tilt);
  //double rH2pp = outerRadius() - 0.5 * lengthEff * sin(tilt);

  double zH2UP = zOuter() + 0.5 * dsDistance * sin(tilt);
  double rH2UP = outerRadius() + 0.5 * dsDistance * cos(tilt);
  //double zH2pUP = zH2UP - 0.5 * lengthEff * cos(tilt);
  double rH2pUP = rH2UP + 0.5 * lengthEff * sin(tilt);
  //double zH2ppUP = zH2UP + 0.5 * lengthEff * cos(tilt);
  //double rH2ppUP = rH2UP - 0.5 * lengthEff * sin(tilt);

  double zH2DOWN = zOuter() - 0.5 * dsDistance * sin(tilt);
  double rH2DOWN = outerRadius() - 0.5 * dsDistance * cos(tilt);
  //double zH2pDOWN = zH2DOWN - 0.5 * lengthEff * cos(tilt);
  //double rH2pDOWN = rH2DOWN + 0.5 * lengthEff * sin(tilt);
  //double zH2ppDOWN = zH2DOWN + 0.5 * lengthEff * cos(tilt);
  //double rH2ppDOWN = rH2DOWN - 0.5 * lengthEff * sin(tilt);


  /*std::cout << "zH2ppUP = " << zH2ppUP << " zH2ppDOWN = " << zH2ppDOWN << std::endl;
  std::cout << "rH2ppUP = " << rH2ppUP <<" rH2ppDOWN = " << rH2ppDOWN << std::endl; 
  std::cout << "atan(rH2ppUP / zH2ppUP) = " << atan(rH2ppUP / zH2ppUP) << std::endl;
  std::cout << "MAX( atan(rH2ppUP / zH2ppUP), atan(rH2ppDOWN / zH2ppDOWN)) = " << MAX( atan(rH2ppUP / zH2ppUP), atan(rH2ppDOWN / zH2ppDOWN)) << std::endl;*/
  //thetaEnd_ = MAX( atan(rH2ppUP / zH2ppUP), atan(rH2ppDOWN / zH2ppDOWN));
  //std::cout << "thetaEnd_ = " << thetaEnd_ << std::endl;





  // MODULE 1 (INNER MODULE)

  thetaInner_ = atan( innerRadius() / zInner() );
  tiltAngleIdealInner_ = 90. - thetaInner_ * 180. / M_PI;
  deltaTiltIdealInner_ = tiltAngle() - tiltAngleIdealInner_;


  //double zH1p = zInner() - 0.5 * lengthEff * cos(tilt);
  //double rH1p = innerRadius() + 0.5 * lengthEff * sin(tilt);
  //double zH1pp = zInner() + 0.5 * lengthEff * cos(tilt);
  //double rH1pp = innerRadius() - 0.5 * lengthEff * sin(tilt);

  double zH1UP = zInner() + 0.5 * dsDistance * sin(tilt);
  double rH1UP = innerRadius() + 0.5 * dsDistance * cos(tilt);
  //double zH1pUP = zH1UP - 0.5 * lengthEff * cos(tilt);
  double rH1pUP = rH1UP + 0.5 * lengthEff * sin(tilt);
  //double zH1ppUP = zH1UP + 0.5 * lengthEff * cos(tilt);
  //double rH1ppUP = rH1UP - 0.5 * lengthEff * sin(tilt);

  double zH1DOWN = zInner() - 0.5 * dsDistance * sin(tilt);
  double rH1DOWN = innerRadius() - 0.5 * dsDistance * cos(tilt);
  //double zH1pDOWN = zH1DOWN - 0.5 * lengthEff * cos(tilt);
  //double rH1pDOWN = rH1DOWN + 0.5 * lengthEff * sin(tilt);
  //double zH1ppDOWN = zH1DOWN + 0.5 * lengthEff * cos(tilt);
  //double rH1ppDOWN = rH1DOWN - 0.5 * lengthEff * sin(tilt);

  //thetaStartInner_ = MIN( atan(rH1pUP / zH1pUP), atan(rH1pDOWN / zH1pDOWN));
  //thetaEndInner_ = MAX( atan(rH1ppUP / zH1ppUP), atan(rH1ppDOWN / zH1ppDOWN));





  // FOR INFO

  //phiOverlapDEG_ = atan(width / (2.* rH2pUP)) + atan(width / (2.* rH1pUP)) - 2. * M_PI / numPhi();

  double T = tan(2.*M_PI / numPhi());
  double A = 1. / (2. * rH1pUP);
  double B = 1. / (2. * rH2pUP);

  double a = T * A * B;
  double b = - (A + B);
  double c = - T;

  double s = (-b - sqrt(b*b - 4*a*c))/(2*a);

  phiOverlap_ = width + s;



  // REAL COORDS (MODULE LENGTH WITH NO Z OVERLAP)

  double zH2pUP_REAL = zH2UP - 0.5 * length * cos(tilt);
  double rH2pUP_REAL = rH2UP + 0.5 * length * sin(tilt);
  double zH2ppUP_REAL = zH2UP + 0.5 * length * cos(tilt);
  double rH2ppUP_REAL = rH2UP - 0.5 * length * sin(tilt);

  double zH2pDOWN_REAL = zH2DOWN - 0.5 * length * cos(tilt);
  double rH2pDOWN_REAL = rH2DOWN + 0.5 * length * sin(tilt);
  double zH2ppDOWN_REAL = zH2DOWN + 0.5 * length * cos(tilt);
  double rH2ppDOWN_REAL = rH2DOWN - 0.5 * length * sin(tilt);

  if ( fabs(rH2pUP_REAL / zH2pUP_REAL) < fabs(rH2pDOWN_REAL / zH2pDOWN_REAL) ) { 
    rStartOuter_REAL_ = rH2pUP_REAL; 
    zStartOuter_REAL_ = zH2pUP_REAL;
  }
  else { 
    rStartOuter_REAL_ = rH2pDOWN_REAL;
    zStartOuter_REAL_ = zH2pDOWN_REAL;
  }
  if ( fabs(rH2ppUP_REAL / zH2ppUP_REAL) > fabs(rH2ppDOWN_REAL / zH2ppDOWN_REAL) ) { 
    rEndOuter_REAL_ = rH2ppUP_REAL; 
    zEndOuter_REAL_ = zH2ppUP_REAL;
  }
  else { 
    rEndOuter_REAL_ = rH2ppDOWN_REAL; 
    zEndOuter_REAL_ = zH2ppDOWN_REAL;
  }


  double zH1pUP_REAL = zH1UP - 0.5 * length * cos(tilt);
  double rH1pUP_REAL = rH1UP + 0.5 * length * sin(tilt);
  double zH1ppUP_REAL = zH1UP + 0.5 * length * cos(tilt);
  double rH1ppUP_REAL = rH1UP - 0.5 * length * sin(tilt);

  double zH1pDOWN_REAL = zH1DOWN - 0.5 * length * cos(tilt);
  double rH1pDOWN_REAL = rH1DOWN + 0.5 * length * sin(tilt);
  double zH1ppDOWN_REAL = zH1DOWN + 0.5 * length * cos(tilt);
  double rH1ppDOWN_REAL = rH1DOWN - 0.5 * length * sin(tilt);

  if ( (rH1pUP_REAL / zH1pUP_REAL) < (rH1pDOWN_REAL / zH1pDOWN_REAL) ) { 
    rStartInner_REAL_ = rH1pUP_REAL; 
    zStartInner_REAL_ = zH1pUP_REAL;
  }
  else { 
    rStartInner_REAL_ = rH1pDOWN_REAL; 
    zStartInner_REAL_ = zH1pDOWN_REAL;
  }
  if ( (rH1ppUP_REAL / zH1ppUP_REAL) > (rH1ppDOWN_REAL / zH1ppDOWN_REAL) ) { 
    rEndInner_REAL_ = rH1ppUP_REAL; 
    zEndInner_REAL_ = zH1ppUP_REAL;
  }
  else { 
    rEndInner_REAL_ = rH1ppDOWN_REAL; 
    zEndInner_REAL_ = zH1ppDOWN_REAL;
  }


}


void TiltedRing::build(double lastThetaEnd) {
  //materialObject_.store(propertyTree());
  //materialObject_.build();


  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();
    buildLeftRight(lastThetaEnd);

  } catch (PathfulException& pe) {
    std::cout << pe.what() << std::endl; // TO DO : should not be necessary to specify pe.what() !! Problem in fullid from capabilities.h ?
    pe.pushPath(fullid(*this)); 
    throw;
  }

  cleanup();
  builtok(true);
}

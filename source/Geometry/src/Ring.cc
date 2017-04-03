#include "Ring.h"
#include "MessageLogger.h"

define_enum_strings(BuildDirection) = { "topdown", "bottomup" };

//
// Constructor - specify unique id & parse geometry config file using boost property tree & read-in module parameters
//
Ring::Ring(int id, BuildDirection direction, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty) :
      m_materialObject(MaterialObject::ROD),
      m_moduleShape           ("moduleShape"           , parsedAndChecked()),
      phiOverlap              ("phiOverlap"            , parsedOnly(), 1.),
      m_buildDirection        ("buildDirection"        , parsedOnly(), direction),
      m_requireOddModsPerSlice("requireOddModsPerSlice", parsedOnly(), false),
      phiSegments             ("phiSegments"           , parsedOnly(), 4),
      m_additionalModules     ("additionalModules"     , parsedOnly(), 0),
      m_alignEdges            ("alignEdges"            , parsedOnly(), true),
      m_ringGap               ("ringGap"               , parsedOnly(), 0.),
      m_smallParity           ("smallParity"           , parsedOnly(), 1),
      smallDelta              ("smallDelta"            , parsedAndChecked()),
      numModules              ("numModules"            , parsedOnly()),
      zRotation               ("zRotation"             , parsedOnly(), 0.),
      ringOuterRadius         ("ringOuterRadius"       , parsedOnly(), -1.),
      ringInnerRadius         ("ringInnerRadius"       , parsedOnly(), -1.)
{
  // Set the geometry config parameters
  this->myid(id);
  this->store(treeProperty);
  if (nodeProperty.count(id)>0) this->store(nodeProperty.at(id));
}

//
// Build all modules -> use either bottom-up or up-down approach
//
void Ring::build(double ringRadius, double ringOffset) {

  m_materialObject.store(propertyTree());
  m_materialObject.build();

  if(m_materialObject.isPopulated()) {
    logUniqueWARNING("Is not possible to define rod material for disks, material ignored.");
  }

  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    // Built either bottom up->down (rectangular) or down >up (wedge-shaped or rectangular)
    if (m_buildDirection() == BOTTOMUP) buildBottomUp(ringRadius);
    else                                buildTopDown(ringRadius);

    // Translate ring to its final position (from average Z pos)
    translateZ(ringOffset);

  }
  catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  cleanup();
  builtok(true);
}

//
// Setup: link lambda functions to various layer related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
//
void Ring::setup() {

  minZ.setup([this]()       { double min = std::numeric_limits<double>::max(); for (const auto& m : m_modules) min = MIN(min, m.minZ()); return min; });
  maxZ.setup([this]()       { double max = -std::numeric_limits<double>::max();for (const auto& m : m_modules) max = MAX(max, m.maxZ()); return max; });
  minZAllMat.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& m : m_modules) min = MIN(min, m.minZAllMat()); return min; });
  maxZAllMat.setup([this]() { double max = -std::numeric_limits<double>::max();for (const auto& m : m_modules) max = MAX(max, m.maxZAllMat()); return max; });
  minR.setup([this]()       { double min = std::numeric_limits<double>::max(); for (const auto& m : m_modules) min = MIN(min, m.minR()); return min; });
  maxR.setup([this]()       { double max = 0;                                  for (const auto& m : m_modules) max = MAX(max, m.maxR()); return max; });

  maxModuleThickness.setup([this]() {
    double max = 0;
    for (const auto& m : m_modules) {
      max = MAX(max, m.thickness());
    }
    return max;
  });
}

//
// Cross-check parameters provided from geometry configuration file
//
inline void Ring::check() {
  PropertyObject::check();
  if (m_moduleShape() == ModuleShape::WEDGE && m_buildDirection() == TOPDOWN) throw PathfulException("Wedge-module rings cannot be built top down");
}

//
// Limit ring geometry by eta cut
//
void Ring::cutAtEta(double eta) { m_modules.erase_if([eta](EndcapModule& m) { return fabs(m.center().Eta()) > eta; }); }

//
// Helper method translating Ring z position by given offset
//
void Ring::translateZ(double zOffset) {
  m_averageZ += zOffset;
  for (auto& m : m_modules) {
    m.translateZ(zOffset);
  }
}

//
// Helper method mirroring the whole Ringc from zPos to -zPos or vice versa
//
void Ring::mirrorZ() {
  m_averageZ *= -1;
  for (auto& m : m_modules) {
    m.mirrorZ();
  }
}

//
// Helper build method -> build modules based on template module
//
void Ring::buildModules(EndcapModule* templ, int numMods, double smallDelta) {

  double alignmentRotation = m_alignEdges() ? 0.5 : 0.;

  for (int i = 0, parity = m_smallParity(); i < numMods; i++, parity *= -1) {

    EndcapModule* mod = GeometryFactory::clone(*templ);
    mod->myid(i+1);
    mod->rotateZ(2*M_PI*(i+alignmentRotation)/numMods); // CUIDADO had a rotation offset of PI/2
    mod->rotateZ(zRotation());
    mod->translateZ(parity*smallDelta); 
    m_modules.push_back(mod);
  }
}

//
// If modules built within a ring using approach bottom-up
//
void Ring::buildBottomUp(double radius) {

  // Update ring radius using ring gap
  double startRadius = radius+m_ringGap();
  if (ringOuterRadius()>0){
    logWARNING("outer radius was set for a bottom-up endcap building. Ignoring ringOuterRadius.");
  }
  if (ringInnerRadius()>0)  {
    if (m_ringGap()!=0) logWARNING("innerRadius and ringGap were both specified. Ignoring ringGap.");
    startRadius = ringInnerRadius();
  }

  int    numMods;
  double modLength;

  EndcapModule* emod = nullptr;

  // Use Wedge-shaped modules
  if (m_moduleShape() == ModuleShape::WEDGE) {

    // Wedge-shaped module -> create & cache parameters
    WedgeModule* wmod = GeometryFactory::make<WedgeModule>();
    wmod->store(propertyTree());

    // Calculate optimal parameters
    auto optimalRingParms = computeOptimalRingParametersWedge(wmod->waferDiameter(), startRadius);
    double alpha          = optimalRingParms.first;
    numMods               = optimalRingParms.second;

    // Build module
    wmod->buildAperture(alpha);
    wmod->buildDistance(startRadius);
    wmod->buildCropDistance(buildCropRadius());
    wmod->build();

    modLength = wmod->length();
    emod = GeometryFactory::make<EndcapModule>(1, wmod, propertyTree());

  }
  // Use rectangular modules
  else {

    // Rectangular module -> create & cache parameters
    RectangularModule* rmod = GeometryFactory::make<RectangularModule>();
    rmod->store(propertyTree());
    rmod->build();

    // Calculate optimal parameters
    auto optimalRingParms = computeOptimalRingParametersRectangle(rmod->width(), startRadius + rmod->length());
    numMods               = optimalRingParms.second;

    modLength = rmod->length();
    emod = GeometryFactory::make<EndcapModule>(1, rmod, propertyTree());
  }

  emod->build();
  emod->translate(XYZVector(startRadius + modLength/2, 0, 0));

  // Update ring inner/outer radius
  ringInnerRadius(radius);
  ringOuterRadius(radius + modLength);

  if (numModules.state()) numMods = numModules();
  else numModules(numMods);
  buildModules(emod, numMods, smallDelta());

  // Delete template module
  delete emod;
}

//
// If modules built within a ring using approach up-down
//
void Ring::buildTopDown(double radius) {

  // Update ring radius using ring gap
  double startRadius = radius-m_ringGap();
  if (ringInnerRadius()>0){
    logWARNING("inner radius was set for a top-down endcap building. Ignoring ringInnerRadius.");
  }
  if (ringOuterRadius()>0)  {
    if (m_ringGap()!=0) logWARNING("outerRadius and ringGap were both specified. Ignoring ringGap.");
    startRadius = ringOuterRadius();
  }

  RectangularModule* rmod = GeometryFactory::make<RectangularModule>();
  rmod->store(propertyTree());

  auto optimalRingParms = computeOptimalRingParametersRectangle(rmod->width(), startRadius);
  int numMods = optimalRingParms.second;

  EndcapModule* emod = GeometryFactory::make<EndcapModule>(1, rmod, propertyTree());
  emod->build();
  emod->translate(XYZVector(startRadius - rmod->length()/2, 0, 0));

  // Update ring inner/outer radius
  ringOuterRadius(radius);
  ringInnerRadius(radius + rmod->length());

  if (numModules.state()) numMods = numModules();
  else numModules(numMods);
  buildModules(emod, numMods, smallDelta());

  // Delte template module
  delete emod;
}

//
// GeometryVisitor pattern -> ring visitable
//
void Ring::accept(GeometryVisitor& v) {
  v.visit(*this);
  for (auto& m : m_modules) { m.accept(v); }
}

//
// GeometryVisitor pattern -> ring visitable (const. option)
//
void Ring::accept(ConstGeometryVisitor& v) const {
  v.visit(*this);
  for (const auto& m : m_modules) { m.accept(v); }
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


double Ring::computeTentativePhiAperture(double moduleWaferDiameter, double minRadius) {
  double r = moduleWaferDiameter/2;

  double l = (minRadius-r);
  double y = pow(r/l, 2);
  double x = solvex(y);

  bool   calculated = false;
  double tempd;

  int i = 0;
  for (; i < c_maxWedgeCalcLoops; i++) {
    l = compute_l(x, y, minRadius);
    y = pow(r/l, 2);
    x = solvex(y);

    tempd = compute_d(x, y, l);
    if (fabs(minRadius - tempd)<1e-03) { //1e-15 {

      calculated = true;
      break;
    }
  }

  // TODO: Just fix, needs to be done properly
  // Starting value if algorithm fails
  double alpha = 2*M_PI/10.;

  if (!calculated) {
    logWARNING("Maximum number of iterations hit while computing wedge geometry, uses default: 10 number of modules in a ring ");
    //std::cout << ">>> " << "Maximum number of iterations hit while computing wedge geometry" << std::endl;
  }
  else alpha = asin(sqrt(x)) * 2;

  return alpha;
}

//
// Helper method computing ring optimal parameters when wedge-shaped modules used
//
std::pair<double, int> Ring::computeOptimalRingParametersWedge(double moduleWaferDiameter, double minRadius) {
  double delta = phiOverlap()/minRadius;// SM: The needed overlap becomes an angle delta by
  //     checking the unsafest point (r=r_min)

  double tentativeAlpha   = computeTentativePhiAperture(moduleWaferDiameter, minRadius) - delta;
  float  tentativeNumMods = 2*M_PI / tentativeAlpha;
  int    optimalNumMods   = (!m_requireOddModsPerSlice()? round(tentativeNumMods/phiSegments() + m_additionalModules()) : roundToOdd(tentativeNumMods/phiSegments()) + m_additionalModules()) * phiSegments();
  float  optimalAlpha     = 2*M_PI/optimalNumMods + delta;

  return std::make_pair(optimalAlpha, optimalNumMods);
}

//
// Helper method computing ring optimal parameters when rectangular-shaped modules used
//
std::pair<double, int> Ring::computeOptimalRingParametersRectangle(double moduleWidth, double maxRadius) {

  double delta            = phiOverlap()/maxRadius;
  double optimalAlpha     = 2*asin(moduleWidth/2. / maxRadius) - delta;
  double tentativeNumMods = 2*M_PI/optimalAlpha;
  int modsPerSlice        = ceil(tentativeNumMods/phiSegments());

  if ((modsPerSlice % 2) == 0 && m_requireOddModsPerSlice()) modsPerSlice++;
  int optimalNumMods = modsPerSlice * phiSegments();

  return std::make_pair(optimalAlpha, optimalNumMods);
}







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


  double zH2p = zOuter() - 0.5 * lengthEff * cos(tilt);
  double rH2p = outerRadius() + 0.5 * lengthEff * sin(tilt);
  //double zH2pp = zOuter() + 0.5 * lengthEff * cos(tilt);
  //double rH2pp = outerRadius() - 0.5 * lengthEff * sin(tilt);

  double zH2UP = zOuter() + 0.5 * dsDistance * sin(tilt);
  double rH2UP = outerRadius() + 0.5 * dsDistance * cos(tilt);
  double zH2pUP = zH2UP - 0.5 * lengthEff * cos(tilt);
  double rH2pUP = rH2UP + 0.5 * lengthEff * sin(tilt);
  //double zH2ppUP = zH2UP + 0.5 * lengthEff * cos(tilt);
  //double rH2ppUP = rH2UP - 0.5 * lengthEff * sin(tilt);

  double zH2DOWN = zOuter() - 0.5 * dsDistance * sin(tilt);
  double rH2DOWN = outerRadius() - 0.5 * dsDistance * cos(tilt);
  double zH2pDOWN = zH2DOWN - 0.5 * lengthEff * cos(tilt);
  double rH2pDOWN = rH2DOWN + 0.5 * lengthEff * sin(tilt);
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


  double zH1p = zInner() - 0.5 * lengthEff * cos(tilt);
  double rH1p = innerRadius() + 0.5 * lengthEff * sin(tilt);
  //double zH1pp = zInner() + 0.5 * lengthEff * cos(tilt);
  //double rH1pp = innerRadius() - 0.5 * lengthEff * sin(tilt);

  double zH1UP = zInner() + 0.5 * dsDistance * sin(tilt);
  double rH1UP = innerRadius() + 0.5 * dsDistance * cos(tilt);
  double zH1pUP = zH1UP - 0.5 * lengthEff * cos(tilt);
  double rH1pUP = rH1UP + 0.5 * lengthEff * sin(tilt);
  //double zH1ppUP = zH1UP + 0.5 * lengthEff * cos(tilt);
  //double rH1ppUP = rH1UP - 0.5 * lengthEff * sin(tilt);

  double zH1DOWN = zInner() - 0.5 * dsDistance * sin(tilt);
  double rH1DOWN = innerRadius() - 0.5 * dsDistance * cos(tilt);
  double zH1pDOWN = zH1DOWN - 0.5 * lengthEff * cos(tilt);
  double rH1pDOWN = rH1DOWN + 0.5 * lengthEff * sin(tilt);
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

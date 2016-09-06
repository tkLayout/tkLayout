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

  minZ.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& m : m_modules) min = MIN(min, m.minZ()); return min; });
  maxZ.setup([this]() { double max = -std::numeric_limits<double>::max();for (const auto& m : m_modules) max = MAX(max, m.maxZ()); return max; });
  minR.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& m : m_modules) min = MIN(min, m.minR()); return min; });
  maxR.setup([this]() { double max = 0;                                  for (const auto& m : m_modules) max = MAX(max, m.maxR()); return max; });

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

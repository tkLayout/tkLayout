#include "Layer.h"
#include "RodPair.h"
#include "MessageLogger.h"
#include "ConversionStation.h"

using std::string;

define_enum_strings(RadiusMode) = { "shrink", "enlarge", "fixed", "auto" };

//
// Constructor - parse geometry config file using boost property tree & read-in Layer parameters
//
Layer::Layer(int id, int barrelNumLayers, bool sameRods, bool barrelMinRFixed, bool barrelMaxRFixed, double barrelRotation,
             const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty) :
 minZ           (string("minZ")              ),
 maxZ           (string("maxZ")              ),
 minR           (string("minR")              ),
 maxR           (string("maxR")              ),
 buildNumModules(       "numModules"         , parsedOnly()),
 outerZ         (       "outerZ"             , parsedOnly()),
 radiusMode     (       "radiusMode"         , parsedAndChecked(), RadiusMode::AUTO),
 requestedAvgRadius(    "radius"             , parsedOnly()),
 avgBuildRadius (       "avgBuildRadius"     , parsedOnly()),
 sameParityRods (       "sameParityRods"     , parsedAndChecked(), true),
 layerRotation  (       "layerRotation"      , parsedOnly()      , 0.),
 tiltedLayerSpecFile(   "tiltedLayerSpecFile", parsedOnly()),
 m_smallDelta   (       "smallDelta"         , parsedAndChecked()),
 m_smallParity  (       "smallParity"        , parsedAndChecked(),-1),
 m_bigDelta     (       "bigDelta"           , parsedAndChecked()),
 m_bigParity    (       "bigParity"          , parsedOnly()      ,-1),
 m_phiOverlap   (       "phiOverlap"         , parsedAndChecked(), 1.),
 m_phiSegments  (       "phiSegments"        , parsedAndChecked(), 4),
 m_useMinMaxRCorrect(   "useMinMaxRCorrect"  , parsedAndChecked(), true),
 m_ringNode     (       "Ring"               , parsedOnly()),
 m_stationsNode (       "Station"            , parsedOnly()),
 m_materialObject(MaterialObject::LAYER),
 m_flangeConversionStation(nullptr),
 m_sameRods(sameRods)
{

  // Set the geometry config parameters
  this->myid(id);
  this->store(treeProperty);
  if (nodeProperty.count(id) > 0) this->store(nodeProperty.at(id));

  // If required the same rods, even/odd rods need to be built starting with the same sign of smallDelta
  if (m_sameRods) sameParityRods(true);

  // Set rotation
  this->layerRotation(layerRotation()+barrelRotation);

  // Set radius mode if not default value set or not defined from barrel level (for backwards 
  // compatibility fix to barrel limits can be switched off using useMinMaxRCorrect variable)
  if ( (myid()==1)               && barrelMinRFixed && m_useMinMaxRCorrect()) this->radiusMode(RadiusMode::FIXED);
  if ( (myid()==barrelNumLayers) && barrelMaxRFixed && m_useMinMaxRCorrect()) this->radiusMode(RadiusMode::FIXED);
}

//
// Limit layer geometry by eta cut
//
void Layer::cutAtEta(double eta)
{
  for (auto& r : m_rods) r.cutAtEta(eta);
  m_rods.erase_if([](const RodPair& r) { return r.numModules() == 0; }); // get rid of rods which have been completely pruned
}

//
// Cross-check parameters provided from geometry configuration file
//
void Layer::check()
{
  PropertyObject::check();

  if      (buildNumModules()>0  &&  outerZ.state()) throw PathfulException("Only one between numModules and outerZ can be specified");
  else if (buildNumModules()==0 && !outerZ.state()) throw PathfulException("At least one between numModules and outerZ must be specified");

  if (m_bigDelta()  <0)            throw PathfulException("Big delta parameter must be positive!");
  if (m_smallDelta()<0)            throw PathfulException("Small delta parameter must be positive!");
  if (m_bigDelta()<m_smallDelta()) throw PathfulException("Big delta parameter is expected to be bigger in size than small delta parameter!");
}

//
// Build recursively individual subdetector systems: rods -> modules
//
void Layer::build(int barrelNumLayers, double barrelMinR, double barrelMaxR) {

  try {
    m_materialObject.store(propertyTree());
    m_materialObject.build();

    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    // Build either straight or tilted geometry
    buildStraight(barrelNumLayers, barrelMinR, barrelMaxR);
    //if (tiltedLayerSpecFile().empty()) buildStraight(barrelNumLayers, barrelMinR, barrelMaxR);
    //else                               buildTilted();

    // Build conversion stations of first or second order
    for (auto& currentStationNode : m_stationsNode) {

      ConversionStation* conversionStation = new ConversionStation();
      conversionStation->store(currentStationNode.second);
      conversionStation->check();
      conversionStation->build();

      if (conversionStation->stationType() == ConversionStation::Type::FLANGE) {

        // Take only first defined flange station
        if (m_flangeConversionStation == nullptr) m_flangeConversionStation = conversionStation;
      }
      else if (conversionStation->stationType() == ConversionStation::Type::SECOND) m_secondConversionStations.push_back(conversionStation);
    }

    cleanup();
    builtok(true);

  }
  catch (PathfulException& pe) {
    pe.pushPath(fullid(*this));
    throw;
  }
}

//
// If straight layer required, build() method internally calls buildStraight()
//
void Layer::buildStraight(int barrelNumLayers, double barrelMinR, double barrelMaxR) {

  //
  // Optimization algorithm

  // Calculate minimum/maximum R boundary taking into account modules layout (+/-bigDelta +-smallDelta +-dsDistance)
  ReadonlyProperty<double, NoDefault> moduleWidth( "width", parsedOnly());
  evaluateProperty(moduleWidth);
  ReadonlyProperty<double, Default> maxDsDistance( "dsDistance", parsedOnly(), 0.0);
  evaluateProperty(maxDsDistance);
  ReadonlyProperty<double, Default> sensorThickness( "sensorThickness", parsedOnly(), 0.1);
  evaluateProperty(sensorThickness);

  double updatedMinR = barrelMinR + m_bigDelta() + m_smallDelta() + maxDsDistance()/2. + sensorThickness()/2.;
  double updatedMaxR = barrelMaxR - m_bigDelta() - m_smallDelta() - maxDsDistance()/2. - sensorThickness()/2.;
         updatedMaxR = sqrt(updatedMaxR*updatedMaxR - moduleWidth()/2.*moduleWidth()/2.); // Need to take the outer envelope, i.e. the module width into account

  // For backwards compatibility with older versions of sofware
  if (!m_useMinMaxRCorrect()) {
    updatedMinR = barrelMinR;
    updatedMaxR = barrelMaxR;
  }

  // Calculate expected radius
  double layerRadius = 0;
  if      ((myid() == 1)              ) layerRadius = updatedMinR;
  else if ((myid() == barrelNumLayers)) layerRadius = updatedMaxR;
  else                                  layerRadius = updatedMinR + (updatedMaxR - updatedMinR)/(barrelNumLayers-1)*(myid()-1);

  // Requested radius out of barrel region -> will use automated algorithm instead
  if (requestedAvgRadius.state() && (requestedAvgRadius()< updatedMinR || requestedAvgRadius()>updatedMaxR)) {

    std::ostringstream message;
    message << "Layer requested to be built on R: " << requestedAvgRadius() <<", which is out of defined barrel region (taking intou account bigDelta & small delta parameters): [";
    message << updatedMinR << ";" << barrelMaxR << "] -> ignoring, will use automated algorithm instead!";
    logWARNING(message);
    this->requestedAvgRadius(layerRadius);
  }
  if (!requestedAvgRadius.state()) requestedAvgRadius(layerRadius);

  //
  // Calculate optimal layer parameters: pessimistic scenario is with bigDelta + smallDelta + dsDistance versus -bigDelta + smallDelta + dsDistance (with -smallDelta it's always better)
  float halfWidthWoOverlap = moduleWidth()/2 - m_phiOverlap()/2;
  float gamma              = atan(halfWidthWoOverlap/(requestedAvgRadius() + m_bigDelta() + m_smallDelta() + maxDsDistance()/2)) + atan(halfWidthWoOverlap/(requestedAvgRadius() - m_bigDelta() + m_smallDelta() + maxDsDistance()/2));
  float modsPerSegment     = 2*M_PI/(gamma * m_phiSegments());

  float optimalRadius;
  int   optimalModsPerSegment;

  switch (radiusMode()) {
  case SHRINK:
    optimalModsPerSegment = floor(modsPerSegment);
    optimalRadius         = calculateOptimalRadius(optimalModsPerSegment*m_phiSegments(), m_bigDelta(), m_smallDelta(), maxDsDistance(), moduleWidth(), m_phiOverlap());
    break;
  case ENLARGE:
    optimalModsPerSegment = ceil(modsPerSegment);
    optimalRadius         = calculateOptimalRadius(optimalModsPerSegment*m_phiSegments(), m_bigDelta(), m_smallDelta(), maxDsDistance(), moduleWidth(), m_phiOverlap());
    break;
  case FIXED:
    optimalModsPerSegment = ceil(modsPerSegment);
    optimalRadius         = requestedAvgRadius();
    break;
  case AUTO: {
    int modsPerSegLo = floor(modsPerSegment);
    int modsPerSegHi = ceil(modsPerSegment);
    float radiusLo   = calculateOptimalRadius(modsPerSegLo*m_phiSegments(), m_bigDelta(), m_smallDelta(), maxDsDistance(), moduleWidth(), m_phiOverlap());
    float radiusHi   = calculateOptimalRadius(modsPerSegHi*m_phiSegments(), m_bigDelta(), m_smallDelta(), maxDsDistance(), moduleWidth(), m_phiOverlap());

    if (fabs(radiusHi - requestedAvgRadius()) < fabs(radiusLo - requestedAvgRadius())) {
      optimalRadius         = radiusHi;
      optimalModsPerSegment = modsPerSegHi;
    } else {
      optimalRadius         = radiusLo;
      optimalModsPerSegment = modsPerSegLo;
    }

    break;
  }
  default:
    throw PathfulException("Invalid value for enum radiusMode");
  }

  //
  // Set rod properties
  avgBuildRadius(optimalRadius);

  int   optimalNumRods = optimalModsPerSegment*m_phiSegments();
  float rodPhiRotation = 2*M_PI/optimalNumRods;

  // Rod pair corresponds to an odd/even rod (shifted by bigDelta & rotated by rod phi rotation angle - given by total number of rods in a layer).
  // The pair stands for a pair of modules in positive/negative Z. Odd/even rod is prototyped as first/second, others are cloned to speed-up building.
  RodPairStraight* firstRod  = nullptr;
  RodPairStraight* secondRod = nullptr;

  // Use these extreme values to calculate optimal positions in Z to have full eta coverage
  // For same rods extremes given (minBarrelR/maxBarrelR) +-bigDelta
  double minRadius = 0.;
  double maxRadius = 0.;
  if (m_sameRods) {

    minRadius = updatedMinR - m_bigDelta();
    maxRadius = updatedMaxR + m_bigDelta();
  }
  // Not same rods -> extreme given by optimal radius +-bigDelta
  else {
    minRadius = optimalRadius - m_bigDelta();
    maxRadius = optimalRadius + m_bigDelta();
  }

  for (int i=1; i<=optimalNumRods; i++) {

    // Set current rod rotation & big parity value
    double rotation    = rodPhiRotation*(i-1) + layerRotation();

    // Alternate bigParity for even versus odd rods
    int    bigParity   = (i%2) ? m_bigParity() : -m_bigParity();

    // Use the same small parity in even rods as for odd rods in case sameParityRods is required
    int    smallParity = (i%2) ? (m_smallParity()) : (sameParityRods() ? m_smallParity() : -m_smallParity());

    // Prototypes of odd rods
    if (i==1) {

      RodPairStraight* rod = nullptr;

      if (buildNumModules() > 0) rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, m_bigDelta(), bigParity, m_smallDelta(), smallParity, buildNumModules(), propertyTree());
      else                       rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, m_bigDelta(), bigParity, m_smallDelta(), smallParity, outerZ(), propertyTree());
      rod->build();
      firstRod = rod;

      m_rods.push_back(rod);
    }
    // Prototype of even rods
    else if (i==2) {

      // If same rods or sameParityRods required, each even/odd rod are the same -> useful from engineering point of view
      if (m_sameRods || sameParityRods()) {

        RodPair* rod = GeometryFactory::clone<RodPairStraight>(*firstRod);
        rod->buildClone(2, 2*m_bigDelta()*bigParity, rodPhiRotation);

        m_rods.push_back(rod);
      }
      else {

        RodPairStraight* rod = nullptr;

        if (buildNumModules() > 0) rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, m_bigDelta(), bigParity, m_smallDelta(), smallParity, buildNumModules(), propertyTree());
        else                       rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, m_bigDelta(), bigParity, m_smallDelta(), smallParity, outerZ(), propertyTree());
        rod->build();
        secondRod = rod;

        m_rods.push_back(rod);
      }

    }
    // Clone prototypes to speed-up building -> need for buildClone() call to update new id and rotate cloned rod by respective angle
    else {

      RodPair*  rod = nullptr;
      double shiftR = 0.0;

      // Build odd rod -> Rotation with respect to first rod, no shift in R
      if ((i%2)==1) {

        rod = GeometryFactory::clone<RodPairStraight>(*firstRod);
        rod->buildClone(i, shiftR, rodPhiRotation*(i-1));

        m_rods.push_back(rod);

      }
      // Build even rod
      else {

        // Clone even to even -> Rotation with respect to second rod, no shift in R
        if (!(m_sameRods || sameParityRods()) ) {

          rod = GeometryFactory::clone<RodPairStraight>(*secondRod);
          rod->buildClone(i, shiftR, rodPhiRotation*(i-2));

        }
        // Clone odd to even -> Rotation with respect to first rod, shift by 2*bigDelta
        else  {
          rod = GeometryFactory::clone<RodPairStraight>(*firstRod);
          rod->buildClone(i, 2*m_bigDelta()*bigParity, rodPhiRotation*(i-1));
        }

        m_rods.push_back(rod);
      }
    }
  } // For NumRods
}

//
// Setup: link lambda functions to various layer related properties (use setup functions for ReadOnly Computable properties)
//
void Layer::setup()
{
  maxZ.setup([&]() { return m_rods.front().maxZ(); });
  minZ.setup([&]() { return m_rods.front().minZ(); });
  maxR.setup([&]() { double max = 0;                                  for (const auto& r : m_rods) { max = MAX(max, r.maxR()); } return max; });
  minR.setup([&]() { double min = std::numeric_limits<double>::max(); for (const auto& r : m_rods) { min = MIN(min, r.minR()); } return min; });
}

//
// GeometryVisitor pattern -> layer visitable
//
void Layer::accept(GeometryVisitor& v)
{
  v.visit(*this);
  for (auto& r : m_rods) { r.accept(v); }
}

//
// GeometryVisitor pattern -> layer visitable (const. option)
//
void Layer::accept(ConstGeometryVisitor& v) const
{
  v.visit(*this);
  for (const auto& r : m_rods) { r.accept(v); }
}

//
// Helper function calculating optimal layer radius for straight option
//
double Layer::calculateOptimalRadius(int numRods, double bigDelta, double smallDelta, double dsDistance, double moduleWidth, double overlap)
{
  double f = (moduleWidth/2) - (overlap/2);
  double R = +bigDelta + smallDelta + dsDistance/2; // One extreme: +bigDelta + smallDelta & +dsDistance (more pessimistic than -smallDelta ...)
  double S = -bigDelta + smallDelta + dsDistance/2; // Another extreme: -bigDelta +smallDelta & +dsDistance (more pessimistic than -smallDelta ...)
  double T = tan(2*M_PI/numRods);

  // Quadratic equation coefficients: 2pi/nRods = atan(f/(r+R)) + atan(f/(r+S)) = atan[(f/(r+R) + f/(r+S))/(1-f*f/(r+R)/r+S)]
  double a = T;
  double b = R*T + S*T - 2*f;
  double c = R*S*T - R*f - S*f - T*f*f;

  return (-b + sqrt(b*b - 4*a*c))/(2*a);
}

//
// If tilted layer required, build() method internally calls buildTilted()
//
//void Layer::buildTilted() {
//  std::ifstream ifs(tiltedLayerSpecFile());
//  if (ifs.fail()) throw PathfulException("Cannot open tilted modules spec file \"" + tiltedLayerSpecFile() + "\"");
//
//  string line;
//  vector<TiltedModuleSpecs> tmspecs1, tmspecs2;
//
//  int numRods = 0;
//
//  while(getline(ifs, line).good()) {
//    if (line.empty()) continue;
//    auto tokens = split<double>(line, " ", false);
//    if (tokens.size() < 7) { logERROR("Failed parsing tilted barrel line: " + line); continue; };
//    TiltedModuleSpecs t1{ tokens[0], tokens[1], tokens[2]*M_PI/180. };
//    TiltedModuleSpecs t2{ tokens[3], tokens[4], tokens[5]*M_PI/180. };
//    if (t1.valid()) tmspecs1.push_back(t1);
//    if (t2.valid()) tmspecs2.push_back(t2);
//    numRods = tokens[6]; // this assumes every row of the spec file has the same value for the last column (num rods in phi)
//  }
//  ifs.close();
//
//  buildNumModules(tmspecs1.size());
//
//  RodTemplate rodTemplate = makeRodTemplate();
//
//  float rodPhiRotation = 2*M_PI/numRods;
//
//  TiltedRodPair* first = GeometryFactory::make<TiltedRodPair>();
//  first->myid(1);
//  first->store(propertyTree());
//  first->build(rodTemplate, tmspecs1);
//  m_rods.push_back(first);
//
//  TiltedRodPair* second = GeometryFactory::make<TiltedRodPair>();
//  second->myid(2);
//  second->store(propertyTree());
//  second->build(rodTemplate, tmspecs2);
//  second->rotateZ(rodPhiRotation);
//  m_rods.push_back(second);
//
//  for (int i = 2; i < numRods; i++) {
//    RodPair* rod = i%2 ? GeometryFactory::clone(*second) : GeometryFactory::clone(*first); // clone rods
//    rod->myid(i+1);
//    rod->rotateZ(rodPhiRotation*(i%2 ? i-1 : i));
//    m_rods.push_back(rod);
//  }
//
//  // computing the layer's place radius as the average of all the modules' radii
//  avgBuildRadius(std::accumulate(tmspecs1.begin(), tmspecs1.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; }));
//  avgBuildRadius(avgBuildRadius() + std::accumulate(tmspecs2.begin(), tmspecs2.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; }));
//  avgBuildRadius(avgBuildRadius() / (tmspecs1.size() + tmspecs2.size()));
//
//}






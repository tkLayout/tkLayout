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
 minZ           (        string("minZ")                        ),
 maxZ           (        string("maxZ")                        ),
 minR           (        string("minR")                        ),
 maxR           (        string("maxR")                        ),
 minRAllMat     (        string("minRAllMat")                  ),
 maxRAllMat     (        string("maxRAllMat")                  ),
 buildNumModules(               "numModules"                   , parsedOnly()),
 buildNumModulesFlat(           "numModulesFlat"               , parsedOnly()),
 buildNumModulesTilted(         "numModulesTilted"             , parsedOnly()),
 outerZ         (               "outerZ"                       , parsedOnly()),
 radiusMode     (               "radiusMode"                   , parsedAndChecked(), RadiusMode::AUTO),
 requestedAvgRadius(            "radius"                       , parsedOnly()),
 avgBuildRadius (               "avgBuildRadius"               , parsedOnly()),
 sameParityRods (               "sameParityRods"               , parsedAndChecked(), true),
 layerRotation  (               "layerRotation"                , parsedOnly()      , 0.),
 isTilted       (               "isTilted"                     , parsedOnly(), false),
 isTiltedAuto   (               "isTiltedAuto"                 , parsedOnly()),
 tiltedLayerSpecFile(           "tiltedLayerSpecFile"          , parsedOnly()),
 smallDelta     (               "smallDelta"                   , parsedAndChecked()),
 m_smallParity  (               "smallParity"                  , parsedAndChecked(),-1),
 bigDelta       (               "bigDelta"                     , parsedAndChecked()),
 m_bigParity    (               "bigParity"                    , parsedOnly()      ,-1),
 phiOverlap     (               "phiOverlap"                   , parsedOnly()),   // used to be default value = 1.
 phiSegments    (               "phiSegments"                  , parsedOnly()),   // used to be default value = 4.
 m_useMinMaxRCorrect(           "useMinMaxRCorrect"            , parsedAndChecked(), true),
 numberRods     (               "numberRods"                   , parsedOnly()),
 m_ringNode     (               "Ring"                         , parsedOnly()),
 m_stationsNode (               "Station"                      , parsedOnly()),
 flatPhiOverlapSmallDeltaMinus( "flatPhiOverlapSmallDeltaMinus", parsedOnly()),
 flatPhiOverlapSmallDeltaPlus(  "flatPhiOverlapSmallDeltaPlus" , parsedOnly()),
 flatThetaEnd(                  "flatThetaEnd"                 , parsedOnly()),
 flatREndInner(                 "flatREndInner"                , parsedOnly()),
 flatREndOuter(                 "flatREndOuter"                , parsedOnly()),
 flatEndModCentreZ(             "flatEndModCenterZ"            , parsedOnly()),
 flatZEnd(                      "flatZEnd"                     , parsedOnly()),
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

  // Set radius mode if not default value set or not defined from barrel level
  if ( (myid()==1)               && barrelMinRFixed) this->radiusMode(RadiusMode::FIXED);
  if ( (myid()==barrelNumLayers) && barrelMaxRFixed) this->radiusMode(RadiusMode::FIXED);
}

//
// Limit layer geometry by eta cut
//
void Layer::cutAtEta(double eta)
{
  for (auto& r : m_flatRods)   r.cutAtEta(eta);
  for (auto& r : m_tiltedRods) r.cutAtEta(eta);

  // Get rid of rods which have been completely pruned
  m_flatRods.erase_if([](const RodPair& r)   { return r.numModules() == 0; });
  m_tiltedRods.erase_if([](const RodPair& r) { return r.numModules() == 0; });
}

//
// Cross-check parameters provided from geometry configuration file
//
void Layer::check()
{
  PropertyObject::check();

  // Non-tilted geometry
  if (!isTilted()) {
    if (buildNumModules()>0  && outerZ.state()) throw PathfulException("Non-tilted geom.: Either numModules or outerZ parameter may be specified. Not both.");
    if (buildNumModules()==0 &&!outerZ.state()) throw PathfulException("Non-tilted geom.: Either numModules or outerZ parameter must be specified!");
    if (!phiOverlap.state())                    throw PathfulException("Flat layer: phiOverlap must be specified.");
    if (!phiSegments.state())                   throw PathfulException("Flat layer: phiSegments must be specified.");

    if (isTiltedAuto.state())      logERROR("Layer " + std::to_string(myid()) + ": doesn't make sense to specify isTiltedAuto. Not being used.");
    if (buildNumModulesFlat()>0)   logERROR("Layer " + std::to_string(myid()) + ": doesn't make sense to specify numModulesFlat. Not being used.");
    if (buildNumModulesTilted()>0) logERROR("Layer " + std::to_string(myid()) + ": doesn't make sense to specify numModulesTilted. Not being used.");

    // Set number flat modules & total number of all modules to the same number (compatibility reasons)
    if (buildNumModules()>0) buildNumModulesFlat(buildNumModules());
  }
  // Tilted geometry
  else {

    //TODO: Check phiOverlap for tilted
    //if (phiOverlap.state())    throw PathfulException("Tilted layer: phiOverlap should not be specified.");
    if (phiSegments.state())   throw PathfulException("Tilted layer: phiSegments should not be specified.");

    // Automatic tilt algorithm applied
    if (isTiltedAuto()) {

      if (buildNumModulesFlat()<0)   throw PathfulException("Tilted layer with automatic placement: numModulesFlat must be specified.");
      if (buildNumModulesTilted()<0) throw PathfulException("Tilted layer with automatic placement: numModulesTilted must be specified.");
      if (outerZ.state())            throw PathfulException("Tilted layer with automatic placement: outerZ can't be defined, expect number of flat & tilted modules.");

      // Check total number of all modules
      if (buildNumModules()>0 && buildNumModules()!=(buildNumModulesFlat() + buildNumModulesTilted())) {
          throw PathfulException("Tilted layer : numModules != (numModulesFlat + numModulesTilted). Specify numModulesFlat and numModulesTilted only, don't specify numModules.");
      }
      else buildNumModules(buildNumModulesFlat() + buildNumModulesTilted());
    }
    // Config file for tilting used instead -> NOT IMPLEMENTED anymore
    else throw PathfulException("Tilted geometry built on configuration file not implemented. Use automatic build instead!");

    // Check that outerZ not defined for tilted part
    if (outerZ.state()) logERROR("Tilted layer: outerZ was specified. Routing of services will be forced to be at Z = outerZ.");
  }

  if (buildNumModules()<0)       throw PathfulException("Number of modules must be positive or zero!");
  if (buildNumModulesFlat()<0)   throw PathfulException("Number of modules in flat part must be positive or zero!");
  if (buildNumModulesTilted()<0) throw PathfulException("Number of modules in tilted part must be positive or zero!");
  if (bigDelta()  <0)            throw PathfulException("Big delta parameter must be positive!");
  if (smallDelta()<0)            throw PathfulException("Small delta parameter must be positive!");
  if (bigDelta()<smallDelta())   throw PathfulException("Big delta parameter is expected to be bigger in size than small delta parameter!");
}


//
// Setup: link lambda functions to various layer related properties (use setup functions for ReadOnly Computable properties)
//
void Layer::setup()
{
  maxFlatZ.setup([&]()       { if (m_flatRods.size()>0) return m_flatRods.front().maxZ(); else return 0.0; });
  minFlatZ.setup([&]()       { if (m_flatRods.size()>0) return m_flatRods.front().minZ(); else return 0.0; });
  maxFlatR.setup([&]()       { double max = 0;                                  for (const auto& r : m_flatRods) { max = MAX(max, r.maxR()); } return max; });
  minFlatR.setup([&]()       { double min = std::numeric_limits<double>::max(); for (const auto& r : m_flatRods) { min = MIN(min, r.minR()); } return min; });
  maxFlatRAllMat.setup([&]() { double max = 0;                                  for (const auto& r : m_flatRods) { max = MAX(max, r.maxRAllMat()); } return max; });
  minFlatRAllMat.setup([&]() { double min = std::numeric_limits<double>::max(); for (const auto& r : m_flatRods) { min = MIN(min, r.minRAllMat()); } return min; });

  maxTiltedZ.setup([&]()       { double max = 0;                                  for (const auto& r : m_tiltedRods) { max = MAX(max, r.maxZ()); } return max; });
  minTiltedZ.setup([&]()       { double min = std::numeric_limits<double>::max(); for (const auto& r : m_tiltedRods) { min = MIN(min, r.minZ()); } return min; });
  maxTiltedR.setup([&]()       { double max = 0;                                  for (const auto& r : m_tiltedRods) { max = MAX(max, r.maxR()); } return max; });
  minTiltedR.setup([&]()       { double min = std::numeric_limits<double>::max(); for (const auto& r : m_tiltedRods) { min = MIN(min, r.minR()); } return min; });
  maxTiltedRAllMat.setup([&]() { double max = 0;                                  for (const auto& r : m_tiltedRods) { max = MAX(max, r.maxRAllMat()); } return max; });
  minTiltedRAllMat.setup([&]() { double min = std::numeric_limits<double>::max(); for (const auto& r : m_tiltedRods) { min = MIN(min, r.minRAllMat()); } return min; });

  maxZ.setup([&]()       { double max = MAX(maxFlatZ(), maxTiltedZ()); return max; });
  minZ.setup([&]()       { double min = MIN(minFlatZ(), minTiltedZ()); return min; });
  maxR.setup([&]()       { double max = MAX(maxFlatR(), maxTiltedR()); return max; });
  minR.setup([&]()       { double min = MIN(minFlatR(), minTiltedR()); return min; });
  maxRAllMat.setup([&]() { double max = MAX(maxFlatRAllMat(), maxTiltedRAllMat()); return max; });
  minRAllMat.setup([&]() { double min = MIN(minFlatRAllMat(), minTiltedRAllMat()); return min; });
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

    // Build either flat layer or flat part of the tilted layer
    buildStraight(barrelNumLayers, barrelMinR, barrelMaxR);

    // Build tilted part of tilted layer
    if (isTilted()) {

      // Calculate geometry properties of straight part needed to build the tilted part
      calculateStraightProperties();

      // Build tilted part
      buildTilted();
    }

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
// If straight layer required (or tilted layer with straight part), build() method internally calls buildStraight()
//
void Layer::buildStraight(int barrelNumLayers, double barrelMinR, double barrelMaxR) {

  //
  // Optimization algorithm

  // Calculate minimum/maximum R boundary taking into account modules layout (+/-bigDelta +-smallDelta +-dsDistance)
  // TODO: Still not taking into account all materials, just the active thickness, not the passive
  ReadonlyProperty<double, NoDefault> moduleWidth( "width", parsedOnly());
  evaluateProperty(moduleWidth);
  ReadonlyProperty<double, Default> maxDsDistance( "dsDistance", parsedOnly(), 0.0);
  evaluateProperty(maxDsDistance);
  ReadonlyProperty<double, Default> sensorThickness( "sensorThickness", parsedOnly(), 0.1);
  evaluateProperty(sensorThickness);

  double updatedMinR = barrelMinR + bigDelta() + smallDelta() + maxDsDistance()/2. + sensorThickness()/2.;
  double updatedMaxR = barrelMaxR - bigDelta() - smallDelta() - maxDsDistance()/2. - sensorThickness()/2.;
         updatedMaxR = sqrt(updatedMaxR*updatedMaxR - moduleWidth()/2.*moduleWidth()/2.); // Need to take the outer envelope, i.e. the module width into account

  // For backwards compatibility with older versions of sofware
  if (!m_useMinMaxRCorrect()) {
    updatedMinR = barrelMinR;
    updatedMaxR = barrelMaxR;
  }

  float optimalRadius;
  int numLayerRods;
  if (!isTilted()) {
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
    float halfWidthWoOverlap = moduleWidth()/2 - phiOverlap()/2;
    float gamma              = atan(halfWidthWoOverlap/(requestedAvgRadius() + bigDelta() + smallDelta() + maxDsDistance()/2)) + atan(halfWidthWoOverlap/(requestedAvgRadius() - bigDelta() + smallDelta() + maxDsDistance()/2));
    float modsPerSegment     = 2*M_PI/(gamma * phiSegments());

    int   optimalModsPerSegment;

    switch (radiusMode()) {
    case SHRINK:
      optimalModsPerSegment = floor(modsPerSegment);
      optimalRadius         = calculateOptimalRadius(optimalModsPerSegment*phiSegments(), bigDelta(), smallDelta(), maxDsDistance(), moduleWidth(), phiOverlap());
      break;
    case ENLARGE:
      optimalModsPerSegment = ceil(modsPerSegment);
      optimalRadius         = calculateOptimalRadius(optimalModsPerSegment*phiSegments(), bigDelta(), smallDelta(), maxDsDistance(), moduleWidth(), phiOverlap());
      break;
    case FIXED:
      optimalModsPerSegment = ceil(modsPerSegment);
      optimalRadius         = requestedAvgRadius();
      break;
    case AUTO: {
      int modsPerSegLo = floor(modsPerSegment);
      int modsPerSegHi = ceil(modsPerSegment);
      float radiusLo   = calculateOptimalRadius(modsPerSegLo*phiSegments(), bigDelta(), smallDelta(), maxDsDistance(), moduleWidth(), phiOverlap());
      float radiusHi   = calculateOptimalRadius(modsPerSegHi*phiSegments(), bigDelta(), smallDelta(), maxDsDistance(), moduleWidth(), phiOverlap());

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
    avgBuildRadiusFlat(optimalRadius);
    avgBuildRadiusTilted(optimalRadius);
    numLayerRods = optimalModsPerSegment*phiSegments();
  }

  // The tilted part is not being optimized using this algorithm
  else {
    optimalRadius = requestedAvgRadius();
    avgBuildRadius(requestedAvgRadius());
    avgBuildRadiusFlat(requestedAvgRadius());
    avgBuildRadiusTilted(requestedAvgRadius());
    numLayerRods = numberRods();
  }

  float rodPhiRotation = 2. * M_PI / numLayerRods;

  // Rod pair corresponds to an odd/even rod (shifted by bigDelta & rotated by rod phi rotation angle - given by total number of rods in a layer).
  // The pair stands for a pair of modules in positive/negative Z. Odd/even rod is prototyped as first/second, others are cloned to speed-up building.
  RodPairStraight* firstRod  = nullptr;
  RodPairStraight* secondRod = nullptr;

  // Use these extreme values to calculate optimal positions in Z to have full eta coverage
  // For same rods extremes given (minBarrelR/maxBarrelR) +-bigDelta
  double minRadius = 0.;
  double maxRadius = 0.;
  if (m_sameRods) {

    minRadius = updatedMinR - bigDelta();
    maxRadius = updatedMaxR + bigDelta();
  }
  // Not same rods -> extreme given by optimal radius +-bigDelta
  else {
    minRadius = optimalRadius - bigDelta();
    maxRadius = optimalRadius + bigDelta();
  }

  for (int i=1; i<=numLayerRods; i++) {

    // Set current rod rotation & big parity value
    double rotation    = rodPhiRotation*(i-1) + layerRotation();

    // Alternate bigParity for even versus odd rods
    int    bigParity   = (i%2) ? m_bigParity() : -m_bigParity();

    // Use the same small parity in even rods as for odd rods in case sameParityRods is required
    int    smallParity;
    if (!isTilted()) smallParity = (i%2) ? (m_smallParity()) : (sameParityRods() ? m_smallParity() : -m_smallParity());
    else             smallParity = pow(-1, buildNumModulesFlat());  // if tilted, smallParity is forced so that central module is at +smallDelta

    // Prototypes of odd rods
    if (i==1) {

      RodPairStraight* rod = nullptr;

      if (buildNumModulesFlat()>0) rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, buildNumModulesFlat(), propertyTree());
      else                         rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, outerZ(), propertyTree());

      rod->build();
      
      firstRod = rod;
      
      m_flatRods.push_back(rod);
    }

    // Prototype of even rods
    else if (i==2) {

      // If same rods or sameParityRods required, each even/odd rod are the same -> useful from engineering point of view
      if (m_sameRods || sameParityRods()) {

        RodPairStraight* rod = GeometryFactory::clone<RodPairStraight>(*firstRod);
        rod->buildClone(2, 2*bigDelta()*bigParity, rodPhiRotation);

        m_flatRods.push_back(rod);
      }
      else {

        RodPairStraight* rod = nullptr;

        if (buildNumModulesFlat()>0) rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, buildNumModulesFlat(), propertyTree());
        else                         rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, outerZ(), propertyTree());
        rod->build();
        secondRod = rod;

        m_flatRods.push_back(rod);
      }

    }
    // Clone prototypes to speed-up building -> need for buildClone() call to update new id and rotate cloned rod by respective angle
    else {

      RodPairStraight*  rod = nullptr;
      double shiftR = 0.0;

      // Build odd rod -> Rotation with respect to first rod, no shift in R
      if ((i%2)==1) {

        rod = GeometryFactory::clone<RodPairStraight>(*firstRod);
        rod->buildClone(i, shiftR, rodPhiRotation*(i-1));

        m_flatRods.push_back(rod);
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
          rod->buildClone(i, 2*bigDelta()*bigParity, rodPhiRotation*(i-1));
        }

        m_flatRods.push_back(rod);
      }
    }
  } // For NumRods
}

//
// If tilted layer required, build() method internally calls buildTilted()
//
void Layer::buildTilted() {

  //
  // Build tilted rods in a template
  RodTemplate rodTemplate;

  //
  // Build tilted rings -> tilted rings duplicate structure of tilted rods (modules are built only once for tilted rods)
  TiltedRing* prevRing = nullptr;
  for (auto iRingTot = (buildNumModulesFlat() + 1); iRingTot < (buildNumModules() + 1); iRingTot++) {

    auto iRing = iRingTot - (buildNumModulesFlat() + 1);

    // Create tilted ring
    TiltedRing* tiltedRing = GeometryFactory::make<TiltedRing>(iRingTot, numberRods(), m_ringNode, propertyTree());

    // Theta, at which the previous ring (or flat part) ends-up
    double lastThetaEnd;
    if (iRing==0) lastThetaEnd = flatThetaEnd();
    else          lastThetaEnd = m_tiltedRings[iRing-1].thetaEndOuter_REAL();

    // Build the tilted ring
    tiltedRing->build(lastThetaEnd);

    // Calculate geometry properties with respect to previous tilted ring or flat part (if first tilted ring)
    tiltedRing->calculateGeometryProperties(prevRing, flatREndInner(), flatREndOuter(), flatEndModCentreZ(), flatZEnd());

    // Update average radius of tilted part
    if (iRing==0) avgBuildRadiusTilted((tiltedRing->innerRadius()+tiltedRing->outerRadius())/2.);
    else          avgBuildRadiusTilted((avgBuildRadiusTilted()*iRing+(tiltedRing->innerRadius()+tiltedRing->outerRadius())/2.)/(iRing+1));

    // Create template rod of tilted modules for each ring
    BarrelModule* module = GeometryFactory::make<BarrelModule>(iRingTot, GeometryFactory::make<RectangularModule>(), m_ringNode, propertyTree());
    module->build();
    rodTemplate.push_back(module);

    // Save the ring
    prevRing = tiltedRing;
    m_tiltedRings.push_back(tiltedRing);
  }

  //
  // Rod template built, i.e. corresponding to 1 ring (module), duplicate & build all ring modules along Z as a tilted rod

  // Build at phi=0+totalLayerRotation -> put at inner radius
  double rotation = 0. + layerRotation();

  bool isOuterRadiusRod = (m_bigParity()==1) ? true : false;
  RodPairTilted* first = GeometryFactory::make<RodPairTilted>(1, rotation, isOuterRadiusRod, propertyTree());
  first->build(rodTemplate, m_tiltedRings, isOuterRadiusRod ? false : true); // If put at inner radius -> flip
  m_tiltedRods.push_back(first);

  // Build at phi=rodPhiRotation+totalLayerRotation -> put at outer radius
  float rodPhiRotation = 2*M_PI/numberRods();
  rotation             = rodPhiRotation + layerRotation();
  isOuterRadiusRod     = (m_bigParity()==1) ? false : true;

  RodPairTilted* second = GeometryFactory::make<RodPairTilted>(2, rotation, isOuterRadiusRod, propertyTree());
  second->build(rodTemplate, m_tiltedRings, isOuterRadiusRod ? false : true); // If put at outer radius -> don't flip
  m_tiltedRods.push_back(second);

  for (auto iRod = 2; iRod<numberRods(); iRod++) {
    RodPair* rod = iRod%2 ? GeometryFactory::clone(*second) : GeometryFactory::clone(*first); // clone rods
    rod->myid(iRod+1);
    rod->rotateZ(rodPhiRotation*(iRod%2 ? iRod-1 : iRod));
    m_tiltedRods.push_back(rod);
  }

  // Compute the place radius for the whole layer as the weighted avarge of tilted & flat part
  int numTilt = buildNumModulesTilted();
  int numFlat = hasStraightCentralModule() ? buildNumModulesFlat()-1 : buildNumModulesFlat(); // Get average radius between inner and outer rod, independently on whether we have even or odd number of modules
  avgBuildRadius( (avgBuildRadiusTilted()*numTilt+avgBuildRadiusFlat()*numFlat) / (numTilt+numFlat) );
}

//
// Flat part has central module -> build procedure started with module positioned at Z=0 (not by edge)
//
bool Layer::hasStraightCentralModule() const {

  bool hasCentralModule = false;

  if (m_flatRods.size() >= 1) {

    const auto& flatRod = dynamic_cast<const RodPairStraight&>(m_flatRods[0]);
    hasCentralModule = flatRod.hasCentralModule();
  }

  return hasCentralModule;
}

//
// Calculate geometry properties of straight part needed to build the tilted part
//
void Layer::calculateStraightProperties() {

  double rInnerMin = std::numeric_limits<double>::max();
  double rOuterMin = std::numeric_limits<double>::max();
  double rInnerMax = 0;
  double rOuterMax = 0;

  //
  // Prepare information about flat part of the layer
  if (buildNumModulesFlat() != 0) {

    if (m_flatRods.size() >= 2) {

      // Calculate min & max radial position of inner & outer rod
      const auto& flatRod1      = dynamic_cast<const RodPairStraight&>(m_flatRods[0]);
      const auto& zPlusModules1 = flatRod1.modules().first;

      int iMod = 0;
      for (const auto& m : zPlusModules1) {

        iMod++;
        (m_bigParity() > 0 ? rOuterMin = MIN(rOuterMin, m.center().Rho() + 0.5*m.dsDistance()) : rInnerMin = MIN(rInnerMin, m.center().Rho() + 0.5*m.dsDistance()));
        (m_bigParity() > 0 ? rOuterMax = MAX(rOuterMax, m.center().Rho() + 0.5*m.dsDistance()) : rInnerMax = MAX(rInnerMax, m.center().Rho() + 0.5*m.dsDistance()));
        if (iMod>=2) break; // Sufficient to have a look at +-smallDelta (first 2 modules)
      }

      const auto& flatRod2      = dynamic_cast<const RodPairStraight&>(m_flatRods[1]);
      const auto& zPlusModules2 = flatRod2.modules().first;

      iMod = 0;
      for (const auto& m : zPlusModules2) {

        iMod++;
        (m_bigParity() > 0 ? rInnerMin = MIN(rInnerMin, m.center().Rho() + 0.5*m.dsDistance()) : rOuterMin = MIN(rOuterMin, m.center().Rho() + 0.5*m.dsDistance()));
        (m_bigParity() > 0 ? rInnerMax = MAX(rInnerMax, m.center().Rho() + 0.5*m.dsDistance()) : rOuterMax = MAX(rOuterMax, m.center().Rho() + 0.5*m.dsDistance()));
        if (iMod>=2) break; // Sufficient to have a look at +-smallDelta (first 2 modules)
      }

      flatThetaEnd((m_bigParity() > 0 ? flatRod1.thetaEnd() : flatRod2.thetaEnd()));
      auto lastMod1     = zPlusModules1.back();
      auto lastMod2     = zPlusModules2.back();
      flatREndInner((m_bigParity() > 0 ? lastMod2.center().Rho() + 0.5* lastMod2.dsDistance() : lastMod1.center().Rho() + 0.5* lastMod1.dsDistance()));
      flatREndOuter((m_bigParity() > 0 ? lastMod1.center().Rho() + 0.5* lastMod1.dsDistance() : lastMod2.center().Rho() + 0.5* lastMod2.dsDistance()));
      flatEndModCentreZ((m_bigParity() > 0 ? lastMod2.center().Z() : lastMod1.center().Z()));
      flatZEnd((m_bigParity() > 0 ? lastMod1.flatMaxZ() : lastMod2.flatMaxZ())); //zEnd_REAL = (m_bigParity() > 0 ? lastMod1.planarMaxZ() : lastMod2.planarMaxZ());  //TAKE CAREEEEEE : REMOVE FLAT PART OVERLAP ?

      // Calculate phiOverlap of the module -> assumed that all modules have the same properties
      double width      = zPlusModules1.size()>0 ? zPlusModules1[0].meanWidth() : 0.0;
      double dsDistance = zPlusModules1.size()>0 ? zPlusModules1[0].dsDistance() : 0.0;
      double T          = tan(2.*M_PI / numberRods());
      double A          = 1. / (2. * rInnerMin);
      double B          = 1. / (2. * rOuterMin);
      double a          = T * A * B;
      double b          = - (A + B);
      double c          = - T;
      double s          = (-b - sqrt(b*b - 4*a*c))/(2*a);
      flatPhiOverlapSmallDeltaMinus(width + s);

      A = 1. / (2. * rInnerMax);
      B = 1. / (2. * rOuterMax);
      a = T * A * B;
      b = - (A + B);
      c = - T;
      s = (-b - sqrt(b*b - 4*a*c))/(2*a);
      flatPhiOverlapSmallDeltaPlus(width + s);

      // Calculate zInner & zOuter
      const auto& minusBigDeltaRod     = dynamic_cast<const RodPairStraight&>(m_bigParity() > 0 ? m_flatRods[1] : m_flatRods[0]);
      const auto& minusBigDeltaModules = minusBigDeltaRod.modules().first;
      const auto& plusBigDeltaRod      = dynamic_cast<const RodPairStraight&>(m_bigParity() > 0 ? m_flatRods[0] : m_flatRods[1]);
      const auto& plusBigDeltaModules  = plusBigDeltaRod.modules().first;

      int i = 0;
      double rStartInner;
      double zStartInner_REAL;
      double rEndInner;
      double zEndInner_REAL;
      int smallParity = minusBigDeltaRod.smallParity();
      for (const auto& m : minusBigDeltaModules) {
        if (i > 0) {
          zStartInner_REAL = m.planarMinZ();
          // Special case where ring has been built going upwards, and with zEndInner_REAL < zError
          if ((zStartInner_REAL < zEndInner_REAL) && (smallParity > 0)) {

            rEndInner -= m.dsDistance();
            rStartInner = m.center().Rho() + 0.5 * m.dsDistance();
          }
          // Standard case
          else {

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
              m_zErrorInnerFlat[i] = fact * (zStartInner_REAL - rStartInner / tan(zErrorInnerAngle));
            }
            else { // Case where (H1pUP, H1ppDOWN) is orthogonal to (Z).
              m_zErrorInnerFlat[i] = fact * zStartInner_REAL;
            }
          }
          else { // Case where consecutive modules are placed at the same r within a rod, for the Pixel for example.
            if (zStartInner_REAL != zEndInner_REAL) { // If consecutive modules (along a rod) are not touching, zError is undefined.
              m_zErrorInnerFlat[i] = std::numeric_limits<double>::quiet_NaN();
            }
            else { // If consecutive modules (along a rod) are touching, zError is infinity.
              m_zErrorInnerFlat[i] = std::numeric_limits<double>::infinity();
            }
          }
        } // i>0

        //zEndInner_REAL = m.planarMaxZ();
        zEndInner_REAL = m.flatMaxZ();
        rEndInner = m.center().Rho() + 0.5 * m.dsDistance();

        i++;
        smallParity = -smallParity;
      }

      i = 0;
      double rStartOuter;
      double zStartOuter_REAL;
      double rEndOuter;
      double zEndOuter_REAL;
      smallParity = plusBigDeltaRod.smallParity();
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
              m_zErrorOuterFlat[i] = fact * (zStartOuter_REAL - rStartOuter / tan(zErrorOuterAngle));
            }
            else {
              m_zErrorOuterFlat[i] = fact * zStartOuter_REAL;
            }
          }
          else {

            if (zStartOuter_REAL != zEndOuter_REAL) {
              m_zErrorOuterFlat[i] = std::numeric_limits<double>::quiet_NaN();
            }
            else {
              m_zErrorOuterFlat[i] = std::numeric_limits<double>::infinity();
            }
          }
        }

        //zEndOuter_REAL = m.planarMaxZ();
        zEndOuter_REAL = m.flatMaxZ();
        rEndOuter = m.center().Rho() + 0.5 * m.dsDistance();

        i++;
        smallParity = -smallParity;
      }
    }
  } // Flat part built?
  else logERROR("Parameters of straight part can't be calculated, because the straight part was not built!!!");
}

//
// GeometryVisitor pattern -> layer visitable
//
void Layer::accept(GeometryVisitor& v)
{
  v.visit(*this);
  for (auto& r : m_flatRods)   { r.accept(v); }
  for (auto& r : m_tiltedRods) { r.accept(v); }
}

//
// GeometryVisitor pattern -> layer visitable (const. option)
//
void Layer::accept(ConstGeometryVisitor& v) const
{
  v.visit(*this);
  for (const auto& r : m_flatRods)   { r.accept(v); }
  for (const auto& r : m_tiltedRods) { r.accept(v); }
}

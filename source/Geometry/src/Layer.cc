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
 minRAllMat     (string("minRAllMat")        ),
 maxRAllMat     (string("maxRAllMat")        ),
 buildNumModules(       "numModules"         , parsedOnly()),
 outerZ         (       "outerZ"             , parsedOnly()),
 radiusMode     (       "radiusMode"         , parsedAndChecked(), RadiusMode::AUTO),
 requestedAvgRadius(    "radius"             , parsedOnly()),
 avgBuildRadius (       "avgBuildRadius"     , parsedOnly()),
 sameParityRods (       "sameParityRods"     , parsedAndChecked(), true),
 layerRotation  (       "layerRotation"      , parsedOnly()      , 0.),
 buildNumModulesFlat(   "numModulesFlat"     , parsedOnly()),
 buildNumModulesTilted( "numModulesTilted"   , parsedOnly()),
 isTilted       (       "isTilted"           , parsedOnly(), false),
 isTiltedAuto   (       "isTiltedAuto"       , parsedOnly()),
 tiltedLayerSpecFile(   "tiltedLayerSpecFile", parsedOnly()),
 smallDelta     (       "smallDelta"         , parsedAndChecked()),
 m_smallParity  (       "smallParity"        , parsedAndChecked(),-1),
 bigDelta       (       "bigDelta"           , parsedAndChecked()),
 m_bigParity    (       "bigParity"          , parsedOnly()      ,-1),
 phiOverlap     (       "phiOverlap"         , parsedAndChecked(), 1.),
 phiSegments    (       "phiSegments"        , parsedAndChecked(), 4),
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

  // Set radius mode if not default value set or not defined from barrel level
  if ( (myid()==1)               && barrelMinRFixed) this->radiusMode(RadiusMode::FIXED);
  if ( (myid()==barrelNumLayers) && barrelMaxRFixed) this->radiusMode(RadiusMode::FIXED);
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

  // Non-tilted geometry
  if (!isTilted()) {

    if (buildNumModules() > 0  && outerZ.state()) throw PathfulException("Only one between numModules and maxZ can be specified");
    if (buildNumModules() == 0 &&!outerZ.state()) throw PathfulException("At least one between numModules and maxZ must be specified");
    if (!phiOverlap.state())                      throw PathfulException("Flat layer: phiOverlap must be specified.");
    if (!phiSegments.state())                     throw PathfulException("Flat layer: phiSegments must be specified.");
    //if (numRods.state())                          throw PathfulException("Flat layer : numRods should not be specified.");

    if (isTiltedAuto.state())          logERROR("Layer " + std::to_string(myid()) + ": doesn't make sense to specify isTiltedAuto. Not being used.");
    if (buildNumModulesFlat.state())   logERROR("Layer " + std::to_string(myid()) + ": doesn't make sense to specify numModulesFlat. Not being used.");
    if (buildNumModulesTilted.state()) logERROR("Layer " + std::to_string(myid()) + ": doesn't make sense to specify numModulesTilted. Not being used.");
  }
  // Tilted geometry
  else {

    if (!isTiltedAuto.state()) throw PathfulException("Tilted layer: isTiltedAuto must be specified.");
    if (phiOverlap.state())    throw PathfulException("Tilted layer: phiOverlap should not be specified.");
    if (phiSegments.state())   throw PathfulException("Tilted layer: phiSegments should not be specified.");

    // Automatic tilt algorithm applied
    if (isTiltedAuto()) {
      if (!buildNumModulesFlat.state())   throw PathfulException("Tilted layer with automatic placement: numModulesFlat must be specified.");
      if (!buildNumModulesTilted.state()) throw PathfulException("Tilted layer with automatic placement: numModulesTilted must be specified.");
      if (buildNumModules() > 0 && buildNumModulesFlat.state() && buildNumModulesTilted.state()) {
        if (buildNumModules() !=  (buildNumModulesFlat() + buildNumModulesTilted())) {
          throw PathfulException("Tilted layer : numModules != (numModulesFlat + numModulesTilted). Specify numModulesFlat and numModulesTilted only, don't specify numModules.");
        }
      }
      //if (!numRods.state()) throw PathfulException("Tilted layer with automatic placement : numRods must be specified.");
    }
    // Config file for tilting used instead
    else {
      // TODO: Not implemented
      throw PathfulException("Tilted geometry built on configuration file not implemented. Use automatic build instead!");

      if (buildNumModulesFlat.state())   logERROR("Tilted layer " + std::to_string(myid()) + ": doesn't make sense to specify numModulesFlat. Not being used.");
      if (buildNumModulesTilted.state()) logERROR("Tilted layer " + std::to_string(myid()) + ": doesn't make sense to specify numModulesTilted. Not being used.");
    }

    if (outerZ.state()) logERROR("Tilted layer: outerZ was specified. Routing of services will be forced to be at Z = outerZ.");
  }

  if (bigDelta()  <0)          throw PathfulException("Big delta parameter must be positive!");
  if (smallDelta()<0)          throw PathfulException("Small delta parameter must be positive!");
  if (bigDelta()<smallDelta()) throw PathfulException("Big delta parameter is expected to be bigger in size than small delta parameter!");
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
    if (!isTilted()) {

      buildStraight(barrelNumLayers, barrelMinR, barrelMaxR);
      buildNumModulesFlat(buildNumModules());
      buildNumModulesTilted(0);
    }
    else buildTilted(barrelNumLayers, barrelMinR, barrelMaxR);

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

  float optimalRadius;
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

  int   optimalNumRods = optimalModsPerSegment*phiSegments();
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

    minRadius = updatedMinR - bigDelta();
    maxRadius = updatedMaxR + bigDelta();
  }
  // Not same rods -> extreme given by optimal radius +-bigDelta
  else {
    minRadius = optimalRadius - bigDelta();
    maxRadius = optimalRadius + bigDelta();
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

      if (buildNumModules() > 0) rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, buildNumModules(), propertyTree());
      else                       rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, outerZ(), propertyTree());
      rod->build();
      firstRod = rod;

      m_rods.push_back(rod);
    }
    // Prototype of even rods
    else if (i==2) {

      // If same rods or sameParityRods required, each even/odd rod are the same -> useful from engineering point of view
      if (m_sameRods || sameParityRods()) {

        RodPair* rod = GeometryFactory::clone<RodPairStraight>(*firstRod);
        rod->buildClone(2, 2*bigDelta()*bigParity, rodPhiRotation);

        m_rods.push_back(rod);
      }
      else {

        RodPairStraight* rod = nullptr;

        if (buildNumModules() > 0) rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, buildNumModules(), propertyTree());
        else                       rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, outerZ(), propertyTree());
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
          rod->buildClone(i, 2*bigDelta()*bigParity, rodPhiRotation*(i-1));
        }

        m_rods.push_back(rod);
      }
    }
  } // For NumRods
}

//
// If tilted layer required, build() method internally calls buildTilted()
//
void Layer::buildTilted(int barrelNumLayers, double barrelMinR, double barrelMaxR) {

//  // Tilted geometry built as described in config file
//  if (!isTiltedAuto()) {
//    //TODO: NOT IMPLEMENTED NOW
//  }
//  // Tilted part built based on automatic algorithm
//  else {
//
//    double flatEndTheta  = M_PI / 2.;
//    double flatInnerEndR = 0;
//    double flatOuterEndR = 0;
//    double flatEndCentreZ= 0;
//    double flatEndMaxZ   = 0;
//
//    double flatInnerMinR = std::numeric_limits<double>::max();
//    double flatOuterMinR = std::numeric_limits<double>::max();
//    double flatInnerMaxR = 0;
//    double flatOuterMaxR = 0;
//    double flatAvgR      = 0;
//
//    double flatPhiOverlapMin = 0;
//    double flatPhiOverlapMax = 0;
//
//    std::vector<double> outerModsR; // List of radii of modules positioned at R + bigDelta
//    std::vector<double> innerModsR; // List of radii of modules positioned at R - bigDelta
//
//    // Build first the flat part of the layer
//    if (buildNumModulesFlat()!=0) {
//
//      // Update first the total number of modules & use this value in buildStraight method to built the flat part
//      buildNumModules(buildNumModulesFlat());
//      buildStraight(barrelNumLayers, barrelMinR, barrelMaxR);
//
//      // Calculate flat part properties
//      if (numRods()>=2) {
//
//        // Calculate min/max R of the flat part & of its last module
//        for (const auto& m : m_rods[0].m_zPlusModules) {
//          (m_bigParity()>0 ? flatOuterMinR = MIN(flatOuterMinR,m.center().Rho()+0.5*m.dsDistance()) : flatInnerMinR = MIN(flatInnerMinR,m.center().Rho()+0.5*m.dsDistance()));
//          (m_bigParity()>0 ? flatOuterMaxR = MAX(flatOuterMaxR,m.center().Rho()+0.5*m.dsDistance()) : flatInnerMaxR = MAX(flatInnerMaxR,m.center().Rho()+0.5*m.dsDistance()));
//        }
//
//        for (const auto& m : m_rods[1].m_zPlusModules) {
//          (m_bigParity()>0 ? flatInnerMinR = MIN(flatInnerMinR,m.center().Rho()+0.5*m.dsDistance()) : flatOuterMinR = MIN(flatOuterMinR,m.center().Rho()+0.5*m.dsDistance()));
//          (m_bigParity()>0 ? flatInnerMaxR = MAX(flatInnerMaxR,m.center().Rho()+0.5*m.dsDistance()) : flatOuterMaxR = MAX(flatOuterMaxR,m.center().Rho()+0.5*m.dsDistance()));
//        }
//
//        // Calculate theta of the last module in the flat part
//        if (m_bigParity()>0) {
//          auto lastMod = m_rods[0].m_zPlusModules.back()
//          flatEndTheta = atan( (lastMod.center().Rho()+0.5*lastMod.dsDistance())/lastMod.planarMaxZ() );
//        }
//        else {
//          auto lastMod = m_rods[1].m_zPlusModules.back()
//          flatEndTheta = atan( (lastMod.center().Rho()+0.5*lastMod.dsDistance())/lastMod.planarMaxZ() );
//        }
//
//        // Calculate R & Z of the last module in the flat part
//        auto lastMod0 = m_rods[0].m_zPlusModules.back();
//        auto lastMod1 = m_rods[1].m_zPlusModules.back();
//        flatInnerEndR = (m_bigParity()>0 ? lastMod1.center().Rho()+0.5*lastMod1.dsDistance() : lastMod0.center().Rho()+0.5* lastMod0.dsDistance());
//        flatOuterEndR = (m_bigParity()>0 ? lastMod0.center().Rho()+0.5*lastMod0.dsDistance() : lastMod1.center().Rho()+0.5* lastMod1.dsDistance());
//        flatEndCentreZ= (m_bigParity()>0 ? lastMod1.center().Z() : lastMod0.center().Z());
//        flatEndMaxZ   = (m_bigParity()>0 ? lastMod0.planarMaxZ() : lastMod1.planarMaxZ());  // TODO: Remove module overlap?
//
//        // Calculate overlaps of flat part in phi
//        double modWidth   = lastMod0.meanWidth(); // All modules in the layer assumed to be the same
//        double dsDistance = lastMod0.dsDistance();
//
//        double T = tan(2.*M_PI/numRods());
//        double A = 0.5/flatInnerMinR;
//        double B = 0.5/flatOuterMinR;
//        double a = T * A * B;
//        double b = - (A + B);
//        double c = - T;
//        double s = (-b - sqrt(b*b - 4*a*c))/(2*a);
//        flatPhiOverlapMin = modWidth + s;
//
//        A = 0.5/flatInnerMaxR;
//        B = 0.5/flatOuterMaxR;
//        a = T * A * B;
//        b = - (A + B);
//        c = - T;
//        s = (-b - sqrt(b*b - 4*a*c))/(2*a);
//        flatPhiOverlapMax = modWidth + s;
//
//        flatAvgR = (flatInnerMinR + flatInnerMaxR + flatOuterMinR + flatOuterMaxR)/4. - 0.5*dsDistance;
//
//  flatRingsGeometryInfo_.calculateFlatRingsGeometryInfo(flatPartRods_, bigParity());
//      }
//      else { logERROR(to_string(flatPartRods_.size()) + " straight rod was built for the whole flat part."); }
//    }
//
//
//    if (tmspecsi.size() != buildNumModulesFlat()) {
//      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " but flat part inner rod has " + to_string(tmspecsi.size()) + " module(s).");
//    }
//    if (tmspecso.size() != buildNumModulesFlat()) {
//      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " but flat part outer rod has " + to_string(tmspecso.size()) + " module(s).");
//    }
//
//
//
//
//    tiltedRingsGeometry_ = makeTiltedRingsTemplate(flatPartThetaEnd);
//
//    if (tiltedRingsGeometry_.size() == buildNumModulesTilted()) {
//      for (int i = 0; i < buildNumModulesTilted(); i++) {
//  int ringNumber = buildNumModulesFlat() + 1 + i;
//  TiltedModuleSpecs ti{ tiltedRingsGeometry_[ringNumber]->innerRadius(), tiltedRingsGeometry_[ringNumber]->zInner(), tiltedRingsGeometry_[ringNumber]->tiltAngle()*M_PI/180. };
//  TiltedModuleSpecs to{ tiltedRingsGeometry_[ringNumber]->outerRadius(), tiltedRingsGeometry_[ringNumber]->zOuter(), tiltedRingsGeometry_[ringNumber]->tiltAngle()*M_PI/180. };
//
//  if (ti.valid()) tmspecsi.push_back(ti);
//  if (to.valid()) tmspecso.push_back(to);
//      }
//      tiltedRingsGeometryInfo_ = TiltedRingsGeometryInfo(buildNumModulesFlat(), flatPartrEndInner, flatPartrEndOuter, flatPartzEnd,  flatPartzEnd_REAL, tiltedRingsGeometry_);
//    }
//    else {
//      logERROR("numModulesTilted = " + to_string(buildNumModulesTilted()) + " but tilted rings geometry template has " + to_string(tiltedRingsGeometry_.size()) + " elements.");
//    }
//
//
//    buildNumModules(buildNumModulesFlat() + buildNumModulesTilted());
//
//    if (tmspecsi.size() != buildNumModules()) {
//      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " and numModulesTilted = " + to_string(buildNumModulesTilted()) + " but tilted rod 1 has " + to_string(tmspecsi.size()) + " module(s) in total.");
//    }
//    if (tmspecso.size() != buildNumModules()) {
//      logERROR("Layer " + to_string(myid()) + " : numModulesFlat = " + to_string(buildNumModulesFlat()) + " and numModulesTilted = " + to_string(buildNumModulesTilted()) + " but tilted rod 2 has " + to_string(tmspecso.size()) + " module(s) in total.");
//    }
//
//  }
//
//  RodTemplate rodTemplate = makeRodTemplate();
//
//  float rodPhiRotation = 2*M_PI/numRods();
//
//  TiltedRodPair* first = GeometryFactory::make<TiltedRodPair>();
//  first->myid(1);
//  first->isOuterRadiusRod(false);
//  first->store(propertyTree());
//  first->build(rodTemplate, tmspecsi, 1);
//  rods_.push_back(first);
//
//  TiltedRodPair* second = GeometryFactory::make<TiltedRodPair>();
//  second->myid(2);
//  second->isOuterRadiusRod(true);
//  second->store(propertyTree());
//  second->build(rodTemplate, tmspecso, 0);
//  second->rotateZ(rodPhiRotation);
//  rods_.push_back(second);
//
//  for (int i = 2; i < numRods(); i++) {
//    RodPair* rod = i%2 ? GeometryFactory::clone(*second) : GeometryFactory::clone(*first); // clone rods
//    rod->myid(i+1);
//    rod->rotateZ(rodPhiRotation*(i%2 ? i-1 : i));
//    rods_.push_back(rod);
//    }
//
//  // computing the layer's place radius as the average of all the modules' radii
//  placeRadius_  = std::accumulate(tmspecsi.begin(), tmspecsi.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
//  placeRadius_ += std::accumulate(tmspecso.begin(), tmspecso.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
//  placeRadius_ /= tmspecsi.size() + tmspecso.size();

}

//
// Setup: link lambda functions to various layer related properties (use setup functions for ReadOnly Computable properties)
//
void Layer::setup()
{
  maxZ.setup([&]()       { if (m_rods.size()>0) return m_rods.front().maxZ(); else return 0.0; });
  minZ.setup([&]()       { if (m_rods.size()>0) return m_rods.front().minZ(); else return 0.0; });
  maxR.setup([&]()       { double max = 0;                                  for (const auto& r : m_rods) { max = MAX(max, r.maxR()); } return max; });
  minR.setup([&]()       { double min = std::numeric_limits<double>::max(); for (const auto& r : m_rods) { min = MIN(min, r.minR()); } return min; });
  maxRAllMat.setup([&]() { double max = 0;                                  for (const auto& r : m_rods) { max = MAX(max, r.maxRAllMat()); } return max; });
  minRAllMat.setup([&]() { double min = std::numeric_limits<double>::max(); for (const auto& r : m_rods) { min = MIN(min, r.minRAllMat()); } return min; });
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






#include "Layer.h"
#include "RodPair.h"
#include "MessageLogger.h"
#include "ConversionStation.h"

using std::string;

define_enum_strings(RadiusMode) = { "shrink", "enlarge", "fixed", "auto" };

FlatRingsGeometryInfo::FlatRingsGeometryInfo() {}

void FlatRingsGeometryInfo::calculateFlatRingsGeometryInfo(std::vector<RodPairStraight*> flatPartRods, int bigParity) {
 
  RodPairStraight* minusBigDeltaRod = (bigParity > 0 ? flatPartRods.at(1) : flatPartRods.front());
  const auto& minusBigDeltaModules = minusBigDeltaRod->modules().first;
  RodPairStraight* plusBigDeltaRod = (bigParity > 0 ? flatPartRods.front() : flatPartRods.at(1));
  const auto& plusBigDeltaModules = plusBigDeltaRod->modules().first;  

  int i = 0;
  double rStartInner;
  double zStartInner_REAL;
  double rEndInner;
  double zEndInner_REAL;
  int smallParity = minusBigDeltaRod->smallParity();
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
  smallParity = plusBigDeltaRod->smallParity();
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
 tiltedLayerSpecFile(   "tiltedLayerSpecFile", parsedOnly()),
 smallDelta     (       "smallDelta"         , parsedAndChecked()),
 m_smallParity  (       "smallParity"        , parsedAndChecked(),-1),
 bigDelta       (       "bigDelta"           , parsedAndChecked()),
 m_bigParity    (       "bigParity"          , parsedOnly()      ,-1),
 phiOverlap     (       "phiOverlap"         , parsedOnly(), 1.),
 phiSegments    (       "phiSegments"        , parsedOnly(), 4),
 m_useMinMaxRCorrect(   "useMinMaxRCorrect"  , parsedAndChecked(), true),
 numberRods        ("numberRods"        , parsedOnly()),
 m_ringNode     (       "Ring"               , parsedOnly()),
 m_stationsNode (       "Station"            , parsedOnly()),
 buildNumModulesFlat("numModulesFlat"     , parsedOnly()),
 buildNumModulesTilted("numModulesTilted"     , parsedOnly()),
 isTilted       ("isTilted"       , parsedOnly(), false),
 isTiltedAuto   ("isTiltedAuto"   , parsedOnly()),
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

  if (!isTilted()) {
    if (buildNumModules() > 0 && outerZ.state()) throw PathfulException("Only one between numModules and outerZ can be specified");
    if (buildNumModules() == 0 && !outerZ.state()) throw PathfulException("At least one between numModules and outerZ must be specified");
    if (!phiOverlap.state()) throw PathfulException("Flat layer : phiOverlap must be specified.");
    if (!phiSegments.state()) throw PathfulException("Flat layer : phiSegments must be specified.");
    if (numberRods.state()) throw PathfulException("Flat layer : numberRods should not be specified.");
    if (isTiltedAuto.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify isTiltedAuto. Not used.");
  }

  if (!isTilted() || (isTilted() && !isTiltedAuto())) {
    if (buildNumModulesFlat.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify numModulesFlat. Not used.");
    if (buildNumModulesTilted.state()) logERROR("Layer " + std::to_string(myid()) + " : doesn't make sense to specify numModulesTilted. Not used.");
  }

  if (isTilted()) {
    if (outerZ.state()) logERROR("Tilted layer : outerZ was specified. Routing of services will be forced to be at Z = outerZ.");
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
      if (!numberRods.state()) throw PathfulException("Tilted layer with automatic placement : numberRods must be specified.");
    }
  }

  if (bigDelta()  <0)            throw PathfulException("Big delta parameter must be positive!");
  if (smallDelta()<0)            throw PathfulException("Small delta parameter must be positive!");
  if (bigDelta()<smallDelta()) throw PathfulException("Big delta parameter is expected to be bigger in size than small delta parameter!");
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



RodTemplate Layer::makeRodTemplate() {
  RodTemplate rodTemplate(buildNumModules() > 0 ? buildNumModules() : (!m_ringNode.empty() ? m_ringNode.rbegin()->first + 1 : 1)); // + 1 to make room for a default constructed module to use when building rods in case the rodTemplate vector doesn't have enough elements
  //std::cout << "rodTemplate.size() = " << rodTemplate.size() << std::endl;
  for (int i = 0; i < rodTemplate.size(); i++) {
    int ringNumber = i+1;
    rodTemplate[i] = std::move(unique_ptr<BarrelModule>(GeometryFactory::make<BarrelModule>(ringNumber, GeometryFactory::make<RectangularModule>(), m_ringNode, propertyTree() )));
    //rodTemplate[i]->store(propertyTree());
    //if (m_ringNode.count(i+1) > 0) rodTemplate[i]->store(m_ringNode.at(i+1));
    rodTemplate[i]->build();
  }
  return rodTemplate;
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
// If straight layer required, build() method internally calls buildStraight()
//
void Layer::buildStraight(int barrelNumLayers, double barrelMinR, double barrelMaxR) {

  if (!isTilted()) {
    if (buildNumModules() > 0 ) buildNumModulesFlat(buildNumModules());
    buildNumModulesTilted(0);
  }
  else { 
    buildNumModules(buildNumModulesFlat());
  }

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
    int    smallParity;
    if (!isTilted()) smallParity = (i%2) ? (m_smallParity()) : (sameParityRods() ? m_smallParity() : -m_smallParity());
    else smallParity = pow(-1, buildNumModulesFlat());  // if tilted, smallParity is forced so that central module is at +smallDelta

    // Prototypes of odd rods
    if (i==1) {

      RodPairStraight* rod = nullptr;

      if (buildNumModules() > 0) rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, buildNumModules(), propertyTree());
      else rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, outerZ(), propertyTree());

      rod->build();
      
      firstRod = rod;

      
      if (!isTilted()) { m_rods.push_back(rod); buildNumModulesFlat(rod->numModulesSide(1)); }
      else { m_flat_rods.push_back(rod); }
    }

    // Prototype of even rods
    else if (i==2) {

      // If same rods or sameParityRods required, each even/odd rod are the same -> useful from engineering point of view
      if (m_sameRods || sameParityRods()) {

        RodPairStraight* rod = GeometryFactory::clone<RodPairStraight>(*firstRod);
        rod->buildClone(2, 2*bigDelta()*bigParity, rodPhiRotation);

        if (!isTilted()) { m_rods.push_back(rod); }
	else { m_flat_rods.push_back(rod); }
      }
      else {

        RodPairStraight* rod = nullptr;

        if (buildNumModules() > 0) rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, buildNumModules(), propertyTree());
        else                       rod = GeometryFactory::make<RodPairStraight>(i, minRadius, maxRadius, optimalRadius, rotation, bigDelta(), bigParity, smallDelta(), smallParity, outerZ(), propertyTree());
        rod->build();
        secondRod = rod;

        if (!isTilted()) { m_rods.push_back(rod); }
	else { m_flat_rods.push_back(rod); }
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

	if (!isTilted()) { m_rods.push_back(rod); }
	else { m_flat_rods.push_back(rod); }

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

        if (!isTilted()) { m_rods.push_back(rod); }
	else { m_flat_rods.push_back(rod); }
      }
    }
  } // For NumRods
}



TiltedRingsTemplate Layer::makeTiltedRingsTemplate(double flatPartThetaEnd) {
  TiltedRingsTemplate tiltedRingsGeometry;

  for (int i = (buildNumModulesFlat() + 1); i < (buildNumModulesFlat() + buildNumModulesTilted() + 1); i++) {

    TiltedRing* tiltedRing = GeometryFactory::make<TiltedRing>();
    tiltedRing->myid(i);
    tiltedRing->store(propertyTree());
    if (m_ringNode.count(i) > 0) tiltedRing->store(m_ringNode.at(i));
    tiltedRing->numPhi(numberRods());

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
      numberRods(tokens[6]); // this assumes every row of the spec file has the same value for the last column (num rods in phi) 
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
    // Warning ! it seems that bigParity is inverted with respect to the tkLayout/dev version

    if (buildNumModulesFlat() != 0) {

      if (m_flat_rods.size() >= 2) {
	RodPairStraight* flatPartRod1 = m_flat_rods.front();
	const auto& zPlusModules1 = flatPartRod1->modules().first;
	for (const auto& m : zPlusModules1) {
	  TiltedModuleSpecs t1{m.center().Rho(), m.center().Z(), 0.0};
	  if (t1.valid()) (m_bigParity() > 0 ? tmspecso.push_back(t1) : tmspecsi.push_back(t1));
	  (m_bigParity() > 0 ? flatPartrOuterSmall = MIN(flatPartrOuterSmall, m.center().Rho() + 0.5*m.dsDistance()) : flatPartrInnerSmall = MIN(flatPartrInnerSmall, m.center().Rho() + 0.5*m.dsDistance()));
	  (m_bigParity() > 0 ? flatPartrOuterBig = MAX(flatPartrOuterBig, m.center().Rho() + 0.5*m.dsDistance()) : flatPartrInnerBig = MAX(flatPartrInnerBig, m.center().Rho() + 0.5*m.dsDistance()));
	}

	RodPairStraight* flatPartRod2 = m_flat_rods.at(1);
	const auto& zPlusModules2 = flatPartRod2->modules().first;
	for (const auto& m : zPlusModules2) {
	  TiltedModuleSpecs t2{m.center().Rho(), m.center().Z(), 0.0};
	  if (t2.valid()) (m_bigParity() > 0 ? tmspecsi.push_back(t2) : tmspecso.push_back(t2));
	  (m_bigParity() > 0 ? flatPartrInnerSmall = MIN(flatPartrInnerSmall, m.center().Rho() + 0.5*m.dsDistance()) : flatPartrOuterSmall = MIN(flatPartrOuterSmall, m.center().Rho() + 0.5*m.dsDistance()));
	  (m_bigParity() > 0 ? flatPartrInnerBig = MAX(flatPartrInnerBig, m.center().Rho() + 0.5*m.dsDistance()) : flatPartrOuterBig = MAX(flatPartrOuterBig, m.center().Rho() + 0.5*m.dsDistance()));
	}

	flatPartThetaEnd = (m_bigParity() > 0 ? flatPartRod1->thetaEnd_REAL() : flatPartRod2->thetaEnd_REAL());
	auto lastMod1 = zPlusModules1.back();
	auto lastMod2 = zPlusModules2.back();	
	flatPartrEndInner = (m_bigParity() > 0 ? lastMod2.center().Rho() + 0.5* lastMod2.dsDistance() : lastMod1.center().Rho() + 0.5* lastMod1.dsDistance());
	flatPartrEndOuter = (m_bigParity() > 0 ? lastMod1.center().Rho() + 0.5* lastMod1.dsDistance() : lastMod2.center().Rho() + 0.5* lastMod2.dsDistance());
	flatPartzEnd = (m_bigParity() > 0 ? lastMod2.center().Z() : lastMod1.center().Z());	
	flatPartzEnd_REAL = (m_bigParity() > 0 ? lastMod1.planarMaxZ() : lastMod2.planarMaxZ());	//TAKE CAREEEEEE : REMOVE FLAT PART OVERLAP ?


	RectangularModule* flatPartrmod = GeometryFactory::make<RectangularModule>();
	flatPartrmod->store(propertyTree());
	if (m_ringNode.count(1) > 0) flatPartrmod->store(m_ringNode.at(1));
	flatPartrmod->build();
	double width = flatPartrmod->width();
	double dsDistance = flatPartrmod->dsDistance();
	double T = tan(2.*M_PI / numberRods());
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

	flatRingsGeometryInfo_.calculateFlatRingsGeometryInfo(m_flat_rods, m_bigParity());


      }
      else { logERROR(to_string(m_flat_rods.size()) + " straight rod was built for the whole flat part."); }     
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

  float rodPhiRotation = 2*M_PI/numberRods();

  //TiltedRodPair* first = GeometryFactory::make<TiltedRodPair>();
  TiltedRodPair* first = GeometryFactory::make<TiltedRodPair>(1, 0., propertyTree());
  //first->myid(1);
  first->isOuterRadiusRod(false);
  //first->store(propertyTree());
  first->build(rodTemplate, tmspecsi, 1);
  m_rods.push_back(first);

  TiltedRodPair* second = GeometryFactory::make<TiltedRodPair>(2, rodPhiRotation, propertyTree());
  //second->myid(2);
  second->isOuterRadiusRod(true);
  //second->store(propertyTree());
  second->build(rodTemplate, tmspecso, 0);
  //second->rotateZ(rodPhiRotation);
  m_rods.push_back(second);

  for (int i = 2; i < numberRods(); i++) {
    RodPair* rod = i%2 ? GeometryFactory::clone(*second) : GeometryFactory::clone(*first); // clone rods
    rod->myid(i+1);
    rod->rotateZ(rodPhiRotation*(i%2 ? i-1 : i));
    m_rods.push_back(rod);
    }

  // computing the layer's place radius as the average of all the modules' radii
  double placeRadius;
  placeRadius  = std::accumulate(tmspecsi.begin(), tmspecsi.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
  placeRadius += std::accumulate(tmspecso.begin(), tmspecso.end(), 0., [](double x, const TiltedModuleSpecs& t) { return x+t.r; });
  placeRadius /= tmspecsi.size() + tmspecso.size();
  avgBuildRadius(placeRadius);

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

    // Build either flat layer, either flat part of tilted layer
    buildStraight(barrelNumLayers, barrelMinR, barrelMaxR);

    // Build tilted part of tilted layer    
    if (isTilted()) buildTilted();


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

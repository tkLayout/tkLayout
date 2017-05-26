#include "RodPair.h"

#include "MessageLogger.h"
#include "SimParms.h"

define_enum_strings(StartZMode) = { "modulecenter", "moduleedge" };

//
//  Constructor - parse geometry config file using boost property tree & read-in Rod parameters -> use number of modules to build rod
//
RodPair::RodPair(int id, double minRadius, double maxRadius, double radius, double rotation, int numModules, const PropertyTree& treeProperty) :
 minZ              (string("minZ")              ),
 maxZ              (string("maxZ")              ),
 minR              (string("minR")              ),
 maxR              (string("maxR")              ),
 minRAllMat        (string("minRAllMat")        ),
 maxRAllMat        (string("maxRAllMat")        ),
 maxModuleThickness(string("maxModuleThickness")),
 m_materialObject(MaterialObject::ROD),
 m_startZMode(             "startZMode"         , parsedOnly(), StartZMode::MODULECENTER),
 m_nModules     (numModules),
 m_outerZ       (0.0),
 m_optimalRadius(radius),
 m_minRadius(    minRadius),
 m_maxRadius(    maxRadius),
 m_rotation(     rotation)
{
  // Set the geometry config parameters
  this->myid(id);
  this->store(treeProperty);
}

//
//  Constructor - parse geometry config file using boost property tree & read-in Rod parameters -> use outerZ to build rod
//
RodPair::RodPair(int id, double minRadius, double maxRadius, double radius, double rotation, double outerZ, const PropertyTree& treeProperty) :
 minZ              (string("minZ")              ),
 maxZ              (string("maxZ")              ),
 minR              (string("minR")              ),
 maxR              (string("maxR")              ),
 minRAllMat        (string("minRAllMat")        ),
 maxRAllMat        (string("maxRAllMat")        ),
 maxModuleThickness(string("maxModuleThickness")),
 m_materialObject(MaterialObject::ROD),
 m_startZMode(             "startZMode"         , parsedOnly(), StartZMode::MODULECENTER),
 m_nModules     (0),
 m_outerZ       (outerZ),
 m_optimalRadius(radius),
 m_minRadius(    minRadius),
 m_maxRadius(    maxRadius),
 m_rotation(     rotation)
{
 // Set the geometry config parameters
 this->myid(id);
 this->store(treeProperty);
}

//
//  Constructor - parse geometry config file using boost property tree & read-in Rod parameters
//
RodPair::RodPair(int id, double rotation, const PropertyTree& treeProperty) :
 minZ              (string("minZ")              ),
 maxZ              (string("maxZ")              ),
 minR              (string("minR")              ),
 maxR              (string("maxR")              ),
 minRAllMat        (string("minRAllMat")        ),
 maxRAllMat        (string("maxRAllMat")        ),
 maxModuleThickness(string("maxModuleThickness")),
 m_materialObject(MaterialObject::ROD),
 m_startZMode(             "startZMode"         , parsedOnly(), StartZMode::MODULECENTER),
 m_nModules     (0),
 m_outerZ       (0),
 m_optimalRadius(0),
 m_minRadius(    0),
 m_maxRadius(    0),
 m_rotation(     rotation)
{
 // Set the geometry config parameters
 this->myid(id);
 this->store(treeProperty);
}

//
// Position newly individual modules if RodPair cloned from a RodPair (i.e. rotate by respective angle and shift in R) and set rod new id
//
void RodPair::buildClone(int id, double shiftR, double rotation)
{
  // Set new id
  this->myid(id);

  // Rotate all modules -> needs to follow after translation!!!
  for (auto& m : m_zPlusModules)  { m.translateR(shiftR); m.rotateZ(rotation); }
  for (auto& m : m_zMinusModules) { m.translateR(shiftR); m.rotateZ(rotation); }
}

//
// Limit rod geometry by eta cut
//
void RodPair::cutAtEta(double eta)
{
  m_zPlusModules.erase_if([eta](BarrelModule& m)  { return fabs(m.center().Eta()) > eta; });
  m_zMinusModules.erase_if([eta](BarrelModule& m) { return fabs(m.center().Eta()) > eta; });
}

void RodPair::rotateZ(double angle) {
  for (auto& m : m_zPlusModules) { m.rotateZ(angle); }
  for (auto& m : m_zMinusModules) { m.rotateZ(angle); }
  //clearComputables();
}

//! GeometryVisitor pattern -> rod visitable
void RodPair::accept(GeometryVisitor& v)
{
  v.visit(*this);
  for (auto& m : m_zPlusModules) { m.accept(v); }
  for (auto& m : m_zMinusModules) { m.accept(v); }
}

//! GeometryVisitor pattern -> rod visitable (const. option)
void RodPair::accept(ConstGeometryVisitor& v) const
{
  v.visit(*this);
  for (const auto& m : m_zPlusModules) { m.accept(v); }
  for (const auto& m : m_zMinusModules) { m.accept(v); }
}

//! Setup: link lambda functions to various rod related properties (use setup functions for ReadOnly Computable properties)
void RodPair::setup() {

  minZ.setup([&]()       { if (m_zMinusModules.size()==0) return std::numeric_limits<double>::max(); else return minget2(m_zMinusModules.begin(), m_zMinusModules.end(), &DetectorModule::minZ); }); // One needs the minZ so don't bother with scanning the zPlus vector
  maxZ.setup([&]()       { if (m_zPlusModules.size()==0)  return std::numeric_limits<double>::min(); else return maxget2(m_zPlusModules.begin() , m_zPlusModules.end() , &DetectorModule::maxZ); });
  minR.setup([&]()       { if (m_zPlusModules.size()==0)  return std::numeric_limits<double>::max(); else return minget2(m_zPlusModules.begin() , m_zPlusModules.end() , &DetectorModule::minR); }); // MinR and maxR can be found just by scanning the zPlus vector, since the rod pair is symmetrical in R
  maxR.setup([&]()       { if (m_zPlusModules.size()==0)  return 0.0                               ; else return maxget2(m_zPlusModules.begin() , m_zPlusModules.end() , &DetectorModule::maxR);  });
  minRAllMat.setup([&]() { if (m_zPlusModules.size()==0)  return std::numeric_limits<double>::max(); else return minget2(m_zPlusModules.begin() , m_zPlusModules.end() , &DetectorModule::minRAllMat); });
  maxRAllMat.setup([&]() { if (m_zPlusModules.size()==0)  return 0.0                               ; else return maxget2(m_zPlusModules.begin() , m_zPlusModules.end() , &DetectorModule::maxRAllMat); });
  maxModuleThickness.setup([&]() { if (m_zPlusModules.size()==0) return 0.0                        ; else return maxget2(m_zPlusModules.begin() , m_zPlusModules.end() , &DetectorModule::thickness); });
}

//
//  Constructor - parse geometry config file using boost property tree & read-in Rod parameters -> use number of modules to build rod
//
RodPairStraight::RodPairStraight(int id, double minRadius, double maxRadius, double radius, double rotation, double bigDelta, int bigParity, double smallDelta, int smallParity, int numModules,
                 const PropertyTree& treeProperty) :
 RodPair(id, minRadius, maxRadius, radius, rotation, numModules, treeProperty),
 forbiddenRange(      "forbiddenRange"      , parsedOnly()),
 zOverlap(            "zOverlap"            , parsedAndChecked() , 1.),
 zError(              "zError"              , parsedAndChecked() , SimParms::getInstance().useLumiRegInGeomBuild() ? SimParms::getInstance().zErrorIP() : 0.0),
 useCompression(      "useCompressionInZ"   , parsedOnly(), true),
 useBalancing(        "useBalancingInZ"     , parsedOnly(), false),
 allowCompressionCuts("allowCompressionCuts", parsedOnly(), true),
 m_ringNode(          "Ring"                , parsedOnly()),
 m_smallDelta(smallDelta),
 m_smallParity(smallParity),
 m_bigDelta(bigDelta),
 m_bigParity(bigParity)
{
  // Set the geometry config parameters (id set through base class)
  this->store(treeProperty);
}

//
//  Constructor - parse geometry config file using boost property tree & read-in Rod parameters -> use outerZ to build rod
//
RodPairStraight::RodPairStraight(int id, double minRadius, double maxRadius, double radius, double rotation, double bigDelta, int bigParity, double smallDelta, int smallParity, double outerZ ,
                 const PropertyTree& treeProperty) :
 RodPair(id, minRadius, maxRadius, radius, rotation, outerZ, treeProperty),
 forbiddenRange(      "forbiddenRange"      , parsedOnly()),
 zOverlap(            "zOverlap"            , parsedAndChecked() , 1.),
 zError(              "zError"              , parsedAndChecked() , SimParms::getInstance().useLumiRegInGeomBuild() ? SimParms::getInstance().zErrorIP() : 0.0),
 useCompression(      "useCompressionInZ"   , parsedOnly(), true),
 useBalancing(        "useBalancingInZ"     , parsedOnly(), false),
 allowCompressionCuts("allowCompressionCuts", parsedOnly(), true),
 m_ringNode(          "Ring"                , parsedOnly()),
 m_smallDelta(smallDelta),
 m_smallParity(smallParity),
 m_bigDelta(bigDelta),
 m_bigParity(bigParity)
{
  // Set the geometry config parameters (id set through base class)
  this->store(treeProperty);
}

//
// Cross-check parameters provided from geometry configuration file
//
void RodPairStraight::check() {

  if (useBalancing() && m_nModules!=0) throw PathfulException("Balancing algorithm can't be used if number of modules specified, instead of outerZ!!!");;
}

//
//Number of iterations allowed when positioning or compressing modules
//
const int   RodPairStraight::c_nIterations = 100;

//
// Safety factor used when compressing modules in Z
//
const float RodPairStraight::c_safetySpaceFactor = 0.1;

//
// Build recursively individual modules
//
void RodPairStraight::build()
{
  try {

    m_materialObject.store(propertyTree());
    m_materialObject.build();

    logINFO(Form("Building %s prototype", fullid(*this).c_str()));
    check();

    // Get module length from a tree
    ReadonlyProperty<double, NoDefault> moduleLength( "length", parsedOnly());
    evaluateProperty(moduleLength);
    ReadonlyProperty<double, Default>   dsDistance("dsDistance", parsedAndChecked(), 0.);
    evaluateProperty(dsDistance);

    //
    // When positioning modules iterate n-times to have +-Z balanced rod, use safety margin to avoid too many iterations
    int    bIter          = 0;   // Balancing algorithm iterates
    double zUnbalance     = 0.0; // Total size of +-Z inbalance/2.
    double startZPos      = 0.0; // Update edge beteween positive & negative modules as balancing algorithm iterates

    while (bIter<c_nIterations) {

      // Unbalance in +-Z Ok, below required limit
      if (bIter!=0 && fabs(zUnbalance)<c_safetySpaceFactor) {
        break;
      }
      else {

        m_zPlusModules.release();
        m_zMinusModules.release();
      }

      // Position variables
      double lastZPos       = 0.0;          // Last ZPos at which module is positioned
      double newZPos        = 0.0;          // New ZPos at which module is positioned
      double lastRPos       = 0.0;          // Last radius using minRodR/maxRodR values -> used to calculate optimal module Z position only
      double newRPos        = 0.0;          // New radius using minRodR/maxRodR values -> used to calculate optimal module Z position only
      double newROptimal    = 0.0;          // Real radius at which module will be positioned in R
      double newDsDistance  = dsDistance(); // ds distance of newly positioned module
      double lastDsDistance = 0.0;          // ds distance of last positioned module
      int    iMod           = 1;            // Current number of built modules

      //
      // Position first module
      int parity = m_smallParity;

      BarrelModule* mod = GeometryFactory::make<BarrelModule>(iMod, GeometryFactory::make<RectangularModule>(), m_ringNode, propertyTree());
      mod->build();

      // Translate
      newRPos   = (parity > 0 ? m_maxRadius + m_smallDelta : m_minRadius - m_smallDelta) - newDsDistance/2; // Ds>0 lower sensor plays a role
      newZPos   = m_startZMode()==StartZMode::MODULECENTER ? 0.0 : startZPos + moduleLength()/2. - zUnbalance;
      //std::cout << "Module Z: " << newZPos << std::endl;

      newROptimal  = m_bigParity > 0 ? m_bigDelta : -m_bigDelta;
      newROptimal += (parity > 0 ? m_optimalRadius + m_smallDelta : m_optimalRadius - m_smallDelta);

      mod->translateR(newROptimal);
      mod->translateZ(newZPos);

      m_zPlusModules.push_back(mod);

      // Update info
      lastZPos       = newZPos;
      lastDsDistance = newDsDistance;
      parity        *= -1;
      iMod++;

      //
      // Position other modules in positive Z

      // To avoid end-less loop
      int    targetMods = m_nModules;
      double targetZ    = m_outerZ;
      int    iIter      = 0;
      while (iMod<=targetMods || (fabs(newZPos)+moduleLength()/2.)<targetZ) {

        iIter++;
        if (iIter>RodPairStraight::RodPairStraight::c_nIterations) logERROR("When positioning modules in positive Z, number of iterations exceeded allowed limit! Quitting!!!");

        BarrelModule* mod = GeometryFactory::make<BarrelModule>(iMod, GeometryFactory::make<RectangularModule>(), m_ringNode, propertyTree());
        mod->build();

        // Translate
        newRPos  = (parity > 0 ? m_maxRadius + m_smallDelta : m_minRadius - m_smallDelta) - newDsDistance/2; // Ds>0 lower sensor plays a role
        lastRPos = (parity > 0 ? m_maxRadius - m_smallDelta : m_minRadius + m_smallDelta) + newDsDistance/2; // Ds>0 upper sensor plays a role
        newZPos  = computeNextZ(lastRPos, newRPos, lastZPos, moduleLength(), lastDsDistance, mod->dsDistance(), BuildDir::RIGHT, parity);

        newROptimal  = m_bigParity > 0 ? m_bigDelta : -m_bigDelta;
        newROptimal += (parity > 0 ? m_optimalRadius + m_smallDelta : m_optimalRadius - m_smallDelta);

        //std::cout << "Module Z: " << newZPos << std::endl;
        mod->translateR(newROptimal);
        mod->translateZ(newZPos);

        m_zPlusModules.push_back(mod);

        // Update info
        lastZPos       = newZPos;
        lastDsDistance = newDsDistance;
        parity        *= -1;
        iMod++;
      }

      //
      // Build modules in negative Z -> balancing can be used only if outerZ defined, not if numModules specified!

      // For module-edge option, number of negative modules given by balancing algorithm -> usable only if outerZ defined, not if the number of modules specified
      if (useBalancing()) targetMods = (m_startZMode()==StartZMode::MODULECENTER) ? m_zPlusModules.size()-1 : std::numeric_limits<int>::max();
      else                targetMods = (m_startZMode()==StartZMode::MODULECENTER) ? m_zPlusModules.size()-1 : m_zPlusModules.size();

      parity     = -m_smallParity;
      iMod       = 1;
      newZPos    = (m_startZMode()==StartZMode::MODULECENTER) ? 0.0 : startZPos + moduleLength()/2. - zUnbalance;
      lastZPos   = newZPos;
      iIter      = 0;
      //std::cout << "TargeMods " << targetMods << std::endl;
      while ((iMod<=targetMods && (fabs(newZPos)+moduleLength()/2.)<targetZ) || (!useBalancing() && iMod<=targetMods)) {

        iIter++;
        if (iIter>RodPairStraight::RodPairStraight::c_nIterations) logERROR("When positioning modules in negative Z, number of iterations exceeded allowed limit! Quitting!!!");

        BarrelModule* mod = GeometryFactory::make<BarrelModule>(iMod, GeometryFactory::make<RectangularModule>(), m_ringNode, propertyTree());
        mod->build();

        // Translate & rotate (rotate after translation!)
        newRPos  = (parity > 0 ? m_maxRadius + m_smallDelta : m_minRadius - m_smallDelta) - newDsDistance/2; // Ds>0 lower sensor plays a role
        lastRPos = (parity > 0 ? m_maxRadius - m_smallDelta : m_minRadius + m_smallDelta) + newDsDistance/2; // Ds>0 upper sensor plays a role
        newZPos  = computeNextZ(lastRPos, newRPos, lastZPos, moduleLength(), lastDsDistance, mod->dsDistance(), BuildDir::LEFT, parity);

        newROptimal  = m_bigParity > 0 ? m_bigDelta : -m_bigDelta;
        newROptimal += (parity > 0 ? m_optimalRadius + m_smallDelta : m_optimalRadius - m_smallDelta);

        //std::cout << "Module Z: " << newZPos << std::endl;
        mod->translateR(newROptimal);
        mod->translateZ(newZPos);

        m_zMinusModules.push_back(mod);

        // Update info
        lastZPos       = newZPos;
        lastDsDistance = newDsDistance;
        parity        *= -1;
        iMod++;
      }

      // Unbalance in +-Z & starting point -> recalculate
      if (m_zPlusModules.size()>0 && m_zMinusModules.size()>0) {

        startZPos      = (m_startZMode()==StartZMode::MODULECENTER) ? 0.0 : m_zPlusModules.front().center().Z() - m_zPlusModules.front().length()/2.;
        double maxZPls = m_zPlusModules.size()>1 ? m_zPlusModules[m_zPlusModules.size()-2].center().Z()+m_zPlusModules[m_zPlusModules.size()-2].length()/2. : 0;
        maxZPls        = MAX(m_zPlusModules[m_zPlusModules.size()-1].center().Z()+m_zPlusModules[m_zPlusModules.size()-1].length()/2., maxZPls);
        double minZMin = m_zMinusModules.size()>1 ? m_zMinusModules[m_zMinusModules.size()-2].center().Z()-m_zMinusModules[m_zMinusModules.size()-2].length()/2. : 0;
        minZMin        = MIN(m_zMinusModules[m_zMinusModules.size()-1].center().Z()-m_zMinusModules[m_zMinusModules.size()-1].length()/2., minZMin);
        zUnbalance     = (maxZPls + minZMin)/2.; // balancing uneven pos/neg stringsdouble

        //std::cout << "Modules: " << m_zPlusModules.size() << " + " << m_zMinusModules.size() << std::endl;
        //if (this->myid()==1) std::cout << ">Balancing +Z x -Z> " << maxZPls << " "<< minZMin << " diff: " << zUnbalance << " " << startZPos << " " << startZPos + moduleLength()/2. - zUnbalance << std::endl;
      }
      else {

        startZPos  = 0.0;
        zUnbalance = 0.0;
      }


      if (!useBalancing()) break; // Do build iteration only once if balancing not needed
      bIter++;                    // Number of allowed iterations in balancing procedure
    } // Balancing

    // Print warning if balancing failed
    if (bIter>=c_nIterations) {

      std::ostringstream tempSS;
      tempSS << "Balancing algorighm in rod pair at avg radius: " << m_optimalRadius << " didn't converge! Layer is skewed";
      tempSS << "Unbalance is " << zUnbalance << " mm";
      logWARNING(tempSS);
    }
    // Print balancing output if balancing performed
    else if (bIter>1) {

      std::ostringstream tempSS;
      tempSS << "Balancing algorithm in rod pair at avg radius " << m_optimalRadius << " converged after " << bIter << " step(s).\n"
             << "   Residual Z unbalance is " << zUnbalance << ".\n"
             << "   Positive string has " << m_zPlusModules.size() << " modules, negative string has " << m_zMinusModules.size() << " modules.\n"
             << "   Z+ rod starts at " << m_zPlusModules.front().center().Z() - m_zPlusModules.front().length()/2.
             << ", Z- rod starts at "  << m_zMinusModules.front().center().Z() + m_zMinusModules.front().length()/2. << ".";
      logINFO(tempSS);
    }

    //
    // Solve colisions in ZPlus
    std::map<int, double> zGuards;
    zGuards[1]      = std::numeric_limits<double>::lowest();
    zGuards[-1]     = std::numeric_limits<double>::lowest();
    double parity   = m_smallParity;
    int iModule     = 1;
    int nCollisions = 0;

    logINFO("Checking for collisions in Z+ rod");

    for (auto& m : m_zPlusModules) {

      double minPhysZ = m.center().Z() - MAX(m.physicalLength(), m.length())/2;
      if (zGuards[parity] > minPhysZ) {

        double offset = zGuards[parity] - minPhysZ;
        m.translateZ(offset);
        nCollisions++;
        logINFO("  Module " + any2str(iModule) + " collides with previous " + (parity > 0 ? "outer" : "inner") + " module. Translated by " + any2str(offset) + " mm");
      }
      zGuards[parity] = m.center().Z() + MAX(m.physicalLength(), m.length())/2;
      parity = -parity;
      iModule++;
    }
    if (nCollisions>0) logWARNING("Some modules in positive Z have been translated to avoid collisions. Check info tab");

    //
    // Solve colisions in ZMinus
    zGuards[1]  = std::numeric_limits<double>::max();
    zGuards[-1] = std::numeric_limits<double>::max();
    parity      = -m_smallParity;
    iModule     = -1;
    nCollisions = 0;

    logINFO("Checking for collisions in Z- rod");

    for (auto& m : m_zMinusModules) {

      double maxPhysZ = m.center().Z() + MAX(m.physicalLength(), m.length())/2;
      if (zGuards[parity] < maxPhysZ) {

        double offset = zGuards[parity] - maxPhysZ;
        m.translateZ(offset);
        nCollisions++;
        logINFO("  Module " + any2str(iModule) + " collides with previous " + (parity > 0 ? "outer" : "inner") + " module. Translated by " + any2str(offset) + " mm");
      }
      zGuards[parity] = m.center().Z() - MAX(m.physicalLength(), m.length())/2;
      parity = -parity;
      iModule--;
    }
    if (nCollisions>0) logWARNING("Some modules in negative Z have been translated to avoid collisions. Check info tab");

    //
    // Compress modules if outerZ defined (to get maximum z position, use reverse operators)
    double curMaxZ = (m_zPlusModules.size()>1) ? MAX(m_zPlusModules.rbegin()->planarMaxZ(),(m_zPlusModules.rbegin()+1)->planarMaxZ()) : (!m_zPlusModules.empty() ? m_zPlusModules.rbegin()->planarMaxZ() : 0.);
    if (useCompression() && (m_nModules==0) && (curMaxZ>m_outerZ)) compressToZ(m_outerZ);

    //
    // Rotate all modules -> needs to follow after translation!!!
    for (auto& m : m_zPlusModules)  { m.rotateZ(m_rotation); }
    for (auto& m : m_zMinusModules) { m.rotateZ(m_rotation); }

    //
    // Set geometry info needed to build tilted part (if required to be built)
    setGeometryInfo();
  }
  catch (PathfulException& pe) {

    pe.pushPath(fullid(*this));
    throw;
  }

  cleanup();
  builtok(true);
}

//
// Helper method calculating module optimal Z position taking into account small/bigDelta, beam spot sizes etc. -> to fully cover eta region
//
double RodPairStraight::computeNextZ(double lastRPos, double newRPos, double lastZPos, double moduleLength, double lastDsDistance, double newDsDistance, BuildDir direction, int parity)
{
  // Update maximum, minimum radius using dsDistance values (Transform given Z positions from central Z to its edge equivalent -> one calculates optimal new Z pos related to edge coverage in eta)
  double lastZ = (direction==BuildDir::RIGHT ?  lastZPos + moduleLength/2. : lastZPos - moduleLength/2.);
  double newZ  = lastZ;

  // Cover beam spot
  double dz = zError();

  // Building in positive Z direction
  if (direction == BuildDir::RIGHT) {

    // For +smallDelta a worse case is +dz, for -smallDelta -dz
    double originZ = parity > 0 ? dz : -dz;

    // Take worse from 2 effects: beam spot size or required overlap in Z
    double newZorigin  = (newZ - zOverlap()) * newRPos/lastRPos;
    double newZshifted = (newZ - originZ) * newRPos/lastRPos + originZ;
    if (SimParms::getInstance().useLumiRegInGeomBuild()) newZ = MIN(newZorigin, newZshifted);
    else                                                 newZ = newZorigin;

    // Transfor back to module central position
    newZ = newZ + moduleLength/2.;
  }

  // Building in negative Z direction
  else {

    // For +smallDelta a worse case is -dz, for -smallDelta +dz
    double originZ = parity > 0 ? -dz : dz;

    // Take worse from 2 effects: beam spot size or required overlap in Z
    double newZorigin  = (newZ + zOverlap()) * newRPos/lastRPos;
    double newZshifted = (newZ - originZ) * newRPos/lastRPos + originZ;
    if (SimParms::getInstance().useLumiRegInGeomBuild()) newZ = MAX(newZorigin, newZshifted);
    else                                                 newZ = newZorigin;

    // Transfor back to module central position
    newZ = newZ - moduleLength/2.;
  }

  // Return centre agaion
  return newZ;
}

//
// Helper method compressing modules in Z direction (symmetrically from minus/plus Z) such as they fit into -outerZ, +outerZ region
//
void RodPairStraight::compressToZ(double zLimit) {

  if (m_zPlusModules.empty()) return;

  logINFO("Layer slated for compression");

  double findMaxZModule = -std::numeric_limits<double>::max();
  double findMinZModule = +std::numeric_limits<double>::max();

  for (auto& m : m_zPlusModules)  if (m.planarMaxZ()>findMaxZModule) findMaxZModule = m.planarMaxZ();
  for (auto& m : m_zMinusModules) if (m.planarMinZ()<findMinZModule) findMinZModule = m.planarMinZ();

  double Deltap =  fabs(zLimit) - findMaxZModule;
  double Deltam = -fabs(zLimit) - findMinZModule;

  logINFO("Rod length cannot exceed " + any2str(zLimit));
  logINFO("  Z+ rod will be compressed by " + any2str(fabs(Deltap)));
  logINFO("  Z- rod will be compressed by " + any2str(fabs(Deltam)));

  // Cut modules beyond allowed maximum Z: outerZ
  if (allowCompressionCuts()) {

    logINFO("Compression algorithm is allowed to cut modules which fall out entirely from the maximum z line");

    //
    // Cut out modules in positive Z
    int zPlusOrigSize = m_zPlusModules.size();
    for (auto it = m_zPlusModules.begin(); it != m_zPlusModules.end();) {

      if (it->planarMinZ()>zLimit) it = m_zPlusModules.erase(it);
      else ++it;
    }
    logINFO("  " + any2str(zPlusOrigSize - m_zPlusModules.size()) + " modules were cut from the Z+ rod");

    if (zPlusOrigSize - m_zPlusModules.size()) {

      findMaxZModule = -std::numeric_limits<double>::max();
      for (auto& m : m_zPlusModules)  if (m.planarMaxZ()>findMaxZModule) findMaxZModule = m.planarMaxZ();

      Deltap = fabs(zLimit) - findMaxZModule;  // we have to use findMaxZModule instead of checking only the last module as there might be inversions at high Z
      logINFO("  Z+ rod now exceeding by " + any2str(fabs(Deltap)));
    }

    //
    // Cut out modules in negative Z
    int zMinusOrigSize = m_zMinusModules.size();
    for (auto it = m_zMinusModules.begin(); it != m_zMinusModules.end();) {

      if (it->planarMaxZ()<-zLimit) it = m_zMinusModules.erase(it);
      else ++it;
    }
    logINFO("  " + any2str(zMinusOrigSize - m_zMinusModules.size()) + " modules were cut from the Z- rod");

    if (zMinusOrigSize - m_zMinusModules.size()) {

      findMinZModule = std::numeric_limits<double>::max();
      for (auto& m : m_zMinusModules) if (m.planarMinZ()<findMinZModule) findMinZModule = m.planarMinZ();

      Deltam = -fabs(zLimit) - findMinZModule;
      logINFO("  Z- rod now exceeding by " + any2str(fabs(Deltam)));
    }
  }

  //
  // Do compression in positive Z
  logINFO("Iterative compression of Z+ rod");
  int iIter=0;

  for (iIter = 0; fabs(Deltap)>c_safetySpaceFactor && iIter<RodPairStraight::c_nIterations; iIter++) {

    double maxCentreZ = -std::numeric_limits<double>::max();
    for (auto& m : m_zPlusModules)  if (m.center().Z()>maxCentreZ) maxCentreZ = m.center().Z();

    // Relative fraction: offset to central Z position of uppermost module in a rod
    double deltap = Deltap/maxCentreZ;
    int    parity = m_smallParity;

    std::map<int, double> zGuards;
    zGuards[1]  = std::numeric_limits<double>::lowest();
    zGuards[-1] = std::numeric_limits<double>::lowest();

    for (auto it = m_zPlusModules.begin(); it<m_zPlusModules.end(); ++it, parity = -parity) {

      // Correction shift for each module
      double translation = deltap*it->center().Z();

      // If first module put to the centre -> avoid any shift
      if (it==m_zPlusModules.begin() && m_startZMode()==StartZMode::MODULECENTER) translation = 0.0;

      // Calculate that new lower Z position of current module so that it doesn't collide with the upper Z position of preceding module
      double minPhysZ    = MIN(it->planarMinZ(), it->center().Z() - it->physicalLength()/2);
      if (minPhysZ + translation<zGuards[parity]) translation = zGuards[parity] - minPhysZ;

      // Shift the module
      it->translateZ(translation);

      // Calculate safe boundary (new upper Z position) to avoid collision with the next module
      double maxPhysZ = MAX(it->planarMaxZ(), it->center().Z() + it->physicalLength()/2);
      zGuards[parity] = maxPhysZ;
    }

    double newMaxZ = -std::numeric_limits<double>::max();
    for (auto& m : m_zPlusModules)  if (m.planarMaxZ()>newMaxZ) newMaxZ = m.planarMaxZ();

    Deltap =  fabs(zLimit) - newMaxZ;
  }
  if (iIter==RodPairStraight::c_nIterations) {
    logINFO(   "Iterative compression didn't terminate after " + any2str(RodPairStraight::c_nIterations) + " iterations. Z+ rod still exceeds by " + any2str(fabs(Deltap)));
    logWARNING("Failed to compress Z+ rod. Still exceeding by " + any2str(fabs(Deltap)) + ". Check info tab.");
  }
  else logINFO("Z+ rod successfully compressed after " + any2str(iIter) + " iterations. Rod now only exceeds by " + any2str(fabs(Deltap)) + " mm.");

  //
  // Do compression in negative Z
  logINFO("Iterative compression of Z- rod");
  iIter=0;

  for (iIter = 0; fabs(Deltam) > c_safetySpaceFactor && iIter < RodPairStraight::c_nIterations; iIter++) {

    double minCentreZ = std::numeric_limits<double>::max();
    for (auto& m : m_zMinusModules)  if (m.center().Z()<minCentreZ) minCentreZ = m.center().Z();

    // Relative fraction: offset to central Z position of lowermost module in a rod
    double deltam = Deltam/minCentreZ; // this can be optimized
    int    parity = -m_smallParity;

    std::map<int, double> zGuards;
    zGuards[1] = std::numeric_limits<double>::max();
    zGuards[-1] = std::numeric_limits<double>::max();

    for (auto it=m_zMinusModules.begin(); it<m_zMinusModules.end(); ++it, parity = -parity) {

      // Correction shift for each module
      double translation = deltam*it->center().Z();

      // Calculate that new upper Z position of current module so that it doesn't collide with the lower Z position of preceding module
      double maxPhysZ = MAX(it->planarMaxZ(), it->center().Z() + it->physicalLength()/2);
      if (maxPhysZ + translation > zGuards[parity]) translation = zGuards[parity] - maxPhysZ;

      // Shift the module
      it->translateZ(translation);

      // Calculate safe boundary (new lower Z position) to avoid collision with the next module
      double minPhysZ = MIN(it->planarMinZ(), it->center().Z() - it->physicalLength()/2);
      zGuards[parity] = minPhysZ;
    }

    double newMinZ = std::numeric_limits<double>::max();
    for (auto& m : m_zMinusModules)  if (m.planarMinZ()<newMinZ) newMinZ = m.planarMinZ();

    Deltam = -fabs(zLimit) - newMinZ;
  }
  if (iIter==RodPairStraight::c_nIterations) {
    logINFO("Iterative compression didn't terminate after " + any2str(RodPairStraight::c_nIterations) + " iterations. Z- rod still exceeds by " + any2str(fabs(Deltam)));
    logWARNING("Failed to compress Z- rod. Still exceeding by " + any2str(fabs(Deltam)) + ". Check info tab.");
  }
  else logINFO("Z- rod successfully compressed after " + any2str(iIter) + " iterations. Rod now only exceeds by " + any2str(fabs(Deltam)) + " mm.");

}

//! Helper method calculating geometry properties needed for building the tilted part of the rod after the straight rod: thetaEnd, phiOverlap, rEndInner, rEndOuter
void RodPairStraight::setGeometryInfo() {

  // Calculate theta end
  if (m_zPlusModules.empty()) { thetaEnd(M_PI/2.); }
  else {

    // Find last module
    auto lastMod = m_zPlusModules.back();

    double dsDistance = lastMod.dsDistance();
    double lastR      = lastMod.center().Rho();

    double rH2ppUP = lastR + 0.5 * dsDistance;  // WARNING !!! FOR THE MOMENT, DOESN T TAKE MODULE WIDTH INTO ACCOUNT, SHOULD BE CHANGED ?

    thetaEnd( atan(rH2ppUP/(lastMod.flatMaxZ())) );
  }
}

//
//  Constructor - parse geometry config file using boost property tree & read-in Rod parameters
//
RodPairTilted::RodPairTilted(int id, double rotation, bool outerRadiusRod, const PropertyTree& treeProperty) :
 RodPair(id, rotation, treeProperty),
 isOuterRadiusRod("isOuterRadiusRod", parsedAndChecked(), outerRadiusRod)
{
  // Set the geometry config parameters (id set through base class)
  this->store(treeProperty);
}

//
// Cross-check parameters provided from geometry configuration file
//
void RodPairTilted::check() {
  PropertyObject::check();

  if (m_startZMode() != StartZMode::MODULECENTER) throw PathfulException("Tilted layer : only startZMode = modulecenter can be specified.");
}

//
// Build tilted rods based on
//
void RodPairTilted::build(const RodTemplate& rodTemplate, const TiltedRings& tiltedRings, bool isFlipped) {

  m_materialObject.store(propertyTree());
  m_materialObject.build();

  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();
    buildModules(m_zPlusModules , rodTemplate, tiltedRings, BuildDir::RIGHT, isFlipped);
    buildModules(m_zMinusModules, rodTemplate, tiltedRings, BuildDir::LEFT , isFlipped);

    //
    // Rotate all modules -> needs to follow after translation!!!
    for (auto& m : m_zPlusModules)  { m.rotateZ(m_rotation); }
    for (auto& m : m_zMinusModules) { m.rotateZ(m_rotation); }
  }
  catch (PathfulException& pe) {

    pe.pushPath(fullid(*this));
    throw;
  }

  cleanup();
  builtok(true);
}

//
// Build modules in given direction (+Z or -Z)
//
void RodPairTilted::buildModules(BarrelModules& modules, const RodTemplate& rodTemplate, const TiltedRings& tiltedRings, BuildDir direction, bool isFlipped) {

  // Build if not empty
  if (!tiltedRings.empty()) {

    int iMod = 0;
    auto iterMod = rodTemplate.begin();
    for (; iMod<tiltedRings.size(); iMod++, ++iterMod) {

      BarrelModule* mod = GeometryFactory::make<BarrelModule>(*iterMod);
      mod->myid(iMod+1);
      mod->side(short(direction));
      mod->tilt(tiltedRings[iMod].tiltAngle() * short(direction));
      if (isOuterRadiusRod()) mod->translateR(tiltedRings[iMod].outerRadius());
      else                    mod->translateR(tiltedRings[iMod].innerRadius());
      mod->flipped(isFlipped);
      if (isOuterRadiusRod()) mod->translateZ(tiltedRings[iMod].zOuter() * short(direction));
      else                    mod->translateZ(tiltedRings[iMod].zInner() * short(direction));

      modules.push_back(mod);
    }
  }
}

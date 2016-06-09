#include "SimParms.h"

#include <MainConfigHandler.h>
#include "global_constants.h"
#include "Visitor.h"
#include "IrradiationMapsManager.h"

//
// Static instance of this class
//
SimParms * SimParms::s_instance = nullptr;

//
// Method constructing static instance of this class
//
SimParms * SimParms::getInstance()
{
  if (s_instance == nullptr) s_instance = new SimParms;
  return s_instance;
}

//
// Private constructor -> needed for singleton pattern
//
SimParms::SimParms() :
 numMinBiasEvents(        "numMinBiasEvents"        , parsedAndChecked()),
 zErrorCollider(          "zErrorCollider"          , parsedAndChecked()),
 rError(                  "rError"                  , parsedAndChecked()),
 useIPConstraint(         "useIPConstraint"         , parsedAndChecked()),
 ptCost(                  "ptCost"                  , parsedAndChecked()),
 stripCost(               "stripCost"               , parsedAndChecked()),
 efficiency(              "efficiency"              , parsedAndChecked()),
 pixelEfficiency(         "pixelEfficiency"         , parsedAndChecked()),
 bunchSpacingNs(          "bunchSpacingNs"          , parsedAndChecked()),
 triggerEtaCut(           "triggerEtaCut"           , parsedOnly(), 2),
 triggerPtCut(            "triggerPtCut"            , parsedOnly(), 1),
 numTriggerTowersEta(     "numTriggerTowersEta"     , parsedOnly(), 1),
 numTriggerTowersPhi(     "numTriggerTowersPhi"     , parsedOnly(), 1),
 timeIntegratedLumi(      "timeIntegratedLumi"      , parsedOnly(), 3000),
 operatingTemp(           "operatingTemp"           , parsedOnly(), -20),
 chargeDepletionVoltage(  "chargeDepletionVoltage"  , parsedOnly(), 600),
 alphaParm(               "alphaParm"               , parsedOnly(), 4e-17),
 referenceTemp(           "referenceTemp"           , parsedOnly(), 20),
 magneticField(           "magneticField"           , parsedOnly(), magnetic_field),
 dipoleMagneticField(     "dipoleMagneticField"     , parsedOnly(), 0.0),
 dipoleDPlResAt10TeV(     "dipoleDPlResAt10TeV"     , parsedOnly(), 0.1),
 dipoleXToX0(             "dipoleXToX0"             , parsedOnly(), 0.1),
 irradiationMapFiles(     "irradiationMapFiles"     , parsedAndChecked()),
 bFieldMapFile(           "bFieldMapFile"           , parsedOnly(), std::string("")),
 chargedMapFile(          "chargedMapFile"          , parsedOnly(), std::string("")),
 chargedMapLowThFile(     "chargedMapLowThFile"     , parsedOnly(), std::string("")),
 chargedMapBOffMatOnFile( "chargedMapBOffMatOnFile" , parsedOnly(), std::string("")),
 chargedMapBOnMatOffFile( "chargedMapBOnMatOffFile" , parsedOnly(), std::string("")),
 chargedMapBOffMatOffFile("chargedMapBOffMatOffFile", parsedOnly(), std::string("")),
 chargedMapBOffTrkOffFile("chargedMapBOffTrkOffFile", parsedOnly(), std::string("")),
 photonsMapFile(          "photonsMapFile"          , parsedOnly(), std::string("")),
 photonsMapLowThFile(     "photonsMapLowThFile"     , parsedOnly(), std::string("")),
 photonsMapBOffMatOnFile( "photonsMapBOffMatOnFile" , parsedOnly(), std::string("")),
 photonsMapBOnMatOffFile( "photonsMapBOnMatOffFile" , parsedOnly(), std::string("")),
 photonsMapBOffMatOffFile("photonsMapBOffMatOffFile", parsedOnly(), std::string("")),
 photonsMapBOffTrkOffFile("photonsMapBOffTrkOffFile", parsedOnly(), std::string("")),
 minTracksEta(            "minTracksEta"            , parsedOnly()),
 maxTracksEta(            "maxTracksEta"            , parsedOnly()),
 taggedTracking(          "TaggedTracking"          , parsedOnly()),
 m_commandLine(""),
 m_geomFile(""),
 m_baseDir(""),
 m_htmlDir(""),
 m_layoutName("Default")
{
  // Read irradiation maps
  m_irradiationMapsManager = new IrradiationMapsManager();
  for (auto file : default_irradiationfiles) {

    std::string path = MainConfigHandler::getInstance().getIrradiationDirectory() + "/" + file;
    irradiationMapFiles.appendString(path);
    m_irradiationMapsManager->addIrradiationMap(path.c_str());
  }
}

//
// Destructor - clear memory
//
SimParms::~SimParms()
{
  if (m_irradiationMapsManager != nullptr) delete m_irradiationMapsManager;
}

//
// Check that sim parameters read-in correctly
//
void SimParms::crosscheck() {
  try {
    check();
    builtok(true);
    cleanup();
  } catch (PathfulException& pe) { pe.pushPath("SimParms"); throw; }
}

//
// SimParms visitable -> implemented accept method to call visitor pattern
//
void SimParms::accept(GeometryVisitor& v) { v.visit(*this); }

//
// SimParms visitable -> implememented const accept method to call visitor pattern
//
void SimParms::accept(ConstGeometryVisitor& v) const { v.visit(*this); }

//
// Set command line options passed over to program to analyze data
//
void SimParms::setCommandLine(int argc, char* argv[])
{
  m_commandLine = argv[1];
  for (auto i=2; i<argc; i++) m_commandLine += std::string(" ") + argv[i];
}

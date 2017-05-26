#include "SimParms.h"

#include <MainConfigHandler.h>
#include "global_constants.h"
#include "Units.h"
#include "Visitor.h"
#include "IrradiationMapsManager.h"

//
// Method constructing static instance of this class
//
SimParms& SimParms::getInstance()
{
  static SimParms s_instance;
  return s_instance;
}

//
// Private constructor -> needed for singleton pattern
//
SimParms::SimParms() :
 numMinBiasEvents(        "numMinBiasEvents"        , parsedAndChecked()),
 zErrorIP(                "zErrorIP"                , parsedAndChecked()),
 rphiErrorIP(             "rphiErrorIP"             , parsedAndChecked()),
 useIPConstraint(         "useIPConstraint"         , parsedAndChecked()),
 useLumiRegInGeomBuild(   "useLumiRegInGeomBuild"   , parsedOnly(), 1),
 useLumiRegInAnalysis(    "useLumiRegInAnalysis"    , parsedOnly(), 0),
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
 magField(                "magField"                , parsedOnly()),
 magFieldZRegions(        "magFieldZRegions"        , parsedOnly()),
 irradiationMapFiles(     "irradiationMapFiles"     , parsedAndChecked()),
 etaRegionRanges(         "etaRegionRanges"         , parsedAndChecked()),
 etaRegionNames(          "etaRegionNames"          , parsedAndChecked()),
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
  m_irradiationMapsManager = std::unique_ptr<IrradiationMapsManager>(new IrradiationMapsManager());
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
  m_irradiationMapsManager.reset();
}

//
// Check that sim parameters read-in correctly & set units
//
void SimParms::crosscheck() {
  try {
    check();
    builtok(true);
    cleanup();
  } catch (PathfulException& pe) { pe.pushPath("SimParms"); throw; }

  // Check that sizes of eta vectors: names & ranges are the same
  if (etaRegionNames.size()!=etaRegionRanges.size()) throw PathfulException("Number of names assigned to individual eta regions don't correspond to number of defined eta regions, check the config file!" , "SimParms");

  // Check that number of magnetic field values correspond to number of magnetic field intervals
  if (magField.size()!=magFieldZRegions.size()) throw PathfulException("Number of magnetic field values doesn't correspond to number of intervals!" , "SimParms");

  // Check that magnetic field non-zero in at least one interval
  bool nonZero = false;
  for (auto i=0; i<magFieldZRegions.size(); i++) {
    if (magField[i]>0) nonZero = true;
  }
  if (nonZero!=true) throw PathfulException("Magnetic field needs to be defined non-zero in at least one Z interval!" , "SimParms");

  // Check that IP errors non-zero in R-Phi & Z if IP constraint will be used
  nonZero = true;
  if (useIPConstraint()) {

    if (zErrorIP()==0)    nonZero = false;
    if (rphiErrorIP()==0) nonZero = false;
  }
  if (!nonZero) throw PathfulException("IP constraint required, but errors on beam spot set to zero in R-Phi/Z!" , "SimParms");

  // Set expected default units
  magField.scaleByUnit(Units::T);
  magFieldZRegions.scaleByUnit(Units::m);

  zErrorIP.scaleByUnit(Units::mm);
  rphiErrorIP.scaleByUnit(Units::mm);

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

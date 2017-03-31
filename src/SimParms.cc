#include "SimParms.hh"
#include "Units.hh"
#include "MainConfigHandler.hh"

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
    numMinBiasEvents("numMinBiasEvents", parsedAndChecked()),
    bunchSpacingNs("bunchSpacingNs", parsedAndChecked()),
    zErrorCollider("zErrorCollider", parsedAndChecked()),
    rphiErrorCollider("rphiErrorCollider", parsedAndChecked()),
    useIPConstraint("useIPConstraint", parsedAndChecked()),
    ptCost("ptCost", parsedAndChecked()),
    stripCost("stripCost", parsedAndChecked()),
    triggerEtaCut("triggerEtaCut", parsedAndChecked()),
    triggerPtCut("triggerPtCut", parsedAndChecked()),
    numTriggerTowersEta("numTriggerTowersEta", parsedAndChecked()),
    numTriggerTowersPhi("numTriggerTowersPhi", parsedAndChecked()),
    timeIntegratedLumi("timeIntegratedLumi", parsedAndChecked()),
    referenceTemp("referenceTemp", parsedAndChecked()),
    alphaParam("alphaParam", parsedAndChecked()),    // radiation-damage coefficient, A/cm
    magField("magneticField", parsedAndChecked()),
    irradiationMapFiles("irradiationMapFiles", parsedAndChecked()),
    minTracksEta("minTracksEta", parsedOnly()),
    maxTracksEta("maxTracksEta", parsedOnly()),
    taggedTracking("TaggedTracking", parsedOnly())
{ }

void SimParms::build() {
  try {
    check();

    std::string irradiationMapDirectory_ = mainConfigHandler::instance().getIrradiationDirectory();
    //iter between irradiation map file names and feed the irradiationMapsManager
    for(std::vector<std::string>::const_iterator iterMapFile = irradiationMapFiles.begin(); iterMapFile != irradiationMapFiles.end(); ++ iterMapFile) {
      std::string fullPath = irradiationMapDirectory_ + "/" + (*iterMapFile);
      irradiationMapsManager_.addIrradiationMap(fullPath.c_str());
    }

    cleanup();
    builtok(true);
  } catch (PathfulException& pe) { pe.pushPath("SimParms"); throw; }

  // Check that IP errors non-zero in R-Phi & Z if IP constraint will be used
  bool nonZero = true;
  if (useIPConstraint()) {

    if (zErrorCollider()==0)    nonZero = false;
    if (rphiErrorCollider()==0) nonZero = false;
  }
  if (!nonZero) throw PathfulException("IP constraint required, but errors on beam spot set to zero in R-Phi/Z!" , "SimParms");


  // Set expected default units
  magField.scaleByUnit(Units::T);

  zErrorCollider.scaleByUnit(Units::mm);
  rphiErrorCollider.scaleByUnit(Units::mm);
}

void SimParms::addIrradiationMapFile(std::string path) {
  irradiationMapFiles.appendString(path);
}


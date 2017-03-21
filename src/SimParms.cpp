#include "SimParms.h"
#include "Units.h"

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
    rError("rError", parsedAndChecked()),
    useIPConstraint("useIPConstraint", parsedAndChecked()),
    ptCost("ptCost", parsedAndChecked()),
    stripCost("stripCost", parsedAndChecked()),
    efficiency("efficiency", parsedAndChecked()),
    pixelEfficiency("pixelEfficiency", parsedAndChecked()),
    triggerEtaCut("triggerEtaCut", parsedAndChecked()),
    triggerPtCut("triggerPtCut", parsedAndChecked()),
    numTriggerTowersEta("numTriggerTowersEta", parsedAndChecked()),
    numTriggerTowersPhi("numTriggerTowersPhi", parsedAndChecked()),
    timeIntegratedLumi("timeIntegratedLumi", parsedAndChecked()),
    referenceTemp("referenceTemp", parsedAndChecked()),
    alphaParam("alphaParam", parsedAndChecked()),    // radiation-damage coefficient, A/cm
    magField("magneticField", parsedAndChecked()),
    pileUp("pileUp", parsedAndChecked()),
    irradiationMapFiles("irradiationMapFiles", parsedAndChecked()),
    minTracksEta("minTracksEta", parsedOnly()),
    maxTracksEta("maxTracksEta", parsedOnly()),
    taggedTracking("TaggedTracking", parsedOnly())
{ }

void SimParms::build() {
  try {
    check();

    //iter between irradiation map file names and feed the irradiationMapsManager
    for(std::vector<std::string>::const_iterator iterMapFile = irradiationMapFiles.begin(); iterMapFile != irradiationMapFiles.end(); ++ iterMapFile) {
      irradiationMapsManager_.addIrradiationMap((*iterMapFile).c_str());
    }

    cleanup();
    builtok(true);
  } catch (PathfulException& pe) { pe.pushPath("SimParms"); throw; }


  // Set expected default units
  magField.scaleByUnit(Units::T);
}

void SimParms::addIrradiationMapFile(std::string path) {
  irradiationMapFiles.appendString(path);
}


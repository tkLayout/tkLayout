#include "SimParms.hh"
#include "Units.hh"
#include "MainConfigHandler.hh"

define_enum_strings(LumiRegShape) = { "ponctual", "flat", "gaussian" };

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

const XYZVector SimParms::getLuminousRegion() const {
  const double dx = lumiRegRError() / sqrt(2.);
  const double dy = dx;

  double dz = 0.;
  if (lumiRegShape() == LumiRegShape::PONCTUAL) dz = 0.;
  else if (lumiRegShape() == LumiRegShape::FLAT) dz = (myDice.Rndm() * 2. - 1.) * lumiRegZError();
  else if (lumiRegShape() == LumiRegShape::GAUSSIAN) dz = myDice.Gaus(0, lumiRegZError());

  return XYZVector(dx, dy, dz); 
}

const XYZVector SimParms::getLuminousRegionInMatBudgetAnalysis() const {
  if (lumiRegShapeInMatBudgetAnalysis() != LumiRegShape::PONCTUAL) logERROR("Non-ponctual IP in Material Budget Analysis not supported.");
  return XYZVector(0., 0., 0.); 
}

void SimParms::addIrradiationMapFile(std::string path) {
  irradiationMapFiles.appendString(path);
}


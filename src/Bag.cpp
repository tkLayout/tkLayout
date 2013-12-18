#include "Bag.h"

const double graphBag::Triggerable     = 0.;
const int graphBag::RhoGraph         = 0x001;
const int graphBag::PhiGraph         = 0x002;
const int graphBag::DGraph           = 0x003;
const int graphBag::CtgthetaGraph    = 0x004;
const int graphBag::Z0Graph          = 0x005;
const int graphBag::PGraph           = 0x006;
const int graphBag::TriggeredGraph   = 0x007;
const int graphBag::IdealGraph       = 0x010;
const int graphBag::RealGraph        = 0x020;
const int graphBag::TriggerGraph     = 0x040;
const int graphBag::StandardGraph    = 0x080;
//const int graphBag::TriggerCorrelationGraph = 0x100;

const int mapBag::efficiencyMap         = 0x001;
const int mapBag::thresholdMap          = 0x002;
const int mapBag::thicknessMap          = 0x004;
const int mapBag::windowMap             = 0x008;
const int mapBag::suggestedSpacingMap   = 0x010;
const int mapBag::suggestedSpacingMapAW = 0x020;
const int mapBag::nominalCutMap      = 0x040;
const int mapBag::irradiatedPowerConsumptionMap = 0x080;
const int mapBag::totalPowerConsumptionMap = 0x100;
const int mapBag::moduleConnectionEtaMap = 0x200;
const int mapBag::moduleConnectionPhiMap = 0x400;
const int mapBag::moduleConnectionEndcapPhiMap = 0x800;

const double profileBag::Triggerable    = 0.;
const int profileBag::TriggeredProfile  = 0x0000007;
const int profileBag::TriggeredFractionProfile  = 0x0000008;
const int profileBag::TriggerPurityProfile = 0x0000010;
const int profileBag::TriggerProfile    = 0x0000040;

// These strings should be different from one another
// Also one should never be a substring of the other
const std::string profileBag::TriggerProfileName = "trigger";
const std::string profileBag::TriggerProfileNameWindow = "windowTrigger";
const std::string profileBag::TurnOnCurveName = "turnOnCurveTrigger";

const double mapBag::dummyMomentum = 0.;

int graphBag::clearTriggerGraphs() {
  return clearGraphs(graphBag::TriggerGraph);
}

int graphBag::clearStandardGraphs() {
  return clearGraphs(graphBag::StandardGraph);
}

int graphBag::clearGraphs(const int& attributeMask) {
  int deleteCounter = 0;
  for (auto it = taggedGraphMap_.begin(); it != taggedGraphMap_.end();) { // first clear tagged graphs with the specified attribute
    if ((it->first.first & attributeMask) == attributeMask) {
      it = taggedGraphMap_.erase(it);
      ++deleteCounter;
    } else ++it;
  }
  std::map<int, std::map<double, TGraph> >::iterator it;
  std::map<int, std::map<double, TGraph> >::iterator nextIt;

  int anAttribute;
  for (it=graphMap_.begin(); it!=graphMap_.end(); ) {
    anAttribute=it->first;
    if ((anAttribute&attributeMask)==attributeMask) {
      nextIt = ((++it)--);
      graphMap_.erase(it);
      it=nextIt;
      ++deleteCounter;
    } else {
      ++it;
    }
  }
  return deleteCounter;
}

int graphBag::buildAttribute(bool ideal, bool isTrigger) {
  int result;
  if (ideal) result = IdealGraph;
  else result = RealGraph;

  if (isTrigger) result |= TriggerGraph;
  else result |= StandardGraph;

  return result;
}

std::map<double, TGraph>& graphBag::getGraphs(const int& attribute) {
  return graphMap_[attribute];
}

std::map<double, TGraph>& graphBag::getTaggedGraphs(int attribute, const string& tag) {
  tagSet_.insert(tag);
  return taggedGraphMap_[std::make_pair(attribute, tag)];
}

std::map<double, TH2D>& mapBag::getMaps(const int& attribute) {
  return mapMap_[attribute];
}

int mapBag::clearMaps(const int& attributeMask) {
  std::map<int, std::map<double, TH2D> >::iterator it;
  std::map<int, std::map<double, TH2D> >::iterator nextIt;

  int deleteCounter = 0;
  int anAttribute;
  for (it=mapMap_.begin(); it!=mapMap_.end(); ) {
    anAttribute=it->first;
    if ((anAttribute&attributeMask)==attributeMask) {
      nextIt = it;
      ++nextIt;
      mapMap_.erase(it);
      it=nextIt;
      ++deleteCounter;
    } else {
      ++it;
    }
  }
  return deleteCounter;
}

int profileBag::clearTriggerProfiles() {
  return clearProfiles(profileBag::TriggerProfile);
}

int profileBag::clearTriggerNamedProfiles() {
  return clearNamedProfiles(profileBag::TriggerProfileName);
}


std::map<double, TProfile>& profileBag::getProfiles(const int& attribute) {
  return profileMap_[attribute];
}

// TODO: this looks like an invitation to use template classes Will
// do as soon as I have time :D (copy-paste worked till now...)
int profileBag::clearProfiles(const int& attributeMask) {
  std::map<int, std::map<double, TProfile> >::iterator it;
  std::map<int, std::map<double, TProfile> >::iterator nextIt;

  int deleteCounter = 0;
  int anAttribute;
  for (it=profileMap_.begin(); it!=profileMap_.end(); ) {
    anAttribute=it->first;
    if ((anAttribute&attributeMask)==attributeMask) {
      nextIt = it;
      ++nextIt;
      profileMap_.erase(it);
      it=nextIt;
      ++deleteCounter;
    } else {
      ++it;
    }
  }
  return deleteCounter;
}

int profileBag::clearNamedProfiles(const std::string& name) {
  std::map<std::string, std::map<double, TProfile> >::iterator it;
  std::map<std::string, std::map<double, TProfile> >::iterator nextIt;

  int deleteCounter = 0;
  std::string anAttribute;
  for (it=namedProfileMap_.begin(); it!=namedProfileMap_.end(); ) {
    anAttribute=it->first;
    if (anAttribute.substr(0, name.size())==name) {
      nextIt = it;
      ++nextIt;
      namedProfileMap_.erase(it);
      it=nextIt;
      ++deleteCounter;
    } else {
      ++it;
    }
  }
  return deleteCounter;    
}

std::vector<std::string> profileBag::getProfileNames(const std::string& name) {
  std::vector<std::string> result;

  std::map<std::string, std::map<double, TProfile> >::iterator it;

  std::string anAttribute;
  for (it=namedProfileMap_.begin(); it!=namedProfileMap_.end(); ++it) {
    anAttribute=it->first;
    if (anAttribute.substr(0, name.size())==name)
      result.push_back(anAttribute);
  }
  return result;    
}


std::map<double, TProfile>& profileBag::getNamedProfiles(const std::string& name) {
  return namedProfileMap_[name];
}


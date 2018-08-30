#include "OuterCabling/OuterBundle.hh"
#include "OuterCabling/OuterCable.hh"


OuterBundle::OuterBundle(const int id, const int stereoBundleId, const Category& type, const std::string subDetectorName, const int layerDiskNumber, const PhiPosition& phiPosition, const bool isPositiveCablingSide, const bool isTiltedPart) :
  type_(type),
  subDetectorName_(subDetectorName),
  layerDiskNumber_(layerDiskNumber),
  phiPosition_(phiPosition),
  isPositiveCablingSide_(isPositiveCablingSide),
  isTiltedPart_(isTiltedPart),
  stereoBundleId_(stereoBundleId)
{
  myid(id);
  plotColor_ = computePlotColor(id, isPositiveCablingSide);
};


const double OuterBundle::minPhi() const { 
  double min = std::numeric_limits<double>::max();
  for (const auto& m : modules_) { min = MIN(min, femodRounded(m->center().Phi(), 2. * M_PI) ); } return min;
}


const double OuterBundle::maxPhi() const { 
  double max = 0.;
  for (const auto& m : modules_) { max = MAX(max, femodRounded(m->center().Phi(), 2. * M_PI) ); } return max;
}

const double OuterBundle::meanPhi() const {
  std::vector<double> modPhis;

  for (const auto& m : modules_) { 
    double phi = femodRounded(m->center().Phi(), 2. * M_PI);
    if (modPhis.size() > 0 && (fabs(modPhis.back() - phi) > M_PI)) {
      if (phi < modPhis.back()) phi += 2.*M_PI;
      else phi -= 2.*M_PI;
    }
    modPhis.push_back(phi);
  } 

  double mean = 0.;
  for (const auto& phi : modPhis) { mean += phi; }
  mean /= numModules();
  return mean;
}


void OuterBundle::moveMinPhiModuleFromOtherBundle(OuterBundle* otherBundle) {
  Container& otherBundleModules = otherBundle->modules();
  const auto minPhiModuleIt = otherBundle->minPhiModule();

  std::move(minPhiModuleIt, minPhiModuleIt + 1, std::back_inserter(modules_));
  otherBundleModules.erase(minPhiModuleIt, minPhiModuleIt + 1);
}


void OuterBundle::moveMaxPhiModuleFromOtherBundle(OuterBundle* otherBundle) {
  Container& otherBundleModules = otherBundle->modules();
  const auto maxPhiModuleIt = otherBundle->maxPhiModule();

  std::move(maxPhiModuleIt, maxPhiModuleIt + 1, std::back_inserter(modules_));
  otherBundleModules.erase(maxPhiModuleIt, maxPhiModuleIt + 1);
}


const std::vector<Module*>::iterator OuterBundle::minPhiModule() {
  const auto modIt = std::min_element(modules_.begin(), modules_.end(), [](Module* a, Module* b) {
	return (femodRounded(a->center().Phi(), 2. * M_PI) <= femodRounded(b->center().Phi(), 2. * M_PI));
      });
  return modIt;
}


const std::vector<Module*>::iterator OuterBundle::maxPhiModule() {
  const auto modIt = std::max_element(modules_.begin(), modules_.end(), [](Module* a, Module* b) {
      return (femodRounded(a->center().Phi(), 2. * M_PI) <= femodRounded(b->center().Phi(), 2. * M_PI));
    });
  return modIt;
}


const int OuterBundle::computePlotColor(const int id, const bool isPositiveCablingSide) const {
  int plotColor = 0;
  int plotId = (isPositiveCablingSide ? id : (id - 20000));
  int plotType = 2 + plotId % 2;  // Barrel : Identifies Flat vs Tilted. Endcap : Identifies PS10GA vs PG10GB vs PS5G vs 2S type.
  int dizaine = plotId / 10;
  int plotPhi = dizaine % 3;  // Barrel : Identifies phiSegmentRef. Endcap : Identifies phiRegionRef.
  plotColor = plotType * 3 + plotPhi;
  return plotColor;
}


const ChannelSection* OuterBundle::opticalChannelSection() const {
  const ChannelSection* opticalSection = nullptr;
  if (cable_) opticalSection = cable_->opticalChannelSection();
  if (!cable_ || !opticalSection) throw PathfulException("cable_ or opticalSection is nullptr");
  return opticalSection;
}


/* Untilted TBPS only: Id of the bundle located on the same Phi, but connected to the tilted modules.
 */
const int OuterBundle::tiltedBundleId() const {
  if (!isBarrelPSFlatPart()) logERROR("Tried to access tiltedBundleId, but not in TBPS flat part.");
  const int bundleId = myid();
  const int tiltedBundleId = bundleId - femod(bundleId, 10);
  return tiltedBundleId;
}


/* Id of the bundle located on the other cabling side, by rotation of 180° around CMS_Y.
 */
const int OuterBundle::stereoBundleId() const {
  if (!isBarrel()) logERROR("Tried to access stereoBundleId in TEDD, where it is not implemented (only implemented for Barrel).");
  return stereoBundleId_;
}


void OuterBundle::setIsPowerRoutedToBarrelLowerSemiNonant(const bool isLower) {
  if (!isBarrel()) logERROR("Tried to set Barrel power routing on a Bundle located in TEDD.");
  isPowerRoutedToBarrelLowerSemiNonant_ = isLower;
}


const bool OuterBundle::isPowerRoutedToBarrelLowerSemiNonant() const {
  if (!isBarrel()) logERROR("Tried to access Barrel power routing from a Bundle located in TEDD.");
  return isPowerRoutedToBarrelLowerSemiNonant_;
}

#include "OuterCabling/Bundle.hh"
#include "OuterCabling/Cable.hh"


Bundle::Bundle(const int id, const int stereoBundleId, const Category& type, const std::string subDetectorName, const int layerDiskNumber, const PhiPosition& phiPosition, const bool isPositiveCablingSide, const bool isTiltedPart) :
  stereoBundleId_(stereoBundleId),
  type_(type),
  subDetectorName_(subDetectorName),
  layerDiskNumber_(layerDiskNumber),
  phiPosition_(phiPosition),
  isPositiveCablingSide_(isPositiveCablingSide),
  isTiltedPart_(isTiltedPart) 
{
  myid(id);
  plotColor_ = computePlotColor(id, isPositiveCablingSide);
};


Bundle::~Bundle() {
  delete cable_;    // TO DO: switch to smart pointers and remove this!
  cable_ = nullptr; 

  delete powerChannelSection_;
  powerChannelSection_ = nullptr;
}


void Bundle::moveMaxPhiModuleFromOtherBundle(Bundle* otherBundle) {
  Container& otherBundleModules = otherBundle->modules();
  auto maxPhiModuleIt = std::max_element(otherBundleModules.begin(), otherBundleModules.end(), [](const Module& a, const Module& b) {
      return (femodRounded(a.center().Phi(), 2. * M_PI) <= femodRounded(b.center().Phi(), 2. * M_PI));
    });

  modules_.transfer(modules_.end(),
		    maxPhiModuleIt,
		    otherBundleModules);
}


void Bundle::moveMinPhiModuleFromOtherBundle(Bundle* otherBundle) {
  Container& otherBundleModules = otherBundle->modules();
  auto minPhiModuleIt = std::min_element(otherBundleModules.begin(), otherBundleModules.end(), [](const Module& a, const Module& b) {
      return (femodRounded(a.center().Phi(), 2. * M_PI) <= femodRounded(b.center().Phi(), 2. * M_PI));
    });

  modules_.transfer(modules_.end(), 
		    minPhiModuleIt,
		    otherBundleModules);
}


const double Bundle::minPhi() const { 
  double min = std::numeric_limits<double>::max();
  for (const auto& m : modules_) { min = MIN(min, femodRounded(m.center().Phi(), 2. * M_PI) ); } return min;
}


const double Bundle::maxPhi() const { 
  double max = 0.;
  for (const auto& m : modules_) { max = MAX(max, femodRounded(m.center().Phi(), 2. * M_PI) ); } return max;
}

const double Bundle::meanPhi() const {
  std::vector<double> modPhis;

  for (const auto& m : modules_) { 
    double phi = femodRounded(m.center().Phi(), 2. * M_PI);
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

Module* Bundle::minPhiModule() const {
  const Module* mod = &(*std::min_element(modules_.begin(), modules_.end(), [](const Module& a, const Module& b) {
	return (femodRounded(a.center().Phi(), 2. * M_PI) <= femodRounded(b.center().Phi(), 2. * M_PI));
      }));
  Module* mod2 = const_cast<Module*>(mod);  // TO DO : Ugly, completely delete this ! actually, PtrSet should be defined and used as a modules container
  return mod2;
}


Module* Bundle::maxPhiModule() const {
  const Module* mod = &(*std::max_element(modules_.begin(), modules_.end(), [](const Module& a, const Module& b) {
	return (femodRounded(a.center().Phi(), 2. * M_PI) <= femodRounded(b.center().Phi(), 2. * M_PI));
      }));
  Module* mod2 = const_cast<Module*>(mod);  // TO DO : Ugly, completely delete this ! actually, PtrSet should be defined and used as a modules container
  return mod2;
}


const int Bundle::computePlotColor(const int id, const bool isPositiveCablingSide) const {
  int plotColor = 0;
  int plotId = (isPositiveCablingSide ? id : (id - 20000));
  int plotType = 2 + plotId % 2;  // Barrel : Identifies Flat vs Tilted. Endcap : Identifies PS10GA vs PG10GB vs PS5G vs 2S type.
  int dizaine = plotId / 10;
  int plotPhi = dizaine % 3;  // Barrel : Identifies phiSegmentRef. Endcap : Identifies phiRegionRef.
  plotColor = plotType * 3 + plotPhi;
  return plotColor;
}


const ChannelSection* Bundle::opticalChannelSection() const {
  const ChannelSection* opticalSection = nullptr;
  if (cable_) opticalSection = cable_->opticalChannelSection();
  if (!cable_ || !opticalSection) throw PathfulException("cable_ or opticalSection is nullptr");
  return opticalSection;
}


/* Untilted TBPS only: Id of the bundle located on the same Phi, but connected to the tilted modules.
 */
const int Bundle::tiltedBundleId() const {
  if (!isBarrelPSFlatPart()) logERROR("Tried to access tiltedBundleId, but not in TBPS flat part.");
  const int bundleId = myid();
  const int tiltedBundleId = bundleId - femod(bundleId, 10);
  return tiltedBundleId;
}


/* Id of the bundle located on the other cabling side, by rotation of 180° around CMS_Y.
 */
const int Bundle::stereoBundleId() const {
  if (!isBarrel()) logERROR("Tried to access stereoBundleId in TEDD, where it is not implemented (only implemented for Barrel).");
  return stereoBundleId_;
}


void Bundle::setIsPowerRoutedToBarrelLowerSemiNonant(const bool isLower) {
  if (!isBarrel()) logERROR("Tried to set Barrel power routing on a Bundle located in TEDD.");
  isPowerRoutedToBarrelLowerSemiNonant_ = isLower;
}


const bool Bundle::isPowerRoutedToBarrelLowerSemiNonant() const {
  if (!isBarrel()) logERROR("Tried to access Barrel power routing from a Bundle located in TEDD.");
  return isPowerRoutedToBarrelLowerSemiNonant_;
}

#include "Cabling/Bundle.hh"
#include "Cabling/Cable.hh"


Bundle::Bundle(const int id, const Category& type, const std::string subDetectorName, const int layerDiskNumber, const PhiPosition& phiPosition, const bool isPositiveCablingSide, const bool isTiltedPart) :
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
}


void Bundle::moveMaxPhiModuleFromOtherBundle(Bundle* otherBundle) {
  Container& otherBundleModules = otherBundle->modules();
  auto maxPhiModuleIt = std::max_element(otherBundleModules.begin(), otherBundleModules.end(), [](const Module& a, const Module& b) {
      return (femod(CoordinateOperations::stereoPhi(a.center()), 2. * M_PI) <= femod(CoordinateOperations::stereoPhi(b.center()), 2. * M_PI));
    });

  modules_.transfer(modules_.end(),
		    maxPhiModuleIt,
		    otherBundleModules);
}

void Bundle::moveMinPhiModuleFromOtherBundle(Bundle* otherBundle) {
  Container& otherBundleModules = otherBundle->modules();
  auto minPhiModuleIt = std::min_element(otherBundleModules.begin(), otherBundleModules.end(), [](const Module& a, const Module& b) {
      return (femod(CoordinateOperations::stereoPhi(a.center()), 2. * M_PI) <= femod(CoordinateOperations::stereoPhi(b.center()), 2. * M_PI));
    });

  modules_.transfer(modules_.end(), 
		    minPhiModuleIt,
		    otherBundleModules);
}

const double Bundle::minPhi() const { 
  double min = std::numeric_limits<double>::max();
  for (const auto& m : modules_) { min = MIN(min, femod(CoordinateOperations::stereoPhi(m.center()), 2. * M_PI) ); } return min;
}

const double Bundle::maxPhi() const { 
  double max = 0.;
  for (const auto& m : modules_) { max = MAX(max, femod(CoordinateOperations::stereoPhi(m.center()), 2. * M_PI) ); } return max;
}

Module* Bundle::minPhiModule() const {
  const Module* mod = &(*std::min_element(modules_.begin(), modules_.end(), [](const Module& a, const Module& b) {
	return (femod(CoordinateOperations::stereoPhi(a.center()), 2. * M_PI) <= femod(CoordinateOperations::stereoPhi(b.center()), 2. * M_PI));
      }));
  Module* mod2 = const_cast<Module*>(mod);  // TO DO : Ugly, completely delete this ! actually, PtrSet should be defined and used as a modules container
  return mod2;
}

Module* Bundle::maxPhiModule() const {
  const Module* mod = &(*std::max_element(modules_.begin(), modules_.end(), [](const Module& a, const Module& b) {
	return (femod(CoordinateOperations::stereoPhi(a.center()), 2. * M_PI) <= femod(CoordinateOperations::stereoPhi(b.center()), 2. * M_PI));
      }));
  Module* mod2 = const_cast<Module*>(mod);  // TO DO : Ugly, completely delete this ! actually, PtrSet should be defined and used as a modules container
  return mod2;
}


const int Bundle::computePlotColor(const int id, const bool isPositiveCablingSide) const {
  int plotColor = 0;
  int plotId = (isPositiveCablingSide ? id : (id - 20000));
  int plotType = 2 + plotId % 2;  // Barrel : Identifies Flat vs Tilted. Endcap : Identifies PS10G vs PG5GA vs PS5GB vs 2S type.
  int dizaine = plotId / 10;
  int plotPhi = dizaine % 3;  // Barrel : Identifies phiSegmentRef. Endcap : Identifies phiRegionRef.
  plotColor = plotType * 3 + plotPhi;
  return plotColor;
}

#include "ITCabling/PowerChain.hh"
#include "Cabling/HvLine.hh"


PowerChain::PowerChain(const int powerChainId, const bool isPositiveZEnd, const bool isPositiveXSide, const std::string subDetectorName, const int layerDiskNumber, const int phiRef, const int ringNumber, const bool isRingInnerEnd) :
  isPositiveZEnd_(isPositiveZEnd),
  isPositiveXSide_(isPositiveXSide),
  subDetectorName_(subDetectorName),
  layerDiskNumber_(layerDiskNumber),
  phiRef_(phiRef),
  ringNumber_(ringNumber),
  isRingInnerEnd_(isRingInnerEnd)
{
  myid(powerChainId);
  plotColor_ = computePlotColor(id, isPositiveZEnd);
};


PowerChain::~PowerChain() {
  delete hvLine_;    // TO DO: switch to smart pointers and remove this!
  hvLine_ = nullptr;
}


/*
const double PowerChain::minPhi() const { 
  double min = std::numeric_limits<double>::max();
  for (const auto& m : modules_) { min = MIN(min, femodRounded(m.center().Phi(), 2. * M_PI) ); } return min;
}


const double PowerChain::maxPhi() const { 
  double max = 0.;
  for (const auto& m : modules_) { max = MAX(max, femodRounded(m.center().Phi(), 2. * M_PI) ); } return max;
}

const double PowerChain::meanPhi() const {
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
*/


const int PowerChain::computePlotColor(const int phiRef, const int quarterRingIndex) const {
  int plotColor = 0;

  const int plotPhi = phiRef % 2;

  int plotId = (isPositiveCablingSide ? id : (id - 20000));
  int plotType = 2 + plotId % 2;  // Barrel : Identifies Flat vs Tilted. Endcap : Identifies PS10GA vs PG10GB vs PS5G vs 2S type.
  int dizaine = plotId / 10;
  int plotPhi = dizaine % 3;  // Barrel : Identifies phiSegmentRef. Endcap : Identifies phiRegionRef.
  plotColor = plotType * 3 + plotPhi;
  return plotColor;
}


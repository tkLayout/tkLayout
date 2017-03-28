#include "WeightDistributionGrid.hh"
#include "MaterialObject.hh"
#include <algorithm>

namespace material {

  WeightDistributionGrid::WeightDistributionGrid(double binDimension) :
    binDimension_(binDimension) {}

  void WeightDistributionGrid::addTotalGrams(double minZ, double minR, double maxZ, double maxR, double length, double surface, const MaterialObject& materialObject) {
    int binIndexZStart = floor(minZ / binDimension_);
    int binIndexZ = binIndexZStart;
    int binIndexR = floor(minR / binDimension_);
    double binMinZ, binMinR, binMaxZ, binMaxR, zMul, rMul;

    do {
      binMinR = floor(binIndexR * binDimension_);
      binMaxR = binMinR + binDimension_;

      do {
        binMinZ = floor(binIndexZ * binDimension_);
        binMaxZ = binMinZ + binDimension_;

        zMul = std::max(binMinZ, minZ) - std::min(binMaxZ, maxZ);
        rMul = std::max(binMinR, minR) - std::min(binMaxR, maxR);

        operator[](std::make_pair(binIndexZ, binIndexR)) += 
          (materialObject.totalGrams(length, surface) * zMul * rMul) / ((maxZ - minZ) * (maxR - minR));

        binIndexZ ++;
      } while (binMaxZ < maxZ);
      
      binIndexR ++;
      binIndexZ = binIndexZStart;
    } while (binMinR < maxR);
  }

  double WeightDistributionGrid::binDimension() const {
    return binDimension_;
  }
}


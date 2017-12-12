#ifndef WEIGHTDISTRIBUTIONGRID_H
#define WEIGHTDISTRIBUTIONGRID_H

#include <map>

namespace material {

  class MaterialObject;
  
  typedef std::map<std::pair<int, int>,double> WeightDistributionGridMapType;

  class WeightDistributionGrid : public WeightDistributionGridMapType {
  public:
    WeightDistributionGrid(double binDimension);
    virtual ~WeightDistributionGrid() {};
    void addTotalGrams(double minZ, double minR, double maxZ, double maxR, double length, double surface, const MaterialObject& materialObject);
    double binDimension() const;
  private:
    double binDimension_;
  };
}

#endif // WEIGHTDISTRIBUTIONGRID_H

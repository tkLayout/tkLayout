#ifndef FLUKAGEOMETRYANALYZER_H
#define FLUKAGEOMETRYANALYZER_H

#include <string>
#include <map>
#include <vector>

#include "Tracker.hh"
#include "MaterialBudget.hh"
#include "InactiveElement.hh"

namespace MaterialBillAnalyzerData {

  typedef std::map<std::string, double> MaterialMap;
  using namespace insur;

  class ServiceElement {
  public: 
    double zmin, zmax, rmin, rmax;
    MaterialMap materialMap;
  };
}

using namespace MaterialBillAnalyzerData;

class MaterialBillAnalyzer {
 private:
  typedef std::map<std::string, MaterialMap> LayerMaterialMap;
  typedef std::vector<ServiceElement> ServicesMaterialVector;
  ServicesMaterialVector servicesMaterialVector_;
  LayerMaterialMap layerMaterialMap_;
  void inspectInactiveElements(const std::vector<InactiveElement>& inactiveElements);
  void inspectModules(std::vector<std::vector<insur::ModuleCap> >& tracker);


 public:
  std::string outputTable;
  void inspectTracker(MaterialBudget&);

};


#endif

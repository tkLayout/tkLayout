#include "AnalyzerVisitors/MaterialBillAnalyzer.hh"

#include "ModuleCap.hh"
#include <iostream>

void MaterialBillAnalyzer::inspectInactiveElements(const std::vector<InactiveElement>& inactiveElements) {
  for (const auto& it : inactiveElements) {
    const std::map<std::string, double>& localMasses = it.getLocalMasses();
    for (auto massIt : localMasses) {
      outputTable += any2str(it.getInnerRadius()) + ", ";
      outputTable += any2str(it.getInnerRadius()+it.getRWidth()) + ", ";
      outputTable += any2str(it.getZOffset()) + ", ";
      outputTable += any2str(it.getZOffset()+it.getZLength()) + ", ";
      outputTable += massIt.first + ", ";
      outputTable += any2str(massIt.second) + "\n";
    }
  }
}

void MaterialBillAnalyzer::inspectModules(std::vector<std::vector<insur::ModuleCap> >& tracker) {
  // loop over layers
  for (auto layerIt : tracker ) {
    // Loop over modules
    for (auto moduleIt : layerIt ) {
      auto myModuleCap = (moduleIt);
      auto myModule = &(myModuleCap.getModule());
      // TODO: put this in a better place
      // (and make a better module typing)
      struct Visitor : public ConstGeometryVisitor {
        std::string id_;
        void visit(const BarrelModule& m) { id_ = m.cntName() + "_L" + any2str(m.layer()); }
        void visit(const EndcapModule& m) { id_ = m.cntName() + "_D" + any2str(m.disk()); }
      };
      Visitor v;
      myModule->accept(v);
      MaterialMap& layerMaterial = layerMaterialMap_[v.id_];
      const std::map<std::string, double>& localMasses = myModuleCap.getLocalMasses();
      for (const auto &it : localMasses)  layerMaterial[it.first]+=it.second;
    }
  }
}

void MaterialBillAnalyzer::inspectTracker(MaterialBudget& mb) {
  outputTable="";

  inspectModules(mb.getBarrelModuleCaps());
  inspectModules(mb.getEndcapModuleCaps());

  outputTable += "material in layers\n";
  outputTable += "layer, material, weight_grams\n";

  for (const auto &it : layerMaterialMap_) {
    for (const auto &layMats : it.second) {
      outputTable += it.first + ", ";
      outputTable += layMats.first+ ", " + any2str(layMats.second) + "\n";
    }
  }

  InactiveSurfaces& is = mb.getInactiveSurfaces();
  std::vector<InactiveElement>& barrelServices = is.getBarrelServices();
  std::vector<InactiveElement>& endcapServices = is.getEndcapServices();
  std::vector<InactiveElement>& supports = is.getSupports();
  outputTable += "other elements\n";
  outputTable += "r_in, r_out, z_in, z_out, material, weight_grams\n";
  inspectInactiveElements(barrelServices);
  inspectInactiveElements(endcapServices);
  inspectInactiveElements(supports);
}


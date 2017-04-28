#ifndef RIBBON_HH
#define RIBBON_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Module.hh"

//#include "CablingVisitable.h"
//#include "CablingVisitor.h"

/*using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;*/

//class Ribbon : public PropertyObject, public Buildable, public Identifiable<int>, public CablingVisitable {
class Ribbon : public PropertyObject, public Identifiable<int> {
  std::string type_;

  //typedef PtrVector<Module> Container;
  typedef std::vector<Module*> Container;
  Container modules_;

  //Property<int, Default> nModulesPerRibbon;

public:
  Ribbon(int id, std::string type) {
    myid(id);
    type_ = type;
  };

  const std::string getType() const { return type_; }

  void addModule(Module& m) { modules_.push_back(&m); }
  const Container& modules() const { return modules_; }

  //void removeModule(int detId) { modules_.erase_if([&](Module* m) { return (m->myDetId() == detId); } ); }

  int numModules() const { return modules_.size(); }

  const double minPhi() const { 
    double min = std::numeric_limits<double>::max();
    for (const auto& m : modules_) { min = MIN(min, femod(m->center().Phi(), 2. * M_PI) ); } return min;
  }

  const double maxPhi() const { 
    double max = 0.;
    for (const auto& m : modules_) { max = MAX(max, femod(m->center().Phi(), 2. * M_PI) ); } return max;
  }

  const Module* minPhiModule() const {
    return *std::min_element(modules_.begin(), modules_.end(), [&](Module* a, Module* b) {
	return (femod(a->center().Phi(), 2. * M_PI) < femod(b->center().Phi(), 2. * M_PI));
      });
  }

  const Module* maxPhiModule() const {
    return *std::max_element(modules_.begin(), modules_.end(), [&](Module* a, Module* b) {
	return (femod(a->center().Phi(), 2. * M_PI) < femod(b->center().Phi(), 2. * M_PI));
      });
  }

  


  /*Ribbon() :
            nModulesPerRibbon      ("nModulesPerRibbon"      , parsedAndChecked(), 6)
  {}
  void setup() {
  }

  Container& modules() { return detectormodules_; }
  const Container& modules() const { return detectormodules_; }
  int nModules() const { return detectormodules_.size(); }
  int maxModules() {return nModulesPerRibbon(); }
  
  void check() override;
  void build();

  void addModule(Module& m) {}

  void accept(CablingVisitor& v) { 
    v.visit(*this); 
    for (Module& m : detectormodules_) { m.accept(v); }
  }
  void accept(ConstCablingVisitor& v) const { 
    v.visit(*this); 
    for (const auto& m : detectormodules_) { m.accept(v); }
    }*/

};



#endif

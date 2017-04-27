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

  std::string getType() { return type_; }

  void addModule(Module& m) { modules_.push_back(&m); }
  const Container& modules() const { return modules_; }






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

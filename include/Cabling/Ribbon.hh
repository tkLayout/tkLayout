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
class Ribbon : public PropertyObject, public Buildable, public Identifiable<int> {
  std::string type_;
  std::string subDetectorName_;
  int layerDiskNumber_;

  double startPhi_;
  double phiRegionWidth_;
  int phiRegionRef_;
  double phiSectorWidth_;
  int phiSectorRef_;


  typedef PtrVector<Module> Container;
  //typedef PtrSet<Module> Container;
  //typedef std::vector<Module*> Container;
  Container modules_;

  //Property<int, Default> nModulesPerRibbon;

public:

  Ribbon(int id, std::string type, std::string subDetectorName, int layerDiskNumber, double startPhi, double phiRegionWidth, int phiRegionRef, const double phiSectorWidth, int phiSectorRef) {
    myid(id);

    type_ = type;
    subDetectorName_ = subDetectorName;
    layerDiskNumber_ = layerDiskNumber;

    startPhi_ = startPhi;
    phiRegionWidth_ = phiRegionWidth;
    phiRegionRef_ = phiRegionRef;
    phiSectorWidth_ = phiSectorWidth;
    phiSectorRef_ = phiSectorRef;


  };

  const std::string type() const { return type_; }
  const std::string subDetectorName() const { return subDetectorName_; }
  const int layerDiskNumber() const { return layerDiskNumber_; }

  const double startPhi() const { return startPhi_; }
  const double phiRegionWidth() const { return phiRegionWidth_; }
  const int phiRegionRef() const { return phiRegionRef_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }



  void addModule(Module* m) { modules_.push_back(m); }
  const Container& modules() const { return modules_; }

  void removeModule(Module* m) {
    int detId = m->myDetId();
    modules_.erase_if([detId](Module& m) { return (m.myDetId() == detId); });
  }

  int numModules() const { return modules_.size(); }

  const double minPhi() const { 
    double min = std::numeric_limits<double>::max();
    for (const auto& m : modules_) { min = MIN(min, femod(m.center().Phi(), 2. * M_PI) ); } return min;
  }

  const double maxPhi() const { 
    double max = 0.;
    for (const auto& m : modules_) { max = MAX(max, femod(m.center().Phi(), 2. * M_PI) ); } return max;
  }

  Module* minPhiModule() const {
    const Module* mod = &(*std::min_element(modules_.begin(), modules_.end(), [](const Module& a, const Module& b) {
	return (femod(a.center().Phi(), 2. * M_PI) <= femod(b.center().Phi(), 2. * M_PI));
	}));
    Module* mod2 = const_cast<Module*>(mod);  // TO DO : Ugly, completely delete this ! actually, PtrSet should be defined and used as a modules container
    return mod2;
  }

  Module* maxPhiModule() const {
    const Module* mod = &(*std::max_element(modules_.begin(), modules_.end(), [](const Module& a, const Module& b) {
	  return (femod(a.center().Phi(), 2. * M_PI) <= femod(b.center().Phi(), 2. * M_PI));
	}));
    Module* mod2 = const_cast<Module*>(mod);  // TO DO : Ugly, completely delete this ! actually, PtrSet should be defined and used as a modules container
    return mod2;
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

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
  std::string subDetectorName_;
  int layerDiskNumber_;

  double startPhi_;
  double phiRegionWidth_;
  int phiRegionRef_;
  double phiSectorWidth_;
  int phiSectorRef_;


  //typedef PtrVector<Module> Container;
  typedef std::vector<Module&> Container;
  Container modules_;

  //Property<int, Default> nModulesPerRibbon;

public:
  /*Property<std::string, AutoDefault> type;
  Property<std::string, AutoDefault> subDetectorName;
  Property<int, AutoDefault> layerDiskNumber;

  Property<double, AutoDefault> startPhi;
  Property<double, AutoDefault> phiRegionWidth;
  Property<int, AutoDefault> phiRegionRef;
  Property<double, AutoDefault> phiSectorWidth;
  Property<int, AutoDefault> phiSectorRef;*/

  Ribbon(int id, std::string type, std::string subDetectorName, int layerDiskNumber, double startPhi, double phiRegionWidth, int phiRegionRef, const double phiSectorWidth, int phiSectorRef) {
    myid(id);
    /*type_ = type;

    startPhi_ = startPhi;
    phiRegionWidth_ = phiRegionWidth;
    phiRegionRef_ = phiRegionRef;
    const double phiSectorWidth_ = phiSectorWidth;
    phiSectorRef_ = phiSectorRef;*/


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



  //const std::string getType() const { return type_; }

  void addModule(Module& m) { modules_.push_back(m); }
  const Container& modules() const { return modules_; }

  void removeModule(Module& m) {
    int detId = m.myDetId();
    modules_.erase(std::remove_if(modules_.begin(), modules_.end(), [&](Module& a) { return (a.myDetId() == detId); } ), modules_.end()); 
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

  Module& minPhiModule() const {
    return *std::min_element(modules_.begin(), modules_.end(), [&](Module& a, Module& b) {
	return (femod(a.center().Phi(), 2. * M_PI) <= femod(b.center().Phi(), 2. * M_PI));
      });
  }

  Module& maxPhiModule() const {
    return *std::max_element(modules_.begin(), modules_.end(), [&](Module& a, Module& b) {
	return (femod(a.center().Phi(), 2. * M_PI) <= femod(b.center().Phi(), 2. * M_PI));
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

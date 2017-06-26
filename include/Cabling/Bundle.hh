#ifndef BUNDLE_HH
#define BUNDLE_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Module.hh"
#include "Cabling/PhiPosition.hh"

//#include "CoordinateOperations.hh"
//#include "CablingVisitable.h"
//#include "CablingVisitor.h"
/*using std::string;
using std::vector;
using std::pair;
using std::unique_ptr;*/

namespace insur { class Cable; }
using insur::Cable;


//class Bundle : public PropertyObject, public Buildable, public Identifiable<int>, public CablingVisitable {
class Bundle : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<Module> Container; 

public:
  Bundle(const int id, const Category& type, const std::string subDetectorName, const int layerDiskNumber, const PhiPosition& phiPosition, const bool isPositiveCablingSide, const bool isTiltedPart);
  ~Bundle();

  // MODULES CONNECTED TO THE BUNDLE.
  const Container& modules() const { return modules_; }
  Container& modules() { return modules_; }
  const int numModules() const { return modules_.size(); }
  void addModule(Module* m) { modules_.push_back(m); }

  // CABLE THE BUNDLE IS CONNECTED TO.
  const Cable* getCable() const { return cable_; }
  void setCable(Cable* cable) { cable_ = cable; }
  
  /*void removeModule(Module* m) {
    int detId = m->myDetId();
    modules_.erase_if([detId](Module& m) { return (m.myDetId() == detId); });
    }*/

  void moveMaxPhiModuleFromOtherBundle(Bundle* otherBundle);
  void moveMinPhiModuleFromOtherBundle(Bundle* otherBundle);

  const double minPhi() const;
  const double maxPhi() const;

  Module* minPhiModule() const;
  Module* maxPhiModule() const;

  const Category& type() const { return type_; }
  const std::string subDetectorName() const { return subDetectorName_; }
  const int layerDiskNumber() const { return layerDiskNumber_; }
  const PhiPosition& phiPosition() const { return phiPosition_; }  
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }
  const bool isTiltedPart() const { return isTiltedPart_; }

  const int plotColor() const { return plotColor_; }


  /*Bundle() :
            nModulesPerBundle      ("nModulesPerBundle"      , parsedAndChecked(), 6)
  {}
  void setup() {
  }

  Container& modules() { return detectormodules_; }
  const Container& modules() const { return detectormodules_; }
  int nModules() const { return detectormodules_.size(); }
  int maxModules() {return nModulesPerBundle(); }
  
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


private:
  const int computePlotColor(const int id, const bool isPositiveCablingSide) const;

  //typedef PtrSet<Module> Container;
  //Property<int, Default> nModulesPerBundle;
  Container modules_;

  Cable* cable_ = nullptr;

  Category type_;
  std::string subDetectorName_;
  int layerDiskNumber_;
  PhiPosition phiPosition_;
  bool isPositiveCablingSide_;
  bool isTiltedPart_;

  int plotColor_;
};



#endif

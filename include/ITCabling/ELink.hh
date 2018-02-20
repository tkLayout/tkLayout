#ifndef ELINK_HH
#define ELINK_HH

#include <vector>
#include <string>

#include "Property.hh"
//#include "Module.hh"
#include "ITCabling/inner_cabling_functions.hh"


//namespace insur { class DetectorModule; }
//using insur::Module;

//namespace insur { class HvLine; }
//using insur::HvLine;


class ELink : public PropertyObject, public Buildable, public Identifiable<int> {

public:
  ELink(const std::string eLinkId);
  //~ELink();

  // MODULE TO WHICH TO ELINK IS CONNECTED
  /*const Module* getModule() const {
    if (!module_) throw PathfulException("module_ is nullptr");
    return module_;
    }*/

  const std::string getName() const { return name_; }

private:
  //Module* module_ = nullptr;

  std::string name_;
};



#endif

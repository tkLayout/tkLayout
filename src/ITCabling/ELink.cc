#include "ITCabling/ELink.hh"
//#include "ITCabling/HvLine.hh"
//#include "Module.hh"


ELink::ELink(const std::string eLinkId) :
  //module_(m),
  name_(eLinkId) 
{};


/*ELink::~ELink() {
  delete module_;    // TO DO: switch to smart pointers and remove this!
  module_ = nullptr;
  }*/


/*void ELink::addModule(Module* m) { 
  modules_.push_back(m);
  }*/

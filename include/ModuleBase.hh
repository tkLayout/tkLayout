#ifndef MODULE_BASE_H
#define MODULE_BASE_H


#include "Visitor.hh"
#include "Property.hh"
#include "Visitable.hh"


//typedef uint32_t ModuleType;


class ModuleBase : public PropertyObject, public Buildable, public Placeable, public Identifiable<int>, public Visitable {
public:
  virtual void build() = 0;
  //virtual ModuleBase* clone() = 0;

};

class ModuleDecorable : public ModuleBase, public Decorable {};

#endif

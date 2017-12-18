#ifndef SUPPORT_H
#define SUPPORT_H

#include "Property.hh"

class Support : public PropertyObject, public Identifiable<int>, public Buildable {
public:
  int index() const { return myid(); }
  Property<double, NoDefault> midZ;

  Support() : midZ("midZ", parsedAndChecked()) {}
  void build() {
    try { 
      check();
    } catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
    builtok(true);
  }

  std::pair<int, double> toPair() const { return { myid(), midZ() }; }
};


#endif

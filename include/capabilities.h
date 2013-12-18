#ifndef CAPABILITIES_H
#define CAPABILITIES_H

#include <typeinfo>
#include <string>

using std::string;


class Buildable {
protected:
  bool built_;
public:
  Buildable() : built_(false) {}
  bool builtok() const { return built_; }
  void builtok(bool state) { built_ = state; }
};

class Placeable {
protected:
  bool placed_;
public:
  Placeable() : placed_(false) {}
  bool placed() const { return placed_; }
  void placed(bool state) { placed_ = state; }
};
/*
typedef string IdentifiableType;

template<class T>
class Identifiable {
  const string base_;
  T::IdType myid_;
public:
  Identifiable() : base_(typeid(T).name()), myid_("NOID") {}
  //template<class U> void myid(U id) { myid_ = any2str(id); }
  void myid(T::IdType id) { myid_ = id; }
  T::IdType myid() const { return myid_; }
  IdentifiableType fullid() const { return base_ + "(" + any2str(myid()) + ")"; }
};
*/
template<class T>
class Identifiable {
  T myid_;
public:
  void myid(const T& id) { myid_ = id; }
  const T& myid() const { return myid_; }
};

template<class T> string fullid(const T& o) { return string(typeid(T).name()) + "(" + any2str(o.myid()) + ")"; }


#endif

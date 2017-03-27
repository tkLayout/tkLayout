#ifndef CAPABILITIES_H
#define CAPABILITIES_H

#include <typeinfo>
#include <string>
#include <bitset>

#include "GeometryFactory.hh"
#include "MessageLogger.hh"

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

class DetIdentifiable {
  uint32_t myDetId_ = 0;
  std::map<int, uint32_t> detIdRef_;

 public:
  uint32_t myDetId() const { return myDetId_; }
  std::bitset<32> myBinaryDetId() const { std::bitset<32> myBinaryDetId(myDetId_); return myBinaryDetId; }
  std::map<int, uint32_t> detIdRef() const { return detIdRef_; }
  // i will create a capabilities.cc and put that in it, sorry in the meantime ;p
  void buildDetId(std::map<int, uint32_t> refs, std::vector<int> schemeShifts) {
    detIdRef_ = refs;
    for (int i = 0; i < schemeShifts.size(); i++) {
      uint32_t ref = refs.at(i);
      int shift = schemeShifts.at(i);
      if (ref <= pow(2, shift) ) {
	myDetId_ <<= shift;
	myDetId_ |= ref;
      }
      else logWARNING("buildDetId : At rank " + any2str(i) + ", ref number " + any2str(ref) + " is reached, while size allocated by the DetId scheme is 2^" + any2str(shift) + ".");
    }
  }
};

template<class T>
class Clonable {
  template<class U> static void conditionalSetup(U* t, typename std::enable_if<std::is_void<decltype(t->setup())>::value>::type* = 0) { t->setup(); } 
  static void conditionalSetup(...) {}
protected:
  Clonable(const Clonable&) {}
public:
  Clonable() {}
  T* clone() const {
    T* t = new T(static_cast<const T&>(*this));
    conditionalSetup(t);
    return t;
  }

  template<class ...U> static T* make(const U&... args) {
    T* t = new T(args...);
    conditionalSetup(t);
    return t;
  }
};



#endif

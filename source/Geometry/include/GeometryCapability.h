#ifndef INCLUDE_GEOMETRY_CAPABILITY_H_
#define INCLUDE_GEOMETRY_CAPABILITY_H_

#include <typeinfo>
#include <string>

#include "StringConverter.h"

/*
 * Set that object has unique identifier that built properly
 */
class Buildable
{

 public:
  Buildable() : m_built(false) {}
  bool    builtok() const     { return m_built; }
  void    builtok(bool state) { m_built = state; }

 protected:
  bool m_built;
}; // Class

/*
 * Set that object has unique identifier that placed in space
 */
class Placeable
{

 public:
  Placeable() : m_placed(false) {}
  bool placed() const     { return m_placed; }
  void placed(bool state) { m_placed = state; }

 protected:
  bool m_placed;
}; // Class

/*
 * Set that object has unique id
 */
template<class T> class Identifiable
{

 public:
  const T& myid() const      { return m_myid; }
  void     myid(const T& id) { m_myid = id; }


 private:
  T m_myid;
}; // Class

/*
 * General method returning object type & its unique id (template object T must be identifiable)
 */
template<class T> std::string fullid(const T& o) { return std::string(typeid(T).name()) + "(" + any2str(o.myid()) + ")"; }

/*
 * Set that object clonable (currently GeometryFactory used instead)
 */
template<class T> class Clonable
{

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

  template<class U> static void conditionalSetup(U* t, typename std::enable_if<std::is_void<decltype(t->setup())>::value>::type* = 0) { t->setup(); }
  static void conditionalSetup(...) {}

 protected:
  Clonable(const Clonable&) {}
}; // Class

#endif /* INCLUDE_GEOMETRY_CAPABILITY_H_ */

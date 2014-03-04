#ifndef PROPERTY_H
#define PROPERTY_H

#include <vector>
#include <map>
#include <string>
#include <list>
#include <set>
#include <memory>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <typeinfo>
#include <typeindex>

#include <boost/property_tree/ptree.hpp>

#include "global_funcs.h"
#include "Decorator.h"
#include "capabilities.h"
#include "StringSet.h"

using std::string;
using std::vector;
using std::map;
using std::list;
using std::shared_ptr;

using boost::property_tree::ptree;


class PathfulException : public std::invalid_argument {
  string path_;
public:
  PathfulException(const string& what, const string& path) : std::invalid_argument(what.c_str()) { pushPath(path); }
  template<class T> PathfulException(const string& what, const T& obj, const string& objid) : std::invalid_argument(what.c_str()) { pushPath(obj, objid); }
  PathfulException(const string& what) : std::invalid_argument(what.c_str()) {}
  void pushPath(const string& p) { path_ = p + (!path_.empty() ? "." + path_ : ""); }
  template<class T> void pushPath(const T& obj, const string& objid) { pushPath(string(typeid(obj).name()) + "(" + objid + ")"); }
  template<class T, class U> void pushPath(const T& obj, const U& objid) { pushPath(string(typeid(obj).name()) + "(" + any2str(objid) + ")"); }
  const string& path() const { return path_; }
  //virtual const char* what() const throw() override { return (pathString() + " : " + std::invalid_argument::what()).c_str(); }
};

struct CheckedPropertyMissing : public PathfulException { CheckedPropertyMissing(const string& objid) : PathfulException("Checked property not set", objid) {} };
struct RequestedUnsetValue : public PathfulException { RequestedUnsetValue() : PathfulException("Tried to get value from an unset property") {} };
struct InvalidPropertyValue : public PathfulException { InvalidPropertyValue(const string& objid, const string& val) : PathfulException("Value '" + val + "' is invalid for property", objid) {} };
struct InvalidComputable : public PathfulException { InvalidComputable() : PathfulException("Computable object is invalid or unset") {} };


template<class T>
struct Stateful { 
  typedef T State;
  virtual State state() const = 0; 
};

struct Parsable : public Stateful<bool> {
  virtual string name() const = 0;
  virtual void fromPtree(const ptree& pt) = 0;
  virtual void fromString(const string& s) = 0;
};

template<typename T>
struct PropertyBase : public Stateful<bool> {
  virtual void operator()(const T&) = 0;
  virtual const T& operator()() const = 0;
  virtual void clear() = 0;
};


template<typename T>
class Default : public PropertyBase<T> {
  T value_;
  const T default_;
public:
  Default(const T& value) : value_(value), default_(value) {}
  void operator()(const T& value) { value_ = value; }
  const T& operator()() const { return value_; }
  bool state() const { return true; }
  void clear() { value_ = default_; }
};

template<typename T>
struct AutoDefault : public Default<T> { AutoDefault() : Default<T>(T()) {} };

template<typename T>
class NoDefault : public PropertyBase<T> {
  T value_;
  bool state_;
public:
  NoDefault() : state_(false) {}
  void operator()(const T& value) { value_ = value; state_ = true; }
  const T& operator()() const { 
    if (!state_) throw RequestedUnsetValue(); 
    return value_; 
  }
  bool state() const { return state_; }
  void clear() { state_ = false; }
};

template<typename T>
class Computable : public PropertyBase<T> {
  typedef std::function<T()> Func;
  Func get;
  mutable T value_;
  mutable bool state_;
public:
  template<class U = Func> Computable(const U& getter /*= []()->T{ throw InvalidComputable(); }*/) : get(Func(getter)), state_(false) {} // default arg commented out due to bug in gcc 4.7.2-5
  Computable() : get([]()->T{ throw InvalidComputable(); }), state_(false) {}
  Computable(const Computable<T>& other) : Computable() { value_ = other.value_; state_ = other.state_; } // Func does not get copied!! needs to be manually setup() again to prevent issues with the captures (capture happens at point of declaration)
  void operator()(const T& value) { value_ = value; state_ = true; }
  const T& operator()() const {
    if (!state_) { value_ = get(); state_ = true; }
    return value_;
  }
  bool state() const { return state_; }
  void clear() { state_ = false; }
  template<class U = Func> void setup(const U& getter) { get = Func(getter); }
};

template<typename T>
class UncachedComputable : public PropertyBase<T> {
  typedef std::function<T()> Func;
  Func get;
  mutable T value_;
  mutable bool state_;
public:
  template<class U = Func> UncachedComputable(const U& getter /*= []()->T{ throw InvalidComputable(); }*/) : get(Func(getter)), state_(false) {}
  UncachedComputable() : get([]()->T{ throw InvalidComputable(); }), state_(false) {}
  UncachedComputable(const UncachedComputable<T>& other) : UncachedComputable() { value_ = other.value_; state_ = other.state_; } // Func does not get copied!! needs to be manually setup() again to prevent issues with the captures (capture happens at point of declaration)
  void operator()(const T& value) { value_ = value; state_ = true; }
  const T& operator()() const { return state_ ? value_ : value_ = get(); }
  bool state() const { return state_; }
  void clear() { state_ = false; }
  template<class U = Func> void setup(const U& getter) { get = Func(getter); }
};

template<typename T>
class Fallback : public PropertyBase<T> {
  typedef PropertyBase<T> Prop;
  const Prop& get;
  T value_;
  bool set_;
public:
  Fallback(const Prop& getter) : get(getter), set_(false) {}
  void operator()(const T& value) { value_ = value; set_ = true; }
  const T& operator()() const { return set_ ? value_ : get(); }
  bool state() const { return get.state(); } 
  void clear() { set_ = false; }
};



typedef std::map<string, Parsable*> PropertyMap;


template<typename T, template<typename> class ValueHolder>
class Property : public PropertyBase<T>, public Parsable {
  ValueHolder<T> valueHolder_;
  const string& name_;
public:
  Property(const string& name, PropertyMap& registrar, const ValueHolder<T>& valueHolder = ValueHolder<T>()) : valueHolder_(valueHolder), name_(StringSet::ref(name)) { registrar[name] = this; }
  Property(const string& name, const ValueHolder<T>& valueHolder = ValueHolder<T>()) : valueHolder_(valueHolder), name_(StringSet::ref(name)) {}
  Property(const ValueHolder<T>& valueHolder = ValueHolder<T>()) : valueHolder_(valueHolder), name_(StringSet::ref("unnamed")) {}
  void operator()(const T& value) { valueHolder_(value); }
  const T& operator()() const { 
    try { return valueHolder_(); } 
    catch(PathfulException& pe) { pe.pushPath(name_); throw; }
  }
  bool state() const { return valueHolder_.state(); }
  void clear() { valueHolder_.clear(); }
  template<class ...U> void setup(const U&... valueHolderArgs) { valueHolder_.setup(valueHolderArgs...); } 
  string name() const { return name_; }
  void fromPtree(const ptree& pt) { valueHolder_(str2any<T>(pt.data())); }
  void fromString(const string& s) { valueHolder_(str2any<T>(s)); }
};

template<typename T, template<typename> class ValueHolder>
class ReadonlyProperty : public Property<T, ValueHolder> {
  void operator()(const T&) {}
public:
  ReadonlyProperty(const string& name, PropertyMap& registrar, const ValueHolder<T>& valueHolder = ValueHolder<T>()) : Property<T, ValueHolder>(name, registrar, valueHolder) {}
  ReadonlyProperty(const string& name, const ValueHolder<T>& valueHolder = ValueHolder<T>()) : Property<T, ValueHolder>(name, valueHolder) {}
  ReadonlyProperty(const ValueHolder<T>& valueHolder = ValueHolder<T>()) : Property<T, ValueHolder>(valueHolder) {}
  ///using Property<T, ValueHolder>::Property; // constructor inheritance not supported by GCC as of version 4.7.2
  using Property<T, ValueHolder>::operator();
  
#define ALLOW_FORCE_SET
#ifdef ALLOW_FORCE_SET
  void force(const T& value) { Property<T, ValueHolder>::operator()(value); } // for testing purposes only
#endif
};


template<typename T, const char Sep = ','>
class PropertyVector : public Parsable {
  std::vector<T> values_;
  const string& name_;
public:
  PropertyVector(const string& name, PropertyMap& registrar, const std::initializer_list<T>& values = {}) : values_(values), name_(StringSet::ref(name)) { registrar[name] = this; }
  PropertyVector(const string& name, const std::initializer_list<T>& values = {}) : values_(values), name_(StringSet::ref(name)) {}
  PropertyVector(const std::initializer_list<T>& values = {}) : values_(values), name_(StringSet::ref("unnamed")) {}
  void operator()(size_t i, const T& value) { values_[i] = value; }
  const T& operator()(size_t i) const { return values_[i]; }
  typename std::vector<T>::const_iterator begin() const { return values_.begin(); }
  typename std::vector<T>::const_iterator end() const { return values_.end(); }
  bool state() const { return !values_.empty(); } 
  void clear() { values_.clear(); }
  string name() const { return name_; }
  void fromPtree(const ptree& pt) { fromString(pt.data()); }
  void fromString(const string& s) { 
    std::vector<T> values = split<T>(trim(s), string(1, Sep), true);
    if (!values[0].empty()) values_.clear();
    else values.erase(values.begin(), values.begin()+1);
    for (const auto& v : values) values_.push_back(v);
  }
  void appendString(const string& s) {
    values_.push_back(trim(s));
  }
}; 


template<typename T>
class PropertyNode : public Parsable, public map<T, ptree> { 
  const string& name_;
public:
  PropertyNode(const string& name, PropertyMap& registrar) : name_(StringSet::ref(name)) { registrar[name] = this; }
  PropertyNode(const string& name) : name_(StringSet::ref(name)) {}
  bool state() const { return !this->empty(); }
  void clear() { map<T, ptree>::clear(); }
  string name() const { return name_; }
  void fromPtree(const ptree& pt) { 
    T key = str2any<T>(pt.data());
    if (this->count(key) == 0) this->insert(make_pair(key, pt)); 
    else for (auto& tel : pt) this->at(key).add_child(tel.first, tel.second);
  }
  void fromString(const string& s) { this->insert(make_pair(str2any<T>(s), ptree())); }
};

template<>
class PropertyNode<int> : public Parsable, public map<int, ptree> {
  const string& name_;
public:
  PropertyNode(const string& name, PropertyMap& registrar) : name_(StringSet::ref(name)) { registrar[name] = this; }
  PropertyNode(const string& name) : name_(StringSet::ref(name)) {}
  bool state() const { return !this->empty(); }
  void clear() { map<int, ptree>::clear(); }
  string name() const { return name_; }
  void fromPtree(const ptree& pt) { 
    vector<int> keys;
    auto tokens = split(pt.data(), ","); // split sequences like 1, 2, 3, 4
    for (auto t : tokens) {
      auto interval = split<int>(t, "-"); // split intervals like 1-5
      if (interval.size() == 1) keys.push_back(interval[0]); // actually not an interval
      else for (int i = interval[0]; i <= interval[1]; i++) keys.push_back(i); // explode the interval
    }
    for (auto k : keys) {
      if (this->count(k) == 0) this->insert(make_pair(k, pt)); 
      else for (auto& tel : pt) this->at(k).add_child(tel.first, tel.second);
    }
  }
  void fromString(const string& s) { this->insert(make_pair(str2any<int>(s), ptree())); }
};


typedef ptree PropertyTree;


class PropertyObject {
  PropertyMap parsedCheckedProperties_, checkedProperties_, parsedProperties_;
  PropertyTree pt_;

  void processProperties(PropertyMap& props) {
    for (auto& propElem : props) {
      auto childRange = pt_.equal_range(propElem.first);
      std::for_each(childRange.first, childRange.second, [&propElem](const ptree::value_type& treeElem) {
        propElem.second->fromPtree(treeElem.second); // takes care of duplicate entries (by overwriting the property value as many times as there are entries with the same key) and of node entries (in that case the PropertyNodes differentiates based on the value)
      });
      pt_.erase(propElem.first);
    }
  }
  void printAll(const PropertyTree& pt) {
    for (auto& p : pt) {
      std::cout << p.first << " = " << p.second.data();
      std::cout << " ... " << p.second.size() << " children";
      std::cout << std::endl;
    }
  }
protected:
  const PropertyTree& propertyTree() const { return pt_; }
  PropertyMap& parsedAndChecked() { return parsedCheckedProperties_; }
  PropertyMap& parsedOnly() { return parsedProperties_; }
  PropertyMap& checkedOnly() { return checkedProperties_; }

public:
  PropertyObject() {}
  virtual void store(const PropertyTree& newpt) {
    if (pt_.empty()) pt_ = newpt;
    else { 
      pt_.data() = newpt.data();
      for (auto& propElem : newpt) { 
        pt_.add_child(propElem.first, propElem.second); 
      } // merging trees in a careless manner, appending children without ever checking if an entry with the same key is already present. the duplicates thus formed will be all grabbed at parsing time by the properties (each duplicate entry overwrites the previous)
    }
//    std::cout << "============ " << pt_.data() << " ===========" << std::endl;
//    printAll(pt_);
    processProperties(parsedCheckedProperties_);
    processProperties(parsedProperties_);
  }
  virtual void check() {
    for (auto& v : parsedCheckedProperties_) {
      if (!v.second->state()) throw CheckedPropertyMissing(v.first); 
    }
    for (auto& v : checkedProperties_) {
      if (!v.second->state()) throw CheckedPropertyMissing(v.first); 
    }
  }

  virtual void cleanup() { pt_.clear(); parsedCheckedProperties_.clear(); parsedProperties_.clear(); }

};


inline ptree getChild(const ptree& pt, const string& name) { return pt.get_child(name, ptree()); }
inline auto getChildRange(const ptree& pt, const string& name) -> decltype(pt.equal_range(name)) { return pt.equal_range(name); } 


std::set<string> preprocessConfiguration(std::istream& is, std::ostream& os, const std::string& istreamid);


#endif

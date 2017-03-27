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

#include "global_funcs.hh"
#include "Decorator.hh"
#include "capabilities.hh"
#include "StringSet.hh"

using std::string;
using std::vector;
using std::map;
using std::pair;
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
  virtual const char* what() const throw() override { return (path() + " : " + std::invalid_argument::what()).c_str(); }
};

struct CheckedPropertyMissing : public PathfulException { CheckedPropertyMissing(const string& objid) : PathfulException("Checked property not set", objid) {} };
struct RequestedUnsetValue : public PathfulException { RequestedUnsetValue() : PathfulException("Tried to get value from an unset property") {} };
struct InvalidPropertyValue : public PathfulException { InvalidPropertyValue(const string& objid) : PathfulException("Value for " + objid + " is invalid") {} };
struct InvalidComputable : public PathfulException { InvalidComputable() : PathfulException("Computable object is invalid or unset") {} };


template<class T>
struct Stateful { 
  typedef T State;
  virtual State state() const = 0; 
};

class Validful {
 public:
  virtual bool valid() const { return true ; }
};

struct Parsable : public Stateful<bool>, public Validful {
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
  
  // Readonly property can't be changed after setting, hence to scale it with used unit, one needs an extra method
  void scaleByUnit(const float& unit) { Property<T, ValueHolder>::operator()(Property<T, ValueHolder>::operator()()*unit); }

#define ALLOW_FORCE_SET
#ifdef ALLOW_FORCE_SET
  void force(const T& value) { Property<T, ValueHolder>::operator()(value); } // for testing purposes only
#endif
};


template<typename T, const char Sep = ','>
class PropertyVector : public Parsable {  // CUIDADO DEPRECATED
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
    string seq = trim(s);
    if (seq.front() != Sep) values_.clear(); // an append is only done when the first character is the separator, in other cases we overwrite (example: if Sep is ',' this is an append: ",X,Y")
    auto values = split<T>(seq, string(1, Sep));
    for (const auto& v : values) values_.push_back(v);
  }
  void appendString(const string& s) {
    values_.push_back(trim(s));
  }
}; 

//! Vector of property variables, which are read-only - property can be read-in (or computed) using data from boost property_tree, can't be (re)set further on
template<typename T, const char Sep = ','>
class ReadonlyPropertyVector : public Parsable {  // CUIDADO DEPRECATED
  std::vector<T> m_values;
  const string&  m_name;

  void operator()(size_t i, const T& value) {}
public:
  ReadonlyPropertyVector(const string& name, PropertyMap& registrar, const std::initializer_list<T>& values = {}) : m_values(values), m_name(StringSet::ref(name)) { registrar[name] = this; }
  ReadonlyPropertyVector(const string& name, const std::initializer_list<T>& values = {}) : m_values(values), m_name(StringSet::ref(name)) {}
  ReadonlyPropertyVector(const std::initializer_list<T>& values = {}) : m_values(values), m_name(StringSet::ref("unnamed")) {}

  const T& operator()(size_t i) const       { return m_values[i]; }
  const T& operator[](size_t i) const       { return m_values[i]; }
  size_t size() const                       { return m_values.size(); }

  typename std::vector<T>::const_iterator begin() const { return m_values.begin(); }
  typename std::vector<T>::const_iterator end()   const { return m_values.end(); }

  // Readonly property can't be changed after setting, hence to scale it with used unit, one needs an extra method
  void scaleByUnit(const float& unit) { for (size_t i=0; i<m_values.size(); i++) m_values[i] *= unit; }

  bool state() const  { return !m_values.empty(); }
  void clear()        { m_values.clear(); }
  string name() const { return m_name; }

  void fromPtree(const ptree& pt)    { fromString(pt.data()); }
  void fromString(const string& s)   {

    string seq = trim(s);
    if (seq.front() != Sep) m_values.clear(); // an append is only done when the first character is the separator, in other cases we overwrite (example: if Sep is ',' this is an append: ",X,Y")
    std::vector<T> values = split<T>(seq, string(1, Sep));
    for (const auto& v : values) m_values.push_back(v);
  }
};

template<typename T, const char Sep = ','>
class MultiProperty : public T, public Parsable {
  const string& name_;
  typedef typename std::decay<decltype(*std::declval<T>().begin())>::type ValueType;
public:
 MultiProperty(const string& name, PropertyMap& registrar, const T& valueHolder = T()) : T(valueHolder), name_(StringSet::ref(name))  { registrar[name] = this; }
 MultiProperty(const string& name, const T& valueHolder = T()) : T(valueHolder), name_(StringSet::ref(name)) {}
 MultiProperty(const T& valueHolder = T()) : T(valueHolder), name_(StringSet::ref("unnamed")) {}
  bool state() const { return !T::empty(); }
  void clear() { T::clear(); }
  //const T& operator()() const { return values_; }
  //T& operator()() const { return values_; }
  //void operator()(const T& values) { values_ = values; }
  //typename T::const_iterator begin() const { return values_.begin(); }
  //typename T::const_iterator end() const { return values_.end(); }
  //typename T::iterator begin() { return values_.begin(); }
  //typename T::iterator end() { return values_.end(); }
  //void clear() { values_.clear(); }
  string name() const { return name_; }
  void fromPtree(const ptree& pt) { fromString(pt.data()); }
  void fromString(const string& s) { 
    string seq = trim(s);
    if (seq.front() != Sep) clear(); // an append is only done when the first character is the separator, in other cases we overwrite (example: if Sep is ',' this is an append: ",X,Y")
    auto values = split<ValueType>(seq, string(1, Sep));
    for (const auto& v : values) T::insert(T::end(), v);
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
    bool weak = pt.data().front() == '_'; // a weak key doesn't cause insertion if a key with the same name is not already present (be mindful of the include order when using weak keys)
    T key = str2any<T>(!weak ? pt.data() : pt.data().substr(1)); // strip leading '_' in case of weak key
    if (this->count(key) == 0 && !weak) this->insert(make_pair(key, pt)); 
    else if (this->count(key) > 0) for (auto& tel : pt) this->at(key).add_child(tel.first, tel.second);
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
    bool weak = pt.data().front() == '_'; // weak key check ('_' can only appear as the leading character, even if the key is an interval and/or a sequence)
    auto tokens = split(!weak ? pt.data() : pt.data().substr(1), ","); // split sequences like 1, 2, 3, 4
    for (auto t : tokens) {
      auto interval = split<int>(t, "-"); // split intervals like 1-5
      if (interval.size() == 1) keys.push_back(interval[0]); // if size == 1, it means no '-' was found and the string was actually not an interval
      else for (int i = interval[0]; i <= interval[1]; i++) keys.push_back(i); // explode the interval
    }
    for (auto k : keys) {
      if (this->count(k) == 0 && !weak) this->insert(make_pair(k, pt)); 
      else if (this->count(k) > 0) for (auto& tel : pt) this->at(k).add_child(tel.first, tel.second);
    }
  }
  void fromString(const string& s) { this->insert(make_pair(str2any<int>(s), ptree())); }
};

template<typename T>
class PropertyNodeUnique : public Parsable, public vector<pair<T, ptree> > {
  const string& name_;
public:
  PropertyNodeUnique(const string& name, PropertyMap& registrar) : name_(StringSet::ref(name)) { registrar[name] = this; }
  PropertyNodeUnique(const string& name) : name_(StringSet::ref(name)) {}
  bool state() const { return !this->empty(); }
  void clear() { vector<pair<T, ptree> >::clear(); }
  string name() const { return name_; }
  void fromPtree(const ptree& pt) { this->push_back(make_pair(pt.data(), pt)); }
  void fromString(const string& s) { this->push_back(make_pair(str2any<T>(s), ptree())); }
};

template <typename T, int numElem,  const char Sep = ','>
class FixedSizeMultiProperty : public MultiProperty<T, Sep> {
  public:
  FixedSizeMultiProperty(const string& name, PropertyMap& registrar) : MultiProperty<T, Sep>(name, registrar) { };
  bool valid() const override {
    int mySize = MultiProperty<T, Sep>::size();
    bool result = (mySize == numElem);
    return result;
  }
};

template <typename T> using RangeProperty = FixedSizeMultiProperty<T, 2, '-'>;

typedef ptree PropertyTree;


class PropertyObject {
  PropertyMap parsedCheckedProperties_, checkedProperties_, parsedProperties_;
  PropertyTree pt_;
  static std::set<string> globalMatchedProperties_;
  static std::set<string> globalUnmatchedProperties_;

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

  void recordMatchedProperties() {
    for (auto& mapel : parsedCheckedProperties_) globalMatchedProperties_.insert(mapel.first);
    for (auto& mapel : parsedProperties_) globalMatchedProperties_.insert(mapel.first);
    for (auto& mapel : checkedProperties_) globalMatchedProperties_.insert(mapel.first);
    for (auto& trel : pt_) globalUnmatchedProperties_.insert(trel.first);
  }

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
    recordMatchedProperties();
  }
  virtual void check() {
    for (auto& v : parsedProperties_) {
      if (v.second->state() && !v.second->valid()) throw InvalidPropertyValue(v.first);
    }
    for (auto& v : parsedCheckedProperties_) {
      if (!v.second->state()) throw CheckedPropertyMissing(v.first); 
      if (v.second->state() && !v.second->valid()) throw InvalidPropertyValue(v.first);
    }
    for (auto& v : checkedProperties_) {
      if (!v.second->state()) throw CheckedPropertyMissing(v.first); 
      if (v.second->state() && !v.second->valid()) throw InvalidPropertyValue(v.first);
    }
  }

  virtual void cleanup() { pt_.clear(); parsedCheckedProperties_.clear(); parsedProperties_.clear(); }
  virtual void cleanupTree() { pt_.clear(); }

  static std::set<string> reportUnmatchedProperties() {
    std::set<string> unmatched;
    std::set_difference(globalUnmatchedProperties_.begin(), globalUnmatchedProperties_.end(),
                        globalMatchedProperties_.begin(), globalMatchedProperties_.end(),
                        std::inserter(unmatched, unmatched.end()));
    return unmatched;
  }

};

inline ptree getChild(const ptree& pt, const string& name) { return pt.get_child(name, ptree()); }
inline auto getChildRange(const ptree& pt, const string& name) -> decltype(pt.equal_range(name)) { return pt.equal_range(name); } 


std::set<string> preprocessConfiguration(std::istream& is, std::ostream& os, const std::string& istreamid);


#endif

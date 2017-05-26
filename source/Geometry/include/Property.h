#ifndef INCLUDE_PROPERTY_H
#define INCLUDE_PROPERTY_H

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

#include "GeometryCapability.h"
#include "string_functions.h"
#include "StringConverter.h"
#include "StringSet.h"

// Used namespaces
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::list;
using std::shared_ptr;

using boost::property_tree::ptree;

// Typedefs
class Parsable;
typedef std::map<string, Parsable*> PropertyMap;
typedef ptree PropertyTree;

inline ptree getChild(const ptree& pt, const string& name) { return pt.get_child(name, ptree()); }
inline auto  getChildRange(const ptree& pt, const string& name) -> decltype(pt.equal_range(name)) { return pt.equal_range(name); }

std::set<string> preprocessConfiguration(std::istream& is, std::ostream& os, const std::string& istreamid);

/*
 * Set of classes used to read in parameters via boost property_tree library. Individual variables read-in as so called properties can be of
 * two given categories:
 *  Property<(int,double,...), type> or ReadonlyProperty<(int,double,...), TYPE),
 *
 *  where TYPE stands for the following:
 *  Default -> initialized by default (default value expected to be set in a constructor) -> overwritten if value provided via PropertyObject store() call using property_tree as input
 *  AutoDefault -> initialized by default (default value expected to be set in a constructor)
 *  NoDefault -> not initialized by default -> value expected to be provided via PropertyObject store() call using property_tree as input
 *  CachedComputable -> value computed through function (lambda function) set by setup() call -> calculated only once (value cached) when its value getter called -> use clear() to uncache
 *  UncachedComputable -> value computed through function (lambda function) set by setup() call -> calculated every time the value getter called
 *
 *  Such variables can be private or public in a class. Their setter method is simply given as variableName(value) and getter method as variableName(). These
 *  variables are automatically
 *
 *  All geometry objects (Trackers -> Barrels, Endcaps ...) must derive from PropertyObject class to be able to read-in all geometry parameters through geometry configuration files.
 */

// ===================================================
// HELPER CLASSES: exceptions, states, ...
// ===================================================

//! Exception class to be used when using boost property_tree
class PathfulException : public std::invalid_argument {
  string m_path;
public:
  PathfulException(const string& what, const string& path) : std::invalid_argument(what.c_str()) { pushPath(path); }
  template<class T> PathfulException(const string& what, const T& obj, const string& objid) : std::invalid_argument(what.c_str()) { pushPath(obj, objid); }
  PathfulException(const string& what) : std::invalid_argument(what.c_str()) {}
  void pushPath(const string& p) { m_path = p + (!m_path.empty() ? "." + m_path : ""); }
  template<class T> void pushPath(const T& obj, const string& objid) { pushPath(string(typeid(obj).name()) + "(" + objid + ")"); }
  template<class T, class U> void pushPath(const T& obj, const U& objid) { pushPath(string(typeid(obj).name()) + "(" + any2str(objid) + ")"); }
  const string& path() const { return m_path; }
  virtual const char* what() const throw() override { return (path() + " : " + std::invalid_argument::what()).c_str(); }
};

struct CheckedPropertyMissing : public PathfulException { CheckedPropertyMissing(const string& objid) : PathfulException("Checked property not set", objid) {} };
struct RequestedUnsetValue    : public PathfulException { RequestedUnsetValue() : PathfulException("Tried to get value from an unset property") {} };
struct InvalidPropertyValue   : public PathfulException { InvalidPropertyValue(const string& objid) : PathfulException("Value for " + objid + " is invalid") {} };
struct InvalidComputable      : public PathfulException { InvalidComputable() : PathfulException("Computable object is invalid or unset") {} };

//! Pure virtual class -> derived objects can have given state
template<class T>
struct Stateful { 
  typedef T State;
  virtual State state() const = 0; 
  virtual ~Stateful() {};
};

//! Pure virtual class -> derived objects wil be valid/non-valid, e.g. property variables are defined/non-defined
class Validful {
 public:
  virtual bool valid() const { return true ; }
  virtual ~Validful() {};
};

//! Pure virtual class -> derived objects will be parsable -> boost property_tree is parsed to get all its data
struct Parsable : public Stateful<bool>, public Validful {
  virtual string name() const = 0;
  virtual void fromPtree(const ptree& pt) = 0;
  virtual void fromString(const string& s) = 0;
};

// ===================================================
// PROPERTY TYPES:
// ===================================================

//! Property base class -> derived objects will have getters, setters and clear method
template<typename T>
struct PropertyBase : public Stateful<bool> {
  virtual void operator()(const T&) = 0;
  virtual const T& operator()() const = 0;
  virtual void clear() = 0;
};

//! Property value initialized by default -> default value expected to be set in a constructor -> overwritten if value provided via PropertyObject store() call using property_tree as input
template<typename T>
class Default : public PropertyBase<T> {
  T       m_value;
  const T m_default;
public:
  Default(const T& value) : m_value(value), m_default(value) {}
  void operator()(const T& value) { m_value = value; }
  const T& operator()() const     { return m_value; }
  bool state() const              { return true; }
  void clear()                    { m_value = m_default; }
};

//! Property value initialized by defaul -> default value expected to be set in a constructor
template<typename T>
struct AutoDefault : public Default<T> { AutoDefault() : Default<T>(T()) {} };

//! Property value not automatically initialized -> value expected to be provided via PropertyObject store() call using property_tree as input
template<typename T>
class NoDefault : public PropertyBase<T> {
  T    m_value;
  bool m_state;
public:
  NoDefault() : m_state(false) {}
  void operator()(const T& value) { m_value = value; m_state = true; }
  const T& operator()() const { 
    if (!m_state) throw RequestedUnsetValue();
    return m_value;
  }
  bool state() const { return m_state; }
  void clear()       { m_state = false; }
};

//! Property value computed computed through function (lambda function) set by setup() call -> calculated only once (value cached) when its value getter called -> use clear() to uncache
template<typename T>
class Computable : public PropertyBase<T> {
  typedef std::function<T()> Func;
  Func get;
  mutable T    m_value;
  mutable bool m_state;
public:
  template<class U = Func> Computable(const U& getter /*= []()->T{ throw InvalidComputable(); }*/) : get(Func(getter)), m_state(false) {} // default arg commented out due to bug in gcc 4.7.2-5
  Computable() : get([]()->T{ throw InvalidComputable(); }), m_state(false) {}
  Computable(const Computable<T>& other) : Computable() { m_value = other.m_value; m_state = other.m_state; } // Func does not get copied!! needs to be manually setup() again to prevent issues with the captures (capture happens at point of declaration)
  void operator()(const T& value) { m_value = value; m_state = true; }
  const T& operator()() const {
    if (!m_state) { m_value = get(); m_state = true; }
    return m_value;
  }
  bool state() const { return m_state; }
  void clear()       { m_state = false; }
  template<class U = Func> void setup(const U& getter) { get = Func(getter); }
};

//! Property value computed computed through function (lambda function) set by setup() call -> calculated every time the value getter called
template<typename T>
class UncachedComputable : public PropertyBase<T> {
  typedef std::function<T()> Func;
  Func get;
  mutable T m_value;
  mutable bool m_state;
public:
  template<class U = Func> UncachedComputable(const U& getter /*= []()->T{ throw InvalidComputable(); }*/) : get(Func(getter)), m_state(false) {}
  UncachedComputable() : get([]()->T{ throw InvalidComputable(); }), m_state(false) {}
  UncachedComputable(const UncachedComputable<T>& other) : UncachedComputable() { m_value = other.m_value; m_state = other.m_state; } // Func does not get copied!! needs to be manually setup() again to prevent issues with the captures (capture happens at point of declaration)
  void operator()(const T& value) { m_value = value; m_state = true; }
  const T& operator()() const { return m_state ? m_value : m_value = get(); }
  bool state() const { return m_state; }
  void clear() { m_state = false; }
  template<class U = Func> void setup(const U& getter) { get = Func(getter); }
};

// ===================================================
// PROPERTIES: variables used with boost property_tree
// ===================================================

//! Property variable - property can be read-in (or computed) using data from boost property_tree or set on the fly
template<typename T, template<typename> class ValueHolder>
class Property : public PropertyBase<T>, public Parsable {
  ValueHolder<T> m_valueHolder;
  const string& m_name;
public:
  Property(const string& name, PropertyMap& registrar, const ValueHolder<T>& valueHolder = ValueHolder<T>()) : m_valueHolder(valueHolder), m_name(StringSet::ref(name)) { registrar[name] = this; }
  Property(const string& name, const ValueHolder<T>& valueHolder = ValueHolder<T>()) : m_valueHolder(valueHolder), m_name(StringSet::ref(name)) {}
  Property(const ValueHolder<T>& valueHolder = ValueHolder<T>()) : m_valueHolder(valueHolder), m_name(StringSet::ref("unnamed")) {}

  void operator()(const T& value) { m_valueHolder(value); }
  const T& operator()() const { 
    try { return m_valueHolder(); }
    catch(PathfulException& pe) { pe.pushPath(m_name); throw; }
  }
  template<class ...U> void setup(const U&... valueHolderArgs) { m_valueHolder.setup(valueHolderArgs...); }

  bool state() const  { return m_valueHolder.state(); }
  void clear()        { m_valueHolder.clear(); }
  string name() const { return m_name; }

  void fromPtree(const ptree& pt)  { m_valueHolder(str2any<T>(pt.data())); }
  void fromString(const string& s) { m_valueHolder(str2any<T>(s)); }
};

//! Property variable (readonly option) - property can be read-in (or computed) using data from boost property_tree, can't be (re)set further on
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

//! Vector of property variables
template<typename T, const char Sep = ','>
class PropertyVector : public Parsable {  // CUIDADO DEPRECATED
  std::vector<T> m_values;
  const string& m_name;
public:
  PropertyVector(const string& name, PropertyMap& registrar, const std::initializer_list<T>& values = {}) : m_values(values), m_name(StringSet::ref(name)) { registrar[name] = this; }
  PropertyVector(const string& name, const std::initializer_list<T>& values = {}) : m_values(values), m_name(StringSet::ref(name)) {}
  PropertyVector(const std::initializer_list<T>& values = {}) : m_values(values), m_name(StringSet::ref("unnamed")) {}

  void operator()(size_t i, const T& value) { m_values[i] = value; }
  const T& operator()(size_t i) const       { return m_values[i]; }
  const T& operator[](size_t i) const       { return m_values[i]; }
  size_t size() const                       { return m_values.size(); }

  typename std::vector<T>::const_iterator begin() const { return m_values.begin(); }
  typename std::vector<T>::const_iterator end()   const { return m_values.end(); }

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

  void appendString(const string& s) { m_values.push_back(trim(s));}
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

//! Map of property variables
template<typename T, const char Sep = ','>
class MultiProperty : public T, public Parsable {
  const string& m_name;
  typedef typename std::decay<decltype(*std::declval<T>().begin())>::type ValueType;
public:
 MultiProperty(const string& name, PropertyMap& registrar, const T& valueHolder = T()) : T(valueHolder), m_name(StringSet::ref(name))  { registrar[name] = this; }
 MultiProperty(const string& name, const T& valueHolder = T()) : T(valueHolder), m_name(StringSet::ref(name)) {}
 MultiProperty(const T& valueHolder = T()) : T(valueHolder), m_name(StringSet::ref("unnamed")) {}
  bool state() const { return !T::empty(); }
  void clear() { T::clear(); }
  //const T& operator()() const { return m_values; }
  //T& operator()() const { return m_values; }
  //void operator()(const T& values) { m_values = values; }
  //typename T::const_iterator begin() const { return m_values.begin(); }
  //typename T::const_iterator end() const { return m_values.end(); }
  //typename T::iterator begin() { return m_values.begin(); }
  //typename T::iterator end() { return m_values.end(); }
  //void clear() { m_values.clear(); }
  string name() const { return m_name; }
  void fromPtree(const ptree& pt) { fromString(pt.data()); }
  void fromString(const string& s) { 
    string seq = trim(s);
    if (seq.front() != Sep) clear(); // an append is only done when the first character is the separator, in other cases we overwrite (example: if Sep is ',' this is an append: ",X,Y")
    auto values = split<ValueType>(seq, string(1, Sep));
    for (const auto& v : values) T::insert(T::end(), v);
  }
};

// ===================================================
// COLLECTIONS OF PROPERTIES: Nodes & Trees (Objects)
// ===================================================

//! Property node (templated): i.e. subtree object
template<typename T>
class PropertyNode : public Parsable, public map<T, ptree> { 
  const string& m_name;
public:
  PropertyNode(const string& name, PropertyMap& registrar) : m_name(StringSet::ref(name)) { registrar[name] = this; }
  PropertyNode(const string& name) : m_name(StringSet::ref(name)) {}
  bool state() const { return !this->empty(); }
  void clear() { map<T, ptree>::clear(); }
  string name() const { return m_name; }
  void fromPtree(const ptree& pt) { 
    bool weak = pt.data().front() == '_'; // a weak key doesn't cause insertion if a key with the same name is not already present (be mindful of the include order when using weak keys)
    T key = str2any<T>(!weak ? pt.data() : pt.data().substr(1)); // strip leading '_' in case of weak key
    if (this->count(key) == 0 && !weak) this->insert(make_pair(key, pt)); 
    else if (this->count(key) > 0) for (auto& tel : pt) this->at(key).add_child(tel.first, tel.second);
  }
  void fromString(const string& s) { this->insert(make_pair(str2any<T>(s), ptree())); }
};

//! Property node (not templated -> using int instead): i.e. subtree object
template<>
class PropertyNode<int> : public Parsable, public map<int, ptree> {
  const string& m_name;
public:
  PropertyNode(const string& name, PropertyMap& registrar) : m_name(StringSet::ref(name)) { registrar[name] = this; }
  PropertyNode(const string& name) : m_name(StringSet::ref(name)) {}
  bool state() const { return !this->empty(); }
  void clear() { map<int, ptree>::clear(); }
  string name() const { return m_name; }
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

//! Property unique node
template<typename T>
class PropertyNodeUnique : public Parsable, public vector<pair<T, ptree> > {
  const string& m_name;
public:
  PropertyNodeUnique(const string& name, PropertyMap& registrar) : m_name(StringSet::ref(name)) { registrar[name] = this; }
  PropertyNodeUnique(const string& name) : m_name(StringSet::ref(name)) {}
  bool state() const { return !this->empty(); }
  void clear() { vector<pair<T, ptree> >::clear(); }
  string name() const { return m_name; }
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

// ===================================================
//
// MAIN PROPERTY CLASS: PropertyObject
//
// ===================================================

//! Property tree object -> all geometry objects reading in configuration from file (in tree format) derive from this class
class PropertyObject {
  PropertyMap parsedCheckedProperties_, checkedProperties_, parsedProperties_;
  PropertyTree m_pTree;
  static std::set<string> globalMatchedProperties_;
  static std::set<string> globalUnmatchedProperties_;

  void processProperties(PropertyMap& props) {
    for (auto& propElem : props) {
      auto childRange = m_pTree.equal_range(propElem.first);
      std::for_each(childRange.first, childRange.second, [&propElem](const ptree::value_type& treeElem) {
        propElem.second->fromPtree(treeElem.second); // takes care of duplicate entries (by overwriting the property value as many times as there are entries with the same key) and of node entries (in that case the PropertyNodes differentiates based on the value)
      });
      m_pTree.erase(propElem.first);
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
  const PropertyTree& propertyTree() const { return m_pTree; }
  PropertyMap& parsedAndChecked() { return parsedCheckedProperties_; }
  PropertyMap& parsedOnly() { return parsedProperties_; }
  PropertyMap& checkedOnly() { return checkedProperties_; }

  void recordMatchedProperties() {
    for (auto& mapel : parsedCheckedProperties_) globalMatchedProperties_.insert(mapel.first);
    for (auto& mapel : parsedProperties_) globalMatchedProperties_.insert(mapel.first);
    for (auto& mapel : checkedProperties_) globalMatchedProperties_.insert(mapel.first);
    for (auto& trel : m_pTree) globalUnmatchedProperties_.insert(trel.first);
  }

public:
  PropertyObject() {}
  virtual ~PropertyObject() {}
  virtual void store(const PropertyTree& newpt) {
    if (m_pTree.empty()) m_pTree = newpt;
    else { 
      m_pTree.data() = newpt.data();
      for (auto& propElem : newpt) { 
        m_pTree.add_child(propElem.first, propElem.second);
      } // merging trees in a careless manner, appending children without ever checking if an entry with the same key is already present. the duplicates thus formed will be all grabbed at parsing time by the properties (each duplicate entry overwrites the previous)
    }
    //std::cout << "============ " << m_pTree.data() << " ===========" << std::endl;
    //printAll(m_pTree);
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
  virtual void evaluateProperty(Parsable& propElem) {

    //printAll(m_pTree);
    auto childRange = m_pTree.equal_range(propElem.name());
    std::for_each(childRange.first, childRange.second, [&propElem](const ptree::value_type& treeElem) {
        propElem.fromPtree(treeElem.second); // takes care of duplicate entries (by overwriting the property value as many times as there are entries with the same key) and of node entries (in that case the PropertyNodes differentiates based on the value)      });
    });
  }
  virtual void cleanup() { m_pTree.clear(); parsedCheckedProperties_.clear(); parsedProperties_.clear(); }
  virtual void cleanupTree() { m_pTree.clear(); }

  static std::set<string> reportUnmatchedProperties() {
    std::set<string> unmatched;
    std::set_difference(globalUnmatchedProperties_.begin(), globalUnmatchedProperties_.end(),
                        globalMatchedProperties_.begin(), globalMatchedProperties_.end(),
                        std::inserter(unmatched, unmatched.end()));
    return unmatched;
  }
}; // Class

#endif /* INCLUDE_PROPERTY_H */

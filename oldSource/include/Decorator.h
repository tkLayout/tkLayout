#include <set>
#include <typeinfo>
#include <typeindex>

#define DECORATOR_RTTI_API


template<class T>
class Decorator {
  T* const decorated_; // the decorated type must be cloneable (i.e. implement a clone method to copy-construct a new instance on the heap)
protected:
  Decorator(T* const decorated) : decorated_(decorated) {}
public:
  typedef T Decorated;
  Decorator(const Decorator<T>& other) : decorated_(((T*)other.decorated_->clone())) {}
  ~Decorator() { delete decorated_; }
  T& decorated() { return *decorated_; }
  const T& decorated() const { return *decorated_; }

#ifdef DECORATOR_RTTI_API
  template<class U> bool is() const { return typeid(*this) == typeid(U) || this->decorated().template is<U>(); }
  bool isAll(std::set<std::type_index>& types) const {
    types.erase(typeid(*this));
    return types.size() > 1 ? this->decorated().isAll(types) : false;
  }
  template<class ...UU> bool isAll() const {
    std::set<std::type_index> types({typeid(UU)...});
    return isAll(types);
  }

  template<class U> const U* as() const { return typeid(*this) == typeid(U) ? static_cast<U* const>(this) : this->decorated().template as<U>(); }
  template<class U> U* as() { return typeid(*this) == typeid(U) ? static_cast<U* const>(this) : this->decorated().template as<U>(); }
#endif
};

class Decorable {
public:
#ifdef DECORATOR_RTTI_API
  template<class U> bool is() const { return typeid(*this) == typeid(U); }
  bool isAll(std::set<std::type_index>& types) const { return types.size() == 1 && types.count(typeid(*this)); }
  template<class U> const U* as() const { return typeid(*this) == typeid(U) ? reinterpret_cast<U* const>(this) : NULL; }
  template<class U> U* as() { return typeid(*this) == typeid(U) ? reinterpret_cast<U* const>(this) : NULL; }
#endif
};


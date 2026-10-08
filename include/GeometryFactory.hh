#ifndef GEOMETRY_FACTORY
#define GEOMETRY_FACTORY

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

class GeometryFactory {
  template<class T> static void conditionalSetup(T* t, typename std::enable_if<std::is_void<decltype(t->setup())>::value>::type* = 0) { t->setup(); } 
  static void conditionalSetup(...) {}
public:
  template<class T> static T* clone(const T& t) { return make<T>(t); }
  template<class T, class ...U> static T* make(const U&... args) {
    T* t = new T(args...);
    conditionalSetup(t);
    return t;
  }
};

struct FactoryCloneAllocator {
  template<class U> static U* allocate_clone(const U& r) { return GeometryFactory::clone(r); } // replaces cloning operation for boost ptr_containers (so that setup() is called where needed)
  template<class U> static void deallocate_clone(const U* r) { boost::delete_clone(r); } // standard boost clone deleter
};

template<class T> using PtrVector = boost::ptr_vector<T, FactoryCloneAllocator>;

template<class K, class V> using PtrMap = boost::ptr_multimap<K, V, std::less<K>, FactoryCloneAllocator>;


#endif

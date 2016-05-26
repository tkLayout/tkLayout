#ifndef INCLUDE_GEOMETRY_FACTORY
#define INCLUDE_GEOMETRY_FACTORY

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

/*
 * Factory class is meant to be used to build individual geometrical objects (barrels, layers, rods, ...)
 * and put them into boost::ptr_vector container, named as PtrVector. (PtrVector deals correctly with the
 * memory management if individual pointers get out of scope, one access vector objects by reference using
 * at() instead of pointers, automatic copy of vector objects is guaranteed if cloned). Use either
 * GeometryFactory::make<GeomObjectType>(constructor parameters) method as new operator, or
 * GeometryFactory::clone<GeomObjectType>(objectInstance) to clone the object. Both methods call setup()
 * method, which is meant to be set by each GeomObjectType classes. That's why the wrapper around
 * boost::ptr_vector allocation has been written.
 */
class GeometryFactory {

 public:

  //! Clone geometrical object: GeometryFactory::clone<GeomObjectType T>(objectInstance t)
  template<class T> static T* clone(const T& t) { return make<T>(t); }

  //! Create new geometrical object and allocate memory on heap: GeometryFactory::make<GeomObjectType T>(T constructor parameters)
  template<class T, class ...U> static T* make(const U&... args) {

    T* t = new T(args...);
    conditionalSetup(t);
    return t;
  }

 private:

  //! Call setup() method for given geometrical object -> called by clone() or make() statical methods
  template<class T> static void conditionalSetup(T* t, typename std::enable_if<std::is_void<decltype(t->setup())>::value>::type* = 0) { t->setup(); }

  static void conditionalSetup(...) {}
};

//! Replaces standard boost clone deleter & cloning operation for boost ptr_containers, so that setup() is called wherever needed
struct FactoryCloneAllocator {

  template<class U> static U* allocate_clone(const U& r) { return GeometryFactory::clone(r); } // replaces cloning operation for boost ptr_containers (so that setup() is called where needed)
  template<class U> static void deallocate_clone(const U* r) { boost::delete_clone(r); }       // standard boost clone deleter
};

//! Vector of pointers using boost ptr_vector memory management
template<class T> using PtrVector = boost::ptr_vector<T, FactoryCloneAllocator>;

//! Map of pointers using boost ptr_vector memory management
template<class K, class V> using PtrMap = boost::ptr_multimap<K, V, std::less<K>, FactoryCloneAllocator>;


#endif /* INCLUDE_GEOMETRY_FACTORY */

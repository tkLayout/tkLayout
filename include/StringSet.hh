#ifndef STRINGSET_H
#define STRINGSET_H

#include <set>
#include <string>

class StringSet {
  std::set<std::string> strings_;
  StringSet() {}
public:
  static StringSet& instance() {
    static StringSet ss;
    return ss;
  }

  static const std::string& ref(const std::string& s) {
    return instance().makeRef(s);
  }

  const std::string& makeRef(const std::string& s) { 
    return *strings_.insert(s).first;
  }

};



#endif

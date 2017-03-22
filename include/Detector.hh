#ifndef TRACKERBUNDLE_HH
#define TRACKERBUNDLE_HH

#include <vector>

class Tracker;

class Detector : public Visitable {
private:
  std::vector<Tracker&> trackers_;
public:
  int addTracker(Tracker&) { trackers_.push_back(tracker); };
  const std::vector<Tracker&>& trackers() { return trackers_; };
  void accept(GeometryVisitor& v) {
    v.visit(*this);
    for (auto& t : trackers_) { t.accept(v); }
  }
  void accept(ConstGeometryVisitor& v) const {
    v.visit(*this);
    for (const auto& t : trackers_) { t.accept(v); }
  }
  void accept(SensorGeometryVisitor& v) {
    v.visit(*this);
    for (auto& t : trackers_) { t.accept(v); }
  }
};

#endif

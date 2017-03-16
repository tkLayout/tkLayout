#ifndef REPORT_HH
#define REPORT_HH

#include <vector>

class Tracker;
class TrackerBundle;
class RootWSite;
class RootWPage;
class RootWContent;

class Report {
public:
  virtual void analyze(Tracker&) { };
  virtual void analyze(TrackerBundle&) { };
  virtual void analyze(const Tracker&) { };
  virtual void analyze(const TrackerBundle&) { };
  virtual void visualize(RootWSite&) { };
  virtual void visualize(RootWPage&) { };
  virtual void visualize(RootWContent&) { };
};

#endif

#ifndef REPORT_HH
#define REPORT_HH

#include <vector>

class TrackerBundle;
class Tracker;
class RootWSite;
class RootWPage;
class RootWContent;

class Report {
public:
  virtual void analyze(TrackerBundle&) { };
  virtual void analyze(Tracker&) { };
  virtual void analyze(const TrackerBundle&) { };
  virtual void analyze(const Tracker&) { };
  virtual void visualizeTo(RootWSite&) { };
  virtual void visualizeTo(RootWPage&) { };
  virtual void visualizeTo(RootWContent&) { };
};

#endif

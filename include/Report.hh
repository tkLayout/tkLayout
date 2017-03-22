#ifndef REPORT_HH
#define REPORT_HH

#include <vector>

class Detector;
class Tracker;
class RootWSite;
class RootWPage;
class RootWContent;

class Report {
public:
  virtual void analyze() { };  
  virtual void analyze(Detector&) { };
  virtual void analyze(Tracker&) { };
  virtual void analyze(const Detector&) { };
  virtual void analyze(const Tracker&) { };
  virtual void visualizeTo(RootWSite&) { };
  virtual void visualizeTo(RootWPage&) { };
  virtual void visualizeTo(RootWContent&) { };
};

#endif

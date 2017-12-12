#ifndef REPORTMODULECOUNT_HH
#define REPORTMODULECOUNT_HH

#include <string>
#include <map>

#include "AnalyzerVisitors/ModuleCount.hh"
#include <Report.hh>
#include "RootWeb.hh"

class DetectorModule;
class Tracker;

class ReportModuleCount : public Report {
private:
  ModuleCountVisitor moduleCounter_;
public:
  void analyze(const Tracker&);
  void visualizeTo(RootWContent&);
};

#endif

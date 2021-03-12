#ifndef REPORTMODULECOUNT_HH
#define REPORTMODULECOUNT_HH

#include <map>
#include <string>

#include "AnalyzerVisitors/ModuleCount.hh"
#include "RootWeb.hh"
#include <Report.hh>

class DetectorModule;
class Tracker;

class ReportModuleCount : public Report {
private:
  ModuleCountVisitor moduleCounter_;

public:
  void analyze(const Tracker &);
  void visualizeTo(RootWContent &);
};

#endif

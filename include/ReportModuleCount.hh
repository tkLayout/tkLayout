#ifndef REPORTMODULECOUNT_HH
#define REPORTMODULECOUNT_HH

#include <string>
#include <map>

#include <Visitor.h>
#include <Report.hh>
#include "rootweb.hh"

class DetectorModule;
class Tracker;

// The visitor counting module per types and subdetector
class ModuleCounterVisitor : public ConstGeometryVisitor {
private:
  bool firstVisit = true;
  std::map<std::string, int> moduleTypes;
  std::map<std::string, int> subDetectors;
  std::map<std::pair<std::string, std::string>, int> count_TypeSub;
  std::string currentSubdetector;
  int currentSubdetectorIndex;
  std::string moduleSummaryType(const DetectorModule& m) const;
  void sortTypesAndDetectors();
public:
  RootWTable* makeTable();
  void preVisit(const Tracker& tracker);
  void visit(const Barrel& barrel);
  void visit(const Endcap& endcap);
  void visit(const DetectorModule& m);
};

class ReportModuleCount : public Report {
private:
  ModuleCounterVisitor moduleCounter_;
public:
  void analyze(const Tracker&);
  void visualize(RootWContent&);
};

#endif

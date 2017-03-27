#ifndef REPORTMODULECOUNT_HH
#define REPORTMODULECOUNT_HH

#include <string>
#include <map>

#include <Visitor.hh>
#include <Report.hh>
#include "RootWeb.hh"

class DetectorModule;
class Tracker;

// The visitor counting module per types and subdetector
class VisitorModuleCount : public ConstGeometryVisitor {
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
  VisitorModuleCount moduleCounter_;
public:
  void analyze(const Tracker&);
  void visualizeTo(RootWContent&);
};

#endif

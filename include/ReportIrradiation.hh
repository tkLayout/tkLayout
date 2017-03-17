#ifndef REPORTIRRADIATION_HH
#define REPORTIRRADIATION_HH

#include <string>
#include <map>

#include <Visitor.h>
#include <Report.hh>
#include "rootweb.hh"
#include "SummaryTable.h"
#include "AnalyzerVisitors/IrradiationPower.h"

class DetectorModule;
class Tracker;

class ReportIrradiation : public Report {
private:
  IrradiationPowerVisitor irradiation_;
  std::map<std::string, SummaryTable> sensorsIrradiationPowerSummary_;
  std::map<std::string, SummaryTable> sensorsIrradiationSummary_;
  void computeIrradiationPowerConsumption(const Tracker&);
  void preparePowerHistograms();

public:
  void analyze(const Tracker&);
  void visualizeTo(RootWSite&);
};

#endif

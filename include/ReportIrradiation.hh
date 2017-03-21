#ifndef REPORTIRRADIATION_HH
#define REPORTIRRADIATION_HH

#include <string>
#include <map>

#include <Visitor.h>
#include <Report.hh>
#include "rootweb.hh"
#include "SummaryTable.h"
#include "VizardTools.hh"
#include "AnalyzerVisitors/IrradiationPower.h"
#include <TH2D.h>
#include "AnalyzerTools.hh"
#include "SimParms.h"

class DetectorModule;
class Tracker;

class ReportIrradiation : public Report, private insur::AnalyzerTools {
private:
  std::map<std::string, SummaryTable> powerSummaries;
  std::map<std::string, SummaryTable> irradiationSummaries;
  SummaryTable irradiationSummaryPerType;
  SummaryTable chipPowerPerType;
  void dumpRadiationTableSummary(RootWPage& myPage, std::map<std::string, SummaryTable>& radiationSummaries, const std::string& title, std::string units);
  std::string createSensorsIrradiationCsv();
  void computeIrradiationPowerConsumption();
  void computeChipPowerConsumptionTable();
  void preparePowerHistograms();
  TH2D sensorsIrradiationPowerMap;
  TH2D totalPowerConsumptionMap;
  Tracker& tracker;
  SimParms& simParms_;
public:
  // TODO: avoid this :-)
  ReportIrradiation(Tracker& t, SimParms& sp) : tracker(t), simParms_(sp) {};
  void analyze();
  void visualizeTo(RootWSite&);
};

#endif

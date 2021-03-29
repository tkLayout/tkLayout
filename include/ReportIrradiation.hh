#ifndef REPORTIRRADIATION_HH
#define REPORTIRRADIATION_HH

#include <map>
#include <string>

#include "AnalyzerTools.hh"
#include "AnalyzerVisitors/IrradiationPower.hh"
#include "RootWeb.hh"
#include "SimParms.hh"
#include "SummaryTable.hh"
#include "VizardTools.hh"
#include <Report.hh>
#include <TH2D.h>
#include <Visitor.hh>

class DetectorModule;
class Tracker;

class ReportIrradiation : public Report, private insur::AnalyzerTools {
private:
  std::map<std::string, SummaryTable> powerSummaries;
  std::map<std::string, SummaryTable> fluenceSummaries;
  SummaryTable fluenceSummaryPerType;
  SummaryTable chipPowerPerType;
  void dumpRadiationTableSummary(
      RootWPage &myPage,
      std::map<std::string, SummaryTable> &radiationSummaries,
      const std::string &title, std::string units);
  std::string createSensorsIrradiationCsv();
  void computeIrradiationPowerConsumption();
  void computeChipPowerConsumptionTable();
  void preparePowerHistograms();
  TH2D sensorsIrradiationPowerMap;
  TH2D totalPowerConsumptionMap;
  Tracker &tracker;

public:
  // TODO: avoid this :-)
  ReportIrradiation(Tracker &t) : tracker(t){};
  void analyze();
  void visualizeTo(RootWSite &);
};

#endif

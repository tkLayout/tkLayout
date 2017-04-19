#ifndef TRIGGERPROCESSORBANDWIDTH_H
#define TRIGGERPROCESSORBANDWIDTH_H

#include <string>
#include <map>
#include <vector>
#include <utility>

#include "TH1.h"
#include "TH2.h"
#include "Math/Point2D.h"
#include "TRandom3.h"
#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"

#include "Tracker.hh"
#include "SimParms.hh"

#include "Visitor.hh"
#include "SummaryTable.hh"

using std::string;
using std::map;
using std::vector;
using std::pair;

namespace AnalyzerHelpers {
  struct Point { double x, y; };
  struct Circle { double x0, y0, r; };
  std::pair<Circle, Circle> findCirclesTwoPoints(const Point& p1, const Point& p2, double r);
  bool isPointInCircle(const Point& p, const Circle& c);
  bool areClockwise(const Point& p1, const Point& p2);

  double calculatePetalAreaMC(const Tracker& tracker, const SimParms& simParms, double crossoverR);
  double calculatePetalAreaModules(const Tracker& tracker, const SimParms& simParms, double crossoverR);
  double calculatePetalCrossover(const Tracker& tracker, const SimParms& simParms);

  bool isModuleInPetal(const DetectorModule& module, double petalPhi, double curvatureR, double crossoverR);
  bool isModuleInCircleSector(const DetectorModule& module, double startPhi, double endPhi);

  bool isModuleInEtaSector(const SimParms& simParms, const Tracker& tracker, const DetectorModule& module, int etaSector); 
  bool isModuleInPhiSector(const SimParms& simParms, const DetectorModule& module, double crossoverR, int phiSector);

}

using namespace AnalyzerHelpers;

class TriggerProcessorBandwidthVisitor : public ConstGeometryVisitor {
  typedef std::map<std::pair<int, int>, int> ProcessorConnections;
  typedef std::map<std::pair<int, int>, double> ProcessorInboundBandwidths;
  typedef std::map<std::pair<int, int>, double> ProcessorInboundStubsPerEvent;
  ProcessorConnections processorConnections_;
  ProcessorInboundBandwidths processorInboundBandwidths_;
  ProcessorInboundStubsPerEvent processorInboundStubsPerEvent_;

  map<string, map<pair<int, int>, double>> &triggerDataBandwidths_, triggerFrequenciesPerEvent_;

  const Tracker* tracker_;
  const SimParms* simParms_;

  int accumulatedLayerOffset_ = 0;
public:
  SummaryTable processorConnectionSummary, processorInboundBandwidthSummary, processorInboundStubPerEventSummary;
  SummaryTable processorCommonConnectionSummary;
  TH1I moduleConnectionsDistribution;
  TH2I processorCommonConnectionMap;
  std::pair<Circle, Circle> sampleTriggerPetal;
  double crossoverR;

  class ModuleConnectionData {
    int phiCpuConnections_, etaCpuConnections_;
    uint32_t detId_;
  public:
    set<std::pair<int, int>> connectedProcessors;
    int phiCpuConnections() const { return phiCpuConnections_; }
    int etaCpuConnections() const { return etaCpuConnections_; }
    uint32_t detId() const { return detId_; }
    int totalCpuConnections() const { return phiCpuConnections_*etaCpuConnections_; }
    void phiCpuConnections(int conn) { phiCpuConnections_ = conn; }
    void etaCpuConnections(int conn) { etaCpuConnections_ = conn; }
    void detId (uint32_t detId) { detId_ = detId; }
    ModuleConnectionData() : phiCpuConnections_(0), etaCpuConnections_(0) {}
  };
  typedef map<const Module*,ModuleConnectionData> ModuleConnectionMap; 
  typedef std::map<std::pair<int, int>, std::set<int> > TriggerSectorMap;

  ModuleConnectionMap moduleConnections;
  TriggerSectorMap sectorMap;

private:
  int numProcEta, numProcPhi;

  double inboundBandwidthTotal = 0.;
  int processorConnectionsTotal = 0;
  double inboundStubsPerEventTotal = 0.;
public:

  TriggerProcessorBandwidthVisitor(map<string, map<pair<int, int>, double>>& triggerDataBandwidths, map<string, map<pair<int, int>, double>>& triggerFrequenciesPerEvent) :
      triggerDataBandwidths_(triggerDataBandwidths),
      triggerFrequenciesPerEvent_(triggerFrequenciesPerEvent)
  {}

  void preVisit();
  void visit(const SimParms& sp);
  void visit(const Tracker& t);
  void visit(const DetectorModule& m);
  void postVisit();
};



#endif

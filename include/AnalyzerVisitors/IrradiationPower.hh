#ifndef IRRADIATIONPOWER_H
#define IRRADIATIONPOWER_H

#include <string>
#include <map>
#include <vector>
#include <utility>

#include "global_constants.hh"
#include "Units.hh"
#include "Tracker.hh"
#include "SimParms.hh"
#include "Visitor.hh"
#include "SummaryTable.hh"

typedef std::tuple<bool, bool, std::string, int, int> ModuleRef;
// Used to identify an irradiated module : all info which matters in that respect.
// bool isBarrel, bool isOuterRadiusRod, std::string barrel/endcap name, int layer/disk number, int ring number

class IrradiationPowerVisitor : public GeometryVisitor {
 private:
  double timeIntegratedLumi_;
  double referenceTemp_;
  double operatingTemp_;
  double alphaParam_;
  double biasVoltage_;
  const IrradiationMapsManager* irradiationMap_;
  std::pair<double, double> getModuleFluenceMeanMax(const IrradiationMapsManager* irradiationMap, const DetectorModule& m);
  const double computeSensorsPower(const double& totalFluence,
				   const double& alphaParam, const double& volume, const double& referenceTemp,
				   const double& operatingTemp, const double& biasVoltage) const;
  bool isBarrel_;
  bool isOuterRadiusRod_;
  std::map<ModuleRef, double> sensorsFluenceMean_;
  std::map<ModuleRef, double> sensorsFluenceMax_;
  std::map<ModuleRef, double> sensorsPowerMean_;
  std::map<ModuleRef, double> sensorsPowerMax_;
  std::map<ModuleRef, int> modulesCounter_;
  std::map<std::string, std::vector<const DetectorModule*> > mapTypeToFluence_;

 public:
  MultiSummaryTable sensorsPowerSummary;
  MultiSummaryTable sensorsFluenceSummary;
  SummaryTable sensorsFluencePerType;
  
  void preVisit();
  void visit(SimParms& sp);
  void visit(Barrel& b);
  void visit(RodPair& r);
  void visit(Endcap& e);
  void visit(DetectorModule& m);
  void postVisit();
};


#endif

#ifndef IRRADIATIONPOWER_H
#define IRRADIATIONPOWER_H

#include <string>
#include <map>
#include <vector>
#include <utility>

#include "global_constants.h"
#include "Tracker.h"
#include "SimParms.h"
#include "Visitor.h"
#include "SummaryTable.h"

class IrradiationPowerVisitor : public GeometryVisitor {
  double timeIntegratedLumi;
  double referenceTemp;
  double operatingTemp;
  double alphaParam;
  double chargeDepletionVoltage;
  const IrradiationMapsManager* irradiationMap_;
public:
  MultiSummaryTable irradiatedPowerConsumptionSummaries;

  void preVisit() {
    irradiatedPowerConsumptionSummaries.clear();   
  }

  void visit(SimParms& sp) {
    timeIntegratedLumi = sp.timeIntegratedLumi();
    referenceTemp    = sp.referenceTemp() + insur::celcius_to_kelvin;   
    alphaParam       = sp.alphaParam();
    irradiationMap_  = &sp.irradiationMapsManager();
  }

  void visit(Barrel& b) {
    irradiatedPowerConsumptionSummaries[b.myid()].setHeader("Layer", "Ring");
    irradiatedPowerConsumptionSummaries[b.myid()].setPrecision(3);        
  }

  void visit(Endcap& e) {
    irradiatedPowerConsumptionSummaries[e.myid()].setHeader("Disk", "Ring");
    irradiatedPowerConsumptionSummaries[e.myid()].setPrecision(3);        
  }

  void visit(DetectorModule& m) {
    operatingTemp    = m.operatingTemp() + insur::celcius_to_kelvin;
    chargeDepletionVoltage = m.chargeDepletionVoltage();
    // <Stefano Mersi>
    // will visit also the modules with z<0, otherwise totals in the summaries will be wrong!
    // if (m.maxZ() < 0) return;
    // </Stefano Mersi>
    double irrxy = 0;
    double irrPoint = 0;
    /*double irrPointCen = 0;
    double irrPoint11 = 0;
    double irrPoint12 = 0;
    double irrPoint21 = 0;
    double irrPoint22 = 0;*/

    auto vertex11 = std::make_pair(m.minZ(), m.minR());
    auto vertex12 = std::make_pair(m.maxZ(), m.minR());
    auto vertex21 = std::make_pair(m.minZ(), m.maxR());
    auto vertex22 = std::make_pair(m.maxZ(), m.maxR());
    
    XYZVector centerVector = m.center();
    auto center = std::make_pair(centerVector.Z(), centerVector.Rho());

    //if (centerVector.Z() < 0) return;
    double volume = 0.;
    for (const auto& s : m.sensors()) volume += s.sensorThickness() * m.area() / 1000.0; // volume is in cm^3

    //calculate irradiation in each vertex (and center) and take the worst
    irrPoint = irradiationMap_->calculateIrradiationPower(center);
    //irrPointCen = irrPoint;
    if(irrPoint > irrxy) irrxy = irrPoint;

    irrPoint = irradiationMap_->calculateIrradiationPower(vertex11);
    //irrPoint11 = irrPoint;
    if(irrPoint > irrxy) irrxy = irrPoint;

    irrPoint = irradiationMap_->calculateIrradiationPower(vertex12);
    //irrPoint12 = irrPoint;
    if(irrPoint > irrxy) irrxy = irrPoint;

    irrPoint = irradiationMap_->calculateIrradiationPower(vertex21);
    //irrPoint21 = irrPoint;
    if(irrPoint > irrxy) irrxy = irrPoint;

    irrPoint = irradiationMap_->calculateIrradiationPower(vertex22);
    //irrPoint22 = irrPoint;
    if(irrPoint > irrxy) irrxy = irrPoint;

    double fluence = irrxy * timeIntegratedLumi; // fluence is in 1MeV-equiv-neutrons/cm^2
    //double fluence = irrxy * timeIntegratedLumi * 1e15 * 80 * 1e-3; // fluence is in 1MeV-equiv-neutrons/cm^2

    // Calculate the leakage current at the reference temperature
    double leakageCurrentAtReferenceTemp = alphaParam * fluence * volume;
    // HAMBURG MODEL
    // Calculate the leakage current at operational temperature
    double leakageCurrent = leakageCurrentAtReferenceTemp * pow(operatingTemp / referenceTemp , 2) * exp(-insur::siliconEffectiveBandGap / (2 * insur::boltzmann_constant) * (1 / operatingTemp - 1 / referenceTemp));

    // Calculate the heat power produced on the silicon sensors due to leakage current
    double irradiatedPowerConsumption = leakageCurrent * chargeDepletionVoltage;

    // Store result
    m.irradiationPower(irradiatedPowerConsumption);
    TableRef tref = m.tableRef();
    irradiatedPowerConsumptionSummaries[tref.table].setCell(tref.row, tref.col, irradiatedPowerConsumption);

    //cout << "mod irr: " << cntName << "," << module->getLayer() << "," << module->getRing() << ";  " << module->getThickness() << "," << center.Rho() << ";  " << volume << "," << fluence << "," << leakCurrentScaled << "," << irradiatedPowerConsumption << endl;
  }
};


#endif

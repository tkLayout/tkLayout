#ifndef IRRADIATIONPOWER_H
#define IRRADIATIONPOWER_H

#include <string>
#include <map>
#include <vector>
#include <utility>

#include "Tracker.h"
#include "SimParms.h"
#include "IrradiationMapsManager.h"

#include "Visitor.h"
#include "SummaryTable.h"
#include "Units.h"

class IrradiationPowerVisitor : public GeometryVisitor {
  double numInvFemtobarns;
  double operatingTemp;
  double chargeDepletionVoltage;
  double alphaParam;
  double referenceTemp;
  const IrradiationMapsManager* irradiationMap_;
public:
  MultiSummaryTable irradiatedPowerConsumptionSummaries;

  void preVisit() {
    irradiatedPowerConsumptionSummaries.clear();   
  }

  void visit(SimParms& sp) {
    numInvFemtobarns = sp.timeIntegratedLumi();
    operatingTemp    = sp.operatingTemp();
    chargeDepletionVoltage = sp.chargeDepletionVoltage();
    alphaParam       = sp.alphaParm();
    referenceTemp    = sp.referenceTemp();
    irradiationMap_  = &sp.irradiationMapsManager();
  }

  void visit(Barrel& b) {
    irradiatedPowerConsumptionSummaries[b.myid()].setHeader("layer", "ring");
    irradiatedPowerConsumptionSummaries[b.myid()].setPrecision(3);        
  }

  void visit(Endcap& e) {
    irradiatedPowerConsumptionSummaries[e.myid()].setHeader("layer", "ring");
    irradiatedPowerConsumptionSummaries[e.myid()].setPrecision(3);        
  }

  void visit(DetectorModule& m) {
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

    double fluence = irrxy/(1./Units::cm2) * numInvFemtobarns; // fluence is printed in 1MeV-equiv-neutrons/cm^2
    //double fluence = irrxy * numInvFemtobarns * 1e15 * 80 * 1e-3; // fluence is in 1MeV-equiv-neutrons/cm^2
    double leakCurrentScaled = alphaParam * fluence * volume * pow((operatingTemp+273.15) / (referenceTemp+273.15), 2) * exp(-1.21/(2*8.617334e-5)*(1/(operatingTemp+273.15)-1/(referenceTemp+273.15))); 
    double irradiatedPowerConsumption = leakCurrentScaled * chargeDepletionVoltage;         
    //cout << "mod irr: " << cntName << "," << module->getLayer() << "," << module->getRing() << ";  " << module->getThickness() << "," << center.Rho() << ";  " << volume << "," << fluence << "," << leakCurrentScaled << "," << irradiatedPowerConsumption << endl;

//    modulePowerConsumptions_[&m] = irradiatedPowerConsumption; // CUIDADO CHECK WHERE IT IS NEEDED
    m.irradiationPower(irradiatedPowerConsumption);

    TableRef tref = m.tableRef();
    irradiatedPowerConsumptionSummaries[tref.table].setCell(tref.row, tref.col, irradiatedPowerConsumption);
  }
};


#endif

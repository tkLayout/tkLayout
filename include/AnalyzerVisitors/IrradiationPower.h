#ifndef IRRADIATIONPOWER_H
#define IRRADIATIONPOWER_H

#include <string>
#include <map>
#include <vector>
#include <utility>

#include "Tracker.h"
#include "SimParms.h"

#include "Visitor.h"
#include "SummaryTable.h"

class IrradiationPowerVisitor : public GeometryVisitor {
  double numInvFemtobarns;
  double operatingTemp;
  double chargeDepletionVoltage;
  double alphaParam;
  double referenceTemp;
  const IrradiationMap* irradiationMap_;
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
    irradiationMap_  = &sp.irradiationMap();
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
    XYZVector center = m.center();
    if (center.Z() < 0) return;
    double volume = 0.;
    for (const auto& s : m.sensors()) volume += s.sensorThickness() * m.area() / 1000.0; // volume is in cm^3
    double x  = center.Z()/25;
    double y  = center.Rho()/25;
    double x1 = floor(x);
    double x2 = ceil(x);
    double y1 = floor(y);
    double y2 = ceil(y);
    double irr11 = irradiationMap_->at(std::make_pair(int(x1), int(y1))); 
    double irr21 = irradiationMap_->at(std::make_pair(int(x2), int(y1)));
    double irr12 = irradiationMap_->at(std::make_pair(int(x1), int(y2)));
    double irr22 = irradiationMap_->at(std::make_pair(int(x2), int(y2)));
    double irrxy = irr11/((x2-x1)*(y2-y1))*(x2-x)*(y2-y) + irr21/((x2-x1)*(y2-y1))*(x-x1)*(y2-y) + irr12/((x2-x1)*(y2-y1))*(x2-x)*(y-y1) + irr22/((x2-x1)*(y2-y1))*(x-x1)*(y-y1); // bilinear interpolation
    double fluence = irrxy * numInvFemtobarns * 1e15 * 80 * 1e-3; // fluence is in 1MeV-equiv-neutrons/cm^2 
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

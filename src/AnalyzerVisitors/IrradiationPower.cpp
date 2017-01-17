#include "AnalyzerVisitors/IrradiationPower.h"

void IrradiationPowerVisitor::preVisit() {
  irradiatedPowerConsumptionSummaries.clear();   
}

void IrradiationPowerVisitor::visit(SimParms& sp) {
  timeIntegratedLumi = sp.timeIntegratedLumi();
  referenceTemp    = sp.referenceTemp() + insur::celcius_to_kelvin;   
  alphaParam       = sp.alphaParam();
  irradiationMap_  = &sp.irradiationMapsManager();
}

void IrradiationPowerVisitor::visit(Barrel& b) {
  irradiatedPowerConsumptionSummaries[b.myid()].setHeader("Layer", "Ring");
  irradiatedPowerConsumptionSummaries[b.myid()].setPrecision(3);        
}

void IrradiationPowerVisitor::visit(Endcap& e) {
  irradiatedPowerConsumptionSummaries[e.myid()].setHeader("Disk", "Ring");
  irradiatedPowerConsumptionSummaries[e.myid()].setPrecision(3);        
}

void IrradiationPowerVisitor::visit(DetectorModule& m) {
  operatingTemp    = m.operatingTemp() + insur::celcius_to_kelvin;
  chargeDepletionVoltage = m.chargeDepletionVoltage();
  double volume = 0.;
  std::vector<double> irradiationValues;
  
  for (const auto& s : m.sensors()) {
    // Calculate total volume occupied by sensors
    volume += s.sensorThickness() * m.area() / 1000.0; // volume in cm^3

    // Calculate irradiation at the center of the sensor
    const std::pair<double,double>& center = std::make_pair(s.center().Z(), s.center().Rho());
    double centerIrradiation = irradiationMap_->calculateIrradiationPower(center);
    irradiationValues.push_back(centerIrradiation);

    // Calculate irradiation at each vertex of the sensor
    for (int i = 0; i < s.envelopePoly().getNumSides(); i++) {
      const std::pair<double,double>& vertex = std::make_pair(s.envelopePoly().getVertex(i).Z(), s.envelopePoly().getVertex(i).Rho());
      double vertexIrradiation = irradiationMap_->calculateIrradiationPower(vertex);
      irradiationValues.push_back(vertexIrradiation);
    }
  }

  // THIS IS TO CALCULATE THE POWER DISSIPATED WITHIN THE SENSORS, DUE TO THE LEAKAGE CURRENT EFFECT

  // For a given module, take the worst fluence value obtained on the sensor(s) (FLUKA simulation)
  // This is a 1 MeV-neutrons-equivalent fluence, for an integrated luminosity = 1 fb-1
  double maxIrradiation = *max_element(irradiationValues.begin(), irradiationValues.end());  // 1MeV-equiv-neutrons / cm^2 / fb-1

  // 1 MeV-neutrons-equivalent fluence
  double fluence = maxIrradiation * timeIntegratedLumi; // 1MeV-equiv-neutrons / cm^2

  // Calculate the leakage current at the reference temperature
  double leakageCurrentAtReferenceTemp = alphaParam * fluence * volume;  // A
  // HAMBURG MODEL
  // Calculate the leakage current at operational temperature
  double leakageCurrent = leakageCurrentAtReferenceTemp * pow(operatingTemp / referenceTemp , 2) * exp(-insur::siliconEffectiveBandGap / (2 * insur::boltzmann_constant) * (1 / operatingTemp - 1 / referenceTemp)); // A

  // Calculate the heat power dissipated within the silicon sensors due to leakage current
  double irradiatedPowerConsumption = leakageCurrent * chargeDepletionVoltage; // W

  // Store result
  m.irradiationPower(irradiatedPowerConsumption);
  TableRef tref = m.tableRef();
  irradiatedPowerConsumptionSummaries[tref.table].setCell(tref.row, tref.col, irradiatedPowerConsumption);
}

#include "AnalyzerVisitors/IrradiationPower.h"

void IrradiationPowerVisitor::preVisit() {
  sensorsIrradiationPowerSummary.clear();   
}

void IrradiationPowerVisitor::visit(SimParms& sp) {
  timeIntegratedLumi_ = sp.timeIntegratedLumi();
  referenceTemp_    = sp.referenceTemp() + insur::celcius_to_kelvin;   
  alphaParam_       = sp.alphaParam();
  irradiationMap_  = &sp.irradiationMapsManager();
}

void IrradiationPowerVisitor::visit(Barrel& b) {
  sensorsIrradiationPowerSummary[b.myid()].setHeader("Layer", "Ring");
  sensorsIrradiationPowerSummary[b.myid()].setPrecision(3);        
}

void IrradiationPowerVisitor::visit(Endcap& e) {
  sensorsIrradiationPowerSummary[e.myid()].setHeader("Disk", "Ring");
  sensorsIrradiationPowerSummary[e.myid()].setPrecision(3);        
}

void IrradiationPowerVisitor::visit(DetectorModule& m) {
  operatingTemp_    = m.operatingTemp() + insur::celcius_to_kelvin;
  biasVoltage_ = m.biasVoltage();
  double volume = 0.;
  std::vector<double> irradiationValues;
  
  for (const auto& s : m.sensors()) {

    // Calculate total volume occupied by sensors
    volume += s.sensorThickness() * m.area() * Units::mm3 / Units::cm3; // convert volume to cm^3

    // Calculate irradiation at the center of the sensor
    const std::pair<double,double>& center = std::make_pair(s.center().Z(), s.center().Rho());
    double centerIrradiation = irradiationMap_->calculateIrradiationPower(center);
    irradiationValues.push_back(centerIrradiation);

    // Calculate irradiation at many vertexes of the sensor
    for (int i = 0; i < s.envelopePoly().getNumSides(); i++) {
      // each vertex
      std::pair<double,double> vertex = std::make_pair(s.envelopePoly().getVertex(i).Z(), s.envelopePoly().getVertex(i).Rho());
      double vertexIrradiation = irradiationMap_->calculateIrradiationPower(vertex);
      irradiationValues.push_back(vertexIrradiation);

      // each middle of 2 consecutive vertexes
      std::pair<double,double> midVertex = std::make_pair(s.envelopeMidPoly().getVertex(i).Z(), s.envelopeMidPoly().getVertex(i).Rho());
      double midVertexIrradiation = irradiationMap_->calculateIrradiationPower(midVertex);
      irradiationValues.push_back(midVertexIrradiation);
    }
  }

  // For a given module, consider all fluence values obtained on the sensor(s) (FLUKA simulation)
  // This is a 1 MeV-neutrons-equivalent fluence, for an integrated luminosity = 1 fb-1
  double sum = std::accumulate(irradiationValues.begin(), irradiationValues.end(), 0.);
  double irradiationMean = sum / irradiationValues.size();                                   // 1MeV-equiv-neutrons / cm^2 / fb-1
  double irradiationMax = *max_element(irradiationValues.begin(), irradiationValues.end());  // 1MeV-equiv-neutrons / cm^2 / fb-1


  // THIS IS TO CALCULATE THE POWER DISSIPATED WITHIN THE SENSORS, DUE TO THE LEAKAGE CURRENT EFFECT
  const double sensorsIrradiationPowerMean = computeSensorsIrradiationPower(irradiationMean, timeIntegratedLumi_, alphaParam_,
								     volume, referenceTemp_, operatingTemp_, biasVoltage_);
  const double sensorsIrradiationPowerMax = computeSensorsIrradiationPower(irradiationMax, timeIntegratedLumi_, alphaParam_,
								      volume, referenceTemp_, operatingTemp_, biasVoltage_);

  // Store results
  m.sensorsIrradiationPowerMean(sensorsIrradiationPowerMean);
  m.sensorsIrradiationPowerMax(sensorsIrradiationPowerMax);

  TableRef tref = m.tableRef();
  std::ostringstream values;
  values.str("");
  values << std::dec << std::fixed << std::setprecision(3) << sensorsIrradiationPowerMean << "," << sensorsIrradiationPowerMax;
  sensorsIrradiationPowerSummary[tref.table].setCell(tref.row, tref.col, values.str());
}


const double IrradiationPowerVisitor::computeSensorsIrradiationPower(const double& irradiation, const double& timeIntegratedLumi,
							       const double& alphaParam, const double& volume, const double& referenceTemp, 
							       const double& operatingTemp, const double& biasVoltage) const {			
  // 1 MeV-neutrons-equivalent fluence
  double fluence = irradiation * timeIntegratedLumi; // 1MeV-equiv-neutrons / cm^2

  // Calculate the leakage current at the reference temperature
  double leakageCurrentAtReferenceTemp = alphaParam * fluence * volume;  // A
  // HAMBURG MODEL
  // Calculate the leakage current at operational temperature
  double leakageCurrent = leakageCurrentAtReferenceTemp * pow(operatingTemp / referenceTemp , 2) * exp(-insur::siliconEffectiveBandGap / (2 * insur::boltzmann_constant) * (1 / operatingTemp - 1 / referenceTemp)); // A

  // Calculate the heat power dissipated within the silicon sensors due to leakage current
  double sensorsIrradiationPower = leakageCurrent * biasVoltage; // W

  return sensorsIrradiationPower;
}

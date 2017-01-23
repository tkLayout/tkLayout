#include "AnalyzerVisitors/IrradiationPower.h"

void IrradiationPowerVisitor::preVisit() {
  sensorsIrradiationPowerSummary.clear();   
}

void IrradiationPowerVisitor::visit(SimParms& sp) {
  timeIntegratedLumi_ = sp.timeIntegratedLumi();
  referenceTemp_    = sp.referenceTemp() + insur::celsius_to_kelvin;   
  alphaParam_       = sp.alphaParam();
  irradiationMap_  = &sp.irradiationMapsManager();
}

void IrradiationPowerVisitor::visit(Barrel& b) {
  isBarrel_ = true;     
}

void IrradiationPowerVisitor::visit(RodPair& r) {
  isOuterRadiusRod_ = r.isOuterRadiusRod();     
}

void IrradiationPowerVisitor::visit(Endcap& e) {
  isBarrel_ = false;
  isOuterRadiusRod_ = false;    
}

void IrradiationPowerVisitor::visit(DetectorModule& m) {
  operatingTemp_    = m.operatingTemp() + insur::celsius_to_kelvin;
  biasVoltage_ = m.biasVoltage();
  double volume = m.totalSensorsVolume() * Units::mm3 / Units::cm3; // Total volume occupied by sensors, converted to cm^3
  
  std::vector<double> irradiationValues;
  for (const auto& s : m.sensors()) {

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

  // STORE RESULTS
  // Results for each module
  m.sensorsIrradiationPowerMean(sensorsIrradiationPowerMean);
  m.sensorsIrradiationPowerMax(sensorsIrradiationPowerMax);

  // Also gather results for all modules of a given type, identified by ModuleRef.
  // This will be used for summary tables.
  TableRef tableRef = m.tableRef();
  ModuleRef moduleRef = std::make_tuple(isBarrel_, isOuterRadiusRod_, tableRef.table, tableRef.row, tableRef.col);
  // mean
  sensorsIrradiationPowerMean_[moduleRef] += sensorsIrradiationPowerMean;
  // max
  if (sensorsIrradiationPowerMax > sensorsIrradiationPowerMax_[moduleRef]) sensorsIrradiationPowerMax_[moduleRef] = sensorsIrradiationPowerMax;
  // counter
  modulesCounter_[moduleRef]++;
}

void IrradiationPowerVisitor::postVisit() {
  // Create summary tables of results.
  // All results are per module category, identified by ModuleRef.
  for (const auto& it : modulesCounter_) {  // for each module category
    if (it.second > 0) {
      // Identifiers of the module category.
      std::string name = std::get<2>(it.first);
      bool isBarrel = std::get<0>(it.first);
      bool isOuterRadiusRod = std::get<1>(it.first);
      if (isBarrel) {
	if (isOuterRadiusRod) name += ", outer radius rods";
	else name += ", inner radius rods";
      }
      int row = std::get<3>(it.first);
      int col = std::get<4>(it.first);

      // Obtain the mean and the max sensorsIrradiationPower, for a given module category.
      double sensorsIrradiationPowerMean = sensorsIrradiationPowerMean_[it.first] / it.second;
      double sensorsIrradiationPowerMax = sensorsIrradiationPowerMax_[it.first];
      std::ostringstream values;
      values.str("");
      values << std::dec << std::fixed << std::setprecision(3) << sensorsIrradiationPowerMean << "," << sensorsIrradiationPowerMax;

      // Store results in summary table.
      if (isBarrel) sensorsIrradiationPowerSummary[name].setHeader("Layer", "Ring");
      else sensorsIrradiationPowerSummary[name].setHeader("Disk", "Ring");
      sensorsIrradiationPowerSummary[name].setPrecision(3);   
      sensorsIrradiationPowerSummary[name].setCell(row, col, values.str());
    }
    else logERROR("Tried to access values from sensorsIrradiationPowerMean_, sensorsIrradiationPowerMax_, but no module in modulesCounter_.");
  }  
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

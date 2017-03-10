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
  isOuterRadiusRod_ = false; // false by default, since no rod here.   
}

void IrradiationPowerVisitor::visit(DetectorModule& m) {
  operatingTemp_    = m.operatingTemp() + insur::celsius_to_kelvin;
  biasVoltage_ = m.biasVoltage();
  double volume = m.totalSensorsVolume() * Units::mm3 / Units::cm3; // Total volume occupied by sensors, converted to cm^3
  
  // A) FOR A GIVEN MODULE, GET THE FLUENCE VALUES (IRRADIATIONS) ON ITS SENSOR(S), FROM THE irradiationmap_.
  // The irradiationMap_ was obtained with FLUKA simulation.
  // The irradiationMap_ values are 1 MeV-neutrons-equivalent fluence, for an integrated luminosity = 1 fb-1 .
  // The values are the mean and the max on different points on the module's sensor(s).
  std::pair<double, double> irradiationMeanMax = getModuleIrradiationMeanMax(irradiationMap_, m);
  double irradiationMean = irradiationMeanMax.first * timeIntegratedLumi_;  // 1MeV-equiv-neutrons / cm^2
  double irradiationMax = irradiationMeanMax.second * timeIntegratedLumi_;  // 1MeV-equiv-neutrons / cm^2


  // B) FOR A GIVEN MODULE, CALCULATE THE POWER DISSIPATED WITHIN THE SENSORS, DUE TO THE LEAKAGE CURRENT EFFECT
  // This use the irradiation on sensors from FLUKA maps, which has just been obtained : irradiationMean, irradiationMax.
  // Many other parameters are used (see below).
  const double sensorsIrradiationPowerMean = computeSensorsIrradiationPower(irradiationMean, alphaParam_,
									    volume, referenceTemp_, operatingTemp_, biasVoltage_);
  const double sensorsIrradiationPowerMax = computeSensorsIrradiationPower(irradiationMax, alphaParam_,
									   volume, referenceTemp_, operatingTemp_, biasVoltage_);

  // C) STORE RESULTS
  // Results for each module
  m.sensorsIrradiationMean(irradiationMean); // 1MeV-equiv-neutrons / cm^2
  m.sensorsIrradiationMax(irradiationMax);   // 1MeV-equiv-neutrons / cm^2
  m.sensorsIrradiationPowerMean(sensorsIrradiationPowerMean);  // W
  m.sensorsIrradiationPowerMax(sensorsIrradiationPowerMax);    // W

  // Also gather results for all modules of a given type, identified by ModuleRef.
  // This will be used for summary tables.
  TableRef tableRef = m.tableRef();
  ModuleRef moduleRef = std::make_tuple(isBarrel_, isOuterRadiusRod_, tableRef.table, tableRef.row, tableRef.col);
  // mean
  sensorsIrradiationPowerMean_[moduleRef] += sensorsIrradiationPowerMean;
  sensorsIrradiationMean_[moduleRef] += irradiationMean;
  // max
  sensorsIrradiationPowerMax_[moduleRef] = MAX(sensorsIrradiationPowerMax_[moduleRef], sensorsIrradiationPowerMax);
  sensorsIrradiationMax_[moduleRef] = MAX(sensorsIrradiationMax_[moduleRef], irradiationMax);
  // counter
  modulesCounter_[moduleRef]++;
}


void IrradiationPowerVisitor::postVisit() {
  // Create summary tables of results. All tables are displayed on website.
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

      // Obtain the mean and the max sensorsIrradiationPower and sensorsIrradiation for a given module category.
      double sensorsIrradiationPowerMean = sensorsIrradiationPowerMean_[it.first] / it.second;
      double sensorsIrradiationPowerMax = sensorsIrradiationPowerMax_[it.first];
      double sensorsIrradiationMean = sensorsIrradiationMean_[it.first] / it.second;
      double sensorsIrradiationMax = sensorsIrradiationMax_[it.first];
      std::ostringstream powerValues, irradiationValues;
      powerValues.str("");
      irradiationValues.str("");
      powerValues << std::dec << std::fixed << std::setprecision(3) << sensorsIrradiationPowerMean << "," << sensorsIrradiationPowerMax;
      //irradiationValues << std::dec << std::scientific << std::setprecision(3) << sensorsIrradiationMean << "," << sensorsIrradiationMax;
      irradiationValues << std::dec << std::scientific << std::setprecision(2) << sensorsIrradiationMax;

      // Store results in the power and irradiation summary tables
      if (isBarrel) {
	sensorsIrradiationPowerSummary[name].setHeader("Layer", "Ring");
	sensorsIrradiationSummary[name].setHeader("Layer", "Ring");
      }
      else {
	sensorsIrradiationPowerSummary[name].setHeader("Disk", "Ring");
	sensorsIrradiationSummary[name].setHeader("Disk", "Ring");
      }
      sensorsIrradiationPowerSummary[name].setPrecision(3);   
      sensorsIrradiationPowerSummary[name].setCell(row, col, powerValues.str());
      sensorsIrradiationSummary[name].setPrecision(3);   
      sensorsIrradiationSummary[name].setCell(row, col, irradiationValues.str());
    }
    else logERROR("Tried to access values from sensorsIrradiationPowerMean_, sensorsIrradiationPowerMax_, but no module in modulesCounter_.");
  }  
}


/** 
    Get, for a given module, the fluence values (irradiations) on its sensor(s), from an irradiationmap.
    Several points on the module's sensor(s) are considered.
    @return : pair (mean, max) of the irradiation values on the different points.
*/
std::pair<double, double> IrradiationPowerVisitor::getModuleIrradiationMeanMax(const IrradiationMapsManager* irradiationMap, const DetectorModule& m) {

  // Get the irradiation from irradiationMap at differents points of the modules's sensor(s).
  std::vector<double> irradiationValues;
  for (const auto& s : m.sensors()) {

    // Center of the sensor
    const std::pair<double,double>& center = std::make_pair(s.center().Z(), s.center().Rho());
    double centerIrradiation = irradiationMap->calculateIrradiationPower(center);
    irradiationValues.push_back(centerIrradiation);

    // Many vertexes of the sensor
    for (int i = 0; i < s.envelopePoly().getNumSides(); i++) {
      // each vertex
      std::pair<double,double> vertex = std::make_pair(s.envelopePoly().getVertex(i).Z(), s.envelopePoly().getVertex(i).Rho());
      double vertexIrradiation = irradiationMap->calculateIrradiationPower(vertex);
      irradiationValues.push_back(vertexIrradiation);

      // each middle of 2 consecutive vertexes
      std::pair<double,double> midVertex = std::make_pair(s.envelopeMidPoly().getVertex(i).Z(), s.envelopeMidPoly().getVertex(i).Rho());
      double midVertexIrradiation = irradiationMap->calculateIrradiationPower(midVertex);
      irradiationValues.push_back(midVertexIrradiation);
    }
  }

  // For a given module, take the average and the max irradiation on all the considered points.
  double sum = std::accumulate(irradiationValues.begin(), irradiationValues.end(), 0.);
  double irradiationMean = sum / irradiationValues.size();                                   // 1MeV-equiv-neutrons / cm^2 / fb-1
  double irradiationMax = *max_element(irradiationValues.begin(), irradiationValues.end());  // 1MeV-equiv-neutrons / cm^2 / fb-1

  std::pair<double, double> irradiationMeanMax = std::make_pair(irradiationMean, irradiationMax);
  return irradiationMeanMax;
}


/** 
    Calculate, for a given module, the power dissipated within the sensors, due to the leakage current effect.
    * @param totalFluence : time-integrated fluence on the module's sensors (1MeV-equiv-neutrons / cm^2)
    * @param alphaParam : radiation-damage constant (A.cm-1)
    * @param volume : total volume occupied by the module's sensor(s) (cm^3)
    * @param referenceTemp : temperature of reference (°C)
    * @param operatingTemp : operating temperature (°C)
    * @param biasVoltage : bias voltage (V)
    */
const double IrradiationPowerVisitor::computeSensorsIrradiationPower(const double& totalFluence, 
								     const double& alphaParam, const double& volume, const double& referenceTemp, 
								     const double& operatingTemp, const double& biasVoltage) const {			
  // Calculate the leakage current at the reference temperature
  double leakageCurrentAtReferenceTemp = alphaParam * totalFluence * volume;  // A
  // HAMBURG MODEL
  // Calculate the leakage current at operational temperature
  double leakageCurrent = leakageCurrentAtReferenceTemp * pow(operatingTemp / referenceTemp , 2) * exp(-insur::siliconEffectiveBandGap / (2 * insur::boltzmann_constant) * (1 / operatingTemp - 1 / referenceTemp)); // A

  // Calculate the heat power dissipated within the silicon sensors due to leakage current
  double sensorsIrradiationPower = leakageCurrent * biasVoltage; // W

  return sensorsIrradiationPower;
}

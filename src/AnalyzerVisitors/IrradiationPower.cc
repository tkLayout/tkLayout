#include "AnalyzerVisitors/IrradiationPower.hh"

void IrradiationPowerVisitor::preVisit() {
  sensorsPowerSummary.clear();
  sensorsFluenceSummary.clear();
  sensorsFluencePerType.clear();
  sensorsDoseSummary.clear();
  sensorsDosePerType.clear();
  lumiInformation="";
  mapInformation="";
}

void IrradiationPowerVisitor::visit(SimParms& sp) {
  timeIntegratedLumi_ = sp.timeIntegratedLumi();
  referenceTemp_    = sp.referenceTemp() + insur::celsius_to_kelvin;
  alphaParam_       = sp.alphaParam();
  irradiationMap_  = &sp.irradiationMapsManager();
  doseMap_         = &sp.doseMapsManager();
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
  std::string modType = m.summaryFullType();
  operatingTemp_    = m.operatingTemp() + insur::celsius_to_kelvin;
  biasVoltage_ = m.biasVoltage();
  double volume = m.totalSensorsVolume() * Units::mm3 / Units::cm3; // Total volume occupied by sensors, converted to cm^3

  // A) FOR A GIVEN MODULE, GET THE FLUENCE VALUES (IRRADIATIONS) ON ITS SENSOR(S), FROM THE irradiationmap_.
  // The irradiationMap_ was obtained with FLUKA simulation.
  // The irradiationMap_ values are 1 MeV-neutrons-equivalent fluence, for an integrated luminosity = 1 fb-1 .
  // The values are the mean and the max on different points on the module's sensor(s).
  std::pair<double, double> irradiationMeanMax = getModuleFluenceMeanMax(irradiationMap_, m);
  double irradiationMean = irradiationMeanMax.first * timeIntegratedLumi_;  // 1MeV-equiv-neutrons / cm^2
  double irradiationMax = irradiationMeanMax.second * timeIntegratedLumi_;  // 1MeV-equiv-neutrons / cm^2
  //Also get the dose
  std::pair<double, double> doseMeanMax = getModuleDoseMeanMax(doseMap_,m);
  double doseMean = doseMeanMax.first * timeIntegratedLumi_;
  double doseMax = doseMeanMax.second * timeIntegratedLumi_;
  if(lumiInformation==""){
    lumiInformation+=std::to_string(timeIntegratedLumi_);
    lumiInformation+=" invfb";
  }

  if(mapInformation==""){
    for(auto irradMaps : irradiationMap_->getIrradiationMaps()){
      int theIdx = irradMaps.getMapName().rfind('/');
      mapInformation+=(irradMaps.getMapName()).substr(theIdx+1);
      mapInformation+=", ";
    }
    for(auto irradMaps : doseMap_->getIrradiationMaps()){
      int theIdx = irradMaps.getMapName().rfind('/');
      mapInformation+=(irradMaps.getMapName()).substr(theIdx+1);
      mapInformation+=", ";
    }
  }

  // B) FOR A GIVEN MODULE, CALCULATE THE POWER DISSIPATED WITHIN THE SENSORS, DUE TO THE LEAKAGE CURRENT EFFECT
  // This use the irradiation on sensors from FLUKA maps, which has just been obtained : irradiationMean, irradiationMax.
  // Many other parameters are used (see below).
  const double sensorsPowerMean = computeSensorsPower(irradiationMean, alphaParam_,
									    volume, referenceTemp_, operatingTemp_, biasVoltage_);
  const double sensorsPowerMax = computeSensorsPower(irradiationMax, alphaParam_,
									   volume, referenceTemp_, operatingTemp_, biasVoltage_);

  // C) STORE RESULTS
  // Results for each module
  //
  m.sensorsIrradiationMean(irradiationMean); // 1MeV-equiv-neutrons / cm^2
  m.sensorsIrradiationMax(irradiationMax);   // 1MeV-equiv-neutrons / cm^2
  m.sensorsDoseMean(doseMean); // 1MeV-equiv-neutrons / cm^2
  m.sensorsDoseMax(doseMax);   // 1MeV-equiv-neutrons / cm^2
  m.sensorsIrradiationPowerMean(sensorsPowerMean);  // W
  m.sensorsIrradiationPowerMax(sensorsPowerMax);    // W
  // Also gather results for all modules of a given type, identified by ModuleRef.
  // This will be used for summary tables.
  TableRef tableRef = m.tableRef();
  ModuleRef moduleRef = std::make_tuple(isBarrel_, isOuterRadiusRod_, tableRef.table, tableRef.row, tableRef.col);
  // mean
  sensorsPowerMean_[moduleRef] += sensorsPowerMean;
  sensorsFluenceMean_[moduleRef] += irradiationMean;
  sensorsDoseMean_[moduleRef] += doseMean;
  // max
  sensorsPowerMax_[moduleRef] = MAX(sensorsPowerMax_[moduleRef], sensorsPowerMax);
  sensorsFluenceMax_[moduleRef] = MAX(sensorsFluenceMax_[moduleRef], irradiationMax);
  sensorsDoseMax_[moduleRef] = MAX(sensorsDoseMax_[moduleRef], doseMax);
  // counter
  modulesCounter_[moduleRef]++;
  // The list of modules per irradiation type
  mapTypeToFluence_[modType].push_back(&m);
}


void IrradiationPowerVisitor::postVisit() {
  // Sort the irradiation and extract the max and 95% percentiles
  sensorsFluencePerType.setCell(0, 0, "Type");
  sensorsFluencePerType.setCell(0, 1, "# mods");
  sensorsFluencePerType.setCell(0, 2, "Max. module irrad (avg) [Hb]");
  sensorsFluencePerType.setCell(0, 3, "95% module irrad (avg) [Hb]");
  sensorsFluencePerType.setCell(0, 4, "z_max [mm]");
  sensorsFluencePerType.setCell(0, 5, "r_max [mm]");
  sensorsDosePerType.setCell(0, 0, "Type");
  sensorsDosePerType.setCell(0, 1, "# mods");
  sensorsDosePerType.setCell(0, 2, "Max. module dose (avg) [Gy]");
  sensorsDosePerType.setCell(0, 3, "95% module dose (avg) [Gy]");
  sensorsDosePerType.setCell(0, 4, "z_max [mm]");
  sensorsDosePerType.setCell(0, 5, "r_max [mm]");
  int iRow=0;
  for (auto& it : mapTypeToFluence_ ) {
    iRow++;
    const std::string& typeName = it.first;
    std::vector<const DetectorModule*>& irrads = it.second;
    std::sort(irrads.begin(), irrads.end(),
	      [] (const DetectorModule* a, const DetectorModule* b) {
		return a->sensorsIrradiationMean() < b->sensorsIrradiationMean();
	      });
    int nModules = irrads.size();
    std::ostringstream irrad_Max("");
    std::ostringstream irrad_95perc("");
    std::ostringstream max_z("");
    std::ostringstream max_r("");
    auto& hottestModule =  irrads.at(nModules-1);
    auto& hottest95Module = irrads.at(ceil(double(nModules)*95/100-1));
    irrad_Max    << std::dec << std::scientific << std::setprecision(2) << hottestModule->sensorsIrradiationMean();
    irrad_95perc << std::dec << std::scientific << std::setprecision(2) <<  hottest95Module->sensorsIrradiationMean();
    max_z << std::dec << std::fixed << std::setprecision(2) << hottestModule->center().Z();
    max_r << std::dec << std::fixed << std::setprecision(2) << hottestModule->center().Rho();
    sensorsFluencePerType.setCell(iRow, 0, typeName);
    sensorsFluencePerType.setCell(iRow, 1, nModules);
    sensorsFluencePerType.setCell(iRow, 2, irrad_Max.str());
    sensorsFluencePerType.setCell(iRow, 3, irrad_95perc.str());
    sensorsFluencePerType.setCell(iRow, 4, max_z.str());
    sensorsFluencePerType.setCell(iRow, 5, max_r.str());
  }
  iRow=0;
  for (auto& it : mapTypeToFluence_ ) {
    iRow++;
    const std::string& typeName = it.first;
    std::vector<const DetectorModule*>& irrads = it.second;
    std::sort(irrads.begin(), irrads.end(),
	      [] (const DetectorModule* a, const DetectorModule* b) {
		return a->sensorsDoseMean() < b->sensorsDoseMean();
	      });
    int nModules = irrads.size();
    std::ostringstream irrad_Max("");
    std::ostringstream irrad_95perc("");
    std::ostringstream max_z("");
    std::ostringstream max_r("");
    auto& hottestModule =  irrads.at(nModules-1);
    auto& hottest95Module = irrads.at(ceil(double(nModules)*95/100-1));
    irrad_Max    << std::dec << std::scientific << std::setprecision(2) << hottestModule->sensorsDoseMean();
    irrad_95perc << std::dec << std::scientific << std::setprecision(2) <<  hottest95Module->sensorsDoseMean();
    max_z << std::dec << std::fixed << std::setprecision(2) << hottestModule->center().Z();
    max_r << std::dec << std::fixed << std::setprecision(2) << hottestModule->center().Rho();
    sensorsDosePerType.setCell(iRow, 0, typeName);
    sensorsDosePerType.setCell(iRow, 1, nModules);
    sensorsDosePerType.setCell(iRow, 2, irrad_Max.str());
    sensorsDosePerType.setCell(iRow, 3, irrad_95perc.str());
    sensorsDosePerType.setCell(iRow, 4, max_z.str());
    sensorsDosePerType.setCell(iRow, 5, max_r.str());
  }
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

      // Obtain the mean and the max sensorsPower and sensorsFluence for a given module category.
      double sensorsPowerMean = sensorsPowerMean_[it.first] / it.second;
      double sensorsFluenceMean = sensorsFluenceMean_[it.first] / it.second;
      double sensorsDoseMean = sensorsDoseMean_[it.first] / it.second;
      std::ostringstream powerValues, irradiationValues, doseValues;
      powerValues.str("");
      irradiationValues.str("");
      doseValues.str("");
      powerValues << std::dec << std::fixed << std::setprecision(3) << sensorsPowerMean;
      irradiationValues << std::dec << std::scientific << std::setprecision(2) << sensorsFluenceMean;
      doseValues << std::dec << std::scientific << std::setprecision(2) << sensorsDoseMean;

      // Store results in the power and irradiation summary tables
      if (isBarrel) {
	sensorsPowerSummary[name].setHeader("Layer", "Ring");
	sensorsFluenceSummary[name].setHeader("Layer", "Ring");
        sensorsDoseSummary[name].setHeader("Layer","Ring");
      }
      else {
	sensorsPowerSummary[name].setHeader("Disk", "Ring");
	sensorsFluenceSummary[name].setHeader("Disk", "Ring");
        sensorsDoseSummary[name].setHeader("Disk","Ring");
      }
      sensorsPowerSummary[name].setPrecision(3);
      sensorsPowerSummary[name].setCell(row, col, powerValues.str());
      sensorsFluenceSummary[name].setPrecision(3);
      sensorsFluenceSummary[name].setCell(row, col, irradiationValues.str());
      sensorsDoseSummary[name].setPrecision(3);
      sensorsDoseSummary[name].setCell(row, col, doseValues.str());
    }
    else logERROR("Tried to access values from sensorsPowerMean_, sensorsPowerMax_, but no module in modulesCounter_.");
  }
}


/**
    Get, for a given module, the fluence values (irradiations) on its sensor(s), from an irradiationmap.
    Several points on the module's sensor(s) are considered.
    @return : pair (mean, max) of the irradiation values on the different points.
*/
std::pair<double, double> IrradiationPowerVisitor::getModuleFluenceMeanMax(const IrradiationMapsManager* irradiationMap, const DetectorModule& m) {

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
    Get, for a given module, the dose values on its sensor(s), from an irradiationmap.
    Several points on the module's sensor(s) are considered.
    @return : pair (mean, max) of the doses at the different points.
*/
std::pair<double, double> IrradiationPowerVisitor::getModuleDoseMeanMax(const IrradiationMapsManager* doseMap, const DetectorModule& m) {

  // Get the dose from doseMap at differents points of the modules's sensor(s).
  std::vector<double> doseValues;
  for (const auto& s : m.sensors()) {

    // Center of the sensor
    const std::pair<double,double>& center = std::make_pair(s.center().Z(), s.center().Rho());
    double centerDose = doseMap->calculateIrradiationPower(center);
    doseValues.push_back(centerDose);

    // Many vertexes of the sensor
    for (int i = 0; i < s.envelopePoly().getNumSides(); i++) {
      // each vertex
      std::pair<double,double> vertex = std::make_pair(s.envelopePoly().getVertex(i).Z(), s.envelopePoly().getVertex(i).Rho());
      double vertexDose = doseMap->calculateIrradiationPower(vertex);
      doseValues.push_back(vertexDose);

      // each middle of 2 consecutive vertexes
      std::pair<double,double> midVertex = std::make_pair(s.envelopeMidPoly().getVertex(i).Z(), s.envelopeMidPoly().getVertex(i).Rho());
      double midVertexDose = doseMap->calculateIrradiationPower(midVertex);
      doseValues.push_back(midVertexDose);
    }
  }

  // For a given module, take the average and the max irradiation on all the considered points.
  double sum = std::accumulate(doseValues.begin(), doseValues.end(), 0.);
  double doseMean = sum / doseValues.size();                                   // Gy / fb-1
  double doseMax = *max_element(doseValues.begin(), doseValues.end());  // Gy / fb-1

  std::pair<double, double> doseMeanMax = std::make_pair(doseMean, doseMax);
  return doseMeanMax;
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
const double IrradiationPowerVisitor::computeSensorsPower(const double& totalFluence,
								     const double& alphaParam, const double& volume, const double& referenceTemp,
								     const double& operatingTemp, const double& biasVoltage) const {			
  // Calculate the leakage current at the reference temperature
  double leakageCurrentAtReferenceTemp = alphaParam * totalFluence * volume;  // A
  // HAMBURG MODEL
  // Calculate the leakage current at operational temperature
  double leakageCurrent = leakageCurrentAtReferenceTemp * pow(operatingTemp / referenceTemp , 2) * exp(-insur::siliconEffectiveBandGap / (2 * insur::boltzmann_constant) * (1 / operatingTemp - 1 / referenceTemp)); // A

  // Calculate the heat power dissipated within the silicon sensors due to leakage current
  double sensorsPower = leakageCurrent * biasVoltage; // W

  return sensorsPower;
}

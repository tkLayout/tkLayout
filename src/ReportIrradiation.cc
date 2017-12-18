#include <ReportIrradiation.hh>
#include <PlotDrawer.hh>

void ReportIrradiation::analyze() {
  computeIrradiationPowerConsumption();
  computeChipPowerConsumptionTable();
  preparePowerHistograms();
}

void ReportIrradiation::preparePowerHistograms() {
  prepareTrackerMap(sensorsIrradiationPowerMap, "sensorsIrradiationPowerMap", "Map of power dissipation in sensors (after irradiation)");
  prepareTrackerMap(totalPowerConsumptionMap, "irradiationPowerConsumptionMap", "Map of power dissipation in modules (after irradiation)");
}

void ReportIrradiation::computeIrradiationPowerConsumption() {
  IrradiationPowerVisitor irradiation_;
  irradiation_.preVisit();
  SimParms::getInstance().accept(irradiation_);
  tracker.accept(irradiation_);
  irradiation_.postVisit();

  powerSummaries = irradiation_.sensorsPowerSummary;
  fluenceSummaries = irradiation_.sensorsFluenceSummary;
  fluenceSummaryPerType = irradiation_.sensorsFluencePerType;
}

void ReportIrradiation::computeChipPowerConsumptionTable() {
  int iType=1;
  struct ModuleTypeVisitor : public ConstGeometryVisitor {
    std::map<std::string, std::vector<std::string> > typeMap;
    void visit(const Module& m) {
      if (!typeMap.count(m.moduleType())) {
	std::vector<std::string>& displayValues = typeMap[m.moduleType()];
	displayValues.push_back(any2str(m.totalPower(), 0));
	displayValues.push_back(any2str(m.powerPerModule(), 0));
	for (const auto& sen : m.sensors()) {
	  displayValues.push_back(any2str(sen.powerPerChannel(), 0));	  
	}
      }
    }
  };
  ModuleTypeVisitor v;
  tracker.accept(v);
  for (auto it = v.typeMap.begin(); it != v.typeMap.end(); ++it) {
    chipPowerPerType.setCell(0, iType,it->first);
    int j=1;
    for (auto& itPowerValue : it->second) {
      if (j==1) chipPowerPerType.setCell(j, 0, "Totap power [mW]");
      else if (j==2) chipPowerPerType.setCell(j, 0, "Power per module [mW]");
      else chipPowerPerType.setCell(j, 0, Form("Power per channel in sensor %d [mW]", j-2));
      chipPowerPerType.setCell(j, iType, itPowerValue);
      j++;
    }
    iType++;
  }
}

std::string ReportIrradiation::createSensorsIrradiationCsv() {
  class TrackerVisitor : public ConstGeometryVisitor {
    std::stringstream output_;
    string sectionName_;
    int layerId_;
    bool isOuterRadiusRod_;
  public:
    void preVisit() {
      output_ << "Section, Layer, Ring, moduleType, dsDistance, isOuterRadiusRod_bool, operatingTemperature_Celsius, biasVoltage_V, meanWidth_mm, length_mm, sensorThickness_mm, sensor(s)Volume(totalPerModule)_mm3, sensorsPowerMean_W, sensorsPowerMax_W, sensorsFluenceMean_Hb, sensorsFluenceMax_Hb" << std::endl;
    }
    void visit(const Barrel& b) { sectionName_ = b.myid(); }
    void visit(const Endcap& e) { sectionName_ = e.myid(); }
    void visit(const Layer& l)  { layerId_ = l.myid(); }
    void visit(const RodPair& r)  { isOuterRadiusRod_ = r.isOuterRadiusRod(); }
    void visit(const Disk& d)  { isOuterRadiusRod_ = false; layerId_ = d.myid(); } // no rod here !
    void visit(const Module& m) {
      output_ << sectionName_ << ", "
	      << layerId_ << ", "
	      << m.moduleRing() << ", "
	      << m.moduleType() << ", "
	      << m.dsDistance() << ", "
	      << isOuterRadiusRod_ << ", "
	      << std::fixed << std::setprecision(6)
	      << m.operatingTemp() << ", "
	      << m.biasVoltage() << ", "
	      << m.meanWidth() << ", "
	      << m.length() << ", "
	      << m.sensorThickness() << ", "
	      << m.totalSensorsVolume() << ", "
	      << std::fixed << std::setprecision(3)
	      << m.sensorsIrradiationPowerMean() << ", "
	      << m.sensorsIrradiationPowerMax() << ", "
	      << m.sensorsIrradiationMean() << ", "
	      << m.sensorsIrradiationMax()	
	      << std::endl;
    }

    std::string output() const { return output_.str(); }
  };

  TrackerVisitor v;
  v.preVisit();
  tracker.accept(v);
  return v.output();
}


void ReportIrradiation::dumpRadiationTableSummary(RootWPage& myPage, std::map<std::string, SummaryTable>& radiationSummaries,
						  const std::string& title, std::string units) {
  for (std::map<std::string, SummaryTable>::iterator it = radiationSummaries.begin(); it != radiationSummaries.end(); ++it) {
    RootWContent& myContent = myPage.addContent(title+std::string(" (") + it->first + ")", false);
    RootWTable* comments = new RootWTable();
    comments->setContent(0, 0, "Values in table represent the mean value per module with fully irradiated sensors ["+units+"]");
    myContent.addItem(comments);
    myContent.addTable().setContent(it->second.getContent());
  }    
} 

void ReportIrradiation::visualizeTo(RootWSite& site) {
  std::string trackerName = tracker.myid();
  std::string pageName = "Irradiation (" + trackerName + ")";
  std::string pageAddress = "irradiation_" + trackerName + ".html";

  RootWPage& myPage = site.addPage(pageName);
  myPage.setAddress(pageAddress);

  // Irradiation on each module type (fine grained)
  RootWContent& summaryContent = myPage.addContent("Irradiation summary per module type");
  RootWTable& summaryTable = summaryContent.addTable();
  summaryTable.setContent(fluenceSummaryPerType.getContent());

  // Power consumption per module tpye (coarse grained)
  RootWContent& chipContent = myPage.addContent("Chip power consumption per module type", false);
  RootWTable& typesTable = chipContent.addTable();
  typesTable.setContent(chipPowerPerType.getContent());
  
  dumpRadiationTableSummary(myPage, powerSummaries, "Power in irradiated sensors", "W");
  dumpRadiationTableSummary(myPage, fluenceSummaries, "Fluence on sensors", "1-MeV-n-eq×cm"+superStart+"-2"+superEnd);

  // Some helper string objects
  ostringstream tempSS;
  std::string tempString;

  struct SensorsIrradiationPower {
    double operator()(const Module& m) { return m.sensorsIrradiationPowerMean(); }  // W
  };

  struct TotalPower {
    double operator()(const Module& m) { return m.sensorsIrradiationPowerMean() + m.totalPower() * Units::mW; }  // W (convert m.totalPower() from mW to W)
  };

  PlotDrawer<YZ, SensorsIrradiationPower, Average> yzSensorsPowerDrawer(0, 0);
  PlotDrawer<YZ, TotalPower, Average> yzTotalPowerDrawer(0, 0);

  yzSensorsPowerDrawer.addModules<CheckType<BARREL | ENDCAP>>(tracker.modules().begin(), tracker.modules().end());
  yzTotalPowerDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end(), BARREL | ENDCAP);

  RootWContent& myContent = myPage.addContent("Power maps", true);

  TCanvas sensorsIrradiationPowerCanvas;
  TCanvas totalPowerCanvas;

  yzSensorsPowerDrawer.drawFrame<HistogramFrameStyle>(sensorsIrradiationPowerCanvas);
  yzSensorsPowerDrawer.drawModules<ContourStyle>(sensorsIrradiationPowerCanvas);


  yzTotalPowerDrawer.drawFrame<HistogramFrameStyle>(totalPowerCanvas);
  yzTotalPowerDrawer.drawModules<ContourStyle>(totalPowerCanvas);

  RootWImage& sensorsIrradiationPowerImage = myContent.addImage(sensorsIrradiationPowerCanvas, insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsIrradiationPowerImage.setComment("Power dissipation in irradiated sensors (due to leakage current) (average per module) (W)");
  sensorsIrradiationPowerImage.setName("sensorsIrradiationPowerMap");
  RootWImage& totalPowerImage = myContent.addImage(totalPowerCanvas, insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  totalPowerImage.setComment("Total power dissipation in irradiated modules (W)");
  totalPowerImage.setName("totalPowerMap");
  
  // Add csv file with sensors irradiation handful info
  RootWContent* filesContent = new RootWContent("power csv files", false);
  myPage.addContent(filesContent);
  RootWTextFile* myTextFile;
  myTextFile = new RootWTextFile(Form("sensorsIrradiation%s.csv", trackerName.c_str()), "Sensors irradiation file");
  myTextFile->addText(createSensorsIrradiationCsv());
  filesContent->addItem(myTextFile);

}



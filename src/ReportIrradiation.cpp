#include <ReportIrradiation.hh>
#include <PlotDrawer.h>

void ReportIrradiation::analyze(const Tracker& tracker) {
  computeIrradiationPowerConsumption(tracker);
  preparePowerHistograms();
}

void ReportIrradiation::preparePowerHistograms() {
  myMapBag.clearMaps(mapBag::sensorsIrradiationPowerMap);
  myMapBag.clearMaps(mapBag::totalPowerConsumptionMap);
  TH2D& sensorsIrradiationPowerMap = myMapBag.getMaps(mapBag::sensorsIrradiationPowerMap)[mapBag::dummyMomentum]; // dummyMomentum is supplied because it is a single map. Multiple maps are indexed like arrays (see above efficiency maps)
  TH2D& totalPowerConsumptionMap = myMapBag.getMaps(mapBag::totalPowerConsumptionMap)[mapBag::dummyMomentum]; // dummyMomentum is supplied because it is a single map. Multiple maps are indexed like arrays (see above efficiency maps)
  prepareTrackerMap(sensorsIrradiationPowerMap, "sensorsIrradiationPowerMap", "Map of power dissipation in sensors (after irradiation)");
  prepareTrackerMap(totalPowerConsumptionMap, "irradiationPowerConsumptionMap", "Map of power dissipation in modules (after irradiation)");
}

void ReportIrradiation::computeIrradiationPowerConsumption(const Tracker& tracker) {
  IrradiationPowerVisitor irradiation_;
  irradiation_.preVisit();
  simParms_->accept(irradiation_);
  tracker.accept(irradiation_);
  irradiation_.postVisit();

  sensorsIrradiationPowerSummary_ = irradiation_.sensorsIrradiationPowerSummary;
  sensorsIrradiationSummary_ = irradiation_.sensorsIrradiationSummary;
}

void ReportIrradiation::visualizeTo(RootWSite& site) {
}

  std::string trackerName = tracker.myid();
  std::string pageName = "Irradiation (" + trackerName + ")";
  std::string pageAddress = "irradiation_" + trackerName + ".html";

  RootWPage& myPage = site.addPage(pageName);
  myPage.setAddress(pageAddress);

  std::map<std::string, SummaryTable>& powerSummaries = a.getSensorsIrradiationPowerSummary();
  std::map<std::string, SummaryTable>& irradiationSummaries = a.getSensorsIrradiationSummary();
  dumpRadiationTableSummary(myPage, powerSummaries, "Power in irradiated sensors", "W");
  dumpRadiationTableSummary(myPage, irradiationSummaries, "Fluence on sensors", "1-MeV-n-eq×cm"+superStart+"-2"+superEnd);

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

  RootWImage& sensorsIrradiationPowerImage = myContent.addImage(sensorsIrradiationPowerCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  sensorsIrradiationPowerImage.setComment("Power dissipation in irradiated sensors (due to leakage current) (average per module) (W)");
  sensorsIrradiationPowerImage.setName("sensorsIrradiationPowerMap");
  RootWImage& totalPowerImage = myContent.addImage(totalPowerCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
  totalPowerImage.setComment("Total power dissipation in irradiated modules (W)");
  totalPowerImage.setName("totalPowerMap");


  // Add csv file with sensors irradiation handful info
  RootWContent* filesContent = new RootWContent("power csv files", false);
  myPage.addContent(filesContent);
  RootWTextFile* myTextFile;
  myTextFile = new RootWTextFile(Form("sensorsIrradiation%s.csv", trackerName.c_str()), "Sensors irradiation file");
  myTextFile->addText(createSensorsIrradiationCsv(tracker));
  filesContent->addItem(myTextFile);

  return true;
}



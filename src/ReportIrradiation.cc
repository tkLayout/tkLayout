#include <ReportIrradiation.hh>
#include <PlotDrawer.hh>
#include <TLegend.h>

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
  doseSummaries = irradiation_.sensorsDoseSummary;
  fluenceSummaryPerType = irradiation_.sensorsFluencePerType;
  doseSummaryPerType = irradiation_.sensorsDosePerType;
  lumiInfo_ = irradiation_.lumiInformation;
  mapNames_ = irradiation_.mapInformation;
  maxFluence_ = irradiation_.maxFluence;
  maxDose_ = irradiation_.maxDose;
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
      output_ << "Module_DetId,Section,Layer,Ring,moduleType,dsDistance,isOuterRadiusRod_bool,r_mm,z_mm,operatingTemperature_Celsius,biasVoltage_V,meanWidth_mm,length_mm,sensorThickness_mm,sensorsVolume_totalPerModule_mm3,sensorsPowerMean_W,sensorsPowerMax_W,sensorsFluenceMean_Hb,sensorsFluenceMax_Hb,sensorsDoseMean_Gy,sensorsDoseMax_Gy" << std::endl;
    }
    void visit(const Barrel& b) { sectionName_ = b.myid(); }
    void visit(const Endcap& e) { sectionName_ = e.myid(); }
    void visit(const Layer& l)  { layerId_ = l.myid(); }
    void visit(const RodPair& r)  { isOuterRadiusRod_ = r.isOuterRadiusRod(); }
    void visit(const Disk& d)  { isOuterRadiusRod_ = false; layerId_ = d.myid(); } // no rod here !
    void visit(const Module& m) {
      output_  << m.myDetId() << ","
	       << sectionName_ << ","
	       << layerId_ << ","
	       << m.moduleRing() << ","
	       << m.moduleType() << ","
	       << m.dsDistance() << ","
	       << isOuterRadiusRod_ << ","
	       << std::fixed << std::setprecision(6)
	       << m.center().Rho() << ","
	       << m.center().Z() << ","
	       << m.operatingTemp() << ","
	       << m.biasVoltage() << ","
	       << m.meanWidth() << ","
	       << m.length() << ","
	       << m.sensorThickness() << ","
	       << m.totalSensorsVolume() << ","
	       << std::fixed << std::setprecision(3)
	       << m.sensorsIrradiationPowerMean() << ","
	       << m.sensorsIrradiationPowerMax() << ","
	       << m.sensorsIrradiationMean() << ","
	       << m.sensorsIrradiationMax() << ","
	       << m.sensorsDoseMean() << ","
	       << m.sensorsDoseMax()
	       << std::endl;
    }

    std::string output() const { return output_.str(); }
  };

  TrackerVisitor v;
  v.preVisit();
  tracker.accept(v);
  return v.output();
}

std::map<std::string,TH1F*> ReportIrradiation::createSensorsIrradiationHistograms(double max_fluence, double max_dose){
  class TrackerVisitor : public ConstGeometryVisitor {
    std::map<std::string,TH1F*> hist_map;
    string sectionName_;
    int layerId_;
    bool isOuterRadiusRod_;
    double theMaxFluence_;
    double theMaxDose_;
  public:
    void setHistogramMaximums(double fluence, double dose){
      theMaxFluence_ = fluence;
      theMaxDose_ = dose;
    }
    void preVisit() {
      for (auto x : {"PS" , "2S", "PS-1.6mm","PS-2.6mm","PS-4.0mm","2S-1.8mm","2S-4.0mm","TEPX" , "TFPX", "TBPX"}) hist_map[std::string("tot_mod_fluence_")+x] = new TH1F(TString::Format("tot_mod_fluence_%s", x),";NIEL fluence[1 MeV n_{eq} / cm^2]; Number of modules",30,0.,1.05*theMaxFluence_);
      for (auto x : {"PS" , "2S", "PS-1.6mm","PS-2.6mm","PS-4.0mm","2S-1.8mm","2S-4.0mm","TEPX" , "TFPX", "TBPX"}) hist_map[std::string("tot_mod_TID_")+x] = new TH1F(TString::Format("tot_mod_TID_%s", x),";Absorbed dose[Gy]; Number of modules",30,0,1.05*theMaxDose_);
    }
    void visit(const Barrel& b) { sectionName_ = b.myid(); }
    void visit(const Endcap& e) { sectionName_ = e.myid(); }
    void visit(const Layer& l)  { layerId_ = l.myid(); }
    void visit(const RodPair& r)  { isOuterRadiusRod_ = r.isOuterRadiusRod(); }
    void visit(const Disk& d)  { isOuterRadiusRod_ = false; layerId_ = d.myid(); } // no rod here !
    void visit(const Module& m) {
      if(sectionName_=="TBPX"){
        hist_map[std::string("tot_mod_fluence_TBPX")]->Fill(m.sensorsIrradiationMean());
        hist_map[std::string("tot_mod_TID_TBPX")]->Fill(m.sensorsDoseMean());
      } else if (sectionName_=="TFPX"){
        hist_map[std::string("tot_mod_fluence_TFPX")]->Fill(m.sensorsIrradiationMean());
        hist_map[std::string("tot_mod_TID_TFPX")]->Fill(m.sensorsDoseMean());
      } else if (sectionName_=="TEPX"){
        hist_map[std::string("tot_mod_fluence_TEPX")]->Fill(m.sensorsIrradiationMean());
        hist_map[std::string("tot_mod_TID_TEPX")]->Fill(m.sensorsDoseMean());
      } else if (m.moduleType()=="pt2S"){
        hist_map[std::string("tot_mod_fluence_2S")]->Fill(m.sensorsIrradiationMean());
        hist_map[std::string("tot_mod_TID_2S")]->Fill(m.sensorsDoseMean());
        if(m.dsDistance()<2.){
          hist_map[std::string("tot_mod_fluence_2S-1.8mm")]->Fill(m.sensorsIrradiationMean());
          hist_map[std::string("tot_mod_TID_2S-1.8mm")]->Fill(m.sensorsDoseMean());
        } else {
          hist_map[std::string("tot_mod_fluence_2S-4.0mm")]->Fill(m.sensorsIrradiationMean());
          hist_map[std::string("tot_mod_TID_2S-4.0mm")]->Fill(m.sensorsDoseMean());
        }
      } else if (m.moduleType()=="ptPS"){
        hist_map[std::string("tot_mod_fluence_PS")]->Fill(m.sensorsIrradiationMean());
        hist_map[std::string("tot_mod_TID_PS")]->Fill(m.sensorsDoseMean());
        if(m.dsDistance()<2.){
          hist_map[std::string("tot_mod_fluence_PS-1.6mm")]->Fill(m.sensorsIrradiationMean());
          hist_map[std::string("tot_mod_TID_PS-1.6mm")]->Fill(m.sensorsDoseMean());
        } else if (m.dsDistance()<3.){
          hist_map[std::string("tot_mod_fluence_PS-2.6mm")]->Fill(m.sensorsIrradiationMean());
          hist_map[std::string("tot_mod_TID_PS-2.6mm")]->Fill(m.sensorsDoseMean());
        } else {
          hist_map[std::string("tot_mod_fluence_PS-4.0mm")]->Fill(m.sensorsIrradiationMean());
          hist_map[std::string("tot_mod_TID_PS-4.0mm")]->Fill(m.sensorsDoseMean());
        }
      }
    }

    std::map<std::string,TH1F*> output()  { 
     for (auto x : {"PS","PS-1.6mm","PS-2.6mm","PS-4.0mm", "2S" , "2S-1.8mm","2S-4.0mm", "TEPX","TFPX","TBPX"}){
       hist_map[std::string("cumulative_fluence_")+x] = dynamic_cast<TH1F*>(hist_map[std::string("tot_mod_fluence_")+x]->GetCumulative(kFALSE));
       hist_map[std::string("cumulative_fluence_")+x]->SetName(TString::Format("cumulative_fluence_%s",x));
       hist_map[std::string("cumulative_fluence_")+x]->GetYaxis()->SetTitle("Fraction of modules");
       hist_map[std::string("cumulative_fluence_")+x]->Scale(1./hist_map[std::string("cumulative_fluence_")+x]->GetBinContent(1));
       hist_map[std::string("cumulative_fluence_")+x]->SetTitle("Fraction of modules with fluence above");
       hist_map[std::string("cumulative_TID_")+x] = dynamic_cast<TH1F*>(hist_map[std::string("tot_mod_TID_")+x]->GetCumulative(kFALSE));
       hist_map[std::string("cumulative_TID_")+x]->SetName(TString::Format("cumulative_TID_%s",x));
       hist_map[std::string("cumulative_TID_")+x]->SetTitle("Fraction of modules with dose above");
       hist_map[std::string("cumulative_TID_")+x]->GetYaxis()->SetTitle("Fraction of modules");
       hist_map[std::string("cumulative_TID_")+x]->Scale(1./hist_map[std::string("cumulative_TID_")+x]->GetBinContent(1));
     }
      std::map<std::string,int>  hist_colours;
      hist_colours["PS"] = 632;
      hist_colours["PS-1.6mm"] = 632;
      hist_colours["PS-2.6mm"] = 910;
      hist_colours["PS-4.0mm"] = 803;
      hist_colours["2S"] = 419;
      hist_colours["2S-1.8mm"] = 419;
      hist_colours["2S-4.0mm"] = 860;
      hist_colours["TEPX"] = 600;
      hist_colours["TFPX"] = 810;
      hist_colours["TBPX"] = 922;
      std::map<std::string,float> bar_offsets;
      bar_offsets["PS"] = 0.2;
      bar_offsets["PS-1.6mm"] = 0.0;
      bar_offsets["PS-2.6mm"] = 0.2;
      bar_offsets["PS-4.0mm"] = 0.4;
      bar_offsets["2S-1.8mm"] = 0.6;
      bar_offsets["2S-4.0mm"] = 0.8;
      bar_offsets["2S"] = 0.4;
      bar_offsets["TEPX"] = 0.2;
      bar_offsets["TFPX"] = 0.4;
      bar_offsets["TBPX"] = 0.6;
      for (auto x : {"PS" ,"PS-1.6mm","PS-2.6mm","PS-4.0mm", "2S","2S-1.8mm","2S-4.0mm", "TEPX" ,"TFPX", "TBPX"}){
        for (auto y : {"tot_mod_fluence_", "tot_mod_TID_", "cumulative_fluence_", "cumulative_TID_"}){
          hist_map[std::string(y)+std::string(x)]->SetLineColor(hist_colours[x]);
          hist_map[std::string(y)+std::string(x)]->SetFillColor(hist_colours[x]);
          hist_map[std::string(y)+std::string(x)]->SetBarWidth(0.2);
          hist_map[std::string(y)+std::string(x)]->SetBarOffset(bar_offsets[x]);
        }
      }
     return hist_map; }
  };

  TrackerVisitor v;
  v.setHistogramMaximums(max_fluence,max_dose);
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

  RootWContent& settingsContent = myPage.addContent("Simulation settings");
  RootWInfo* lumInfo;
  lumInfo = new RootWInfo("Integrated luminosity");
  lumInfo->setValue(lumiInfo_);
  settingsContent.addItem(lumInfo);
  RootWInfo* mapNamesInfo;
  mapNamesInfo = new RootWInfo("Dose and irradiation map names");
  mapNamesInfo->setValue(mapNames_);
  settingsContent.addItem(mapNamesInfo);
  RootWInfo* notesInfo;
  notesInfo = new RootWInfo("Note");
  std::string notesString_ = "The NIEL(irrad) and TID(dose) total radiation is computed for each module at the 4 corners of the module,"
    " and at the central point of each side.<br/>\n"
    "The average value is computed for each module, and the highest of thse is reported under \"Max. module irrad(dose)\", while the 95-th percentile is reported under \"95% module irrad(dose)\"<br/>\n"
    "The dose (NIEL or TID) values are obtained from a linear interpolation of the nearest points in the (r,z) plane in the radiation map quoted here above."
    "If more than one map covers the same point, the one with the finer grid is used (typically, the Inner Tracker is covered by a finer grid).\n";
  notesInfo->setValue(notesString_);
  settingsContent.addItem(notesInfo);
  notesInfo = new RootWInfo("Note");
  notesInfo->setValue("the total integrated luminosity used here (" + lumiInfo_ + ") is used for scaling, but it is not necessarily the most up-to-date estimate.\n");
  settingsContent.addItem(notesInfo);
  notesInfo = new RootWInfo("Note");
  notesInfo->setValue("the estimate of bias current is based on a linear scaling: all relevant parameters are under the <a href=\"info.html\">info tab</a>\n");
  settingsContent.addItem(notesInfo);  



  // Irradiation on each module type (fine grained)
  RootWContent& summaryContent = myPage.addContent("Irradiation summary per module type");
  RootWTable& summaryTable = summaryContent.addTable();
  summaryTable.setContent(fluenceSummaryPerType.getContent());

  // Dose on each module type (fine grained)
  RootWContent& doseSummaryContent = myPage.addContent("Dose summary per module type");
  RootWTable& doseSummaryTable = doseSummaryContent.addTable();
  doseSummaryTable.setContent(doseSummaryPerType.getContent());


  // Power consumption per module tpye (coarse grained)
  RootWContent& chipContent = myPage.addContent("Chip power consumption per module type", false);
  RootWTable& typesTable = chipContent.addTable();
  typesTable.setContent(chipPowerPerType.getContent());
  
  dumpRadiationTableSummary(myPage, powerSummaries, "Power in irradiated sensors", "W");
  dumpRadiationTableSummary(myPage, fluenceSummaries, "Niel fluence on sensors", "1-MeV-n-eq×cm"+superStart+"-2"+superEnd);
  dumpRadiationTableSummary(myPage, doseSummaries, "Dose on sensors", "Gy");

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

  std::unique_ptr<TCanvas> sensorsIrradiationPowerCanvas(new TCanvas());
  std::unique_ptr<TCanvas> totalPowerCanvas(new TCanvas());

  yzSensorsPowerDrawer.drawFrame<HistogramFrameStyle>(*sensorsIrradiationPowerCanvas.get());
  yzSensorsPowerDrawer.drawModules<ContourStyle>(*sensorsIrradiationPowerCanvas.get());


  yzTotalPowerDrawer.drawFrame<HistogramFrameStyle>(*totalPowerCanvas.get());
  yzTotalPowerDrawer.drawModules<ContourStyle>(*totalPowerCanvas.get());

  RootWImage& sensorsIrradiationPowerImage = myContent.addImage(std::move(sensorsIrradiationPowerCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsIrradiationPowerImage.setComment("Power dissipation in irradiated sensors (due to leakage current) (average per module) (W)");
  sensorsIrradiationPowerImage.setName("sensorsIrradiationPowerMap");
  RootWImage& totalPowerImage = myContent.addImage(std::move(totalPowerCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  totalPowerImage.setComment("Total power dissipation in irradiated modules (W)");
  totalPowerImage.setName("totalPowerMap");

   
  RootWContent& myupdContent = myPage.addContent("Number of modules per dose and fluence figures", true);

  std::map<std::string,TH1F*> histos = createSensorsIrradiationHistograms(maxFluence_, maxDose_);
  std::unique_ptr<TCanvas> irradiationFullCanvas(new TCanvas());
  irradiationFullCanvas->cd();
  TLegend* legFullCanvas = new TLegend(0.7,0.7,0.9,0.9);
  legFullCanvas->SetFillStyle(0);
  int drawnhists=0;
  for (auto x : {"2S" , "PS", "TEPX" , "TFPX", "TBPX"}){
    if(histos[std::string("tot_mod_fluence_")+x]->GetEntries() > 0 ){
      if(drawnhists==0){
        histos[std::string("tot_mod_fluence_")+x]->GetYaxis()->SetRangeUser(0,1.5*histos[std::string("tot_mod_fluence_")+x]->GetMaximum());
        histos[std::string("tot_mod_fluence_")+x]->DrawCopy("B");
        drawnhists+=1;
      } else {
        histos[std::string("tot_mod_fluence_")+x]->DrawCopy("BSAME");
      }
      legFullCanvas->AddEntry(histos[std::string("tot_mod_fluence_")+x],x,"F");
    }
  }
  if(drawnhists>0){
    legFullCanvas->Draw("SAME");
  }


  std::unique_ptr<TCanvas> irradiationFullCumulCanvas(new TCanvas());
  irradiationFullCumulCanvas->cd();
  TLegend* legFullCumulCanvas = new TLegend(0.7,0.7,0.9,0.9);
  legFullCumulCanvas->SetFillStyle(0);
  drawnhists=0;
  for (auto x : {"2S" , "PS", "TEPX" , "TFPX", "TBPX"}){
    if(histos[std::string("tot_mod_fluence_")+x]->GetEntries() > 0 ){
      if(drawnhists==0){
        histos[std::string("cumulative_fluence_")+x]->DrawCopy("L");
        drawnhists+=1;
      } else {
        histos[std::string("cumulative_fluence_")+x]->DrawCopy("LSAME");
      }
      legFullCumulCanvas->AddEntry(histos[std::string("cumulative_fluence_")+x],x,"L");
    }
  }
  if(drawnhists>0){
    legFullCumulCanvas->Draw("SAME");
  }


  std::unique_ptr<TCanvas> irradiationSplitCanvas(new TCanvas());
  irradiationSplitCanvas->cd();
  TLegend* legSplitCanvas = new TLegend(0.7,0.7,0.9,0.9);
  legSplitCanvas->SetFillStyle(0);
  drawnhists=0;
  for (auto x : {"2S-1.8mm","2S-4.0mm","PS-1.6mm", "PS-2.6mm" ,"PS-4.0mm"}){
    if(histos[std::string("tot_mod_fluence_")+x]->GetEntries() > 0 ){
      if(drawnhists==0){
        histos[std::string("tot_mod_fluence_")+x]->GetYaxis()->SetRangeUser(0,1.5*histos[std::string("tot_mod_fluence_")+x]->GetMaximum());
        histos[std::string("tot_mod_fluence_")+x]->DrawCopy("B");
        drawnhists+=1;
      } else {
        histos[std::string("tot_mod_fluence_")+x]->DrawCopy("BSAME");
      }
      legSplitCanvas->AddEntry(histos[std::string("tot_mod_fluence_")+x],x,"F");
    }
  }
  if(drawnhists>0){
    legSplitCanvas->Draw("SAME");
  }

  std::unique_ptr<TCanvas> irradiationSplitCumulCanvas(new TCanvas());
  irradiationSplitCumulCanvas->cd();
  TLegend* legSplitCumulCanvas = new TLegend(0.7,0.7,0.9,0.9);
  legSplitCumulCanvas->SetFillStyle(0);
  drawnhists=0;
  for (auto x : {"2S-1.8mm","2S-4.0mm","PS-1.6mm","PS-2.6mm", "PS-4.0mm"}){
    if(histos[std::string("tot_mod_fluence_")+x]->GetEntries() > 0 ){
      if(drawnhists==0){
        histos[std::string("cumulative_fluence_")+x]->DrawCopy("L");
        drawnhists+=1;
      } else {
        histos[std::string("cumulative_fluence_")+x]->DrawCopy("LSAME");
      }
      legSplitCumulCanvas->AddEntry(histos[std::string("cumulative_fluence_")+x],x,"L");
    }
  }
  if(drawnhists>0){
    legSplitCumulCanvas->Draw("SAME");
  }



  std::unique_ptr<TCanvas> doseFullCanvas(new TCanvas());
  doseFullCanvas->cd();
  TLegend* legDoseFullCanvas = new TLegend(0.7,0.7,0.9,0.9);
  legDoseFullCanvas->SetFillStyle(0);
  drawnhists=0;
  for (auto x : {"2S" , "PS", "TEPX" , "TFPX", "TBPX"}){
    if(histos[std::string("tot_mod_TID_")+x]->GetEntries() > 0 ){
      if(drawnhists==0){
        histos[std::string("tot_mod_TID_")+x]->GetYaxis()->SetRangeUser(0,1.5*histos[std::string("tot_mod_TID_")+x]->GetMaximum());
        histos[std::string("tot_mod_TID_")+x]->DrawCopy("B");
        drawnhists+=1;
      } else {
       histos[std::string("tot_mod_TID_")+x]->DrawCopy("BSAME");
      }
      legDoseFullCanvas->AddEntry(histos[std::string("tot_mod_TID_")+x],x,"F");
    }
  }
  if(drawnhists>0){
    legDoseFullCanvas->Draw("SAME");
  }

  std::unique_ptr<TCanvas> doseFullCumulCanvas(new TCanvas());
  doseFullCumulCanvas->cd();
  TLegend* legDoseFullCumulCanvas = new TLegend(0.7,0.7,0.9,0.9);
  legDoseFullCumulCanvas->SetFillStyle(0);
  drawnhists=0;
  for (auto x : {"2S" , "PS", "TEPX" , "TFPX", "TBPX"}){
    if(histos[std::string("tot_mod_TID_")+x]->GetEntries() > 0 ){
      if(drawnhists==0){
        histos[std::string("cumulative_TID_")+x]->DrawCopy("L");
        drawnhists+=1;
      } else {
       histos[std::string("cumulative_TID_")+x]->DrawCopy("LSAME");
      }
      legDoseFullCumulCanvas->AddEntry(histos[std::string("cumulative_TID_")+x],x,"L");
    }
  }
  if(drawnhists>0){
    legDoseFullCumulCanvas->Draw("SAME");
  }

  std::unique_ptr<TCanvas> doseSplitCanvas(new TCanvas());
  doseSplitCanvas->cd();
  TLegend* legDoseSplitCanvas = new TLegend(0.7,0.7,0.9,0.9);
  legDoseSplitCanvas->SetFillStyle(0);
  drawnhists=0;
  for (auto x : {"2S-1.8mm","2S-4.0mm","PS-1.6mm","PS-2.6mm","PS-4.0mm"}){
    if(histos[std::string("tot_mod_TID_")+x]->GetEntries() > 0 ){
      if(drawnhists==0){
        histos[std::string("tot_mod_TID_")+x]->GetYaxis()->SetRangeUser(0,1.5*histos[std::string("tot_mod_TID_")+x]->GetMaximum());
        histos[std::string("tot_mod_TID_")+x]->DrawCopy("B");
        drawnhists+=1;
      } else {
       histos[std::string("tot_mod_TID_")+x]->DrawCopy("BSAME");
      }
      legDoseSplitCanvas->AddEntry(histos[std::string("tot_mod_TID_")+x],x,"F");
    }
  }
  legDoseSplitCanvas->Draw("SAME");


  std::unique_ptr<TCanvas> doseSplitCumulCanvas(new TCanvas());
  doseSplitCumulCanvas->cd();
  TLegend* legDoseSplitCumulCanvas = new TLegend(0.7,0.7,0.9,0.9);
  legDoseSplitCumulCanvas->SetFillStyle(0);
  drawnhists=0;
  for (auto x : {"2S-1.8mm","2S-4.0mm","PS-1.6mm","PS-2.6mm" , "PS-4.0mm"}){
    if(histos[std::string("tot_mod_TID_")+x]->GetEntries() > 0 ){
      if(drawnhists==0){
        histos[std::string("cumulative_TID_")+x]->DrawCopy("C");
        drawnhists+=1;
      } else {
       histos[std::string("cumulative_TID_")+x]->DrawCopy("CSAME");
      }
      legDoseSplitCumulCanvas->AddEntry(histos[std::string("cumulative_TID_")+x],x,"L");
    }
  }
  legDoseSplitCumulCanvas->Draw("SAME");

  RootWImage& sensorsIrradFullHistogram = myupdContent.addImage(std::move(irradiationFullCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsIrradFullHistogram.setName("sensorsIrradiationFull");

  RootWImage& sensorsIrradFullCumulHistogram = myupdContent.addImage(std::move(irradiationFullCumulCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsIrradFullCumulHistogram.setName("sensorsIrradiationFullCumulative");

  RootWImage& sensorsIrradSplitHistogram = myupdContent.addImage(std::move(irradiationSplitCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsIrradSplitHistogram.setName("sensorsIrradiationSplit");

  RootWImage& sensorsIrradSplitCumulHistogram = myupdContent.addImage(std::move(irradiationSplitCumulCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsIrradSplitCumulHistogram.setName("sensorsIrradiationSplitCumulative");


   
  RootWImage& sensorsDoseFullHistogram = myupdContent.addImage(std::move(doseFullCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsDoseFullHistogram.setName("sensorsDoseFull");

  RootWImage& sensorsDoseFullCumulHistogram = myupdContent.addImage(std::move(doseFullCumulCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsDoseFullCumulHistogram.setName("sensorsDoseFullCumul");

   RootWImage& sensorsDoseSplitHistogram = myupdContent.addImage(std::move(doseSplitCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsDoseSplitHistogram.setName("sensorsDoseSplit");

  RootWImage& sensorsDoseSplitCumulHistogram = myupdContent.addImage(std::move(doseSplitCumulCanvas), insur::vis_std_canvas_sizeX, insur::vis_min_canvas_sizeY);
  sensorsDoseSplitCumulHistogram.setName("sensorsDoseSplitCumul");

 
  // Add csv file with sensors irradiation handful info
  RootWContent* filesContent = new RootWContent("power csv files", false);
  myPage.addContent(filesContent);
  RootWTextFile* myTextFile;
  myTextFile = new RootWTextFile(Form("sensorsIrradiation%s.csv", trackerName.c_str()), "Sensors irradiation file");
  myTextFile->addText(createSensorsIrradiationCsv());
  filesContent->addItem(myTextFile);

}



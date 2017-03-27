#ifndef ANALYZERVISITOR_H
#define ANALYZERVISITOR_H

#include <string>
#include <map>
#include <vector>
#include <utility>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "Tracker.hh"
#include "Barrel.hh"
#include "Endcap.hh"
#include "Layer.hh"
#include "Disk.hh"
#include "RodPair.hh"
#include "Ring.hh"
#include "Module.hh"
#include "SimParms.hh"
#include "PtErrorAdapter.hh"
#include "MessageLogger.hh"

#include "Visitor.hh"

#include "Bag.hh"
#include "SummaryTable.hh"

#include "AnalyzerVisitors/TriggerProcessorBandwidth.hh"
#include "AnalyzerVisitors/IrradiationPower.hh"
#include "AnalyzerVisitors/Bandwidth.hh"
#include "AnalyzerVisitors/TriggerDistanceTuningPlots.hh"
#include "AnalyzerVisitors/TriggerFrequency.hh"

using std::string;
using std::map;
using std::vector;
using std::pair;




namespace AnalyzerHelpers {

  void drawModuleOnMap(const Module& m, double val, TH2D& map, TH2D& counter);
  void drawModuleOnMap(const Module& m, double val, TH2D& map);

}



class TriggerEfficiencyMapVisitor : public ConstGeometryVisitor {
  double myPt_;
  TH2D& myMap_;
  TH2D* counter_;
public:
  TriggerEfficiencyMapVisitor(TH2D& map, double pt) : myMap_(map), myPt_(pt) { counter_ = (TH2D*)map.Clone(); }

  void visit(const DetectorModule& aModule) {
    // returns immediately if module is not pt enabled
    if (aModule.sensorLayout() != PT) return;
    double myValue = PtErrorAdapter(aModule).getTriggerProbability(myPt_);
    if (myValue>=0) AnalyzerHelpers::drawModuleOnMap(aModule, myValue, myMap_, *counter_);
  }

  void postVisit() {
    for (int i=1; i<=myMap_.GetNbinsX(); ++i)
      for (int j=1; j<=myMap_.GetNbinsY(); ++j)
        if (counter_->GetBinContent(i,j)!=0)
          myMap_.SetBinContent(i,j, myMap_.GetBinContent(i,j) / counter_->GetBinContent(i,j));
    // ... and get rid of the counter
  }

  ~TriggerEfficiencyMapVisitor() { delete counter_; }
};





class PtThresholdMapVisitor : public ConstGeometryVisitor {
  double myPt_;
  TH2D& myMap_;
  TH2D* counter_;
public:
  PtThresholdMapVisitor(TH2D& map, double pt) : myMap_(map), myPt_(pt) { counter_ = (TH2D*)map.Clone(); }

  void visit(const DetectorModule& aModule) {
    if (aModule.sensorLayout() != PT) return;
    double myValue = PtErrorAdapter(aModule).getPtThreshold(myPt_);
    if (myValue >= 0) AnalyzerHelpers::drawModuleOnMap(aModule, myValue, myMap_, *counter_);
  }

  void postVisit() {
    for (int i=1; i<=myMap_.GetNbinsX(); ++i)
      for (int j=1; j<=myMap_.GetNbinsY(); ++j)
        if (counter_->GetBinContent(i,j)!=0)
          myMap_.SetBinContent(i,j, myMap_.GetBinContent(i,j) / counter_->GetBinContent(i,j));
    // ... and get rid of the counter
  }

  ~PtThresholdMapVisitor() { delete counter_; }
};





class SpacingCutVisitor : public ConstGeometryVisitor {
  TH2D& suggestedSpacingMap_;
  TH2D& suggestedSpacingMapAW_;
  TH2D& nominalCutMap_;
  TH2D *counter_, *counterSpacing_, *counterSpacingAW_;
  ModuleOptimalSpacings& moduleOptimalSpacings_;
public:
  SpacingCutVisitor(TH2D& suggestedSpacingMap, TH2D& suggestedSpacingMapAW, TH2D& nominalCutMap, ModuleOptimalSpacings& moduleOptimalSpacings) : 
      suggestedSpacingMap_(suggestedSpacingMap), suggestedSpacingMapAW_(suggestedSpacingMapAW), nominalCutMap_(nominalCutMap), moduleOptimalSpacings_(moduleOptimalSpacings) {
        counter_ = (TH2D*)suggestedSpacingMap.Clone();
        counterSpacing_ = (TH2D*)suggestedSpacingMap.Clone();
        counterSpacingAW_ = (TH2D*)suggestedSpacingMapAW.Clone();
  }

  void visit(const DetectorModule& aModule) {
    if (aModule.sensorLayout() != PT) return;
    double mySuggestedSpacing = moduleOptimalSpacings_[&aModule][5]; // TODO: put this 5 in a configuration of some sort
    double mySuggestedSpacingAW = moduleOptimalSpacings_[&aModule][aModule.triggerWindow()];
    double nominalCut = PtErrorAdapter(aModule).getPtCut();

    AnalyzerHelpers::drawModuleOnMap(aModule, nominalCut, nominalCutMap_);
    if (mySuggestedSpacing != 0) AnalyzerHelpers::drawModuleOnMap(aModule, mySuggestedSpacing, suggestedSpacingMap_, *counterSpacing_);
    if (mySuggestedSpacingAW != 0) AnalyzerHelpers::drawModuleOnMap(aModule, mySuggestedSpacingAW, suggestedSpacingMapAW_, *counterSpacingAW_);

  }

  void postVisit() {
    for (int i=1; i<=suggestedSpacingMap_.GetNbinsX(); ++i) {
      for (int j=1; j<=suggestedSpacingMap_.GetNbinsY(); ++j) {
        if (counter_->GetBinContent(i,j)!=0) {
          nominalCutMap_.SetBinContent(i,j, nominalCutMap_.GetBinContent(i,j) / counter_->GetBinContent(i,j));
          if ((suggestedSpacingMap_.GetBinContent(i,j)/counterSpacing_->GetBinContent(i,j))>50) {
            std::cout << "debug: for bin " << i << ", " << j << " suggestedSpacing is " << suggestedSpacingMap_.GetBinContent(i,j)
              << " and counter is " << counterSpacing_->GetBinContent(i,j) << std::endl;
          }
          suggestedSpacingMap_.SetBinContent(i,j, suggestedSpacingMap_.GetBinContent(i,j) / counterSpacing_->GetBinContent(i,j));
          suggestedSpacingMapAW_.SetBinContent(i,j, suggestedSpacingMapAW_.GetBinContent(i,j) / counterSpacingAW_->GetBinContent(i,j));
        }
      }
    }
  }

  ~SpacingCutVisitor() {
    delete counter_;
    delete counterSpacing_;
    delete counterSpacingAW_;
  }

};




#endif

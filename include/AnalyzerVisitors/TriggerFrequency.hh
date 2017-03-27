#ifndef TRIGGERFREQUENCY_H
#define TRIGGERFREQUENCY_H

#include <string>
#include <map>
#include <vector>
#include <utility>

#include <TH1.h>

#include "PtErrorAdapter.hh"
#include "Tracker.hh"
#include "SimParms.hh"

#include "Visitor.hh"
#include "SummaryTable.hh"

class TriggerFrequencyVisitor : public ConstGeometryVisitor {
  typedef std::map<std::pair<std::string, int>, TH1D*> StubRateHistos;

  std::map<std::string, std::map<std::pair<int,int>, int>>   triggerFrequencyCounts_;
  std::map<std::string, std::map<std::pair<int,int>, double>>  triggerFrequencyAverageTrue_, 
                                                               triggerFrequencyInterestingParticleTrue_, 
                                                               triggerFrequencyAverageFake_, 
                                                               triggerFrequencyAverageMisfiltered_,
                                                               triggerFrequencyAverageCombinatorial_,
                                                               triggerDataBandwidths_; // trigger frequency by module in Z and R, averaged over Phi
  StubRateHistos totalStubRateHistos_, trueStubRateHistos_;

  int nbins_;
  double bunchSpacingNs_, nMB_, interestingPt_;

  void setupSummaries(const string& cntName) {
    triggerFrequencyTrueSummaries[cntName].setHeader("Layer", "Ring");
    triggerFrequencyFakeSummaries[cntName].setHeader("Layer", "Ring");
    triggerFrequencyInterestingSummaries[cntName].setHeader("Layer", "Ring");
    triggerFrequencyMisfilteredSummaries[cntName].setHeader("Layer", "Ring");
    triggerFrequencyCombinatorialSummaries[cntName].setHeader("Layer", "Ring");
    triggerRateSummaries[cntName].setHeader("Layer", "Ring");
    triggerEfficiencySummaries[cntName].setHeader("Layer", "Ring");
    triggerPuritySummaries[cntName].setHeader("Layer", "Ring");
    triggerDataBandwidthSummaries[cntName].setHeader("Layer", "Ring");
    triggerFrequencyTrueSummaries[cntName].setPrecision(3);
    triggerFrequencyFakeSummaries[cntName].setPrecision(3);
    triggerFrequencyInterestingSummaries[cntName].setPrecision(3);
    triggerFrequencyInterestingSummaries[cntName].setPrecision(3);
    triggerFrequencyMisfilteredSummaries[cntName].setPrecision(3);
    triggerRateSummaries[cntName].setPrecision(3);
    triggerEfficiencySummaries[cntName].setPrecision(3);
    triggerPuritySummaries[cntName].setPrecision(3);
    triggerDataBandwidthSummaries[cntName].setPrecision(3);
  }
public:
  std::map<std::string, std::map<std::pair<int,int>, double>> triggerFrequenciesPerEvent;
  MultiSummaryTable triggerFrequencyTrueSummaries, 
                    triggerFrequencyFakeSummaries, 
                    triggerFrequencyInterestingSummaries, 
                    triggerFrequencyMisfilteredSummaries,
                    triggerFrequencyCombinatorialSummaries,
                    triggerRateSummaries, 
                    triggerEfficiencySummaries,
                    triggerPuritySummaries, 
                    triggerDataBandwidthSummaries,
                    stripOccupancySummaries,
                    hitOccupancySummaries;

  void visit(const SimParms& sp) {
    bunchSpacingNs_ = sp.bunchSpacingNs();
    nMB_ = sp.numMinBiasEvents();
    interestingPt_ = sp.triggerPtCut();
  }

  void visit(const Barrel& b) { setupSummaries(b.myid()); }

  void visit(const Layer& l) {
    nbins_ = l.numModulesPerRod();
  }

  void visit(const Endcap& e) { setupSummaries(e.myid()); }

  void visit(const Disk& d) {
    nbins_ = d.numRings();
  }

  void visit(const DetectorModule& module) {

    XYZVector center = module.center();
    // TODO: check this too
    if ((center.Z()<0) || module.posRef().phi > 2/*(center.Phi()<0) || (center.Phi()>M_PI/2)*/ || (module.dsDistance()==0.0)) return;

    TH1D* currentTotalHisto;
    TH1D* currentTrueHisto;

    string table = module.tableRef().table;
    int row = module.tableRef().row;
    int col = module.tableRef().col;

    PtErrorAdapter pterr(module);

    if (totalStubRateHistos_.count(std::make_pair(table, row)) == 0) {
      currentTotalHisto = new TH1D(("totalStubsPerEventHisto" + table + any2str(row)).c_str(), ";Modules;MHz/cm^2", nbins_, 0.5, nbins_+0.5);
      currentTrueHisto = new TH1D(("trueStubsPerEventHisto" + table + any2str(row)).c_str(), ";Modules;MHz/cm^2", nbins_, 0.5, nbins_+0.5); 
      totalStubRateHistos_[std::make_pair(table, row)] = currentTotalHisto; 
      trueStubRateHistos_[std::make_pair(table, row)] = currentTrueHisto; 
    } else {
      currentTotalHisto = totalStubRateHistos_[std::make_pair(table, row)]; 
      currentTrueHisto = trueStubRateHistos_[std::make_pair(table, row)]; 
    }


    int curCnt = triggerFrequencyCounts_[table][std::make_pair(row, col)]++;
    double curAvgTrue = triggerFrequencyAverageTrue_[table][std::make_pair(row, col)];
    double curAvgInteresting = triggerFrequencyInterestingParticleTrue_[table][std::make_pair(row, col)];
    //double curAvgFake = triggerFrequencyAverageFake_[table][std::make_pair(row, col)];
    double curAvgMisfiltered = triggerFrequencyAverageMisfiltered_[table][std::make_pair(row, col)];
    double curAvgCombinatorial = triggerFrequencyAverageCombinatorial_[table][std::make_pair(row, col)];

    //curAvgTrue  = curAvgTrue + (module->getTriggerFrequencyTruePerEvent()*tracker.getNMB() - curAvgTrue)/(curCnt+1);
    //curAvgFake  = curAvgFake + (module->getTriggerFrequencyFakePerEvent()*pow(tracker.getNMB(),2) - curAvgFake)/(curCnt+1); // triggerFrequencyFake scales with the square of Nmb!

    double highPtParticlesRate = pterr.getParticleFrequencyPerEventAbove(interestingPt_)*nMB_;
    double trueStubRate = pterr.getTriggerFrequencyTruePerEventAbove(interestingPt_)*nMB_; // highPtParticlesRate * triggerEfficiency
    double misfilteredStubRate = pterr.getTriggerFrequencyTruePerEventBelow(interestingPt_)*nMB_; // low-Pt particles improperly considered to be high-pT, due to pT measurement errors, for which we form stubs
    double combinatorialStubRate = pterr.getTriggerFrequencyFakePerEvent()*pow(nMB_,2); // stubs due to occupancy combinatorics - i.e. random pixels/strips turned on in the upper and lower sensors caused by separate tracks or secondaries which happen to fall within the trigger window
    double fakeStubRate = misfilteredStubRate + combinatorialStubRate; // combinatoricStubRate scales with the square of Nmb, while misfilteredStubRate scales linearly with Nmb

    curAvgTrue  += (trueStubRate - curAvgTrue)/(curCnt+1);
    curAvgInteresting += (highPtParticlesRate - curAvgInteresting)/(curCnt+1);
    curAvgMisfiltered  += (misfilteredStubRate - curAvgMisfiltered)/(curCnt+1);
    curAvgCombinatorial += (combinatorialStubRate - curAvgCombinatorial)/(curCnt+1);

    double curAvgFake = curAvgMisfiltered + curAvgCombinatorial;
    double curAvgTotal = curAvgTrue + curAvgFake;

    triggerFrequencyAverageTrue_[table][std::make_pair(row, col)] = curAvgTrue;            
    triggerFrequencyInterestingParticleTrue_[table][std::make_pair(row, col)] = curAvgInteresting;    
    triggerFrequencyAverageFake_[table][std::make_pair(row, col)] = curAvgFake;    
    triggerFrequencyAverageMisfiltered_[table][std::make_pair(row, col)] = curAvgMisfiltered;
    triggerFrequencyAverageCombinatorial_[table][std::make_pair(row, col)] = curAvgCombinatorial;

    int triggerDataHeaderBits  = module.numTriggerDataHeaderBits();
    int triggerDataPayloadBits = module.numTriggerDataPayloadBits();
    double triggerDataBandwidth = (triggerDataHeaderBits + curAvgTotal*triggerDataPayloadBits) / (bunchSpacingNs_); // GIGABIT/second
    triggerDataBandwidths_[table][std::make_pair(row, col)] = triggerDataBandwidth;
    triggerFrequenciesPerEvent[table][std::make_pair(row, col)] = curAvgTotal;


    //                currentTotalGraph->SetPoint(module->getRing()-1, module->getRing(), curAvgTotal*(1000/tracker.getBunchSpacingNs())*(100/module->getArea()));
    //                currentTrueGraph->SetPoint(module->getRing()-1, module->getRing(), curAvgTrue*(1000/tracker.getBunchSpacingNs())*(100/module->getArea()));

    currentTotalHisto->SetBinContent(col, curAvgTotal*(1000/bunchSpacingNs_)*(100/module.area()));
    currentTrueHisto->SetBinContent(col, curAvgTrue*(1000/bunchSpacingNs_)*(100/module.area()));

    triggerFrequencyTrueSummaries[table].setCell(row, col, curAvgTrue);
    triggerFrequencyInterestingSummaries[table].setCell(row, col, curAvgInteresting);
    triggerFrequencyFakeSummaries[table].setCell(row, col, curAvgFake);
    triggerFrequencyMisfilteredSummaries[table].setCell(row, col, curAvgMisfiltered);
    triggerFrequencyCombinatorialSummaries[table].setCell(row, col, curAvgCombinatorial);
    triggerRateSummaries[table].setCell(row, col, curAvgTotal);             
    triggerEfficiencySummaries[table].setCell(row, col, curAvgTrue/curAvgInteresting);                
    triggerPuritySummaries[table].setCell(row, col, curAvgTrue/(curAvgTrue+curAvgFake));                
    triggerDataBandwidthSummaries[table].setCell(row, col, triggerDataBandwidth);

    stripOccupancySummaries[table].setCell(row, col, module.stripOccupancyPerEvent()*nMB_*100);
    hitOccupancySummaries[table].setCell(row, col, module.hitOccupancyPerEvent()*nMB_*100);

  }

};


#endif

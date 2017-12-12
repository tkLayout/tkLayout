#ifndef BANDWIDTH_H
#define BANDWIDTH_H

#include <string>
#include <map>
#include <vector>
#include <utility>

#include "TH1.h"

#include "Tracker.hh"
#include "SimParms.hh"

#include "Visitor.hh"
#include "SummaryTable.hh"

class BandwidthVisitor : public ConstGeometryVisitor {
  TH1D &chanHitDistribution_, &bandwidthDistribution_, &bandwidthDistributionSparsified_;

  double nMB_;
public:
  BandwidthVisitor(TH1D& chanHitDistribution, TH1D& bandwidthDistribution, TH1D& bandwidthDistributionSparsified) :
      chanHitDistribution_(chanHitDistribution),
      bandwidthDistribution_(bandwidthDistribution),
      bandwidthDistributionSparsified_(bandwidthDistributionSparsified)
  {}

  void preVisit() {
    chanHitDistribution_.Reset();
    bandwidthDistribution_.Reset();
    bandwidthDistributionSparsified_.Reset();
    chanHitDistribution_.SetNameTitle("NHitChannels", "Number of hit channels;Hit Channels;Modules");
    bandwidthDistribution_.SetNameTitle("BandWidthDist", "Module Needed Bandwidth;Bandwidth (bps);Modules");
    bandwidthDistributionSparsified_.SetNameTitle("BandWidthDistSp", "Module Needed Bandwidth (sparsified);Bandwidth (bps);Modules");
    chanHitDistribution_.SetBins(200, 0., 400);
    bandwidthDistribution_.SetBins(100, 0., 6E+8);
    bandwidthDistributionSparsified_.SetBins(100, 0., 6E+8);
    bandwidthDistribution_.SetLineColor(kBlack);
    bandwidthDistributionSparsified_.SetLineColor(kRed);
  }

  void visit(const SimParms& sp) {
    nMB_ = sp.numMinBiasEvents();
  }

  void visit(const DetectorModule& m) {
    if (m.sensors().back().type() == SensorType::Strip) {
      for (auto s : m.sensors()) {
        double occupancy = m.hitOccupancyPerEvent();
        double hitChannels = occupancy * nMB_ * s.numChannels();
        chanHitDistribution_.Fill(hitChannels);
        int nChips = s.totalROCs();

        // Binary unsparsified (bps)
        bandwidthDistribution_.Fill((16*nChips + s.numChannels())*100E3);

        int spHdr = m.numSparsifiedHeaderBits();
        int spPay = m.numSparsifiedPayloadBits();      

        bandwidthDistributionSparsified_.Fill(((spHdr*nChips)+(hitChannels*spPay))*100E3);
      }
    }
  }
};


#endif

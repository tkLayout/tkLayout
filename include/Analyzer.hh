// File:   Analyzer.h
// Author: ndemaio
//
// Created on May 19, 2009, 1:58 PM
//

/**
 * @file Analyzer.h
 * @brief This class takes care of analysing the material budget
 */

#ifndef ANALYZER_H
#define ANALYZER_H

#define MY_RANDOM_SEED 0xcaffe

#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <ModuleCap.hh>
#include <InactiveElement.hh>
#include <InactiveSurfaces.hh>
#include <MaterialBudget.hh>
#include <TCanvas.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <Math/Vector3D.h>

#include <global_funcs.hh>

#include "TRandom3.h"
#include "Module.hh"
#include "AnalyzerVisitor.hh"
#include "Bag.hh"
#include "SummaryTable.hh"
#include "TagMaker.hh"
#include "Hit.hh"
#include "TrackNew.hh"

#include <TFile.h>
#include <TProfile.h>
#include <TF1.h>
#include <TAxis.h>
#include <TCanvas.h>

#include <AnalyzerTools.hh>

// Forward declaration
class TProfile;

namespace insur {

  /**
   * A warning that may occur during processing
   */
  static const std::string msg_module_warning = "Warning: tracker module with undefined subdetector type found.";

  /**
   * Two comparison functions for <i>std::pair<int, int></i> entries.
   */
  bool compareIntPairFirst(std::pair<int, int> p, std::pair<int, int> q);
  bool compareIntPairSecond(std::pair<int, int> p, std::pair<int, int> q);
  /**
   * @class Analyzer
   * @brief This class analyses the properties of a given <i>MaterialBudget</i> instance with respect to eta.
   *
   * It simulates a series of tracks that start at the origin (z = 0), maintain a fixed value of PI / 2 for phi and cover
   * an eta range from 0 to the maximal eta found in the provided geometry. Each volume hit by a track contributes
   * its radiation and interaction lengths to a grand total for that track. Those grand totals, recorded by eta, are
   * stored in a series of histograms that give a complete profile of the expected interaction of the tracker itself with
   * the particles that pass through it.
   */

  typedef std::map<std::pair<std::string, int>, TH1D*> StubRateHistos;
  typedef std::vector<Module*> ModuleVector;
  typedef std::vector<Layer*> LayerVector;
  typedef TriggerProcessorBandwidthVisitor::ModuleConnectionMap ModuleConnectionMap;
  typedef TriggerProcessorBandwidthVisitor::TriggerSectorMap TriggerSectorMap;

  // TODO:
  // Move this to track.hh?
  typedef std::vector<Track> TrackCollection;
  //typedef double TrackCollectionKey;
  typedef std::map<int, TrackNewCollection> TrackNewCollectionMap;


  class Analyzer : private AnalyzerTools {
  public:
    Analyzer();
    virtual ~Analyzer() {}
    std::map<std::string, TH1D*>& getHistoActiveComponentsR() { return rComponents; }
    std::map<std::string, TH1D*>& getHistoActiveComponentsI() { return iComponents; }
    std::map<std::string, TH1D*>& getHistoBeamPipeR() { return rComponentsBeamPipe; }
    std::map<std::string, TH1D*>& getHistoBeamPipeI() { return iComponentsBeamPipe; }
    std::map<std::string, TH1D*>& getHistoPixelIntersticeR() { return rComponentsPixelInterstice; }
    std::map<std::string, TH1D*>& getHistoPixelIntersticeI() { return iComponentsPixelInterstice; }
    std::map<std::string, TH1D*>& getHistoPixelTrackingVolumeR() { return rComponentsPixelTrackingVolume; }
    std::map<std::string, TH1D*>& getHistoPixelTrackingVolumeI() { return iComponentsPixelTrackingVolume; }
    std::map<std::string, TH1D*>& getHistoIntersticeR() { return rComponentsInterstice; }
    std::map<std::string, TH1D*>& getHistoIntersticeI() { return iComponentsInterstice; }
    std::map<std::string, TH1D*>& getHistoOuterTrackingVolumeR() { return rComponentsOuterTrackingVolume; }
    std::map<std::string, TH1D*>& getHistoOuterTrackingVolumeI() { return iComponentsOuterTrackingVolume; }
    //std::map<std::string, TH1D*>& getHistoTotalTrackingVolumeR() { return rComponentsTotalTrackingVolume; }
    //std::map<std::string, TH1D*>& getHistoTotalTrackingVolumeI() { return iComponentsTotalTrackingVolume; }
    //std::vector<std::string>& getComponentsTrackingVolume() { return componentsTotalTrackingVolumeOrder; }
    TH1D& getHistoModulesBarrelsR() { return ractivebarrel; }
    TH1D& getHistoModulesBarrelsI() { return iactivebarrel; }
    TH1D& getHistoModulesEndcapsR() { return ractiveendcap; }
    TH1D& getHistoModulesEndcapsI() { return iactiveendcap; }
    TH1D& getHistoServicesBarrelsR() { return rserfbarrel; }
    TH1D& getHistoServicesBarrelsI() { return iserfbarrel; }
    TH1D& getHistoServicesEndcapsR() { return rserfendcap; }
    TH1D& getHistoServicesEndcapsI() { return iserfendcap; }
    TH1D& getHistoSupportsBarrelsR() { return rlazybarrel; }
    TH1D& getHistoSupportsBarrelsI() { return ilazybarrel; }
    TH1D& getHistoSupportsEndcapsR() { return rlazyendcap; }
    TH1D& getHistoSupportsEndcapsI() { return ilazyendcap; }
    TH1D& getHistoSupportsBarrelTubesR() { return rlazybtube; }
    TH1D& getHistoSupportsBarrelTubesI() { return ilazybtube; }
    TH1D& getHistoSupportsTubesR() { return rlazytube; }
    TH1D& getHistoSupportsTubesI() { return ilazytube; }
    TH1D& getHistoSupportsUserDefinedR() { return rlazyuserdef; }
    TH1D& getHistoSupportsUserDefinedI() { return ilazyuserdef; }
    TH1D& getHistoBarrelsAllR() { return rbarrelall; }
    TH1D& getHistoBarrelsAllI() { return ibarrelall; }
    TH1D& getHistoEndcapsAllR() { return rendcapall; }
    TH1D& getHistoEndcapsAllI() { return iendcapall; }
    TH1D& getHistoModulesAllR() { return ractiveall; }
    TH1D& getHistoModulesAllI() { return iactiveall; }
    TH1D& getHistoServicesAllR() { return rserfall; }
    TH1D& getHistoServicesAllI() { return iserfall; }
    TH1D& getHistoSupportsAllR() { return rlazyall; }
    TH1D& getHistoSupportsAllI() { return ilazyall; }
    std::map<std::string, TH1D*>& getHistoServicesDetailsR() { return rComponentsServicesDetails; }
    std::map<std::string, TH1D*>& getHistoServicesDetailsI() { return iComponentsServicesDetails; }
    TH1D& getHistoGlobalR() { return rglobal; }
    TH1D& getHistoGlobalI() { return iglobal; }
    TH2D& getHistoIsoR() { return isor; }
    TH2D& getHistoIsoI() { return isoi; }
    TH2D& getHistoMapRadiation();
    TH2D& getHistoMapInteraction();
    TH1D& getHistoOptimalSpacing(bool actualWindow);
    //std::vector<Track>& getTracks() { return tv; } // useless ?! remove !
    std::map<int, TGraph>& getRhoGraphs(bool ideal, bool isTrigger);
    std::map<int, TGraph>& getPhiGraphs(bool ideal, bool isTrigger);
    std::map<int, TGraph>& getDGraphs(bool ideal, bool isTrigger);
    std::map<int, TGraph>& getCtgThetaGraphs(bool ideal, bool isTrigger);
    std::map<int, TGraph>& getZ0Graphs(bool ideal, bool isTrigger);
    std::map<int, TGraph>& getPGraphs(bool ideal, bool isTrigger);
    std::map<int, TGraph>& getLGraphs(bool ideal, bool isTrigger);
    std::map<int, TGraph>& getBetaGraphs(bool ideal, bool isTrigger);
    std::map<int, TGraph>& getOmegaGraphs(bool ideal, bool isTrigger);
    GraphBag& getGraphBag() { return myGraphBag; }
    mapBag& getMapBag() { return myMapBag; }
    profileBag& getProfileBag() { return myProfileBag; }
    std::map<int, TGraphErrors>& getSpacingTuningGraphs() { return spacingTuningGraphs; }
    std::map<int, TGraphErrors>& getSpacingTuningGraphsBad() { return spacingTuningGraphsBad; }
    TH1D& getSpacingTuningFrame() { return spacingTuningFrame; }
    const double& getTriggerRangeLowLimit(const std::string& typeName ) { return triggerRangeLowLimit[typeName] ; }
    const double& getTriggerRangeHighLimit(const std::string& typeName ) { return triggerRangeHighLimit[typeName] ; }
    /*virtual*/ void analyzeMaterialBudget(MaterialBudget& mb, const std::vector<double>& momenta, int etaSteps = 50, MaterialBudget* pm = NULL);
    void computeTriggerProcessorsBandwidth(Tracker& tracker);
    void analyzeTaggedTracking(MaterialBudget& mb,
                               const std::vector<double>& momenta,
                               const std::vector<double>& triggerMomenta,
                               const std::vector<double>& thresholdProbabilities,
                               bool isPixel,
                               bool& debugResolution,
                               int etaSteps = 50,
                               MaterialBudget* pm = nullptr);
    bool checkFile(const std::string& fileName, const std::string& filePath);
    bool isTripletFromDifLayers(TrackNew& track, int iHit, bool propagOutIn);
    bool analyzePatterReco(MaterialBudget& mb, mainConfigHandler& mainConfig, int etaSteps = 50, MaterialBudget* pm = nullptr);
    std::vector<TProfile*> hisPatternRecoInOutPt;//! InOut approach - tracker: Bkg contamination probability accumulated across eta for set of pT
    std::vector<TProfile*> hisPatternRecoInOutP; //! InOut approach - inner tracker: Bkg contamination probability accumulated across eta for set of pT
    std::map<std::string, std::vector<TProfile*>> hisPtHitDProjInOut;     //!< InOut approach: D0 projection @ ith+3 measurement plane at given eta for set of pt
    std::map<std::string, std::vector<TProfile*>> hisPHitDProjInOut;      //!< InOut approach: D0 projection @ ith+3 measurement plane at given eta for set of p
    std::map<std::string, std::vector<TProfile*>> hisPtHitZProjInOut;     //!< InOut approach: Z0 projection @ ith+3 measurement plane at given eta for set of pt
    std::map<std::string, std::vector<TProfile*>> hisPHitZProjInOut;      //!< InOut approach: Z0 projection @ ith+3 measurement plane at given eta for set of p
    std::map<std::string, std::vector<TProfile*>> hisPtHitProbContamInOut;//!< InOut approach: Calculated occupancy from flux & pile-up @ ith+3 measurement plane at given eta for set of pt
    std::map<std::string, std::vector<TProfile*>> hisPHitProbContamInOut; //!< InOut approach: Calculated occupancy from flux & pile-up @ ith+3 measurement plane at given eta for set of p
    std::vector<TProfile*> hisPatternRecoOutInPt;//! OutIn approach - tracker: Bkg contamination probability accumulated across eta for set of pT
    std::vector<TProfile*> hisPatternRecoOutInP; //! OutIn approach - inner tracker: Bkg contamination probability accumulated across eta for set of pT
    std::map<std::string, std::vector<TProfile*>> hisPtHitDProjOutIn;     //!< OutIn approach: D0 projection @ ith+3 measurement plane at given eta for set of pt
    std::map<std::string, std::vector<TProfile*>> hisPHitDProjOutIn;      //!< OutIn approach: D0 projection @ ith+3 measurement plane at given eta for set of p
    std::map<std::string, std::vector<TProfile*>> hisPtHitZProjOutIn;     //!< OutIn approach: Z0 projection @ ith+3 measurement plane at given eta for set of pt
    std::map<std::string, std::vector<TProfile*>> hisPHitZProjOutIn;      //!< OutIn approach: Z0 projection @ ith+3 measurement plane at given eta for set of p
    std::map<std::string, std::vector<TProfile*>> hisPtHitProbContamOutIn;//!< OutIn approach: Calculated occupancy from flux & pile-up @ ith+3 measurement plane at given eta for set of pt
    std::map<std::string, std::vector<TProfile*>> hisPHitProbContamOutIn; //!< OutIn approach: Calculated occupancy from flux & pile-up @ ith+3 measurement plane at given eta for set of p
    
    virtual void analyzeTriggerEfficiency(Tracker& tracker,
                                          const std::vector<double>& triggerMomenta,
                                          const std::vector<double>& thresholdProbabilities,
                                          int etaSteps = 50);
    void createTriggerDistanceTuningPlots(Tracker& tracker, const std::vector<double>& triggerMomenta);
    void analyzeGeometry(Tracker& tracker, int nTracks = 1000);
    void computeBandwidth(Tracker& tracker);
    void computeTriggerFrequency(Tracker& tracker);
    void analyzePower(Tracker& tracker);
    void createGeometryLite(Tracker& tracker);
    TH2D& getMapPhiEta() { return mapPhiEta; }
    TCanvas& getEtaProfileCanvas() {return etaProfileCanvas; }
    TH1D& getHitDistribution() {return hitDistribution; }
    TProfile& getTotalEtaProfile() {return totalEtaProfile; }
    TProfile& getTotalEtaProfileSensors() {return totalEtaProfileSensors; }
    TProfile& getTotalEtaProfileStubs() {return totalEtaProfileStubs; }
    TProfile& getTotalEtaProfileLayers() {return totalEtaProfileLayers; }
    TGraph& getPowerDensity() {return powerDensity;};
    std::vector<TProfile>& getTypeEtaProfiles() {return typeEtaProfile; }
    std::vector<TProfile>& getTypeEtaProfilesSensors() {return typeEtaProfileSensors; }
    std::vector<TProfile>& getTypeEtaProfilesStubs() {return typeEtaProfileStubs; }
    std::map<std::string, TProfile>& getLayerEtaCoverageProfiles() {return layerEtaCoverageProfile;}
    std::map<std::string, TProfile>& getLayerEtaCoverageProfilesStubs() {return layerEtaCoverageProfileStubs; }
    std::map<std::string, std::map<std::string, TH1I*>>& getStubEfficiencyCoverageProfiles() { return stubEfficiencyCoverageProfiles_; } // map of maps: inner map has momenta as keys
    std::vector<TObject> getSavingVector();
    TCanvas* getGeomLite() {if (geomLiteCreated) return geomLite; else return NULL; };
    TCanvas* getGeomLiteXY() {if (geomLiteXYCreated) return geomLiteXY; else return NULL; };
    TCanvas* getGeomLiteYZ() {if (geomLiteYZCreated) return geomLiteYZ; else return NULL; };
    TCanvas* getGeomLiteEC() {if (geomLiteECCreated) return geomLiteEC; else return NULL; };
    TH1D& getChanHitDistribution() { return chanHitDistribution; };
    TH1D& getBandwidthDistribution() { return bandwidthDistribution; };
    TH1D& getBandwidthDistributionSparsified() { return bandwidthDistributionSparsified; }
    TH1I& getModuleConnectionsDistribution() { return moduleConnectionsDistribution; }
    const ModuleConnectionMap& getModuleConnectionMap() const { return moduleConnections_; }
    int getGeometryTracksUsed() {return geometryTracksUsed; }
    int getMaterialTracksUsed() {return materialTracksUsed; }
    // Hadrons
    TGraph& getHadronTotalHitsGraph() {return hadronTotalHitsGraph;};
    TGraph& getHadronAverageHitsGraph() {return hadronAverageHitsGraph;};

    StubRateHistos& getTotalStubRateHistos() { return totalStubRateHistos_; }
    StubRateHistos& getTrueStubRateHistos() { return trueStubRateHistos_; }

    std::vector<double>& getHadronNeededHitsFraction() {return hadronNeededHitsFraction;};
    std::vector<TGraph>& getHadronGoodTracksFraction() { return hadronGoodTracksFraction; };


    static std::vector<double> average(TGraph& myGraph, std::vector<double> cuts);

    static const double ZeroHitsRequired;
    static const double OneHitRequired;

    void computeWeightSummary(MaterialBudget& mb);
    std::map<std::string, SummaryTable>& getBarrelWeightSummary() { return barrelWeights;};
    std::map<std::string, SummaryTable>& getEndcapWeightSummary() { return endcapWeights;};
    std::map<std::string, SummaryTable>& getBarrelWeightComponentSummary() { return barrelComponentWeights;};
    std::map<std::string, SummaryTable>& getEndcapWeightComponentSummary() { return endcapComponentWeights;};
    std::map<std::string, double>& getTypeWeigth() { return typeWeight; };
    std::map<std::string, double>& getTagWeigth() { return tagWeight; };
    std::map<std::string, TH2D>& getParametrizedResolutionLocalXBarrelMap() {return parametrizedResolutionLocalXBarrelMap; }
    std::map<std::string, TH2D>& getParametrizedResolutionLocalXEndcapsMap() { return parametrizedResolutionLocalXEndcapsMap; }
    std::map<std::string, TH2D>& getParametrizedResolutionLocalYBarrelMap() { return parametrizedResolutionLocalYBarrelMap; }
    std::map<std::string, TH2D>& getParametrizedResolutionLocalYEndcapsMap() { return parametrizedResolutionLocalYEndcapsMap; }
    std::map<std::string, TH1D>& getParametrizedResolutionLocalXBarrelDistribution() {return parametrizedResolutionLocalXBarrelDistribution; }
    std::map<std::string, TH1D>& getParametrizedResolutionLocalXEndcapsDistribution() { return parametrizedResolutionLocalXEndcapsDistribution; }
    std::map<std::string, TH1D>& getParametrizedResolutionLocalYBarrelDistribution() { return parametrizedResolutionLocalYBarrelDistribution; }
    std::map<std::string, TH1D>& getParametrizedResolutionLocalYEndcapsDistribution() { return parametrizedResolutionLocalYEndcapsDistribution; }
    std::map<std::string, SummaryTable>& getTriggerFrequencyTrueSummaries() { return triggerFrequencyTrueSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerFrequencyInterestingSummaries() { return triggerFrequencyInterestingSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerFrequencyFakeSummaries() { return triggerFrequencyFakeSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerFrequencyMisfilteredSummaries() { return triggerFrequencyMisfilteredSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerFrequencyCombinatorialSummaries() { return triggerFrequencyCombinatorialSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerRateSummaries() { return triggerRateSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerEfficiencySummaries() { return triggerEfficiencySummaries_; }
    std::map<std::string, SummaryTable>& getTriggerPuritySummaries() { return triggerPuritySummaries_; }
    std::map<std::string, SummaryTable>& getTriggerDataBandwidthSummaries() { return triggerDataBandwidthSummaries_; }
    
    double getTriggerPetalCrossoverR() const { return triggerPetalCrossoverR_; }
    const std::pair<Circle, Circle>& getSampleTriggerPetal() const { return sampleTriggerPetal_; }

    std::map<std::string, SummaryTable>& getStripOccupancySummaries() { return stripOccupancySummaries_; }
    std::map<std::string, SummaryTable>& getHitOccupancySummaries() { return hitOccupancySummaries_; }

    SummaryTable& getProcessorConnectionSummary() { return processorConnectionSummary_; }
    SummaryTable& getProcessorCommonConnectionSummary() { return processorCommonConnectionSummary_; }
    TH2I& getProcessorCommonConnectionMap() { return processorCommonConnectionMap_; }
    std::map<std::string, SummaryTable>& getModuleConnectionSummaries() { return moduleConnectionSummaries_; }
    SummaryTable& getProcessorInboundBandwidthSummary() { return processorInboundBandwidthSummary_; }
    SummaryTable& getProcessorInboundStubPerEventSummary() { return processorInboundStubPerEventSummary_; }

    const TriggerSectorMap& getTriggerSectorMap() const { return triggerSectorMap_; }

    inline double getEtaMaxMaterial() const { return insur::geom_max_eta_coverage;}
    inline double getEtaMaxGeometry() const { return insur::geom_max_eta_coverage;}
    inline double getEtaMaxTracker()  const { return insur::geom_max_eta_coverage;}
    inline double getEtaMaxTrigger()  const { return insur::geom_max_eta_coverage;}

    const std::string & getBillOfMaterials() { return billOfMaterials_ ; }
  protected:
    /**
     * @struct Cell
     * @brief This struct contains the cumulative radiation and interaction lengths over the area described by r and z.
     * @param rlength The cumulative radiation length
     * @param ilength The cumulative interaction length
     * @param rmin The minimal radius value of the cell
     * @param rmax The maximal radius value of the cell
     * @param etamin The minimal eta value of the cell
     * @param etamax The maximal eta value of the cell
     */
    struct Cell { double rlength; double ilength; double rmin; double rmax; double etamin; double etamax; };
    std::vector<std::vector<Cell> > cells;
    TH1D ractivebarrel, ractiveendcap, rserfbarrel, rserfendcap, rlazybarrel, rlazyendcap, rlazybtube, rlazytube, rlazyuserdef;
    TH1D iactivebarrel, iactiveendcap, iserfbarrel, iserfendcap, ilazybarrel, ilazyendcap, ilazybtube, ilazytube, ilazyuserdef;
    TH1D rbarrelall, rendcapall, ractiveall, rserfall, rlazyall;
    TH1D ibarrelall, iendcapall, iactiveall, iserfall, ilazyall;
    std::map<std::string, TH1D*> rComponentsServicesDetails, iComponentsServicesDetails;
    TH1D rglobal, iglobal;

    std::map<std::string, TH1D*> rComponents, iComponents;
    std::map<std::string, TH1D*> rComponentsBeamPipe, iComponentsBeamPipe;
    std::map<std::string, TH1D*> rComponentsPixelInterstice, iComponentsPixelInterstice;
    std::map<std::string, TH1D*> rComponentsPixelTrackingVolume, iComponentsPixelTrackingVolume;
    std::map<std::string, TH1D*> rComponentsInterstice, iComponentsInterstice;
    std::map<std::string, TH1D*> rComponentsOuterTrackingVolume, iComponentsOuterTrackingVolume;
    TH2D isor, isoi;
    TH2D mapRadiation, mapInteraction;
    TH2I mapRadiationCount, mapInteractionCount;
    TH2D mapRadiationCalib, mapInteractionCalib;
    TH2D mapPhiEta;
    TCanvas etaProfileCanvas;
    TCanvas* geomLite; bool geomLiteCreated;
    TCanvas* geomLiteXY; bool geomLiteXYCreated;
    TCanvas* geomLiteYZ; bool geomLiteYZCreated;
    TCanvas* geomLiteEC; bool geomLiteECCreated;
    TH1D chanHitDistribution;
    TH1D bandwidthDistribution;
    TH1D bandwidthDistributionSparsified;
    TH1D optimalSpacingDistribution;
    TH1D optimalSpacingDistributionAW;

    TH1I moduleConnectionsDistribution;

    std::map<std::string, SummaryTable> barrelWeights;
    std::map<std::string, SummaryTable> endcapWeights;
    std::map<std::string, SummaryTable> barrelComponentWeights;
    std::map<std::string, SummaryTable> endcapComponentWeights;
    std::map<std::string, double> typeWeight;
    std::map<std::string, double> tagWeight;

    std::map<std::string, TH2D> parametrizedResolutionLocalXBarrelMap;
    std::map<std::string, TH2D> parametrizedResolutionLocalXEndcapsMap;
    std::map<std::string, TH2D> parametrizedResolutionLocalYBarrelMap;
    std::map<std::string, TH2D> parametrizedResolutionLocalYEndcapsMap; 
    std::map<std::string, TH1D> parametrizedResolutionLocalXBarrelDistribution;
    std::map<std::string, TH1D> parametrizedResolutionLocalXEndcapsDistribution;
    std::map<std::string, TH1D> parametrizedResolutionLocalYBarrelDistribution;
    std::map<std::string, TH1D> parametrizedResolutionLocalYEndcapsDistribution;

    std::map<std::string, std::map<std::pair<int, int>, double> > triggerDataBandwidths_;
    std::map<std::string, std::map<std::pair<int, int>, double> > triggerFrequenciesPerEvent_;
    std::map<std::string, SummaryTable> triggerFrequencyTrueSummaries_, triggerFrequencyFakeSummaries_, triggerFrequencyMisfilteredSummaries_, triggerFrequencyCombinatorialSummaries_, triggerFrequencyInterestingSummaries_;
    std::map<std::string, SummaryTable> triggerRateSummaries_, triggerEfficiencySummaries_, triggerPuritySummaries_;
    std::map<std::string, SummaryTable> triggerDataBandwidthSummaries_;

    std::map<std::string, SummaryTable> stripOccupancySummaries_;
    std::map<std::string, SummaryTable> hitOccupancySummaries_;

    SummaryTable processorConnectionSummary_;
    SummaryTable processorCommonConnectionSummary_;
    TH2I processorCommonConnectionMap_;
    std::map<std::string, SummaryTable> moduleConnectionSummaries_;
    SummaryTable processorInboundBandwidthSummary_;
    SummaryTable processorInboundStubPerEventSummary_;

    double triggerPetalCrossoverR_;
    std::pair<Circle, Circle> sampleTriggerPetal_;

    ModuleConnectionMap moduleConnections_;
    TriggerSectorMap triggerSectorMap_;
    
    TH1D hitDistribution;
    GraphBag myGraphBag;
    mapBag myMapBag;
    profileBag myProfileBag;
    std::map<int, TGraphErrors> spacingTuningGraphs; // TODO: find a way to communicate the limits, not their plots!
    std::map<int, TGraphErrors> spacingTuningGraphsBad; // TODO: find a way to communicate the limits, not their plots!
    ModuleOptimalSpacings moduleOptimalSpacings;
    TH1D spacingTuningFrame;
    std::map<std::string, double> triggerRangeLowLimit;
    std::map<std::string, double> triggerRangeHighLimit;

    // Hadrons
    TGraph hadronTotalHitsGraph;
    TGraph hadronAverageHitsGraph;
    std::vector<double> hadronNeededHitsFraction;
    std::vector<TGraph> hadronGoodTracksFraction;


    StubRateHistos totalStubRateHistos_;
    StubRateHistos trueStubRateHistos_;

    TGraph powerDensity;
    TProfile totalEtaProfile, totalEtaProfileSensors, totalEtaProfileStubs, totalEtaProfileLayers;
    std::vector<TProfile> typeEtaProfile, typeEtaProfileSensors, typeEtaProfileStubs;
    std::map<std::string, TProfile> layerEtaCoverageProfile, layerEtaCoverageProfileStubs;

    std::map<std::string, std::map<std::string, TH1I*>> stubEfficiencyCoverageProfiles_;

    std::vector<TObject> savingGeometryV; // Vector of ROOT objects to be saved
    std::vector<TObject> savingMaterialV; // Vector of ROOT objects to be saved

    Material findAllHits(MaterialBudget& mb, MaterialBudget* pm, TrackNew& track);


    void computeDetailedWeights(std::vector<std::vector<ModuleCap> >& tracker, std::map<std::string, SummaryTable>& weightTables, bool byMaterial);
    virtual Material analyzeModules(std::vector<std::vector<ModuleCap> >& tr, double eta, double theta, double phi, Track& t, 
                                    std::map<std::string, Material>& sumComponentsRI, bool isPixel = false);

    int findHitsModules(Tracker& tracker, double z0, double eta, double theta, double phi, Track& t);

    virtual Material findHitsModules(std::vector<std::vector<ModuleCap> >& tr, TrackNew& t, bool isPixel = false);
    virtual Material findHitsModuleLayer(std::vector<ModuleCap>& layer, TrackNew& t, bool isPixel = false);

    virtual Material findModuleLayerRI(std::vector<ModuleCap>& layer, double eta, double theta, double phi, Track& t, 
                                       std::map<std::string, Material>& sumComponentsRI, bool isPixel = false);
    virtual Material analyzeInactiveSurfaces(std::vector<InactiveElement>& elements, double eta, double theta, 
                                             Track& t, std::map<std::string, Material>& sumServicesComponentsRI, MaterialProperties::Category cat = MaterialProperties::no_cat, bool isPixel = false);
    virtual Material findHitsInactiveSurfaces(std::vector<InactiveElement>& elements, TrackNew& t, bool isPixel = false);

    void clearGraphsPt(int graphAttributes, const std::string& aTag);
    void clearGraphsP(int graphAttributes, const std::string& aTag);
    void calculateGraphsConstPt(const int& aMomentum,
                                const TrackNewCollection& aTrackCollection,
                                int graphAttributes,
                                const string& graphTag);
    void calculateGraphsConstP(const int& aMomentum,
                               const TrackNewCollection& aTrackCollection,
                               int graphAttributes,
                               const string& graphTag);
    void calculateParametrizedResolutionPlots(std::map<std::string, TrackNewCollectionMap>& taggedTrackPtCollectionMap);
    void fillTriggerEfficiencyGraphs(const Tracker& tracker,
                                     const std::vector<double>& triggerMomenta,
                                     const std::vector<Track>& trackVector);
    void fillTriggerPerformanceMaps(Tracker& tracker);
    //void fillPowerMap(Tracker& tracker);
    void clearMaterialBudgetHistograms();
    void prepareTriggerPerformanceHistograms(const int& nTracks, const double& etaMax, const vector<double>& triggerMomenta, const vector<double>& thresholdProbabilities);
    void prepareTriggerProcessorHistograms();
    void clearGeometryHistograms();
    void clearCells();
    void setHistogramBinsBoundaries(int bins, double min, double max);
    void setCellBoundaries(int bins, double minr, double maxr, double minz, double maxz);
    void fillCell(double r, double eta, double theta, Material mat);
    void fillMapRT(const double& r, const double& theta, const Material& mat);
    void fillMapRZ(const double& r, const double& z, const Material& mat);
    void transformEtaToZ();
    double findXThreshold(const TProfile& aProfile, const double& yThreshold, const bool& goForward );
    std::pair<double, double> computeMinMaxTracksEta(const Tracker& t) const;
  private:
    // A random number generator
    TRandom3 myDice; 
    int findCellIndexR(double r);
    int findCellIndexEta(double eta);
    int createResetCounters(Tracker& tracker, std::map <std::string, int> &modTypes);
    std::pair <XYZVector, double > shootDirection(double minEta, double maxEta);
    std::vector<std::pair<Module*, HitType>> trackHit(const XYZVector& origin, const XYZVector& direction, Tracker::Modules& properModules);
    void resetTypeCounter(std::map<std::string, int> &modTypes);
    double diffclock(clock_t clock1, clock_t clock2);
    Color_t colorPicker(std::string);
    std::map<std::string, Color_t> colorPickMap;
    Color_t lastPickedColor;
    int geometryTracksUsed;
    int materialTracksUsed;
    void fillAvailableSpacing(Tracker& tracker, std::vector<double>& spacingOptions);
    static constexpr double maximum_n_planes = 13.;

    bool isModuleInEtaSector(const Tracker& tracker, const Module* module, int etaSector) const;
    bool isModuleInPhiSector(const Tracker& tracker, const Module* module, int phiSector) const;

    static int bsCounter;
    
    std::string billOfMaterials_;
  };
}
#endif  /* _ANALYZER_H */


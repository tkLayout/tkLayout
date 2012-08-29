// 
// File:   Analyzer.h
// Author: ndemaio
//
// Created on May 19, 2009, 1:58 PM
//

/**
 * @file Analyzer.h
 * @brief This class takes care of analysing the material budget
 */

#ifndef _ANALYZER_H
#define _ANALYZER_H

#define MY_RANDOM_SEED 0xcaffe

#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <hit.hh>
#include <module.hh>
#include <ModuleCap.h>
#include <InactiveElement.h>
#include <InactiveSurfaces.h>
#include <MaterialBudget.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <global_funcs.h>

#include "TRandom3.h"

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
   * @class SummaryTable
   * @brief A generic object to build summary tables
   */

  class SummaryTable {
  public:
    SummaryTable() : numRows_(0), numColumns_(0), rowOffset_(0), columnOffset_(0), precision_(-1), summaryCellPosition_(0,0), summaryLabelPosition_(0,0) {};
    void setHeader(std::string rowHeader, std::string columnHeader, int rowOffset = 0, int columnOffset = 0) { // has to be called before filling the table with content or the row and column numbering will not be correctly set
      rowOffset_ = rowOffset; columnOffset_ = columnOffset;
      summaryTable[make_pair(0,0)] = columnHeader + " &rarr;<br>" + rowHeader + " &darr;";
    }
    void setPrecision(int precision) { precision_ = precision; } // has to be called before filling the table or conversions from floating point won't have the desired precision

    template<typename T> void setCell(const int row, const int column, const T& content) { setCell(row, column, any2str(content, precision_)); }
    template<typename T> void setSummaryCell(std::string label, const T& content) { setSummaryCell(label, any2str(content, precision_)); }

    std::string getCell(int row, int column) { return summaryTable[make_pair(row,column)];} // this actually alters the map if the cell's not there = DANGEROUS

    bool hasCell(int row, int column) const { return summaryTable.count(make_pair(row,column)); }  // tests whether a cell has already been inserted = SAFE
    bool hasSummaryCell() const { return summaryCellPosition_ > std::make_pair(0, 0); }

    std::map<std::pair<int, int>, std::string>& getContent() { return summaryTable; }

    void clear() { summaryTable.clear(); }
  private:
    std::map<std::pair<int, int>, std::string> summaryTable;
    int numRows_, numColumns_;
    int rowOffset_, columnOffset_; // from which number rows and columns headers should start
    int precision_; // precision to convert floating point numbers with
    std::pair<int, int> summaryCellPosition_, summaryLabelPosition_;
  };

  /**
   * @class graphBag
   * @brief A bag of graphs sorted by variable, scope and track's pt
   */
  class graphBag {
  public:
    static const double Triggerable;
    static const int RhoGraph;
    static const int PhiGraph;
    static const int DGraph;
    static const int CtgthetaGraph;
    static const int Z0Graph;
    static const int PGraph;
    static const int TriggeredGraph;
    static const int IdealGraph;
    static const int RealGraph;
    static const int TriggerGraph;
    static const int StandardGraph;
    std::map<double, TGraph>& getGraphs(const int& attribute);
    int clearTriggerGraphs();
    int clearStandardGraphs();
    static int buildAttribute(bool ideal, bool isTrigger);
    //static std::pair<double, double> splitMomenta(double momentum);
    //static double joinMomenta(double momentum1, double momentum2);
  private:
    std::map<int, std::map<double, TGraph> > graphMap_;
    int clearGraphs(const int& attributeMask);
  };

  /**
   * @class mapBag
   * @brief A bag of graphs sorted by variable, scope and track's pt
   */
  class mapBag {
  public:
    static const int efficiencyMap;
    static const int thresholdMap;
    static const int thicknessMap;
    static const int windowMap;
    static const int suggestedSpacingMap;
    static const int suggestedSpacingMapAW;
    static const int nominalCutMap;
    static const int irradiatedPowerConsumptionMap;
    static const int totalPowerConsumptionMap;
    static const int moduleConnectionEtaMap;
    static const int moduleConnectionPhiMap;
    static const int moduleConnectionEndcapPhiMap;
    static const double dummyMomentum;
    std::map<double, TH2D>& getMaps(const int& attribute);
    int clearMaps(const int& attributeMask);
  private:
    std::map<int, std::map<double, TH2D> > mapMap_;
  };

  /**
   * @class profileBag
   * @brief A bag of TProfiles sorted by a double variable and scope
   */
  class profileBag {
  public:
    static const double Triggerable;
    static const int TriggeredProfile;
    static const int TriggerProfile;
    static const int TriggeredFractionProfile;
    static const std::string TriggerProfileName;
    static const std::string TriggerProfileNameWindow;
    static const std::string TurnOnCurveName;
    std::map<double, TProfile>& getProfiles(const int& attribute);
    int clearTriggerProfiles();
    int clearTriggerNamedProfiles();
    std::map<double, TProfile>& getNamedProfiles(const std::string& name);
    std::vector<std::string> getProfileNames(const std::string& name);
  private:
    int clearProfiles(const int& attributeMask);
    int clearNamedProfiles(const std::string& name);
    std::map<int, std::map<double, TProfile> > profileMap_;
    std::map<std::string, std::map<double, TProfile> > namedProfileMap_;
  };     


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

  class Analyzer {
  public:
    Analyzer();
    virtual ~Analyzer() {}
    std::map<std::string, TH1D*>& getHistoActiveComponentsR() { return rComponents; }
    std::map<std::string, TH1D*>& getHistoActiveComponentsI() { return iComponents; }
    TH1D& getHistoModulesBarrelsR() { return ractivebarrel; }
    TH1D& getHistoModulesBarrelsI() { return iactivebarrel; }
    TH1D& getHistoModulesEndcapsR() {return ractiveendcap; }
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
    TH1D& getHistoExtraServicesR() { return rextraservices; }
    TH1D& getHistoExtraServicesI() { return iextraservices; }
    TH1D& getHistoExtraSupportsR() { return rextrasupports; }
    TH1D& getHistoExtraSupportsI() { return iextrasupports; }
    TH1D& getHistoGlobalR() { return rglobal; }
    TH1D& getHistoGlobalI() { return iglobal; }
    TH2D& getHistoIsoR() { return isor; }
    TH2D& getHistoIsoI() { return isoi; }
    TH2D& getHistoMapRadiation();
    TH2D& getHistoMapInteraction();
    TH1D& getHistoOptimalSpacing(bool actualWindow);
    //std::vector<Track>& getTracks() { return tv; } // useless ?! remove !
    std::map<double, TGraph>& getRhoGraphs(bool ideal, bool isTrigger);
    std::map<double, TGraph>& getPhiGraphs(bool ideal, bool isTrigger);
    std::map<double, TGraph>& getDGraphs(bool ideal, bool isTrigger);
    std::map<double, TGraph>& getCtgThetaGraphs(bool ideal, bool isTrigger);
    std::map<double, TGraph>& getZ0Graphs(bool ideal, bool isTrigger);
    std::map<double, TGraph>& getPGraphs(bool ideal, bool isTrigger);
    graphBag& getGraphBag() { return myGraphBag; }
    mapBag& getMapBag() { return myMapBag; }
    profileBag& getProfileBag() { return myProfileBag; }
    std::map<int, TGraphErrors>& getSpacingTuningGraphs() { return spacingTuningGraphs; }
    std::map<int, TGraphErrors>& getSpacingTuningGraphsBad() { return spacingTuningGraphsBad; }
    TH1D& getSpacingTuningFrame() { return spacingTuningFrame; }
    const double& getTriggerRangeLowLimit(const std::string& typeName ) { return triggerRangeLowLimit[typeName] ; }
    const double& getTriggerRangeHighLimit(const std::string& typeName ) { return triggerRangeHighLimit[typeName] ; }
    virtual void analyzeMaterialBudget(MaterialBudget& mb, const std::vector<double>& momenta, int etaSteps = 50, MaterialBudget* pm = NULL, bool computeResolution = false);
    //virtual void analyzeMaterialBudgetTrigger(MaterialBudget& mb, std::vector<double>& momenta, int etaSteps = 50, MaterialBudget* pm = NULL);
    void computeTriggerProcessorsBandwidth(Tracker& tracker);
    virtual void analyzeTrigger(MaterialBudget& mb,
                                const std::vector<double>& momenta,
                                const std::vector<double>& triggerMomenta,
                                const std::vector<double>& thresholdProbabilities,
                                int etaSteps = 50, MaterialBudget* pm = NULL);
    virtual void analyzeTriggerEfficiency(Tracker& tracker,
                                          const std::vector<double>& triggerMomenta,
                                          const std::vector<double>& thresholdProbabilities,
                                          int etaSteps = 50);
    void createTriggerDistanceTuningPlots(Tracker& tracker, const std::vector<double>& triggerMomenta);
    void analyzeGeometry(Tracker& tracker, int nTracks = 1000); // TODO: why virtual?
    void computeBandwidth(Tracker& tracker);
    void computeTriggerFrequency(Tracker& tracker);
    void computeIrradiatedPowerConsumption(Tracker& tracker);
    void analyzePower(Tracker& tracker);
    void createGeometryLite(Tracker& tracker);
    TH2D& getMapPhiEta() { return mapPhiEta; }
    TCanvas& getEtaProfileCanvas() {return etaProfileCanvas;};
    TH1D& getHitDistribution() {return hitDistribution;};
    TProfile& getTotalEtaProfile() {return totalEtaProfile;};
    TGraph& getPowerDensity() {return powerDensity;};
    std::vector<TProfile>& getTypeEtaProfiles() {return typeEtaProfile;};
    std::vector<TObject> getSavingVector();
    TCanvas* getGeomLite() {if (geomLiteCreated) return geomLite; else return NULL; };
    TCanvas* getGeomLiteXY() {if (geomLiteXYCreated) return geomLiteXY; else return NULL; };
    TCanvas* getGeomLiteYZ() {if (geomLiteYZCreated) return geomLiteYZ; else return NULL; };
    TCanvas* getGeomLiteEC() {if (geomLiteECCreated) return geomLiteEC; else return NULL; };
    TH1D& getChanHitDistribution() { return chanHitDistribution; };
    TH1D& getBandwidthDistribution() { return bandwidthDistribution; };
    TH1D& getBandwidthDistributionSparsified() { return bandwidthDistributionSparsified; }
    TH1I& getModuleConnectionsDistribution() { return moduleConnectionsDistribution; }
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
    std::map<std::string, SummaryTable>& getTriggerFrequencyTrueSummaries() { return triggerFrequencyTrueSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerFrequencyInterestingSummaries() { return triggerFrequencyInterestingSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerFrequencyFakeSummaries() { return triggerFrequencyFakeSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerRateSummaries() { return triggerRateSummaries_; }
    std::map<std::string, SummaryTable>& getTriggerEfficiencySummaries() { return triggerEfficiencySummaries_; }
    std::map<std::string, SummaryTable>& getTriggerPuritySummaries() { return triggerPuritySummaries_; }
    std::map<std::string, SummaryTable>& getTriggerDataBandwidthSummaries() { return triggerDataBandwidthSummaries_; }
    std::map<std::string, SummaryTable>& getIrradiatedPowerConsumptionSummaries() { return irradiatedPowerConsumptionSummaries_; }

    SummaryTable& getProcessorConnectionSummary() { return processorConnectionSummary_; }
    std::map<std::string, SummaryTable>& getModuleConnectionSummaries() { return moduleConnectionSummaries_; }
    SummaryTable& getProcessorInboundBandwidthSummary() { return processorInboundBandwidthSummary_; }
    SummaryTable& getProcessorInboundStubPerEventSummary() { return processorInboundStubPerEventSummary_; }
   
    

    // double getEtaMaxMaterial() { return etaMaxMaterial; } 
    double getEtaMaxMaterial() { return getEtaMaxTracking(); }
    double getEtaMaxGeometry() { return etaMaxGeometry; } 
    double getEtaMaxTracking();
    double getEtaMaxTrigger();
    // void setEtaMaxMaterial(const double& newValue) { etaMaxMaterial = newValue; } 
    void setEtaMaxGeometry(const double& newValue) { etaMaxGeometry = newValue; } 
    void addCut(const std::string& cutName, const double& trackingCut, const double& triggerCut);
    const std::vector<double>& getTrackingCuts() { return trackingCuts; }
    const std::vector<double>& getTriggerCuts() { return triggerCuts; }
    const std::vector<std::string>& getCutNames() { return cutNames; }

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
    TH1D rextraservices, rextrasupports;
    TH1D iextraservices, iextrasupports;
    TH1D rglobal, iglobal;

    std::map<std::string, TH1D*> rComponents, iComponents;

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

    std::map<std::string, std::map<std::pair<int, int>, double> > triggerDataBandwidths_;
    std::map<std::string, std::map<std::pair<int, int>, double> > triggerFrequenciesPerEvent_;
    std::map<std::string, SummaryTable> triggerFrequencyTrueSummaries_, triggerFrequencyFakeSummaries_, triggerFrequencyInterestingSummaries_;
    std::map<std::string, SummaryTable> triggerRateSummaries_, triggerEfficiencySummaries_, triggerPuritySummaries_;
    std::map<std::string, SummaryTable> triggerDataBandwidthSummaries_;
    std::map<std::string, SummaryTable> irradiatedPowerConsumptionSummaries_;

    SummaryTable processorConnectionSummary_;
    std::map<std::string, SummaryTable> moduleConnectionSummaries_;
    SummaryTable processorInboundBandwidthSummary_;
    SummaryTable processorInboundStubPerEventSummary_;

    TH1D hitDistribution;
    graphBag myGraphBag;
    mapBag myMapBag;
    profileBag myProfileBag;
    std::map<int, TGraphErrors> spacingTuningGraphs; // TODO: find a way to communicate the limits, not their plots!
    std::map<int, TGraphErrors> spacingTuningGraphsBad; // TODO: find a way to communicate the limits, not their plots!
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
    TProfile totalEtaProfile;
    std::vector<TProfile> typeEtaProfile;

    std::vector<TObject> savingGeometryV; // Vector of ROOT objects to be saved
    std::vector<TObject> savingMaterialV; // Vector of ROOT objects to be saved

    Material findAllHits(MaterialBudget& mb, MaterialBudget* pm, 
                         double& eta, double& theta, double& phi, Track& track);


    void computeDetailedWeights(std::vector<std::vector<ModuleCap> >& tracker, std::map<std::string, SummaryTable>& weightTables, bool byMaterial);
    virtual Material analyzeModules(std::vector<std::vector<ModuleCap> >& tr, double eta, double theta, double phi, Track& t, 
                                    std::map<std::string, Material>& sumComponentsRI, bool isPixel = false);

    int findHitsModules(Tracker& tracker, double z0, double eta, double theta, double phi, Track& t);

    virtual Material findHitsModules(std::vector<std::vector<ModuleCap> >& tr,
                                     double eta, double theta, double phi, Track& t, bool isPixel = false);
    virtual Material findHitsModuleLayer(std::vector<ModuleCap>& layer, double eta, double theta, double phi, Track& t, bool isPixel = false);

    virtual Material findModuleLayerRI(std::vector<ModuleCap>& layer, double eta, double theta, double phi, Track& t, 
                                       std::map<std::string, Material>& sumComponentsRI, bool isPixel = false);
    virtual Material analyzeInactiveSurfaces(std::vector<InactiveElement>& elements, double eta, double theta, 
                                             Track& t, MaterialProperties::Category cat = MaterialProperties::no_cat, bool isPixel = false);
    virtual Material findHitsInactiveSurfaces(std::vector<InactiveElement>& elements, double eta, double theta,
                                              Track& t, bool isPixel = false);

    void calculateGraphs(const std::vector<double>& p,
                         const std::vector<Track>& trackVector,
                         int graphAttributes);
    void fillTriggerEfficiencyGraphs(const std::vector<double>& triggerMomenta,
                                     const std::vector<Track>& trackVector);
    void fillTriggerPerformanceMaps(Tracker& tracker);
    void fillPowerMap(Tracker& tracker);
    void clearMaterialBudgetHistograms();
    void prepareTriggerPerformanceHistograms(const int& nTracks, const double& etaMax, const vector<double>& triggerMomenta, const vector<double>& thresholdProbabilities);
    void preparePowerHistograms();
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
  private:
    // A random number generator
    TRandom3 myDice; 
    int findCellIndexR(double r);
    int findCellIndexEta(double eta);
    int createResetCounters(Tracker& tracker, std::map <std::string, int> &modTypes);
    std::pair <XYZVector, double > shootDirection(double minEta, double maxEta);
    ModuleVector trackHit(const XYZVector& origin, const XYZVector& direction, ModuleVector* properModules);
    void resetTypeCounter(std::map<std::string, int> &modTypes);
    double diffclock(clock_t clock1, clock_t clock2);
    Color_t colorPicker(std::string);
    std::map<std::string, Color_t> colorPickMap;
    Color_t lastPickedColor;
    int geometryTracksUsed;
    int materialTracksUsed;
    void prepareTrackerMap(TH2D& myMap, const std::string& name, const std::string& title);
    void prepareRadialTrackerMap(TH2D& myMap, const std::string& name, const std::string& title);
    void fillAvailableSpacing(Tracker& tracker, std::vector<double>& spacingOptions);
    static const double maximum_n_planes = 13;

    bool isModuleInEtaSector(const Tracker& tracker, const Module* module, int etaSector) const;
    bool isModuleInPhiSector(const Tracker& tracker, const Module* module, int phiSector) const;

    /*
     * Eta values to show results
     */
    //double etaMaxMaterial;
    double etaMaxGeometry;

    std::vector<std::string> cutNames;
    std::vector<double> trackingCuts;
    std::vector<double> triggerCuts;

    static int bsCounter;
  };
}
#endif  /* _ANALYZER_H */


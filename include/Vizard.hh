//
// File:   Vizard.h
// Author: ndemaio
//
// Created on November 20, 2008, 12:41 PM
//

/**
 * @file Vizard.h
 * @brief This is the header file for the core visualisation class
 */

#ifndef _VIZARD_H
#define	_VIZARD_H
// Standard includes
#include <fstream>
#include <sstream>
// ROOT objects
#include <TGeoManager.h>
#include <TGeoMedium.h>
#include <TGeoMaterial.h>
#include <TGeoVolume.h>
#include <TGeoArb8.h>
#include <TGeoMatrix.h>
#include <TFile.h>
#include <THStack.h>
#include <TStyle.h>
#include <TText.h>
#include <TLatex.h>
#include <TColor.h>
#include <TLine.h>
#include <TArrow.h>
#include <TEllipse.h>
#include <TView.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaletteAxis.h>
#include <TProfile.h>
#include <TPolyLine3D.h>
#include <TArc.h>
// Program constants
#include <global_constants.hh>
// Custom objects

#include "MaterialTab.hh"

#include <Tracker.hh>
#include "OuterCabling/OuterCablingMap.hh"
#include "InnerCabling/InnerCablingMap.hh"
#include <Analyzer.hh>
#include <TagMaker.hh>




#include <InactiveSurfaces.hh>
#include "Module.hh"
#include <RootWeb.hh>
#include <vector>
#include <set>
#include <Palette.hh>

#include <PlotDrawer.hh>
#include <AnalyzerVisitors/GeometricInfo.hh>
#include "VizardTools.hh"


using namespace material;


namespace material {
  class WeightDistributionGrid;
}

using namespace boost::accumulators;
using material::WeightDistributionGrid;
//namespace insur { class CablingMap; }
//using insur::CablingMap;

namespace insur {

  class graphIndex {
  public:
    graphIndex() {} ;
    graphIndex(const graphIndex& ref) {
      ideal = ref.ideal;
      p = ref.p;
      name = ref.name;
    };
    std::string name;
    bool ideal;
    double p;
    bool operator<( const graphIndex& other ) const { // returns true if this > other
      if (this->name!=other.name) return ((this->name)>other.name);
      if ((!(this->ideal))&&(other.ideal)) return false;
      if ((this->ideal)&&(!other.ideal)) return true;
      return ((this->p)>other.p);
    }
  };


  typedef std::map<std::string, double> WeightsPerComponent;
  typedef std::map<std::string, WeightsPerComponent> WeightsPerMechanicalCategory;
  typedef std::map<std::string, WeightsPerMechanicalCategory> WeightsPerSubdetector;


  /**
   * @class Vizard
   * @brief This class bundles a number of output functions for different parts and stages of the material budget buildup.
   *
   * It provides one to write a simplified geometry of active and inactive surfaces to a <i>ROOT</i> file, another to
   * save the neighbour relations between different inactive surfaces and layers/discs to a text file, and a third to print
   * the radiation and interaction length histograms to an image and embed that in HTML after a tracker layout has
   * been analysed.
   *
   * A function to write the neighbour relations to a DOT file instead of the quick and dirty internal format that is
   * used now is planned but not implemented yet.
   */
  class Vizard {
  public:
    Vizard();
    virtual ~Vizard();
    bool needsReset() { return geometry_created; }
    void buildVisualization(Tracker& am, InactiveSurfaces& is, bool simplified);
    void display(std::string rootfilename = "");
    void display(Tracker& am, InactiveSurfaces& is, std::string rootfilename = "", bool simplified = true);
    void writeNeighbourGraph(InactiveSurfaces& is);
    void writeNeighbourGraph(InactiveSurfaces& is, std::string outfile);
    void writeNeighbourGraph(InactiveSurfaces& is, std::ostream& outstream);
    void dotGraph(InactiveSurfaces& is, std::string outfile); // temporary, does nothing yet

    // TODO: all these functions should check if the corresponding data is present
    // and return true or false, depending if they created the output or not
    void histogramSummary(Analyzer& a, MaterialBudget& materialBudget, RootWSite& site);
    void histogramSummary(Analyzer& a, MaterialBudget& materialBudget, RootWSite& site, std::string alternativeName);
    void totalMaterialSummary(Analyzer& analyzer, Analyzer& pixelAnalyzer, RootWSite& site);
    void weigthSummary(Analyzer& a, MaterialBudget& materialBudget, WeightDistributionGrid& weightGrid, RootWSite& site, std::string alternativeName);
    bool geometrySummary(Analyzer& a, Tracker& tracker, InactiveSurfaces* inactive, RootWSite& site, bool& debugResolution, std::string alternativeName = "");
    bool outerCablingSummary(Analyzer& a, Tracker& tracker, RootWSite& site);
    bool innerCablingSummary(Analyzer& a, Tracker& tracker, RootWSite& site);
    bool bandwidthSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site);
    bool triggerProcessorsSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site);
    bool errorSummary(Analyzer& a, RootWSite& site, std::string additionalTag, bool isTrigger);
    bool taggedErrorSummary(Analyzer& a, RootWSite& site);
    bool patternRecoSummary(Analyzer& a, mainConfigHandler& mainConfig, RootWSite& site);
    bool triggerSummary(Analyzer& a, Tracker& tracker, RootWSite& site, bool extended);
    bool neighbourGraphSummary(InactiveSurfaces& is, RootWSite& site);
    WeightsPerSubdetector computeDetailedWeights(MaterialBudget& mb, RootWPage& page);
    bool additionalInfoSite(const std::string& settingsfile,
                            Analyzer& analyzer, Analyzer& pixelAnalyzer, Tracker& tracker, RootWSite& site);
    bool makeLogPage(RootWSite& site);
    void setCommandLine(std::string commandLine) { commandLine_ = commandLine; };
    void createXmlSite(RootWSite& site,std::string xmldir,std::string layoutdir);

  protected:
    TGeoManager* gm;
    TGeoVolume* top;
    TGeoVolumeAssembly* active;
    TGeoVolumeAssembly* inactive;
    TGeoVolumeAssembly* barrels;
    TGeoVolumeAssembly* endcaps;
    TGeoVolumeAssembly* services;
    TGeoVolumeAssembly* supports;
    TGeoMedium* medvac;
    TGeoMedium* medact;
    TGeoMedium* medserf;
    TGeoMedium* medlazy;
    TGeoMaterial* matvac;
    TGeoMaterial* matact;
    TGeoMaterial* matserf;
    TGeoMaterial* matlazy;
  private:
    TProfile* totalEtaProfileSensors_ = 0, *totalEtaProfileSensorsPixel_ = 0;
    TProfile* totalEtaProfileLayers_ = 0, *totalEtaProfileLayersPixel_ = 0;
    std::map<std::string, TH1D*> rCompsPixelTrackingVolume_, iCompsPixelTrackingVolume_;
    std::map<std::string, TH1D*> rServicesDetailsPixelTrackingVolume_, iServicesDetailsPixelTrackingVolume_;
  bool geometry_created;
    std::string commandLine_;
    int detailedModules(std::vector<Layer*>* layers,
                        TGeoVolume* v, TGeoCombiTrans* t, TGeoVolumeAssembly* a, int counter);
    TGeoCombiTrans* modulePlacement(Module* m, TGeoVolume* v);
    double averageHistogramValues(TH1D& histo, double cutoff);
    double averageHistogramValues(TH1D& histo, double cutoffStart, double cutoffEnd);

    void createSummaryCanvas(double maxZ, double maxRho, Analyzer& analyzer, std::unique_ptr<TCanvas> &YZCanvas, std::unique_ptr<TCanvas> &XYCanvas, std::unique_ptr<TCanvas> &XYCanvasEC);
    void createSummaryCanvasNicer(Tracker& tracker, std::unique_ptr<TCanvas> &YZCanvas, std::unique_ptr<TCanvas> &YZCanvasBarrel, std::unique_ptr<TCanvas> &XYCanvas, std::vector<std::unique_ptr<TCanvas> > &XYCanvasEC);

    // OT CABLING
    void createOuterCablingPlotsBundles(const Tracker& tracker, std::unique_ptr<TCanvas> &YZCanvas, std::unique_ptr<TCanvas> &XYCanvas, std::unique_ptr<TCanvas> &XYNegCanvas, 
					       std::vector<std::unique_ptr<TCanvas> > &XYPosBundlesDisks, std::vector<std::unique_ptr<TCanvas> > &XYPosBundlesDiskSurfaces,
					       std::vector<std::unique_ptr<TCanvas> > &XYNegBundlesDisks, std::vector<std::unique_ptr<TCanvas> > &XYNegBundlesDiskSurfaces);
    void createOuterCablingPlotsDTCs(Tracker& tracker, std::unique_ptr<TCanvas> &YZCanvas, std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYNegFlatCanvas, 
					    std::unique_ptr<TCanvas> &XYCanvas, std::unique_ptr<TCanvas> &XYFlatCanvas, std::vector<std::unique_ptr<TCanvas> > &XYCanvasEC);
    void createOuterCablingPlotsServicesChannelsOptical(Tracker& tracker, const OuterCablingMap* myCablingMap, std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYNegFlatCanvas, std::unique_ptr<TCanvas> &XYCanvas, std::unique_ptr<TCanvas> &XYFlatCanvas, std::vector<std::unique_ptr<TCanvas> > &XYCanvasEC);
    void createOuterCablingPlotsServicesChannelsPower(Tracker& tracker, const OuterCablingMap* myCablingMap,
						     std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYNegFlatCanvas, std::unique_ptr<TCanvas> &XYCanvas, std::unique_ptr<TCanvas> &XYFlatCanvas, 
						     std::vector<std::unique_ptr<TCanvas> > &XYCanvasesDisk, std::vector<std::unique_ptr<TCanvas> > &XYNegCanvasesDisk);
    RootWTable* opticalServicesChannels(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide, const ChannelSlot requestedSlot = ChannelSlot::UNKNOWN);
    void analyzeOpticalServicesChannels(const OuterCablingMap* myCablingMap, std::map<int, std::vector<int> > &cablesPerChannel, std::map<int, int> &psBundlesPerChannel, std::map<int, int> &ssBundlesPerChannel, const bool isPositiveCablingSide, const ChannelSlot requestedSlot = ChannelSlot::UNKNOWN);
    RootWTable* createOpticalServicesChannelTable(const std::map<int, std::vector<int> > &cablesPerChannel, const std::map<int, int> &psBundlesPerChannel, const std::map<int, int> &ssBundlesPerChannel, const bool isPositiveCablingSide, const ChannelSlot requestedSlot = ChannelSlot::UNKNOWN);
    RootWTable* powerServicesChannels(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide, const std::vector<ChannelSlot>& slots);
    void analyzePowerServicesChannels(const OuterCablingMap* myCablingMap, std::map<int, int> &psBundlesPerChannel, std::map<int, int> &ssBundlesPerChannel, const bool isPositiveCablingSide, const ChannelSlot requestedSlot = ChannelSlot::UNKNOWN);
    void createPowerServicesChannelTable(RootWTable* channelsTable, const std::map<int, int> &psBundlesPerChannel, const std::map<int, int> &ssBundlesPerChannel, const bool isPositiveCablingSide, const ChannelSlot requestedSlot = ChannelSlot::UNKNOWN);

    // IT CABLING
    void createInnerCablingPlotsPowerChains(const Tracker& tracker, 
						   std::vector<std::unique_ptr<TCanvas> > &ZPhiLayerPlots,
						   std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYCentralCanvas, std::unique_ptr<TCanvas> &XYCanvas,
						   std::vector<std::unique_ptr<TCanvas> > &XYPosPowerChainsDiskSurfaces);
    void createInnerCablingPlotsGBTs(const Tracker& tracker,
						 std::vector<std::unique_ptr<TCanvas> > &ZPhiLayerPlots,
						 std::vector<std::unique_ptr<TCanvas> > &XYPosGBTsDiskSurfaces);
    void createInnerCablingPlotsBundles(const Tracker& tracker,
						    std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYPosCanvas,
						    std::vector<std::unique_ptr<TCanvas> > &XYPosBundlesDisks);
    void createInnerCablingPlotsDTCs(const Tracker& tracker,
						 std::unique_ptr<TCanvas> &RZCanvas,
						 std::unique_ptr<TCanvas> &XYPosCanvas,
						 std::vector<std::unique_ptr<TCanvas> > &XYPosDTCsDisks);
    void computeInnerCablingCount(const InnerCablingMap* myInnerCablingMap,
				  int& numSensorsOneXSide, int& numSensorsPlusXSidePlusZEnd, int& numSensorsPlusXSideMinusZEnd,
				  int& numPowerChainsOneXSide, int& numPowerChainsPlusXSidePlusZEnd, int& numPowerChainsPlusXSideMinusZEnd,
				  int& numELinksOneXSide, int& numELinksPlusXSidePlusZEnd, int& numELinksPlusXSideMinusZEnd,
				  int& numBundlesOneXSide, int& numBundlesPlusXSidePlusZEnd, int& numBundlesPlusXSideMinusZEnd,
				  int& numGBTsOneXSide, int& numGBTsPlusXSidePlusZEnd, int& numGBTsPlusXSideMinusZEnd,
				  int& numDTCsOneXSide, int& numDTCsPlusXSidePlusZEnd, int& numDTCsPlusXSideMinusZEnd) const;

    enum {ViewSectionXY=3, ViewSectionYZ=1, ViewSectionXZ=2};
    void drawEtaTicks(double maxL, double maxR, double tickDistance, double tickLength, double textDistance, Style_t labelFont, Float_t labelSize,
                      double etaStep, double etaMax, double etaLongLine);
    void drawTicks(Analyzer& a, TView* myView, double maxL, double maxR, int noAxis=1, double spacing = 100., Option_t* option = "same"); // shold become obsolete
    void drawGrid(double maxL, double maxR, int noAxis=1, double spacing = 100., Option_t* option = "same"); // shold become obsolete
    bool drawEtaProfilesAny(TProfile& totalEtaProfile, std::vector<TProfile>& etaProfiles, bool total=true); // generic business logic called by hit or stub version with appropriate parameters
    bool drawEtaProfiles(TCanvas& myCanvas, Analyzer& analyzer); // for hits
    bool drawEtaProfiles(TVirtualPad& myPad, Analyzer& analyzer); // for hits
    bool drawEtaProfilesSensors(TCanvas& myCanvas, Analyzer& analyzer, bool total=true);
    bool drawEtaProfilesSensors(TVirtualPad& myPad, Analyzer& analyzer, bool total=true);
    bool drawEtaProfilesStubs(TCanvas& myCanvas, Analyzer& analyzer);
    bool drawEtaProfilesStubs(TVirtualPad& myPad, Analyzer& analyzer);
    bool drawTracksDistributionPerNumberOfStubs(TCanvas& myCanvas, Analyzer& analyzer, const bool isPixelTracker); // plots for each number of stubs
    bool drawTracksDistributionPerNumberOfStubs(TVirtualPad& myPad, Analyzer& analyzer, const bool isPixelTracker); // plots for each number of stubs
    bool drawEtaProfilesLayers(TCanvas& myCanvas, Analyzer& analyzer);
    bool drawEtaProfilesLayers(TVirtualPad& myPad, Analyzer& analyzer);

    // COVERAGE PER LAYER: HITS AND STUBS
    bool drawHitCoveragePerLayer(RootWPage& myPage, Analyzer& analyzer, const bool isPixelTracker);
    bool drawStubCoveragePerLayer(RootWPage& myPage, Analyzer& analyzer, const bool isPixelTracker);
    bool drawStubWith3HitsCoveragePerLayer(RootWPage& myPage, Analyzer& analyzer);
    bool drawCoveragePerlayer(RootWPage& myPage, const bool isPixelTracker, const std::string type, 
			      std::map<std::string, TProfile>& coveragePerLayer,			      
			      std::map<std::string, CoveragePerNumberOfHits>& coveragePerLayerDetails,
			      std::map<std::string, std::map<int, double> >& stubWith3HitsCountPerDiskAndRing);
    void drawCoveragePerlayerDetails(CoveragePerNumberOfHits& coveragePerLayerDetails, const std::string type, TLegend* layerLegend, const int plotMaxNumberOfHits);

    int momentumColor(int iMomentum);
    void closeGraph(TGraph& myGraph);

    double getDrawAreaZ(const Tracker& tracker) const { return tracker.maxZ()*1.1; }
    double getDrawAreaR(const Tracker& tracker) const { return tracker.maxR()*1.1; }
    double getDrawAreaX(const Tracker& tracker) const { return tracker.maxR()*1.1; }
    double getDrawAreaY(const Tracker& tracker) const { return tracker.maxR()*1.1; }

    void fillPlotMap(std::string& plotName,
                     std::map<graphIndex, TGraph*>& myPlotMap,
                     Analyzer *a,
                     std::map<int, TGraph>& (Analyzer::*retriveFunction)(bool, bool), bool isTrigger);

    void fillTaggedPlotMap(GraphBag& gb,
                           const string& plotName,
                           int graphType,
                           const string& tag,
                           std::map<graphIndex, TGraph*>& myPlotMap);
    std::string triggerSectorMapCsv_;
    std::string moduleConnectionsCsv_;

    void createTriggerSectorMapCsv(const TriggerSectorMap& tsm);
    void createModuleConnectionsCsv(const ModuleConnectionMap& moduleConnections);
    std::string createAllModulesCsv(const Tracker& t, bool& withHeader);
    std::string createBarrelModulesCsv(const Tracker& t);
    std::string createEndcapModulesCsv(const Tracker& t);
    std::string createModulesDetIdListCsv();
    std::string createSensorsDetIdListCsv();

    std::string createChemicalElementsCsv();
    std::string createChemicalMixturesCsv(const bool hasChemicalFormula);

    std::string createModulesToDTCsCsv(const Tracker& t, const bool isPositiveCablingSide);
    std::string createDTCsToModulesCsv(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide);
    std::string createCMSSWOuterTrackerCablingMapCsv(const Tracker& tracker);
    std::string createBundlesToEndcapModulesCsv(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide);
    std::string countBundlesToEndcapModulesCombinations(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide);

    std::string createInnerTrackerModulesToDTCsCsv(const Tracker& tracker);
    std::string createInnerTrackerDTCsToModulesCsv(const InnerCablingMap* myInnerCablingMap) ;

    TProfile* newProfile(TH1D* sourceHistogram, double xlow, double xup, int desiredNBins = 0);
    TProfile& newProfile(const TGraph& sourceGraph, double xlow, double xup, int nrebin = 1, int nBins = 0);
    TProfile& newProfile_timesSin(const TGraph& sourceGraph, double xlow, double xup, int nrebin = 1, int nBins = 0);
    void stackHistos(std::vector<std::pair<std::string, TH1D*>>& histoMap, RootWTable*& myTable, int& index, THStack*& totalStack, THStack*& myStack, TLegend*& legend, bool& isRadiation);
    void stackHistos(std::map<std::string, TH1D*>& histoMap, RootWTable*& myTable, int& index, THStack*& totalStack, THStack*& myStack, TLegend*& legend, bool& isRadiation);

    // int getNiceColor(unsigned int plotIndex);
    std::vector<Tracker*> trackers_;
    std::vector<MaterialBudget*> materialBudgets_;
    std::unique_ptr<TCanvas> drawFullLayoutRZ();
    std::unique_ptr<TCanvas> drawFullLayoutServicesRZ();
    std::unique_ptr<TCanvas> drawFullLayoutBarrelXY();

    void plotAndPrintVolumeMaterials(WeightsPerSubdetector& totalWeights, std::stringstream& allVolumesStream, std::stringstream& modulesStream, 
				     const std::map<LocalElement, double, ElementNameCompare>& allMasses, 
				     const double z1, const double z2, const double r1, const double r2, const double rl, const double il,
				     std::map<std::string, int>& subdetectorColors, const std::vector<int>& allColors, int& colorIndex,
				     const bool isModule, const int serviceId = 0, const double serviceLength = 0., 
				     const Module* detectorModule = nullptr, const bool printModulesCsv = false);
    void plotVolumeBox(const std::string subdetectorName, 
		       std::map<std::string, int>& subdetectorColors, const std::vector<int>& allColors, int& colorIndex,
		       const bool isEmpty, 
		       const double z1, const double z2, const double r1, const double r2, const bool isFilled = true);
    const int computeSubdetectorColor(const std::string subdetectorName, 
				      std::map<std::string, int>& subdetectorColors, const std::vector<int>& allColors, int& colorIndex,
				      const bool isEmpty);

    void drawCircle(double radius, bool full, int color=kBlack);
    void drawPhiSectorsBoundaries(const double phiSectorWidth, const bool isRotatedY180 = false);
    void drawFrameOfReference(const bool isRotatedY180, const double scalingFactor = 1.);
    const std::pair<double, double> computeInnerCablingPlotsScalingFactors(const Tracker& tracker);
    const std::pair<double, double> computeInnerCablingPlotsMaxRadii(const Tracker& tracker);
    void computeServicesChannelsLegend(TLegend* legend, const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide, const bool isPowerCabling);
  };



}
#endif	/* _VIZARD_H */

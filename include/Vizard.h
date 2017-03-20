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
#include <TView.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaletteAxis.h>
#include <TProfile.h>
#include <TPolyLine3D.h>
#include <TArc.h>
// Program constants
#include <global_constants.h>
// Custom objects
#include <SimParms.h>
#include <Tracker.h>
#include <Analyzer.h>
#include <TagMaker.h>

#include <InactiveSurfaces.h>
#include "Module.h"
#include <rootweb.hh>
#include <vector>
#include <set>
#include <Palette.h>

#include <PlotDrawer.h>
#include "VizardTools.hh"

namespace material {
  class WeightDistributionGrid;
}

using namespace boost::accumulators;
using material::WeightDistributionGrid;

namespace insur {

  /*
   * Assorted messages that may pop up
   */
  static const std::string msg_uninitialised = "Vizard::buildVisualization(am, is) needs to be called first to build the visual geometry objects.";
  static const std::string root_wrong = "Something went wrong creating output file. Existing geometry was not written to file.";
  static const std::string graph_wrong = "File stream reported error state: neighbour graph not written to file.";
  static const std::string exc_badalloc_graph = "Error: caught bad_alloc exception in Vizard::writeNeighbourGraph(). ";
  static const std::string graph_nowrite = "Neighbour graph was not written to file.";

  //clearStart="<tt>";
  //clearEnd="</tt>";

  // gStyle stuff
  static const int style_grid = 3;
  static const int materialNBins = 300;

  // Colors for plot background and such
  static const int color_plot_background       = kWhite;
  static const int color_pad_background        = kGray;
  static const int color_grid                  = kGreen-10;
  static const int color_hard_grid             = kGray;
  static const std::vector<string> color_names = {"Black","BrightBlue","Red","BrightGreen","Yellow","Pink","Aqua","Green","Blue"};


  // Pads to plot the tracker ortho views
  static const unsigned int padYZ = 1;
  static const unsigned int padXY = 2;
  static const unsigned int padProfile = 3;
  static const unsigned int padEC = 4;

  // Formatting parameters
  static const int coordPrecision = 3;
  static const int zOverlapPrecision = 5;
  static const int anglePrecision = 1;
  static const int areaPrecision = 1;
  static const int occupancyPrecision = 1;
  static const int rphiResolutionPrecision = 0;
  static const int rphiResolutionRmsePrecision = 1;
  static const int pitchPrecision = 0;
  static const int stripLengthPrecision = 1;
  static const int millionChannelPrecision = 2;
  static const int totalPowerPrecision = 2;
  static const int modulePowerPrecision = 0;
  static const int powerPrecision = 1;
  static const int costPrecision  = 1;
  static const int powerPerUnitPrecision = 2;
  static const int costPerUnitPrecision  = 1;
  static const int minimumBiasPrecision = 0;
  static const int luminosityPrecision = 0;
  static const int alphaParamPrecision = 4;
  static const int weightPrecision = 0;

  static const int tableResolutionPrecisionHigh = 5;
  static const int tableResolutionPrecisionStd  = 3;

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
    void histogramSummary(Analyzer& a, MaterialBudget& materialBudget, bool debugServices, RootWSite& site);
    void histogramSummary(Analyzer& a, MaterialBudget& materialBudget, bool debugServices, RootWSite& site, std::string alternativeName);
    void weigthSummart(Analyzer& a, WeightDistributionGrid& weightGrid, RootWSite& site, std::string alternativeName);
    bool geometrySummary(Analyzer& a, Tracker& tracker, SimParms& simparms, InactiveSurfaces* inactive, RootWSite& site, bool& debugResolution, std::string alternativeName = "");
    bool bandwidthSummary(Analyzer& analyzer, Tracker& tracker, SimParms& simparms, RootWSite& site);
    bool triggerProcessorsSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site);
    bool errorSummary(Analyzer& a, RootWSite& site, std::string additionalTag, bool isTrigger);
    bool taggedErrorSummary(Analyzer& a, RootWSite& site);
    bool triggerSummary(Analyzer& a, Tracker& tracker, RootWSite& site, bool extended);
    bool neighbourGraphSummary(InactiveSurfaces& is, RootWSite& site);
    void drawInactiveSurfacesSummary(MaterialBudget& mb, RootWPage& page);
    bool additionalInfoSite(const std::string& settingsfile,
                            Analyzer& analyzer, Analyzer& pixelAnalyzer, Tracker& tracker, SimParms& simparms, RootWSite& site);
    bool makeLogPage(RootWSite& site);
    void setCommandLine(std::string commandLine) { commandLine_ = commandLine; }
      bool createXmlSite(RootWSite& site,std::string xmldir,std::string layoutdir);

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
  bool geometry_created;
    std::string commandLine_;
    int detailedModules(std::vector<Layer*>* layers,
                        TGeoVolume* v, TGeoCombiTrans* t, TGeoVolumeAssembly* a, int counter);
    TGeoCombiTrans* modulePlacement(Module* m, TGeoVolume* v);
    double averageHistogramValues(TH1D& histo, double cutoff);
    double averageHistogramValues(TH1D& histo, double cutoffStart, double cutoffEnd);

    void createSummaryCanvas(double maxZ, double maxRho, Analyzer& analyzer, TCanvas *&YZCanvas, TCanvas *&XYCanvas, TCanvas *&XYCanvasEC);
    void createSummaryCanvasNicer(Tracker& tracker, TCanvas *&YZCanvas, TCanvas *&YZCanvasBarrel, TCanvas *&XYCanvas, std::vector<TCanvas*> &XYCanvasEC);

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
    bool drawEtaProfilesLayers(TCanvas& myCanvas, Analyzer& analyzer);
    bool drawEtaProfilesLayers(TVirtualPad& myPad, Analyzer& analyzer);
    bool drawEtaCoverageAny(RootWPage& myPage, std::map<std::string, TProfile>& layerEtaCoverage, const std::string& type); // generic business logic called by hit or stub version
    bool drawEtaCoverage(RootWPage& myPage, Analyzer& analyzer); // for hits
    bool drawEtaCoverageStubs(RootWPage& myPage, Analyzer& analyzer);
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

    TProfile* newProfile(TH1D* sourceHistogram, double xlow, double xup, int desiredNBins = 0);
    TProfile& newProfile(const TGraph& sourceGraph, double xlow, double xup, int nrebin = 1, int nBins = 0);
    TProfile& newProfile_timesSin(const TGraph& sourceGraph, double xlow, double xup, int nrebin = 1, int nBins = 0);
    void stackHistos(std::vector<std::pair<std::string, TH1D*>>& histoMap, RootWTable*& myTable, int& index, THStack*& totalStack, THStack*& myStack, TLegend*& legend, bool& isRadiation);
    void stackHistos(std::map<std::string, TH1D*>& histoMap, RootWTable*& myTable, int& index, THStack*& totalStack, THStack*& myStack, TLegend*& legend, bool& isRadiation);
    // int getNiceColor(unsigned int plotIndex);
    std::vector<Tracker*> trackers_;
    TCanvas* drawFullLayoutRZ();
    TCanvas* drawFullLayoutBarrelXY();

    void drawCircle(double radius, bool full, int color=kBlack);

  };



}
#endif	/* _VIZARD_H */

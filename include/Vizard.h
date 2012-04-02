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
// Program constants
#include <global_constants.h>
// Custom objects
#include <tracker.hh>
#include <Analyzer.h>
#include <InactiveSurfaces.h>
#include <rootweb.hh>
#include <vector>
#include <set>
#include <Palette.h>

namespace insur {

    /*
     * Assorted messages that may pop up
     */
    static const std::string msg_uninitialised = "Vizard::buildVisualization(am, is) needs to be called first to build the visual geometry objects.";
    static const std::string root_wrong = "Something went wrong creating output file. Existing geometry was not written to file.";
    static const std::string graph_wrong = "File stream reported error state: neighbour graph not written to file.";
    static const std::string exc_badalloc_graph = "Error: caught bad_alloc exception in Vizard::writeNeighbourGraph(). ";
    static const std::string graph_nowrite = "Neighbour graph was not written to file.";

    // Some strings for the html formatting
    static const std::string subStart = "<sub>";      // These only should be needed
    static const std::string subEnd = "</sub>";
    static const std::string superStart = "<sup>";
    static const std::string superEnd = "</sup>";
    static const std::string smallStart = "<small>";
    static const std::string smallEnd = "</small>";
    static const std::string emphStart="<b>";
    static const std::string emphEnd="</b>";
    //clearStart="<tt>";
    //clearEnd="</tt>";

    // gStyle stuff
    static const int style_grid = 3;

    // Colors for plot background and such
    static const int color_plot_background = kWhite;
    static const int color_pad_background = kGray;
    static const int color_grid = kGreen-10;
    static const int color_hard_grid = kGray;

    // Pads to plot the tracker ortho views
    static const unsigned int padYZ = 1;
    static const unsigned int padXY = 2;
    static const unsigned int padProfile = 3;
    static const unsigned int padEC = 4;

    // Formatting parameters
    static const int coordPrecision = 0;
    static const int areaPrecision = 1;
    static const int occupancyPrecision = 1;
    static const int rphiResolutionPrecision = 0;
    static const int pitchPrecision = 0;
    static const int stripLengthPrecision = 1;
    static const int millionChannelPrecision = 2;
    static const int powerPrecision = 1;
    static const int costPrecision  = 1;
    static const int powerPerUnitPrecision = 2;
    static const int costPerUnitPrecision  = 1;
    static const int minimumBiasPrecision = 0;
    static const int weightPrecision = 0;

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
        void histogramSummary(Analyzer& a, std::string outfilename);

	// TODO: all these functions should check if the corresponding data is present
	// and return true or false, depending if they created the output or not
	void histogramSummary(Analyzer& a, RootWSite& site);
	void histogramSummary(Analyzer& a, RootWSite& site, std::string alternativeName);
	void weigthSummart(Analyzer& a, RootWSite& site, std::string alternativeName);
	bool geometrySummary(Analyzer& a, Tracker& tracker, RootWSite& site);
	bool bandwidthSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site);
	bool irradiatedPowerSummary(Analyzer& a, RootWSite& site);
	bool errorSummary(Analyzer& a, RootWSite& site, std::string additionalTag, bool isTrigger);
	bool triggerSummary(Analyzer& a, RootWSite& site, bool extended);
	bool neighbourGraphSummary(InactiveSurfaces& is, RootWSite& site); 
	bool additionalInfoSite(const std::string& geomfile, const std::string& settingsfile,
				const std::string& matfile, const std::string& pixmatfile,
				bool defaultMaterial, bool defaultPixelMaterial,
				Analyzer& analyzer, Tracker& tracker, RootWSite& site);
	bool makeLogPage(RootWSite& site);
	std::string getSummaryString();
	std::string getSummaryLabelString();

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
        bool geometry_created;
        int detailedModules(std::vector<Layer*>* layers,
                TGeoVolume* v, TGeoCombiTrans* t, TGeoVolumeAssembly* a, int counter);
        TGeoCombiTrans* modulePlacement(Module* m, TGeoVolume* v);
        double averageHistogramValues(TH1D& histo, double cutoff);
        double averageHistogramValues(TH1D& histo, double cutoffStart, double cutoffEnd);

	void createSummaryCanvas(double maxZ, double maxRho, Analyzer& analyzer, TCanvas *&YZCanvas, TCanvas *&XYCanvas, TCanvas *&XYCanvasEC);
	void createSummaryCanvasNice(Tracker& tracker, TCanvas *&YZCanvas, TCanvas *&XYCanvas, TCanvas *&XYCanvasEC);
	void createColorPlotCanvas(Tracker& tracker, int plotVariable, TCanvas *&RZCanvas);

	enum {ViewSectionXY=3, ViewSectionYZ=1, ViewSectionXZ=2};
	void drawEtaTicks(double maxL, double maxR, double tickDistance, double tickLength, double textDistance, Style_t labelFont, Float_t labelSize,
			  double etaStep, double etaMax, double etaLongLine);
	void drawTicks(TView* myView, double maxL, double maxR, int noAxis=1, double spacing = 100., Option_t* option = "same"); // shold become obsolete
	void drawGrid(double maxL, double maxR, int noAxis=1, double spacing = 100., Option_t* option = "same"); // shold become obsolete
	bool drawEtaProfiles(TCanvas& myCanvas, Analyzer& analyzer);
	bool drawEtaProfiles(TVirtualPad& myPad, Analyzer& analyzer);
        int momentumColor(int iMomentum);
	void closeGraph(TGraph& myGraph);

	void fillPlotMap(std::string& plotName, 
			 std::map<graphIndex, TGraph*>& myPlotMap,
			 Analyzer *a,
			 std::map<double, TGraph>& (Analyzer::*retriveFunction)(bool, bool), bool isTrigger);
        std::string summaryCsv_;
        std::string summaryCsvLabels_;
	std::string occupancyCsv_;
	void setSummaryString(std::string);
	void addSummaryElement(std::string element, bool first = false);
	void setSummaryLabelString(std::string);
	void addSummaryLabelElement(std::string element, bool first = false);
	void addSummaryElement(double element, bool first = false);

	void setOccupancyString(std::string newString) { occupancyCsv_ = newString; };
	void addOccupancyElement(double element);
	void addOccupancyElement(std::string element);
	void addOccupancyEOL();

	// int getNiceColor(unsigned int plotIndex);
    };

    //    /**
    //     * @class moduleGraphicPosition
    //     *
    //     * @brief a class to hold a simplifies (2D) representation of the
    //     * module position in a selected plane with the purpose of drawing
    //     * it on a canvas. The apporximation allows to avoid duplication
    //     * of drawing (which is heavy for vector graphics)
    //     */
    //    class moduleGraphicPosition {
    //    private:df
    //      // Approximate to the tenth of millimiter
    //      static const double fraction=10;
    //      double round(const double& x);
    //    public:
    //      double x[4];
    //      double y[4];
    //      Color_t color;
    //      void setRZ(const Module&);
    //      void setRPhi(const Module&);
    //      bool operator <(const moduleGraphicPosition& b) const;
    //    };

    /**
     * @class linePosition
     *
     * @brief a class to hold a the position and color of a line
     */
    class linePosition {
    private:
      int x[2];
      int y[2];
      static int round(const double& x);
    public:
      // Member variables
      // Approximate to the tenth of millimiter
      static const double fraction=10;
      Color_t color;
      // Member functions
      bool operator <(const linePosition& b) const;
      linePosition& operator=(const linePosition& lp);
      void draw(Width_t lwidth) const;
      void setCoordinates(const double& x1, const double& y1, const double& x2, const double& y2);
      void setPointCoordinates(const int& iPoint, const double& x, const double& y);
      bool isTrivial();
    };

    /**
     * @class moduleValue
     *
     * @brief a class to hold values for the modules to be drawn
     */
    class moduleValue {
    private:
      double nValues;
      double sumValues;
    public:
      moduleValue();
      void fill(const double& value);
      double getAverage() const;
    };

    /**
     * @class linePlacer
     *
     * @brief a class to place lines in the right position and color,
     * depending on the module
     */
    class linePlacer {
    private:
      std::set<linePosition> linePositions;
      std::map<linePosition, moduleValue> lineValues;
      TPaletteAxis *myPalette;
      static const Width_t default_linewidth = 2;
    public:
      linePlacer();
      void setPalette(TPaletteAxis* newPalette) { myPalette = newPalette ; }
      void setRZ(const Module&);
      void setRPhi(const Module&);
      void setRZValue(const Module&, const double& value);
      void setRPhiValue(const Module&, const double& value);
      void drawLines();
      void drawLineValues();
      void drawLines(Width_t lwidth);
      void drawLineValues(Width_t lwidth);
      void setMinMax(TH2D& myFrame);
    };

    

}
#endif	/* _VIZARD_H */


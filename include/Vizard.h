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
#include <TColor.h>
#include <TView.h>
#include <TLegend.h>
#include <TGraph.h>
// Program constants
#include <global_constants.h>
// Custom objects
#include <tracker.hh>
#include <Analyzer.h>
#include <InactiveSurfaces.h>
#include <rootweb.hh>
#include <vector>

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
    static const int pitchPrecision = 0;
    static const int stripLengthPrecision = 1;
    static const int millionChannelPrecision = 2;
    static const int powerPrecision = 1;
    static const int costPrecision  = 1;
    static const int powerPerUnitPrecision = 2;
    static const int costPerUnitPrecision  = 1;
    static const int minimumBiasPrecision = 0;
    
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
        void dotGraph(InactiveSurfaces& is, std::string outfile); // temporary, does nothing yet
        void histogramSummary(Analyzer& a, std::string outfilename);
#ifdef USING_ROOTWEB
	// TODO: all these functions should check if the corresponding data is present
	// and return true or false, depending if they created the output or not
	void histogramSummary(Analyzer& a, RootWSite& site);
	bool geometrySummary(Analyzer& a, Tracker& tracker, RootWSite& site);
	bool bandwidthSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site);
        bool errorSummary(Analyzer& a, RootWSite& site);
	bool additionalInfoSite(std::string& geomfile, std::string& settingsfile, std::string& matfile, Analyzer& analyzer, Tracker& tracker, RootWSite& site);
	bool makeLogPage(RootWSite& site);
#endif
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
	void createSummaryCanvas(double maxZ, double maxRho, Analyzer& analyzer, TCanvas *&summaryCanvas, TCanvas *&YZCanvas);
	enum {ViewSectionXY=3, ViewSectionYZ=1, ViewSectionXZ=2};
	void drawTicks(TView* myView, double maxL, double maxR, int noAxis=1, double spacing = 100., Option_t* option = "same"); // shold become obsolete
	void drawGrid(double maxL, double maxR, int noAxis=1, double spacing = 100., Option_t* option = "same"); // shold become obsolete
    };
}
#endif	/* _VIZARD_H */


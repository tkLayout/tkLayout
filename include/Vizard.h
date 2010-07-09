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

#include <fstream>
#include <sstream>
#include <TGeoManager.h>
#include <TGeoMedium.h>
#include <TGeoMaterial.h>
#include <TGeoVolume.h>
#include <TGeoArb8.h>
#include <TGeoMatrix.h>
#include <TFile.h>
#include <THStack.h>
#include <global_constants.h>
#include <tracker.hh>
#include <Analyzer.h>
#include <InactiveSurfaces.h>
#include <rootweb.hh>

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
        void buildVisualization(Tracker& am, InactiveSurfaces& is, bool simplified);
        void display(std::string rootfilename = "");
        void display(Tracker& am, InactiveSurfaces& is, std::string rootfilename = "", bool simplified = true);
        void writeNeighbourGraph(InactiveSurfaces& is);
        void writeNeighbourGraph(InactiveSurfaces& is, std::string outfile);
        void dotGraph(InactiveSurfaces& is, std::string outfile); // temporary, does nothing yet
        void histogramSummary(Analyzer& a, std::string outfilename);
	void histogramSummary(Analyzer& a, RootWSite& site);
	bool geometrySummary(Analyzer& a, Tracker& tracker, RootWSite& site);
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
    };
}
#endif	/* _VIZARD_H */


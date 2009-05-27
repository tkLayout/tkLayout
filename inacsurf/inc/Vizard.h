// 
// File:   Vizard.h
// Author: ndemaio
//
// Created on November 20, 2008, 12:41 PM
//

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
#include <global_constants.h>
#include <tracker.hh>
#include<Analyzer.h>
#include <InactiveSurfaces.h>
namespace insur {
    // TODO: move messages to static strings
    
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
        TGeoCombiTrans* modulePlacement(Module* m, TGeoVolume* v);
    };
}
#endif	/* _VIZARD_H */


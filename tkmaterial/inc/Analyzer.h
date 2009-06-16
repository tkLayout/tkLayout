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
#define	_ANALYZER_H

#include<string>
#include<vector>
#include<module.hh>
#include<ModuleCap.h>
#include<InactiveElement.h>
#include<InactiveSurfaces.h>
#include <MaterialBudget.h>
namespace insur {
    /**
     * A warning that may occur during processing
     */
    static const std::string msg_module_warning = "Warning: tracker module with undefined subdetector type found.";
    
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
    class Analyzer {
    public:
        Analyzer() {}
        virtual ~Analyzer() {}
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
        TH1D& getHistoGlobalR() { return rglobal; }
        TH1D& getHistoGlobalI() { return iglobal; }
        TH2D& getHistoIsoR() { return isor; }
        TH2D& getHistoIsoI() { return isoi; }
        virtual void analyzeMaterialBudget(MaterialBudget& mb, int etaSteps = 50);
    protected:
        static const int r_bins = 25;
        static const int z_bins = 50;
        TH1D ractivebarrel, ractiveendcap, rserfbarrel, rserfendcap, rlazybarrel, rlazyendcap, rlazytube, rlazyuserdef;
        TH1D iactivebarrel, iactiveendcap, iserfbarrel, iserfendcap, ilazybarrel, ilazyendcap, ilazytube, ilazyuserdef;
        TH1D rbarrelall, rendcapall, ractiveall, rserfall, rlazyall;
        TH1D ibarrelall, iendcapall, iactiveall, iserfall, ilazyall;
        TH1D rglobal, iglobal;
        TH2D isor, isoi;
        virtual std::pair<double, double> analyzeModules(std::vector<std::vector<ModuleCap> >& tr, double eta, double theta);
        virtual std::pair<double, double> findModuleLayerRI(std::vector<ModuleCap>& layer, double theta, double phi);
        virtual std::pair<double, double> analyzeInactiveSurfaces(std::vector<InactiveElement>& elements, double eta,
                                                                                                       double theta, MaterialProperties::Category cat = MaterialProperties::no_cat);
        void clearHistograms();
        void setHistogramBinsBoundaries(int bins, double min, double max);
    private:
    };
}
#endif	/* _ANALYZER_H */


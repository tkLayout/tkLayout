// 
// File:   Extractor.h
// Author: ndemaio
//
// Created on January 31, 2010, 3:33 AM
//

/**
 * @file Extractor.h
 * @brief This is the header file for the class that analyses a full tracker and material budget for export
 */

#ifndef _EXTRACTOR_H
#define	_EXTRACTOR_H

#include <tk2CMSSW_datatypes.h>
#include <tk2CMSSW_strings.h>
#include <set>
#include <cmath>
#include <sstream>
#include <tracker.hh>
#include <MaterialTable.h>
#include <MaterialBudget.h>

namespace insur {
    /**
     * @class Extractor
     * @brief This class bundles the analysis functions that prepare an existing material budget and table for output to CMSSW XML.
     *
     * The only public function of the class receives the material budget and table that make up the input. The output goes into an instance
     * of a struct - provided as a reference - that bundles a series of vectors of internal datatypes. ons of internal data types. The analysis
     * results will be stored in those, ready to be formatted and written to file.
     */
    class Extractor {
    public:
        void analyse(MaterialTable& mt, MaterialBudget& mb, CMSSWBundle& d, bool wt = false);
    protected:
        void analyseElements(MaterialTable&mattab, std::vector<Element>& elems);
        void analyseBarrelContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
                                                                           std::vector<std::pair<double, double> >& down);
        void analyseEndcapContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
                                                                             std::vector<std::pair<double, double> >& down);
        void analyseLayers(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& bc, Tracker& tr, std::vector<Composite>& c,
                                        std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<AlgoInfo>& a,
                                        std::vector<Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, bool wt = false);
        void analyseDiscs(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& ec, Tracker& tr, std::vector<Composite>& c,
                                      std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<AlgoInfo>& a,
                                      std::vector<Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, bool wt = false);
        void analyseBarrelServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                                                    std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt = false);
        void analyseEndcapServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                                                      std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt = false);
        void analyseSupports(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                                            std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt = false);
    private:
        Composite createComposite(std::string name, double density, MaterialProperties& mp, bool nosensors = false);
        std::vector<ModuleCap>::iterator findPartnerModule(std::vector<ModuleCap>::iterator i,
                                                                                                std::vector<ModuleCap>::iterator g, int ponrod, bool find_first = false);
        double findDeltaR(std::vector<Module*>::iterator start, std::vector<Module*>::iterator stop, double middle);
        double findDeltaZ(std::vector<Module*>::iterator start, std::vector<Module*>::iterator stop, double middle);
        int findSpecParIndex(std::vector<SpecParInfo>& specs, std::string label);
        double calculateSensorThickness(ModuleCap& mc, MaterialTable& mt);
        std::string stringParam(std::string name, std::string value);
        std::string numericParam(std::string name, std::string value);
        std::string vectorParam(double x, double y, double z);
        double compositeDensity(ModuleCap& mc, bool nosensors = false);
        double compositeDensity(InactiveElement& ie);
        double fromRim(double r, double w);
        int Z(double x0, double A);
    };
}
#endif	/* _EXTRACTOR_H */


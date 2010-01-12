// 
// File:   tk2CMSSW.h
// Author: ndemaio
//
// Created on August 28, 2009, 5:23 PM
//

/**
 * @file tk2CMSSW.h
 * @brief This is the header file for the CMSSW XML translator class
 */

#ifndef _TK2CMSSW_H
#define	_TK2CMSSW_H

#include <tk2CMSSW_strings.h>
#include <set>
#include <cmath>
#include <fstream>
#include <sstream>
#include <MaterialTable.h>
#include <MaterialBudget.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

namespace bfs = boost::filesystem;
namespace insur {
    class tk2CMSSW {
    public:
        tk2CMSSW() {}
        virtual ~tk2CMSSW() {}
        void translate(MaterialTable& mt, MaterialBudget& mb, std::string outsubdir = "");
    protected:
        enum CompType { wt, vl, ap };
        enum ShapeType { bx, tb, tp, pc };
        struct Rotation {
            std::string name;
            double phix;
            double phiy;
            double phiz;
            double thetax;
            double thetay;
            double thetaz;
        };
        struct Translation {
            double dx;
            double dy;
            double dz;
        };
        struct Element {
            std::string tag;
            double density;
            int atomic_number;
            double atomic_weight;
        };
        struct Composite {
            std::string name;
            double density;
            CompType method;
            std::vector<std::pair<std::string, double> > elements;
        };
        struct LogicalInfo {
            std::string name_tag;
            std::string shape_tag;
            std::string material_tag;
        };
        struct ShapeInfo {
            ShapeType type;
            std::string name_tag;
            double dx;
            double dxx;
            double dy;
            double dz;
            double rmin;
            double rmax;
            std::vector<std::pair<double, double> > rzup;
            std::vector<std::pair<double, double> > rzdown;
        };
        struct PosInfo {
            std::string parent_tag;
            std::string child_tag;
            int copy;
            Rotation rot;
            Translation trans;
        };
        struct AlgoInfo {
            std::string name;
            std::string parent;
            std::vector<std::string> parameters;
        };
        struct RingInfo {
            std::string name;
            std::string childname;
            bool fw;
            int modules;
            double rin;
            double rout;
            double rmid;
            double mthk;
            double phi;
        };
        struct SpecParInfo {
            std::string name;
            std::pair<std::string, std::string> parameter;
            std::vector<std::string> partselectors;
        };
        std::vector<Element> elements;
        std::vector<Composite> composites;
        std::vector<LogicalInfo> logic;
        std::vector<ShapeInfo> shapes;
        std::vector<PosInfo> positions;
        std::vector<AlgoInfo> algos;
        std::vector<SpecParInfo> specs;
        void materialSection(std::string name, std::vector<Element>& e, std::vector<Composite>& c, std::ostringstream& stream);
        void logicalPartSection(std::vector<LogicalInfo>& l, std::string label,  std::ostringstream& stream);
        void solidSection(std::vector<ShapeInfo>& s, std::string label, std::ostringstream& stream);
        void posPartSection(std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::string label, std::ostringstream& stream);
        void specParSection(std::vector<SpecParInfo>& t, std::string label, std::ostringstream& stream);
        void algorithm(std::string name, std::string parent, std::vector<std::string>& params, std::ostringstream& stream);
        void elementaryMaterial(std::string tag, double density, int a_number, double a_weight, std::ostringstream& stream);
        void compositeMaterial(std::string name, double density, CompType method,
                                               std::vector<std::pair<std::string, double> >& es, std::ostringstream& stream);
        void logicalPart(std::string name, std::string solid, std::string material, std::ostringstream& stream);
        void box(std::string name, double dx, double dy, double dz, std::ostringstream& stream);
        void trapezoid(std::string name, double dx, double dxx, double dy, double dz, std::ostringstream& stream);
        void tubs(std::string name, double rmin, double rmax, double dz, std::ostringstream& stream);
        void polycone(std::string name, std::vector<std::pair<double, double> >& rzu,
                               std::vector<std::pair<double, double> >& rzd, std::ostringstream& stream);
        void posPart(std::string parent, std::string child, Rotation& rot, Translation& trans, int copy, std::ostringstream& stream);
        void rotation(std::string name, double phix, double phiy, double phiz,
                                                          double thetax, double thetay, double thetaz, std::ostringstream& stream);
        void translation(double x, double y, double z, std::ostringstream& stream);
        std::string stringParam(std::string name, std::string value);
        std::string numericParam(std::string name, std::string value);
        std::string vectorParam(double x, double y, double z);
        void specPar(std::string name, std::pair<std::string, std::string> param, std::vector<std::string>& partsel, std::ostringstream& stream);
    private:
        void analyse(MaterialTable& mt, MaterialBudget& mb, std::vector<Element>& e,
                             std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                             std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::vector<SpecParInfo>& t);
        void analyseElements(MaterialTable&mattab, std::vector<Element>& elems);
        void analyseBarrelContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
                                                                           std::vector<std::pair<double, double> >& down);
        void analyseForwardEndcapContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
                                                                                           std::vector<std::pair<double, double> >& down);
        void analyseBackwardEndcapContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
                                                                                              std::vector<std::pair<double, double> >& down);
        void analyseEndcapContainer(std::vector<Module*>::iterator i, std::vector<Module*>::iterator g,
                                                         std::vector<std::pair<double, double> >& up, std::vector<std::pair<double, double> >& down);
        void analyseLayers(std::vector<std::vector<ModuleCap> >& bc, Tracker& tr,
                                        std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                                        std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::vector<SpecParInfo>& t);
        void analyseDiscs(std::vector<std::vector<ModuleCap> >& ec, Tracker& tr,
                                        std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                                        std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::vector<SpecParInfo>& t);
        void analyseBarrelServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                                                    std::vector<PosInfo>& p, std::vector<SpecParInfo>& t);
        void analyseEndcapServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                                                      std::vector<PosInfo>& p, std::vector<SpecParInfo>& t);
        void analyseSupports(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                                            std::vector<PosInfo>& p, std::vector<SpecParInfo>& t);
        tk2CMSSW::Composite createComposite(std::string name, double density, MaterialProperties& mp);
        std::vector<ModuleCap>::iterator findPartnerModule(std::vector<ModuleCap>::iterator i,
                                                                                                std::vector<ModuleCap>::iterator g, int ponrod, bool find_first = false);
        double findDeltaR(std::vector<Module*>::iterator start, std::vector<Module*>::iterator stop, double middle);
        double findDeltaZ(std::vector<Module*>::iterator start, std::vector<Module*>::iterator stop, double middle);
        int findSpecParIndex(std::vector<SpecParInfo>& specs, std::string label);
        double fromRim(double r, double w);
        double compositeDensity(InactiveElement& ie);
        double compositeDensity(ModuleCap& mc);
        int Z(double x0, double A);
        void print();
    };
}
#endif	/* _TK2CMSSW_H */


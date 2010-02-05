// 
// File:   XMLWriter.h
// Author: ndemaio
//
// Created on January 31, 2010, 3:33 AM
//

/**
 * @file XMLWriter.h
 * @brief This is the header file for the class that writes a set of extracted tracker data to CMSSW XML
 */

#ifndef _XMLWRITER_H
#define	_XMLWRITER_H

#include <tk2CMSSW_datatypes.h>
#include <tk2CMSSW_strings.h>
#include <iostream>
#include <sstream>
#include <fstream>

namespace insur {
    class XMLWriter {
    public:
        void writeXML(std::vector<Element> e, std::vector<Composite> c, std::vector<LogicalInfo> l, std::vector<ShapeInfo> s,
                                 std::vector<PosInfo> p, std::vector<AlgoInfo> a, std::vector<SpecParInfo> t, std::ostringstream& gb,
                                 std::ostringstream& tb, std::ostringstream& pb, std::ostringstream& sb, std::ostringstream& mb);
    protected:
        void prodcuts(std::vector<SpecParInfo>& t, std::ostringstream& stream);
        void trackersens(std::vector<SpecParInfo>& t, std::ostringstream& stream);
        void recomaterial(std::vector<SpecParInfo>& t, std::ostringstream& stream);
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
        void specPar(std::string name, std::pair<std::string, std::string> param, std::vector<std::string>& partsel, std::ostringstream& stream);
    private:
        std::vector<std::pair<std::string, std::vector<std::string> > >& buildPaths(std::vector<SpecParInfo>& specs,
                                                                                                               std::vector<std::pair<std::string, std::vector<std::string> > >& blocks);
        bool endcapsInTopology(std::vector<SpecParInfo>& specs);
        int findNumericPrefixSize(std::string s);
        int findEntry(std::vector<SpecParInfo>& specs, std::string name);
        std::vector<std::pair<std::string, std::vector<std::string> > >::iterator findEntry(std::string name,
                                                                                                                          std::vector<std::pair<std::string, std::vector<std::string> > >& data);
    };
}
#endif	/* _XMLWRITER_H */


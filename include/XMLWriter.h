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
#include <global_constants.h>
#include <iostream>
#include <sstream>
#include <fstream>

namespace insur {
    /**
     * @class XMLWriter
     * @brief This class bundles the output functions that write tracker information from a series of internal collections to CMSSW XML files.
     *
     * The public functions of the class produce the files <i>pixbar.xml, pixfwd.xml, tracker.xml, trackerStructureTopology.xml,
     * trackerProdCuts.xml, trackersens.xml</i> and <i>trackerRecoMaterial.xml</i> from a series of data collections and
     * writes them to the provided output streams. The file <i>tracker.xml</i> is created from scratch while the others use skeleton
     * files that they add to.
     */
    class XMLWriter {
    public:
        void pixbar(std::vector<ShapeInfo>& s, std::ifstream& in, std::ofstream& out);
        void pixfwd(std::vector<ShapeInfo>& s, std::ifstream& in, std::ofstream& out);
        void tracker(CMSSWBundle& d, std::ofstream& out, std::istream& trackerVolumeTemplate, bool wt = false);
        void topology(std::vector<SpecParInfo>& t, std::ifstream& in, std::ofstream& out);
        void prodcuts(std::vector<SpecParInfo>& t, std::ifstream& in, std::ofstream& out);
        void trackersens(std::vector<SpecParInfo>& t, std::ifstream& in, std::ofstream& out);
        void recomaterial(std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, std::ifstream& in, std::ofstream& out, bool wt = false);
        void setExtendedHeader(const std::string& header) { extendedHeader_ = header; }
        void setSimpleHeader(const std::string& header) { simpleHeader_ = header; }
        const std::string& getExtendedHeader() const { return extendedHeader_; }
        const std::string& getSimpleHeader() const { return simpleHeader_; }
    protected:
        void trackerLogicalVolume(std::ostringstream& stream, std::istream& instream); // takes the stream containing the tracker logical volume template and outputs it to the outstream
        void materialSection(std::string name, std::vector<Element>& e, std::vector<Composite>& c, std::ostringstream& stream);
        void rotationSection(std::map<std::string,Rotation>& r, std::string label, std::ostringstream& stream);
        void logicalPartSection(std::vector<LogicalInfo>& l, std::string label,  std::ostringstream& stream, bool wt = false);
        void solidSection(std::vector<ShapeInfo>& s, std::vector<ShapeOperationInfo>& so, std::string label, std::ostringstream& stream, std::istream& trackerVolumeTemplate, bool notobtid, bool wt = false);
        void posPartSection(std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::string label, std::ostringstream& stream);
        void specParSection(std::vector<SpecParInfo>& t, std::string label, std::ostringstream& stream);
        void algorithm(std::string name, std::string parent, std::vector<std::string>& params, std::ostringstream& stream);
        void elementaryMaterial(std::string tag, double density, int a_number, double a_weight, std::ostringstream& stream);
        void compositeMaterial(std::string name, double density, CompType method,
                                               std::vector<std::pair<std::string, double> >& es, std::ostringstream& stream);
        void logicalPart(std::string name, std::string solid, std::string material, std::ostringstream& stream);
        void box(std::string name, double dx, double dy, double dz, std::ostringstream& stream);
        void trapezoid(std::string name, double dx, double dxx, double dy, double dyy, double dz, std::ostringstream& stream);
        void tubs(std::string name, double rmin, double rmax, double dz, std::ostringstream& stream);
	void cone(std::string name, double rmin1, double rmax1, double rmin2, double rmax2, double dz, std::ostringstream& stream);
        void polycone(std::string name, std::vector<std::pair<double, double> >& rzu,
                               std::vector<std::pair<double, double> >& rzd, std::ostringstream& stream);
	void sunion(std::string name, std::string rSolid1, std::string rSolid2, std::ostringstream& stream);
	void sintersection(std::string name, std::string rSolid1, std::string rSolid2, std::ostringstream& stream);
        void posPart(std::string parent, std::string child, std::string rotref, Translation& trans, int copy, std::ostringstream& stream);
        void rotation(std::string name, double thetax, double phix, double thetay, double phiy,
                                                          double thetaz, double phiz, std::ostringstream& stream);
        void translation(double x, double y, double z, std::ostringstream& stream);
        void specPar(std::string name, std::pair<std::string, std::string> param, std::vector<std::string>& partsel, std::ostringstream& stream);
        void specPar1(std::string name, std::pair<std::string, std::string> param, std::vector<std::string>& partsel, std::ostringstream& stream);
        void specParROC(std::vector<std::string>& partsel, std::vector<ModuleROCInfo>& minfo, std::pair<std::string, std::string> param, std::ofstream& stream);
    private:
        std::vector<PathInfo>& buildPaths(std::vector<SpecParInfo>& specs, std::vector<PathInfo>& blocks, bool wt = false);
        bool endcapsInTopology(std::vector<SpecParInfo>& specs);
        int findNumericPrefixSize(std::string s);
        int findEntry(std::vector<SpecParInfo>& specs, std::string name);
        std::vector<PathInfo>::iterator findEntry(std::string name, std::vector<PathInfo>& data);
        std::string extendedHeader_, simpleHeader_; // headers containing generation information which are inserted after the preamble
    };
}
#endif	/* _XMLWRITER_H */


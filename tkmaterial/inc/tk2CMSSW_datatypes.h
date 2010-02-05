#ifndef _TK2CMSSW_DATATYPES_H
#define _TK2CMSSW_DATATYPES_H

/**
 * @file tk2CMSSW_datatypes.h
 * @brief This header file lists the internal structs and enumerations needed to translate a full tracker and material budget to CMSSW XML
 */

#include <string>
#include <vector>

namespace insur {
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
}
#endif /* _TK2CMSSW_DATATYPES_H */

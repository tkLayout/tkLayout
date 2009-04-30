//
// File:   MatCalc.h
// Author: ndemaio
//
// Created on October 24, 2008, 11:24 AM
//

/**
 * @file MatCalc.h
 * @brief
 */

#ifndef _MATCALC_H
#define	_MATCALC_H

#include <list>
#include <string>
#include <stdexcept>
#include <ModuleCap.h>
#include <MaterialTable.h>
#include <InactiveElement.h>
namespace insur {
    static const std::string err_no_such_type = "Error: the requested type is not on the list.";
    static const std::string err_no_material = "Error: no material with the specified tag was found.";
    static const std::string err_no_service = "Error: no service material with the specified properties was found.";
    static const std::string err_no_support = "Error: no support material with the specified properties was found.";
    static const std::string err_unknown_type = "Error: unknown module type.";
    static const std::string err_up_general = "Error updating parameter entry.";
    static const std::string err_conversion = "Error: material unit not recognised during conversion.";
    static const std::string err_matadd_weird = "Something weird happened when trying to add an entry to one of the vectors for material parameters...";
    static const std::string msg_negative_area = "Warning: module surface is negative.";
    static const std::string msg_abort = "Aborting function.";
    /**
     * @class MatCalc
     * @brief
     */
    class MatCalc {
    public:
        enum Modtype { un_mod, rphi, stereo, pt };
        enum Matunit { gr, mm3, mm, grpm };
        MatCalc() { init_done = false; }
        virtual ~MatCalc() {}
        bool initDone();
        void initDone(bool yes);
        int getStripsAcross(Modtype type); // throws exception
        int getSegmentsAlong(Modtype type); // throws exception
        std::pair<int, int> getDefaultDimensions(Modtype type); // throws exception
        void addTypeInfo(Modtype type, int strips, int segments);
        void updateTypeInfoStrips(Modtype type, int strips);
        void updateTypeInfoSegments(Modtype type, int segments);
        void addModuleParameters(std::string tag, Modtype type,
                double A, Matunit uA, double B, Matunit uB, double C, Matunit uC, double D, Matunit uD, bool local);
        void addServiceParameters(std::string tag, double Q, Matunit uQ);
        void addServiceParameters(std::string tagIn, double In, Matunit uIn, std::string tagOut, double Out, Matunit uOut, bool local);
        void addSupportParameters(std::string tag, double M, Matunit uM, MaterialProperties::Category cM);
        void clearTypeVector();
        void clearModVectors();
        void clearModVector(Modtype type);
        void copyContents(Modtype source, Modtype dest);
        void appendContents(Modtype source, Modtype dest);
        // TODO: add vector operations for services and supports
        bool typeRegistered(Modtype type);
        uint registeredTypes();
        MaterialTable& getMaterialTable();
        // TODO: add more functions as necessary
        virtual bool calculateBarrelMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps);
        virtual bool calculateEndcapMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps);
        virtual bool calculateBarrelServiceMaterials(
            std::vector<std::vector<ModuleCap> >& barrelcaps,
                std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices);
        virtual bool calculateEndcapServiceMaterials(
            std::vector<std::vector<ModuleCap> >& endcapcaps,
                std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices);
        virtual bool calculateSupportMaterials(std::vector<InactiveElement>& supports);
        void printInternals();
    protected:
        struct TypeInfo {
            Modtype type;
            int strips_across;
            int segments_along;
        };
        struct SingleMod {
            std::string tag;
            double A, B, C, D;
            Matunit uA, uB, uC, uD;
            bool is_local;
        };
        struct SingleSerLocal {
            std::string tag;
            double Q;
            Matunit uQ;
        };
        struct SingleSerExit {
            std::string tagIn;
            double In;
            Matunit uIn;
            std::string tagOut;
            double Out;
            Matunit uOut;
            bool is_local;
        };
        struct SingleSup {
            std::string tag;
            double M;
            Matunit uM;
            MaterialProperties::Category cM;
        };
        struct MatInfo {
            std::vector<TypeInfo> typeinfo;
            std::vector<SingleMod> modinforphi;
            std::vector<SingleMod> modinfostereo;
            std::vector<SingleMod> modinfopt;
            std::vector<SingleSerLocal> serlocalinfo;
            std::vector<SingleSerExit> serexitinfo;
            std::vector<SingleSup> supinfo;
        };
        bool init_done;
        MaterialTable mt;
        MatInfo internals;
        std::vector<SingleMod>& getModVector(Modtype type); // throws exception
        TypeInfo& getTypeInfoByType(Modtype type); // throws exception
        SingleMod& getSingleMod(std::string tag, Modtype type, Matunit uA, Matunit uB, Matunit uC, Matunit uD, bool local); // throws exception
        SingleSerLocal& getSingleSer(std::string tag, Matunit u); // throws exception
        SingleSerExit& getSingleSer(std::string tag1, std::string tag2, Matunit u1, Matunit u2, bool local); // throws exception
        SingleSup& getSingleSup(std::string tag, Matunit uM, MaterialProperties::Category cM); // throws exception
    private:
        bool entryExists(Modtype type);
        bool entryExists(std::string tag, Modtype type, Matunit uA, Matunit uB, Matunit uC, Matunit uD, bool local);
        bool entryExists(std:: string tag, Matunit uQ);
        bool entryExists(std::string tag1, std::string tag2, Matunit uIn, Matunit uOut, bool local);
        bool entryExists(std::string tag, Matunit uM, MaterialProperties::Category cM);
        int findRods(std::vector<std::vector<ModuleCap> >& caps, int layer);
        double convert(double value, Matunit unit, double densityorlength, double surface = 0); // throws exception
        void adjacentDifferentCategory(std::vector<ModuleCap>& source, InactiveElement& dest, int r, double l, double s);
        void adjacentSameCategory(InactiveElement& source, InactiveElement& dest);
    };
}
#endif	/* _MATCALC_H */


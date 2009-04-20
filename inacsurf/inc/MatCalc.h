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

#include <string>
#include <stdexcept>
#include <ModuleCap.h>
#include <MaterialTable.h>
#include <InactiveElement.h>
namespace insur {
    static const std::string err_no_such_type = "Error: the requested type is not on the list.";
    static const std::string err_no_material = "Error: no material with the specified tag was found.";
    static const std::string err_no_service = "Error: no service with the specified tag was found.";
    static const std::string err_unknown_type = "Error: unknown module type.";
    static const std::string err_up_general = "Error updating parameter entry.";
    static const std::string err_matadd_weird = "Something weird happened when trying to add an entry to one of the vectors for material parameters...";
    /**
     * @class MatCalc
     * @brief
     */
    class MatCalc {
    public:
        enum Modtype { rphi, stereo, pt };
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
        // TODO: add more functions as necessary
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
    private:
        bool entryExists(Modtype type);
        bool entryExists(std::string tag, Modtype type, Matunit uA, Matunit uB, Matunit uC, Matunit uD, bool local);
        bool entryExists(std:: string tag, Matunit uQ);
        bool entryExists(std::string tag1, std::string tag2, Matunit uIn, Matunit uOut, bool local);
        int findRods(std::vector<std::vector<ModuleCap> >& caps, int layer);
        double convert(double value, Matunit unit, double densityorlength, double surface = 0);
    };
}
#endif	/* _MATCALC_H */


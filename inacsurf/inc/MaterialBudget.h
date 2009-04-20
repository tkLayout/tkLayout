//
// File:   MaterialBudget.h
// Author: ndemaio
//
// Created on March 10, 2009, 6:40 PM
//

/**
 * @file MaterialBudget.h
 * @brief
 */

#ifndef _MATERIALBUDGET_H
#define	_MATERIALBUDGET_H

#include <vector>
#include <tracker.hh>
#include <InactiveSurfaces.h>
#include <ModuleCap.h>
#include <MatCalc.h>
namespace insur {
    static const std::string err_materials_supports = "Error calculating materials for supports.";
    static const std::string err_materials_bservices = "Error calculating materials for barrel services.";
    static const std::string err_materials_eservices = "Error calculating materials for endcap services.";
    static const std::string err_materials_bmodules = "Error calculating materials for barrel modules.";
    static const std::string err_materials_emodules = "Error calculating materials for endcap modules.";
    
    class MaterialBudget {
    public:
        MaterialBudget(Tracker& tr, InactiveSurfaces& is);
        virtual ~MaterialBudget();
        void materialsAll(MatCalc& calc);
        void print();
    protected:
        Tracker* tracker;
        InactiveSurfaces* inactive;
        std::vector<std::vector<ModuleCap> > capsbarrelmods, capsendmods;
        void materialsSupports(MatCalc& calc);
        void materialsServices(MatCalc& calc);
        void materialsModules(MatCalc& calc);
    private:
        MaterialBudget();
        MaterialBudget(const MaterialBudget& budget);
    };
}
#endif	/* _MATERIALBUDGET_H */


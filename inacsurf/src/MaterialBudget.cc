/**
 * @file MaterialBudget.cc
 * @brief
 */

#include <MaterialBudget.h>
namespace insur {
    MaterialBudget::MaterialBudget(Tracker& tr, InactiveSurfaces& is) {
        tracker = &tr;
        inactive = &is;
        std::vector<Layer*>::const_iterator liter = tracker->getBarrelLayers()->begin();
        std::vector<Layer*>::const_iterator lend = tracker->getBarrelLayers()->end();
        std::vector<Module*>::const_iterator miter;
        std::vector<Module*>::const_iterator mend;
        // barrel layers
        while (liter != lend) {
            std::vector<ModuleCap> tmp;
            capsbarrelmods.push_back(tmp);
            miter = (*liter)->getModuleVector()->begin();
            mend = (*liter)->getModuleVector()->end();
            while (miter != mend) {
                ModuleCap cap(*(*miter));
                capsbarrelmods.back().push_back(cap);
                miter++;
            }
            liter++;
        }
        // endcap layers
        liter = tracker->getEndcapLayers()->begin();
        lend = tracker->getEndcapLayers()->end();
        while (liter != lend) {
            std::vector<ModuleCap> tmp;
            capsendmods.push_back(tmp);
            miter = (*liter)->getModuleVector()->begin();
            mend = (*liter)->getModuleVector()->end();
            while (miter != mend) {
                ModuleCap cap(*(*miter));
                capsendmods.back().push_back(cap);
                miter++;
            }
            liter++;
        }
    }
    
    MaterialBudget::~MaterialBudget() {}
    
    void MaterialBudget::materialsAll(MatCalc& calc) {
        materialsModules(calc);
        materialsServices(calc);
        materialsSupports(calc);
    }
    
    void MaterialBudget::materialsSupports(MatCalc& calc) {
        if (!calc.calculateSupportMaterials(inactive->getSupports()))
            std::cout << err_materials_supports << std::endl;
    }
    
    void MaterialBudget::materialsServices(MatCalc& calc) {
        if (calc.calculateBarrelServiceMaterials(capsbarrelmods, inactive->getBarrelServices())) {
            if (!calc.calculateEndcapServiceMaterials(capsendmods, inactive->getEndcapServices()))
                std::cout << err_materials_eservices << std::endl;
        }
        else std::cout << err_materials_bservices << std::endl;
    }
    void MaterialBudget::materialsModules(MatCalc& calc) {
        if (calc.calculateBarrelMaterials(capsbarrelmods)) {
            if (!calc.calculateEndcapMaterials(capsendmods))
                std::cout << err_materials_emodules << std::endl;
        }
        else std::cout << err_materials_bmodules << std::endl;
    }
    
    void MaterialBudget::print() {
        // TODO: print statements
    }
}

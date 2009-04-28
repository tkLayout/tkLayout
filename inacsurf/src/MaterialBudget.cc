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
                cap.setCategory(MaterialProperties::b_mod);
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
                cap.setCategory(MaterialProperties::e_mod);
                capsendmods.back().push_back(cap);
                miter++;
            }
            liter++;
        }
    }
    
    MaterialBudget::~MaterialBudget() {}
    
    Tracker& MaterialBudget::getTracker() { return *tracker; }
    
    InactiveSurfaces& MaterialBudget::getInactiveSurfaces() { return *inactive; }
    
    std::vector<std::vector<ModuleCap> >& MaterialBudget::getBarrelModuleCaps() { return capsbarrelmods; }
    
    std::vector<std::vector<ModuleCap> >& MaterialBudget::getEndcapModuleCaps() { return capsendmods; }
    
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
        if (calc.calculateBarrelServiceMaterials(capsbarrelmods, inactive->getBarrelServices(), inactive->getEndcapServices())) {
            if (!calc.calculateEndcapServiceMaterials(capsendmods, inactive->getBarrelServices(), inactive->getEndcapServices()))
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
        std::cout << "-----Material Budget Internal State-----" << std::endl;
        int a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
        for (uint i = 0; i < inactive->getSupports().size(); i++) {
            if (inactive->getSupportPart(i).isVertical()) a++;
            else b++;
        }
        for (uint i = 0; i < tracker->getBarrelLayers()->size(); i++) c = c + tracker->getBarrelLayers()->at(i)->getModuleVector()->size();
        for (uint i = 0; i < capsbarrelmods.size(); i++) d = d + capsbarrelmods.at(i).size();
        for (uint i = 0; i < tracker->getEndcapLayers()->size(); i++) e = e + tracker->getEndcapLayers()->at(i)->getModuleVector()->size();
        for (uint i = 0; i < capsendmods.size(); i++) f = f + capsendmods.at(i).size();
        std::cout << "Tracker: " << c << " barrel modules in " << tracker->getBarrelLayers()->size() << " layers." << std::endl;
        std::cout << "MaterialBudget: " << d << " barrel modcaps in " << capsbarrelmods.size() << " vectors." << std::endl;
        std::cout << "Tracker: " << e << " endcap modules in " << tracker->getEndcapLayers()->size() << " discs." << std::endl;
        std::cout << "MaterialBudget: " << f << " endcap modcaps in " << capsendmods.size() << " vectors." << std::endl;
        std::cout << "InactiveServices: " << inactive->getBarrelServices().size() << " barrel services and ";
        std::cout << inactive->getEndcapServices().size() << " endcap services." << std::endl;
        std::cout << "InactiveSurfaces: " << inactive->getSupports().size() << " supports in total, "  << a << " in barrels, " << b << " in endcaps." << std::endl;
        int am = 0, bm = 0, cm = 0, dm = 0, em = 0, fm = 0;
        double ar = 0, br = 0, cr = 0, dr = 0, er = 0, fr = 0, ai = 0, bi = 0, ci = 0, di = 0, ei = 0, fi = 0;
        for (uint i = 0; i < capsbarrelmods.size(); i++) {
            for (uint j = 0; j < capsbarrelmods.at(i).size(); j++) {
                am = am + capsbarrelmods.at(i).at(j).getTotalMass();
                ar = ar + capsbarrelmods.at(i).at(j).getRadiationLength();
                ai = ai + capsbarrelmods.at(i).at(j).getInteractionLength();
            }
        }
        std::cout << "Total module mass in barrels: " << am << std::endl;
        std::cout << "Average module mass in barrels: " << (am / (double)d) << std::endl;
        std::cout << "Average module radiation length in barrels: " << (ar / (double)d) << std::endl;
        std::cout << "Average module interaction length in barrels: " << (ai / (double)d) << std::endl;
        for (uint i = 0; i < capsendmods.size(); i++) {
            for (uint j = 0; j < capsendmods.at(i).size(); j++) {
                bm = bm + capsendmods.at(i).at(j).getTotalMass();
                br = br + capsendmods.at(i).at(j).getRadiationLength();
                bi = bi + capsendmods.at(i).at(j).getInteractionLength();
            }
        }
        std::cout << "Total module mass in endcaps: " << bm << std::endl;
        std::cout << "Average module mass in endcaps: " << (bm / (double)f) << std::endl;
        std::cout << "Average module radiation length in endcaps: " << (br / (double)f) << std::endl;
        std::cout << "Average module interaction length in endcaps: " << (bi / (double)f) << std::endl;
        for (uint i = 0; i < inactive->getBarrelServices().size(); i++) {
            cm = cm + inactive->getBarrelServicePart(i).getTotalMass();
            cr = cr + inactive->getBarrelServicePart(i).getRadiationLength();
            ci = ci + inactive->getBarrelServicePart(i).getInteractionLength();
        }
        std::cout << "Total service mass in barrels: " << cm << std::endl;
        std::cout << "Average service mass in barrels: " << (cm / (double)(inactive->getBarrelServices().size())) << std::endl;
        std::cout << "Average service radiation length in barrels: " << (cr / (double)(inactive->getBarrelServices().size())) << std::endl;
        std::cout << "Average service interaction length in barrels: " << (ci / (double)(inactive->getBarrelServices().size())) << std::endl;
        for (uint i = 0; i < inactive->getEndcapServices().size(); i++) {
            dm = dm + inactive->getEndcapServicePart(i).getTotalMass();
            dr = dr + inactive->getEndcapServicePart(i).getRadiationLength();
            di = di + inactive->getEndcapServicePart(i).getInteractionLength();
        }
        std::cout << "Total service mass in endcaps: " << dm << std::endl;
        std::cout << "Average service mass in endcaps: " << (dm / (double)(inactive->getEndcapServices().size())) << std::endl;
        std::cout << "Average service radiation length in endcaps: " << (dr / (double)(inactive->getEndcapServices().size())) << std::endl;
        std::cout << "Average service interaction length in endcaps: " << (di / (double)(inactive->getEndcapServices().size())) << std::endl;
        for (uint i = 0; i < inactive->getSupports().size(); i++) {
            if (inactive->getSupportPart(i).isVertical()) {
                em = em + inactive->getSupportPart(i).getTotalMass();
                er = er + inactive->getSupportPart(i).getRadiationLength();
                ei = ei + inactive->getSupportPart(i).getInteractionLength();
            }
            else {
                fm = fm + inactive->getSupportPart(i).getTotalMass();
                fr = fr + inactive->getSupportPart(i).getRadiationLength();
                fi = fi + inactive->getSupportPart(i).getInteractionLength();
            }
        }
        std::cout << "Total support mass in barrels: " << em << std::endl;
        std::cout << "Average support mass in barrels: " << (em / (double)a) << std::endl;
        std::cout << "Average support radiation length in barrels: " << (er / (double)a) << std::endl;
        std::cout << "Average support interaction length in barrels: " << (ei / (double)a) << std::endl;
        std::cout << "Total support mass in endcaps: " << fm << std::endl;
        std::cout << "Average support mass in endcaps: " << (fm / (double)b) << std::endl;
        std::cout << "Average support radiation length in endcaps: " << (fr / (double)b) << std::endl;
        std::cout << "Average support interaction length in endcaps: " << (fi / (double)b) << std::endl;
        std::cout << "Total mass in barrels: " << (am + cm + em) << std::endl;
        std::cout << "Total mass in endcaps: " << (bm + dm + fm) << std::endl;
        std::cout << "Average radiation length in barrels: " << (((ar / (double)d) + (cr / (double)(inactive->getBarrelServices().size())) + (er / (double)a)) / 3.0) << std::endl;
        std::cout << "Average radiation length in endcaps: " << (((br / (double)f) + (dr / (double)(inactive->getEndcapServices().size())) + (fr / (double) b)) / 3.0) << std::endl;
        std::cout << "Average interaction length in barrels: " << (((ai / (double)d) + (ci / (double)(inactive->getBarrelServices().size())) + (ei / (double)a)) / 3.0) << std::endl;
        std::cout << "Average interaction length in endcaps: " << (((bi / (double)f) + (di / (double)(inactive->getEndcapServices().size())) + (fi / (double)b)) / 3.0) << std::endl;
        std::cout << std::endl;
    }
}

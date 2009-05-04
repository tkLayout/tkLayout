#include <MatCalcDummy.h>
namespace insur {
    MatCalcDummy::MatCalcDummy() { init(); }
    
    bool MatCalcDummy::calculateBarrelMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps) {
        for (unsigned int i = 0; i < barrelcaps.size(); i++) {
            for (unsigned int j = 0; j < barrelcaps.at(i).size(); j++) {
                if (barrelcaps.at(i).at(j).getSurface() <= 0) {
                    std::cout << "MatCalcDummy::calculateBarrelMaterials(): Warning: surface of barrelcap reported as " << barrelcaps.at(i).at(j).getSurface() << std::endl;
                }
                barrelcaps.at(i).at(j).addLocalMass(ta, ab);
                barrelcaps.at(i).at(j).calculateTotalMass();
                barrelcaps.at(i).at(j).calculateRadiationLength(mt);
                barrelcaps.at(i).at(j).calculateInteractionLength(mt);
            }
        }
        return true;
    }
    
    bool MatCalcDummy::calculateEndcapMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps) {
        for (unsigned int i = 0; i < endcapcaps.size(); i++) {
            for (unsigned int j = 0; j < endcapcaps.at(i).size(); j++) {
                if (endcapcaps.at(i).at(j).getSurface() <= 0) {
                    std::cout << "MatCalcDummy::calculateEndcapMaterials(): Warning: surface of endcapcap reported as " << endcapcaps.at(i).at(j).getSurface() << std::endl;
                }
                endcapcaps.at(i).at(j).addLocalMass(ta, ae);
                endcapcaps.at(i).at(j).calculateTotalMass();
                endcapcaps.at(i).at(j).calculateRadiationLength(mt);
                endcapcaps.at(i).at(j).calculateInteractionLength(mt);
            }
        }
        return true;
    }
    
    bool MatCalcDummy::calculateBarrelServiceMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps,
                            std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices) {
        for (unsigned int i = 0; i < barrelservices.size(); i++) {
            barrelservices.at(i).addLocalMass(ts, sb);
            barrelservices.at(i).calculateTotalMass();
            barrelservices.at(i).calculateRadiationLength(mt);
            barrelservices.at(i).calculateInteractionLength(mt);
        }
        return true;
    }
    
    bool MatCalcDummy::calculateEndcapServiceMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps,
                            std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices) {
        for (unsigned int i = 0; i < endcapservices.size(); i++) {
            endcapservices.at(i).addLocalMass(ts, se);
            endcapservices.at(i).calculateTotalMass();
            endcapservices.at(i).calculateRadiationLength(mt);
            endcapservices.at(i).calculateInteractionLength(mt);
        }
        return true;
    }
    
    bool MatCalcDummy::calculateSupportMaterials(std::vector<InactiveElement>& supports) {
        for (unsigned int i = 0; i < supports.size(); i++) {
            if (supports.at(i).isVertical()) supports.at(i).addLocalMass(tl, lb);
            else supports.at(i).addLocalMass(tl, le);
            supports.at(i).calculateTotalMass();
            supports.at(i).calculateRadiationLength(mt);
            supports.at(i).calculateInteractionLength(mt);
        }
        return true;
    }
    
    void MatCalcDummy::init() {
        ta = "Si_dummy";
        ts = "Cu_dummy";
        tl = "CF_dummy";
        ab = 100.0;
        ae = 200.0;
        sb = 500.0;
        se = 1000.0;
        lb = 2500.0;
        le = 2000.0;
        mt.addMaterial(ta, d_silicon, 1, 1.5);
        mt.addMaterial(ts, d_copper, 2, 2.5);
        mt.addMaterial(tl, d_carbon, 3, 3.5);
    }
}

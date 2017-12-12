/**
 * @file MatCalcDummy.cc
 * @brief This is the implementation of the derivative material assignment test class
 */

#include <MatCalcDummy.hh>
namespace insur {
    /**
     * The constructor simply calls the init function...
     */
    MatCalcDummy::MatCalcDummy() { init(); }
    
    /**
     * The function to assign materials to the barrel modules simply puts the same static amount in each one. It overrides the
     * more complex version of the parent class with something quick and dirty.
     * @param barrelcaps A reference to the collection of <i>ModuleCap</i> objects that sit on top of the barrel modules
     * @return True, because the process is always considered successful - even if one of the module surfaces was reported as negative
     */
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
    
    
    /**
     * The function to assign materials to the endcap modules simply puts the same static amount in each one. It overrides
     * the more complex version of the parent class with something quick and dirty.
     * @param endcapcaps A reference to the collection of <i>ModuleCap</i> objects that sit on top of the endcap modules
     * @return True, because the process is always considered successful - even if one of the module surfaces was reported as negative
     */
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
    
    /**
     * The function to assign materials to the barrel services simply puts the same static amount in each one. It overrides
     * the more complex version of the parent class with something quick and dirty.
     * @param barrelcaps A reference to the collection of <i>ModuleCap</i> that sits on top of the barrel modules; unused
     * @param barrelservices A reference to the collection of barrel services
     * @param endcapservices A reference to the collection of endcap services; unused
     * @return True, because the process is always considered successful
     */
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
    
    /**
     * The function to assign materials to the endcap services simply puts the same static amount in each one. It overrides
     * the more complex version of the parent class with something quick and dirty.
     * @param endcapcaps A reference to the collection of <i>ModuleCap</i> that sits on top of the endcap modules; unused
     * @param barrelservices A reference to the collection of barrel services; unused
     * @param endcapservices A reference to the collection of endcap services
     * @return True, because the process is always considered successful
     */
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
    
    /**
     * The function to assign materials to the supports simply puts the same static amount in each one. It overrides the 
     * more complex version of the parent class with something quick and dirty.
     * @param supports A reference to the collection of support parts
     * @return True, because the process is always considered successful
     */
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
    
    /**
     * This function initialises the material variables and the material table to the dummy values necessary to assign
     * test materials and calculate the radiation and interaction lengths from them.
     */
    void MatCalcDummy::init() {
        ta = "Si_dummy";
        ts = "Cu_dummy";
        tl = "CF_dummy";
        // arbitrary numbers
        ab = 10.0;
        ae = 10.0;
        sb = 500.0;
        se = 500.0;
        lb = 2000.0;
        le = 2000.0;
        mt.addMaterial(ta, mat_d_silicon, 1, 1.5);
        mt.addMaterial(ts, mat_d_copper, 2, 2.5);
        mt.addMaterial(tl, mat_d_carbon, 3, 3.5);
    }
}

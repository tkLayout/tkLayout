/**
 * @file MaterialProperties.cc
 * @brief This is the implementation for a base class for elements that have an influence on the material budget of a tracker layout
 */

#include <MaterialProperties.hh>
#include<MaterialTab.hh>

RILength& RILength::operator+=(const RILength &a) {
  interaction += a.interaction;
  radiation += a.radiation;
  return *this;
}

const RILength RILength::operator+(const RILength &other) const {
  RILength result = *this;     // Make a copy of myself.  Same as Material result(*this);
  result += other;            // Use += to add other to the copy.
  return result;              // All done!
}



namespace insur {
    /*-----public functions-----*/
    /**
     * The constructor sets a few defaults. The flags for the initialisation status of the material vectors are
     * set to false, the element category to none, i.e. unidentified, and the numeric values for the sums of masses
     * and the radiation and interaction lengths to -1 to indicate that they are all uninitialised.
     */
    MaterialProperties::MaterialProperties() {
        msl_set = false;
        trck = true;
        cat = no_cat;
        total_mass = 0;
        local_mass = 0;
        r_length = 0;
        i_length = 0;
    }
    
    /**
     * Get the category of this element.
     * @return The category identifier as defined in the enumeration <i>Category</i>
     */
    MaterialProperties::Category MaterialProperties::getCategory() { return cat; }
    
    /**
     * Set the category of this element.
     * @param c The category identifier as defined in the enumeration <i>Category</i>
     */
    void MaterialProperties::setCategory(Category c) { cat = c; }
    
    /**
     * Get the surface of this element that is relevant for tracking. This is a virtual function that is meant to be
     * implemented by the subclasses since the material properties as collected in this class know nothing about
     * the geometry of the element they belong to. As soon as a subclass contains the relevant information, it is
     * expected to provide a meaningful value that this function can return.
     * @return A constant value of <i>-1</i>
     */
    double MaterialProperties::getSurface() const { return -1; }

    double MaterialProperties::getLength() const { return -1; }
    
    /**
     * Get the local mass of one of the materials, as identified by its name, that make up the element.
     * If the material does not appear on the list, the function throws an exception.
     * @param tag The name of the material
     * @return The mass of the requested material
     */
    double MaterialProperties::getLocalMass(std::string tag) { // throws exception
        if (!localmasses.count(tag)) throw std::runtime_error("MaterialProperties::getLocalMass(std::string): " + err_local_mass + ": " + tag);
        return localmasses.at(tag);
    }

    /**
     * Get the local mass of one of the components, as identified by its name, that make up the element.
     * If the component does not appear on the list, the function throws an exception.
     * @param tag The name of the component
     * @return The mass of the requested component
     */
    double MaterialProperties::getLocalMassComp(std::string comp) { // throws exception
        if (!localmassesComp.count(comp)) throw std::runtime_error("MaterialProperties::getLocalMass(std::string): " + err_local_mass + ": " + comp);
        return localmassesComp.at(comp);
    }
    
    
    const std::map<std::string, double>& MaterialProperties::getLocalMasses() const { return localmasses; }
    const std::map<std::string, double>& MaterialProperties::getLocalMassesComp() const { return localmassesComp; }

    /**
     * Add the local mass for a material, as specified by its tag, to the internal list.
     * If the given material is already listed with a mass value, the new value is added
     * to the existing one.
     * @param tag The name of the material
     * @param ms The mass value
     */
  void MaterialProperties::addLocalMass(std::string tag, double ms) {
        msl_set = true;
        localmasses[tag] += ms;
    }

    /**
     * Add the local mass for a material, as specified by its tag, to the internal list.
     * Also keeps track of the originating component
     * If the given material is already listed with a mass value, the new value is added
     * to the existing one.
     * @param tag The name of the material
     * @param comp The name of the component
     * @param ms The mass value
     */
  void MaterialProperties::addLocalMass(std::string tag, std::string comp, double ms, int minZ) {
        msl_set = true;
        localmasses[tag] += ms;
        localmassesComp[getSubName(comp)] += ms;
        localCompMats[comp][tag] += ms; 
    }
    
    /**
     * Get the number of registered local masses for the materials found in the inactive element.
     * @return The size of the internal mass vector
     */
    unsigned int MaterialProperties::localMassCount() { return localmasses.size(); }
    
    /**
     * Get the number of registered local masses for the components found in the inactive element.
     * @return The size of the internal mass vector
     */
    unsigned int MaterialProperties::localMassCompCount() { return localmassesComp.size(); }
    
    /**
     * Reset the state of the internal mass vector to empty, discarding all entries.
     */
    void MaterialProperties::clearMassVectors() {
        localmasses.clear();
        localmassesComp.clear();
        localCompMats.clear();
    }
    
    /**
     * Copy the entire mass vector to another instance of <i>MaterialProperties</i>
     * @param mp The destination object
     */
    void MaterialProperties::copyMassVectors(MaterialProperties& mp) {
      mp.clearMassVectors(); //TODO: why?!?!?!?!?!
        //for (unsigned int i = 0; i < localMassCount(); i++) mp.addLocalMass(localmasses.at(i));
        //for (unsigned int i = 0; i < localMassCompCount(); i++) mp.addLocalMassComp(localmassesComp.at(i));
        for (std::map<std::string, std::map<std::string, double> >::iterator compit = localCompMats.begin(); compit != localCompMats.end(); ++compit)
            for (std::map<std::string, double>::iterator matit = compit->second.begin(); matit != compit->second.end(); ++matit)
                mp.addLocalMass(matit->first, compit->first, matit->second);
    }
    
    /**
     * Get the cumulative mass of the inactive element.
     * @return The overall mass, taking into account all registered materials; -1 if the value has not yet been computed
     */
    double MaterialProperties::getTotalMass() const { return total_mass; }
    
    /**
     * Get the local mass of the inactive element.
     * @return The overall mass, taking into account all registered materials; -1 if the value has not yet been computed
     */
    double MaterialProperties::getLocalMass() { return local_mass; }
    
    /**
     * Get the radiation length of the inactive element.
     * @return The overall radiation length, taking into account all registered materials; -1 if the value has not yet been computed
     */
    double MaterialProperties::getRadiationLength() { return r_length; }
    

    const std::map<std::string, RILength>& MaterialProperties::getComponentsRI() const { return componentsRI; } // CUIDADO: I know it parts with the old API but it's so much more practical this way

    /**
     * Get the intraction length of the inactive element.
     * @return The overall radiation length, taking into account all registered materials; -1 if the value has not yet been computed
     */
    double MaterialProperties::getInteractionLength() { return i_length; }

    RILength MaterialProperties::getMaterialLengths() {
      RILength myMaterials;
      myMaterials.interaction = i_length;
      myMaterials.radiation = r_length;
      return myMaterials;
    }
    
    /**
     * Calculate the overall mass of the inactive element from the internal mass vectors,
     * if the one of the mass vectors has at least one element in it, and an offset value.
     * @param offset A starting value for the calculation
     */
    void MaterialProperties::calculateTotalMass(double offset) {
        calculateLocalMass(offset);
        if (msl_set) total_mass = local_mass + offset;
        else total_mass = offset;
    }
    
    /**
     * Calculate the overall local mass of the inactive element from the internal mass vector,
     * if that mass vectors has at least one element in it, and an offset value.
     * @param offset A starting value for the calculation
     */
    void MaterialProperties::calculateLocalMass(double offset) {
        if (msl_set) {
            local_mass = offset;
            for (std::map<std::string, double>::iterator it = localmasses.begin(); it != localmasses.end(); ++it) {
                local_mass += it->second;
            }
        }
    }
    
    /**
     * Calculate the overall radiation length of the inactive element from the material table and the material vectors,
     * if they have at least one element in them, and an offset value.
     * @param materials A reference to the external material table
     * @param offset A starting value for the calculation
     */
    void MaterialProperties::calculateRadiationLength(MaterialTable& materials, double offset) {
        if (getSurface() > 0) {
            r_length = offset;
            if (msl_set) {
                // local mass loop
                for (std::map<std::string, double>::iterator it = localmasses.begin(); it != localmasses.end(); ++it) {
                    r_length += it->second / (materials.getMaterial(it->first).rlength * getSurface() / 100.0);
                }
                for (std::map<std::string, std::map<std::string, double> >::iterator cit = localCompMats.begin(); cit != localCompMats.end(); ++cit) {
                    for (std::map<std::string, double>::iterator mit = cit->second.begin(); mit != cit->second.end(); ++mit) {
                        componentsRI[getSuperName(cit->first)].radiation += mit->second / (materials.getMaterial(mit->first).rlength * getSurface() / 100.0);
                    }
                }
            }
        }
    }
    
    /**
     * Calculate the overall interaction length of the inactive element from the material table and the material vectors,
     * if they have at least one element in them, and an offset value.
     * @param materials A reference to the external material table
     * @param offset A starting value for the calculation
     */
    void MaterialProperties::calculateInteractionLength(MaterialTable& materials, double offset) {
        if (getSurface() > 0) {
            i_length = offset;
            if (msl_set) {
                // local mass loop
                for (std::map<std::string, double>::iterator it = localmasses.begin(); it != localmasses.end(); ++it) {
                    i_length += it->second / (materials.getMaterial(it->first).ilength * getSurface() / 100.0);
                }
                    
                for (std::map<std::string, std::map<std::string, double> >::iterator cit = localCompMats.begin(); cit != localCompMats.end(); ++cit) {
                    for (std::map<std::string, double>::iterator mit = cit->second.begin(); mit != cit->second.end(); ++mit) {
                        componentsRI[getSuperName(cit->first)].interaction += mit->second / (materials.getMaterial(mit->first).ilength * getSurface() / 100.0);
                    }
                }
            }
        }
    }

  // Versions with new material tab definition
    void MaterialProperties::calculateRadiationLength(double offset) {
      const material::MaterialTab& materialTab = material::MaterialTab::instance();
 
        if (getSurface() > 0) {
            r_length = offset;
            if (msl_set) {
                // local mass loop
                for (std::map<std::string, double>::iterator it = localmasses.begin(); it != localmasses.end(); ++it) {
                    r_length += it->second / (materialTab.radiationLength(it->first) * getSurface() / 100.0);
                }
                for (std::map<std::string, std::map<std::string, double> >::iterator cit = localCompMats.begin(); cit != localCompMats.end(); ++cit) {
                    for (std::map<std::string, double>::iterator mit = cit->second.begin(); mit != cit->second.end(); ++mit) {
                        componentsRI[getSuperName(cit->first)].radiation += mit->second / (materialTab.radiationLength(mit->first) * getSurface() / 100.0);
                    }
                }
            }
        }
    }
    
    void MaterialProperties::calculateInteractionLength(double offset) {
      const material::MaterialTab& materialTab =  material::MaterialTab::instance();

        if (getSurface() > 0) {
            i_length = offset;
            if (msl_set) {
                // local mass loop
                for (std::map<std::string, double>::iterator it = localmasses.begin(); it != localmasses.end(); ++it) {
                    i_length += it->second / (materialTab.interactionLength(it->first) * getSurface() / 100.0);
                }
                    
                for (std::map<std::string, std::map<std::string, double> >::iterator cit = localCompMats.begin(); cit != localCompMats.end(); ++cit) {
                    for (std::map<std::string, double>::iterator mit = cit->second.begin(); mit != cit->second.end(); ++mit) {
                        componentsRI[getSuperName(cit->first)].interaction += mit->second / (materialTab.interactionLength(mit->first) * getSurface() / 100.0);
                    }
                }
            }
        }
    }

    
    /**
     * Find out if the volume is relevant for tracking during analysis.
     * @return True if the material properties of this volume matter for the tracker analysis, false otherwise
     */
    bool MaterialProperties::track() { return trck; }
    
    /**
     * Set this volume's relevance for tracking during analysis.
     * @param tracking_on The flag that turns tracking of material properties for the volume on or off
     */
    void MaterialProperties::track(bool tracking_on) { trck = tracking_on; }
    
    /**
     * Print the material properties.
     */
    void MaterialProperties::print() {
        std::cout << "Material properties (current state)" << std::endl;
        std::cout << "localmasses: vector with " << localmasses.size() << " elements." << std::endl;
        int i = 0;
        for (std::map<std::string, double>::const_iterator it = localmasses.begin(); it != localmasses.end(); ++it)
            std::cout << "Material " << i++ << " (material, mass): (" << it->first << ", " << it->second << ")" << std::endl;
        i = 0;

        std::cout << "total_mass = " << total_mass << std::endl;
        std::cout << "local_mass = " << local_mass << std::endl;
        std::cout << "r_length = " << r_length << std::endl;
        std::cout << "i_length = " <<i_length  << std::endl;
    }
    


    /*-----protected-----*/

    std::string MaterialProperties::getSuperName(std::string name) const {
        std::stringstream ss(name);
        std::pair<std::string, std::string> split;
        std::getline(ss, split.first, '_');
        std::getline(ss, split.second, '_');
        return !split.second.empty() ? split.second : split.first;
    }

    std::string MaterialProperties::getSubName(std::string name) const {
        std::stringstream ss(name);
        std::pair<std::string, std::string> split;
        std::getline(ss, split.first, '_');
        std::getline(ss, split.second, '_');
        return split.first;
    }

define_enum_strings(MaterialProperties::Category) = { "Nocat", "Bmod", "Emod", "Bser", "Eser", "Bsup", "Esup", "Osup", "Tsup", "Usup" };
}

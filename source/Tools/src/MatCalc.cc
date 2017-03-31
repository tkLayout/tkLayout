/**
 * @file MatCalc.cc
 * @brief This class implements the algorithms that allocate materials to the various types of structures in the tracker
 */

#include <MatCalc.h>

    // public
    /**
     * This function returns the initialisation status.
     * @return True if the internal data structures have been initialised from a config file, false otherwise
     */
    bool MatCalc::initDone() { return init_done; }
    
    /**
     * This function sets the initialisation status.
     * @param yes True if the internal data structures contain information, false otherwise
     */
    void MatCalc::initDone(bool yes) { init_done = yes; }
    
    /**
     * This function resets the internal data structures.
     */
    void MatCalc::reset() {
        internals.modinfo.clear();
        internals.serlocalinfo.clear();
        internals.serexitinfo.clear();
        internals.supinfo.clear();
        init_done = false;
    }
    
    /**
     * Getter for the default number of strips of a given type. If the requested type is unknown
     * or not on the list, an exception is thrown.
     * @param type The module type for the requested value
     * @return The default number of strips across a module
     */
    int MatCalc::getStripsAcross(std::string type) {
        try {
            TypeInfo& info = getTypeInfoByType(type);
            return info.strips_across;
        }
        catch(std::range_error& re) {
            std::cerr << "MatCalc::getStripsAcross(): " << re.what() << std::endl;
            throw;
        }
    }
    
    /**
     * Getter for the default number of segments of a given type. If the requested type is unknown
     * or not on the list, an exception is thrown.
     * @param type The module type for the requested value
     * @return The default number of segments along a module
     */
    int MatCalc::getSegmentsAlong(std::string type) {
        try {
            TypeInfo& info = getTypeInfoByType(type);
            return info.segments_along;
        }
        catch(std::range_error& re) {
            std::cerr << "MatCalc::getSegmentsAlong(): " << re.what() << std::endl;
            throw;
        }
    }
    
    /**
     * Compact getter for the default number of strips and segments of a given type. If the
     * requested type is unknown or not on the list, an exception is thrown.
     * @param type The module type for the requested value
     * @return The default number of strips across and segments along a module, bundled into a <i>std::pair</i>
     */
    std::pair<int, int> MatCalc::getDefaultDimensions(std::string type) {
        try {
            std::pair<int, int> result;
            TypeInfo& info = getTypeInfoByType(type);
            result.first = info.strips_across;
            result.second = info.segments_along;
            return result;
        }
        catch(std::range_error& re) {
            std::cerr << "MatCalc::getDefaultDimensions(): " << re.what() << std::endl;
            throw;
        }
    }
    
    /**
     * Adding a module type and its default dimensions is done in here.
     * @param type The type of the module that needs to be added
     * @param strips The default number of strips across the module
     * @param segments The default number of segments along the module
     */
    void MatCalc::addTypeInfo(std::string type, int strips, int segments) {
        if (!entryExists(type)) {
            TypeInfo info;
            info.type = type;
            info.strips_across = strips;
            info.segments_along = segments;
            std::vector<SingleMod> tmp;
            std::pair<TypeInfo, std::vector<SingleMod> > tpair(info, tmp);
            internals.modinfo.push_back(tpair);
        }
    }
    
    /**
     * Updating the default number of strips on a given module type is done here.
     * @param type The type of the module that needs to be changed
     * @param strips The new default number of strips across a module
     */
    void MatCalc::updateTypeInfoStrips(std::string type, int strips) {
        if (entryExists(type)) {
            try {
                TypeInfo& info = getTypeInfoByType(type);
                info.strips_across = strips;
            }
            catch (std::range_error& re) {
                std::cerr << "MatCalc::updateTypeInfoStrips(): " << err_up_general << std::endl;
            }
        }
    }
    
    /**
     * Updating the default number of segments on a given module type is done here.
     * @param type The type of the module that needs to be changed
     * @param segments The new default number of segments along a module
     */
    void MatCalc::updateTypeInfoSegments(std::string type, int segments) {
        if (entryExists(type)) {
            try {
                TypeInfo& info = getTypeInfoByType(type);
                info.segments_along = segments;
            }
            catch (std::range_error& re) {
                std::cerr << "MatCalc::updateTypeInfoSegments(): " << err_up_general << std::endl;
            }
        }
    }
    
    /**
     * Add a parameter entry to the module info vector corresponding to <i>type</i> by copying the supplied
     * values. Add the new values to the old ones if a material with the given parameters already appears on
     * the list.
     * @param tag A unique material identifier
     * @param type The type of module the material parameters refer to
     * @param A The coefficient for the layout-dependent and travelling portion
     * @param uA The unit of parameter <i>A</i>
     * @param B The coefficient for the layout-dependent and localised portion
     *@param uB The unit of parameter <i>B</i>
     * @param C The coefficient for the per-module and travelling portion
     *@param uC The unit of parameter <i>C</i>
     * @param D The coefficient for the per-module and localised portion
     *@param uD The unit of parameter <i>D</i>
     * @param local A flag indicating whether the amounts as a whole will be treated as local or as exiting the module later
     */
  void MatCalc::addModuleParameters(std::string tag, std::string type, std::string comp,
            double A, Matunit uA, double B, Matunit uB, double C, Matunit uC, double D, Matunit uD, bool local) {
        try {
	  if (!entryExists(tag, type, comp, uA, uB, uC, uD, local)) {
                // Create new entry if there is none
                SingleMod mod;
                mod.tag = tag;
		mod.comp = comp;
                mod.A = A;
                mod.uA = uA;
                mod.B = B;
                mod.uB = uB;
                mod.C = C;
                mod.uC = uC;
                mod.D = D;
                mod.uD = uD;
                mod.is_local = local;
                std::vector<SingleMod>& vect = getModVector(type);
                vect.push_back(mod);
            }
            else {
                // Add new amounts to existing entry otherwise
	      SingleMod& mod = getSingleMod(tag, type, comp, uA, uB, uC, uD, local);
                mod.A = mod.A + A;
                mod.B = mod.B + B;
                mod.C = mod.C + C;
                mod.D = mod.D + D;
            }
        }
        catch(std::range_error& re) {
            std::cerr << "MatCalc::addModuleParameters(): " << err_matadd_weird << std::endl;
            std::cerr << re.what() << std::endl;
        }
    }
    
    /**
     * Add a parameter entry to the local service info vector.
     * @param tag A unique material identifier
     * @param Q The amount of material to be added to the service
     * @param uQ The unit of parameter <i>Q</i>
     */
    void MatCalc::addServiceParameters(std::string tag, double Q, Matunit uQ) {
        if (!entryExists(tag, uQ)) {
            // Create new entry if there is none
            SingleSerLocal serl;
            serl.tag = tag;
            serl.Q = Q;
            serl.uQ = uQ;
            internals.serlocalinfo.push_back(serl);
        }
        else {
            // Add new amount to existing entry otherwise
            SingleSerLocal& serl = getSingleSer(tag, uQ);
            serl.Q = serl.Q + Q;
        }
    }
    
    /**
     * Add a parameter entry to the info vector mapping transitions from modules to services.
     * @param tagIn A unique identifier for the module material
     * @param In The amount of module material that is to be converted
     * @param uIn The unit of parameter <i>In</i>
     * @param tagOut A unique identifier for the service material
     * @param Out The amount of service material that comes out of the conversion
     * @param uOut The unit of parameter <i>Out</i>
     * @param local A flag indicating whether the resulting amount of service material will be treated as local or exiting
     */
    void MatCalc::addServiceParameters(std::string tagIn, double In, Matunit uIn,
            std::string tagOut, double Out, Matunit uOut, bool local) {
        if (!entryExists(tagIn, tagOut, uIn, uOut, local)) {
            // Create new entry if there is none
            SingleSerExit ser;
            ser.tagIn = tagIn;
            ser.In = In;
            ser.uIn = uIn;
            ser.tagOut = tagOut;
            ser.Out = Out;
            ser.uOut = uOut;
            ser.is_local = local;
            internals.serexitinfo.push_back(ser);
        }
        else {
            // Add new amount to existing entry otherwise
            SingleSerExit& ser = getSingleSer(tagIn, tagOut, uIn, uOut, local);
            ser.Out = ser.Out + ser.In / In * Out;
        }
    }
    
    /**
     * Add a parameter entry to the supports info vector.
     * @param tag A unique material identifier
     * @param M The amount of material to be added to the support structure
     * @param uM The unit of parameter <i>M</i>
     * @param cM The category of supports that this amount of material belongs to
     */
    void MatCalc::addSupportParameters(std::string tag, double M, Matunit uM, MaterialProperties::Category cM) {
        if (!entryExists(tag, uM, cM)) {
            // Create new entry if there is none
            SingleSup sup;
            sup.tag = tag;
            sup.M = M;
            sup.uM = uM;
            sup.cM = cM;
            internals.supinfo.push_back(sup);
        }
        else {
            // Add new amount to existing entry otherwise
            SingleSup& sup = getSingleSup(tag, uM, cM);
            sup.M = sup.M + M;
        }
    }
    
    
    /**
     * Reset the module information vector
     */
    void MatCalc::clearModVectors() { internals.modinfo.clear(); }
    
    /**
     * Reset a module info vector of a given type.
     * @param type The type indentifier of the vector that needs to be cleared
     */
    void MatCalc::clearModVector(std::string type) {
        try {
            std::vector<SingleMod>& vect = getModVector(type);
            vect.clear();
        }
        catch(std::range_error& re) {
            std::cerr << "MatCalc::clearModVector(): " << re.what() << std::endl;
        }
    }
    
    /**
     * Copy the contents of one module info vector into another.
     * @param source The type of the source vector
     * @param dest The type of the destination vector
     */
    void MatCalc::copyContents(std::string source, std::string dest) {
        clearModVector(dest);
        appendContents(source, dest);
    }
    
    /**
     * Append the contents of one module info vector to another.
     * @param source The type of the source vector
     * @param dest The type of the destination vector
     */
    void MatCalc::appendContents(std::string source, std::string dest) {
        if (source.compare(dest) != 0) {
            try {
                std::vector<SingleMod>& in = getModVector(source);
                std::vector<SingleMod>& out = getModVector(dest);
                std::copy(in.begin(), in.end(), std::back_inserter(out));
            }
            catch(std::range_error& re) {
                std::cerr << "MatCalc::appendContents(): " << re.what() << std::endl;
            }
        }
    }
    
    /**
     * Check if there is internal information about a a given module type or not.
     * @param type The identifier used to query the collection of type info entries
     * @return True if the given module type appears in the <i>typeinfo</i> vector, false otherwise
     */
    bool MatCalc::typeRegistered(std::string type) {
        std::vector<std::pair<TypeInfo, std::vector<SingleMod> > >::const_iterator iter = internals.modinfo.begin();
        std::vector<std::pair<TypeInfo, std::vector<SingleMod> > >::const_iterator guard = internals.modinfo.end();
        while (iter != guard) {
            if (type.compare(iter->first.type) == 0) return true;
            else iter++;
        }
        return false;
    }
    
    /**
     * This is a convenince function that returns the number of registered module types.
     * @return The number of module types with registered default dimensions
     */
    unsigned int MatCalc::registeredTypes() { return internals.modinfo.size(); }
    
    /**
     * Getter for the internal global material table. The material table always exists but
     * it may be empty if the class has not been initialised.
     * @return A reference to the global material table
     */
    MaterialTable& MatCalc::getMaterialTable() { return mt; }
    
    /**
     * This is the core function that assigns materials to barrel modules once the calculator has been initialised.
     * It loops through the elements of the provided vector of vectors and calculates the relevant material mix for
     * the element based on its position along a rod. If a material is declared with a unit other than grammes, the
     * resulting mass is calculated from the value and unit of the material combined with the geometry and size
     * of the element before it is added to the total.
     * @param barrelcaps The collection mapping to the barrel modules that need to have a material mix assigned to them
     * @return True if there were no errors during processing, false otherwise
     */
    bool MatCalc::calculateBarrelMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps) { // sorry, but this code is a POS
      // layer loop
      for (unsigned int i = 0; i < barrelcaps.size(); i++) {
        if (barrelcaps.at(i).size() > 0) {
          try {
            int rindex = 0;
            int maxRing = 0; // only considers Z+
            std::vector<double> stripseg_scalars;
            std::vector<std::string> mtypes;
            std::vector<std::list<int> > modinrings;
            // module loop for ring types and multipliers for strips and segments
            for (unsigned int j = 0; j < barrelcaps.at(i).size(); j++) {
              // ring index of current module
              rindex = barrelcaps.at(i).at(j).getModule().uniRef().ring;
              maxRing = barrelcaps.at(i).at(j).getModule().uniRef().side > 0 ? MAX(barrelcaps.at(i).at(j).getModule().uniRef().ring, maxRing) : maxRing;
              // collect ring types
              if ((int)mtypes.size() < rindex) {
                while ((int)mtypes.size() < rindex) mtypes.push_back("");
              }
              if (mtypes.at(rindex - 1).empty()) mtypes.at(rindex - 1) = barrelcaps.at(i).at(j).getModule().moduleType();
              // collect multipliers for strips and segments
              if ((int)stripseg_scalars.size() < rindex) {
                while ((int)stripseg_scalars.size() < rindex) stripseg_scalars.push_back(0.0);
              }
              if (stripseg_scalars.at(rindex - 1) == 0.0) {
                if (!mtypes.at(rindex - 1).empty()) {
                  stripseg_scalars.at(rindex - 1) = (double)barrelcaps.at(i).at(j).getModule().outerSensor().numStripsAcross();
                  stripseg_scalars.at(rindex - 1) = stripseg_scalars.at(rindex - 1) / (double)getStripsAcross(mtypes.at(rindex - 1));
                  stripseg_scalars.at(rindex - 1) = stripseg_scalars.at(rindex - 1) * ((double)barrelcaps.at(i).at(j).getModule().outerSensor().numSegments()); // CUIDADO: as of now only the outer sensor is considered for the scaling
                  stripseg_scalars.at(rindex - 1) = stripseg_scalars.at(rindex - 1) / (double)getSegmentsAlong(mtypes.at(rindex - 1));
                }
              }
              // record module indices per ring
              if ((int)modinrings.size() < rindex) {
                while ((int)modinrings.size() < rindex) {
                  std::list<int> tmp;
                  modinrings.push_back(tmp);
                }
              }
              modinrings.at(rindex - 1).push_back(j);
            }
            rindex = modinrings.size(); // CUIDADO Tentative fix, rindex is max of ring index found in modules

            // ring loop
            for (int j = 0; j < /*maxRing*/ rindex; j++) { // CUIDADO rindex WTF!?!?!?
              if (!modinrings.at(j).empty()) {
                double A, B, C, D;
                double density, surface, length;
                std::list<int>::iterator first = modinrings.at(j).begin();
                std::list<int>::iterator start = modinrings.at(j).begin(); start++;
                std::list<int>::iterator stop = modinrings.at(j).end();
                surface = barrelcaps.at(i).at(*first).getSurface();
                if (surface < 0) {
                  std::cerr << msg_negative_area << " Barrel module in layer " << i << ", position " << j;
                  std::cerr << " with index " << *first << " within the layer. " << msg_abort << std::endl;
                  return false;
                }
                else {
                  length = barrelcaps.at(i).at(*first).getModule().length();

#define TILTED_HOTFIX
#ifdef TILTED_HOTFIX
                  // hot fix to scale the services materials for the non contiguous modules in the tilted barrel
                  //if (barrelcaps.at(i).at(modinrings.at(j).front()).getModule().tiltAngle() != 0.) {  // CUIDADO Hot fixing the hot fix -- now applied to untilted modules too!
                    if (j == 0) { // the following code crashes for a layer with a single module, but c'mon, when is that gonna happen?
                      double zMin = barrelcaps.at(i).at(modinrings.at(j).front()).getModule().center().Z();
                      double zThis = barrelcaps.at(i).at(modinrings.at(j).front()).getModule().minZ();
                      double zNext = barrelcaps.at(i).at(modinrings.at(j+1).front()).getModule().center().Z();
                      length = (zThis - zMin) + (zNext - zThis)/2.;
                    } else if (j < modinrings.size()-1 &&
                        (barrelcaps.at(i).at(modinrings.at(j).front()).getModule().side() == // CUIDADO we check if we switch sides (due to asym barrel), 
                         barrelcaps.at(i).at(modinrings.at(j+1).front()).getModule().side())) { // jumping to the farthest z- mod (nasty bug in orig code)
                      double zPrev = barrelcaps.at(i).at(modinrings.at(j-1).front()).getModule().center().Z();
                      double zNext = barrelcaps.at(i).at(modinrings.at(j+1).front()).getModule().center().Z();
                      length = (zNext - zPrev) / 2; // override length with half the distance between modules, so as to scale material defined as grams per metres over a bigger length
                    } else if (barrelcaps.at(i).at(modinrings.at(j).front()).getModule().side() == // CUIDADO we check if we've switched sides
                               barrelcaps.at(i).at(modinrings.at(j-1).front()).getModule().side()) { // jumping to the farthest z- mod (nasty bug in orig code)
                      double zPrev = barrelcaps.at(i).at(modinrings.at(j-1).front()).getModule().center().Z();
                      double zThis = barrelcaps.at(i).at(modinrings.at(j).front()).getModule().center().Z();
                      double zMax = barrelcaps.at(i).at(modinrings.at(j).front()).getModule().maxZ();
                      length = (zMax - zThis) + (zThis - zPrev)/2.;
                    }
                  //}
#endif
                  std::vector<SingleMod>& vect = getModVector(mtypes.at(j));
                  std::vector<SingleMod>::const_iterator guard = vect.end();
                  std::vector<SingleMod>::const_iterator iter;
                  // materials loop
                  for (iter = vect.begin(); iter != guard; iter++) {
                    // measurement for unit conversion
                    density = mt.getMaterial(iter->tag).density;
                    // unit conversion per parameter (internal unit is grammes)
                    if (iter->uB == grpm) B = convert(iter->B, iter->uB, length);
                    else B = convert(iter->B, iter->uB, density, surface);
                    if (iter->uD == grpm) D = convert(iter->D, iter->uD, length);
                    else D = convert(iter->D, iter->uD, density, surface);
                    // parameter scaling
                    B = B * stripseg_scalars.at(j);
                    // save converted and scaled material
                    if (iter->is_local) barrelcaps.at(i).at(*first).addLocalMass(iter->tag, iter->comp, B + D); // CUIDADO MIXED UP j and *first?? was j
                    //else barrelcaps.at(i).at(*first).addExitingMass(iter->tag, iter->comp, B + D); // was j
                  }
                  //accumulation of travelling parameters for outer rings
                  if (j > 0) {
                    // inner rings accumulation loop
                    for (int k = 0; k < j; k++) {
                      surface = barrelcaps.at(i).at(modinrings.at(k).front()).getSurface();
                      if (surface < 0) {
                        std::cerr << msg_negative_area << " Barrel module in layer " << i << ", position " << k;
                        std::cerr << " with index " << modinrings.at(k).front() << " within the layer. " << msg_abort << std::endl;
                        return false;
                      }
                      else {
                        length = barrelcaps.at(i).at(modinrings.at(k).front()).getModule().length();
#ifdef TILTED_HOTFIX
                        // hot fix to scale the services materials for the non contiguous modules in the tilted barrel
                        // if (barrelcaps.at(i).at(modinrings.at(j).front()).getModule().tiltAngle() != 0.) {
                          if (j < modinrings.size()-1 && // we never accumulate material for the first module, so we don't need to take care of the case like in the firt block of hotfix
                              (barrelcaps.at(i).at(modinrings.at(j).front()).getModule().side() == // CUIDADO we check if we switch sides (due to asym barrel), 
                               barrelcaps.at(i).at(modinrings.at(j+1).front()).getModule().side())) { // jumping to the farthest z- mod (nasty bug in orig code)
                            double zPrev = barrelcaps.at(i).at(modinrings.at(j-1).front()).getModule().center().Z();
                            double zNext = barrelcaps.at(i).at(modinrings.at(j+1).front()).getModule().center().Z();
                            length = (zNext - zPrev) / 2; // override length with half the distance between modules, so as to scale material defined as grams per metres over a bigger length
                          } else if (barrelcaps.at(i).at(modinrings.at(j).front()).getModule().side() == // CUIDADO we check if we've switched sides
                                     barrelcaps.at(i).at(modinrings.at(j-1).front()).getModule().side()) { // jumping to the farthest z- mod (nasty bug in orig code)
                            double zPrev = barrelcaps.at(i).at(modinrings.at(j-1).front()).getModule().center().Z();
                            double zThis = barrelcaps.at(i).at(modinrings.at(j).front()).getModule().center().Z();
                            double zMax = barrelcaps.at(i).at(modinrings.at(j).front()).getModule().maxZ();
                            length = (zMax - zThis) + (zThis - zPrev)/2.;
                          }
                        // }
#endif
                        std::vector<SingleMod>& vect = getModVector(mtypes.at(k));
                        std::vector<SingleMod>::const_iterator iter, guard = vect.end();
                        for (iter = vect.begin(); iter != guard; iter++) {
                          // measurement for unit conversion
                          density = mt.getMaterial(iter->tag).density;
                          // unit conversion per parameter (internal unit is grammes)
                          if (iter->uA == grpm) A = convert(iter->A, iter->uA, length);
                          else A = convert(iter->A, iter->uA, density, surface);
                          if (iter->uC == grpm) C = convert(iter->C, iter->uC, length);
                          else C = convert(iter->C, iter->uC, density, surface);
                          // parameter scaling
                          A = A * stripseg_scalars.at(k);
                          // save converted and scaled material
                          if (iter->is_local) barrelcaps.at(i).at(*first).addLocalMass(iter->tag, iter->comp, A + C); // was j
                          //else barrelcaps.at(i).at(*first).addExitingMass(iter->tag, iter->comp, A + C); // was j
                        }
                      }
                    }
                  }
                  // use recorded materials to calculate material properties for every module in current ring
                  barrelcaps.at(i).at(*first).calculateTotalMass();
                  barrelcaps.at(i).at(*first).calculateRadiationLength(mt);
                  barrelcaps.at(i).at(*first).calculateInteractionLength(mt);
                  while (start != stop) {
                    barrelcaps.at(i).at(*first).copyMassVectors(barrelcaps.at(i).at(*start));
                    barrelcaps.at(i).at(*start).calculateTotalMass();
                    barrelcaps.at(i).at(*start).calculateRadiationLength(mt);
                    barrelcaps.at(i).at(*start).calculateInteractionLength(mt);
                    start++;
                  }
                }
              }
            }
          }
          catch(std::range_error& re) {
            std::cerr << "MatCalc::calculateBarrelMaterials(): " << re.what() << " " << msg_abort << std::endl;
            return false;
          }
        }
      }
      return true;
    }

    /**
     * This is the core function that assigns materials to endcap modules once the calculator has been initialised.
     * It loops through the elements of the provided vector of vectors and calculates the relevant material mix for
     * the element based on the ring and disc it belongs to. If a material is declared with a unit other than grammes, the
     * resulting mass is calculated from the value and unit of the material combined with the geometry and size
     * of the element before it is added to the total.
     * @param endcapcaps The collection mapping to the endcap modules that need to have a material mix assigned to them
     * @return True if there were no errors during processing, false otherwise
     */
    bool MatCalc::calculateEndcapMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps) {
      // disc loop
      for (unsigned int i = 0; i < endcapcaps.size(); i++) {
        if (endcapcaps.at(i).size() > 0) {
          try {
            int rindex;
            std::vector<int> mods;
            std::vector<double> stripseg_scalars;
            std::vector<std::string> mtypes;
            std::vector<std::list<int> > modinrings;
            // module loop for ring types and multipliers for strips and segments
            for (unsigned int j = 0; j < endcapcaps.at(i).size(); j++) {
              // sum up the number of modules per ring
              rindex = endcapcaps.at(i).at(j).getModule().uniRef().ring;
              if ((int)mods.size() < rindex) {
                while ((int)mods.size() < rindex) mods.push_back(0);
              }
              mods.at(rindex - 1) = mods.at(rindex - 1) + 1;
              // collect ring types
              if ((int)mtypes.size() < rindex) {
                while ((int)mtypes.size() < rindex) mtypes.push_back("");
              }
              if (mtypes.at(rindex - 1).empty()) mtypes.at(rindex - 1) = endcapcaps.at(i).at(j).getModule().moduleType();
              // collect multipliers for strips and segments
              if ((int)stripseg_scalars.size() < rindex) {
                while ((int)stripseg_scalars.size() < rindex) stripseg_scalars.push_back(0.0);
              }
              if (stripseg_scalars.at(rindex - 1) == 0.0) {
                if (!mtypes.at(rindex - 1).empty()) {
                  stripseg_scalars.at(rindex - 1) = (double)endcapcaps.at(i).at(j).getModule().outerSensor().numStripsAcross(); // CUIDADO as it this now, only the outer sensor is considered for the scaling
                  stripseg_scalars.at(rindex - 1) = stripseg_scalars.at(rindex - 1) / (double)getStripsAcross(mtypes.at(rindex - 1));
                  stripseg_scalars.at(rindex - 1) = stripseg_scalars.at(rindex - 1) * ((double)endcapcaps.at(i).at(j).getModule().outerSensor().numSegments());
                  stripseg_scalars.at(rindex - 1) = stripseg_scalars.at(rindex - 1) / (double)getSegmentsAlong(mtypes.at(rindex - 1));
                }
              }
              // record module indices per ring
              if ((int)modinrings.size() < rindex) {
                while ((int)modinrings.size() < rindex) {
                  std::list<int> tmp;
                  modinrings.push_back(tmp);
                }
              }
              modinrings.at(rindex - 1).push_back(j);
            }
            rindex = modinrings.size();
            // ring loop
            for (int j = 0; j < rindex; j++) {
              if (!modinrings.at(j).empty()) {
                double A, B, C, D;
                double density, surface, length;
                std::list<int>::iterator first = modinrings.at(j).begin();
                std::list<int>::iterator start = modinrings.at(j).begin(); start++;
                std::list<int>::iterator stop = modinrings.at(j).end();
                surface = endcapcaps.at(i).at(*first).getSurface();
                if (surface < 0) {
                  std::cerr << msg_negative_area << " Endcap module in disc " << i << ", ring " << j;
                  std::cerr << " with index " << *first << " within the disc. " << msg_abort << std::endl;
                  return false;
                }
                // calculation of static parameters for all rings
                else {
                  length = endcapcaps.at(i).at(*first).getModule().length();
                  std::vector<SingleMod>& vect = getModVector(mtypes.at(j));
                  std::vector<SingleMod>::const_iterator iter, guard = vect.end();
                  // materials loop
                  for (iter = vect.begin(); iter != guard; iter++) {
                    // unit conversion per parameter (internal unit is grammes)
                    density = mt.getMaterial(iter->tag).density;
                    if (iter->uB == grpm) B = convert(iter->B, iter->uB, length);
                    else B = convert(iter->B, iter->uB, density, surface);
                    if (iter->uD == grpm) D = convert(iter->D, iter->uD, length);
                    else D = convert(iter->D, iter->uD, density, surface);
                    // parameter scaling
                    B = B * stripseg_scalars.at(j);
                    // save converted and scaled material
                    if (iter->is_local) endcapcaps.at(i).at(*first).addLocalMass(iter->tag, iter->comp, B + D);
                    //else endcapcaps.at(i).at(*first).addExitingMass(iter->tag, iter->comp, B + D);
                  }
                  // accumulation of travelling parameters for outer rings
                  if (j > 0) {
                    // inner rings accumulation loop
                    for (int k = 0; k < j; k++) {
                      if (mods.at(k) > 0) {
                        surface = endcapcaps.at(i).at(modinrings.at(k).front()).getSurface();
                        if (surface < 0) {
                          std::cerr << msg_negative_area << " Endcap module in disc " << i << ", ring " << k;
                          std::cerr << " with index " << modinrings.at(k).front() << " within the disc. " << msg_abort << std::endl;
                          return false;
                        }
                        else {
                          length = endcapcaps.at(i).at(modinrings.at(k).front()).getModule().length();
                          std::vector<SingleMod>& vect = getModVector(mtypes.at(k));
                          std::vector<SingleMod>::const_iterator guard = vect.end();
                          std::vector<SingleMod>::const_iterator iter;
                          // materials loop
                          for (iter = vect.begin(); iter != guard; iter++) {
                            // unit conversion per parameter (internal unit is grammes)
                            density = mt.getMaterial(iter->tag).density;
                            if (iter->uA == grpm) A = convert(iter->A, iter->uA, length);
                            else A = convert(iter->A, iter->uA, density, surface);
                            if (iter->uC == grpm) C = convert(iter->C, iter->uC, length);
                            else C = convert(iter->C, iter->uC, density, surface);
                            // parameter scaling
                            A = A * (double)mods.at(k) / (double)mods.at(j) * stripseg_scalars.at(k);
                            C = C * (double)mods.at(k) / (double)mods.at(j);
                            // save converted and scaled material
                            if (iter->is_local) endcapcaps.at(i).at(*first).addLocalMass(iter->tag, iter->comp, A + C);
                            //else endcapcaps.at(i).at(*first).addExitingMass(iter->tag, iter->comp, A + C);
                          }
                        }
                      }
                    }
                  }
                  // use recorded materials to calculate material properties for every module in current ring
                  endcapcaps.at(i).at(*first).calculateTotalMass();
                  endcapcaps.at(i).at(*first).calculateRadiationLength(mt);
                  endcapcaps.at(i).at(*first).calculateInteractionLength(mt);
                  while (start != stop) {
                    endcapcaps.at(i).at(*first).copyMassVectors(endcapcaps.at(i).at(*start));
                    endcapcaps.at(i).at(*start).calculateTotalMass();
                    endcapcaps.at(i).at(*start).calculateRadiationLength(mt);
                    endcapcaps.at(i).at(*start).calculateInteractionLength(mt);
                    start++;
                  }
                }
              }
            }
          }
          catch (std::range_error& re) {
            std::cerr << "MatCalc::calculateEndcapMaterials(): " << re.what() << " " << msg_abort << std::endl;
            return false;
          }
        }
      }
      return true;
    }

    /**
     * This is the core function that assigns materials to barrel services once the calculator has been initialised.
     * It loops through the elements of the provided vector and calculates the relevant material mix for the element
     * based on the feeder and neighbour volumes. If a material is declared with a unit other than grammes, the
     * resulting mass is calculated from the value and unit of the material combined with the geometry and size
     * of the element before it is added to the total.
     * @param barrelcaps The collection mapping to the barrel modules that may act as feeder volumes to the service
     * @param barrelservices The collection of barrel services that need to have a material mix assigned to them
     * @param endcapservices The collection of endcap service volumes - they may act as neighbour volumes to the barrel services
     * @return True if there were no errors during processing, false otherwise
     */
    bool MatCalc::calculateBarrelServiceMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps,
                                                  std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices) {
      int feeder, neighbour;
      InactiveElement::InType ftype, ntype;
      double length, surface;
      // main loop over barrel service volumes
      for(unsigned int i = 0; i < barrelservices.size(); i++) {
//        try {
          // collect information about feeder and neighbour volumes
          ftype = barrelservices.at(i).getFeederType();
          feeder = barrelservices.at(i).getFeederIndex();
          ntype = barrelservices.at(i).getNeighbourType();
          neighbour = barrelservices.at(i).getNeighbourIndex();
          // measurements for unit conversion
          if (barrelservices.at(i).isVertical()) length = barrelservices.at(i).getRWidth();
          else length = barrelservices.at(i).getZLength();
          surface = barrelservices.at(i).getSurface();
          // handle incoming materials from feeder volume
          switch(ftype) {
          case InactiveElement::no_in : break;
          case InactiveElement::tracker : adjacentDifferentCategory(barrelcaps.at(feeder), barrelservices.at(i),
                                                                    findBarrelRods(barrelcaps, feeder), length, surface);
                                          break;
          case InactiveElement::barrel : adjacentSameCategory(barrelservices.at(feeder), barrelservices.at(i));
                                         break;
          case InactiveElement::endcap : adjacentSameCategory(endcapservices.at(feeder), barrelservices.at(i));
          }
          // handle incoming materials from neighbour volume
          switch(ntype) {
          case InactiveElement::no_in : break;
          case InactiveElement::tracker : adjacentDifferentCategory(barrelcaps.at(neighbour), barrelservices.at(i),
                                                                    findBarrelRods(barrelcaps, neighbour), length, surface);
                                          break;
          case InactiveElement::barrel : adjacentSameCategory(barrelservices.at(neighbour), barrelservices.at(i));
                                         break;
          case InactiveElement::endcap : adjacentSameCategory(endcapservices.at(neighbour), barrelservices.at(i));
          }
          // use recorded materials to calculate material properties
          barrelservices.at(i).calculateTotalMass();
          barrelservices.at(i).calculateRadiationLength(mt);
          barrelservices.at(i).calculateInteractionLength(mt);
//        }
//        catch(std::runtime_error& re) {
//          std::cerr << "MatCalc::calculateBarrelServiceMaterials(): " << re.what() << " " << msg_abort << std::endl;
//          return false;
//        }
//        catch(std::exception& e) {
//          std::cerr << "MatCalc::calculateBarrelServiceMaterials(): " << e.what() << " " << msg_abort << std::endl;
//          return false;
//        }
      }
      return true;
    }

    /**
     * This is the core function that assigns materials to endcap services once the calculator has been initialised.
     * It loops through the elements of the provided vector and calculates the relevant material mix for the element
     * based on the feeder and neighbour volumes. If a material is declared with a unit other than grammes, the
     * resulting mass is calculated from the value and unit of the material combined with the geometry and size
     * of the element before it is added to the total.
     * @param endcapcaps The collection mapping to the endcap modules that may act as feeder volumes to the service
     * @param barrelservices The collection of barrel service volumes - they may act as neighbour volumes to the endcap services
     * @param endcapservices The collection of endcap services that need to have a material mix assigned to them
     * @return True if there were no errors during processing, false otherwise
     */
    bool MatCalc::calculateEndcapServiceMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps,
                                                  std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices) {
      int feeder, neighbour;
      InactiveElement::InType ftype, ntype;
      double length, surface;
      // main loop over endcap service volumes
      for (unsigned int i = 0; i < endcapservices.size(); i++) {
        //try {
          // collect information about feeder and neighbour volumes
          ftype = endcapservices.at(i).getFeederType();
          feeder = endcapservices.at(i).getFeederIndex();
          ntype = endcapservices.at(i).getNeighbourType();
          neighbour = endcapservices.at(i).getNeighbourIndex();
          // measurements for unit conversion
          if (endcapservices.at(i).isVertical()) length = endcapservices.at(i).getRWidth();
          else length = endcapservices.at(i).getZLength();
          surface = endcapservices.at(i).getSurface();
          // handle incoming materials from feeder volumes
          switch(ftype) {
          case InactiveElement::no_in : break;
          case InactiveElement::tracker : adjacentDifferentCategory(endcapcaps.at(feeder), endcapservices.at(i),
                                                                    findEndcapRods(endcapcaps, feeder), length, surface);
                                          break;
          case InactiveElement::barrel : adjacentSameCategory(barrelservices.at(feeder), endcapservices.at(i));
                                         break;
          case InactiveElement::endcap : adjacentSameCategory(endcapservices.at(feeder), endcapservices.at(i));
          }
          // handle incoming materials from neighbour volumes
          switch(ntype) {
          case InactiveElement::no_in : break;
          case InactiveElement::tracker : adjacentDifferentCategory(endcapcaps.at(neighbour), endcapservices.at(i),
                                                                    findEndcapRods(endcapcaps, neighbour), length, surface);
                                          break;
          case InactiveElement::barrel : adjacentSameCategory(barrelservices.at(neighbour), endcapservices.at(i));
                                         break;
          case InactiveElement::endcap : adjacentSameCategory(endcapservices.at(neighbour), endcapservices.at(i));
          }
          // use recorded materials to calculate material properties
          endcapservices.at(i).calculateTotalMass();
          endcapservices.at(i).calculateRadiationLength(mt);
          endcapservices.at(i).calculateInteractionLength(mt);
       // }
       // catch(std::runtime_error& re) {
       //   std::cerr << "MatCalc::calculateEndcapServiceMaterials(): " << re.what() << " " << msg_abort << std::endl;
       //   return false;
       // }
       // catch(std::exception& e) {
       //   std::cerr << "MatCalc::calculateEndcapServiceMaterials(): " << e.what() << " " << msg_abort << std::endl;
       //   return false;
       // }
      }
      return true;
    }

    /**
     * This is the core function that assigns materials to support structures once the calculator has been initialised.
     * It loops through the elements of the provided vector, checks the category of the support and calculates the
     * relevant material mix for the element. If a material is declared with a unit other than grammes, the resulting
     * mass is calculated from the value and unit of the material combined with the geometry and size of the element
     * before it is added to the total.
     * @param supports The vector of support parts that need to have a material mix assigned to them
     * @return True if there were no errors during processing, false otherwise
     */
    bool MatCalc::calculateSupportMaterials(std::vector<InactiveElement>& supports) {
      double length, surface;
      try {
        // main loop over support volumes
        for (unsigned int i = 0; i < supports.size(); i++) {
          // measurements for unit conversion
          surface = supports.at(i).getSurface();
          if (supports.at(i).isVertical()) length = supports.at(i).getRWidth();
          else length = supports.at(i).getZLength();
          std::vector<SingleSup>::const_iterator iter, guard = internals.supinfo.end();
          // materials loop
          for (iter = internals.supinfo.begin(); iter != guard; iter++) {
            if (iter->cM == supports.at(i).getCategory()) {
              // unit conversion (internal unit is grammes)
              double M;
              if (iter->uM == grpm) M = convert(iter->M, iter->uM, length);
              else M = convert(iter->M, iter->uM, mt.getMaterial(iter->tag).density, surface);
              // save converted and scaled material
              supports.at(i).addLocalMass(iter->tag, M);
            }
          }
          // use recorded materials to calculate material properties
          supports.at(i).calculateTotalMass();
          supports.at(i).calculateRadiationLength(mt);
          supports.at(i).calculateInteractionLength(mt);
        }
      }
      catch(std::exception& e) {
        std::cerr << "MatCalc::calculateSupportMaterials(): " << e.what() << " " << msg_abort << std::endl;
        return false;
      }
      return true;
    }

    /**
     * This function prints a summary of the internal state of the material calculator.
     */
    void MatCalc::printInternals() {
      if (init_done) {
        std::cout << "-----MatCalc parsed parameters-----" << std::endl << std::endl;
        std::cout << "Number of module types found: " << internals.modinfo.size() << std::endl;
        std::vector<std::pair<TypeInfo, std::vector<SingleMod> > >::const_iterator i, g = internals.modinfo.end();
        for (i = internals.modinfo.begin(); i != g; i++) {
          std::cout << "Type " << i->first.type << ": " << i->first.strips_across << " strips, " << i->first.segments_along << " segments." << std::endl;
          std::cout << "Number of type " << i->first.type << " materials found: " << i->second.size() << std::endl;
          for (std::vector<SingleMod>::const_iterator iter = i->second.begin(); iter != i->second.end(); iter++) {
            std::cout << iter->tag << ": " << (iter->is_local ? "local" : "exiting") << " material with A = " << iter->A << " (unit ";
            std::cout << iter->uA << "), B = " << iter->B << " (unit " << iter->uB << "), C = " << iter->C << " (unit " << iter->uC;
            std::cout << "), D = " << iter->D << " (unit " << iter->uD << ")." << std::endl;
          }
        }
        std::cout << std::endl << "Number of parsed entries for local services: " << internals.serlocalinfo.size() << std::endl;
        for (std::vector<SingleSerLocal>::const_iterator iter = internals.serlocalinfo.begin(); iter != internals.serlocalinfo.end(); iter++) {
          std::cout << iter->tag << ": Q = " << iter->Q << " (unit " << iter->uQ << ")." << std::endl;
        }
        std::cout << std::endl << "Number of parsed entries for exiting services: " << internals.serexitinfo.size() << std::endl;
        for (std::vector<SingleSerExit>::const_iterator iter = internals.serexitinfo.begin(); iter != internals.serexitinfo.end(); iter++) {
          std::cout << iter->In << " (unit " << iter->uIn << ") " << iter->tagIn << " --> " << iter->Out << " (unit " << iter->uOut;
          std::cout << ") " << iter->tagOut << " (" << (iter->is_local ? "local" : "not local") << " to service)." << std::endl;
        }
        std::cout << std::endl << "Number of support materials found: " << internals.supinfo.size() << std::endl;
        for (std::vector<SingleSup>::const_iterator iter = internals.supinfo.begin(); iter != internals.supinfo.end(); iter++) {
          std::cout << iter->tag << ": M = " << iter->M << " (unit " << iter->uM << ", category " << iter->cM << ")." << std::endl;
        }
      }
      else std::cout << "No parameters have been parsed into the internal data structures so far." << std::endl;
    }

    // protected
    /**
     * Getter for the module info vector of the given type. If that type is unknown or not on the list,
     * an exception is thrown.
     * @param type The module type of the requested info vector
     * @return A reference to the requested module info vector
     */
    std::vector<MatCalc::SingleMod>& MatCalc::getModVector(std::string type) { // throws exception
      std::vector<std::pair<TypeInfo, std::vector<SingleMod> > >::iterator iter = internals.modinfo.begin();
      std::vector<std::pair<TypeInfo, std::vector<SingleMod> > >::iterator guard = internals.modinfo.end();
      while (iter != guard) {
        if (type.compare(iter->first.type) == 0) return iter->second;
        iter++;
      }
      throw std::range_error("Exception in MatCalc::getModVector(): " + err_unknown_type);
    }

    /**
     * Getter for the complete type information of a given module type. If that type is unknown or not on
     * the list, an exception is thrown.
     * @param type The module type of the requested vector entry
     * @return A reference to the requested module type information
     */
    MatCalc::TypeInfo& MatCalc::getTypeInfoByType(std::string type) { // throws exception
      std::vector<std::pair<TypeInfo, std::vector<SingleMod> > >::iterator iter = internals.modinfo.begin();
      std::vector<std::pair<TypeInfo, std::vector<SingleMod> > >::iterator guard = internals.modinfo.end();
      while (iter != guard) {
        if (type.compare(iter->first.type) == 0) return iter->first;
        iter++;
      }
      throw std::range_error("Exception in MatCalc::getTypeInfoTyType(): " + err_no_such_type);
    }

    /**
     * Getter for a single material entry from one of the module info vectors. If the given parameter combination
     * does not exist on the list, an exception is thrown.
     * @param tag A unique material identifier
     * @param type The module type that the material should be applied to
     * @param uA The unit of parameter <i>A</i>
     * @param uB The unit of parameter <i>B</i>
     * @param uC The unit of parameter <i>C</i>
     * @param uD The unit of parameter <i>D</i>
     * @param local A flag indicating whether the material in question is local or exiting
     * @return A reference to the requested entry
     */
    MatCalc::SingleMod& MatCalc::getSingleMod(std::string tag, std::string type, std::string comp,
                                              Matunit uA, Matunit uB, Matunit uC, Matunit uD, bool local) { // throws exception
      std::vector<SingleMod>& vect = getModVector(type);
      std::vector<SingleMod>::iterator iter = vect.begin();
      while (iter != vect.end()) {
        if ((tag.compare(iter->tag) == 0) && (iter->is_local == local)) {
          if (comp.compare(iter->comp) == 0)
            if ((iter->uA == uA) && (iter->uB == uB) && (iter->uC == uC) && (iter->uD == uD)) return *iter;
        }
        iter++;
      }
      throw std::range_error("Exception in MatCalc::getSingleMod(): " + err_no_material);
    }

    /**
     * Getter for a single entry from the info vector for materials local to the layer/service boundary. If
     * the given parameter combination does not exist on the list, an exception is thrown.
     * @param tag A unique material identifier
     * @param u The unit of the requested material
     * @return A reference to the requested entry
     */
    MatCalc::SingleSerLocal& MatCalc::getSingleSer(std::string tag, Matunit u) { // throws exception
      std::vector<SingleSerLocal>::iterator iter = internals.serlocalinfo.begin();
      std::vector<SingleSerLocal>::iterator guard = internals.serlocalinfo.end();
      while (iter != guard) {
        if ((tag.compare(iter->tag) == 0) && (iter->uQ == u)) return *iter;
        iter++;
      }
      throw std::range_error("Exception in MatCalc::getSingleSer(2): " + err_no_service);
    }

    /**
     * Getter for a single entry from the info vector for material mappings at the layer/service boundary. If
     * the given parameter combination does not exist on the list, an exception is thrown.
     * @param tag1 A unique material identifier for the source material
     * @param tag2 A unique material identifier for the destination material
     * @param u1 The unit of the source material
     * @param u2 The unit of the destination material
     * @param local A flag indicating whether the destination material stays local to the service or exits it
     * @return A reference to the requested entry
     */
    MatCalc::SingleSerExit& MatCalc::getSingleSer(std::string tag1, std::string tag2, Matunit u1, Matunit u2, bool local) { // throws exception
      std::vector<SingleSerExit>::iterator iter = internals.serexitinfo.begin();
      std::vector<SingleSerExit>::iterator guard = internals.serexitinfo.end();
      while (iter != guard) {
        if ((tag1.compare(iter->tagIn) == 0) && (tag2.compare(iter->tagOut) == 0)) {
          if ((iter->uIn == u1) && (iter->uOut == u2)) {
            if (iter->is_local == local) return *iter;
          }
        }
        iter++;
      }
      throw std::range_error("Exception in MatCalc::getSingleSer(5): " + err_no_service);
    }

    /**
     * Getter for a single entry from the support parts info vector. If the given parameter combination does
     * not exist, an exception is thrown.
     * @param tag A unique material identifier
     * @param uM The unit of the requested material
     * @param cM The category of support parts that the material belongs to
     * @return A reference to the requested entry
     */
    MatCalc::SingleSup& MatCalc::getSingleSup(std::string tag, Matunit uM, MaterialProperties::Category cM) { // throws exception
      std::vector<SingleSup>::iterator iter = internals.supinfo.begin();
      std::vector<SingleSup>::iterator guard = internals.supinfo.end();
      while (iter != guard) {
        if (tag.compare(iter->tag) == 0) {
          if ((iter->uM == uM) && (iter->cM == cM)) return *iter;
        }
        iter++;
      }
      throw std::range_error("Exception in MatCalc::getSingleSup(): " + err_no_support);
    }

    // private
    /**
     * Check if there is an entry for a given type in the internal type info list.
     * @param type The requested module type
     * @return True if such an entry exists, false otherwise
     */
    bool MatCalc::entryExists(std::string type) {
      std::vector<std::pair<TypeInfo, std::vector<SingleMod> > >::iterator iter = internals.modinfo.begin();
      std::vector<std::pair<TypeInfo, std::vector<SingleMod> > >::iterator guard = internals.modinfo.end();
      while (iter != guard) {
        if (type.compare(iter->first.type) == 0) return true;
        iter++;
      }
      return false;
    }

    /**
     * Check if there is an entry of the given specifics in one of  the internal lists of module materials.
     * @param tag A unique material identifier
     * @param type The module type that the material should be applied to
     * @param uA The unit of parameter <i>A</i>
     * @param uB The unit of parameter <i>B</i>
     * @param uC The unit of parameter <i>C</i>
     * @param uD The unit of parameter <i>D</i>
     * @param local A flag indicating whether the material in question is local or exiting
     * @return True if such an entry exists, false otherwise
     */
    bool MatCalc::entryExists(std::string tag, std::string type, std::string comp, Matunit uA, Matunit uB, Matunit uC, Matunit uD, bool local) {
      try {
        std::vector<SingleMod>& vect = getModVector(type);
        std::vector<SingleMod>::const_iterator iter = vect.begin();
        std::vector<SingleMod>::const_iterator guard = vect.end();
        while (iter != guard) {
          bool found = (tag.compare(iter->tag) == 0) && (iter->is_local == local);
          found = found && (comp.compare(iter->comp) == 0);
          found = found && (iter->uA == uA) && (iter->uB == uB) && (iter->uC == uC) && (iter->uD == uD);
          if (found) return true;
          iter++;
        }
      }
      catch (std::range_error& re) {
        std::cerr << "MatCalc::entryExists(7): " << re.what() << std::endl;
        return false;
      }
      return false;
    }

    /**
     * Check if there is an entry of the given specifics in the internal list of local materials at the layer/service boundary.
     * @param tag A unique material identifier
     * @param uQ The unit of the material in question
     * @return True if such an entry exists, false otherwise
     */
    bool MatCalc::entryExists(std::string tag, Matunit uQ) {
      std::vector<SingleSerLocal>::iterator iter = internals.serlocalinfo.begin();
      std::vector<SingleSerLocal>::iterator guard = internals.serlocalinfo.end();
      while (iter != guard) {
        if ((tag.compare(iter->tag) == 0) && (iter->uQ == uQ)) return true;
        iter++;
      }
      return false;
    }

    /**
     * Check if there is an entry of the given specifics in the internal list of material mappings for the layer/service boundary.
     * @param tag1 A unique identifier for the source material that the mapping starts from
     * @param tag2 A unique identifier for the destination material that the source material is converted to
     * @param u1 The unit of the source material
     * @param u2 The unit of the destination material
     * @param local A flag indicating whether the destination material stays local to the service or exits it
     * @return True if such an entry exists, false otherwise
     */
    bool MatCalc::entryExists(std::string tag1, std::string tag2, Matunit u1, Matunit u2, bool local) {
      std::vector<SingleSerExit>::iterator iter = internals.serexitinfo.begin();
      std::vector<SingleSerExit>::iterator guard = internals.serexitinfo.end();
      while (iter != guard) {
        bool found = (tag1.compare(iter->tagIn) == 0) && (tag2.compare(iter->tagOut) == 0);
        found = found && (iter->uIn == u1) && (iter->uOut == u2) && (iter->is_local == local);
        if (found) return true;
        iter++;
      }
      return false;
    }

    /**
     * Check if there is an entry of the given specifics in the internal list of support materials.
     * @param tag A unique material identifier
     * @param uM The unit of the requested material
     * @param cM The category of support parts that the material belongs to
     * @return True if such an entry exists, false otherwise
     */
    bool MatCalc::entryExists(std::string tag, Matunit uM, MaterialProperties::Category cM) {
      std::vector<SingleSup>::iterator iter = internals.supinfo.begin();
      std::vector<SingleSup>::iterator guard = internals.supinfo.end();
      while (iter != guard) {
        if ((tag.compare(iter->tag) == 0) && (iter->uM == uM) && (iter->cM == cM)) return true;
        iter++;
      }
      return false;
    }

    /**
     * This convenience function finds the number of rods in a layer.
     * @param caps The collection of <i>ModuleCap</i> objects that maps to a series of layers in a tracker
     * @param layer The layer under investigation
     * @return The number of rods in the given layer
     */
    int MatCalc::findBarrelRods(std::vector<std::vector<ModuleCap> >& caps, int layer) {
      return findEndcapRods(caps, layer) / 2;
    }

    /**
     * This convenience function finds the number of modules in the last ring of a disc.
     * @param caps The collection of <i>ModuleCap</i> objects that maps to a series of discs in a tracker
     * @param layer The disc under investigation
     * @return The number of modules in the outmost ring of the given disc
     */
    int MatCalc::findEndcapRods(std::vector<std::vector<ModuleCap> >& caps, int layer) {
      int res = 0;
      int index = 1;
      if ((layer >= 0) && (layer < (int)caps.size())) {
        for (unsigned int i = 0; i < caps.at(layer).size(); i++) {
          if (caps.at(layer).at(i).getModule().uniRef().ring > index) {
            index = caps.at(layer).at(i).getModule().uniRef().ring;
            res = 1;
          }
          else if (caps.at(layer).at(i).getModule().uniRef().ring == index) res++;
        }
      }
      return res;
    }

    /**
     * This is a convenience function that converts an amount of material from one of the units available to the
     * config file to the internal unit of grams. Since this normally requires additional information about the
     * geometry of the volume or about the material, those parameters must be supplied as well. If the caller
     * attempts to convert from an unknown unit, an exception is thrown.
     * @param value The amount of material that is to be converted
     * @param unit The unit that the conversion starts from
     * @param densityorlength The density of the material if the unit is in mm or mm3, the traversed length if it is in g/m
     * @param surface The surface of the tracker volume
     * @return The material mass after conversion
     */
    double MatCalc::convert(double value, Matunit unit, double densityorlength, double surface) { // throws exception
      switch(unit) {
      case gr : return value;
      case mm3 : return densityorlength * value / 1000.0;
      case mm : return densityorlength * surface * value / 1000.0;
      case grpm : return densityorlength * value / 1000.0;
      default : throw std::range_error("Exception in MatCalc::convert(): " + err_conversion);
      }
    }

    /**
     * This function propagates materials from a layer or disc to the first adjacent service.
     * @param source The layer or disc vector that the materials come out of
     * @param dest The service that is currently being processed
     * @param r The number of rods within the source layer
     * @param l The length that e.g. a cable would have to travel to cross along the service volume
     * @param s The surface of the service volume
     */
    void MatCalc::adjacentDifferentCategory(std::vector<ModuleCap>& source, InactiveElement& dest, int r, double l, double s) {
      // S-labelled service materials
      std::vector<SingleSerLocal>::const_iterator liter, lguard = internals.serlocalinfo.end();
      // materials loop
      for (liter = internals.serlocalinfo.begin(); liter != lguard; liter++) {
        // unit conversion (internal unit is grammes)
        double Q;
        if (liter->uQ == grpm) Q = convert(liter->Q, liter->uQ, l);
        else Q = convert(liter->Q, liter->uQ, mt.getMaterial(liter->tag).density, s);
        // save converted and scaled material
        dest.addLocalMass(liter->tag, r * Q);
      }
      // D-labelled service materials
      int modsonrod = 0, lastmod = 0;
      // loop to find information about contributing source modules
      for (unsigned int j = 0; j < source.size(); j++) {
        if (modsonrod < source.at(j).getModule().uniRef().ring) {
          // modsonrod finds the number of modules along a rod
          modsonrod = source.at(j).getModule().uniRef().ring;
          // lastmod finds the index of a sample module at the end of a rod
          lastmod = j;
        }
      }
      std::vector<SingleSerExit>::const_iterator eiter, eguard = internals.serexitinfo.end();
      // materials loop
      for (eiter = internals.serexitinfo.begin(); eiter != eguard; eiter++) {
        try {
          // unit conversion per parameter (internal unit is grammes)
          double In, Out;
          if (eiter->uIn == grpm) In = convert(eiter->In, eiter->uIn, source.at(lastmod).getModule().length());
          else In = convert(eiter->In, eiter->uIn, mt.getMaterial(eiter->tagIn).density, source.at(lastmod).getSurface());
          if (eiter->uOut == grpm) Out = convert(eiter->Out, eiter->uOut, l);
          else Out = convert(eiter->Out, eiter->uOut, mt.getMaterial(eiter->tagOut).density, s);
          Out = Out * source.at(lastmod).getExitingMass(eiter->tagIn) / In;
          if (eiter->is_local) dest.addLocalMass(eiter->tagOut, r * Out);
          // save converted and scaled material
          else dest.addExitingMass(eiter->tagOut, r * Out);
        }
        catch (std::runtime_error& re) {
          std::cerr << "MatCalc::adjacentDifferentCategory(): " << err_no_material << msg_ignore_tag << eiter->tagIn << "." << std::endl;
          std::cerr << "MatCalc::adjacentDifferentCategory(): exception was: " << re.what() << std::endl;
        }
      }
    }

    /**
     * This function propagates the materials of a service to the next.
     * @param source The neighbour volume of the current service
     * @param dest The service that is currently being processed
     */
    void MatCalc::adjacentSameCategory(InactiveElement& source, InactiveElement& dest) {
      double tmp, scalar;
      // find the multiplier to convert from the source surface to the destination surface
      scalar = dest.getSurface() / source.getSurface();
      for (std::map<std::string, double>::const_iterator it = source.getExitingMasses().begin(); it != source.getExitingMasses().end(); ++it) {
        // scale incoming material with respect to the new geometry
        tmp = scalar * it->second;
        // save scaled material
        dest.addExitingMass(it->first, tmp);
      }
    }

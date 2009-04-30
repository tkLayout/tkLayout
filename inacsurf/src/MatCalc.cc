/**
 * @file MatCalc.cc
 * @brief
 */

#include <MatCalc.h>
namespace insur {
    // public
    bool MatCalc::initDone() { return init_done; }
    
    void MatCalc::initDone(bool yes) { init_done = yes; }
    
    int MatCalc::getStripsAcross(Modtype type) {
        try {
            TypeInfo& info = getTypeInfoByType(type);
            return info.strips_across;
        }
        catch(std::range_error& re) {
            std::cerr << re.what() << std::endl;
            throw;
        }
    }
    
    int MatCalc::getSegmentsAlong(Modtype type) {
        try {
            TypeInfo& info = getTypeInfoByType(type);
            return info.segments_along;
        }
        catch(std::range_error& re) {
            std::cerr << re.what() << std::endl;
            throw;
        }
    }
    
    std::pair<int, int> MatCalc::getDefaultDimensions(Modtype type) {
        try {
            std::pair<int, int> result;
            TypeInfo& info = getTypeInfoByType(type);
            result.first = info.strips_across;
            result.second = info.segments_along;
            return result;
        }
        catch(std::range_error& re) {
            std::cerr << re.what() << std::endl;
            throw;
        }
    }
    
    void MatCalc::addTypeInfo(Modtype type, int strips, int segments) {
        if (!entryExists(type)) {
            TypeInfo info;
            info.type = type;
            info.strips_across = strips;
            info.segments_along = segments;
            internals.typeinfo.push_back(info);
        }
    }
    
    void MatCalc::updateTypeInfoStrips(Modtype type, int strips) {
        if (entryExists(type)) {
            try {
                TypeInfo& info = getTypeInfoByType(type);
                info.strips_across = strips;
            }
            catch (std::range_error re) {
                std::cerr << err_up_general << std::endl;
            }
        }
    }
    
    void MatCalc::updateTypeInfoSegments(Modtype type, int segments) {
        if (entryExists(type)) {
            try {
                TypeInfo& info = getTypeInfoByType(type);
                info.segments_along = segments;
            }
            catch (std::range_error re) {
                std::cerr << err_up_general << std::endl;
            }
        }
    }
    
    /**
     * Add a parameter entry to the vector corresponding to <i>type</i> by copying the supplied values.
     * Add the new values to the old ones if a material with the given parameters already appears on the list.
     * @param tag A unique material identifier
     * @param type The type of module the material parameters refer to
     * @param A The coefficient for the layout-dependent and travelling portion
     * @param B The coefficient for the layout-dependent and localised portion
     * @param C The coefficient for the per-module and travelling portion
     * @param D The coefficient for the per-module and localised portion
     * @param local A flag indicating whether the amounts as a whole will be treated as local or as exiting the module later
     */
    void MatCalc::addModuleParameters(std::string tag, Modtype type,
            double A, Matunit uA, double B, Matunit uB, double C, Matunit uC, double D, Matunit uD, bool local) {
        try {
            if (!entryExists(tag, type, uA, uB, uC, uD, local)) {
                SingleMod mod;
                mod.tag = tag;
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
                SingleMod& mod = getSingleMod(tag, type, uA, uB, uC, uD, local);
                mod.A = mod.A + A;
                mod.B = mod.B + B;
                mod.C = mod.C + C;
                mod.D = mod.D + D;
            }
        }
        catch(std::range_error& re) {
            std::cerr << err_matadd_weird << std::endl;
            std::cerr << re.what() << std::endl;
        }
    }
    
    void MatCalc::addServiceParameters(std::string tag, double Q, Matunit uQ) {
        if (!entryExists(tag, uQ)) {
            SingleSerLocal serl;
            serl.tag = tag;
            serl.Q = Q;
            serl.uQ = uQ;
            internals.serlocalinfo.push_back(serl);
        }
        else {
            SingleSerLocal& serl = getSingleSer(tag, uQ);
            serl.Q = serl.Q + Q;
        }
    }
    
    void MatCalc::addServiceParameters(std::string tagIn, double In, Matunit uIn,
            std::string tagOut, double Out, Matunit uOut, bool local) {
        if (!entryExists(tagIn, tagOut, uIn, uOut, local)) {
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
            SingleSerExit& ser = getSingleSer(tagIn, tagOut, uIn, uOut, local);
            ser.In = ser.In + In;
            ser.Out = ser.Out + Out;
        }
    }
    
    void MatCalc::addSupportParameters(std::string tag, double M, Matunit uM, MaterialProperties::Category cM) {
        if (!entryExists(tag, uM, cM)) {
            SingleSup sup;
            sup.tag = tag;
            sup.M = M;
            sup.uM = uM;
            sup.cM = cM;
            internals.supinfo.push_back(sup);
        }
        else {
            SingleSup& sup = getSingleSup(tag, uM, cM);
            sup.M = sup.M + M;
        }
    }
    
    void MatCalc::clearTypeVector() { internals.typeinfo.clear(); }
    
    void MatCalc::clearModVectors() {
        internals.modinforphi.clear();
        internals.modinfostereo.clear();
        internals.modinfopt.clear();
    }
    
    void MatCalc::clearModVector(Modtype type) {
        try {
            std::vector<SingleMod>& vect = getModVector(type);
            vect.clear();
        }
        catch(std::range_error& re) {
            std::cerr << re.what() << std::endl;
        }
    }
    
    void MatCalc::copyContents(Modtype source, Modtype dest) {
        clearModVector(dest);
        appendContents(source, dest);
    }
    
    void MatCalc::appendContents(Modtype source, Modtype dest) {
        if (source != dest) {
            try {
                std::vector<SingleMod>& in = getModVector(source);
                std::vector<SingleMod>& out = getModVector(dest);
                std::copy(in.begin(), in.end(), std::back_inserter(out));
            }
            catch(std::range_error& re) {
                std::cout << re.what() << std::endl;
            }
        }
    }
    
    bool MatCalc::typeRegistered(Modtype type) {
        std::vector<TypeInfo>::const_iterator iter = internals.typeinfo.begin();
        std::vector<TypeInfo>::const_iterator guard = internals.typeinfo.end();
        while (iter != guard) {
            if (iter->type == type) return true;
            else iter++;
        }
        return false;
    }
    
    uint MatCalc::registeredTypes() { return internals.typeinfo.size(); }
    
    MaterialTable& MatCalc::getMaterialTable() { return mt; }
    
    bool MatCalc::calculateBarrelMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps) {
        for (uint i = 0; i < barrelcaps.size(); i++) {
            if (barrelcaps.at(i).size() > 0) {
                Modtype mtype;
                if (barrelcaps.at(i).at(0).getModule().getType().compare(type_rphi) == 0) mtype = rphi;
                else if(barrelcaps.at(i).at(0).getModule().getType().compare(type_stereo) == 0) mtype = stereo;
                else if(barrelcaps.at(i).at(0).getModule().getType().compare(type_pt) == 0) mtype = pt;
                else {
                    std::cerr << err_unknown_type << " " << msg_abort << std::endl;
                    return false;
                }
                try {
                    double stripseg_scalar = (double)barrelcaps.at(i).at(0).getModule().getNStripsAcross() / (double)getStripsAcross(mtype);
                    stripseg_scalar = stripseg_scalar * (double)barrelcaps.at(i).at(0).getModule().getNSegments() / (double)getSegmentsAlong(mtype);
                    for (uint j = 0; j < barrelcaps.at(i).size(); j++) {
                        std::vector<SingleMod>& vect = getModVector(mtype);
                        std::vector<SingleMod>::const_iterator guard = vect.end();
                        std::vector<SingleMod>::const_iterator iter;
                        int index = barrelcaps.at(i).at(j).getModule().getRing(); // assuming physicists' counting scheme: starting at 1!!!
                        for (iter = vect.begin(); iter != guard; iter++) {
                            double A, B, C, D;
                            double density, surface, length;
                            density = mt.getMaterial(iter->tag).density;
                            surface = barrelcaps.at(i).at(j).getSurface();
                            length = barrelcaps.at(i).at(j).getModule().getHeight();
                            if (iter->uA == grpm) A = convert(iter->A, iter->uA, length);
                            else A = convert(iter->A, iter->uA, density, surface);
                            if (iter->uB == grpm) B = convert(iter->B, iter->uB, length);
                            else B = convert(iter->B, iter->uB, density, surface);
                            if (iter->uC == grpm) C = convert(iter->C, iter->uC, length);
                            else C = convert(iter->C, iter->uC, density, surface);
                            if (iter->uD == grpm) D = convert(iter->D, iter->uD, length);
                            else D = convert(iter->D, iter->uD, density, surface);
                            A = A * stripseg_scalar * (double)(index - 1);
                            B = B * stripseg_scalar;
                            C = C * (double)(index - 1);
                            if (iter->is_local) barrelcaps.at(i).at(j).addLocalMass(iter->tag, A + B + C + D);
                            else barrelcaps.at(i).at(j).addExitingMass(iter->tag, A + B + C + D);
                        }
                        barrelcaps.at(i).at(j).calculateTotalMass();
                        barrelcaps.at(i).at(j).calculateRadiationLength(mt);
                        barrelcaps.at(i).at(j).calculateInteractionLength(mt);
                    }
                }
                catch(std::range_error re) {
                    std::cerr << re.what() << " " << msg_abort << std::endl;
                    return false;
                }
            }
        }
        return true;
    }
    
    bool MatCalc::calculateEndcapMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps) {
        for (uint i = 0; i < endcapcaps.size(); i++) {
            if (endcapcaps.at(i).size() > 0) {
                try {
                    int rindex;
                    std::vector<int> mods;
                    std::vector<double> stripseg_scalars;
                    std::vector<Modtype> mtypes;
                    std::vector<std::list<int> > modinrings;
                    for (uint j = 0; j < endcapcaps.at(i).size(); j++) {
                        rindex = endcapcaps.at(i).at(j).getModule().getRing();
                        if ((int)mods.size() < rindex) {
                            while ((int)mods.size() < rindex) mods.push_back(0);
                        }
                        mods.at(rindex - 1) = mods.at(rindex - 1) + 1;
                        if ((int)mtypes.size() < rindex) {
                            while ((int)mtypes.size() < rindex) mtypes.push_back(un_mod);
                        }
                        if (mtypes.at(rindex - 1) == un_mod) {
                            if (endcapcaps.at(i).at(0).getModule().getType().compare(type_rphi) == 0) mtypes.at(rindex - 1) = rphi;
                            else if(endcapcaps.at(i).at(0).getModule().getType().compare(type_stereo) == 0) mtypes.at(rindex - 1) = stereo;
                            else if(endcapcaps.at(i).at(0).getModule().getType().compare(type_pt) == 0) mtypes.at(rindex - 1) = pt;
                            else {
                                std::cerr << err_unknown_type << " Encountered type value '" << endcapcaps.at(i).at(j).getModule().getType() << "'" << std::endl;
                                return false;
                            }
                        }
                        if ((int)stripseg_scalars.size() < rindex) {
                            while ((int)stripseg_scalars.size() < rindex) stripseg_scalars.push_back(0.0);
                        }
                        if (stripseg_scalars.at(rindex - 1) == 0.0) {
                            if (mtypes.at(rindex - 1) != un_mod) {
                                stripseg_scalars.at(rindex - 1) = (double)endcapcaps.at(i).at(j).getModule().getNStripsAcross();
                                stripseg_scalars.at(rindex - 1) = stripseg_scalars.at(rindex - 1) / (double)getStripsAcross(mtypes.at(rindex - 1));
                                stripseg_scalars.at(rindex - 1) = stripseg_scalars.at(rindex - 1) * (double)endcapcaps.at(i).at(j).getModule().getNSegments();
                                stripseg_scalars.at(rindex - 1) = stripseg_scalars.at(rindex - 1) / (double)getSegmentsAlong(mtypes.at(rindex - 1));
                            }
                        }
                        if ((int)modinrings.size() < rindex) {
                            while ((int)modinrings.size() < rindex) {
                                std::list<int> tmp;
                                modinrings.push_back(tmp);
                            }
                        }
                        modinrings.at(rindex - 1).push_back(j);
                    }
                    rindex = modinrings.size();
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
                            else {
                                length = endcapcaps.at(i).at(*first).getModule().getHeight();
                                std::vector<SingleMod>& vect = getModVector(mtypes.at(j));
                                std::vector<SingleMod>::const_iterator guard = vect.end();
                                std::vector<SingleMod>::const_iterator iter;
                                for (iter = vect.begin(); iter != guard; iter++) {
                                    density = mt.getMaterial(iter->tag).density;
                                    if (iter->uB == grpm) B = convert(iter->B, iter->uB, length);
                                    else B = convert(iter->B, iter->uB, density, surface);
                                    if (iter->uD == grpm) D = convert(iter->D, iter->uD, length);
                                    else D = convert(iter->D, iter->uD, density, surface);
                                    B = B * stripseg_scalars.at(j);
                                    if (iter->is_local) endcapcaps.at(i).at(*first).addLocalMass(iter->tag, B + D);
                                    else endcapcaps.at(i).at(*first).addExitingMass(iter->tag, B + D);
                                }
                                if (j > 0) {
                                    for (int k = 0; k < j; k++) {
                                        if (mods.at(k) > 0) {
                                            surface = endcapcaps.at(i).at(modinrings.at(k).front()).getSurface();
                                            if (surface < 0) {
                                                std::cerr << msg_negative_area << " Endcap module in disc " << i << ", ring " << k;
                                                std::cerr << " with index " << *first << " within the disc. " << msg_abort << std::endl;
                                                return false;
                                            }
                                            else {
                                                length = endcapcaps.at(i).at(modinrings.at(k).front()).getModule().getHeight();
                                                std::vector<SingleMod>& vect = getModVector(mtypes.at(j));
                                                std::vector<SingleMod>::const_iterator guard = vect.end();
                                                std::vector<SingleMod>::const_iterator iter;
                                                for (iter = vect.begin(); iter != guard; iter++) {
                                                    density = mt.getMaterial(iter->tag).density;
                                                    if (iter->uA == grpm) A = convert(iter->A, iter->uA, length);
                                                    else A = convert(iter->A, iter->uA, density, surface);
                                                    if (iter->uC == grpm) C = convert(iter->C, iter->uC, length);
                                                    else C = convert(iter->C, iter->uC, density, surface);
                                                    A = A * (double)mods.at(k) / (double)mods.at(j) * stripseg_scalars.at(k);
                                                    C = C * (double)mods.at(k) / (double)mods.at(j);
                                                    if (iter->is_local) endcapcaps.at(i).at(*first).addLocalMass(iter->tag, A + C);
                                                    else endcapcaps.at(i).at(*first).addExitingMass(iter->tag, A + C);
                                                }
                                            }
                                        }
                                    }
                                }
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
                catch (std::range_error re) {
                    std::cerr << re.what() << " " << msg_abort << std::endl;
                    return false;
                }
            }
        }
        return true;
    }
    
    bool MatCalc::calculateBarrelServiceMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps,
            std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices) {
        int feeder, neighbour;
        InactiveElement::InType ftype, ntype;
        double length, surface;
        for(uint i = 0; i < barrelservices.size(); i++) {
            try {
                ftype = barrelservices.at(i).getFeederType();
                feeder = barrelservices.at(i).getFeederIndex();
                ntype = barrelservices.at(i).getNeighbourType();
                neighbour = barrelservices.at(i).getNeighbourIndex();
                if (barrelservices.at(i).isVertical()) length = barrelservices.at(i).getRWidth();
                else length = barrelservices.at(i).getZLength();
                surface = barrelservices.at(i).getSurface();
                // feeder
                switch(ftype) {
                    case InactiveElement::no_in : break;
                    case InactiveElement::tracker : adjacentDifferentCategory(barrelcaps.at(feeder), barrelservices.at(i), findRods(barrelcaps, feeder), length, surface);
                    break;
                    case InactiveElement::barrel : adjacentSameCategory(barrelservices.at(feeder), barrelservices.at(i));
                    break;
                    case InactiveElement::endcap : adjacentSameCategory(endcapservices.at(feeder), barrelservices.at(i));
                }
                // neighbour
                switch(ntype) {
                    case InactiveElement::no_in : break;
                    case InactiveElement::tracker : adjacentDifferentCategory(barrelcaps.at(neighbour), barrelservices.at(i), findRods(barrelcaps, neighbour), length, surface);
                    break;
                    case InactiveElement::barrel : adjacentSameCategory(barrelservices.at(neighbour), barrelservices.at(i));
                    break;
                    case InactiveElement::endcap : adjacentSameCategory(endcapservices.at(neighbour), barrelservices.at(i));
                }
                barrelservices.at(i).calculateTotalMass();
                barrelservices.at(i).calculateRadiationLength(mt);
                barrelservices.at(i).calculateInteractionLength(mt);
            }
            catch(std::runtime_error re) {
                std::cerr << re.what() << " " << msg_abort << std::endl;
                return false;
            }
            catch(std::exception e) {
                std::cerr << e.what() << " " << msg_abort << std::endl;
                return false;
            }
        }
        return true;
    }
    
    bool MatCalc::calculateEndcapServiceMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps,
            std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices) {
        int feeder, neighbour;
        InactiveElement::InType ftype, ntype;
        double length, surface;
        for (uint i = 0; i < endcapservices.size(); i++) {
            try {
                ftype = endcapservices.at(i).getFeederType();
                feeder = endcapservices.at(i).getFeederIndex();
                ntype = endcapservices.at(i).getNeighbourType();
                neighbour = endcapservices.at(i).getNeighbourIndex();
                if (endcapservices.at(i).isVertical()) length = endcapservices.at(i).getRWidth();
                else length = endcapservices.at(i).getZLength();
                surface = endcapservices.at(i).getSurface();
                // feeder
                switch(ftype) {
                    case InactiveElement::no_in : break;
                    case InactiveElement::tracker : adjacentDifferentCategory(endcapcaps.at(feeder), endcapservices.at(i), findRods(endcapcaps, feeder), length, surface);
                    break;
                    case InactiveElement::barrel : adjacentSameCategory(barrelservices.at(feeder), endcapservices.at(i));
                    break;
                    case InactiveElement::endcap : adjacentSameCategory(endcapservices.at(feeder), endcapservices.at(i));
                }
                // neighbour
                switch(ntype) {
                    case InactiveElement::no_in : break;
                    case InactiveElement::tracker : adjacentDifferentCategory(endcapcaps.at(neighbour), endcapservices.at(i), findRods(endcapcaps, neighbour), length, surface);
                    break;
                    case InactiveElement::barrel : adjacentSameCategory(barrelservices.at(neighbour), endcapservices.at(i));
                    break;
                    case InactiveElement::endcap : adjacentSameCategory(endcapservices.at(neighbour), endcapservices.at(i));
                }
                endcapservices.at(i).calculateTotalMass();
                endcapservices.at(i).calculateRadiationLength(mt);
                endcapservices.at(i).calculateInteractionLength(mt);
            }
            catch(std::runtime_error re) {
                std::cerr << re.what() << " " << msg_abort << std::endl;
                return false;
            }
            catch(std::exception e) {
                std::cerr << e.what() << " " << msg_abort << std::endl;
                return false;
            }
        }
        return true;
    }
    
    bool MatCalc::calculateSupportMaterials(std::vector<InactiveElement>& supports) {
        double length, surface;
        try {
            for (uint i = 0; i < supports.size(); i++) {
                surface = supports.at(i).getSurface();
                if (supports.at(i).isVertical()) length = supports.at(i).getRWidth();
                else length = supports.at(i).getZLength();
                std::vector<SingleSup>::const_iterator iter, guard = internals.supinfo.end();
                for (iter = internals.supinfo.begin(); iter != guard; iter++) {
                    if (iter->cM == supports.at(i).getCategory()) {
                        int M;
                        if (iter->uM == grpm) M = convert(iter->M, iter->uM, length);
                        else M = convert(iter->M, iter->uM, mt.getMaterial(iter->tag).density, surface);
                        supports.at(i).addLocalMass(iter->tag, M);
                    }
                }
                supports.at(i).calculateTotalMass();
                supports.at(i).calculateRadiationLength(mt);
                supports.at(i).calculateInteractionLength(mt);
            }
        }
        catch(std::exception e) {
            std::cerr << e.what() << " " << msg_abort << std::endl;
            return false;
        }
        return true;
    }
    
    void MatCalc::printInternals() {
        if (init_done) {
            std::cout << "-----MatCalc parsed parameters-----" << std::endl << std::endl;
            std::cout << "Number of module types found: " << internals.typeinfo.size() << std::endl;
            for (std::vector<TypeInfo>::const_iterator iter = internals.typeinfo.begin(); iter != internals.typeinfo.end(); iter++) {
                std::cout << "Type " << iter->type << ": " << iter->strips_across << " strips across, ";
                std::cout << iter->segments_along << " segments along." << std::endl;
            }
            std::cout << std::endl << "Number of rphi materials found: " << internals.modinforphi.size() << std::endl;
            for (std::vector<SingleMod>::const_iterator iter = internals.modinforphi.begin(); iter != internals.modinforphi.end(); iter++) {
                std::cout << iter->tag << ": " << (iter->is_local ? "local" : "exiting") << " material with A = " << iter->A << " (unit ";
                std::cout << iter->uA << "), B = " << iter->B << " (unit " << iter->uB << "), C = " << iter->C << " (unit " << iter->uC;
                std::cout << "), D = " << iter->D << " (unit " << iter->uD << ")." << std::endl;
            }
            std::cout << std::endl << "Number of stereo materials found: " << internals.modinfostereo.size() << std::endl;
            for (std::vector<SingleMod>::const_iterator iter = internals.modinfostereo.begin(); iter != internals.modinfostereo.end(); iter++) {
                std::cout << iter->tag << ": " << (iter->is_local ? "local" : "exiting") << " material with A = " << iter->A << " (unit ";
                std::cout << iter->uA << "), B = " << iter->B << " (unit " << iter->uB << "), C = " << iter->C << " (unit " << iter->uC;
                std::cout << "), D = " << iter->D << " (unit " << iter->uD << ")." << std::endl;
            }
            std::cout << std::endl << "Number of pt materials found: " << internals.modinfopt.size() << std::endl;
            for (std::vector<SingleMod>::const_iterator iter = internals.modinfopt.begin(); iter != internals.modinfopt.end(); iter++) {
                std::cout << iter->tag << ": " << (iter->is_local ? "local" : "exiting") << " material with A = " << iter->A << " (unit ";
                std::cout << iter->uA << "), B = " << iter->B << " (unit " << iter->uB << "), C = " << iter->C << " (unit " << iter->uC;
                std::cout << "), D = " << iter->D << " (unit " << iter->uD << ")." << std::endl;
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
    std::vector<MatCalc::SingleMod>& MatCalc::getModVector(Modtype type) { // throws exception
        switch(type) {
            case rphi : return internals.modinforphi;
            case stereo : return internals.modinfostereo;
            case pt : return internals.modinfopt;
            default : throw std::range_error(err_unknown_type);
        }
    }
    
    MatCalc::TypeInfo& MatCalc::getTypeInfoByType(Modtype type) { // throws exception
        std::vector<TypeInfo>::iterator iter = internals.typeinfo.begin();
        std::vector<TypeInfo>::iterator guard = internals.typeinfo.end();
        while (iter != guard) {
            if (iter->type == type) return *iter;
            iter++;
        }
        throw std::range_error(err_no_such_type);
    }
    
    MatCalc::SingleMod& MatCalc::getSingleMod(std::string tag, Modtype type, Matunit uA, Matunit uB, Matunit uC, Matunit uD, bool local) { // throws exception
        std::vector<SingleMod>& vect = getModVector(type);
        std::vector<SingleMod>::iterator iter = vect.begin();
        while (iter != vect.end()) {
            if ((tag.compare(iter->tag) == 0) && (iter->is_local == local)) {
                if ((iter->uA == uA) && (iter->uB == uB) && (iter->uC == uC) && (iter->uD == uD)) return *iter;
            }
            iter++;
        }
        throw std::range_error(err_no_material);
    }
    
    MatCalc::SingleSerLocal& MatCalc::getSingleSer(std::string tag, Matunit u) { // throws exception
        std::vector<SingleSerLocal>::iterator iter = internals.serlocalinfo.begin();
        std::vector<SingleSerLocal>::iterator guard = internals.serlocalinfo.end();
        while (iter != guard) {
            if ((tag.compare(iter->tag) == 0) && (iter->uQ == u)) return *iter;
            iter++;
        }
        throw std::range_error(err_no_service);
    }
    
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
        throw std::range_error(err_no_service);
    }
    
    MatCalc::SingleSup& MatCalc::getSingleSup(std::string tag, Matunit uM, MaterialProperties::Category cM) { // throws exception
        std::vector<SingleSup>::iterator iter = internals.supinfo.begin();
        std::vector<SingleSup>::iterator guard = internals.supinfo.end();
        while (iter != guard) {
            if (tag.compare(iter->tag) == 0) {
                if ((iter->uM == uM) && (iter->cM == cM)) return *iter;
            }
            iter++;
        }
        throw std::range_error(err_no_support);
    }
    
    // private
    bool MatCalc::entryExists(Modtype type) {
        std::vector<TypeInfo>::iterator iter = internals.typeinfo.begin();
        std::vector<TypeInfo>::iterator guard = internals.typeinfo.end();
        while (iter != guard) {
            if (iter->type == type) return true;
            iter++;
        }
        return false;
    }
    
    bool MatCalc::entryExists(std::string tag, Modtype type, Matunit uA, Matunit uB, Matunit uC, Matunit uD, bool local) {
        try {
            std::vector<SingleMod>& vect = getModVector(type);
            std::vector<SingleMod>::const_iterator iter = vect.begin();
            std::vector<SingleMod>::const_iterator guard = vect.end();
            while (iter != guard) {
                bool found = (tag.compare(iter->tag) == 0) && (iter->is_local == local);
                found = found && (iter->uA == uA) && (iter->uB == uB) && (iter->uC == uC) && (iter->uD == uD);
                if (found) return true;
                iter++;
            }
        }
        catch (std::range_error& re) {
            std::cerr << re.what() << std::endl;
            return false;
        }
        return false;
    }
    
    bool MatCalc::entryExists(std::string tag, Matunit uQ) {
        std::vector<SingleSerLocal>::iterator iter = internals.serlocalinfo.begin();
        std::vector<SingleSerLocal>::iterator guard = internals.serlocalinfo.end();
        while (iter != guard) {
            if ((tag.compare(iter->tag) == 0) && (iter->uQ == uQ)) return true;
            iter++;
        }
        return false;
    }
    
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
    
    bool MatCalc::entryExists(std::string tag, Matunit uM, MaterialProperties::Category cM) {
        std::vector<SingleSup>::iterator iter = internals.supinfo.begin();
        std::vector<SingleSup>::iterator guard = internals.supinfo.end();
        while (iter != guard) {
            if ((tag.compare(iter->tag) == 0) && (iter->uM == uM) && (iter->cM == cM)) return true;
            iter++;
        }
        return false;
    }
    
    int MatCalc::findRods(std::vector<std::vector<ModuleCap> >& caps, int layer) {
        int res = 0;
        int index = 1;
        if (layer < (int)caps.size()) {
            for (uint i = 0; i < caps.at(layer).size(); i++) {
                if (caps.at(layer).at(i).getModule().getRing() > index) {
                    index = caps.at(layer).at(i).getModule().getRing();
                    res = 1;
                }
                else if (caps.at(layer).at(i).getModule().getRing() == index) res++;
            }
        }
        return res;
    }
    
    double MatCalc::convert(double value, Matunit unit, double densityorlength, double surface) { // throws exception
        switch(unit) {
            case gr : return value;
            case mm3 : return densityorlength * value;
            case mm : return densityorlength * surface * value;
            case grpm : return densityorlength * value;
            default : throw std::range_error(err_conversion);
        }
    }
    
    void MatCalc::adjacentDifferentCategory(std::vector<ModuleCap>& source, InactiveElement& dest, int r, double l, double s) {
        // S-labelled service materials
        std::vector<SingleSerLocal>::const_iterator liter, lguard = internals.serlocalinfo.end();
        for (liter = internals.serlocalinfo.begin(); liter != lguard; liter++) {
            int Q;
            if (liter->uQ == grpm) Q = convert(liter->Q, liter->uQ, l);
            else Q = convert(liter->Q, liter->uQ, mt.getMaterial(liter->tag).density, s);
            dest.addLocalMass(liter->tag, r * Q);
        }
        // D-labelled service materials
        int modsonrod = 0, lastmod = 0;
        for (uint j = 0; j < source.size(); j++) {
            if (modsonrod < source.at(j).getModule().getRing()) {
                modsonrod = source.at(j).getModule().getRing();
                lastmod = j;
            }
        }
        std::vector<SingleSerExit>::const_iterator eiter, eguard = internals.serexitinfo.end();
        for (eiter = internals.serexitinfo.begin(); eiter != eguard; eiter++) {
            int In, Out;
            if (eiter->uIn == grpm) In = convert(eiter->In, eiter->uIn, l);
            else In = convert(eiter->In, eiter->uIn, mt.getMaterial(eiter->tagIn).density, s);
            if (eiter->uOut == grpm) Out = convert(eiter->Out, eiter->uOut, l);
            else Out = convert(eiter->Out, eiter->uOut, mt.getMaterial(eiter->tagOut).density, s);
            Out = (double)Out * source.at(lastmod).getExitingMass(eiter->tagIn) / (double)In;
            if (eiter->is_local) dest.addLocalMass(eiter->tagOut, r * Out);
            else dest.addExitingMass(eiter->tagOut, r * Out);
        }
    }
    
    void MatCalc::adjacentSameCategory(InactiveElement& source, InactiveElement& dest) {
        double tmp;
        for (uint j = 0; j < source.exitingMassCount(); j++) {
            tmp = source.getExitingMass(j);
            if (source.isVertical()) tmp = tmp / source.getRWidth();
            else tmp = tmp / source.getZLength();
            if (dest.isVertical()) tmp = tmp * dest.getRWidth();
            else tmp = tmp * dest.getZLength();
            dest.addExitingMass(source.getExitingTag(j), tmp);
        }
    }
}

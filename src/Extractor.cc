/**
 * @file Extractor.cc
 * @brief This class extracts the tracker geometry and material properties from an existing setup and groups them in preparation for translation to CMSSW XML
 */

#include <Extractor.h>
namespace insur {
    //public
    /**
     * This is the public analysis function that extracts the information that is necessary to convert a given material budget to a
     * series of CMSSW XML files. None of this information is written to file, though. Instead, it is stored in a custom struct that
     * consists of a series of vectors listing different types of information chunks. Those can then be used later to write the specific
     * XML blocks in the output files. The function itself does bookkeeping of the input and output data, not actual analysis. It
     * delegates that to a series of internal custom functions, making sure that each of them gets the correct parameters to run.
     * @param mt A reference to the global material table; used as input
     * @param mb A reference to the material budget that is to be analysed; used as input
     * @param d A reference to the bundle of vectors that will contain the information extracted during analysis; used for output
     */
    void Extractor::analyse(MaterialTable& mt, MaterialBudget& mb, CMSSWBundle& d, bool wt) {
        std::cout << "Starting analysis..." << std::endl;
        Tracker& tr = mb.getTracker();
        InactiveSurfaces& is = mb.getInactiveSurfaces();
        std::vector<std::vector<ModuleCap> >& bc = mb.getBarrelModuleCaps();
        std::vector<std::vector<ModuleCap> >& ec = mb.getEndcapModuleCaps();
        std::vector<Element>& e = d.elements;
        std::vector<Composite>& c = d.composites;
        std::vector<LogicalInfo>& l = d.logic;
        std::vector<ShapeInfo>& s = d.shapes;
        std::vector<PosInfo>& p = d.positions;
        std::vector<AlgoInfo>& a = d.algos;
        std::vector<Rotation>& r = d.rots;
        std::vector<SpecParInfo>& t = d.specs;
        std::vector<RILengthInfo>& ri = d.lrilength;
        //int layer;
        e.clear();
        c.clear();
        l.clear();
        s.clear();
        p.clear();
        a.clear();
        r.clear();
        t.clear();
        ri.clear();
        // inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        AlgoInfo alg;
        Rotation rot;
        SpecParInfo spec;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        pos.trans.dz = 0.0;
        pos.rotref = "";
        // initialise rotation list with Harry's tilt mod
        rot.name = xml_barrel_tilt;
        rot.thetax = 90.0;
        rot.phix = 270.0;
        rot.thetay = 180.0;
        rot.phiy = 0.0;
        rot.thetaz = 90.0;
        rot.phiz = 0.0;
        r.push_back(rot);
        // define top-level barrel volume container (polycone)
        if (!wt) {
            shape.type = pc;
            shape.name_tag = xml_tob;
            analyseBarrelContainer(tr, shape.rzup, shape.rzdown);
            s.push_back(shape);
            std::cout << "Barrel container done." << std::endl;
            // define top-level endcap volume containers (polycone)
            analyseEndcapContainer(tr, shape.rzup, shape.rzdown);
            if (!(shape.rzup.empty() || shape.rzdown.empty())) {
                shape.name_tag = xml_tid;
                s.push_back(shape);
            }
            std::cout << "Endcap container done." << std::endl;
        }
        // translate entries in mt to elementary materials
        analyseElements(mt, e);
        std::cout << "Elementary materials done." << std::endl;
        // analyse barrel
        analyseLayers(mt, bc, tr, c, l, s, p, a, r, t, ri, wt);
        std::cout << "Barrel layers done." << std::endl;
        // analyse endcaps
        analyseDiscs(mt, ec, tr, c, l, s, p, a, r, t, ri, wt);
        std::cout << "Endcap discs done." << std::endl;
        // barrel services
        analyseBarrelServices(is, c, l, s, p, t);
        std::cout << "Barrel services done." << std::endl;
        // endcap services
        analyseEndcapServices(is, c, l, s, p, t);
        std::cout << "Endcap services done." << std::endl;
        // supports
        analyseSupports(is, c, l, s, p, t);
        std::cout << "Support structures done." << std::endl;
        std::cout << "Analysis done." << std::endl;
    }
    
    //protected
    /**
     * This is one of the smaller analysis functions that provide the core functionality of this class. It goes through
     * the global material table, treating each entry as an elementary material and copying or converting the information
     * for that material (density, standard radiation length, standard interaction length) to the properties (density,
     * atomic weight and atomic number) used by CMSSW to describe elementary materials.
     * @param mattab A reference to the global material table; used as input
     * @param elems A reference to the collection of elementary material information; used as output
     */
    void Extractor::analyseElements(MaterialTable&mattab, std::vector<Element>& elems) {
        for (unsigned int i = 0; i < mattab.rowCount(); i++) {
            Element e;
            MaterialRow& r = mattab.getMaterial(i);
            e.tag = r.tag;
            e.density = r.density;
            e.atomic_weight = pow((r.ilength / 35.), 3); // magic!
            e.atomic_number = Z(r.rlength, e.atomic_weight);
            elems.push_back(e);
        }
    }
    
    /**
     * This is one of the smaller analysis functions that provide the core functionality of this class. Its main purpose is that
     * of extracting a series of (r, z) points that will be used later to extend the polycone volume enclosing the entire pixel
     * and tracker barrels. Since those points need to be in order around the enclosing polygon, they are grouped into two
     * different vectors, one for z- and one for z+. Because of the way the points are extracted, mirroring points in z- and
     * z+ will be extracted bottom-to-top and placed in vectors <i>up</i> and <i>down</i>, respectively. The idea is
     * that these two vectors should be traversed in opposite directions later: <i>up</i> from first to last element, <i>down</i>
     * from last to first.
     * @param t A reference to the collection of topology information
     * @param up A reference to a vector listing polygon points by increasing radius
     * @param down A reference to a vector listing polygon points by decreasing radius
     */
    void Extractor::analyseBarrelContainer(Tracker& t, std::vector<std::pair<double, double> >& up,
            std::vector<std::pair<double, double> >& down) {
        bool is_short, previous_short = false;
        std::pair<double, double> rz;
        double rmax = 0.0, zmax = 0.0, zmin = 0.0;
        unsigned int layer, n_of_layers = t.getBarrelLayers()->size();
        std::vector<Layer*>* bl = t.getBarrelLayers();
        up.clear();
        down.clear();
        for (layer = 0; layer < n_of_layers; layer++) {
            is_short = (bl->at(layer)->getMinZ() > 0) || (bl->at(layer)->getMaxZ() < 0);
            if (is_short) {
                //short layer on z- side
                if (bl->at(layer)->getMaxZ() < 0) {
                    //indices 0, 1
                    if ((layer == 0) || ((layer == 1) && (previous_short))) {
                        rz.first = bl->at(layer)->getMinRho();
                        rz.second = bl->at(layer)->getMinZ();
                        up.push_back(rz);
                    }
                    else {
                        //new barrel reached
                        if (bl->at(layer)->getMinZ() != zmin) {
                            //new layer sticks out compared to old layer
                            if (bl->at(layer)->getMinZ() > zmin) rz.first = bl->at(layer)->getMinRho();
                            //old layer sticks out compared to new layer
                            else rz.first = rmax;
                            rz.second = zmin;
                            up.push_back(rz);
                            rz.second = bl->at(layer)->getMinZ();
                            up.push_back(rz);
                        }
                    }
                    //indices size - 2, size - 1
                    if ((layer == n_of_layers - 1) || ((layer == n_of_layers - 2) && (previous_short))) {
                        rz.first = bl->at(layer)->getMaxRho();
                        rz.second = bl->at(layer)->getMinZ();
                        up.push_back(rz);
                    }
                }
                //short layer on z+ side
                else {
                    //indices 0, 1
                    if ((layer == 0) || ((layer == 1) && (previous_short))) {
                        rz.first = bl->at(layer)->getMinRho();
                        rz.second = bl->at(layer)->getMaxZ();
                        down.push_back(rz);
                    }
                    else {
                        //new barrel reached
                        if (bl->at(layer)->getMaxZ() != zmax) {
                            //new layer sticks out compared to old layer
                            if (bl->at(layer)->getMaxZ() > zmax) rz.first = bl->at(layer)->getMinRho();
                            //old layer sticks out compared to new layer
                            else rz.first = rmax;
                            rz.second = zmax;
                            down.push_back(rz);
                            rz.second = bl->at(layer)->getMaxZ();
                            down.push_back(rz);
                        }
                    }
                    //indices size - 2, size - 1
                    if ((layer == n_of_layers - 1) || ((layer == n_of_layers - 2) && (previous_short))) {
                        rz.first = bl->at(layer)->getMaxRho();
                        rz.second = bl->at(layer)->getMinZ();
                        down.push_back(rz);
                    }
                }
            }
            //regular layer across z=0
            else {
                //index 0
                if (layer == 0) {
                    rz.first = bl->at(layer)->getMinRho();
                    rz.second = bl->at(layer)->getMinZ();
                    up.push_back(rz);
                    rz.second = bl->at(layer)->getMaxZ();
                    down.push_back(rz);
                }
                else {
                    //new barrel reached
                    if (bl->at(layer)->getMaxZ() != zmax) {
                        //new layer sticks out compared to old layer
                        if (bl->at(layer)->getMaxZ() > zmax) rz.first = bl->at(layer)->getMinRho();
                        //old layer sticks out compared to new layer
                        else rz.first = rmax;
                        rz.second = zmin;
                        up.push_back(rz);
                        rz.second = zmax;
                        down.push_back(rz);
                        rz.second = bl->at(layer)->getMinZ();
                        up.push_back(rz);
                        rz.second = bl->at(layer)->getMaxZ();
                        down.push_back(rz);
                    }
                }
                //index size - 1
                if (layer == n_of_layers - 1) {
                    rz.first = bl->at(layer)->getMaxRho();
                    rz.second = bl->at(layer)->getMinZ();
                    up.push_back(rz);
                    rz.second = bl->at(layer)->getMaxZ();
                    down.push_back(rz);
                }
            }
            rmax = bl->at(layer)->getMaxRho();
            if (bl->at(layer)->getMinZ() < 0) zmin = bl->at(layer)->getMinZ();
            if (bl->at(layer)->getMaxZ() > 0) zmax = bl->at(layer)->getMaxZ();
            previous_short = is_short;
        }
    }
    
    /**
     *This is one of the smaller analysis functions that provide the core functionality of this class. Its main purpose is that
     * of extracting a series of (r, z) points that will be used later to extend the polycone volume enclosing one of the pixel
     * and tracker endcaps, namely those in z+. Since those points need to be in order around the enclosing polygon, they
     * are grouped into two different vectors, one for for those lying to the left of an imaginary line vertically bisecting the
     * endcaps, the other for those lying to the right of it. Because of the way the points are extracted, mirroring points to
     * the left and right will be extracted bottom-to-top and placed in vectors <i>up</i> and <i>down</i>, respectively.
     * The idea is that these two vectors should be traversed in opposite directions later: <i>up</i> from first to last element,
     * <i>down</i> from last to first.
     * @param t A reference to the collection of topology information
     * @param up A reference to a vector listing polygon points by increasing radius
     * @param down A reference to a vector listing polygon points by decreasing radius
     */
    void Extractor::analyseEndcapContainer(Tracker& t,
            std::vector<std::pair<double, double> >& up, std::vector<std::pair<double, double> >& down) {
        int first, last;
        std::pair<double, double> rz;
        double rmin = 0.0, rmax = 0.0, zmax = 0.0;
        std::vector<Layer*>* el = t.getEndcapLayers();
        first = 0;
        last = el->size();
        while (first < last) {
            if (el->at(first)->getMaxZ() > 0) break;
            first++;
        }
        up.clear();
        down.clear();
        for (int i = first; i < last; i++) {
            // special treatment for first disc
            if (i == first) {
                rmin = el->at(i)->getMinRho();
                rmax = el->at(i)->getMaxRho();
                rz.first = rmax;
                rz.second = el->at(i)->getMinZ() - xml_z_pixfwd;
                up.push_back(rz);
                rz.first = rmin;
                down.push_back(rz);
            }
            // disc beyond the first
            else {
                // endcap change larger->smaller
                if (rmax > el->at(i)->getMaxRho()) {
                    rz.second = zmax - xml_z_pixfwd;
                    rz.first = rmax;
                    up.push_back(rz);
                    rz.first = rmin;
                    down.push_back(rz);
                    rmax = el->at(i)->getMaxRho();
                    rmin = el->at(i)->getMinRho();
                    rz.first = rmax;
                    up.push_back(rz);
                    rz.first = rmin;
                    down.push_back(rz);
                }
                // endcap change smaller->larger
                if (rmax < el->at(i)->getMaxRho()) {
                    rz.second = el->at(i)->getMinZ() - xml_z_pixfwd;
                    rz.first = rmax;
                    up.push_back(rz);
                    rz.first = rmin;
                    down.push_back(rz);
                    rmax = el->at(i)->getMaxRho();
                    rmin = el->at(i)->getMinRho();
                    rz.first = rmax;
                    up.push_back(rz);
                    rz.first = rmin;
                    down.push_back(rz);
                }
            }
            zmax = el->at(i)->getMaxZ();
            // special treatment for last disc
            if (i == last - 1) {
                rz.first = rmax;
                rz.second = zmax - xml_z_pixfwd;
                up.push_back(rz);
                rz.first = rmin;
                down.push_back(rz);
            }
        }
    }
    
    /**
     * This is one of the two main analysis functions that provide the core functionality of this class. It examines the barrel layers
     * and the modules within, extracting a great range of different pieces of information from the geometry layout. These are shapes
     * for individual modules, but also for their enclosing volumes, divided into rods and the layers themselves. They form hierarchies
     * of volumes, one inside the other, and are placed within their parent volumes according to individual placement rules or algorithms,
     *  sometimes using globally defined rotations. They also include some topology, such as which volumes contain the active surfaces,
     * and how those active surfaces are subdivided and connected to the readout electronics. Last but not least, overall radiation and
     * interaction lengths for each layer are calculated and stored; those are used as approximative values for certain CMSSW functions
     * later on. Short layers are treated according to some specialised rules, but also processed in here.
     * @param mt A reference to the global material table; used as input
     * @param bc A reference to the collection of material properties of the barrel modules; used as input
     * @param tr A reference to the tracker object; used as input
     * @param c A reference to the collection of composite material information; used for output
     * @param l A reference to the collection of volume hierarchy information; used for output
     * @param s A reference to the collection of shape parameters; used for output
     * @param p A reference to the collection of volume positionings; used for output
     * @param a A reference to the collection of algorithm calls and their parameters; used for output
     * @param r A reference to the collection of rotations; used for output
     * @param t A reference to the collection of topology information; used for output
     * @param ri A reference to the collection of overall radiation and interaction lengths per layer or disc; used for output
     */
    void Extractor::analyseLayers(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& bc, Tracker& tr,
            std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<PosInfo>& p,
            std::vector<AlgoInfo>& a, std::vector<Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, bool wt) {
        int layer;
        std::string nspace;
        if (wt) nspace = xml_newfileident;
        else nspace = xml_fileident;
        std::vector<std::vector<ModuleCap> >::iterator oiter, oguard;
        std::vector<ModuleCap>::iterator iiter, iguard;
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        AlgoInfo alg;
        Rotation rot;
        SpecParInfo rocdims, lspec, rspec, mspec;
        RILengthInfo ril;
        shape.dyy = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        pos.trans.dz = 0.0;
        rot.phix = 0.0;
        rot.phiy = 0.0;
        rot.phiz = 0.0;
        rot.thetax = 0.0;
        rot.thetay = 0.0;
        rot.thetaz = 0.0;
        lspec.name = xml_subdet_layer + xml_par_tail;
        lspec.parameter.first = xml_tkddd_structure;
        lspec.parameter.second = xml_det_layer;
        rspec.name = xml_subdet_rod + xml_par_tail;
        rspec.parameter.first = xml_tkddd_structure;
        rspec.parameter.second = xml_det_rod;
        mspec.name = xml_subdet_tobdet + xml_par_tail;
        mspec.parameter.first = xml_tkddd_structure;
        mspec.parameter.second = xml_det_tobdet;
        ril.barrel = true;
        ril.index = 0;
        // b_mod: one composite for every module position on rod
        // s and l: one entry for every module position on rod (box), one for every layer (tube), rods TBD
        // p: one entry for every layer (two for short layers), two modules, one wafer and active for each ring on rod
        // a: rods within layer (twice in case of a short layer)
        layer = 1;
        alg.name = xml_tobalgo;
        oguard = bc.end();
        // barrel caps layer loop
        for (oiter = bc.begin(); oiter != oguard; oiter++) {
            double rmin = tr.getBarrelLayers()->at(layer - 1)->getMinRho();
            rmin = rmin - tr.getBarrelLayers()->at(layer - 1)->getMaxModuleThickness() / 2.0;
            double rmax = tr.getBarrelLayers()->at(layer - 1)->getMaxRho();
            rmax = rmax + tr.getBarrelLayers()->at(layer - 1)->getMaxModuleThickness() / 2.0;
            double zmin = tr.getBarrelLayers()->at(layer - 1)->getMinZ();
            double zmax = tr.getBarrelLayers()->at(layer - 1)->getMaxZ();
            double deltar = findDeltaR(tr.getBarrelLayers()->at(layer - 1)->getModuleVector()->begin(),
                    tr.getBarrelLayers()->at(layer - 1)->getModuleVector()->end(), (rmin + rmax) / 2.0);
            double ds, dt = 0.0;
            double rtotal = 0.0, itotal = 0.0;
            int count = 0;
            if (deltar == 0.0) continue;
            bool is_short = (zmax < 0.0) || (zmin > 0.0);
            bool is_relevant = !is_short || (zmin > 0.0);
            if (is_relevant) {
                shape.type = bx;
                shape.rmin = 0.0;
                shape.rmax = 0.0;
                ril.index = layer;
                std::set<int> rings;
                std::ostringstream lname, rname, pconverter;
                lname << xml_layer << layer;
                rname << xml_rod << layer;
                iguard = oiter->end();
                // module caps loop
                for (iiter = oiter->begin(); iiter != iguard; iiter++) {
                    if (rings.find(iiter->getModule().getRing()) == rings.end()) {
                        std::vector<ModuleCap>::iterator partner;
                        std::ostringstream matname, shapename, specname;
                        // module composite material
                        matname << xml_base_actcomp << "L" << layer << "P" << iiter->getModule().getRing();
                        c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));
                        // module box
                        shapename << iiter->getModule().getRing() << lname.str();
                        shape.name_tag = xml_barrel_module + shapename.str();
                        shape.dx = iiter->getModule().getArea() / iiter->getModule().getHeight() / 2.0;
                        shape.dy = iiter->getModule().getHeight() / 2.0;
                        shape.dz = iiter->getModule().getModuleThickness() / 2.0;
                        s.push_back(shape);
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = nspace + ":" + logic.name_tag;
                        logic.material_tag = nspace + ":" + matname.str();
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        pos.rotref = nspace + ":" + xml_barrel_tilt;
                        if ((iiter->getModule().getMeanPoint().Rho() > (rmax - deltar / 2.0))
                                || ((iiter->getModule().getMeanPoint().Rho() < ((rmin + rmax) / 2.0))
                                && (iiter->getModule().getMeanPoint().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dz;
                        else pos.trans.dx = shape.dz - deltar / 2.0;
                        if (is_short) {
                            pos.parent_tag = nspace + ":" + rname.str() + xml_plus;
                            pos.trans.dz = iiter->getModule().getMinZ() - ((zmax + zmin) / 2.0) + shape.dy;
                            p.push_back(pos);
                            pos.parent_tag = nspace + ":" + rname.str() + xml_minus;
                            pos.trans.dz = -pos.trans.dz;
                            pos.copy = 2;
                            p.push_back(pos);
                            pos.copy = 1;
                        }
                        else {
                            pos.parent_tag = nspace + ":" + rname.str();
                            partner = findPartnerModule(iiter, iguard, iiter->getModule().getRing());
                            if (iiter->getModule().getMeanPoint().Z() > 0) {
                                pos.trans.dz = iiter->getModule().getMaxZ() - shape.dy;
                                p.push_back(pos);
                                if (partner != iguard) {
                                    if ((partner->getModule().getMeanPoint().Rho() > (rmax - deltar / 2.0))
                                            || ((partner->getModule().getMeanPoint().Rho() < ((rmin + rmax) / 2.0))
                                            && (partner->getModule().getMeanPoint().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dz;
                                    else pos.trans.dx = shape.dz - deltar / 2.0;
                                    pos.trans.dz = partner->getModule().getMaxZ() - shape.dy;
                                    pos.copy = 2;
                                    p.push_back(pos);
                                    pos.copy = 1;
                                }
                            }
                            else {
                                pos.trans.dz = iiter->getModule().getMaxZ() - shape.dy;
                                pos.copy = 2;
                                p.push_back(pos);
                                pos.copy = 1;
                                if (partner != iguard) {
                                    if ((partner->getModule().getMeanPoint().Rho() > (rmax - deltar / 2.0))
                                            || ((partner->getModule().getMeanPoint().Rho() < ((rmin + rmax) / 2.0))
                                            && (partner->getModule().getMeanPoint().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dz;
                                    else pos.trans.dx = shape.dz - deltar / 2.0;
                                    pos.trans.dz = partner->getModule().getMaxZ() - shape.dy;
                                    p.push_back(pos);
                                }
                            }
                        }
                        pos.rotref = "";
                        // wafer
                        shape.name_tag = xml_barrel_module + shapename.str() + xml_base_waf;
                        shape.dz = calculateSensorThickness(*iiter, mt) / 2.0;
                        if (iiter->getModule().getNFaces() == 2) shape.dz = shape.dz / 2.0;
                        s.push_back(shape);
                        pos.parent_tag = logic.shape_tag;
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = nspace + ":" + logic.name_tag;
                        logic.material_tag = xml_material_air;
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        pos.trans.dx = 0.0;
                        pos.trans.dz = shape.dz - iiter->getModule().getModuleThickness() / 2.0;
                        p.push_back(pos);
                        if (iiter->getModule().getNFaces() == 2) {
                            pos.trans.dz = pos.trans.dz + 2 * shape.dz + iiter->getModule().getStereoDistance();
                            pos.copy = 2;
                            if (iiter->getModule().getStereoRotation() != 0) {
                                rot.name = type_stereo + xml_barrel_module + shapename.str();
                                rot.thetax = 90.0;
                                rot.phix = iiter->getModule().getStereoRotation() / M_PI * 180;
                                rot.thetay = 90.0;
                                rot.phiy = 90.0 + iiter->getModule().getStereoRotation() / M_PI * 180;
                                r.push_back(rot);
                                pos.rotref = nspace + ":" + rot.name;
                            }
                            p.push_back(pos);
                            pos.rotref.clear();
                            rot.name.clear();
                            rot.thetax = 0.0;
                            rot.phix = 0.0;
                            rot.thetay = 0.0;
                            rot.phiy = 0.0;
                            pos.copy = 1;
                        }
                        // active surface
                        shape.name_tag = xml_barrel_module + shapename.str() + xml_base_act;
                        s.push_back(shape);
                        pos.parent_tag = logic.shape_tag;
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = nspace + ":" + logic.name_tag;
                        logic.material_tag = nspace + ":" + xml_sensor_silicon;
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        pos.trans.dz = 0.0;
                        p.push_back(pos);
                        // topology
                        mspec.partselectors.push_back(logic.name_tag);
                        specname << xml_roc_x << xml_par_tail << (iiter->getModule().getNStripsAcross() / xml_roc_rows);
                        int id = findSpecParIndex(t, specname.str());
                        if (id >= 0) t.at(id).partselectors.push_back(logic.name_tag);
                        else {
                            rocdims.partselectors.clear();
                            rocdims.name = specname.str();
                            rocdims.parameter.first = xml_roc_x;
                            specname.str("");
                            specname << (iiter->getModule().getNStripsAcross() / xml_roc_rows);
                            rocdims.parameter.second = specname.str();
                            rocdims.partselectors.push_back(logic.name_tag);
                            t.push_back(rocdims);
                        }
                        specname.str("");
                        specname << xml_roc_y << xml_par_tail << iiter->getModule().getNSegments();
                        id = findSpecParIndex(t, specname.str());
                        if (id >= 0) t.at(id).partselectors.push_back(logic.name_tag);
                        else {
                            rocdims.partselectors.clear();
                            rocdims.name = specname.str();
                            rocdims.parameter.first = xml_roc_y;
                            specname.str("");
                            specname << iiter->getModule().getNSegments();
                            rocdims.parameter.second = specname.str();
                            rocdims.partselectors.push_back(logic.name_tag);
                            t.push_back(rocdims);
                        }
                        // material properties
                        rtotal = rtotal + iiter->getRadiationLength();
                        itotal = itotal + iiter->getInteractionLength();
                        count++;
                        rings.insert(iiter->getModule().getRing());
                        dt = iiter->getModule().getModuleThickness();
                    }
                }
                if (count > 0) {
                    ril.rlength = rtotal / (double)count;
                    ril.ilength = itotal / (double)count;
                    ri.push_back(ril);
                }
                // rod(s)
                shape.name_tag = rname.str();
                if (is_short) shape.name_tag = shape.name_tag + xml_plus;
                shape.dx = deltar / 2.0;
                if (is_short) shape.dz = (zmax - zmin) / 2.0;
                else shape.dz = zmax;
                s.push_back(shape);
                logic.name_tag = shape.name_tag;
                logic.shape_tag = nspace + ":" + logic.name_tag;
                logic.material_tag = xml_material_air;
                l.push_back(logic);
                rspec.partselectors.push_back(logic.name_tag);
                pconverter << logic.shape_tag;
                if (is_short) {
                    shape.name_tag = rname.str() + xml_minus;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = nspace + ":" + logic.name_tag;
                    l.push_back(logic);
                    rspec.partselectors.push_back(logic.name_tag);
                }
                ds = fromRim(rmax, shape.dy);
                // extra containers for short layers
                if (wt && is_short) {
                    shape.type = tb;
                    shape.dx = 0.0;
                    shape.dy = 0.0;
                    shape.rmin = rmin;
                    shape.rmax = rmax;
                    shape.name_tag = lname.str() + xml_plus;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = nspace + ":" + logic.name_tag;
                    l.push_back(logic);
                    pos.trans.dx = 0.0;
                    pos.trans.dz = zmin + (zmax - zmin) / 2.0;
                    pos.parent_tag = nspace + ":" + lname.str();
                    pos.child_tag = logic.shape_tag;
                    p.push_back(pos);
                    shape.name_tag = lname.str() + xml_minus;
                    s.push_back(shape);
                    logic.name_tag = shape.name_tag;
                    logic.shape_tag = nspace + ":" + logic.name_tag;
                    l.push_back(logic);
                    pos.trans.dz = -pos.trans.dz;
                    pos.child_tag = logic.shape_tag;
                    p.push_back(pos);
                }
                // layer
                shape.type = tb;
                shape.dx = 0.0;
                shape.dy = 0.0;
                pos.trans.dx = 0.0;
                pos.trans.dz = 0.0;
                shape.name_tag = lname.str();
                shape.rmin = rmin;
                shape.rmax = rmax;
                shape.dz = zmax;
                s.push_back(shape);
                logic.name_tag = shape.name_tag;
                logic.shape_tag = nspace + ":" + logic.name_tag;
                l.push_back(logic);
                pos.parent_tag = xml_pixbarident + ":" + xml_pixbar;
                pos.child_tag = logic.shape_tag;
                p.push_back(pos);
                lspec.partselectors.push_back(logic.name_tag);
                // rods in layer algorithm(s)
                alg.parent = logic.shape_tag;
                if (wt && is_short) alg.parent = alg.parent + xml_plus;
                alg.parameters.push_back(stringParam(xml_childparam, pconverter.str()));
                pconverter.str("");
                pconverter << (tr.getBarrelLayers()->at(layer - 1)->getTilt() + 90) << "*deg";
                alg.parameters.push_back(numericParam(xml_tilt, pconverter.str()));
                pconverter.str("");
                pconverter << tr.getBarrelLayers()->at(layer - 1)->getStartAngle();
                alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
                pconverter.str("");
                alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
                pconverter << (rmin + deltar / 2.0) << "*mm";
                alg.parameters.push_back(numericParam(xml_radiusin, pconverter.str()));
                pconverter.str("");
                pconverter << (rmax - ds - deltar / 2.0 - 2.0 * dt) << "*mm";
                alg.parameters.push_back(numericParam(xml_radiusout, pconverter.str()));
                pconverter.str("");
                if (!wt && is_short) {
                    pconverter << (zmin + (zmax - zmin) / 2.0) << "*mm";
                    alg.parameters.push_back(numericParam(xml_zposition, pconverter.str()));
                    pconverter.str("");
                }
                else alg.parameters.push_back(numericParam(xml_zposition, "0.0*mm"));
                pconverter << static_cast<BarrelLayer*>(tr.getBarrelLayers()->at(layer - 1))->getRods();
                alg.parameters.push_back(numericParam(xml_number, pconverter.str()));
                alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
                alg.parameters.push_back(numericParam(xml_incrcopyno, "1"));
                a.push_back(alg);
                // extras for short layers
                if (is_short) {
                    pconverter.str("");
                    if (wt) alg.parent = logic.shape_tag + xml_minus;
                    pconverter << nspace << ":" << rname.str() << xml_minus;
                    alg.parameters.front() = stringParam(xml_childparam, pconverter.str());
                    pconverter.str("");
                    if (wt) pconverter << "0.0*mm";
                    else pconverter << -(zmin + (zmax - zmin) / 2.0) << "*mm";
                    alg.parameters.at(6) = numericParam(xml_zposition, pconverter.str());
                    a.push_back(alg);
                }
                alg.parameters.clear();
            }
            layer++;
        }
        if (!lspec.partselectors.empty()) t.push_back(lspec);
        if (!rspec.partselectors.empty()) t.push_back(rspec);
        if (!mspec.partselectors.empty()) t.push_back(mspec);
    }
    
    /**
     * This is one of the two main analysis functions that provide the core functionality of this class. It examines the endcap discs in z+
     * and the rings and modules within, extracting a great range of different pieces of information from the geometry layout. These
     * are shapes for individual modules, but also for their enclosing volumes, divided into rings and then discs. They form hierarchies
     * of volumes, one inside the other, and are placed within their parent volumes according to individual placement rules or algorithms,
     *  sometimes using globally defined rotations. They also include some topology, such as which volumes contain the active surfaces,
     * and how those active surfaces are subdivided and connected to the readout electronics. Last but not least, overall radiation and
     * interaction lengths for each disc are calculated and stored; those are used as approximative values for certain CMSSW functions
     * later on.
     * @param mt A reference to the global material table; used as input
     * @param ec A reference to the collection of material properties of the endcap modules; used as input
     * @param tr A reference to the tracker object; used as input
     * @param c A reference to the collection of composite material information; used for output
     * @param l A reference to the collection of volume hierarchy information; used for output
     * @param s A reference to the collection of shape parameters; used for output
     * @param p A reference to the collection of volume positionings; used for output
     * @param a A reference to the collection of algorithm calls and their parameters; used for output
     * @param r A reference to the collection of rotations; used for output
     * @param t A reference to the collection of topology information; used for output
     * @param ri A reference to the collection of overall radiation and interaction lengths per layer or disc; used for output
     */
    void Extractor::analyseDiscs(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& ec, Tracker& tr,
            std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<PosInfo>& p,
            std::vector<AlgoInfo>& a, std::vector<Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, bool wt) {
        int layer;
        std::string nspace;
        if (wt) nspace = xml_newfileident;
        else nspace = xml_fileident;
        std::vector<std::vector<ModuleCap> >::iterator oiter, oguard;
        std::vector<ModuleCap>::iterator iiter, iguard;
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        AlgoInfo alg;
        Rotation rot;
        SpecParInfo rocdims, dspec, rspec, mspec;
        RILengthInfo ril;
        shape.dyy = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        pos.trans.dz = 0.0;
        rot.phix = 0.0;
        rot.phiy = 0.0;
        rot.phiz = 0.0;
        rot.thetax = 0.0;
        rot.thetay = 0.0;
        rot.thetaz = 0.0;
        dspec.name = xml_subdet_wheel + xml_par_tail;
        dspec.parameter.first = xml_tkddd_structure;
        dspec.parameter.second = xml_det_wheel;
        rspec.name = xml_subdet_ring + xml_par_tail;
        rspec.parameter.first = xml_tkddd_structure;
        rspec.parameter.second = xml_det_ring;
        mspec.name = xml_subdet_tiddet + xml_par_tail;
        mspec.parameter.first = xml_tkddd_structure;
        mspec.parameter.second = xml_det_tiddet;
        ril.barrel = false;
        ril.index = 0;
        // e_mod: one composite for every ring
        // s and l: one entry for every ring module, one for every ring, one for every disc
        // p: one entry for every disc, one for every ring, one module, wafer and active per ring
        // a: two per ring with modules inside ring
        layer = 1;
        alg.name = xml_ecalgo;
        oguard = ec.end();
        // endcap caps layer loop
        for (oiter = ec.begin(); oiter != oguard; oiter++) {
            if (tr.getEndcapLayers()->at(layer - 1)->getMinZ() > 0) {
                ril.index = layer;
                std::set<int> ridx;
                std::map<int, RingInfo> rinfo;
                double rmin = tr.getEndcapLayers()->at(layer - 1)->getMinRho();
                double rmax = tr.getEndcapLayers()->at(layer - 1)->getMaxRho();
                double zmax = tr.getEndcapLayers()->at(layer - 1)->getMaxZ();
                zmax = zmax + tr.getEndcapLayers()->at(layer - 1)->getMaxModuleThickness() / 2.0;
                double zmin = tr.getEndcapLayers()->at(layer - 1)->getMinZ();
                zmin = zmin - tr.getEndcapLayers()->at(layer - 1)->getMaxModuleThickness() / 2.0;
                std::ostringstream dname, pconverter;
                double rtotal = 0.0, itotal = 0.0;
                int count = 0;
                dname << xml_disc << layer;
                shape.type = tp;
                shape.rmin = 0.0;
                shape.rmax = 0.0;
                pos.trans.dz = 0.0;
                iguard = oiter->end();
                // endcap module caps loop
                for (iiter = oiter->begin(); iiter != iguard; iiter++) {
                    // new ring
                    if (ridx.find(iiter->getModule().getRing()) == ridx.end()) {
                        ridx.insert(iiter->getModule().getRing());
                        std::ostringstream matname, rname, mname, specname;
                        matname << xml_base_actcomp << "D" << layer << "R" << iiter->getModule().getRing();
                        c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));
                        rname << xml_ring << iiter->getModule().getRing() << dname.str();
                        mname << xml_endcap_module << iiter->getModule().getRing() << dname.str();
                        // collect ring info
                        RingInfo rinf;
                        rinf.name = rname.str();
                        rinf.childname = mname.str();
                        rinf.fw = (iiter->getModule().getMeanPoint().Z() < (zmin + zmax) / 2.0);
                        rinf.modules = static_cast<EndcapLayer*>(tr.getEndcapLayers()->at(layer - 1))->getModulesOnRing().at(iiter->getModule().getRing() - 1);
                        rinf.rin = iiter->getModule().getMinRho();
                        rinf.rout = iiter->getModule().getMaxRho();
                        rinf.rmid = iiter->getModule().getMeanPoint().Rho();
                        rinf.mthk = iiter->getModule().getModuleThickness();
                        rinf.phi = iiter->getModule().getMeanPoint().Phi();
                        rinfo.insert(std::pair<int, RingInfo>(iiter->getModule().getRing(), rinf));
                        // module trapezoid
                        shape.name_tag = mname.str();
                        shape.dx = iiter->getModule().getHeight() / 2.0;
                        shape.dy = static_cast<EndcapModule&>(iiter->getModule()).getWidthLo() / 2.0;
                        shape.dyy = static_cast<EndcapModule&>(iiter->getModule()).getWidthHi() / 2.0;
                        shape.dz = iiter->getModule().getModuleThickness() / 2.0;
                        s.push_back(shape);
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = nspace + ":" + logic.name_tag;
                        logic.material_tag = nspace + ":" + matname.str();
                        l.push_back(logic);
                        // wafer
                        pos.parent_tag = logic.shape_tag;
                        shape.name_tag = mname.str() + xml_base_waf;
                        shape.dz = calculateSensorThickness(*iiter, mt) / 2.0;
                        if (iiter->getModule().getNFaces() == 2) shape.dz = shape.dz / 2.0;
                        s.push_back(shape);
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = nspace + ":" + logic.name_tag;
                        logic.material_tag = xml_material_air;
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        if (iiter->getModule().getMaxZ() > 0) pos.trans.dz = shape.dz - iiter->getModule().getModuleThickness() / 2.0;
                        else pos.trans.dz = iiter->getModule().getModuleThickness() / 2.0 - shape.dz;
                        p.push_back(pos);
                        if (iiter->getModule().getNFaces() == 2) {
                            if (iiter->getModule().getMaxZ() > 0) pos.trans.dz = pos.trans.dz + 2 * shape.dz + iiter->getModule().getStereoDistance();
                            else pos.trans.dz = pos.trans.dz - 2 * shape.dz - iiter->getModule().getStereoDistance();
                            pos.copy = 2;
                            if (iiter->getModule().getStereoRotation() != 0) {
                                rot.name = type_stereo + xml_endcap_module + mname.str();
                                rot.thetax = 90.0;
                                rot.phix = iiter->getModule().getStereoRotation() / M_PI * 180;
                                rot.thetay = 90.0;
                                rot.phiy = 90.0 + iiter->getModule().getStereoRotation() / M_PI * 180;
                                r.push_back(rot);
                                pos.rotref = nspace + ":" + rot.name;
                            }
                            p.push_back(pos);
                            pos.rotref.clear();
                            rot.name.clear();
                            rot.thetax = 0.0;
                            rot.phix = 0.0;
                            rot.thetay = 0.0;
                            rot.phiy = 0.0;
                            pos.copy = 1;
                        }
                        // active surface
                        pos.parent_tag = logic.shape_tag;
                        shape.name_tag = mname.str() + xml_base_act;
                        s.push_back(shape);
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = nspace + ":" + logic.name_tag;
                        logic.material_tag = nspace + ":" + xml_sensor_silicon;
                        l.push_back(logic);
                        pos.child_tag = logic.shape_tag;
                        pos.trans.dz = 0.0;
                        p.push_back(pos);
                        // topology
                        mspec.partselectors.push_back(logic.name_tag);
                        specname << xml_roc_x << xml_par_tail << (iiter->getModule().getNStripsAcross() / xml_roc_rows);
                        int id = findSpecParIndex(t, specname.str());
                        if (id >= 0) t.at(id).partselectors.push_back(logic.name_tag);
                        else {
                            rocdims.partselectors.clear();
                            rocdims.name = specname.str();
                            rocdims.parameter.first = xml_roc_x;
                            specname.str("");
                            specname << (iiter->getModule().getNStripsAcross() / xml_roc_rows);
                            rocdims.parameter.second = specname.str();
                            rocdims.partselectors.push_back(logic.name_tag);
                            t.push_back(rocdims);
                        }
                        specname.str("");
                        specname << xml_roc_y << xml_par_tail << iiter->getModule().getNSegments();
                        id = findSpecParIndex(t, specname.str());
                        if (id >= 0) t.at(id).partselectors.push_back(logic.name_tag);
                        else {
                            rocdims.partselectors.clear();
                            rocdims.name = specname.str();
                            rocdims.parameter.first = xml_roc_y;
                            specname.str("");
                            specname << iiter->getModule().getNSegments();
                            rocdims.parameter.second = specname.str();
                            rocdims.partselectors.push_back(logic.name_tag);
                            t.push_back(rocdims);
                        }
                        // material properties
                        rtotal = rtotal + iiter->getRadiationLength();
                        itotal = itotal + iiter->getInteractionLength();
                        count++;
                    }
                }
                if (count > 0) {
                    ril.rlength = rtotal / (double)count;
                    ril.ilength = itotal / (double)count;
                    ri.push_back(ril);
                }
                // rings
                shape.type = tb;
                shape.dx = 0.0;
                shape.dy = 0.0;
                shape.dyy = 0.0;
                shape.dz = findDeltaZ(tr.getEndcapLayers()->at(layer - 1)->getModuleVector()->begin(),
                        tr.getEndcapLayers()->at(layer - 1)->getModuleVector()->end(), (zmin + zmax) / 2.0) / 2.0;
                std::set<int>::const_iterator siter, sguard = ridx.end();
                for (siter = ridx.begin(); siter != sguard; siter++) {
                    if (rinfo[*siter].modules > 0) {
                        shape.name_tag = rinfo[*siter].name;
                        shape.rmin = rinfo[*siter].rin;
                        shape.rmax = rinfo[*siter].rout;
                        s.push_back(shape);
                        logic.name_tag = shape.name_tag;
                        logic.shape_tag = nspace + ":" + logic.name_tag;
                        logic.material_tag = xml_material_air;
                        l.push_back(logic);
                        pos.parent_tag = nspace + ":" + dname.str() + xml_plus;
                        pos.child_tag = logic.shape_tag;
                        if (rinfo[*siter].fw) pos.trans.dz = (zmin - zmax) / 2.0 + shape.dz;
                        else pos.trans.dz = (zmax - zmin) / 2.0 - shape.dz;
                        p.push_back(pos);
                        pos.parent_tag = nspace + ":" + dname.str() + xml_minus;
                        p.push_back(pos);
                        rspec.partselectors.push_back(logic.name_tag);
                        alg.parent = logic.shape_tag;
                        alg.parameters.push_back(stringParam(xml_childparam, nspace + ":" + rinfo[*siter].childname));
                        pconverter << (rinfo[*siter].modules / 2);
                        alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
                        pconverter.str("");
                        alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
                        alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
                        alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
                        pconverter << rinfo[*siter].phi;
                        alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
                        pconverter.str("");
                        pconverter << rinfo[*siter].rmid;
                        alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
                        pconverter.str("");
                        alg.parameters.push_back(vectorParam(0, 0, shape.dz - rinfo[*siter].mthk / 2.0));
                        a.push_back(alg);
                        alg.parameters.clear();
                        alg.parameters.push_back(stringParam(xml_childparam, nspace + ":" + rinfo[*siter].childname));
                        pconverter << (rinfo[*siter].modules / 2);
                        alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
                        pconverter.str("");
                        alg.parameters.push_back(numericParam(xml_startcopyno, "2"));
                        alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
                        alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
                        pconverter << (rinfo[*siter].phi + 2 * PI / (double)(rinfo[*siter].modules));
                        alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
                        pconverter.str("");
                        pconverter << rinfo[*siter].rmid;
                        alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
                        pconverter.str("");
                        alg.parameters.push_back(vectorParam(0, 0, rinfo[*siter].mthk / 2.0 - shape.dz));
                        a.push_back(alg);
                        alg.parameters.clear();
                    }
                }
                //disc
                shape.name_tag = dname.str();
                shape.rmin = rmin;
                shape.rmax = rmax;
                shape.dz = (zmax - zmin) / 2.0;
                s.push_back(shape);
                logic.name_tag = shape.name_tag + xml_plus;
                logic.shape_tag = nspace + ":" + shape.name_tag;
                logic.material_tag = xml_material_air;
                l.push_back(logic);
                pos.parent_tag = xml_pixfwdident + ":" + xml_pixfwd_plus;
                pos.child_tag = nspace + ":" + logic.name_tag;
                pos.trans.dz = (zmax + zmin) / 2.0 - xml_z_pixfwd;
                p.push_back(pos);
                dspec.partselectors.push_back(logic.name_tag);
                logic.name_tag = shape.name_tag + xml_minus;
                l.push_back(logic);
                pos.parent_tag = xml_pixfwdident + ":" + xml_pixfwd_minus;
                pos.child_tag = nspace + ":" + logic.name_tag;
                p.push_back(pos);
                dspec.partselectors.push_back(logic.name_tag);
            }
            layer++;
        }
        if (!dspec.partselectors.empty()) t.push_back(dspec);
        if (!rspec.partselectors.empty()) t.push_back(rspec);
        if (!mspec.partselectors.empty()) t.push_back(mspec);
    }
    
    /**
     * This is one of the smaller analysis functions that provide the core functionality of this class. It does a number of things:
     * it creates a composite material information struct for each barrel service, and it adds the remaining information about
     * the volume to the collections of hierarchy, shape, position and topology information. All of these future CMSSW XML
     * blocks are given unique names based on the properties of the barrel service they came from.
     * @param is A reference to the collection of inactive surfaces that the barrel services are a part of; used as input
     * @param c A reference to the collection of composite material information; used for output
     * @param l A reference to the collection of volume hierarchy information; used for output
     * @param s A reference to the collection of shape parameters; used for output
     * @param p A reference to the collection of volume positionings; used for output
     * @param t A reference to the collection of topology information; used for output
     */
    void Extractor::analyseBarrelServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
            std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt) {
        std::string nspace;
        if (wt) nspace = xml_newfileident;
        else nspace = xml_fileident;
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        shape.type = tb;
        shape.dx = 0.0;
        shape.dy = 0.0;
        shape.dyy = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        // b_ser: one composite for every service volume on the z+ side
        // s, l and p: one entry per service volume
        std::vector<InactiveElement>::iterator iter, guard;
        std::vector<InactiveElement>& bs = is.getBarrelServices();
        guard = bs.end();
        for (iter = bs.begin(); iter != guard; iter++) {
            std::ostringstream matname, shapename;
            matname << xml_base_serfcomp << iter->getCategory() << "R" << (int)(iter->getInnerRadius()) << "dZ" << (int)(iter->getZLength());
            shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(iter->getZOffset());
            if ((iter->getZOffset() + iter->getZLength()) > 0) c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
            shape.name_tag = shapename.str();
            shape.dz = iter->getZLength() / 2.0;
            shape.rmin = iter->getInnerRadius();
            shape.rmax = shape.rmin + iter->getRWidth();
            s.push_back(shape);
            logic.name_tag = shapename.str();
            logic.shape_tag = nspace + ":" + shapename.str();
            logic.material_tag = nspace + ":" + matname.str();
            l.push_back(logic);
            pos.parent_tag = nspace + ":" + xml_tracker;
            pos.child_tag = logic.shape_tag;
            pos.trans.dz = iter->getZOffset() + shape.dz;
            p.push_back(pos);
        }
    }
    
    /**
     * This is one of the smaller analysis functions that provide the core functionality of this class. It does a number of things:
     * it creates a composite material information struct for each endcap service, and it adds the remaining information about
     * the volume to the collections of hierarchy, shape, position and topology information. All of these future CMSSW XML
     * blocks are given unique names based on the properties of the encap service they came from.
     * @param is A reference to the collection of inactive surfaces that the endcap services are a part of; used as input
     * @param c A reference to the collection of composite material information; used for output
     * @param l A reference to the collection of volume hierarchy information; used for output
     * @param s A reference to the collection of shape parameters; used for output
     * @param p A reference to the collection of volume positionings; used for output
     * @param t A reference to the collection of topology information; used for output
     */
    void Extractor::analyseEndcapServices(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
            std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt) {
        std::string nspace;
        if (wt) nspace = xml_newfileident;
        else nspace = xml_fileident;
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        shape.type = tb;
        shape.dx = 0.0;
        shape.dy = 0.0;
        shape.dyy = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        // e_ser: one composite for every service volume on the z+ side
        // s, l and p: one entry per service volume
        std::vector<InactiveElement>::iterator iter, guard;
        std::vector<InactiveElement>& es = is.getEndcapServices();
        guard = es.end();
        for (iter = es.begin(); iter != guard; iter++) {
            std::ostringstream matname, shapename;
            matname << xml_base_serfcomp << iter->getCategory() << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
            shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
            if ((iter->getZOffset() + iter->getZLength()) > 0) c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
            shape.name_tag = shapename.str();
            shape.dz = iter->getZLength() / 2.0;
            shape.rmin = iter->getInnerRadius();
            shape.rmax = shape.rmin + iter->getRWidth();
            s.push_back(shape);
            logic.name_tag = shapename.str();
            logic.shape_tag = nspace + ":" + shapename.str();
            logic.material_tag = nspace + ":" + matname.str();
            l.push_back(logic);
            pos.parent_tag = nspace + ":" + xml_tracker;
            pos.child_tag = logic.shape_tag;
            pos.trans.dz = iter->getZOffset() + shape.dz;
            p.push_back(pos);
        }
    }
    
    /**
     * This is one of the smaller analysis functions that provide the core functionality of this class. It does a number of things:
     * it creates a composite material information struct for each support volume, and it adds the remaining information about
     * that volume to the collections of hierarchy, shape, position and topology information. All of these future CMSSW XML
     * blocks are given unique names based on the properties of the support structures they came from.
     * @param is A reference to the collection of inactive surfaces that the supports are a part of; used as input
     * @param c A reference to the collection of composite material information; used for output
     * @param l A reference to the collection of volume hierarchy information; used for output
     * @param s A reference to the collection of shape parameters; used for output
     * @param p A reference to the collection of volume positionings; used for output
     * @param t A reference to the collection of topology information; used for output
     */
    void Extractor::analyseSupports(InactiveSurfaces& is, std::vector<Composite>& c, std::vector<LogicalInfo>& l,
            std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt) {
        std::string nspace;
        if (wt) nspace = xml_newfileident;
        else nspace = xml_fileident;
        // container inits
        ShapeInfo shape;
        LogicalInfo logic;
        PosInfo pos;
        shape.type = tb;
        shape.dx = 0.0;
        shape.dy = 0.0;
        shape.dyy = 0.0;
        pos.copy = 1;
        pos.trans.dx = 0.0;
        pos.trans.dy = 0.0;
        // b_sup, e_sup, o_sup, t_sup, u_sup: one composite per category
        // l, s and p: one entry per support part
        std::set<MaterialProperties::Category> found;
        std::set<MaterialProperties::Category>::iterator fres;
        std::vector<InactiveElement>::iterator iter, guard;
        std::vector<InactiveElement>& sp = is.getSupports();
        guard = sp.end();
        // support volume loop
        for (iter = sp.begin(); iter != guard; iter++) {
            std::ostringstream matname, shapename;
            matname << xml_base_lazycomp << iter->getCategory();
            shapename << xml_base_lazy << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset()));
            fres = found.find(iter->getCategory());
            if (fres == found.end()) {
                c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
                found.insert(iter->getCategory());
            }
            shape.name_tag = shapename.str();
            shape.dz = iter->getZLength() / 2.0;
            shape.rmin = iter->getInnerRadius();
            shape.rmax = shape.rmin + iter->getRWidth();
            s.push_back(shape);
            logic.name_tag = shapename.str();
            logic.shape_tag = nspace + ":" + shapename.str();
            logic.material_tag = nspace + ":" + matname.str();
            l.push_back(logic);
            pos.parent_tag = nspace + ":" + xml_tracker;
            pos.child_tag = logic.shape_tag;
            if ((iter->getCategory() == MaterialProperties::o_sup) ||
                    (iter->getCategory() == MaterialProperties::t_sup)) pos.trans.dz = 0.0;
            else pos.trans.dz = iter->getZOffset() + shape.dz;
            p.push_back(pos);
        }
    }
    
    //private
    /**
     * Bundle the information for a composite material from a list of components into an instance of a <i>Composite</i> struct.
     * The function includes an option to skip the sensor silicon in the list of elementary materials, omitting it completely when
     * assembling the information for the <i>Composite</i> struct. This is useful for active volumes, where the sensor silicon
     * is assigned to a separate volume during the tranlation to CMSSW XML. On the other hand, an inactive volume will want
     * to include all of its material components - regardless of what they're labelled.
     * @param name The name of the new composite material
     * @param density The overall density of the new composite material
     * @param mp A reference to the material properties object that provides the list of elementary materials
     * @param nosensors A flag to include or exclude sensor silicon from the material list; true if it is excluded, false otherwise
     * @return A new instance of a <i>Composite</i> struct that bundles the material information for further processing
     */
    Composite Extractor::createComposite(std::string name, double density, MaterialProperties& mp, bool nosensors) {
        Composite comp;
        comp.name = name;
        comp.density = density;
        comp.method = wt;
        double m = 0.0;
        for (unsigned int i = 0; i < mp.localMassCount(); i++) {
            if (!nosensors || (mp.getLocalTag(i).compare(xml_sensor_silicon) != 0)) {
                std::pair<std::string, double> p;
                p.first = mp.getLocalTag(i);
                p.second = mp.getLocalMass(i);
                comp.elements.push_back(p);
                m = m + mp.getLocalMass(i);
            }
        }
        for (unsigned int i = 0; i < mp.exitingMassCount(); i++) {
            if (!nosensors || (mp.getExitingTag(i).compare(xml_sensor_silicon) != 0)) {
                std::pair<std::string, double> p;
                p.first = mp.getExitingTag(i);
                p.second = mp.getExitingMass(i);
                bool found = false;
                std::vector<std::pair<std::string, double> >::iterator iter, guard = comp.elements.end();
                for (iter = comp.elements.begin(); iter != guard; iter++) {
                    if (iter->first == p.first) {
                        found = true;
                        break;
                    }
                }
                if (found) iter->second = iter->second + p.second;
                else comp.elements.push_back(p);
                m = m + mp.getExitingMass(i);
            }
        }
        for (unsigned int i = 0; i < comp.elements.size(); i++)
            comp.elements.at(i).second = comp.elements.at(i).second / m;
        return comp;
    }
    
    /**
     * Find the partner module of a given one in a layer, i.e. a module that is on the same rod but on the opposite side of z=0.
     * @param i An iterator pointing to the start of the search range
     * @param g An iterator pointing to one past the end of the search range
     * @param ponrod The position along the rod of the original module
     * @param find_first A flag indicating whether to stop the search at the first module with the desired position, regardless of which side of z=0 it is on; default is false
     * @return An iterator pointing to the partner module, or to one past the end of the range if no partner is found
     */
    std::vector<ModuleCap>::iterator Extractor::findPartnerModule(std::vector<ModuleCap>::iterator i,
            std::vector<ModuleCap>::iterator g, int ponrod, bool find_first) {
        std::vector<ModuleCap>::iterator res = i;
        if (i != g) {
            bool plus = false;
            if (!find_first) plus = i->getModule().getMeanPoint().Z() > 0;
            while (res != g) {
                if (res->getModule().getRing() == ponrod) {
                    if (find_first) break;
                    else {
                        if((plus && (res->getModule().getMeanPoint().Z() < 0))
                                || (!plus && (res->getModule().getMeanPoint().Z() > 0))) break;
                    }
                }
                res++;
            }
        }
        return res;
    }
    
    /**
     * Find the total width in r of the volume enclosing a layer.
     * @param start An iterator pointing to the start of a range of modules in a vector
     * @param stop An iterator pointing to one past the end of a range of modules in a vector
     * @param middle The midpoint of the layer radius range
     * @return The thickness of the layer, or zero if the layer is degenerate
     */
    double Extractor::findDeltaR(std::vector<Module*>::iterator start,
            std::vector<Module*>::iterator stop, double middle) {
        std::vector<Module*>::iterator iter, mod1, mod2;
        double dr = 0.0;
        iter = start;
        mod1 = stop;
        mod2 = stop;
        for (iter = start; iter != stop; iter++) {
            if ((*iter)->getMeanPoint().Rho() > middle) {
                mod1 = iter;
                break;
            }
        }
        for (iter = mod1; iter != stop; iter++) {
            if ((*iter)->getMeanPoint().Rho() > middle) {
                if ((*iter)->getMeanPoint().Rho() < (*mod1)->getMeanPoint().Rho()) {
                    mod2 = iter;
                    break;
                }
                else if (!((*iter)->getMeanPoint().Rho() == (*mod1)->getMeanPoint().Rho())) {
                    mod2 = mod1;
                    mod1 = iter;
                    break;
                }
            }
        }
        dr = (*mod1)->getMinRho() - (*mod2)->getMinRho() + (*mod1)->getModuleThickness();
        return dr;
    }
    
    /**
     * Find half the width in z of the volume enclosing an endcap ring.
     * @param start An iterator pointing to the start of a range of modules in a vector
     * @param stop An iterator pointing to one past the end of a range of modules in a vector
     * @param middle The midpoint in z of the disc that the ring in question belongs to
     * @return Half the thickness of the ring volume, or zero if the module collection is degenerate
     */
    double Extractor::findDeltaZ(std::vector<Module*>::iterator start,
            std::vector<Module*>::iterator stop, double middle) {
        std::vector<Module*>::iterator iter, mod1, mod2;
        double dz = 0.0;
        iter = start;
        mod1 = stop;
        mod2 = stop;
        for (iter = start; iter != stop; iter++) {
            if ((*iter)->getMinZ() > middle) {
                mod1 = iter;
                break;
            }
        }
        for (iter = mod1; iter != stop; iter++) {
            if ((*iter)->getMinZ() > middle) {
                if ((*iter)->getMinZ() < (*mod1)->getMinZ()) {
                    mod2 = iter;
                    break;
                }
                else if (!((*iter)->getMinZ() == (*mod1)->getMinZ())) {
                    mod2 = mod1;
                    mod1 = iter;
                    break;
                }
            }
        }
        dz = (*mod1)->getMaxZ() - (*mod2)->getMinZ();
        return dz;
    }
    
    /**
     * Given the name of a <i>SpecPar</i> block, find the corresponding index in the vector of <i>SpecParInfo</i> structs.
     * @param specs The list of <i>SpecParInfo</i> structs collecting topological information for the tracker
     * @param label The name of the requested <i>SpecPar</i> block
     * @return The index of the requested struct, or <i>-1</i> if the name was not found in the collection
     */
    int Extractor::findSpecParIndex(std::vector<SpecParInfo>& specs, std::string label) {
        int idx = 0, size = (int)(specs.size());
        while (idx < size) {
            if (specs.at(idx).name.compare(label) == 0) return idx;
            idx++;
        }
        return -1;
    }
    
    /**
     * Calculate the thickness of the sensor material in a module from the amount of sensor silicon <i>SenSi</i>
     * and the dimensions of the module.
     * @param mc A reference to the module materials
     * @param mt A reference to the global list of material properties
     * @return The resulting sensor thickness
     */
    double Extractor::calculateSensorThickness(ModuleCap& mc, MaterialTable& mt) {
        double t = 0.0;
        double m = 0.0, d = 0.0;
        for (unsigned int i = 0; i < mc.localMassCount(); i++) {
            if (mc.getLocalTag(i).compare(xml_sensor_silicon) == 0) m = m + mc.getLocalMass(i);
        }
        for (unsigned int i = 0; i < mc.exitingMassCount(); i++) {
            if (mc.getExitingTag(i).compare(xml_sensor_silicon) == 0) m = m + mc.getExitingMass(i);
        }
        try { d = mt.getMaterial(xml_sensor_silicon).density; }
        catch (std::exception& e) { return 0.0; }
        t = 1000 * m / (d * mc.getSurface());
        return t;
    }
    
    /**
     * Pre-format a named string type parameter as a CMSSW XML string.
     * @param name The name of the string parameter
     * @param value The value of the string parameter, i.e. the string itself
     * @return The formatted XML tag
     */
    std::string Extractor::stringParam(std::string name, std::string value) {
        std::string res;
        res = xml_algorithm_string + name + xml_algorithm_value + value + xml_general_endline;
        return res;
    }
    
    /**
     * Pre-format a named numeric parameter as a CMSSW XML string.
     * @param name The name of the parameter
     * @param value The value of the parameter in string form
     * @return The formatted numeric parameter tag
     */
    std::string Extractor::numericParam(std::string name, std::string value) {
        std::string res;
        res = xml_algorithm_numeric + name + xml_algorithm_value + value + xml_general_endline;
        return res;
    }
    
    /**
     * Pre-format a 3D vector parameter as a CMSSW XML string.
     * @param x The coordinate along x
     * @param y The coordinate along y
     * @param z The coordinate along z
     * @return The formatted vector parameter tag
     */
    std::string Extractor::vectorParam(double x, double y, double z) {
        std::ostringstream res;
        res << xml_algorithm_vector_open << x << "," << y << "," << z << xml_algorithm_vector_close;
        return res.str();
    }
    
    /**
     * Calculate the composite density of the material mix in a module. A boolean flag controls whether the sensor silicon
     * <i>SenSi</i> should be added to or excluded from the overall mix.
     * @param mc A reference to the module materials
     * @param nosensors Then sensor silicon flag; true if that material should be excluded, false otherwise
     * @return The calculated overall density in <i>g/cm3</i>
     */
    double Extractor::compositeDensity(ModuleCap& mc, bool nosensors) {
        double d = mc.getSurface() * mc.getModule().getModuleThickness();
        if (nosensors) {
            double m = 0.0;
            for (unsigned int i = 0; i < mc.localMassCount(); i++) {
                if (mc.getLocalTag(i).compare(xml_sensor_silicon) != 0) m = m + mc.getLocalMass(i);
            }
            for (unsigned int i = 0; i < mc.exitingMassCount(); i++) {
                if (mc.getExitingTag(i).compare(xml_sensor_silicon) != 0) m = m + mc.getExitingMass(i);
            }
            d = 1000 * m / d;
        }
        else d = 1000 * mc.getTotalMass() / d;
        return d;
    }
    
    /**
     * Compute the overall density of the materials in an inactive element.
     * @param ie A reference to the inactive surface
     * @return The calculated material density in <i>g/cm3</i>
     */
    double Extractor::compositeDensity(InactiveElement& ie) {
        double d = ie.getRWidth() + ie.getInnerRadius();
        d = d * d - ie.getInnerRadius() * ie.getInnerRadius();
        d = 1000 * ie.getTotalMass() / (PI * ie.getZLength() * d);
        return d;
    }
    
    /**
     * Calculate the radial distance of an outer rod surface from the outer limit of its layer.
     * @param r The outer radius of the layer
     * @param w Half the width of the rod
     * @return The distance of the rod surface from the layer arc
     */
    double Extractor::fromRim(double r, double w) {
        double s = asin(w / r);
        s = 1 - cos(s);
        s = s * r;
        return s;
    }
    
    /**
     * Calculate the atomic number of an elementary material from its radiation length and atomic weight.
     * @param x0 The radiation length
     * @param A The atomic weight
     * @return The atomic number
     */
    int Extractor::Z(double x0, double A) {
      // TODO: understand this: why do we need to
      // get an integer value as output?
      double d = 4 - 4 * (1.0 - 181.0 * A / x0);
      if (d > 0) return int(floor((sqrt(d) - 2.0) / 2.0 + 0.5));
      else return -1;
    }
}

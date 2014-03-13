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

    // These are ALL vectors!
    e.clear(); // Element
    c.clear(); // Composite
    l.clear(); // LogicalInfo
    s.clear(); // ShapeInfo
    p.clear(); // PosInfo
    a.clear(); // AlgoInfo
    r.clear(); // Rotation
    t.clear(); // SpecParInfo
    ri.clear(); // RILengthInfo

    // Initialization

    // Default positioning
    PosInfo pos;        
    pos.copy = 1;
    pos.trans.dx = 0.0;
    pos.trans.dy = 0.0;
    pos.trans.dz = 0.0;
    pos.rotref = "";

    // Initialise rotation list with Harry's tilt mod
    Rotation rot;
    rot.name = xml_barrel_tilt;
    rot.thetax = 90.0;
    rot.phix = 270.0;
    rot.thetay = 180.0;
    rot.phiy = 0.0;
    rot.thetaz = 90.0;
    rot.phiz = 0.0;
    r.push_back(rot);

    // EndcapRot is not needed any more
    /*
       rot.name = xml_endcap_rot;
       rot.thetax = 90.0;
       rot.phix = 90.0;
       rot.thetay = 90.0;
       rot.phiy = 180.0;
       rot.thetaz = 0.0;
       rot.phiz = 0.0;
       r.push_back(rot);
       */

    // These seem not to be needed here
    //LogicalInfo logic;
    //AlgoInfo alg;
    //SpecParInfo spec;

    // Define top-level barrel and endcap volume containers (polycone)
    // This just fill the polycone profiles of the two volumes
    ShapeInfo shape;
    if (!wt) {
      shape.type = pc;
      shape.name_tag = xml_tob;

      // Barrel
      analyseBarrelContainer(tr, shape.rzup, shape.rzdown);
      s.push_back(shape);
      std::cout << "Barrel container done." << std::endl;

      // Endcap
      analyseEndcapContainer(tr, shape.rzup, shape.rzdown);
      if (!(shape.rzup.empty() || shape.rzdown.empty())) {
        shape.name_tag = xml_tid;
        s.push_back(shape);
      }
      std::cout << "Endcap container done." << std::endl;
    }

    // Translate entries in mt to elementary materials
    analyseElements(mt, e);
    std::cout << "Elementary materials done." << std::endl;
    // Analyse barrel
    analyseLayers(mt, bc, tr, c, l, s, p, a, r, t, ri, wt);
    std::cout << "Barrel layers done." << std::endl;
    // Analyse endcaps
    analyseDiscs(mt, ec, tr, c, l, s, p, a, r, t, ri, wt);
    std::cout << "Endcap discs done." << std::endl;
    // Barrel services
    analyseBarrelServices(is, c, l, s, p, t);
    std::cout << "Barrel services done." << std::endl;
    // Endcap services
    analyseEndcapServices(is, c, l, s, p, t);
    std::cout << "Endcap services done." << std::endl;
    // Supports
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
    up.clear();
    down.clear();

    LayerAggregator lagg;
    t.accept(lagg);
    auto bl = lagg.getBarrelLayers();
    int n_of_layers = bl->size();
    for (int layer = 0; layer < n_of_layers; layer++) {
      is_short = (bl->at(layer)->minZ() > 0) || (bl->at(layer)->maxZ() < 0);
      if (is_short) {
        //short layer on z- side
        if (bl->at(layer)->maxZ() < 0) {
          //indices 0, 1
          if ((layer == 0) || ((layer == 1) && (previous_short))) {
            rz.first = bl->at(layer)->minR();
            rz.second = bl->at(layer)->minZ();
            up.push_back(rz);
          }
          else {
            //new barrel reached
            if (bl->at(layer)->minZ() != zmin) {
              //new layer sticks out compared to old layer
              if (bl->at(layer)->minZ() > zmin) rz.first = bl->at(layer)->minR();
              //old layer sticks out compared to new layer
              else rz.first = rmax;
              rz.second = zmin;
              up.push_back(rz);
              rz.second = bl->at(layer)->minZ();
              up.push_back(rz);
            }
          }
          //indices size - 2, size - 1
          if ((layer == n_of_layers - 1) || ((layer == n_of_layers - 2) && (previous_short))) {
            rz.first = bl->at(layer)->maxR();
            rz.second = bl->at(layer)->minZ();
            up.push_back(rz);
          }
        }
        //short layer on z+ side
        else {
          //indices 0, 1
          if ((layer == 0) || ((layer == 1) && (previous_short))) {
            rz.first = bl->at(layer)->minR();
            rz.second = bl->at(layer)->maxZ();
            down.push_back(rz);
          }
          else {
            //new barrel reached
            if (bl->at(layer)->maxZ() != zmax) {
              //new layer sticks out compared to old layer
              if (bl->at(layer)->maxZ() > zmax) rz.first = bl->at(layer)->minR();
              //old layer sticks out compared to new layer
              else rz.first = rmax;
              rz.second = zmax;
              down.push_back(rz);
              rz.second = bl->at(layer)->maxZ();
              down.push_back(rz);
            }
          }
          //indices size - 2, size - 1
          if ((layer == n_of_layers - 1) || ((layer == n_of_layers - 2) && (previous_short))) {
            rz.first = bl->at(layer)->maxR();
            rz.second = bl->at(layer)->minZ();
            down.push_back(rz);
          }
        }
      }
      //regular layer across z=0
      else {
        //index 0
        if (layer == 0) {
          rz.first = bl->at(layer)->minR();
          rz.second = bl->at(layer)->minZ();
          up.push_back(rz);
          rz.second = bl->at(layer)->maxZ();
          down.push_back(rz);
        }
        else {
          //new barrel reached
          if (bl->at(layer)->maxZ() != zmax) {
            //new layer sticks out compared to old layer
            if (bl->at(layer)->maxZ() > zmax) rz.first = bl->at(layer)->minR();
            //old layer sticks out compared to new layer
            else rz.first = rmax;
            rz.second = zmin;
            up.push_back(rz);
            rz.second = zmax;
            down.push_back(rz);
            rz.second = bl->at(layer)->minZ();
            up.push_back(rz);
            rz.second = bl->at(layer)->maxZ();
            down.push_back(rz);
          }
        }
        //index size - 1
        if (layer == n_of_layers - 1) {
          rz.first = bl->at(layer)->maxR();
          rz.second = bl->at(layer)->minZ();
          up.push_back(rz);
          rz.second = bl->at(layer)->maxZ();
          down.push_back(rz);
        }
      }
      rmax = bl->at(layer)->maxR();
      if (bl->at(layer)->minZ() < 0) zmin = bl->at(layer)->minZ();
      if (bl->at(layer)->maxZ() > 0) zmax = bl->at(layer)->maxZ();
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
    LayerAggregator lagg; t.accept(lagg);
    auto el = lagg.getEndcapLayers();
    first = 0;
    last = el->size();
    while (first < last) {
      if (el->at(first)->maxZ() > 0) break;
      first++;
    }
    up.clear();
    down.clear();
    for (int i = first; i < last; i++) {
      // special treatment for first disc
      if (i == first) {
        rmin = el->at(i)->minR();
        rmax = el->at(i)->maxR();
        rz.first = rmax;
        rz.second = el->at(i)->minZ() - xml_z_pixfwd;
        up.push_back(rz);
        rz.first = rmin;
        down.push_back(rz);
      }
      // disc beyond the first
      else {
        // endcap change larger->smaller
        if (rmax > el->at(i)->maxR()) {
          rz.second = zmax - xml_z_pixfwd;
          rz.first = rmax;
          up.push_back(rz);
          rz.first = rmin;
          down.push_back(rz);
          rmax = el->at(i)->maxR();
          rmin = el->at(i)->minR();
          rz.first = rmax;
          up.push_back(rz);
          rz.first = rmin;
          down.push_back(rz);
        }
        // endcap change smaller->larger
        if (rmax < el->at(i)->maxR()) {
          rz.second = el->at(i)->minZ() - xml_z_pixfwd;
          rz.first = rmax;
          up.push_back(rz);
          rz.first = rmin;
          down.push_back(rz);
          rmax = el->at(i)->maxR();
          rmin = el->at(i)->minR();
          rz.first = rmax;
          up.push_back(rz);
          rz.first = rmin;
          down.push_back(rz);
        }
      }
      zmax = el->at(i)->maxZ();
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

    // Container inits
    ShapeInfo shape;
    shape.dyy = 0.0;

    LogicalInfo logic;

    PosInfo pos;
    pos.copy = 1;
    pos.trans.dx = 0.0;
    pos.trans.dy = 0.0;
    pos.trans.dz = 0.0;

    AlgoInfo alg;

    Rotation rot;
    rot.phix = 0.0;
    rot.phiy = 0.0;
    rot.phiz = 0.0;
    rot.thetax = 0.0;
    rot.thetay = 0.0;
    rot.thetaz = 0.0;

    ModuleROCInfo minfo;
    ModuleROCInfo minfo_zero={}; 
    SpecParInfo rocdims, lspec, rspec, mspec;
    // Layer
    lspec.name = xml_subdet_layer + xml_par_tail;
    lspec.parameter.first = xml_tkddd_structure;
    lspec.parameter.second = xml_det_layer;
    // Rod
    rspec.name = xml_subdet_rod + xml_par_tail;
    rspec.parameter.first = xml_tkddd_structure;
    rspec.parameter.second = xml_det_rod;
    // Module
    mspec.name = xml_subdet_tobdet + xml_par_tail;
    mspec.parameter.first = xml_tkddd_structure;
    mspec.parameter.second = xml_det_tobdet;

    // name goes into trackerStructureTopology
    // parameter.first goes into trackerStructureTopology
    // parameter.second goes into trackerStructureTopology?
    //   >> well, as long as the values are TOBLayer, TOBRod, TOBDet
    //   >> they do not go into any file!

    RILengthInfo ril;
    ril.barrel = true;
    ril.index = 0;

    // b_mod: one composite for every module position on rod
    // s and l: one entry for every module position on rod (box), one for every layer (tube), rods TBD
    // p: one entry for every layer (two for short layers), two modules, one wafer and active for each ring on rod
    // a: rods within layer (twice in case of a short layer)
    layer = 1;
    alg.name = xml_tobalgo;
    oguard = bc.end();

    LayerAggregator lagg;
    tr.accept(lagg);

    // barrel caps layer loop
    for (oiter = bc.begin(); oiter != oguard; oiter++) {
      struct ThicknessVisitor : public ConstGeometryVisitor {
        double max = 0;
        void visit(const Module& m) { max = MAX(max, m.thickness()); }
      };
      ThicknessVisitor v;
      lagg.getBarrelLayers()->at(layer-1)->accept(v);

      double rmin = lagg.getBarrelLayers()->at(layer - 1)->minR();
      rmin = rmin - v.max / 2.0;
      double rmax = lagg.getBarrelLayers()->at(layer - 1)->maxR();
      rmax = rmax + v.max / 2.0;
      double zmin = lagg.getBarrelLayers()->at(layer - 1)->minZ();
      double zmax = lagg.getBarrelLayers()->at(layer - 1)->maxZ();
      // CUIDADO maxR and minR should take into account dsDistance and sensorThickness (check if it's been already fixed in DetectorModule!!)
      double deltar = rmax - rmin; //findDeltaR(lagg.getBarrelLayers()->at(layer - 1)->getModuleVector()->begin(),
      //           lagg.getBarrelLayers()->at(layer - 1)->getModuleVector()->end(), (rmin + rmax) / 2.0);

      double rodThickness = lagg.getBarrelLayers()->at(layer - 1)->rods().at(0).thickness();
      double ds, dt = 0.0;
      double rtotal = 0.0, itotal = 0.0;

      int count = 0;

      if (deltar == 0.0) continue;

      bool is_short = (zmax < 0.0) || (zmin > 0.0);
      bool is_relevant = !is_short || (zmin > 0.0);

      if (is_relevant) {

        shape.type = bx; // box
        shape.rmin = 0.0;
        shape.rmax = 0.0;

        ril.index = layer;

        std::set<int> rings;
        std::ostringstream lname, rname, pconverter;

        lname << xml_layer << layer;
        rname << xml_rod << layer;

        // lname and rname are Layer1 and Rod1 and go into all files

        iguard = oiter->end();

        // module caps loop
        for (iiter = oiter->begin(); iiter != iguard; iiter++) {
          int modRing = iiter->getModule().uniRef().ring;
          if (rings.find(modRing) == rings.end()) {

            // This is the Barrel Case
            std::vector<ModuleCap>::iterator partner;
            std::ostringstream matname, shapename, specname;

            // module composite material
            matname << xml_base_actcomp << "L" << layer << "P" << modRing; 
            c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));

            // module box
            shapename << modRing << lname.str();
            shape.name_tag = xml_barrel_module + shapename.str();

            shape.dx = iiter->getModule().area() / iiter->getModule().length() / 2.0;
            shape.dy = iiter->getModule().length() / 2.0;
            shape.dz = iiter->getModule().thickness() / 2.0;
            s.push_back(shape);

            logic.name_tag = shape.name_tag;
            logic.shape_tag = nspace + ":" + logic.name_tag;
            logic.material_tag = nspace + ":" + matname.str();
            l.push_back(logic);

            // name_tag is BModule1Layer1 and it goes into all files

            pos.child_tag = logic.shape_tag;
            pos.rotref = nspace + ":" + xml_barrel_tilt;

            if ((iiter->getModule().center().Rho() > (rmax - deltar / 2.0))
                || ((iiter->getModule().center().Rho() < ((rmin + rmax) / 2.0))
                    && (iiter->getModule().center().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dz;
            else pos.trans.dx = shape.dz - deltar / 2.0;

            if (is_short) {
              pos.parent_tag = nspace + ":" + rname.str() + xml_plus;
              pos.trans.dz = iiter->getModule().minZ() - ((zmax + zmin) / 2.0) + shape.dy;
              p.push_back(pos);
              pos.parent_tag = nspace + ":" + rname.str() + xml_minus;
              pos.trans.dz = -pos.trans.dz;
              pos.copy = 2; // This is a copy of the BModule (FW/BW barrel half)
              p.push_back(pos);
              pos.copy = 1;
            } else {
              pos.parent_tag = nspace + ":" + rname.str();
              partner = findPartnerModule(iiter, iguard, modRing);
              if (iiter->getModule().uniRef().side > 0) { 
                pos.trans.dz = iiter->getModule().maxZ() - shape.dy;
                p.push_back(pos);
                if (partner != iguard) {
                  if ((partner->getModule().center().Rho() > (rmax - deltar / 2.0))
                      || ((partner->getModule().center().Rho() < ((rmin + rmax) / 2.0))
                          && (partner->getModule().center().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dz;
                  else pos.trans.dx = shape.dz - deltar / 2.0;
                  pos.trans.dz = partner->getModule().maxZ() - shape.dy;
                  pos.copy = 2; // This is a copy of the BModule (FW/BW barrel half)
                  p.push_back(pos);
                  pos.copy = 1;
                }
              } else {
                pos.trans.dz = iiter->getModule().maxZ() - shape.dy;
                pos.copy = 2; // This is a copy of the BModule (FW/BW barrel half)
                p.push_back(pos);
                pos.copy = 1;
                if (partner != iguard) {
                  if ((partner->getModule().center().Rho() > (rmax - deltar / 2.0))
                      || ((partner->getModule().center().Rho() < ((rmin + rmax) / 2.0))
                          && (partner->getModule().center().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dz;
                  else pos.trans.dx = shape.dz - deltar / 2.0;
                  pos.trans.dz = partner->getModule().maxZ() - shape.dy;
                  p.push_back(pos);
                }
              }
            }
            pos.rotref = "";

            // wafer
            string xml_base_inout = "";
            if (iiter->getModule().numSensors() == 2) xml_base_inout = xml_base_inner;

            shape.name_tag = xml_barrel_module + shapename.str() + xml_base_inout + xml_base_waf;
            shape.dz = iiter->getModule().sensorThickness() / 2.0; //CUIDADO WAS calculateSensorThickness(*iiter, mt) / 2.0;
            //if (iiter->getModule().numSensors() == 2) shape.dz = shape.dz / 2.0; // CUIDADO calcSensThick returned 2x what getSensThick returns, it means that now one-sided sensors are half as thick if not compensated for in the config files
            s.push_back(shape);

            pos.parent_tag = logic.shape_tag;

            logic.name_tag = shape.name_tag;
            logic.shape_tag = nspace + ":" + logic.name_tag;
            logic.material_tag = xml_material_air;
            l.push_back(logic);

            pos.child_tag = logic.shape_tag;
            pos.trans.dx = 0.0;
            //pos.trans.dz = shape.dz - iiter->getModule().tmoduleThickness() / 2.0;
            pos.trans.dz = /*shape.dz*/ - iiter->getModule().dsDistance() / 2.0; 
            p.push_back(pos);

            if (iiter->getModule().numSensors() == 2) {

              xml_base_inout = xml_base_outer;

              shape.name_tag = xml_barrel_module + shapename.str() + xml_base_inout + xml_base_waf;
              s.push_back(shape);

              logic.name_tag = shape.name_tag;
              logic.shape_tag = nspace + ":" + logic.name_tag;
              l.push_back(logic);

              pos.child_tag = logic.shape_tag;
              pos.trans.dz = pos.trans.dz + /*2 * shape.dz +*/ iiter->getModule().dsDistance();  // CUIDADO: was with 2*shape.dz, but why???
              //pos.copy = 2;

              if (iiter->getModule().stereoRotation() != 0) {
                rot.name = type_stereo + xml_barrel_module + shapename.str();
                rot.thetax = 90.0;
                rot.phix = iiter->getModule().stereoRotation() / M_PI * 180;
                rot.thetay = 90.0;
                rot.phiy = 90.0 + iiter->getModule().stereoRotation() / M_PI * 180;
                r.push_back(rot);
                pos.rotref = nspace + ":" + rot.name;
              }
              p.push_back(pos);

              // Now reset
              pos.rotref.clear();
              rot.name.clear();
              rot.thetax = 0.0;
              rot.phix = 0.0;
              rot.thetay = 0.0;
              rot.phiy = 0.0;
              pos.copy = 1;
            }



            // active surface
            xml_base_inout = "";
            if (iiter->getModule().numSensors() == 2) xml_base_inout = xml_base_inner;

            shape.name_tag = xml_barrel_module + shapename.str() + xml_base_inout + xml_base_act;
            s.push_back(shape);

            //pos.parent_tag = logic.shape_tag;
            pos.parent_tag = nspace + ":" + xml_barrel_module + shapename.str() + xml_base_inout + xml_base_waf;

            logic.name_tag = shape.name_tag;
            logic.shape_tag = nspace + ":" + logic.name_tag;
            logic.material_tag = nspace + ":" + xml_sensor_silicon;
            l.push_back(logic);

            pos.child_tag = logic.shape_tag;
            pos.trans.dz = 0.0;
            p.push_back(pos);

            // topology
            mspec.partselectors.push_back(logic.name_tag);

            minfo.name		= iiter->getModule().moduleType();
            minfo.rocrows	= any2str<int>(iiter->getModule().innerSensor().numROCRows());  // in case of single sensor module innerSensor() and outerSensor() point to the same sensor
            minfo.roccols	= any2str<int>(iiter->getModule().innerSensor().numROCCols());
            minfo.rocx		= any2str<int>(iiter->getModule().innerSensor().numROCX());
            minfo.rocy		= any2str<int>(iiter->getModule().innerSensor().numROCY());

            mspec.moduletypes.push_back(minfo);

            if (iiter->getModule().numSensors() == 2) { 

              // active surface
              xml_base_inout = xml_base_outer;

              shape.name_tag = xml_barrel_module + shapename.str() + xml_base_inout + xml_base_act;
              s.push_back(shape);

              pos.parent_tag = nspace + ":" + xml_barrel_module + shapename.str() + xml_base_inout + xml_base_waf;

              logic.name_tag = shape.name_tag;
              logic.shape_tag = nspace + ":" + logic.name_tag;
              l.push_back(logic);

              pos.child_tag = logic.shape_tag;
              p.push_back(pos);

              // topology
              mspec.partselectors.push_back(logic.name_tag);

              minfo.rocrows	= any2str<int>(iiter->getModule().outerSensor().numROCRows());
              minfo.roccols	= any2str<int>(iiter->getModule().outerSensor().numROCCols());
              minfo.rocx		= any2str<int>(iiter->getModule().outerSensor().numROCX());
              minfo.rocy		= any2str<int>(iiter->getModule().outerSensor().numROCY());

              mspec.moduletypes.push_back(minfo);
            } // End of replica for Pt-modules

            // material properties
            rtotal = rtotal + iiter->getRadiationLength();
            itotal = itotal + iiter->getInteractionLength();
            count++;
            rings.insert(modRing);
            dt = iiter->getModule().thickness();
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
        shape.dy = shape.dx; // CUIDADO EXPERIMENTAL - Rods must have the dx of modules as dy because modules are afterwards rotated with the Harry's tilt mod 
        shape.dx = rodThickness / 2.0;
        if (is_short) shape.dz = (zmax - zmin) / 2.0;
        else shape.dz = zmax;
        s.push_back(shape);
        logic.name_tag = shape.name_tag;
        logic.shape_tag = nspace + ":" + logic.name_tag;
        logic.material_tag = xml_material_air;
        l.push_back(logic);
        rspec.partselectors.push_back(logic.name_tag);
        rspec.moduletypes.push_back(minfo_zero);
        //rspec.moduletypes.push_back("");
        pconverter << logic.shape_tag;
        if (is_short) {
          shape.name_tag = rname.str() + xml_minus;
          s.push_back(shape);
          logic.name_tag = shape.name_tag;
          logic.shape_tag = nspace + ":" + logic.name_tag;
          l.push_back(logic);
          rspec.partselectors.push_back(logic.name_tag);
          //rspec.moduletypes.push_back("");
          rspec.moduletypes.push_back(minfo_zero);
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
        //lspec.moduletypes.push_back("");
        lspec.moduletypes.push_back(minfo_zero);
        // rods in layer algorithm(s)
        alg.parent = logic.shape_tag;
        if (wt && is_short) alg.parent = alg.parent + xml_plus;
        alg.parameters.push_back(stringParam(xml_childparam, pconverter.str()));
        pconverter.str("");
        pconverter << (lagg.getBarrelLayers()->at(layer - 1)->tilt() + 90) << "*deg";
        alg.parameters.push_back(numericParam(xml_tilt, pconverter.str()));
        pconverter.str("");
        pconverter << lagg.getBarrelLayers()->at(layer - 1)->startAngle();
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
        pconverter << lagg.getBarrelLayers()->at(layer - 1)->numRods();
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

    // Container inits
    ShapeInfo shape;
    shape.dyy = 0.0;

    LogicalInfo logic;

    PosInfo pos;
    pos.copy = 1;
    pos.trans.dx = 0.0;
    pos.trans.dy = 0.0;
    pos.trans.dz = 0.0;

    AlgoInfo alg;

    Rotation rot;
    rot.phix = 0.0;
    rot.phiy = 0.0;
    rot.phiz = 0.0;
    rot.thetax = 0.0;
    rot.thetay = 0.0;
    rot.thetaz = 0.0;

    ModuleROCInfo minfo;
    ModuleROCInfo minfo_zero={}; 
    SpecParInfo rocdims, dspec, rspec, mspec;
    // Disk
    dspec.name = xml_subdet_wheel + xml_par_tail;
    dspec.parameter.first = xml_tkddd_structure;
    dspec.parameter.second = xml_det_wheel;
    // Ring
    rspec.name = xml_subdet_ring + xml_par_tail;
    rspec.parameter.first = xml_tkddd_structure;
    rspec.parameter.second = xml_det_ring;
    // Module
    mspec.name = xml_subdet_tiddet + xml_par_tail;
    mspec.parameter.first = xml_tkddd_structure;
    mspec.parameter.second = xml_det_tiddet;

    RILengthInfo ril;
    ril.barrel = false;
    ril.index = 0;

    // e_mod: one composite for every ring
    // s and l: one entry for every ring module, one for every ring, one for every disc
    // p: one entry for every disc, one for every ring, one module, wafer and active per ring
    // a: two per ring with modules inside ring
    layer = 1;
    alg.name = xml_ecalgo;
    oguard = ec.end();

    LayerAggregator lagg;
    tr.accept(lagg);

    // endcap caps layer loop
    for (oiter = ec.begin(); oiter != oguard; oiter++) {
      if (lagg.getEndcapLayers()->at(layer - 1)->minZ() > 0) {

        ril.index = layer;
        std::set<int> ridx;
        std::map<int, RingInfo> rinfo;

        struct ThicknessVisitor : public ConstGeometryVisitor {
          double max = 0;
          void visit(const Module& m) { max = MAX(max, m.thickness()); }
        };
        ThicknessVisitor v;
        lagg.getEndcapLayers()->at(layer-1)->accept(v);
		// CUIDADO refactor using Disk methods for min/max Z,R
        double rmin = lagg.getEndcapLayers()->at(layer - 1)->minR();
        double rmax = lagg.getEndcapLayers()->at(layer - 1)->maxR();
        double zmax = lagg.getEndcapLayers()->at(layer - 1)->maxZ();
        zmax = zmax + v.max / 2.0; 
        double zmin = lagg.getEndcapLayers()->at(layer - 1)->minZ();
        zmin = zmin - v.max / 2.0;

        std::ostringstream dname, pconverter;

        double rtotal = 0.0, itotal = 0.0;
        int count = 0;
        dname << xml_disc << layer;

        //shape.type = tp;
        shape.rmin = 0.0;
        shape.rmax = 0.0;
        pos.trans.dz = 0.0;
        iguard = oiter->end();

        // endcap module caps loop
        for (iiter = oiter->begin(); iiter != iguard; iiter++) {
          int modRing = iiter->getModule().uniRef().ring;
          // new ring
          if (ridx.find(modRing) == ridx.end()) {

            // This is the Barrel Case
            ridx.insert(modRing);
            std::ostringstream matname, rname, mname, specname;

            // module composite material
            matname << xml_base_actcomp << "D" << layer << "R" << modRing;
            c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));

            rname << xml_ring << modRing << dname.str();
            mname << xml_endcap_module << modRing << dname.str();

            // collect ring info
            RingInfo rinf;
            rinf.name = rname.str();
            rinf.childname = mname.str();
            rinf.fw = (iiter->getModule().center().Z() < (zmin + zmax) / 2.0);
            rinf.modules = lagg.getEndcapLayers()->at(layer - 1)->rings().at(modRing-1).numModules();
            rinf.rin = iiter->getModule().minR();
            rinf.rout = iiter->getModule().maxR();
            rinf.rmid = iiter->getModule().center().Rho();
            rinf.mthk = iiter->getModule().thickness();
            rinf.phi = iiter->getModule().center().Phi();
            rinfo.insert(std::pair<int, RingInfo>(modRing, rinf));

            // module trapezoid

            shape.type = iiter->getModule().shape() == RECTANGULAR ? bx : tp;

            shape.name_tag = mname.str();
            shape.dx = iiter->getModule().minWidth() / 2.0;
            shape.dxx = iiter->getModule().maxWidth() / 2.0;
            shape.dy = iiter->getModule().length() / 2.0;
            shape.dyy = iiter->getModule().length() / 2.0;
            //shape.dx = iiter->getModule().length() / 2.0;
            //shape.dy = iiter->getModule().minWidth() / 2.0;
            //shape.dyy = iiter->getModule().maxWidth() / 2.0;
            shape.dz = iiter->getModule().thickness() / 2.0;
            s.push_back(shape);

            logic.name_tag = shape.name_tag;
            logic.shape_tag = nspace + ":" + logic.name_tag;
            logic.material_tag = nspace + ":" + matname.str();
            l.push_back(logic);



            // wafer -- same x and y size of parent shape, but different thickness
            string xml_base_inout = "";
            if (iiter->getModule().numSensors() == 2) xml_base_inout = xml_base_inner;


            pos.parent_tag = logic.shape_tag;

            shape.name_tag = mname.str() + xml_base_inout+ xml_base_waf;
            shape.dz = iiter->getModule().sensorThickness() / 2.0; // CUIDADO WAS calculateSensorThickness(*iiter, mt) / 2.0;
            //if (iiter->getModule().numSensors() == 2) shape.dz = shape.dz / 2.0; // CUIDADO calcSensThick returned 2x what getSensThick returns, it means that now one-sided sensors are half as thick if not compensated for in the config files
            s.push_back(shape);

            logic.name_tag = shape.name_tag;
            logic.shape_tag = nspace + ":" + logic.name_tag;
            logic.material_tag = xml_material_air;
            l.push_back(logic);

            pos.child_tag = logic.shape_tag;

            if (iiter->getModule().maxZ() > 0) pos.trans.dz = /*shape.dz*/ - iiter->getModule().dsDistance() / 2.0; // CUIDADO WAS getModule().moduleThickness()
            else pos.trans.dz = iiter->getModule().dsDistance() / 2.0 /*- shape.dz*/; // DITTO HERE
            p.push_back(pos);
            if (iiter->getModule().numSensors() == 2) {

              xml_base_inout = xml_base_outer;

              //pos.parent_tag = logic.shape_tag;

              shape.name_tag = mname.str() + xml_base_inout+ xml_base_waf;
              s.push_back(shape);

              logic.name_tag = shape.name_tag;
              logic.shape_tag = nspace + ":" + logic.name_tag;
              l.push_back(logic);

              pos.child_tag = logic.shape_tag;

              if (iiter->getModule().maxZ() > 0) pos.trans.dz = /*pos.trans.dz + 2 * shape.dz +*/  iiter->getModule().dsDistance() / 2.0; // CUIDADO removed pos.trans.dz + 2*shape.dz, added / 2.0
              else pos.trans.dz = /* pos.trans.dz - 2 * shape.dz -*/ - iiter->getModule().dsDistance() / 2.0;
              //pos.copy = 2;
              if (iiter->getModule().stereoRotation() != 0) {
                rot.name = type_stereo + xml_endcap_module + mname.str();
                rot.thetax = 90.0;
                rot.phix = iiter->getModule().stereoRotation() / M_PI * 180;
                rot.thetay = 90.0;
                rot.phiy = 90.0 + iiter->getModule().stereoRotation() / M_PI * 180;
                r.push_back(rot);
                pos.rotref = nspace + ":" + rot.name;
              }

              p.push_back(pos);

              // Now reset
              pos.rotref.clear();
              rot.name.clear();
              rot.thetax = 0.0;
              rot.phix = 0.0;
              rot.thetay = 0.0;
              rot.phiy = 0.0;
              pos.copy = 1;
            }



            // active surface
            xml_base_inout = "";
            if (iiter->getModule().numSensors() == 2) xml_base_inout = xml_base_inner;

            //pos.parent_tag = logic.shape_tag;
            pos.parent_tag = nspace + ":" + mname.str() + xml_base_inout + xml_base_waf;

            shape.name_tag = mname.str() + xml_base_inout + xml_base_act;
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

            minfo.name		= iiter->getModule().moduleType();
            minfo.rocrows	= any2str<int>(iiter->getModule().innerSensor().numROCRows());
            minfo.roccols	= any2str<int>(iiter->getModule().innerSensor().numROCCols());
            minfo.rocx		= any2str<int>(iiter->getModule().innerSensor().numROCX());
            minfo.rocy		= any2str<int>(iiter->getModule().innerSensor().numROCY());

            mspec.moduletypes.push_back(minfo);

            if (iiter->getModule().numSensors() == 2) {

              // active surface
              xml_base_inout = xml_base_outer;

              //pos.parent_tag = logic.shape_tag;
              pos.parent_tag = nspace + ":" + mname.str() + xml_base_inout + xml_base_waf;

              shape.name_tag = mname.str() + xml_base_inout + xml_base_act;
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

              minfo.rocrows	= any2str<int>(iiter->getModule().outerSensor().numROCRows());
              minfo.roccols	= any2str<int>(iiter->getModule().outerSensor().numROCCols());
              minfo.rocx		= any2str<int>(iiter->getModule().outerSensor().numROCX());
              minfo.rocy		= any2str<int>(iiter->getModule().outerSensor().numROCY());

              mspec.moduletypes.push_back(minfo);
              //mspec.moduletypes.push_back(iiter->getModule().getType());

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
        shape.dz = (zmax - zmin) / 2.0; //findDeltaZ(lagg.getEndcapLayers()->at(layer - 1)->getModuleVector()->begin(), // CUIDADO what the hell is this??
        //lagg.getEndcapLayers()->at(layer - 1)->getModuleVector()->end(), (zmin + zmax) / 2.0) / 2.0;

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

            pos.parent_tag = nspace + ":" + dname.str(); // CUIDADO ended with: + xml_plus;
            pos.child_tag = logic.shape_tag;

            if (rinfo[*siter].fw) pos.trans.dz = (zmin - zmax) / 2.0 + shape.dz;
            else pos.trans.dz = (zmax - zmin) / 2.0 - shape.dz;
            p.push_back(pos);
            //pos.parent_tag = nspace + ":" + dname.str(); // CUIDADO ended with: + xml_minus;
            //p.push_back(pos);

            rspec.partselectors.push_back(logic.name_tag);
            rspec.moduletypes.push_back(minfo_zero);

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

        logic.name_tag = shape.name_tag; // CUIDADO ended with + xml_plus;
        //logic.extra = xml_plus;
        logic.shape_tag = nspace + ":" + shape.name_tag;
        logic.material_tag = xml_material_air;
        l.push_back(logic);

        pos.parent_tag = xml_pixfwdident + ":" + xml_pixfwd;
        pos.child_tag = nspace + ":" + logic.name_tag;
        pos.trans.dz = (zmax + zmin) / 2.0 - xml_z_pixfwd;
        p.push_back(pos);

        dspec.partselectors.push_back(logic.name_tag);
        dspec.moduletypes.push_back(minfo_zero);
        dspec.partextras.push_back(logic.extra);
        //   logic.name_tag = shape.name_tag; // CUIDADO ended with + xml_minus;
        //   logic.extra = xml_minus;
        //   l.push_back(logic);
        //   pos.parent_tag = xml_pixfwdident + ":" + xml_pixfwd;
        //   pos.child_tag = nspace + ":" + logic.name_tag;
        //   p.push_back(pos);
        //dspec.partselectors.push_back(logic.name_tag); // CUIDADO dspec still needs to be duplicated for minus discs (I think)
        //dspec.partextras.push_back(logic.extra);
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
      if ((iter->getZOffset() + iter->getZLength()) > 0) { // This is necessary because of replication of Forward volumes!
        c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
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
    for (std::map<std::string, double>::const_iterator it = mp.getLocalMasses().begin(); it != mp.getLocalMasses().end(); ++it) {
      if (!nosensors || (it->first.compare(xml_sensor_silicon) != 0)) {
        //    std::pair<std::string, double> p;
        //    p.first = mp.getLocalTag(i);
        //    p.second = mp.getLocalMass(i);
        comp.elements.push_back(*it);
        //    m = m + mp.getLocalMass(i);
        m += it->second;
      }
    }
    for (std::map<std::string, double>::const_iterator it = mp.getExitingMasses().begin(); it != mp.getExitingMasses().end(); ++it) {
      if (!nosensors || (it->first.compare(xml_sensor_silicon) != 0)) {
        std::pair<std::string, double> p = *it;
        //    p.first = mp.getExitingTag(i);
        //    p.second = mp.getExitingMass(i);
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
        m += it->second; // mp.getExitingMass(i);
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
      if (!find_first) plus = i->getModule().uniRef().side > 0; //i->getModule().center().Z() > 0;
      while (res != g) {
        if (res->getModule().uniRef().ring == ponrod) {
          if (find_first) break;
          else {
            if((plus && res->getModule().uniRef().side < 0 /*(res->getModule().center().Z() < 0)*/)
               || (!plus && res->getModule().uniRef().side > 0 /*(res->getModule().center().Z() > 0)*/)) break;
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
      if ((*iter)->center().Rho() > middle) {
        mod1 = iter;
        break;
      }
    }
    for (iter = mod1; iter != stop; iter++) {
      if ((*iter)->center().Rho() > middle) {
        if ((*iter)->center().Rho() < (*mod1)->center().Rho()) {
          mod2 = iter;
          break;
        }
        else if (!((*iter)->center().Rho() == (*mod1)->center().Rho())) {
          mod2 = mod1;
          mod1 = iter;
          break;
        }
      }
    }
    dr = (*mod1)->minR() - (*mod2)->minR() + (*mod1)->thickness();
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
      if ((*iter)->minZ() > middle) {
        mod1 = iter;
        break;
      }
    }
    for (iter = mod1; iter != stop; iter++) {
      if ((*iter)->minZ() > middle) {
        if ((*iter)->minZ() < (*mod1)->minZ()) {
          mod2 = iter;
          break;
        }
        else if (!((*iter)->minZ() == (*mod1)->minZ())) {
          mod2 = mod1;
          mod1 = iter;
          break;
        }
      }
    }
    dz = (*mod1)->maxZ() - (*mod2)->minZ();
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
    for (std::map<std::string, double>::const_iterator it = mc.getLocalMasses().begin(); it != mc.getLocalMasses().end(); ++it) {
      if (it->first.compare(xml_sensor_silicon) == 0) m += it->second;
    }
    for (std::map<std::string, double>::const_iterator it = mc.getExitingMasses().begin(); it != mc.getExitingMasses().end(); ++it) {
      if (it->first.compare(xml_sensor_silicon) == 0) m += it->second;
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
    double d = mc.getSurface() * mc.getModule().thickness();
    if (nosensors) {
      double m = 0.0;
      for (std::map<std::string, double>::const_iterator it = mc.getLocalMasses().begin(); it != mc.getLocalMasses().end(); ++it) {
        if (it->first.compare(xml_sensor_silicon) != 0) m += it->second;
      }
      for (std::map<std::string, double>::const_iterator it = mc.getExitingMasses().begin(); it != mc.getExitingMasses().end(); ++it) {
        if (it->first.compare(xml_sensor_silicon) != 0) m += it->second;
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

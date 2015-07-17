/**
 * @file Extractor.cc
 * @brief This class extracts the tracker geometry and material properties from an existing setup and groups them in preparation for translation to CMSSW XML
 */


#define __USE_MODULECOMPLEX__
//#define __DEBUGPRINT__
//#define __FLIPSENSORS_OUT__
//#define __FLIPSENSORS_IN__

#include <Extractor.h>
#ifdef __USE_MODULECOMPLEX__
#include <cstdlib>
#endif
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
    std::map<std::string,Rotation>& r = d.rots;
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
    rot.name = xml_default_mod_rot;
    rot.thetax = 90.0;
    rot.phix = 90.0;
    rot.thetay = 0.0;
    rot.phiy = 0.0;
    rot.thetaz = 90.0;
    rot.phiz = 0.0;
    r.insert(std::pair<const std::string,Rotation>(rot.name,rot));

    // Flip module (Fix Y axis)
    rot.name = xml_flip_mod_rot;
    rot.thetax = 90.0;
    rot.phix = 180.0;
    rot.thetay = 90.0;
    rot.phiy = 90.0;
    rot.thetaz = 180.0;
    rot.phiz = 0.0;
    r.insert(std::pair<const std::string,Rotation>(rot.name,rot));


#if defined(__FLIPSENSORS_IN__) || defined(__FLIPSENSORS_OUT__)
    rot.name = rot_sensor_tag;
#if 1  // Fix Y axis case
    rot.thetax = 90.0;
    rot.phix   = 180.0;
    rot.thetay = 90.0;
    rot.phiy   = 90.0;
#else  // Fix X axis case
    rot.thetax = 90.0;
    rot.phix   = 0.0;
    rot.thetay = 90.0;
    rot.phiy   = 270.0;
#endif
    rot.thetaz = 180.0;
    rot.phiz   = 0.0;
    r.insert(std::pair<const std::string,Rotation>(rot.name,rot));
#endif

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
    analyseLayers(mt, tr, c, l, s, p, a, r, t, ri, wt);
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
  void Extractor::analyseLayers(MaterialTable& mt/*, std::vector<std::vector<ModuleCap> >& bc*/, Tracker& tr,
                                std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<PosInfo>& p,
                                std::vector<AlgoInfo>& a, std::map<std::string,Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, bool wt) {
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

    LayerAggregator lagg;
    tr.accept(lagg);
    lagg.postVisit();
    std::vector<std::vector<ModuleCap> >& bc = lagg.getBarrelCap();
    oguard = bc.end();

    // barrel caps layer loop
    for (oiter = bc.begin(); oiter != oguard; oiter++) {
      struct ThicknessVisitor : public ConstGeometryVisitor {
        double max = 0;
        void visit(const Module& m) { max = MAX(max, m.thickness()); }
      };
      ThicknessVisitor v;
      lagg.getBarrelLayers()->at(layer-1)->accept(v);
      
      double rmin = lagg.getBarrelLayers()->at(layer - 1)->minR();
      //rmin = rmin - v.max / 2.0; // no need to manually add/subtract the thickness, it is already taken into account by the new volumetric minR/maxR functions
      double rmax = lagg.getBarrelLayers()->at(layer - 1)->maxR();
      //rmax = rmax + v.max / 2.0;
      double zmin = lagg.getBarrelLayers()->at(layer - 1)->minZ();
      double zmax = lagg.getBarrelLayers()->at(layer - 1)->maxZ();
      double flatzmax;
      double rodThickness = lagg.getBarrelLayers()->at(layer - 1)->rodThickness();
      double deltar = rodThickness; //rmax - rmin; //findDeltaR(lagg.getBarrelLayers()->at(layer - 1)->getModuleVector()->begin(),
      // lagg.getBarrelLayers()->at(layer - 1)->getModuleVector()->end(), (rmin + rmax) / 2.0);

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

	std::set<int> ridx;
        std::map<int, BTiltedRingInfo> rinfoplus;
        std::map<int, BTiltedRingInfo> rinfominus;
        std::ostringstream lname, rodname, pconverter;

        lname << xml_layer << layer;
        rodname << xml_rod << layer;

        // lname and rodname are Layer1 and Rod1 and go into all files

        iguard = oiter->end();

        // module caps loop
        for (iiter = oiter->begin(); iiter != iguard; iiter++) {
          
	  if (iiter->getModule().uniRef().side > 0 && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)){

	    int modRing = iiter->getModule().uniRef().ring;
	    // Tilt angle of the module
	    double tiltAngle = 0;
	    if (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty()) {
	      tiltAngle = iiter->getModule().tiltAngle() * 180 / M_PI;
	    }
	    
	    //std::cout << "iiter->getModule().uniRef().phi = " << iiter->getModule().uniRef().phi << " iiter->getModule().center().Rho() = " << iiter->getModule().center().Rho() << " iiter->getModule().center().Z() = " << iiter->getModule().center().Z() << " iiter->getModule().flipped() = " << iiter->getModule().flipped() << std::endl;

            std::ostringstream mname;
            mname << xml_barrel_module << modRing << lname.str();
#ifdef __USE_MODULECOMPLEX__ 
            std::string parentName = mname.str(); 
            //if (iiter->getModule().flipped()) parentName = mname.str() + xml_flipped; 
            ModuleComplex modcomplex(mname.str(),parentName,*iiter);
            modcomplex.buildSubVolumes();
#endif
	    if (iiter->getModule().uniRef().phi == 1) {
	      
	      ridx.insert(modRing);
	      //if (ridx.find(modRing) == ridx.end()) {

            // This is the Barrel Case
            std::vector<ModuleCap>::iterator partner;
            std::ostringstream ringname, matname, specname;

#ifndef __USE_MODULECOMPLEX__ 
            // module composite material
            matname << xml_base_actcomp << "L" << layer << "P" << modRing;
            c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));
#endif

	    ringname << xml_ring << modRing << lname.str();

            // module box
            shape.name_tag = mname.str();
#ifndef __USE_MODULECOMPLEX__
            shape.dx = iiter->getModule().area() / iiter->getModule().length() / 2.0;
            shape.dy = iiter->getModule().length() / 2.0;
            shape.dz = iiter->getModule().thickness() / 2.0;
#else
            shape.dx = modcomplex.getExpandedModuleWidth()/2.0;
            shape.dy = modcomplex.getExpandedModuleLength()/2.0;
            shape.dz = modcomplex.getExpandedModuleThickness()/2.0;
#endif
	    

	    if (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0)) {
	      // collect ring info
	      BTiltedRingInfo rinf;
	      rinf.name = ringname.str() + xml_plus;
	      rinf.childname = mname.str() + xml_tilted + xml_plus;
	      rinf.flippedchildname = mname.str() + xml_tilted + xml_plus + xml_flipped;
	      rinf.inner_flipped = iiter->getModule().flipped();
	      rinf.r1 = iiter->getModule().center().Rho();
	      rinf.z1 = iiter->getModule().center().Z();
	      rinf.modules = lagg.getBarrelLayers()->at(layer - 1)->numRods();
#ifndef __USE_MODULECOMPLEX__
	      rinf.mdx = shape.dx;
	      rinf.mdy = shape.dy;
	      rinf.mdz = shape.dz;
#else
	      rinf.rmin = modcomplex.getRmin();
	      rinf.zmin = modcomplex.getZmin();
#endif
	      rinf.tiltAngle = tiltAngle;
	      rinf.phi = iiter->getModule().uniRef().phi;
	      rinfoplus.insert(std::pair<int, BTiltedRingInfo>(modRing, rinf));

	      rinf.name = ringname.str() + xml_minus;
	      rinf.childname = mname.str() + xml_tilted + xml_minus;
	      rinf.flippedchildname = mname.str() + xml_tilted + xml_minus + xml_flipped;
	      rinf.z1 = - iiter->getModule().center().Z();
	      rinfominus.insert(std::pair<int, BTiltedRingInfo>(modRing, rinf));

	      // module boxes
	      shape.name_tag = mname.str() + xml_tilted + xml_plus + xml_flipped;
	      s.push_back(shape);
	      shape.name_tag = mname.str() + xml_tilted + xml_plus;
	      s.push_back(shape);
	      shape.name_tag = mname.str() + xml_tilted + xml_minus + xml_flipped;
	      s.push_back(shape);
	      shape.name_tag = mname.str() + xml_tilted + xml_minus;
	      s.push_back(shape);
	    }
	    if (iiter->getModule().flipped() || (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0))) {
	      shape.name_tag = mname.str() + xml_flipped;
	      s.push_back(shape);
	    }
	    shape.name_tag = mname.str();
	    s.push_back(shape);

#ifdef __USE_MODULECOMPLEX__ 
            // Get it back for sensors
            shape.dx = iiter->getModule().area() / iiter->getModule().length() / 2.0;
            shape.dy = iiter->getModule().length() / 2.0;
            shape.dz = iiter->getModule().thickness() / 2.0;
#endif


#if 0
            logic.material_tag = nspace + ":" + matname.str();
#else
            logic.material_tag = xml_material_air;
#endif
	    if (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0)) {
	      logic.name_tag = mname.str() + xml_tilted + xml_plus + xml_flipped;
	      logic.shape_tag = nspace + ":" + logic.name_tag ;
	      l.push_back(logic);
	      logic.name_tag = mname.str() + xml_tilted + xml_plus;
	      logic.shape_tag = nspace + ":" + logic.name_tag;
	      l.push_back(logic);
	      logic.name_tag = mname.str() + xml_tilted + xml_minus + xml_flipped;
	      logic.shape_tag = nspace + ":" + logic.name_tag;
	      l.push_back(logic);
	      logic.name_tag = mname.str() + xml_tilted + xml_minus;
	      logic.shape_tag = nspace + ":" + logic.name_tag;
	      l.push_back(logic);
	    }
	    if (iiter->getModule().flipped() || (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0))) {
	      logic.name_tag = mname.str() + xml_flipped;
	      logic.shape_tag = nspace + ":" + logic.name_tag;
	      l.push_back(logic);
	    }
            logic.name_tag = mname.str();
            logic.shape_tag = nspace + ":" + logic.name_tag;
	    l.push_back(logic);

            // name_tag is BModule1Layer1 and it goes into all files
	    if (!iiter->getModule().flipped()) pos.child_tag = nspace + ":" + mname.str();
            else pos.child_tag = nspace + ":" + mname.str() + xml_flipped;
	    
	    std::string xml_tilted_mod_rot;
            if (lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty()) {pos.rotref = nspace + ":" + xml_default_mod_rot;}
	    else {
	      if (tiltAngle == 0) {pos.rotref = nspace + ":" + xml_default_mod_rot;}
	      else {
		pconverter << xml_default_mod_rot << "_tilt_" << tiltAngle;
		xml_tilted_mod_rot = pconverter.str();
		pconverter.str("");
	      
		rot.name = xml_tilted_mod_rot + "_" + xml_plus;
		if (r.count(rot.name) == 0 ) {
		  rot.thetax = 90.0;
		  rot.phix = 90.0;
		  rot.thetay = tiltAngle;
		  rot.phiy = 180.0;
		  rot.thetaz = 90.0 - tiltAngle;
		  rot.phiz = 0.0;
		  r.insert(std::pair<const std::string,Rotation>(rot.name,rot));
		}
		  
		rot.name = xml_tilted_mod_rot + "_" + xml_minus;
		if (r.count(rot.name) == 0 ) {
		  rot.thetax = 90.0;
		  rot.phix = 90.0;
		  rot.thetay = tiltAngle;
		  rot.phiy = 0.0;
		  rot.thetaz = 90.0 + tiltAngle;
		  rot.phiz = 0.0;
		  r.insert(std::pair<const std::string,Rotation>(rot.name,rot));
		}
	      }
	    }
	    //if (logic.shape_tag == "tracker:BModule5Layer1") { pos.rotref = nspace + ":" + xml_barrel_tilt_47; }

	    // TODO: comment this lime
            if ((iiter->getModule().center().Rho() > (rmax - deltar / 2.0))
                || ((iiter->getModule().center().Rho() < ((rmin + rmax) / 2.0))
                    && (iiter->getModule().center().Rho() > (rmin + deltar / 2.0)))) {
	      pos.trans.dx = deltar / 2.0 - shape.dz; 
	    }
	    else pos.trans.dx = shape.dz - deltar / 2.0;
	    // debug cout << the values here for a known rod

            if (is_short) {
              pos.parent_tag = nspace + ":" + rodname.str() + xml_plus;
              pos.trans.dz = iiter->getModule().minZ() - ((zmax + zmin) / 2.0) + shape.dy;
	      if (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0)) {
		pos.trans.dx = 0;
		pos.trans.dz = 0;
		pos.rotref = nspace + ":" + xml_tilted_mod_rot + "_" + xml_plus;
		pos.parent_tag =  nspace + ":" + mname.str() + xml_tilted + xml_plus + xml_flipped;
		pos.child_tag = nspace + ":" + mname.str() + xml_flipped;
		p.push_back(pos);
		pos.parent_tag = nspace + ":" + mname.str() + xml_tilted + xml_plus;
		pos.child_tag = nspace + ":" + mname.str();
	      }
	      p.push_back(pos);
              pos.parent_tag = nspace + ":" + rodname.str() + xml_minus;
              pos.trans.dz = -pos.trans.dz;
	      pos.copy = 2; // This is a copy of the BModule (FW/BW barrel half)
	      if (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0)) {
		pos.trans.dx = 0;
		pos.trans.dz = 0;
		pos.rotref = nspace + ":" + xml_tilted_mod_rot + "_" + xml_minus;
		pos.parent_tag = nspace + ":" + mname.str() + xml_tilted + xml_minus + xml_flipped;
		pos.child_tag = nspace + ":" + mname.str() + xml_flipped;
		p.push_back(pos);
		pos.parent_tag = nspace + ":" + mname.str() + xml_tilted + xml_minus;
		pos.child_tag = nspace + ":" + mname.str();
	      }
	      p.push_back(pos);
	      pos.copy = 1;
	      if (iiter->getModule().flipped() || (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0))) {
		pos.trans.dx = 0;
		pos.trans.dz = 0;
		pos.rotref = nspace + ":" + xml_flip_mod_rot;
		pos.parent_tag = nspace + ":" + mname.str() + xml_flipped;
		pos.child_tag = nspace + ":" + mname.str();
		p.push_back(pos);
	      }
            }
	    
	    else {
              pos.parent_tag = nspace + ":" + rodname.str();
              partner = findPartnerModule(iiter, iguard, modRing);
	      pos.trans.dz = iiter->getModule().maxZ() - shape.dy;
	      if (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0)) {
		pos.trans.dx = 0;
		pos.trans.dz = 0;
		pos.rotref = nspace + ":" + xml_tilted_mod_rot + "_" + xml_plus;
		pos.parent_tag = nspace + ":" + mname.str() + xml_tilted + xml_plus + xml_flipped;
		pos.child_tag = nspace + ":" + mname.str() + xml_flipped;
		p.push_back(pos);
		pos.parent_tag = nspace + ":" + mname.str() + xml_tilted + xml_plus;
		pos.child_tag = nspace + ":" + mname.str();
	      }
	      p.push_back(pos);
 
	      if (partner != iguard) {
		if ((partner->getModule().center().Rho() > (rmax - deltar / 2.0))
		    || ((partner->getModule().center().Rho() < ((rmin + rmax) / 2.0))
			&& (partner->getModule().center().Rho() > (rmin + deltar / 2.0)))) pos.trans.dx = deltar / 2.0 - shape.dz;
		else pos.trans.dx = shape.dz - deltar / 2.0;
		pos.trans.dz = partner->getModule().maxZ() - shape.dy;
		pos.copy = 2; // This is a copy of the BModule (FW/BW barrel half)
		if (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0)) {
		  pos.trans.dx = 0;
		  pos.trans.dz = 0;
		  pos.rotref = nspace + ":" + xml_tilted_mod_rot + "_" + xml_minus;
		  pos.parent_tag = nspace + ":" + mname.str() + xml_tilted + xml_minus + xml_flipped;
		  pos.child_tag = nspace + ":" + mname.str() + xml_flipped;
		  p.push_back(pos);
		  pos.parent_tag = nspace + ":" + mname.str() + xml_tilted + xml_minus;
		  pos.child_tag = nspace + ":" + mname.str();
		}
		p.push_back(pos);
		pos.copy = 1;
	      }

	      if (iiter->getModule().flipped() || (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty() && (tiltAngle != 0))) {
		pos.trans.dx = 0;
		pos.trans.dz = 0;
		pos.rotref = nspace + ":" + xml_flip_mod_rot;
		pos.parent_tag = nspace + ":" + mname.str() + xml_flipped;
		pos.child_tag = nspace + ":" + mname.str();
		p.push_back(pos);
	      }
            }
            pos.rotref = "";

            // wafer
            string xml_base_inout = "";
            if (iiter->getModule().numSensors() == 2) xml_base_inout = xml_base_inner;

            shape.name_tag = mname.str() + xml_base_inout + xml_base_waf;
            shape.dz = iiter->getModule().sensorThickness() / 2.0;
            //if (iiter->getModule().numSensors() == 2) shape.dz = shape.dz / 2.0; // CUIDADO calcSensThick returned 2x what getSensThick returns, it means that now one-sided sensors are half as thick if not compensated for in the config files
            s.push_back(shape);

            pos.parent_tag = nspace + ":" + mname.str();

            logic.name_tag = mname.str() + xml_base_inout + xml_base_waf;
            logic.shape_tag = nspace + ":" + logic.name_tag;
            logic.material_tag = xml_material_air;
            l.push_back(logic);

            pos.child_tag = nspace + ":" + mname.str() + xml_base_inout + xml_base_waf;
            pos.trans.dx = 0.0;
            pos.trans.dz = /*shape.dz*/ - iiter->getModule().dsDistance() / 2.0; 
            p.push_back(pos);

            if (iiter->getModule().numSensors() == 2) {

              xml_base_inout = xml_base_outer;

              shape.name_tag = mname.str() + xml_base_inout + xml_base_waf;
              s.push_back(shape);

              logic.name_tag = mname.str() + xml_base_inout + xml_base_waf;
              logic.shape_tag = nspace + ":" + logic.name_tag;
              l.push_back(logic);

              pos.child_tag = nspace + ":" + mname.str() + xml_base_inout + xml_base_waf;
              pos.trans.dz = pos.trans.dz + /*2 * shape.dz +*/ iiter->getModule().dsDistance();  // CUIDADO: was with 2*shape.dz, but why???
              //pos.copy = 2;

              if (iiter->getModule().stereoRotation() != 0) {
                rot.name = type_stereo + mname.str();
                rot.thetax = 90.0;
                rot.phix = iiter->getModule().stereoRotation() / M_PI * 180;
                rot.thetay = 90.0;
                rot.phiy = 90.0 + iiter->getModule().stereoRotation() / M_PI * 180;
                r.insert(std::pair<const std::string,Rotation>(rot.name,rot));
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

            shape.name_tag = mname.str() + xml_base_inout + xml_base_act;
            s.push_back(shape);

            //pos.parent_tag = logic.shape_tag;
            pos.parent_tag = nspace + ":" + mname.str() + xml_base_inout + xml_base_waf;

            logic.name_tag = mname.str() + xml_base_inout + xml_base_act;
            logic.shape_tag = nspace + ":" + logic.name_tag;
            logic.material_tag = nspace + ":" + xml_sensor_silicon;
            l.push_back(logic);

            pos.child_tag = nspace + ":" + mname.str() + xml_base_inout + xml_base_act;
            pos.trans.dz = 0.0;
#ifdef __FLIPSENSORS_IN__ // Flip INNER sensors
            pos.rotref = nspace + ":" + rot_sensor_tag;
#endif
            p.push_back(pos);

            // topology
            mspec.partselectors.push_back(mname.str() + xml_base_inout + xml_base_act);

            minfo.name		= iiter->getModule().moduleType();
            minfo.rocrows	= any2str<int>(iiter->getModule().innerSensor().numROCRows());  // in case of single sensor module innerSensor() and outerSensor() point to the same sensor
            minfo.roccols	= any2str<int>(iiter->getModule().innerSensor().numROCCols());
            minfo.rocx		= any2str<int>(iiter->getModule().innerSensor().numROCX());
            minfo.rocy		= any2str<int>(iiter->getModule().innerSensor().numROCY());

            mspec.moduletypes.push_back(minfo);

            if (iiter->getModule().numSensors() == 2) { 

              // active surface
              xml_base_inout = xml_base_outer;

              shape.name_tag = mname.str() + xml_base_inout + xml_base_act;
              s.push_back(shape);

              pos.parent_tag = nspace + ":" + mname.str() + xml_base_inout + xml_base_waf;

              logic.name_tag = mname.str() + xml_base_inout + xml_base_act;
              logic.shape_tag = nspace + ":" + logic.name_tag;
              l.push_back(logic);

              pos.child_tag = nspace + ":" + mname.str() + xml_base_inout + xml_base_act;
#ifdef __FLIPSENSORS_OUT__ // Flip OUTER sensors
              pos.rotref = nspace + ":" + rot_sensor_tag;
#endif
              p.push_back(pos);

              // topology
              mspec.partselectors.push_back(mname.str() + xml_base_inout + xml_base_act);

              minfo.rocrows	= any2str<int>(iiter->getModule().outerSensor().numROCRows());
              minfo.roccols	= any2str<int>(iiter->getModule().outerSensor().numROCCols());
              minfo.rocx	= any2str<int>(iiter->getModule().outerSensor().numROCX());
              minfo.rocy	= any2str<int>(iiter->getModule().outerSensor().numROCY());

              mspec.moduletypes.push_back(minfo);
#ifdef __USE_MODULECOMPLEX__ 
              modcomplex.addMaterialInfo(c);
              modcomplex.addShapeInfo(s);
              modcomplex.addLogicInfo(l);
              modcomplex.addPositionInfo(p);
#ifdef __DEBUGPRINT__
              modcomplex.print();
#endif
#endif
            } // End of replica for Pt-modules

            // material properties
            rtotal = rtotal + iiter->getRadiationLength();
            itotal = itotal + iiter->getInteractionLength();
            count++;
            dt = iiter->getModule().thickness();
          }

	    if (iiter->getModule().uniRef().phi == 2 && !lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty()) {
	      if (tiltAngle == 0) { 
		flatzmax = iiter->getModule().center().Z() + iiter->getModule().thickness() / 2.0;
	      }
	      else {
		std::map<int,BTiltedRingInfo>::iterator it;
		it = rinfoplus.find(modRing);
		if (it != rinfoplus.end()) {
		  it->second.r2 = iiter->getModule().center().Rho();
		  it->second.z2 = iiter->getModule().center().Z();
	          it->second.rmax = modcomplex.getRmax();
	          it->second.zmax = modcomplex.getZmax();
		}
		it = rinfominus.find(modRing);
		if (it != rinfominus.end()) {
		  it->second.r2 = iiter->getModule().center().Rho();
		  it->second.z2 = - iiter->getModule().center().Z();
	          it->second.rmax = modcomplex.getRmax();
	          it->second.zmax = modcomplex.getZmax();
		}
	      }
	    }

	  }
	}
        if (count > 0) {
          ril.rlength = rtotal / (double)count;
          ril.ilength = itotal / (double)count;
          ri.push_back(ril);
        }

	ds = fromRim(rmax, shape.dx);

        // rod(s)
        shape.name_tag = rodname.str();
        if (is_short) shape.name_tag = shape.name_tag + xml_plus;
        shape.dy = shape.dx;
        shape.dx = rodThickness / 2.0;
        if (is_short) shape.dz = (zmax - zmin) / 2.0;
        else shape.dz = zmax;
	if (!lagg.getBarrelLayers()->at(layer - 1)->tiltedLayerSpecFile().empty()) { shape.dz = flatzmax; }
        s.push_back(shape);
        logic.name_tag = rodname.str();
        logic.shape_tag = nspace + ":" + logic.name_tag;
        logic.material_tag = xml_material_air;
        l.push_back(logic);
        rspec.partselectors.push_back(rodname.str());
        rspec.moduletypes.push_back(minfo_zero);
        //rspec.moduletypes.push_back("");
        if (is_short) {
          shape.name_tag = rodname.str() + xml_minus;
          s.push_back(shape);
          logic.name_tag = rodname.str() + xml_minus;
          logic.shape_tag = nspace + ":" + logic.name_tag;
          l.push_back(logic);
          rspec.partselectors.push_back(rodname.str() + xml_minus);
          //rspec.moduletypes.push_back("");
          rspec.moduletypes.push_back(minfo_zero);
        }
        // extra containers for short layers
        if (wt && is_short) {
          shape.type = tb;
          shape.dx = 0.0;
          shape.dy = 0.0;
          shape.rmin = rmin;
          shape.rmax = rmax;
          shape.name_tag = lname.str() + xml_plus;
          s.push_back(shape);
          logic.name_tag = lname.str() + xml_plus;
          logic.shape_tag = nspace + ":" + logic.name_tag;
          l.push_back(logic);
          pos.trans.dx = 0.0;
          pos.trans.dz = zmin + (zmax - zmin) / 2.0;
          pos.parent_tag = nspace + ":" + lname.str();
          pos.child_tag = nspace + ":" + lname.str() + xml_plus;
          p.push_back(pos);
          shape.name_tag = lname.str() + xml_minus;
          s.push_back(shape);
          logic.name_tag = lname.str() + xml_minus;
          logic.shape_tag = nspace + ":" + logic.name_tag;
          l.push_back(logic);
          pos.trans.dz = - pos.trans.dz;
          pos.child_tag =  nspace + ":" + lname.str() + xml_minus;
          p.push_back(pos);
        }
        // rods in layer algorithm(s)
	alg.name = xml_tobalgo;
        alg.parent = nspace + ":" + lname.str();
        if (wt && is_short) alg.parent = alg.parent + xml_plus;
	pconverter <<  nspace + ":" + rodname.str();
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
	pconverter.str("");
        alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
        alg.parameters.push_back(numericParam(xml_incrcopyno, "1"));
        a.push_back(alg);
        // extras for short layers
        if (is_short) {
          if (wt) alg.parent = nspace + ":" + lname.str() + xml_minus;
          pconverter << nspace << ":" << rodname.str() << xml_minus;
          alg.parameters.front() = stringParam(xml_childparam, pconverter.str());
          pconverter.str("");
          if (wt) pconverter << "0.0*mm";
          else pconverter << -(zmin + (zmax - zmin) / 2.0) << "*mm";
          alg.parameters.at(6) = numericParam(xml_zposition, pconverter.str());
	  pconverter.str("");
          a.push_back(alg);
        }
        alg.parameters.clear();

        // tilted rings
	if ( !rinfoplus.empty() || !rinfominus.empty() ) {

	  std::map<std::string, std::map<int,BTiltedRingInfo>> rinfototal;
	  if ( !rinfoplus.empty() ) { rinfototal.insert({"rinfoplus", rinfoplus}); }
	  if ( !rinfominus.empty() ) { rinfototal.insert({"rinfominus", rinfominus}); }

	  for (auto const &rinfoside : rinfototal) {
	    shape.type = tb;
	    shape.dx = 0.0;
	    shape.dy = 0.0;
	    shape.dyy = 0.0;
	    pos.trans.dx = 0;
	    pos.trans.dy = 0;
	    pos.trans.dz = 0;

	    for (auto const &ringinfo : rinfoside.second) {
	      auto const& rinfo = ringinfo.second;
	      if (rinfo.modules > 0) {
		shape.name_tag = rinfo.name;
#ifndef __USE_MODULECOMPLEX__ 
		shape.rmin = rinfo.r1 - sqrt( pow(rinfo.mdy , 2) + pow(rinfo.mdz , 2)) * sin(rinfo.tiltAngle * M_PI / 180 + atan(rinfo.mdz / rinfo.mdy));
		shape.rmax = rinfo.r2 + sqrt( pow(rinfo.mdy , 2) + pow(rinfo.mdz , 2)) * sin(rinfo.tiltAngle * M_PI / 180 + atan(rinfo.mdz / rinfo.mdy));
		shape.dz = abs(rinfo.z2 - rinfo.z1) / 2 + sqrt( pow(rinfo.mdy , 2) + pow(rinfo.mdz , 2)) * cos(rinfo.tiltAngle * M_PI / 180 - atan(rinfo.mdz / rinfo.mdy));
#else
                shape.rmin = rinfo.rmin;
                shape.rmax = rinfo.rmax;
                shape.dz = (rinfo.zmax-rinfo.zmin)/2.;
#endif
		s.push_back(shape);

		logic.name_tag = rinfo.name;
		logic.shape_tag = nspace + ":" + logic.name_tag;
		logic.material_tag = xml_material_air;
		l.push_back(logic);
	      
		pos.parent_tag = nspace + ":" + lname.str();
		pos.child_tag = nspace + ":" + rinfo.name;
	      
		//if (rinfo.fw) pos.trans.dz = (zmin - zmax) / 2.0 + shape.dz;
		//else pos.trans.dz = (zmax - zmin) / 2.0 - shape.dz;
		pos.trans.dz = (rinfo.z1 + rinfo.z2) / 2; 
		p.push_back(pos);
	      
		rspec.partselectors.push_back(rinfo.name);
		//rspec.moduletypes.push_back(minfo_zero);
	      
		alg.name = xml_tilted_rings_algo;
		alg.parent = nspace + ":" + rinfo.name;
		if (rinfo.inner_flipped) {
		  alg.parameters.push_back(stringParam(xml_childparam, nspace + ":" + rinfo.flippedchildname));
		} else {
		  alg.parameters.push_back(stringParam(xml_childparam, nspace + ":" + rinfo.childname));
		}
		pconverter << (rinfo.modules / 2);
		alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
		pconverter.str("");
		alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
		alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
		alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
		pconverter << 2 * M_PI / (double)(rinfo.modules) * (rinfo.phi - 1);
		alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
		pconverter.str("");
		pconverter << rinfo.r1;
		alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
		pconverter.str("");
		alg.parameters.push_back(vectorParam(0, 0, (rinfo.z1 - rinfo.z2) / 2.0));
		a.push_back(alg);
		alg.parameters.clear();
	      
		alg.name = xml_tilted_rings_algo;
		alg.parent = logic.shape_tag;
		if (rinfo.inner_flipped) {
		  alg.parameters.push_back(stringParam(xml_childparam, nspace + ":" + rinfo.childname));
		} else {
		  alg.parameters.push_back(stringParam(xml_childparam, nspace + ":" + rinfo.flippedchildname));
		}
		pconverter << (rinfo.modules / 2);
		alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
		pconverter.str("");
		alg.parameters.push_back(numericParam(xml_startcopyno, "2"));
		alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
		alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
		pconverter << 2 * M_PI / (double)(rinfo.modules) * (rinfo.phi);
		alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
		pconverter.str("");
		pconverter << rinfo.r2;
		alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
		pconverter.str("");
		alg.parameters.push_back(vectorParam(0, 0, (rinfo.z2 - rinfo.z1) / 2.0));
		a.push_back(alg);
		alg.parameters.clear();
	      }
	    }
	  }
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
        logic.name_tag = lname.str();
        logic.shape_tag = nspace + ":" + logic.name_tag;
        l.push_back(logic);
        pos.parent_tag = xml_pixbarident + ":" + xml_2OTbar;
        pos.child_tag = nspace + ":" + lname.str();
        p.push_back(pos);
        lspec.partselectors.push_back(lname.str());
        //lspec.moduletypes.push_back("");
        lspec.moduletypes.push_back(minfo_zero);
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
                               std::vector<AlgoInfo>& a, std::map<std::string,Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, bool wt) {
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
        std::map<int, ERingInfo> rinfo;

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
        //zmax = zmax + v.max / 2.0; 
        double zmin = lagg.getEndcapLayers()->at(layer - 1)->minZ();
        //zmin = zmin - v.max / 2.0;
        double ringThickness = lagg.getEndcapLayers()->at(layer - 1)->maxRingThickness(); // all the ring volumes will have the same thickness (equal to the thickest ring in the disk)  
        double diskThickness = lagg.getEndcapLayers()->at(layer - 1)->thickness();

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

#ifndef __USE_MODULECOMPLEX__
            // module composite material
            matname << xml_base_actcomp << "D" << layer << "R" << modRing;
            c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));
#endif

            rname << xml_ring << modRing << dname.str();
            mname << xml_endcap_module << modRing << dname.str();

#ifdef __USE_MODULECOMPLEX__
            if (iiter->getModule().flipped()) { // flip case
              // We will add an additional hierarchy to rotate this module.
              // Here we will take a naming scheme e.g.:
              // "EModule1Disc6Flipped" for the outer shell,
              // "EModule1Disc6" to be added contents.
              shape.name_tag = mname.str(); 

              // Define the outer shell
              logic.name_tag = mname.str();
              logic.shape_tag = nspace + ":" + mname.str() + xml_flipped;
              logic.material_tag = xml_material_air;
              pos.parent_tag = nspace + ":" + mname.str() + xml_flipped;
              pos.child_tag = nspace + ":" + shape.name_tag;
              pos.rotref = nspace + ":" + xml_flip_mod_rot;
              l.push_back(logic);
              p.push_back(pos);
              pos.rotref.clear();
            } else {
              shape.name_tag = mname.str();
            }
            std::string parentName = shape.name_tag; // e.g. EModule1Disc6 / EModule1Disc6Flipped
            ModuleComplex modcomplex(mname.str(),parentName,*iiter);
            modcomplex.buildSubVolumes();

#endif

            // collect ring info
            ERingInfo rinf;
            rinf.name = rname.str();
            rinf.childname = mname.str();
            rinf.fw = (iiter->getModule().center().Z() < (zmin + zmax) / 2.0);
            //rinf.modules = lagg.getEndcapLayers()->at(layer - 1)->rings().at(modRing-1).numModules();
            //rinf.modules = lagg.getEndcapLayers()->at(layer - 1)->rings().at(modRing).numModules();
            rinf.modules = lagg.getEndcapLayers()->at(layer - 1)->ringsMap().at(modRing)->numModules();

#ifndef __USE_MODULECOMPLEX__
            rinf.rin = iiter->getModule().minR();
            rinf.rout = iiter->getModule().maxR();
#else
            rinf.rin  = modcomplex.getRmin();
            rinf.rout = modcomplex.getRmax();
#endif
            rinf.rmid = iiter->getModule().center().Rho();
            rinf.mthk = iiter->getModule().thickness();
            rinf.phi = iiter->getModule().center().Phi();
            rinfo.insert(std::pair<int, ERingInfo>(modRing, rinf));

            // module trapezoid

            shape.type = iiter->getModule().shape() == RECTANGULAR ? bx : tp;

            shape.name_tag = mname.str();
#ifndef __USE_MODULECOMPLEX__
            shape.dx = iiter->getModule().minWidth() / 2.0;
            shape.dxx = iiter->getModule().maxWidth() / 2.0;
            shape.dy = iiter->getModule().length() / 2.0;
            shape.dyy = iiter->getModule().length() / 2.0;
            shape.dz = iiter->getModule().thickness() / 2.0;
#else       
            if (shape.type==bx) {
              shape.dx = modcomplex.getExpandedModuleWidth()/2.0;
              shape.dy = modcomplex.getExpandedModuleLength()/2.0;
              shape.dz = modcomplex.getExpandedModuleThickness()/2.0;
            } else { // obsolete?
              shape.dx = iiter->getModule().minWidth() / 2.0 + iiter->getModule().serviceHybridWidth();
              shape.dxx = iiter->getModule().maxWidth() / 2.0 + iiter->getModule().serviceHybridWidth();
              shape.dy = iiter->getModule().length() / 2.0 + iiter->getModule().frontEndHybridWidth();
              shape.dyy = iiter->getModule().length() / 2.0 + iiter->getModule().frontEndHybridWidth();
              shape.dz = iiter->getModule().thickness() / 2.0 + iiter->getModule().supportPlateThickness();
            }
#endif
            s.push_back(shape);
#ifdef __USE_MODULECOMPLEX__ 
            // Get it back for sensors
            shape.dx = iiter->getModule().minWidth() / 2.0;
            shape.dxx = iiter->getModule().maxWidth() / 2.0;
            shape.dy = iiter->getModule().length() / 2.0;
            shape.dyy = iiter->getModule().length() / 2.0;
            shape.dz = iiter->getModule().thickness() / 2.0;
#endif

            logic.name_tag = shape.name_tag;
            logic.shape_tag = nspace + ":" + logic.name_tag;
#if 0
            logic.material_tag = nspace + ":" + matname.str();
#else
            logic.material_tag = xml_material_air;
#endif
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
                r.insert(std::pair<const std::string,Rotation>(rot.name,rot));
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
#ifdef __FLIPSENSORS_IN__ // Flip INNER sensors
            pos.rotref = nspace + ":" + rot_sensor_tag;
#endif
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
#ifdef __FLIPSENSORS_OUT__ // Flip OUTER sensors
              pos.rotref = nspace + ":" + rot_sensor_tag;
#endif
              p.push_back(pos);

              // topology
              mspec.partselectors.push_back(logic.name_tag);

              minfo.rocrows	= any2str<int>(iiter->getModule().outerSensor().numROCRows());
              minfo.roccols	= any2str<int>(iiter->getModule().outerSensor().numROCCols());
              minfo.rocx		= any2str<int>(iiter->getModule().outerSensor().numROCX());
              minfo.rocy		= any2str<int>(iiter->getModule().outerSensor().numROCY());

              mspec.moduletypes.push_back(minfo);
              //mspec.moduletypes.push_back(iiter->getModule().getType());
#ifdef __USE_MODULECOMPLEX__ 
              modcomplex.addMaterialInfo(c);
              modcomplex.addShapeInfo(s);
              modcomplex.addLogicInfo(l);
              modcomplex.addPositionInfo(p);
#ifdef __DEBUGPRINT__
              modcomplex.print();
#endif
#endif
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
        shape.dz = ringThickness / 2.0; //findDeltaZ(lagg.getEndcapLayers()->at(layer - 1)->getModuleVector()->begin(), // CUIDADO what the hell is this??
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
            pconverter << (rinfo[*siter].phi + 2 * M_PI / (double)(rinfo[*siter].modules));
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
        shape.dz = diskThickness/2.0; //(zmax - zmin) / 2.0;
        s.push_back(shape);

        logic.name_tag = shape.name_tag; // CUIDADO ended with + xml_plus;
        //logic.extra = xml_plus;
        logic.shape_tag = nspace + ":" + shape.name_tag;
        logic.material_tag = xml_material_air;
        l.push_back(logic);

        pos.parent_tag = xml_pixfwdident + ":" + xml_2OTfwd;
        pos.child_tag = nspace + ":" + logic.name_tag;
        pos.trans.dz = (zmax + zmin) / 2.0 - xml_z_pixfwd;
        p.push_back(pos);

        dspec.partselectors.push_back(logic.name_tag);
        dspec.moduletypes.push_back(minfo_zero);
        dspec.partextras.push_back(logic.extra);
        //   logic.name_tag = shape.name_tag; // CUIDADO ended with + xml_minus;
        //   logic.extra = xml_minus;
        //   l.push_back(logic);
        //   pos.parent_tag = xml_pixfwdident + ":" + xml_2OTfwd;
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
#if 1
    int  previousInnerRadius = -1;
#endif
    for (iter = bs.begin(); iter != guard; iter++) {
#if 1
      if ( (int)(iter->getZOffset()) == 0 ) {
        if ( previousInnerRadius == (int)(iter->getInnerRadius()) ) continue;  
        else previousInnerRadius = (int)(iter->getInnerRadius());
      }
#endif
      std::ostringstream matname, shapename;
#if 0
      matname << xml_base_serfcomp << iter->getCategory() << "R" << (int)(iter->getInnerRadius()) << "dZ" << (int)(iter->getZLength());
      shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(iter->getZOffset());
      if ((iter->getZOffset() + iter->getZLength()) > 0) c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
#else
      matname << xml_base_serfcomp << iter->getCategory() << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
      shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
      if ((iter->getZOffset() + iter->getZLength()) > 0 ) {
        if ( iter->getLocalMasses().size() ) {
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

          pos.parent_tag = xml_pixbarident + ":" + xml_2OTbar; //xml_tracker;
          pos.child_tag = logic.shape_tag;
          pos.trans.dz = iter->getZOffset() + shape.dz;
          p.push_back(pos);
        
        } else {
          std::stringstream msg;
          msg << shapename.str() << " is not exported to XML because it is empty." << std::ends;
          logWARNING( msg.str() ); 
        }
      }
#endif
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
#if 0
      if ((iter->getZOffset() + iter->getZLength()) > 0) { // This is necessary because of replication of Forward volumes!
        c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
#else
      if ( (iter->getZOffset() + iter->getZLength()) > 0 ) { 
        if ( iter->getLocalMasses().size() ) {
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

          pos.parent_tag = xml_pixfwdident + ":" + xml_2OTfwd; // xml_tracker;
          pos.child_tag = logic.shape_tag;
          pos.trans.dz = iter->getZOffset() + shape.dz;
          p.push_back(pos);
        }
        else {
          std::stringstream msg;
          msg << shapename.str() << " is not exported to XML because it is empty." << std::ends;
          logWARNING( msg.str() ); 
        }
#endif
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
#if 0
      shapename << xml_base_lazy /*<< any2str(iter->getCategory()) */<< "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset()));
#else
      shapename << xml_base_lazy /*<< any2str(iter->getCategory()) */<< "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(iter->getZLength() / 2.0 + iter->getZOffset());
#endif

      fres = found.find(iter->getCategory());
#if 0
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

      switch (iter->getCategory()) {
      case MaterialProperties::b_sup:
      case MaterialProperties::t_sup:
      case MaterialProperties::u_sup:
      case MaterialProperties::o_sup:
          pos.parent_tag = xml_pixbarident + ":" + xml_2OTbar;
          break;
      case MaterialProperties::e_sup:
          pos.parent_tag = xml_pixfwdident + ":" + xml_2OTfwd;
          break;
      default:
          pos.parent_tag = nspace + ":" + xml_tracker;
      }
      pos.child_tag = logic.shape_tag;
      if ((iter->getCategory() == MaterialProperties::o_sup) ||
          (iter->getCategory() == MaterialProperties::t_sup)) pos.trans.dz = 0.0;
      else pos.trans.dz = iter->getZOffset() + shape.dz;
      p.push_back(pos);
#else
      if (fres == found.end() && iter->getLocalMasses().size() ) { 
        c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
        found.insert(iter->getCategory());

        shape.name_tag = shapename.str();
        shape.dz = iter->getZLength() / 2.0;
        shape.rmin = iter->getInnerRadius();
        shape.rmax = shape.rmin + iter->getRWidth();
        s.push_back(shape);

        logic.name_tag = shapename.str();
        logic.shape_tag = nspace + ":" + shapename.str();
        logic.material_tag = nspace + ":" + matname.str();
        l.push_back(logic);

        switch (iter->getCategory()) {
        case MaterialProperties::b_sup:
        case MaterialProperties::t_sup:
        case MaterialProperties::u_sup:
        case MaterialProperties::o_sup:
            pos.parent_tag = xml_pixbarident + ":" + xml_2OTbar;
            break;
        case MaterialProperties::e_sup:
            pos.parent_tag = xml_pixfwdident + ":" + xml_2OTfwd;
            break;
        default:
            pos.parent_tag = nspace + ":" + xml_tracker;
        }
        pos.child_tag = logic.shape_tag;
        if ((iter->getCategory() == MaterialProperties::o_sup) ||
            (iter->getCategory() == MaterialProperties::t_sup)) pos.trans.dz = 0.0;
        else pos.trans.dz = iter->getZOffset() + shape.dz;
        p.push_back(pos);
      }
#endif
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
    d = 1000 * ie.getTotalMass() / (M_PI * ie.getZLength() * d);
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

#ifdef __USE_MODULECOMPLEX__ 
  // These value should be consistent with 
  // the configuration file
  const int ModuleComplex::HybridFBLR_0  = 0; // Front + Back + Right + Left
  const int ModuleComplex::InnerSensor   = 1; 
  const int ModuleComplex::OuterSensor   = 2; 
  const int ModuleComplex::HybridFront   = 3; 
  const int ModuleComplex::HybridBack    = 4; 
  const int ModuleComplex::HybridLeft    = 5; 
  const int ModuleComplex::HybridRight   = 6; 
  const int ModuleComplex::HybridBetween = 7; 
  const int ModuleComplex::SupportPlate  = 8; // Support Plate 
  const int ModuleComplex::nTypes        = 9; 
  // extras
  const int ModuleComplex::HybridFB        = 34; 
  const int ModuleComplex::HybridLR        = 56; 
  const int ModuleComplex::HybridFBLR_3456 = 3456; 

  const double ModuleComplex::kmm3Tocm3 = 1e-3; 

  ModuleComplex::ModuleComplex(std::string moduleName,
                               std::string parentName,
                               ModuleCap&  modcap        ) : moduleId(moduleName),
                                                             parentId(parentName),
                                                             modulecap(modcap),
                                                             module(modcap.getModule()),
                                                             modWidth(module.area()/module.length()),
                                                             modLength(module.length()),
                                                             modThickness(module.thickness()),
                                                             sensorThickness(module.sensorThickness()),
                                                             sensorDistance(module.dsDistance()),
                                                             frontEndHybridWidth(module.frontEndHybridWidth()),
                                                             serviceHybridWidth(module.serviceHybridWidth()),
                                                             hybridThickness(module.hybridThickness()),
                                                             supportPlateThickness(module.supportPlateThickness()),
                                                             hybridTotalMass(0.),
                                                             hybridTotalVolume_mm3(-1.),
                                                             hybridFrontAndBackVolume_mm3(-1.),
                                                             hybridLeftAndRightVolume_mm3(-1.),
                                                             moduleMassWithoutSensors_expected(0.),
                                                             expandedModWidth(modWidth+2*serviceHybridWidth),
                                                             expandedModLength(modLength+2*frontEndHybridWidth),
                                                             expandedModThickness(sensorDistance+2*(supportPlateThickness+sensorThickness)),
                                                             center(module.center()),
                                                             normal(module.normal()),
                                                             prefix_xmlfile("tracker:"),
                                                             prefix_material("hybridcomposite") {
  }

  ModuleComplex::~ModuleComplex() {
    std::vector<Volume*>::const_iterator vit;
    for ( vit = volumes.begin(); vit != volumes.end(); vit++ ) delete *vit;
  }

  void ModuleComplex::buildSubVolumes() {
  //  Top View
  //  ------------------------------
  //  |            L(5)            |  
  //  |----------------------------|     y
  //  |     |                |     |     ^
  //  |B(4) |     Between    | F(3)|     |
  //  |     |       (7)      |     |     +----> x
  //  |----------------------------|
  //  |            R(6)            |     
  //  ------------------------------     
  //                                            z
  //  Side View                                 ^
  //         ---------------- OuterSensor(2)    |
  //  ====== ================ ====== Hybrids    +----> x
  //         ---------------- InnerSensor(1)
  //  ============================== 
  //          SupportPlate(8)                      
  //
  //  R(6) and L(5) are Front-End Hybrids
  //  B(4) and F(3) are Service Hybdrids
  //
  //
    Volume* vol[nTypes];
    //Unused pointers
    vol[HybridFBLR_0] = 0;
    vol[InnerSensor]  = 0;
    vol[OuterSensor]  = 0;

    double dx = serviceHybridWidth;              
    double dy = modLength; 
    double dz = hybridThickness;  
    double posx = (modWidth+serviceHybridWidth)/2.;
    double posy = 0.;
    double posz = 0.;
    // Hybrid FrontSide Volume
    vol[HybridFront] = new Volume(moduleId+"FSide",HybridFront,parentId,dx,dy,dz,posx,posy,posz);

    posx = -(modWidth+serviceHybridWidth)/2.;
    posy = 0.;
    posz = 0.;
    // Hybrid BackSide Volume
    vol[HybridBack] = new Volume(moduleId+"BSide",HybridBack,parentId,dx,dy,dz,posx,posy,posz);

    dx = modWidth+2*serviceHybridWidth;  
    dy = frontEndHybridWidth;
    posx = 0.;
    posy = (modLength+frontEndHybridWidth)/2.;
    posz = 0.;
    // Hybrid LeftSide Volume
    vol[HybridLeft] = new Volume(moduleId+"LSide",HybridLeft,parentId,dx,dy,dz,posx,posy,posz);

    posx = 0.;
    posy = -(modLength+frontEndHybridWidth)/2.;
    posz = 0.;
    // Hybrid RightSide Volume
    vol[HybridRight] = new Volume(moduleId+"RSide",HybridRight,parentId,dx,dy,dz,posx,posy,posz);

    dx = modWidth; 
    dy = modLength; 
    posx = 0.;
    posy = 0.;
    posz = 0.;
    // Hybrid Between Volume
    vol[HybridBetween] = new Volume(moduleId+"Between",HybridBetween,parentId,dx,dy,dz,posx,posy,posz);

    dx = expandedModWidth;  
    dy = expandedModLength; 
    dz = supportPlateThickness; 
    posx = 0.;
    posy = 0.;
    posz = - ( ( sensorDistance + supportPlateThickness )/2. + sensorThickness ); 
    // SupportPlate
    vol[SupportPlate] = new Volume(moduleId+"SupportPlate",SupportPlate,parentId,dx,dy,dz,posx,posy,posz);

    // ===========================================================
    // Finding Rmin/Rmax/Zmin/Zmax with considering hybrid volumes
    // ===========================================================
    //
    // Module polygon
    //   top view
    //   v1                v2
    //    *---------------*
    //    |       ^ my    |
    //    |       |   mx  |
    //    |       *------>|
    //    |     center    |
    //    |               |
    //    *---------------*
    //   v0                v3
    //  (v4)
    //
    //   side view
    //    ----------------- top
    //    ----------------- bottom
    //
    vector<double> rv; // radius list from which we will find min/max.
    vector<double> zv; // z list from which we will find min/max.

    // mx: (v2+v3)/2, my: (v1+v2)/2
    XYZVector mx = 0.5*( module.basePoly().getVertex(2) + module.basePoly().getVertex(3) ) - center ;
    XYZVector my = 0.5*( module.basePoly().getVertex(1) + module.basePoly().getVertex(2) ) - center ;

    // new vertexes after expasion due to hybrid volumes
    const int npoints = 5; // v0,v1,v2,v3,v4(=v0)
    XYZVector v[npoints-1]; 
    v[0] = module.center() - (expandedModWidth/modWidth)*mx - (expandedModLength/modLength)*my; 
    v[1] = module.center() - (expandedModWidth/modWidth)*mx + (expandedModLength/modLength)*my; 
    v[2] = module.center() + (expandedModWidth/modWidth)*mx + (expandedModLength/modLength)*my; 
    v[3] = module.center() + (expandedModWidth/modWidth)*mx - (expandedModLength/modLength)*my; 

    // Fill all the vertex candidates (8 points)
    XYZVector v_top[npoints];    // module top surface
    XYZVector v_bottom[npoints]; // module bottom surface
    for (int ip = 0; ip < npoints-1; ip++) {
      v_top[ip]    = v[ip] + 0.5*modThickness*normal;
      v_bottom[ip] = v[ip] - 0.5*modThickness*normal;

        // for debuging
        vertex.push_back(v_top[ip]);
        vertex.push_back(v_bottom[ip]);

      zv.push_back(v_top[ip].Z());
      zv.push_back(v_bottom[ip].Z());
      v_top[ip].SetZ(0.);    // projection to xy plain.
      v_bottom[ip].SetZ(0.); // projection to xy plain.
      rv.push_back(v_top[ip].R());
      rv.push_back(v_bottom[ip].R());
    }
    // copy v0 as v4 for convenience
    v_top[npoints-1] = v_top[0];
    v_bottom[npoints-1] = v_bottom[0];

    // Fill possible side point candidates (8 sides)
    // No need to consider the 4 sides connecting top-bottom e.g. the side made from v0_top and v0_bottom.
    for (int ip = 0; ip < npoints-1; ip++) {
      // top surface
      XYZVector sp_top; // midpoint of each side
      sp_top = v_top[ip]-v_top[ip+1];
      if ( sp_top.Dot(v_top[ip]) * sp_top.Dot(v_top[ip+1]) < 0 ) { // minimum is a point divided internally 
        rv.push_back((v_top[ip]-sp_top.Unit().Dot(v_top[ip])*sp_top.Unit()).R());
      }
      // bottom surface
      XYZVector sp_bottom; // midpoint of each side
      sp_bottom = v_bottom[ip]-v_bottom[ip+1];
      if ( sp_bottom.Dot(v_bottom[ip]) * sp_bottom.Dot(v_bottom[ip+1]) < 0 ) { // minimum is a point divided internally 
        rv.push_back((v_bottom[ip]-sp_bottom.Unit().Dot(v_bottom[ip])*sp_bottom.Unit()).R());
      }
    }
    
    // Find min and max
    sort(rv.begin(),rv.end());
    sort(zv.begin(),zv.end());
    rmin = rv.at(0);
    rmax = rv.back();
    zmin = zv.at(0);
    zmax = zv.back();



    // Material assignment

    ElementsVector matElements = module.getLocalElements();
    ElementsVector::const_iterator meit;
    for (meit = matElements.begin(); meit != matElements.end(); meit++) {

       const MaterialObject::Element* el = *meit;

       // We skip in the case of ...
       // Filter sensor material and unexpected targetVolumes.
       if ( el->componentName() == "Sensor"     ||
            el->componentName() == "PS Sensors" ||
            el->componentName() == "2S Sensors"    ) {
          continue; // We will not handle sensors in this class 
       } else if ( el->targetVolume() == InnerSensor ||
                   el->targetVolume() == OuterSensor   ) { // Unexpected targetVolume ID 
         std::cerr << "!!!! ERROR !!!! : Found unexpected targetVolume." << std::endl;
         std::cerr << "targetVolume " << el->targetVolume() << " is only for sensors. Exit." << std::endl;
         std::exit(1);
       } else if ( el->targetVolume() >= nTypes   &&
                   el->targetVolume() != HybridFB &&     
                   el->targetVolume() != HybridLR &&     
                   el->targetVolume() != HybridFBLR_3456  ) {  
         std::cerr << "!!!! ERROR !!!! : Found unexpected targetVolume." << std::endl;
         std::cerr << "targetVolume " << el->targetVolume() << " is not supported. Exit." << std::endl;
         std::exit(1);
       }

       moduleMassWithoutSensors_expected += el->quantityInGrams(module);

       if ( el->targetVolume() == HybridFront   ||
            el->targetVolume() == HybridBack    ||
            el->targetVolume() == HybridLeft    ||
            el->targetVolume() == HybridRight   ||
            el->targetVolume() == HybridBetween ||
            el->targetVolume() == SupportPlate     ) {
          vol[el->targetVolume()]->addMaterial(el->elementName(),el->quantityInGrams(module));
          vol[el->targetVolume()]->addMass(el->quantityInGrams(module));
       } else if ( el->targetVolume() == HybridFB ) { 
          if (hybridFrontAndBackVolume_mm3 < 0) { // Need only once
            hybridFrontAndBackVolume_mm3 = vol[HybridFront]->getVolume()
                                         + vol[HybridBack]->getVolume();
          }
          vol[HybridFront]->addMaterial(el->elementName(),el->quantityInGrams(module));
          vol[HybridBack]->addMaterial(el->elementName(),el->quantityInGrams(module));
          vol[HybridFront]->addMass(el->quantityInGrams(module)*vol[HybridFront]->getVolume()/hybridFrontAndBackVolume_mm3);
          vol[HybridBack]->addMass(el->quantityInGrams(module)*vol[HybridBack]->getVolume()/hybridFrontAndBackVolume_mm3);
       } else if ( el->targetVolume() == HybridLR ) { 
          if (hybridLeftAndRightVolume_mm3 < 0) { // Need only once
            hybridLeftAndRightVolume_mm3 = vol[HybridLeft]->getVolume()
                                         + vol[HybridRight]->getVolume();
          }
          vol[HybridLeft]->addMaterial(el->elementName(),el->quantityInGrams(module));
          vol[HybridRight]->addMaterial(el->elementName(),el->quantityInGrams(module));
          vol[HybridLeft]->addMass(el->quantityInGrams(module)*vol[HybridLeft]->getVolume()/hybridLeftAndRightVolume_mm3);
          vol[HybridRight]->addMass(el->quantityInGrams(module)*vol[HybridRight]->getVolume()/hybridLeftAndRightVolume_mm3);
       } else if ( el->targetVolume() == HybridFBLR_0 || el->targetVolume() == HybridFBLR_3456 ) { // Uniformly Distribute
          vol[HybridFront]->addMaterial(el->elementName(),el->quantityInGrams(module));
          vol[HybridBack]->addMaterial(el->elementName(),el->quantityInGrams(module));
          vol[HybridLeft]->addMaterial(el->elementName(),el->quantityInGrams(module));
          vol[HybridRight]->addMaterial(el->elementName(),el->quantityInGrams(module));

          if (hybridTotalVolume_mm3 < 0) { // Need only once
            hybridTotalVolume_mm3 = vol[HybridFront]->getVolume()
                                  + vol[HybridBack]->getVolume()
                                  + vol[HybridLeft]->getVolume()
                                  + vol[HybridRight]->getVolume();
          }

          // Uniform density distribution and consistent with total mass
          vol[HybridFront]->addMass(el->quantityInGrams(module)*vol[HybridFront]->getVolume()/hybridTotalVolume_mm3); 
          vol[HybridBack]->addMass(el->quantityInGrams(module)*vol[HybridBack]->getVolume()/hybridTotalVolume_mm3);   
          vol[HybridLeft]->addMass(el->quantityInGrams(module)*vol[HybridLeft]->getVolume()/hybridTotalVolume_mm3);   
          vol[HybridRight]->addMass(el->quantityInGrams(module)*vol[HybridRight]->getVolume()/hybridTotalVolume_mm3);
       }
    }

    volumes.push_back(vol[HybridFront]);
    volumes.push_back(vol[HybridBack]);
    volumes.push_back(vol[HybridLeft]);
    volumes.push_back(vol[HybridRight]);
    volumes.push_back(vol[HybridBetween]);
    volumes.push_back(vol[SupportPlate]);

  }

  void ModuleComplex::addShapeInfo(std::vector<ShapeInfo>& vec) {
    ShapeInfo ele;
    ele.type = bx; // Box 
    std::vector<Volume*>::const_iterator vit;
    for ( vit = volumes.begin(); vit != volumes.end(); vit++ ) {
      if ( !((*vit)->getDensity()>0.) ) continue; 
      ele.name_tag = (*vit)->getName();
      ele.dx = (*vit)->getDx()/2.; // half length
      ele.dy = (*vit)->getDy()/2.; // half length
      ele.dz = (*vit)->getDz()/2.; // half length
      vec.push_back(ele);
    }
  }

  void ModuleComplex::addLogicInfo(std::vector<LogicalInfo>& vec) {
    LogicalInfo ele;
    std::vector<Volume*>::const_iterator vit;
    for ( vit = volumes.begin(); vit != volumes.end(); vit++ ) {
      if ( !((*vit)->getDensity()>0.) ) continue; 
      ele.name_tag     = (*vit)->getName();
      ele.shape_tag    = prefix_xmlfile + (*vit)->getName(); 
      ele.material_tag = prefix_xmlfile + prefix_material + (*vit)->getName();
      vec.push_back(ele);
    }
  }

  void ModuleComplex::addPositionInfo(std::vector<PosInfo>& vec) {
    PosInfo ele;
    ele.copy = 1;
    std::vector<Volume*>::const_iterator vit;
    for ( vit = volumes.begin(); vit != volumes.end(); vit++ ) {
      if ( !((*vit)->getDensity()>0.) ) continue; 
      ele.parent_tag = prefix_xmlfile + (*vit)->getParentName();
      ele.child_tag  = prefix_xmlfile + (*vit)->getName();
      ele.trans.dx   = (*vit)->getX();
      ele.trans.dy   = (*vit)->getY();
      ele.trans.dz   = (*vit)->getZ();
      vec.push_back(ele);
    }
  }

  void ModuleComplex::addMaterialInfo(std::vector<Composite>& vec) {
    std::vector<Volume*>::const_iterator vit;
    for ( vit = volumes.begin(); vit != volumes.end(); vit++ ) {
      if ( !((*vit)->getDensity()>0.) ) continue; 
      Composite comp;
      comp.name    = prefix_material + (*vit)->getName();
      comp.density = (*vit)->getDensity();
      comp.method  = wt;

      double m = 0.0;
      for (std::map<std::string, double>::const_iterator it = (*vit)->getMaterialList().begin(); 
                                                         it != (*vit)->getMaterialList().end(); ++it) {
          comp.elements.push_back(*it);
          m += it->second;
      }
   
      for (unsigned int i = 0; i < comp.elements.size(); i++)
        comp.elements.at(i).second = comp.elements.at(i).second / m;
      
      vec.push_back(comp);
    }
  }


  void ModuleComplex::print() const {
    std::cout << "ModuleComplex::print():" << std::endl;
    std::cout << "  Module Name:" << moduleId << std::endl;
    std::cout << "  Geometry Information:" << std::endl;
    std::cout << "    center postion : (" << center.X() << ","
                                          << center.Y() << "," 
                                          << center.Z() << ")" << std::endl;
    std::cout << "    normal vector  : (" << normal.X() << ","
                                          << normal.Y() << "," 
                                          << normal.Z() << ")" << std::endl;
    std::cout << "    module width     : " << expandedModWidth     << std::endl; 
    std::cout << "    module length    : " << expandedModLength    << std::endl; 
    std::cout << "    module thickness : " << expandedModThickness << std::endl; 
    std::cout << "    vertex points    : ";
    for (unsigned int i = 0; i < vertex.size(); i++) {
      if (i!=vertex.size()-1) 
        std::cout << "(" << vertex[i].X() << "," << vertex[i].Y() << "," << vertex[i].Z() << "), ";
      else 
        std::cout << "(" << vertex[i].X() << "," << vertex[i].Y() << "," << vertex[i].Z() << ")" << std::endl;
    }
    std::vector<Volume*>::const_iterator vit;
    double moduleTotalMass = 0.;
    for ( vit = volumes.begin(); vit != volumes.end(); vit++ ) {
      (*vit)->print();
      moduleTotalMass += (*vit)->getMass();
    }
    std::cerr << "  Module Total Mass = " << moduleTotalMass 
              << " (" << moduleMassWithoutSensors_expected << " is expected.)" << std::endl; 
  }
#endif
}

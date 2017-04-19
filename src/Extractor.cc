/**
 * @file Extractor.cc
 * @brief This class extracts the tracker geometry and material properties from an existing setup and groups them in preparation for translation to CMSSW XML
 */


//#define __DEBUGPRINT__
//#define __FLIPSENSORS_OUT__
//#define __FLIPSENSORS_IN__

#include <Extractor.hh>
#include <cstdlib>

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
  void Extractor::analyse(MaterialTable& mt, MaterialBudget& mb, XmlTags& trackerXmlTags, CMSSWBundle& d, bool wt) {

    std::cout << "Starting analysis..." << std::endl;

    Tracker& tr = mb.getTracker();
    InactiveSurfaces& is = mb.getInactiveSurfaces();
    bool isPixelTracker = tr.isPixelTracker();

    std::vector<std::vector<ModuleCap> >& bc = mb.getBarrelModuleCaps();
    std::vector<std::vector<ModuleCap> >& ec = mb.getEndcapModuleCaps();

    std::vector<Element>& e = d.elements;
    std::vector<Composite>& c = d.composites;
    std::vector<LogicalInfo>& l = d.logic;
    std::vector<ShapeInfo>& s = d.shapes;
    std::vector<ShapeOperationInfo>& so = d.shapeOps;
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
    so.clear(); // ShapeOperationInfo
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
    // This rotation places an unflipped module within a rod
    Rotation rot;
    rot.name = xml_places_unflipped_mod_in_rod;
    rot.thetax = 90.0;
    rot.phix = 90.0;
    rot.thetay = 0.0;
    rot.phiy = 0.0;
    rot.thetaz = 90.0;
    rot.phiz = 0.0;
    r.insert(std::pair<const std::string,Rotation>(rot.name,rot));

    // This rotation places a flipped module within a rod
    rot.name = xml_places_flipped_mod_in_rod;
    rot.thetax = 90.0;
    rot.phix = 270.0;
    rot.thetay = 0.0;
    rot.phiy = 0.0;
    rot.thetaz = 90.0;
    rot.phiz = 180.0;
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
    ShapeOperationInfo shapeOp;
    if (!wt) {
      shape.type = pc;
      shape.name_tag = xml_tob;

      // Barrel
      analyseBarrelContainer(tr, trackerXmlTags, shape.rzup, shape.rzdown);
      s.push_back(shape);
      std::cout << "Barrel container done." << std::endl;

      // Endcap
      analyseEndcapContainer(ec, tr, trackerXmlTags, shape.rzup, shape.rzdown);
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
    analyseLayers(mt, tr, trackerXmlTags, c, l, s, so, p, a, r, t, ri, wt);
    std::cout << "Barrel layers done." << std::endl;
    // Analyse endcaps
    analyseDiscs(mt, ec, tr, trackerXmlTags, c, l, s, so, p, a, r, t, ri, wt);
    std::cout << "Endcap discs done." << std::endl;
    // Analyse services
    analyseServices(is, isPixelTracker, trackerXmlTags, c, l, s, p, t);
    std::cout << "Barrel services done." << std::endl;
    // Analyse supports
    analyseSupports(is, isPixelTracker, trackerXmlTags, c, l, s, p, t);
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
      std::pair<double, int> AZ = getAZ(r.rlength, r.ilength);
      e.atomic_weight = AZ.first;
      e.atomic_number = AZ.second;
      // Z and A are calculated from radiation length and nuclear interaction lengths.
      // THIS IS BECAUSE RADIATION LENGTH AND NUCLEAR INTERACTION LENGTH CANNOT BE TRANSMITED DIRECTLY TO CMSSW.
      // Hence, Z and A (and density) only are transmitted to CMSSW.
      // On CMSSW side, radiation lengths and nuclear interaction lengths will be recomputed from this Z, A, and density info.
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
  void Extractor::analyseBarrelContainer(Tracker& t, XmlTags& trackerXmlTags, std::vector<std::pair<double, double> >& up,
                                         std::vector<std::pair<double, double> >& down) {
    
    bool isPixelTracker = t.isPixelTracker();

    std::pair<double, double> rz;
    double rmax = 0.0, zmax = 0.0, zmin = 0.0;
    up.clear();
    down.clear();
    std::vector<std::vector<ModuleCap> >::iterator oiter;
    std::vector<ModuleCap>::iterator iiter;

    LayerAggregator lagg;
    t.accept(lagg);
    auto bl = lagg.getBarrelLayers();
    
    int layer = 1;
    int n_of_layers = bl->size();

    lagg.postVisit();
    std::vector<std::vector<ModuleCap> >& bc = lagg.getBarrelCap();
    
    for (oiter = bc.begin(); oiter != bc.end(); oiter++) {
      double lrmin = std::numeric_limits<double>::max();
      double lrmax = 0;
      double lzmin = 0;
      double lzmax = 0; 
 
      for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	if (iiter->getModule().uniRef().side > 0 && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)){
	  int modRing = iiter->getModule().uniRef().ring;
	  // layer name
	  std::ostringstream lname;
	  lname << xml_layer << layer; // e.g. Layer1
	  // module name
	  std::ostringstream mname;
	  mname << xml_barrel_module << modRing << lname.str(); // .e.g. BModule1Layer1
	  // parent module name
	  std::string parentName = mname.str();
	  // build module volumes, with hybrids taken into account
	  ModuleComplex modcomplex(mname.str(),parentName,*iiter);
	  modcomplex.buildSubVolumes();
	  lrmin = MIN(lrmin, modcomplex.getRmin());
	  lrmax = MAX(lrmax, modcomplex.getRmax());	    
	  lzmax = MAX(lzmax, modcomplex.getZmax());
	  lzmin = -lzmax;
	}
      }

      if (layer == 1) {
	rz.first = lrmin;
	rz.second = lzmin;
	up.push_back(rz);
	rz.second = lzmax;
	down.push_back(rz);
      }
      else {
	//new barrel reached
	if (lzmax != zmax) {
	  //new layer sticks out compared to old layer
	  if (lzmax > zmax) rz.first = lrmin;
	  //old layer sticks out compared to new layer
	  else rz.first = rmax;
	  rz.second = zmin;
	  up.push_back(rz);
	  rz.second = zmax;
	  down.push_back(rz);
	  rz.second = lzmin;
	  up.push_back(rz);
	  rz.second = lzmax;
	  down.push_back(rz);
	}
      }
      //index size - 1
      if (layer == n_of_layers) {
	rz.first = lrmax;
	rz.second = lzmin;
	up.push_back(rz);
	rz.second = lzmax;
	down.push_back(rz);
      }
      rmax = lrmax;
      if (lzmin < 0) zmin = lzmin;
      if (lzmax > 0) zmax = lzmax;

      layer++;
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
  void Extractor::analyseEndcapContainer(std::vector<std::vector<ModuleCap> >& ec, Tracker& t, XmlTags& trackerXmlTags,
                                         std::vector<std::pair<double, double> >& up, std::vector<std::pair<double, double> >& down) {
    
    bool isPixelTracker = t.isPixelTracker();

    int first, last;
    std::pair<double, double> rz;
    double rmin = 0.0, rmax = 0.0, zmax = 0.0;
    up.clear();
    down.clear();
    std::vector<std::vector<ModuleCap> >::iterator oiter;
    std::vector<ModuleCap>::iterator iiter;  

    LayerAggregator lagg; 
    t.accept(lagg);

    auto el = lagg.getEndcapLayers();
    int layer = 1;
    int n_of_layers = el->size();
    bool hasfirst = false;

    //lagg.postVisit();   
    //std::vector<std::vector<ModuleCap> >& ec = lagg.getEndcapCap();
    
    for (oiter = ec.begin(); oiter != ec.end(); oiter++) {
      std::set<int> ridx;
      double lrmin = std::numeric_limits<double>::max();
      double lrmax = 0;
      double lzmin = std::numeric_limits<double>::max();
      double lzmax = 0; 
 
      for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	int modRing = iiter->getModule().uniRef().ring;
	if (ridx.find(modRing) == ridx.end()) {
	  ridx.insert(modRing);
	  // disk name
	  std::ostringstream dname;
	  dname << xml_disc << layer; // e.g. Disc6
	  // module name
	  std::ostringstream mname;
	  mname << xml_endcap_module << modRing << dname.str(); // e.g. EModule1Disc6
	  // parent module name  
	  std::string parentName = mname.str();
	  // build module volumes, with hybrids taken into account
	  ModuleComplex modcomplex(mname.str(),parentName,*iiter);
	  modcomplex.buildSubVolumes();
	  lrmin = MIN(lrmin, modcomplex.getRmin());
	  lrmax = MAX(lrmax, modcomplex.getRmax());	  
	  lzmin = MIN(lzmin, modcomplex.getZmin());  
	  lzmax = MAX(lzmax, modcomplex.getZmax());
	}
      }

      if ((lzmax > 0) && (!hasfirst)) {
	first = layer;
	hasfirst = true;
      }

      if (layer >= first) {
	if (layer == first) {
	  rmin = lrmin;
	  rmax = lrmax;
	  rz.first = rmax;
	  rz.second = lzmin - xml_z_pixfwd;
	  up.push_back(rz);
	  rz.first = rmin;
	  down.push_back(rz);
	}
	// disc beyond the first
	else {
	  // endcap change larger->smaller
	  if (rmax > lrmax) {
	    rz.second = zmax - xml_z_pixfwd;
	    rz.first = rmax;
	    up.push_back(rz);
	    rz.first = rmin;
	    down.push_back(rz);
	    rmax = lrmax;
	    rmin = lrmin;
	    rz.first = rmax;
	    up.push_back(rz);
	    rz.first = rmin;
	    down.push_back(rz);
	  }
	  // endcap change smaller->larger
	  if (rmax < lrmax) {
	    rz.second = lzmin - xml_z_pixfwd;
	    rz.first = rmax;
	    up.push_back(rz);
	    rz.first = rmin;
	    down.push_back(rz);
	    rmax = lrmax;
	    rmin = lrmin;
	    rz.first = rmax;
	    up.push_back(rz);
	    rz.first = rmin;
	    down.push_back(rz);
	  }
	}
	zmax = lzmax;
	// special treatment for last disc
	if (layer == n_of_layers) {
	  rz.first = rmax;
	  rz.second = zmax - xml_z_pixfwd;
	  up.push_back(rz);
	  rz.first = rmin;
	  down.push_back(rz);
	}
      }
      layer++;
    }
  }

  /**
   * This is one of the two main analysis functions that provide the core functionality of this class. 
   * It examines the barrel layers and the modules within, extracting a great range of different pieces of information 
   * from the geometry layout. The volumes considered are layers, rods, potentially tilted rings, modules (with wafer, 
   * active surfaces, hybrids, support plate). They form hierarchies of volumes, one inside the other. 
   * Output information are volume hierarchy, material, shapes, positioning (potential use of algorithm and rotations). 
   * They is also some topology, such as which volumes contain the active surfaces, and how those active surfaces are subdivided
   * and connected to the readout electronics. Last but not least, overall radiation and interaction lengths for each layer are 
   * calculated and stored; those are used as approximative values for certain CMSSW functions later on.
   * @param mt A reference to the global material table; used as input
   * @param bc A reference to the collection of material properties of the barrel modules; used as input
   * @param tr A reference to the tracker object; used as input
   * @param c A reference to the collection of composite material information; used for output
   * @param l A reference to the collection of volume hierarchy information; used for output
   * @param s A reference to the collection of shape parameters; used for output
   * @param so A reference to the collection of operations on physical volumes
   * @param p A reference to the collection of volume positionings; used for output
   * @param a A reference to the collection of algorithm calls and their parameters; used for output
   * @param r A reference to the collection of rotations; used for output
   * @param t A reference to the collection of topology information; used for output
   * @param ri A reference to the collection of overall radiation and interaction lengths per layer or disc; used for output
   */
  void Extractor::analyseLayers(MaterialTable& mt/*, std::vector<std::vector<ModuleCap> >& bc*/, Tracker& tr, XmlTags& trackerXmlTags,
                                std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<ShapeOperationInfo>& so, 
				std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::map<std::string,Rotation>& r, std::vector<SpecParInfo>& t, 
				std::vector<RILengthInfo>& ri, bool wt) {

    bool isPixelTracker = tr.isPixelTracker();

    // Container inits
    ShapeInfo shape;
    shape.dyy = 0.0;

    ShapeOperationInfo shapeOp;

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
    SpecParInfo rocdims, lspec, rspec, srspec, trspec, sspec, mspec;    
    // Layer
    lspec.name = trackerXmlTags.topo_layer_name + xml_par_tail;
    lspec.parameter.first = xml_tkddd_structure;
    lspec.parameter.second = trackerXmlTags.topo_layer_value;
    // Rod (straight or tilted)
    rspec.name = xml_subdet_straight_or_tilted_rod + xml_par_tail;
    rspec.parameter.first = xml_tkddd_structure;
    rspec.parameter.second = xml_det_straight_or_tilted_rod;
    // Straight Rod
    srspec.name = trackerXmlTags.topo_straight_rod_name + xml_par_tail;
    srspec.parameter.first = xml_tkddd_structure;
    srspec.parameter.second = trackerXmlTags.topo_straight_rod_value;
    // Tilted Ring (if any)
    trspec.name = trackerXmlTags.topo_tilted_ring_name + xml_par_tail;
    trspec.parameter.first = xml_tkddd_structure;
    trspec.parameter.second = trackerXmlTags.topo_tilted_ring_value;
    // Module stack
    sspec.name = trackerXmlTags.topo_bmodule_name + xml_par_tail;
    sspec.parameter.first = xml_tkddd_structure;
    sspec.parameter.second =  trackerXmlTags.topo_bmodule_value;
    // Module detectors
    mspec.name = xml_subdet_tobdet + xml_par_tail;
    mspec.parameter.first = xml_tkddd_structure;
    mspec.parameter.second = xml_det_tobdet;


    // material properties
    RILengthInfo ril;
    ril.barrel = true;
    ril.index = 0;


    // aggregate information about the modules
    LayerAggregator lagg;
    tr.accept(lagg);
    lagg.postVisit();
    std::vector<std::vector<ModuleCap> >& bc = lagg.getBarrelCap();

    std::vector<std::vector<ModuleCap> >::iterator oiter;
    std::vector<ModuleCap>::iterator iiter;
    
    int layer = 1;
    std::map<int, double> rodThickness;
    std::map<int, double> rodWidth;

    // LOOP ON LAYERS
    for (oiter = bc.begin(); oiter != bc.end(); oiter++) {
      
      // is the layer tilted ?
      bool isTilted = lagg.getBarrelLayers()->at(layer - 1)->isTilted();
     
      // TO DO : NOT THAT EXTREMA OF MODULES WITH HYBRIDS HAVE BEEN INTRODUCED INTO DETECTORMODULE (USED TO BE IN MODULE COMPLEX CLASS ONLY)
      // ALL THIS SHOULD BE DELETED AND REPLACED BY module->minRwithHybrids() GETTERS.

      // Calculate geometrical extrema of rod (straight layer), or of rod part + tilted ring (tilted layer)
      // straight layer : x and y extrema of rod
      double xmin = std::numeric_limits<double>::max();
      double xmax = 0;
      double ymin = std::numeric_limits<double>::max();
      double ymax = 0;
      // straight or tilted layer : z and r extrema of rod (straight layer), or of {rod part + tilted ring} (tilted layer)
      //double zmin = std::numeric_limits<double>::max();
      double zmax = 0;
      double rmin = std::numeric_limits<double>::max();
      double rmax = 0;
      // tilted layer : x, y, z, and r extrema of rod part
      double flatPartMinX = std::numeric_limits<double>::max();
      double flatPartMaxX = 0;
      double flatPartMinY = std::numeric_limits<double>::max();
      double flatPartMaxY = 0;
      double flatPartMaxZ = 0;
      double flatPartMinR = std::numeric_limits<double>::max();
      double flatPartMaxR = 0;
      // tilted layer : z extrema of one before last module of flat part
      int flatPartNumModules = lagg.getBarrelLayers()->at(layer - 1)->buildNumModulesFlat();
      double flatPartOneBeforeLastModuleMaxZ = 0;
      // straight or tilted layer : radii of rods (straight layer) or of rod parts (tilted layer)
      double RadiusIn = 0;
      double RadiusOut = 0;
      // loop on module caps
      for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	// only positive side, and modules with uniref phi == 1 or 2
	if (iiter->getModule().uniRef().side > 0 && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)) {
	  int modRing = iiter->getModule().uniRef().ring;
	  // layer name
	  std::ostringstream lname;
	  lname << xml_layer << layer; // e.g. Layer1
	  // module name
	  std::ostringstream mname;
	  mname << xml_barrel_module << modRing << lname.str(); //.e.g. BModule1Layer1
	  // parent module name
	  std::string parentName = mname.str();
	  // build module volumes, with hybrids taken into account
	  ModuleComplex modcomplex(mname.str(),parentName,*iiter);
	  modcomplex.buildSubVolumes();
	  if (iiter->getModule().uniRef().phi == 1) {
	    rodThickness.insert( {layer,  modcomplex.getExpandedModuleThickness() / 2.} );
	    rodWidth.insert( {layer, modcomplex.getExpandedModuleWidth() / 2.} );
	    xmin = MIN(xmin, modcomplex.getXmin());
	    xmax = MAX(xmax, modcomplex.getXmax());
	    ymin = MIN(ymin, modcomplex.getYmin());
	    ymax = MAX(ymax, modcomplex.getYmax());
	    // tilted layer : rod part
	    if (isTilted && iiter->getModule().tiltAngle() == 0) {
	      flatPartMinX = MIN(flatPartMinX, modcomplex.getXmin());
	      flatPartMaxX = MAX(flatPartMaxX, modcomplex.getXmax());
	      flatPartMinY = MIN(flatPartMinY, modcomplex.getYmin());
	      flatPartMaxY = MAX(flatPartMaxY, modcomplex.getYmax());
	    }
	  }
	  // for z and r, uniref phi == 2 has to be taken into account too
	  // (because different from uniref phi == 1 in case of tilted layer)
	  zmax = MAX(zmax, modcomplex.getZmax());
	  //zmin = -zmax;
	  rmin = MIN(rmin, modcomplex.getRmin());
	  rmax = MAX(rmax, modcomplex.getRmax());
	  // tilted layer : rod part
	  if (isTilted && (iiter->getModule().tiltAngle() == 0)) {
	    flatPartMaxZ = MAX(flatPartMaxZ, modcomplex.getZmax());
	    flatPartMinR = MIN(flatPartMinR, modcomplex.getRmin());
	    flatPartMaxR = MAX(flatPartMaxR, modcomplex.getRmax());
	    if (flatPartNumModules >= 2 && modRing == (flatPartNumModules - 1)) { flatPartOneBeforeLastModuleMaxZ = MAX(flatPartOneBeforeLastModuleMaxZ, modcomplex.getZmax()); }
	  }
	  // both modRings 1 and 2 have to be taken into account because of small delta
	  if (iiter->getModule().uniRef().phi == 1 && (modRing == 1 || modRing == 2)) { RadiusIn = RadiusIn + iiter->getModule().center().Rho() / 2; }
	  if (iiter->getModule().uniRef().phi == 2 && (modRing == 1 || modRing == 2)) { RadiusOut = RadiusOut + iiter->getModule().center().Rho() / 2; }
	  //std::cout << 
	}
      }

      if ((rmax - rmin) == 0.0) continue;


      shape.type = bx; // box
      shape.rmin = 0.0;
      shape.rmax = 0.0;


      // for material properties
      double rtotal = 0.0, itotal = 0.0;     
      int count = 0;
      ril.index = layer;
      

      std::ostringstream lname, rodname, rodNextPhiName, pconverter;
      lname << xml_layer << layer; // e.g. Layer1
      if (!isPixelTracker) rodname << xml_rod << layer; // e.g.Rod1
      else {
	rodname << xml_rod << xml_flipped << layer; // e.g.RodFlipped1
	rodNextPhiName << xml_rod << xml_unflipped << layer; // e.g.RodUnflipped1
      }

      double rodStartPhiAngle, rodNextPhiStartPhiAngle;

      // information on tilted rings, indexed by ring number
      std::map<int, BTiltedRingInfo> rinfoplus; // positive-z side
      std::map<int, BTiltedRingInfo> rinfominus; // negative-z side


      // LOOP ON MODULE CAPS 
      for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
        
	// ONLY POSITIVE SIDE, AND MODULES WITH UNIREF PHI == 1 OR 2
	if (iiter->getModule().uniRef().side > 0 && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)) {

	  // ring number (position on rod, or tilted ring number)
	  int modRing = iiter->getModule().uniRef().ring; 
    
	  // tilt angle of the module
	  double tiltAngle = 0;
	  if (isTilted) { tiltAngle = iiter->getModule().tiltAngle() * 180. / M_PI; }
	    
	  //std::cout << "iiter->getModule().uniRef().phi = " << iiter->getModule().uniRef().phi << " iiter->getModule().posRef().phi = " << iiter->getModule().posRef().phi << "iiter->getModule().center().Phi() = " << iiter->getModule().center().Phi() * 180. / M_PI << " iiter->getModule().center().Rho() = " << iiter->getModule().center().Rho() << " iiter->getModule().center().X() = " << iiter->getModule().center().X() << " iiter->getModule().center().Y() = " << iiter->getModule().center().Y() << " iiter->getModule().center().Z() = " << iiter->getModule().center().Z() << " tiltAngle = " << tiltAngle << " iiter->getModule().flipped() = " << iiter->getModule().flipped() << " iiter->getModule().moduleType() = " << iiter->getModule().moduleType() << std::endl;

	  // module name
	  std::ostringstream mname;
	  mname << xml_barrel_module << modRing << lname.str(); // e.g. BModule1Layer1

	  // parent module name
	  std::string parentName = mname.str();

	  std::vector<ModuleCap>::iterator partner;

	  // build module volumes, with hybrids taken into account
	  ModuleComplex modcomplex(mname.str(),parentName,*iiter);
	  modcomplex.buildSubVolumes();
#ifdef __DEBUGPRINT__
	  modcomplex.print();
#endif

	  // ROD 1 (STRAIGHT LAYER), OR ROD 1 + MODULES WITH UNIREF PHI == 1 OF THE TILTED RINGS (TILTED LAYER)
	  if (iiter->getModule().uniRef().phi == 1) {           

	    rodStartPhiAngle = iiter->getModule().center().Phi();

            std::ostringstream ringname;
	    ringname << xml_ring << modRing << lname.str();


	    // MODULE

	    // For SolidSection in tracker.xml : module's box shape
            shape.name_tag = mname.str();
            //shape.dx = iiter->getModule().area() / iiter->getModule().length() / 2.0;
            //shape.dy = iiter->getModule().length() / 2.0;
            //shape.dz = iiter->getModule().thickness() / 2.0;
            shape.dx = modcomplex.getExpandedModuleWidth()/2.0;
            shape.dy = modcomplex.getExpandedModuleLength()/2.0;
            shape.dz = modcomplex.getExpandedModuleThickness()/2.0;
	    s.push_back(shape);
	    

	    // For LogicalPartSection in tracker.xml : module's material
            //logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
            logic.material_tag = xml_material_air;
            logic.name_tag = mname.str();
            logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
	    l.push_back(logic);	    
	    // module composite material
            //matname << xml_base_actcomp << "L" << layer << "P" << modRing;
            //c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));
            

	    // For PosPart section in tracker.xml : module's positions in rod (straight layer) or rod part (tilted layer)
            if (!isTilted || (isTilted && (tiltAngle == 0))) {
	      pos.parent_tag = trackerXmlTags.nspace + ":" + rodname.str();
	      pos.child_tag = trackerXmlTags.nspace + ":" + mname.str();
	      partner = findPartnerModule(iiter, oiter->end(), modRing);

	      pos.trans.dx = iiter->getModule().center().Rho() - RadiusIn;
	      pos.trans.dz = iiter->getModule().center().Z();
	      if (!iiter->getModule().flipped()) { pos.rotref = trackerXmlTags.nspace + ":" + xml_places_unflipped_mod_in_rod; }
	      else { pos.rotref = trackerXmlTags.nspace + ":" + xml_places_flipped_mod_in_rod; }
	      p.push_back(pos);
	      
	      // This is a copy of the BModule (FW/BW barrel half)
	      if (partner != oiter->end()) {
		pos.trans.dx = partner->getModule().center().Rho() - RadiusIn;
		pos.trans.dz = partner->getModule().center().Z();
		if (!partner->getModule().flipped()) { pos.rotref = trackerXmlTags.nspace + ":" + xml_places_unflipped_mod_in_rod; }
		else { pos.rotref = trackerXmlTags.nspace + ":" + xml_places_flipped_mod_in_rod; }
		pos.copy = 2; 
		p.push_back(pos);
		pos.copy = 1;
	      }
	      pos.rotref = "";
	    }
	  }

	  if (isPixelTracker && iiter->getModule().uniRef().phi == 2) {
	    if (!isTilted || (isTilted && (tiltAngle == 0))) {
	      rodNextPhiStartPhiAngle = iiter->getModule().center().Phi();
	      pos.parent_tag = trackerXmlTags.nspace + ":" + rodNextPhiName.str();
	      pos.child_tag = trackerXmlTags.nspace + ":" + mname.str();
	      partner = findPartnerModule(iiter, oiter->end(), modRing);

	      pos.trans.dx = iiter->getModule().center().Rho() - RadiusOut;
	      pos.trans.dz = iiter->getModule().center().Z();
	      if (!iiter->getModule().flipped()) { pos.rotref = trackerXmlTags.nspace + ":" + xml_places_unflipped_mod_in_rod; }
	      else { pos.rotref = trackerXmlTags.nspace + ":" + xml_places_flipped_mod_in_rod; }
	      p.push_back(pos);
	      
	      // This is a copy of the BModule (FW/BW barrel half)
	      if (partner != oiter->end()) {
		pos.trans.dx = partner->getModule().center().Rho() - RadiusOut;
		pos.trans.dz = partner->getModule().center().Z();
		if (!partner->getModule().flipped()) { pos.rotref = trackerXmlTags.nspace + ":" + xml_places_unflipped_mod_in_rod; }
		else { pos.rotref = trackerXmlTags.nspace + ":" + xml_places_flipped_mod_in_rod; }
		pos.copy = 2; 
		p.push_back(pos);
		pos.copy = 1;
	      }
	      pos.rotref = "";
	    }
	  }

	  if (iiter->getModule().uniRef().phi == 1) {

	    std::ostringstream ringname;
	    ringname << xml_ring << modRing << lname.str();

	    // Topology
	    sspec.partselectors.push_back(mname.str());
	    sspec.moduletypes.push_back(minfo_zero);


            // WAFER
            string xml_base_lowerupper = "";
            if (iiter->getModule().numSensors() == 2) {
	      xml_base_lowerupper = xml_base_lower;
	      shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_waf;
	    }	    
	    else {
	      if (iiter->getModule().isPixelModule()) shape.name_tag = mname.str() + xml_PX + xml_base_waf;
	      else if (iiter->getModule().isTimingModule()) shape.name_tag = mname.str() + xml_timing + xml_base_waf;
	      else { std::cerr << "Wafer : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	    }

	    // SolidSection
            shape.name_tag = shape.name_tag;
	    shape.dx = iiter->getModule().area() / iiter->getModule().length() / 2.0;
	    shape.dy = iiter->getModule().length() / 2.0;
            shape.dz = iiter->getModule().sensorThickness() / 2.0;
            s.push_back(shape);   

	    // LogicalPartSection
            logic.name_tag = shape.name_tag;
            logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
            logic.material_tag = xml_material_air;
            l.push_back(logic);

	    // PosPart section
            pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str();
            pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
            pos.trans.dx = 0.0;
            pos.trans.dz = /*shape.dz*/ - iiter->getModule().dsDistance() / 2.0; 
            p.push_back(pos);

            if (iiter->getModule().numSensors() == 2) {

              xml_base_lowerupper = xml_base_upper;

	      // SolidSection
              shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_waf;
              s.push_back(shape);

	      // LogicalPartSection
              logic.name_tag = shape.name_tag;
              logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
              l.push_back(logic);

	      // PosPart section
              pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
              pos.trans.dz = pos.trans.dz + /*2 * shape.dz +*/ iiter->getModule().dsDistance();  // CUIDADO: was with 2*shape.dz, but why???
              //pos.copy = 2;

              if (iiter->getModule().stereoRotation() != 0) {
                rot.name = type_stereo + mname.str();
                rot.thetax = 90.0;
                rot.phix = iiter->getModule().stereoRotation() / M_PI * 180.;
                rot.thetay = 90.0;
                rot.phiy = 90.0 + iiter->getModule().stereoRotation() / M_PI * 180.;
                r.insert(std::pair<const std::string,Rotation>(rot.name,rot));
                pos.rotref = trackerXmlTags.nspace + ":" + rot.name;
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


            // ACTIVE SURFACE
            xml_base_lowerupper = "";
            if (iiter->getModule().numSensors() == 2) xml_base_lowerupper = xml_base_lower;

	    if (iiter->getModule().moduleType() == "ptPS") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_ps + xml_base_pixel + xml_base_act;
	    else if (iiter->getModule().moduleType() == "pt2S") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_2s+ xml_base_act;
	    else if (iiter->getModule().isTimingModule()) shape.name_tag = mname.str() + xml_timing + xml_base_act;
	    else if (iiter->getModule().isPixelModule()) shape.name_tag = mname.str() + xml_PX + xml_base_Act;
	    else { std::cerr << "Active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }

	    // SolidSection
	    shape.dx = iiter->getModule().area() / iiter->getModule().length() / 2.0;
	    shape.dy = iiter->getModule().length() / 2.0;
            shape.dz = iiter->getModule().sensorThickness() / 2.0;
            s.push_back(shape);   

	    // LogicalPartSection
            logic.name_tag = shape.name_tag;
            logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
            logic.material_tag = xml_fileident + ":" + xml_tkLayout_material + xml_sensor_silicon;
            l.push_back(logic);

	    // PosPart section
	    if (iiter->getModule().numSensors() == 2) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_base_lowerupper + xml_base_waf;
	    else {
	      if (iiter->getModule().isTimingModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_timing + xml_base_waf;
	      else if (iiter->getModule().isPixelModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_PX + xml_base_waf;
	      else { std::cerr << "Positioning active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	    }
            pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
            pos.trans.dz = 0.0;
#ifdef __FLIPSENSORS_IN__ // Flip INNER sensors
            pos.rotref = trackerXmlTags.nspace + ":" + rot_sensor_tag;
#endif
            p.push_back(pos);

            // Topology
	    minfo.name		= iiter->getModule().moduleType();
            mspec.partselectors.push_back(shape.name_tag);
            minfo.rocrows	= any2str<int>(iiter->getModule().innerSensor().numROCRows());  // in case of single sensor module innerSensor() and outerSensor() point to the same sensor
            minfo.roccols	= any2str<int>(iiter->getModule().innerSensor().numROCCols());
            minfo.rocx		= any2str<int>(iiter->getModule().innerSensor().numROCX());
            minfo.rocy		= any2str<int>(iiter->getModule().innerSensor().numROCY());
            mspec.moduletypes.push_back(minfo);

            if (iiter->getModule().numSensors() == 2) { 

              xml_base_lowerupper = xml_base_upper;

	      // SolidSection
	      if (iiter->getModule().moduleType() == "ptPS") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_ps + xml_base_strip + xml_base_act;
	      else if (iiter->getModule().moduleType() == "pt2S") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_2s+ xml_base_act;
	      else { std::cerr << "Active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
              s.push_back(shape);

	      // LogicalPartSection
              logic.name_tag = shape.name_tag;
              logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
              l.push_back(logic);

	      // PosPart section
	      pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_base_lowerupper + xml_base_waf;
              pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
#ifdef __FLIPSENSORS_OUT__ // Flip OUTER sensors
              pos.rotref = trackerXmlTags.nspace + ":" + rot_sensor_tag;
#endif
              p.push_back(pos);

              // Topology
              mspec.partselectors.push_back(shape.name_tag);
              minfo.rocrows	= any2str<int>(iiter->getModule().outerSensor().numROCRows());
              minfo.roccols	= any2str<int>(iiter->getModule().outerSensor().numROCCols());
              minfo.rocx	= any2str<int>(iiter->getModule().outerSensor().numROCX());
              minfo.rocy	= any2str<int>(iiter->getModule().outerSensor().numROCY());
              mspec.moduletypes.push_back(minfo);
	    } // End of replica for Pt-modules

	    modcomplex.addMaterialInfo(c);
	    modcomplex.addShapeInfo(s);
	    modcomplex.addLogicInfo(l);
	    modcomplex.addPositionInfo(p);
#ifdef __DEBUGPRINT__
	    modcomplex.print();
#endif
	    

	    // collect tilted ring info
	    if (isTilted && (tiltAngle != 0)) {
	      BTiltedRingInfo rinf;
	      // ring on positive-z side
	      rinf.name = ringname.str() + xml_plus;	      
	      rinf.childname = mname.str();
	      rinf.isZPlus = 1;
	      rinf.tiltAngle = tiltAngle;
	      rinf.bw_flipped = iiter->getModule().flipped();
	      rinf.modules = lagg.getBarrelLayers()->at(layer - 1)->numRods();
	      rinf.r1 = iiter->getModule().center().Rho();
	      rinf.z1 = iiter->getModule().center().Z();
	      rinf.startPhiAngle1 = iiter->getModule().center().Phi();     
	      rinf.rmin = modcomplex.getRmin();
	      rinf.zmin = modcomplex.getZmin();
	      rinf.rminatzmin = modcomplex.getRminatZmin();      
	      rinfoplus.insert(std::pair<int, BTiltedRingInfo>(modRing, rinf));

	      // same ring on negative-z side
	      rinf.name = ringname.str() + xml_minus;
	      rinf.isZPlus = 0;
	      rinf.z1 = - iiter->getModule().center().Z();
	      rinfominus.insert(std::pair<int, BTiltedRingInfo>(modRing, rinf));
	    }


            // material properties
            rtotal = rtotal + iiter->getRadiationLength();
            itotal = itotal + iiter->getInteractionLength();
            count++;
	    //double dt = modcomplex.getExpandedModuleThickness();
	  }


	  // ONLY MODULES WITH UNIREF PHI == 2 OF THE TILTED RINGS (TILTED LAYER)
	  if (isTilted && (iiter->getModule().uniRef().phi == 2)) {
	    std::map<int,BTiltedRingInfo>::iterator it;
	    // fill the info of the z-positive ring with matching ring number
	    it = rinfoplus.find(modRing);
	    if (it != rinfoplus.end()) {
	      it->second.fw_flipped = iiter->getModule().flipped();
	      it->second.r2 = iiter->getModule().center().Rho();
	      it->second.z2 = iiter->getModule().center().Z();
	      it->second.startPhiAngle2 = iiter->getModule().center().Phi();
	      it->second.rmax = modcomplex.getRmax();
	      it->second.zmax = modcomplex.getZmax();
	      it->second.rmaxatzmax = modcomplex.getRmaxatZmax();
	    }
	    // fill the info of the z-negative ring with matching ring number
	    it = rinfominus.find(modRing);
	    if (it != rinfominus.end()) {
	      it->second.fw_flipped = iiter->getModule().flipped();
	      it->second.r2 = iiter->getModule().center().Rho();
	      it->second.z2 = - iiter->getModule().center().Z();
	      it->second.startPhiAngle2 = iiter->getModule().center().Phi();
	      it->second.rmax = modcomplex.getRmax();
	      it->second.zmax = modcomplex.getZmax();
	      it->second.rmaxatzmax = modcomplex.getRmaxatZmax();
	    }
	  }
	}
      }
      // material properties
      if (count > 0) {
	ril.rlength = rtotal / (double)count;
	ril.ilength = itotal / (double)count;
	ri.push_back(ril);
      }


      // rod(s)
      shape.name_tag = rodname.str();
      shape.dx = (ymax - ymin) / 2 + xml_epsilon;
      if (isTilted) shape.name_tag = rodname.str() + "Full";
      if (isTilted) shape.dx = (flatPartMaxY - flatPartMinY) / 2 + xml_epsilon;
      if (isPixelTracker) shape.dx = rodThickness.at(layer) + xml_epsilon;   // TO DO : IMPLEMENT SMALL DELTA FOR THE TILTED !!!! OTHERWISE, THIS VALUE WILL BE FALSE FOR A TILTED PIXEL TRACKER EXPORT !!!
      shape.dy = (xmax - xmin) / 2 + xml_epsilon;
      if (isTilted) shape.dy = (flatPartMaxX - flatPartMinX) / 2 + xml_epsilon;
      if (isPixelTracker) shape.dy = rodWidth.at(layer) + xml_epsilon;
      shape.dz = zmax + xml_epsilon;
      if (isTilted) shape.dz = flatPartMaxZ + xml_epsilon;
      s.push_back(shape);

      // Subtraction of an air volume from the flat part rod container volume, to avoid collision with first tilted ring
      if (isTilted && flatPartNumModules >= 2) {
	shape.name_tag = rodname.str() + "Air";
	shape.dx = shape.dx / 2.0;
	shape.dy = shape.dy + xml_epsilon;
	shape.dz = (shape.dz - flatPartOneBeforeLastModuleMaxZ) / 2.;
	s.push_back(shape);
	
	shapeOp.name_tag = rodname.str() + "SubtractionIntermediate";
	shapeOp.type = substract;
	shapeOp.rSolid1 = rodname.str() + "Full";
	shapeOp.rSolid2 = rodname.str() + "Air";
	shapeOp.trans.dx = shape.dx + xml_epsilon;
	shapeOp.trans.dy = 0.;
	shapeOp.trans.dz = flatPartMaxZ + xml_epsilon - shape.dz + xml_epsilon;
	so.push_back(shapeOp);

	shapeOp.name_tag = rodname.str();
	shapeOp.type = substract;
	shapeOp.rSolid1 = rodname.str() + "SubtractionIntermediate";
	shapeOp.rSolid2 = rodname.str() + "Air";
	shapeOp.trans.dx = shape.dx + xml_epsilon;
	shapeOp.trans.dy = 0.;
	shapeOp.trans.dz = -shapeOp.trans.dz;
	so.push_back(shapeOp);
      }

      logic.name_tag = rodname.str();
      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
      logic.material_tag = xml_material_air;
      l.push_back(logic);

      if (isPixelTracker) {
	shape.name_tag = rodNextPhiName.str();
	s.push_back(shape);
	logic.name_tag = rodNextPhiName.str();
	logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
	l.push_back(logic);
	rspec.partselectors.push_back(rodNextPhiName.str());
	srspec.partselectors.push_back(rodNextPhiName.str());
      }

      rspec.partselectors.push_back(rodname.str());
      rspec.moduletypes.push_back(minfo_zero);
      srspec.partselectors.push_back(rodname.str());
      srspec.moduletypes.push_back(minfo_zero);
        

      
      // rods in layer algorithm(s)
      if (!isPixelTracker) {
	alg.name = xml_phialt_algo;
	alg.parent = trackerXmlTags.nspace + ":" + lname.str();
	pconverter <<  trackerXmlTags.nspace + ":" + rodname.str();
	alg.parameters.push_back(stringParam(xml_childparam, pconverter.str()));
	alg.parameters.push_back(numericParam(xml_tilt, "90*deg")); // This "tilt" here has nothing to do with the tilt angle of a tilted TBPS.
	// It is an angle used internally by PhiAltAlgo to shift in Phi the startAngle ( in (X,Y) plane).
	// 90 deg corresponds to no shift. SHOULD NOT BE MODIFIED!!
	pconverter.str("");
	pconverter << rodStartPhiAngle * 180. / M_PI << "*deg";
	alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
	pconverter.str("");
	alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
	pconverter << RadiusIn << "*mm";
	alg.parameters.push_back(numericParam(xml_radiusin, pconverter.str()));
	pconverter.str("");
	pconverter << RadiusOut << "*mm";
	alg.parameters.push_back(numericParam(xml_radiusout, pconverter.str()));
	pconverter.str("");
	alg.parameters.push_back(numericParam(xml_zposition, "0.0*mm"));
	pconverter << lagg.getBarrelLayers()->at(layer - 1)->numRods();
	alg.parameters.push_back(numericParam(xml_number, pconverter.str()));
	pconverter.str("");
	alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
	alg.parameters.push_back(numericParam(xml_incrcopyno, "1"));
	a.push_back(alg);
	alg.parameters.clear();
      }
      else {
	alg.name = xml_angular_algo;
	alg.parent = trackerXmlTags.nspace + ":" + lname.str();
	pconverter <<  trackerXmlTags.nspace + ":" + rodname.str();
	alg.parameters.push_back(stringParam(xml_childparam, pconverter.str()));
	pconverter.str("");
	pconverter << rodStartPhiAngle * 180. / M_PI << "*deg";
	alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
	pconverter.str("");
	alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
	pconverter << RadiusIn << "*mm";
	alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
	pconverter.str("");
	alg.parameters.push_back(vectorParam(0., 0., 0.));
	pconverter << lagg.getBarrelLayers()->at(layer - 1)->numRods() / 2;
	alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
	pconverter.str("");
	alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
	alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
	a.push_back(alg);
	alg.parameters.clear();

	alg.name = xml_angular_algo;
	alg.parent = trackerXmlTags.nspace + ":" + lname.str();
	pconverter <<  trackerXmlTags.nspace + ":" + rodNextPhiName.str();
	alg.parameters.push_back(stringParam(xml_childparam, pconverter.str()));
	pconverter.str("");
	pconverter << rodNextPhiStartPhiAngle * 180. / M_PI << "*deg";
	alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
	pconverter.str("");
	alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
	pconverter << RadiusOut << "*mm";
	alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
	pconverter.str("");
	alg.parameters.push_back(vectorParam(0., 0., 0.));
	pconverter << lagg.getBarrelLayers()->at(layer - 1)->numRods() / 2;
	alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
	pconverter.str("");
	alg.parameters.push_back(numericParam(xml_startcopyno, "2"));
	alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
	a.push_back(alg);
	alg.parameters.clear();
      }

      // reset
      shape.dx = 0.0;
      shape.dy = 0.0;
      shape.dyy = 0.0;	  
      pos.trans.dx = 0;
      pos.trans.dy = 0;
      pos.trans.dz = 0;

      // tilted rings
      if ( !rinfoplus.empty() || !rinfominus.empty() ) {

	std::map<std::string, std::map<int,BTiltedRingInfo>> rinfototal;
	if ( !rinfoplus.empty() ) { rinfototal.insert({"rinfoplus", rinfoplus}); }
	if ( !rinfominus.empty() ) { rinfototal.insert({"rinfominus", rinfominus}); }

	for (auto const &rinfoside : rinfototal) {

	  for (auto const &ringinfo : rinfoside.second) {
	    auto const& rinfo = ringinfo.second;
	    if (rinfo.modules > 0) {

	      // reset
	      shape.rmin = 0.0;
	      shape.rmax = 0.0;

	      // section of cone
	      shape.name_tag = rinfo.name + "Cone";
	      shape.type = co;
	      shape.dz = (rinfo.zmax - rinfo.zmin) / 2. + xml_epsilon;
	      if (rinfo.isZPlus) {
		shape.rmin1 = rinfo.rminatzmin - xml_epsilon * tan(rinfo.tiltAngle * M_PI / 180.0);
		shape.rmax1 = rinfo.rmaxatzmax + 2 * shape.dz * tan(rinfo.tiltAngle * M_PI / 180.0) + xml_epsilon * tan(rinfo.tiltAngle * M_PI / 180.0);
		shape.rmin2 = rinfo.rminatzmin - 2 * shape.dz * tan(rinfo.tiltAngle * M_PI / 180.0) - xml_epsilon * tan(rinfo.tiltAngle * M_PI / 180.0);
		shape.rmax2 = rinfo.rmaxatzmax + xml_epsilon * tan(rinfo.tiltAngle * M_PI / 180.0);
	      }
	      else {
		shape.rmin1 = rinfo.rminatzmin - 2 * shape.dz * tan(rinfo.tiltAngle * M_PI / 180.0) - xml_epsilon * tan(rinfo.tiltAngle * M_PI / 180.0);
		shape.rmax1 = rinfo.rmaxatzmax + xml_epsilon * tan(rinfo.tiltAngle * M_PI / 180.0);
		shape.rmin2 = rinfo.rminatzmin - xml_epsilon * tan(rinfo.tiltAngle * M_PI / 180.0);
		shape.rmax2 = rinfo.rmaxatzmax + 2 * shape.dz * tan(rinfo.tiltAngle * M_PI / 180.0) + xml_epsilon * tan(rinfo.tiltAngle * M_PI / 180.0);
	      }
	      s.push_back(shape);

	      // reset
	      shape.rmin1 = 0.0;
	      shape.rmax1 = 0.0;
	      shape.rmin2 = 0.0;
	      shape.rmax2 = 0.0;

	      // section of tub
	      shape.type = tb;
	      shape.name_tag = rinfo.name + "Tub";
	      shape.dz = (rinfo.zmax - rinfo.zmin) / 2. + xml_epsilon;
	      shape.rmin = rinfo.rmin - xml_epsilon;
	      shape.rmax = rinfo.rmax + xml_epsilon;
	      s.push_back(shape);

	      // intersection of sections of cone and tub
	      // Please note that the layer's dimensions rely on the fact this intersection is made,
	      // so that layer's extrema are ~rmin and ~rmax
	      shapeOp.name_tag = rinfo.name;
	      shapeOp.type = intersec;
	      shapeOp.rSolid1 = rinfo.name + "Cone";
	      shapeOp.rSolid2 = rinfo.name + "Tub";
	      so.push_back(shapeOp);

	      logic.name_tag = rinfo.name;
	      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
	      logic.material_tag = xml_material_air;
	      l.push_back(logic);
	      
	      pos.parent_tag = trackerXmlTags.nspace + ":" + lname.str();
	      pos.child_tag = trackerXmlTags.nspace + ":" + rinfo.name;
	      pos.trans.dz = (rinfo.z1 + rinfo.z2) / 2.0; 
	      p.push_back(pos);
	      
	      rspec.partselectors.push_back(rinfo.name);
	      //rspec.moduletypes.push_back(minfo_zero);
	      trspec.partselectors.push_back(rinfo.name);
	      //trspec.moduletypes.push_back(minfo_zero);
	      
	      // backward part of the ring
	      alg.name = xml_trackerring_algo;
	      alg.parent = trackerXmlTags.nspace + ":" + rinfo.name;
	      alg.parameters.push_back(stringParam(xml_childparam, trackerXmlTags.nspace + ":" + rinfo.childname));
	      pconverter << (rinfo.modules / 2);
	      alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
	      pconverter.str("");
	      alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
	      alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
	      alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
	      pconverter << rinfo.startPhiAngle1 * 180. / M_PI << "*deg";
	      alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
	      pconverter.str("");
	      pconverter << rinfo.r1;
	      alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
	      pconverter.str("");
	      alg.parameters.push_back(vectorParam(0, 0, (rinfo.z1 - rinfo.z2) / 2.0));	      
	      pconverter << rinfo.isZPlus;
	      alg.parameters.push_back(numericParam(xml_iszplus, pconverter.str()));
	      pconverter.str("");
	      pconverter << rinfo.tiltAngle << "*deg";
	      alg.parameters.push_back(numericParam(xml_tiltangle, pconverter.str()));
	      pconverter.str("");
	      pconverter << rinfo.bw_flipped;
	      alg.parameters.push_back(numericParam(xml_isflipped, pconverter.str()));
	      pconverter.str("");
	      a.push_back(alg);
	      alg.parameters.clear();
	      
	      // forward part of the ring
	      alg.name =  xml_trackerring_algo;
	      alg.parent = trackerXmlTags.nspace + ":" + rinfo.name;
	      alg.parameters.push_back(stringParam(xml_childparam, trackerXmlTags.nspace + ":" + rinfo.childname));
	      pconverter << (rinfo.modules / 2);
	      alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
	      pconverter.str("");
	      alg.parameters.push_back(numericParam(xml_startcopyno, "2"));
	      alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
	      alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
	      pconverter << rinfo.startPhiAngle2 * 180. / M_PI << "*deg";
	      alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
	      pconverter.str("");
	      pconverter << rinfo.r2;
	      alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
	      pconverter.str("");
	      alg.parameters.push_back(vectorParam(0, 0, (rinfo.z2 - rinfo.z1) / 2.0));
	      pconverter << rinfo.isZPlus;
	      alg.parameters.push_back(numericParam(xml_iszplus, pconverter.str()));
	      pconverter.str("");
	      pconverter << rinfo.tiltAngle << "*deg";
	      alg.parameters.push_back(numericParam(xml_tiltangle, pconverter.str()));
	      pconverter.str("");
	      pconverter << rinfo.fw_flipped;
	      alg.parameters.push_back(numericParam(xml_isflipped, pconverter.str()));
	      pconverter.str("");
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
      shape.rmin = rmin - 2 * xml_epsilon;
      shape.rmax = rmax + 2 * xml_epsilon;
      shape.dz = zmax + 2 * xml_epsilon;
      s.push_back(shape);
      logic.name_tag = lname.str();
      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
      l.push_back(logic);
      pos.parent_tag = xml_pixbarident + ":" + trackerXmlTags.bar;
      pos.child_tag = trackerXmlTags.nspace + ":" + lname.str();
      p.push_back(pos);
      lspec.partselectors.push_back(lname.str());
      //lspec.moduletypes.push_back("");
      lspec.moduletypes.push_back(minfo_zero);

      layer++;
    }
    if (!lspec.partselectors.empty()) t.push_back(lspec);
    if (!rspec.partselectors.empty()) t.push_back(rspec); 
    if (!srspec.partselectors.empty()) t.push_back(srspec);
    if (!trspec.partselectors.empty()) t.push_back(trspec);
    if (!sspec.partselectors.empty()) t.push_back(sspec);
    if (!mspec.partselectors.empty()) t.push_back(mspec);
  }
  
  /**
   * This is one of the two main analysis functions that provide the core functionality of this class. It examines the endcap discs in z+
   * and the rings and modules within, extracting a great range of different pieces of information from the geometry layout. These
   * are shapes for individual modules, but also for their enclosing volumes, divided into rings and then discs. They form hierarchies
   * of volumes, one inside the other.
   * Output information are volume hierarchy, material, shapes, positioning (potential use of algorithm and rotations). 
   * They is also some topology, such as which volumes contain the active surfaces, and how those active surfaces are subdivided
   * and connected to the readout electronics. Last but not least, overall radiation and interaction lengths for each layer are 
   * calculated and stored; those are used as approximative values for certain CMSSW functions later on.
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
  void Extractor::analyseDiscs(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& ec, Tracker& tr, XmlTags& trackerXmlTags,
                               std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<ShapeOperationInfo>& so,
			       std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::map<std::string,Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, bool wt) {

    bool isPixelTracker = tr.isPixelTracker();

    // Container inits
    ShapeInfo shape;
    shape.dyy = 0.0;

    ShapeOperationInfo shapeOp;

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
    SpecParInfo rocdims, dspec, rspec, sspec, mspec;
    // Disc
    dspec.name = trackerXmlTags.topo_disc_name + xml_par_tail;
    dspec.parameter.first = xml_tkddd_structure;
    dspec.parameter.second = trackerXmlTags.topo_disc_value;
    // Ring
    rspec.name = trackerXmlTags.topo_ring_name + xml_par_tail;
    rspec.parameter.first = xml_tkddd_structure;
    rspec.parameter.second = trackerXmlTags.topo_ring_value;
    // Module stack
    sspec.name = trackerXmlTags.topo_emodule_name + xml_par_tail;
    sspec.parameter.first = xml_tkddd_structure;
    sspec.parameter.second = trackerXmlTags.topo_emodule_value;
    // Module detectors
    mspec.name = xml_subdet_tiddet + xml_par_tail;
    mspec.parameter.first = xml_tkddd_structure;
    mspec.parameter.second = xml_det_tiddet;


    // material properties
    RILengthInfo ril;
    ril.barrel = false;
    ril.index = 0;  


    LayerAggregator lagg;
    tr.accept(lagg);

    std::vector<std::vector<ModuleCap> >::iterator oiter;
    std::vector<ModuleCap>::iterator iiter;

    int layer = 1;
    int discNumber = 1;


    // LOOP ON DISKS
    for (oiter = ec.begin(); oiter != ec.end(); oiter++) {

      if (lagg.getEndcapLayers()->at(layer - 1)->minZ() > 0) {
   
	std::set<int> ringsIndexes; // VERY UGLY !! TO DO : IMPLEMENT ringsIndexes() IN CLASS DISK
	for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	  int modRing = iiter->getModule().uniRef().ring;
	  if (ringsIndexes.find(modRing) == ringsIndexes.end()) ringsIndexes.insert(modRing);
	}

	// Calculate z extrema of the disk, and diskThickness
	// r extrema of disk and ring
	double rmin = std::numeric_limits<double>::max();
	double rmax = 0;
	// z extrema of disk
	double zmin = std::numeric_limits<double>::max();
	double zmax = 0;
	// z extrema of ring
	std::map<int,double> ringzmin, ringzmax;
	for (const auto& i : ringsIndexes) {
	  ringzmin.insert( {i, std::numeric_limits<double>::max()} );
	  ringzmax.insert( {i, 0.} );
	}

	// loop on module caps
	for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	  if (iiter->getModule().uniRef().side > 0 && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)) {
	    int modRing = iiter->getModule().uniRef().ring;
	    //disk name
	    std::ostringstream dname;
	    dname << xml_disc << discNumber; // e.g. Disc6
	    // module name
	    std::ostringstream mname;
	    mname << xml_endcap_module << modRing << dname.str(); // e.g. EModule1Disc6
	    // parent module name
	    std::string parentName = mname.str();
	    // build module volumes, with hybrids taken into account
	    ModuleComplex modcomplex(mname.str(),parentName,*iiter);
	    modcomplex.buildSubVolumes();
	    rmin = MIN(rmin, modcomplex.getRmin());
	    rmax = MAX(rmax, modcomplex.getRmax());
	    zmin = MIN(zmin, modcomplex.getZmin());
	    zmax = MAX(zmax, modcomplex.getZmax());
	    ringzmin.at(modRing) = MIN(ringzmin.at(modRing), modcomplex.getZmin());  
	    ringzmax.at(modRing) = MAX(ringzmax.at(modRing), modcomplex.getZmax());
	  }
	}

	double diskZ = 0;
	if (ringsIndexes.size() < 2) std::cout << "!!!!!!Disk with less than 2 rings, unexpected" << std::endl;
	else {
	  int firstRingIndex = *(ringsIndexes.begin());
	  int secondRingIndex = *ringsIndexes.begin() + 1;
	  diskZ = (ringzmin.at(firstRingIndex) + ringzmax.at(firstRingIndex) + ringzmin.at(secondRingIndex) + ringzmax.at(secondRingIndex)) / 4.;
	}
	
	//double diskThickness = zmax - zmin;
	double diskThickness = 2. * MAX(fabs(zmin - diskZ), fabs(zmax - diskZ));

	//shape.type = tp;
        shape.rmin = 0.0;
        shape.rmax = 0.0;
        pos.trans.dz = 0.0;
	shapeOp.trans.dz = 0.0;

	// for material properties
        double rtotal = 0.0, itotal = 0.0;
        int count = 0;
	ril.index = discNumber;


	//if (zmin > 0) {	
        std::ostringstream dname, pconverter;
	//disk name
        dname << xml_disc << discNumber; // e.g. Disc6

        std::map<int, ERingInfo> rinfo;
	std::set<int> ridx;


        // LOOP ON MODULE CAPS
        for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	  if (iiter->getModule().uniRef().side > 0 && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)) {
	    // ring number
	    int modRing = iiter->getModule().uniRef().ring;

	    //std::cout << "iiter->getModule().uniRef().phi = " << iiter->getModule().uniRef().phi << " iiter->getModule().uniRef().ring = " << iiter->getModule().uniRef().ring << "iiter->getModule().center().Phi() = " << iiter->getModule().center().Phi() * 180. / M_PI << " iiter->getModule().center().Rho() = " << iiter->getModule().center().Rho() << " iiter->getModule().center().X() = " << iiter->getModule().center().X() << " iiter->getModule().center().Y() = " << iiter->getModule().center().Y() << " iiter->getModule().center().Z() = " << iiter->getModule().center().Z() << " iiter->getModule().flipped() = " << iiter->getModule().flipped() << " iiter->getModule().moduleType() = " << iiter->getModule().moduleType() << std::endl;


	    if (iiter->getModule().uniRef().phi == 1) {
	      // new ring
	      //if (ridx.find(modRing) == ridx.end()) {
	      ridx.insert(modRing);

	      std::ostringstream matname, rname, mname, specname;
	      // ring name
	      rname << xml_ring << modRing << dname.str(); // e.g. Ring1Disc6
	      // module name
	      mname << xml_endcap_module << modRing << dname.str(); // e.g. EModule1Disc6
 
	      // parent module name
	      std::string parentName = mname.str();

	      // build module volumes, with hybrids taken into account
	      ModuleComplex modcomplex(mname.str(),parentName,*iiter);
	      modcomplex.buildSubVolumes();          
#ifdef __DEBUGPRINT__
	      modcomplex.print();
#endif


	      // MODULE

	      // module box
	      shape.name_tag = mname.str();
	      shape.type = iiter->getModule().shape() == RECTANGULAR ? bx : tp;
	      //shape.dx = iiter->getModule().minWidth() / 2.0;
	      //shape.dxx = iiter->getModule().maxWidth() / 2.0;
	      //shape.dy = iiter->getModule().length() / 2.0;
	      //shape.dyy = iiter->getModule().length() / 2.0;
	      //shape.dz = iiter->getModule().thickness() / 2.0;    
	      if (shape.type==bx) {
		shape.dx = modcomplex.getExpandedModuleWidth()/2.0;
		shape.dy = modcomplex.getExpandedModuleLength()/2.0;
		shape.dz = modcomplex.getExpandedModuleThickness()/2.0;
	      } else { // obsolete !
		shape.dx = iiter->getModule().minWidth() / 2.0 + iiter->getModule().serviceHybridWidth();
		shape.dxx = iiter->getModule().maxWidth() / 2.0 + iiter->getModule().serviceHybridWidth();
		shape.dy = iiter->getModule().length() / 2.0 + iiter->getModule().frontEndHybridWidth();
		shape.dyy = iiter->getModule().length() / 2.0 + iiter->getModule().frontEndHybridWidth();
		shape.dz = iiter->getModule().thickness() / 2.0 + iiter->getModule().supportPlateThickness();
	      }
	      s.push_back(shape);

	      // Get it back for sensors
	      shape.dx = iiter->getModule().minWidth() / 2.0;
	      shape.dxx = iiter->getModule().maxWidth() / 2.0;
	      shape.dy = iiter->getModule().length() / 2.0;
	      shape.dyy = iiter->getModule().length() / 2.0;
	      shape.dz = iiter->getModule().thickness() / 2.0;

	      logic.name_tag = mname.str();
	      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;

	      //logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
	      logic.material_tag = xml_material_air;
	      l.push_back(logic);
	      // module composite material
	      //matname << xml_base_actcomp << "D" << discNumber << "R" << modRing;
	      //c.push_back(createComposite(matname.str(), compositeDensity(*iiter, true), *iiter, true));

	    //Topology
	    sspec.partselectors.push_back(mname.str());
            sspec.moduletypes.push_back(minfo_zero);



	      // WAFER -- same x and y size of parent shape, but different thickness
	    string xml_base_lowerupper = "";
            if (iiter->getModule().numSensors() == 2) {
	      xml_base_lowerupper = xml_base_lower;
	      shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_waf;
	    }
	    else {
	      if (iiter->getModule().isPixelModule()) shape.name_tag = mname.str() + xml_PX + xml_base_waf;
	      else if (iiter->getModule().isTimingModule()) shape.name_tag = mname.str() + xml_timing + xml_base_waf;
	      else { std::cerr << "Wafer : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	    }     

	      shape.dz = iiter->getModule().sensorThickness() / 2.0; // CUIDADO WAS calculateSensorThickness(*iiter, mt) / 2.0;
	      //if (iiter->getModule().numSensors() == 2) shape.dz = shape.dz / 2.0; // CUIDADO calcSensThick returned 2x what getSensThick returns, it means that now one-sided sensors are half as thick if not compensated for in the config files
	      s.push_back(shape);

	      logic.name_tag = shape.name_tag;
	      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
	      logic.material_tag = xml_material_air;
	      l.push_back(logic);

	      pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str();
	      pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
	      if (iiter->getModule().uniRef().side > 0) pos.trans.dz = /*shape.dz*/ - iiter->getModule().dsDistance() / 2.0; // CUIDADO WAS getModule().moduleThickness()
	      else pos.trans.dz = iiter->getModule().dsDistance() / 2.0 /*- shape.dz*/; // DITTO HERE
	      p.push_back(pos);

	      if (iiter->getModule().numSensors() == 2) {

		xml_base_lowerupper = xml_base_upper;

		//pos.parent_tag = logic.shape_tag;

		shape.name_tag = mname.str() + xml_base_lowerupper+ xml_base_waf;
		s.push_back(shape);

		logic.name_tag = shape.name_tag;
		logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
		l.push_back(logic);

		pos.child_tag = logic.shape_tag;

		if (iiter->getModule().uniRef().side > 0) pos.trans.dz = /*pos.trans.dz + 2 * shape.dz +*/  iiter->getModule().dsDistance() / 2.0; // CUIDADO removed pos.trans.dz + 2*shape.dz, added / 2.0
		else pos.trans.dz = /* pos.trans.dz - 2 * shape.dz -*/ - iiter->getModule().dsDistance() / 2.0;
		//pos.copy = 2;
		if (iiter->getModule().stereoRotation() != 0) {
		  rot.name = type_stereo + xml_endcap_module + mname.str();
		  rot.thetax = 90.0;
		  rot.phix = iiter->getModule().stereoRotation() / M_PI * 180.;
		  rot.thetay = 90.0;
		  rot.phiy = 90.0 + iiter->getModule().stereoRotation() / M_PI * 180.;
		  r.insert(std::pair<const std::string,Rotation>(rot.name,rot));
		  pos.rotref = trackerXmlTags.nspace + ":" + rot.name;
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


	      // ACTIVE SURFACE
	      xml_base_lowerupper = "";
	      if (iiter->getModule().numSensors() == 2) xml_base_lowerupper = xml_base_lower; 
	      
	      if (iiter->getModule().moduleType() == "ptPS") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_ps + xml_base_pixel + xml_base_act;
	      else if (iiter->getModule().moduleType() == "pt2S") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_2s+ xml_base_act;
	      else if (iiter->getModule().isTimingModule()) shape.name_tag = mname.str() + xml_timing + xml_base_act;
	      else if (iiter->getModule().isPixelModule()) shape.name_tag = mname.str() + xml_PX + xml_base_Act;
	      else { std::cerr << "Active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	      s.push_back(shape);

	      logic.name_tag = shape.name_tag;
	      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
	      logic.material_tag = xml_fileident + ":" + xml_tkLayout_material + xml_sensor_silicon;
	      l.push_back(logic);

	      if (iiter->getModule().numSensors() == 2) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_base_lowerupper + xml_base_waf;
	      else {
		if (iiter->getModule().isTimingModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_timing + xml_base_waf;
		else if (iiter->getModule().isPixelModule())  pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_PX + xml_base_waf;
		else { std::cerr << "Positioning active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	      }

	      pos.child_tag = logic.shape_tag;
	      pos.trans.dz = 0.0;
#ifdef __FLIPSENSORS_IN__ // Flip INNER sensors
	      pos.rotref = trackerXmlTags.nspace + ":" + rot_sensor_tag;
#endif
	      p.push_back(pos);

	      // Topology
	      mspec.partselectors.push_back(logic.name_tag);
	      minfo.name		= iiter->getModule().moduleType();
	      minfo.rocrows	= any2str<int>(iiter->getModule().innerSensor().numROCRows());
	      minfo.roccols	= any2str<int>(iiter->getModule().innerSensor().numROCCols());
	      minfo.rocx		= any2str<int>(iiter->getModule().innerSensor().numROCX());
	      minfo.rocy		= any2str<int>(iiter->getModule().innerSensor().numROCY());
	      mspec.moduletypes.push_back(minfo);

	      if (iiter->getModule().numSensors() == 2) {

		xml_base_lowerupper = xml_base_upper;

		if (iiter->getModule().moduleType() == "ptPS") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_ps + xml_base_strip + xml_base_act;
		else if (iiter->getModule().moduleType() == "pt2S") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_2s+ xml_base_act;
		else { std::cerr << "Active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
		s.push_back(shape);

		logic.name_tag = shape.name_tag;
		logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
		logic.material_tag = xml_fileident + ":" + xml_tkLayout_material + xml_sensor_silicon;
		l.push_back(logic);

		pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_base_lowerupper + xml_base_waf;
		pos.child_tag = logic.shape_tag;
		pos.trans.dz = 0.0;
#ifdef __FLIPSENSORS_OUT__ // Flip OUTER sensors
		pos.rotref = trackerXmlTags.nspace + ":" + rot_sensor_tag;
#endif
		p.push_back(pos);

		// Topology
		mspec.partselectors.push_back(logic.name_tag);
		minfo.rocrows	= any2str<int>(iiter->getModule().outerSensor().numROCRows());
		minfo.roccols	= any2str<int>(iiter->getModule().outerSensor().numROCCols());
		minfo.rocx		= any2str<int>(iiter->getModule().outerSensor().numROCX());
		minfo.rocy		= any2str<int>(iiter->getModule().outerSensor().numROCY());
		mspec.moduletypes.push_back(minfo);
	      }

	      modcomplex.addMaterialInfo(c);
	      modcomplex.addShapeInfo(s);
	      modcomplex.addLogicInfo(l);
	      modcomplex.addPositionInfo(p);
#ifdef __DEBUGPRINT__
	      modcomplex.print();
#endif
	      


	      // collect ring info
	      ERingInfo rinf;
	      rinf.name = rname.str();
	      rinf.childname = mname.str();
	      rinf.fw = (iiter->getModule().center().Z() > (zmin + zmax) / 2.0);
	      rinf.isZPlus = iiter->getModule().uniRef().side;
	      rinf.fw_flipped = iiter->getModule().flipped();	      
	      rinf.modules = lagg.getEndcapLayers()->at(layer - 1)->ringsMap().at(modRing)->numModules();
	      rinf.mthk = modcomplex.getExpandedModuleThickness();
	      rinf.rmin  = modcomplex.getRmin();
	      rinf.rmid = iiter->getModule().center().Rho();
	      rinf.rmax = modcomplex.getRmax();
	      rinf.zmin = ringzmin.at(modRing);
	      rinf.zmax = ringzmax.at(modRing);
	      rinf.zfw = iiter->getModule().center().Z();
	      rinf.startPhiAnglefw = iiter->getModule().center().Phi();
	      rinfo.insert(std::pair<int, ERingInfo>(modRing, rinf));


	      // material properties
	      rtotal = rtotal + iiter->getRadiationLength();
	      itotal = itotal + iiter->getInteractionLength();
	      count++;
	    }

	    if (iiter->getModule().uniRef().phi == 2) {
	      std::map<int,ERingInfo>::iterator it;
	      // fill the info of the z-backward part of the ring with matching ring number
	      it = rinfo.find(modRing);
	      if (it != rinfo.end()) {
		it->second.zbw = iiter->getModule().center().Z();
		it->second.startPhiAnglebw = iiter->getModule().center().Phi();
	      }
	    }
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
        //findDeltaZ(lagg.getEndcapLayers()->at(layer - 1)->getModuleVector()->begin(), // CUIDADO what the hell is this??
        //lagg.getEndcapLayers()->at(layer - 1)->getModuleVector()->end(), (zmin + zmax) / 2.0) / 2.0;

        std::set<int>::const_iterator siter, sguard = ridx.end();
        for (siter = ridx.begin(); siter != sguard; siter++) {
          if (rinfo[*siter].modules > 0) {

            shape.name_tag = rinfo[*siter].name;
            shape.rmin = rinfo[*siter].rmin - xml_epsilon;
            shape.rmax = rinfo[*siter].rmax + xml_epsilon;
	    shape.dz = (rinfo[*siter].zmax - rinfo[*siter].zmin) / 2.0 + xml_epsilon;
            s.push_back(shape);

            logic.name_tag = shape.name_tag;
            logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
            logic.material_tag = xml_material_air;
            l.push_back(logic);

            pos.parent_tag = trackerXmlTags.nspace + ":" + dname.str(); // CUIDADO ended with: + xml_plus;
            pos.child_tag = logic.shape_tag;

	    pos.trans.dz = (rinfo[*siter].zmin + rinfo[*siter].zmax) / 2.0 - diskZ;
            p.push_back(pos);
            //pos.parent_tag = trackerXmlTags.nspace + ":" + dname.str(); // CUIDADO ended with: + xml_minus;
            //p.push_back(pos);

            rspec.partselectors.push_back(logic.name_tag);
            rspec.moduletypes.push_back(minfo_zero);

	    // forward part of the ring
	    alg.name = xml_trackerring_algo;
            alg.parent = logic.shape_tag;
            alg.parameters.push_back(stringParam(xml_childparam, trackerXmlTags.nspace + ":" + rinfo[*siter].childname));
            pconverter << (rinfo[*siter].modules / 2);
            alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
            pconverter.str("");
            alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
            alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
            alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
            pconverter << rinfo[*siter].startPhiAnglefw * 180. / M_PI << "*deg";
            alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
            pconverter.str("");
            pconverter << rinfo[*siter].rmid;
            alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
            pconverter.str("");
	    alg.parameters.push_back(vectorParam(0, 0, rinfo[*siter].zfw - (rinfo[*siter].zmin + rinfo[*siter].zmax) / 2.0));
	    pconverter << rinfo[*siter].isZPlus;
	    alg.parameters.push_back(numericParam(xml_iszplus, pconverter.str()));
	    pconverter.str("");
	    alg.parameters.push_back(numericParam(xml_tiltangle, "90*deg"));
	    pconverter << rinfo[*siter].fw_flipped;
	    alg.parameters.push_back(numericParam(xml_isflipped, pconverter.str()));
	    pconverter.str("");
            a.push_back(alg);
            alg.parameters.clear();

	    // backward part of the ring
	    alg.name = xml_trackerring_algo;
            alg.parameters.push_back(stringParam(xml_childparam, trackerXmlTags.nspace + ":" + rinfo[*siter].childname));
            pconverter << (rinfo[*siter].modules / 2);
            alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
            pconverter.str("");
            alg.parameters.push_back(numericParam(xml_startcopyno, "2"));
            alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
            alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
            pconverter << rinfo[*siter].startPhiAnglebw * 180. / M_PI << "*deg";
            alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
            pconverter.str("");
            pconverter << rinfo[*siter].rmid;
            alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
            pconverter.str("");
	    alg.parameters.push_back(vectorParam(0, 0, rinfo[*siter].zbw - (rinfo[*siter].zmin + rinfo[*siter].zmax) / 2.0));
	    pconverter << rinfo[*siter].isZPlus;
	    alg.parameters.push_back(numericParam(xml_iszplus, pconverter.str()));
	    pconverter.str("");
	    alg.parameters.push_back(numericParam(xml_tiltangle, "90*deg"));
	    pconverter << !rinfo[*siter].fw_flipped;
	    alg.parameters.push_back(numericParam(xml_isflipped, pconverter.str()));
	    pconverter.str("");
            a.push_back(alg);
            alg.parameters.clear();
          }
        }

        //disc
        shape.name_tag = dname.str();
        shape.rmin = rmin - 2 * xml_epsilon;
        shape.rmax = rmax + 2 * xml_epsilon;
        shape.dz = diskThickness / 2.0 + 2 * xml_epsilon; //(zmax - zmin) / 2.0;
        s.push_back(shape);

	shape.name_tag = dname.str();
        logic.name_tag = shape.name_tag; // CUIDADO ended with + xml_plus;
        logic.shape_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
        logic.material_tag = xml_material_air;
        l.push_back(logic);

        pos.parent_tag = xml_pixfwdident + ":" + trackerXmlTags.fwd;
        pos.child_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
        //pos.trans.dz = (zmax + zmin) / 2.0 - xml_z_pixfwd;
	pos.trans.dz = diskZ - xml_z_pixfwd;
        p.push_back(pos);

        dspec.partselectors.push_back(logic.name_tag);
        dspec.moduletypes.push_back(minfo_zero);
        dspec.partextras.push_back(logic.extra);
        //   logic.name_tag = shape.name_tag; // CUIDADO ended with + xml_minus;
        //   logic.extra = xml_minus;
        //   l.push_back(logic);
        //   pos.parent_tag = xml_pixfwdident + ":" + trackerXmlTags.fwd;
        //   pos.child_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
        //   p.push_back(pos);
        //dspec.partselectors.push_back(logic.name_tag); // CUIDADO dspec still needs to be duplicated for minus discs (I think)
        //dspec.partextras.push_back(logic.extra);
	discNumber++;
      }
      layer++;
    }
    if (!dspec.partselectors.empty()) t.push_back(dspec);
    if (!rspec.partselectors.empty()) t.push_back(rspec);
    if (!sspec.partselectors.empty()) t.push_back(sspec);
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
  void Extractor::analyseServices(InactiveSurfaces& is, bool& isPixelTracker, XmlTags& trackerXmlTags,
					std::vector<Composite>& c, std::vector<LogicalInfo>& l,
                                        std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt) {
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
    double serviceBarrelRMin = std::numeric_limits<double>::max();
    double serviceBarrelRMax = 0;
    double serviceEndcapsRMin = std::numeric_limits<double>::max();
    double serviceEndcapsRMax = 0;
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
      matname << xml_base_serfcomp << "R" << (int)(iter->getInnerRadius()) << "dZ" << (int)(iter->getZLength());
      shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(iter->getZOffset());
      if ((iter->getZOffset() + iter->getZLength()) > 0) c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));
#else
      matname << xml_base_serfcomp << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
      shapename << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0));
      if ((iter->getZOffset() + iter->getZLength()) > 0 ) {
        if ( iter->getLocalMasses().size() ) {
          c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));

	  // TO DO : CALCULATION OF OUTERMOST SHAPES BOUNDARIES
	  double startEndcaps;
	  if (!isPixelTracker) startEndcaps = 1250.;
	  //else startEndcaps = 300.; // PIXEL 1_1_1
	  else startEndcaps = 227.;   // PIXEL 4_0_2_1
          
	  // BARREL services
	  if ((iter->getZOffset() + iter->getZLength() / 2.0) < startEndcaps ) {
	    shape.name_tag = shapename.str();
	    shape.dz = iter->getZLength() / 2.0;
	    shape.rmin = iter->getInnerRadius();
	    shape.rmax = shape.rmin + iter->getRWidth();
	    s.push_back(shape);
	    if (shape.rmin < serviceBarrelRMin) serviceBarrelRMin = shape.rmin;
	    if (shape.rmax > serviceBarrelRMax) serviceBarrelRMax = shape.rmax;

	    logic.name_tag = shapename.str();
	    logic.shape_tag = trackerXmlTags.nspace + ":" + shapename.str();
	    logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
	    l.push_back(logic);

	    pos.parent_tag = xml_pixbarident + ":" + trackerXmlTags.bar; //xml_tracker;
	    pos.child_tag = logic.shape_tag;
	    pos.trans.dz = iter->getZOffset() + shape.dz;
	    p.push_back(pos);
	    pos.copy = 2;
	    pos.trans.dz = -pos.trans.dz;
	    pos.rotref = trackerXmlTags.nspace + ":" + xml_flip_mod_rot;
	    p.push_back(pos);
	  }

	  // ENDCAPS services
	  else {
	    // cut in 2 the services that belong to both Barrel and Endcaps mother volumes
	    if (iter->getZOffset() < startEndcaps) {
	      std::ostringstream shapenameBarrel, shapenameEndcaps;
	      shapenameBarrel << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0)) << "BarrelPart";
	      shapenameEndcaps << xml_base_serf << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset() + iter->getZLength() / 2.0)) << "EndcapsPart";

	      // Barrel part
	      shape.name_tag = shapenameBarrel.str();
	      shape.dz = (startEndcaps - iter->getZOffset()) / 2.0 - xml_epsilon;
	      shape.rmin = iter->getInnerRadius();
	      shape.rmax = shape.rmin + iter->getRWidth();
	      s.push_back(shape);
	      if (shape.rmin < serviceBarrelRMin) serviceBarrelRMin = shape.rmin;
	      if (shape.rmax > serviceBarrelRMax) serviceBarrelRMax = shape.rmax;

	      logic.name_tag = shapenameBarrel.str();
	      logic.shape_tag = trackerXmlTags.nspace + ":" + shapenameBarrel.str();
	      logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
	      l.push_back(logic);

	      pos.parent_tag = xml_pixbarident + ":" + trackerXmlTags.bar; //xml_tracker;
	      pos.child_tag = logic.shape_tag;
	      pos.trans.dz = iter->getZOffset() + shape.dz + xml_epsilon;
	      p.push_back(pos);
	      pos.copy = 2;
	      pos.trans.dz = -pos.trans.dz;
	      pos.rotref = trackerXmlTags.nspace + ":" + xml_flip_mod_rot;
	      p.push_back(pos);
	      

	      pos.copy = 1;
	      pos.rotref.clear();

	      // Endcaps part
	      shape.name_tag = shapenameEndcaps.str();
	      shape.dz = (iter->getZOffset() + iter->getZLength() - startEndcaps) / 2.0 - xml_epsilon;
	      shape.rmin = iter->getInnerRadius();
	      shape.rmax = shape.rmin + iter->getRWidth();
	      s.push_back(shape);    
	      if (shape.rmin < serviceEndcapsRMin) serviceEndcapsRMin = shape.rmin;
	      if (shape.rmax > serviceEndcapsRMax) serviceEndcapsRMax = shape.rmax;

	      logic.name_tag = shapenameEndcaps.str();
	      logic.shape_tag = trackerXmlTags.nspace + ":" + shapenameEndcaps.str();
	      logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
	      l.push_back(logic);

	      pos.parent_tag = xml_pixfwdident + ":" + trackerXmlTags.fwd; // xml_tracker;
	      pos.child_tag = logic.shape_tag;
	      pos.trans.dz = startEndcaps + shape.dz + xml_epsilon - xml_z_pixfwd;
	      p.push_back(pos);


	    }

	    // ENDCAPS-only services
	    else {
	      shape.name_tag = shapename.str();
	      shape.dz = iter->getZLength() / 2.0;
	      shape.rmin = iter->getInnerRadius();
	      shape.rmax = shape.rmin + iter->getRWidth();
	      s.push_back(shape);
	      if (shape.rmin < serviceEndcapsRMin) serviceEndcapsRMin = shape.rmin;
	      if (shape.rmax > serviceEndcapsRMax) serviceEndcapsRMax = shape.rmax;

	      logic.name_tag = shapename.str();
	      logic.shape_tag = trackerXmlTags.nspace + ":" + shapename.str();
	      logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
	      l.push_back(logic);

	      pos.parent_tag = xml_pixfwdident + ":" + trackerXmlTags.fwd; // xml_tracker;
	      pos.child_tag = logic.shape_tag;
	      pos.trans.dz = iter->getZOffset() + shape.dz - xml_z_pixfwd;
	      p.push_back(pos);
	    }
	  }



	  pos.copy = 1;
	  pos.rotref.clear();

        } 
	else {
          std::stringstream msg;
          msg << shapename.str() << " is not exported to XML because it is empty." << std::ends;
          logWARNING( msg.str() ); 
        }
      }
#endif
    }
    // DEBUG EXTREMA
    /*std::cout << "serviceBarrelRMin = " << serviceBarrelRMin << std::endl;
    std::cout << "serviceBarrelRMax = " << serviceBarrelRMax << std::endl;
    std::cout << "serviceEndcapsRMin = " << serviceEndcapsRMin << std::endl;
    std::cout << "serviceEndcapsRMax = " << serviceEndcapsRMax << std::endl;*/
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
  void Extractor::analyseSupports(InactiveSurfaces& is, bool& isPixelTracker, XmlTags& trackerXmlTags,
				  std::vector<Composite>& c, std::vector<LogicalInfo>& l,
                                  std::vector<ShapeInfo>& s, std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt) {
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

    double supportBarrelRMin = std::numeric_limits<double>::max();
    double supportBarrelRMax = 0;
    double supportEndcapsRMin = std::numeric_limits<double>::max();
    double supportEndcapsRMax = 0;
    // support volume loop
    for (iter = sp.begin(); iter != guard; iter++) {
      std::ostringstream matname, shapename;
      matname << xml_base_lazycomp << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(iter->getZLength() / 2.0 + iter->getZOffset());
#if 0
      shapename << xml_base_lazy /*<< any2str(iter->getCategory()) */<< "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(fabs(iter->getZOffset()));
#else
      shapename << xml_base_lazy /*<< any2str(iter->getCategory()) */<< "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(iter->getZLength() / 2.0 + iter->getZOffset());
#endif

      if ((iter->getZOffset() + iter->getZLength()) > 0 ) {
        if ( iter->getLocalMasses().size() ) {
          c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter));

	  // TO DO : CALCULATION OF OUTERMOST SHAPES BOUNDARIES
	  double startEndcaps;
	  if (!isPixelTracker) startEndcaps = 1250.;
	  //else startEndcaps = 300.; // PIXEL 1_1_1
	  else startEndcaps = 227.;   // PIXEL 4_0_2_1          

	  // BARREL supports
	  if (iter->getZOffset() < startEndcaps ) {
	    shape.name_tag = shapename.str();
	    shape.dz = iter->getZLength() / 2.0;
	    shape.rmin = iter->getInnerRadius();
	    shape.rmax = shape.rmin + iter->getRWidth();
	    s.push_back(shape);
	    if (shape.rmin < supportBarrelRMin) supportBarrelRMin = shape.rmin;
	    if (shape.rmax > supportBarrelRMax) supportBarrelRMax = shape.rmax;

	    logic.name_tag = shapename.str();
	    logic.shape_tag = trackerXmlTags.nspace + ":" + shapename.str();
	    logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
	    l.push_back(logic);

	    pos.parent_tag = xml_pixbarident + ":" + trackerXmlTags.bar; //xml_tracker;
	    pos.child_tag = logic.shape_tag;
	    pos.trans.dz = iter->getZOffset() + shape.dz;
	    p.push_back(pos);
	    pos.copy = 2;
	    pos.trans.dz = -pos.trans.dz;
	    pos.rotref = trackerXmlTags.nspace + ":" + xml_flip_mod_rot;
	    p.push_back(pos);
	  }

	  // ENDCAPS supports
	  else {
	      shape.name_tag = shapename.str();
	      shape.dz = iter->getZLength() / 2.0;
	      shape.rmin = iter->getInnerRadius();
	      shape.rmax = shape.rmin + iter->getRWidth();
	      s.push_back(shape);
	      if (shape.rmin < supportEndcapsRMin) supportEndcapsRMin = shape.rmin;
	      if (shape.rmax > supportEndcapsRMax) supportEndcapsRMax = shape.rmax;

	      logic.name_tag = shapename.str();
	      logic.shape_tag = trackerXmlTags.nspace + ":" + shapename.str();
	      logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
	      l.push_back(logic);

	      pos.parent_tag = xml_pixfwdident + ":" + trackerXmlTags.fwd; // xml_tracker;
	      pos.child_tag = logic.shape_tag;
	      pos.trans.dz = iter->getZOffset() + shape.dz - xml_z_pixfwd;
	      p.push_back(pos);
	  }


	  pos.copy = 1;
	  pos.rotref.clear();

        } 
	else {
          std::stringstream msg;
          msg << shapename.str() << " is not exported to XML because it is empty." << std::ends;
          logWARNING( msg.str() ); 
        }
      }

    }
    // DEBUG EXTREMA
    /*std::cout << "supportBarrelRMin = " << supportBarrelRMin << std::endl;
    std::cout << "supportBarrelRMax = " << supportBarrelRMax << std::endl;
    std::cout << "supportEndcapsRMin = " << supportEndcapsRMin << std::endl;
    std::cout << "supportEndcapsRMax = " << supportEndcapsRMax << std::endl;*/
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
    for (const auto& it : mp.getLocalMasses()) {   
      if (!nosensors || (it.first.compare(xml_sensor_silicon) != 0)) { // SenSi element is not treated here. 
	if (comp.elements.find(it.first) != comp.elements.end()) { // should not be necessary, but better safe than sorry
	  throw PathfulException("Tried to insert several times " + it.first + " in component " + comp.name);
        } 
	else { 
	  comp.elements.insert(it); // add element to the composite
	  m += it.second; // calculate total mass of the composite
	}      
      }
    }
    // element info is the ratio of the mass of the element by the mass of the composite
    for (auto& elem : comp.elements) {
      elem.second /= m;
    }
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
  // Obsolete and not used.
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
  // Obsolete and not used.
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
   * Calculate estimated atomic number and atomic weight from radiation and nuclear interaction length.
   * @param radiationLength
   * @param interactionLength
   * @return Estimated atomic number and atomic weight.
   */
  std::pair<double, int> Extractor::getAZ(double radiationLength, double interactionLength) {
    const double A_a = 31.645;
    const double A_b = 11.57238;
    double A_estimate = pow((interactionLength-A_b)/A_a,3); // Wait, A value will be corrected in a few lines.
    int Z = findClosest_Z(A_estimate, radiationLength); // Now, we have Z !

    // Correct A to obtain the right radLength, but leave 1/3 of the
    // correction aside, given the relative dependencies
    double radEstimate = compute_X0(Z, A_estimate);   
    A_estimate /= pow(radEstimate/radiationLength, 2/3.); // Now, we have A !

    // Bonus : calculate estimated radiation length and interaction lengths that one can get from the provided Z and A values.
    radEstimate = compute_X0(Z, A_estimate);  
    double intEstimate = pow(A_estimate, 1./3.)*A_a+A_b;
    // On CMSSW side, radiation lengths and nuclear interaction lengths will be recomputed from Z, A and density.
    // Now, a question is whether the computed radiation and interaction lengths on CMMSW side (similar to radEstimate and intEstimate) are closed to the values we initially had at hand (double radiationLength and interactionLength) !
    // The estimated errors can be calculated by :
    // Error on radiation length : (radEstimate-radiationLength)/radiationLength
    // Error on nuclear interaction length : (intEstimate-interactionLength)/interactionLength

    std::pair<double, int> AZ = std::make_pair(A_estimate, Z);
    return AZ;
  }

  /**
   * Calculate the closest atomic number corresponding to a given estimated atomic weight and radiation length.
   * @param pA Estimated atomic weight
   * @param radiationLength
   * @return Closest atomic number.
   */
  int Extractor::findClosest_Z(double pA, double radiationLength) {
    double distance = fabs(radiationLength - compute_X0(1, pA));
    int closest_Z = 1;
    for (int i=2; i<85; ++i) {
      double aDistance = fabs(radiationLength-compute_X0(i, pA));
      if (aDistance<distance) distance = aDistance , closest_Z=i;
    }
    return closest_Z;
  }

  /**
   * Calculate the estimated radiation length associated to an estimated atomic number Z and estimated atomic weight A.
   */
  double Extractor::compute_X0(int pZ, double pA) {
    double Z = pZ;
    const double alpha = 1/137.035999139;
    double a = alpha*Z;
    double a_2 = a*a;
    double a_3 = a*a*a;
    double a_4 = a_2*a_2;
    double a_6 = a_3*a_3;
    double f_Z = a_2*(1/(1+a_2)+0.20206-0.0369*a_3+0.0083*a_4-0.002*a_6);
    double L_Z, L_prime_Z;
    switch (pZ) {
    case 1:
      // H
      L_Z       = 5.31;
      L_prime_Z = 6.144;
      break;
    case 2:
      // He
      L_Z       = 4.79;
      L_prime_Z = 5.621;
      break;
    case 3:
      // Li
      L_Z       = 4.74;
      L_prime_Z = 5.805;
      break;
    case 4:
      // Be:
      L_Z       = 4.71;
      L_prime_Z = 5.924;
      break;
    default:
      // Z >4
      L_Z = log(184.15) - 1/3. * log(Z);
      L_prime_Z = log(1194) - 2/3. * log(Z);
    }
    
    // X_0
    return 716.408 * pA / (Z*Z*(L_Z-f_Z) + Z*L_prime_Z);
  }

 

  const double ModuleComplex::kmm3Tocm3 = 1e-3;

// These value should be consistent with 
// the configuration file
           // OUTER TRACKER MODULE
  const int ModuleComplex::HybridFBLR_0  = 0; // Front + Back + Right + Left
  const int ModuleComplex::InnerSensor   = 1; 
  const int ModuleComplex::OuterSensor   = 2; 
  const int ModuleComplex::HybridFront   = 3; 
  const int ModuleComplex::HybridBack    = 4; 
  const int ModuleComplex::HybridLeft    = 5; 
  const int ModuleComplex::HybridRight   = 6; 
  const int ModuleComplex::HybridBetween = 7; 
  const int ModuleComplex::SupportPlate  = 8; // Support Plate
  // extras
  const int ModuleComplex::HybridFB        = 34; 
  const int ModuleComplex::HybridLR        = 56; 
  const int ModuleComplex::HybridFBLR_3456 = 3456;

           // PIXEL MODULE
  const int ModuleComplex::PixelModuleNull   = 0;
  const int ModuleComplex::PixelModuleHybrid   = 1; 
  const int ModuleComplex::PixelModuleSensor   = 2; 
  const int ModuleComplex::PixelModuleChip     = 3; 

  

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
                                                             chipThickness(module.chipThickness()),
                                                             hybridTotalMass(0.),
                                                             hybridTotalVolume_mm3(-1.),
                                                             hybridFrontAndBackVolume_mm3(-1.),
                                                             hybridLeftAndRightVolume_mm3(-1.),
                                                             moduleMassWithoutSensors_expected(0.),
                                                             expandedModWidth(modWidth+2*serviceHybridWidth),
                                                             expandedModLength(modLength+2*frontEndHybridWidth),
                                                             center(module.center()),
                                                             normal(module.normal()),
                                                             prefix_material(xml_hybrid_comp) {
    if (!module.isPixelModule()) {
      expandedModThickness = sensorDistance + 2.0 * (supportPlateThickness + sensorThickness);   
      prefix_xmlfile = xml_fileident + ":";
      nTypes = 9;
    }
    else {
      expandedModThickness = sensorThickness + 2.0 * MAX(chipThickness, hybridThickness);
      prefix_xmlfile = xml_PX_fileident + ":";
      nTypes = 4;
    }
  }

  ModuleComplex::~ModuleComplex() {
    std::vector<Volume*>::const_iterator vit;
    for ( vit = volumes.begin(); vit != volumes.end(); vit++ ) delete *vit;
  }

  void ModuleComplex::buildSubVolumes() {
    Volume* vol[nTypes];
    if (!module.isPixelModule()) {
      //                                                   OUTER TRACKER MODULE
      //
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
      //  R(6) and L(5) are Front-End Hybrids.
      //  B(4) and F(3) are Service Hybdrids.
      //
      //  SupportPlate(8) thickness is of course null for 2S modules
    
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
    }

    else {
      //                                                      PIXEL MODULE
      //
      //  Top View 
      //        ------------------           y
      //        |                |           ^
      //        |     Hybrid     |           |
      //        |       (1)      |           +----> x
      //        ------------------    
      //                                             z
      //                                             ^
      //  Side View                                  |
      //         ================ Hybrid  (1)        +----> x
      //         ---------------- Sensor  (2)
      //         ================ Chip    (3)
      //
      // Chip(3) volume can contain Bumps and any other material for simplification.

      //Unused pointers
      vol[PixelModuleNull] = 0;

      double dx = modWidth;              
      double dy = modLength; 
      double dz = hybridThickness;  
      double posx = 0.;
      double posy = 0.;
      double posz = sensorThickness / 2. + hybridThickness / 2.;

      // Hybrid Volume (Top Inactive)
      vol[PixelModuleHybrid] = new Volume(moduleId + "Hybrid", PixelModuleHybrid, parentId, dx, dy, dz, posx, posy, posz);

      dz = chipThickness;
      posz = - sensorThickness / 2. - chipThickness / 2.; 
      // Chip Volume (Bottom Inactive)
      vol[PixelModuleChip] = new Volume(moduleId + "Chip", PixelModuleChip, parentId, dx, dy, dz, posx, posy, posz);
    }


    // =========================================================================================================
    // Finding Xmin/Xmax/Ymin/Ymax/Zmin/Zmax/Rmin/Rmax/RminatZmin/RmaxatZmax, taking hybrid volumes into account
    // =========================================================================================================
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
    vector<double> xv; // x list (in global frame of reference) from which we will find min/max.
    vector<double> yv; // y list (in global frame of reference) from which we will find min/max.
    vector<double> zv; // z (in global frame of reference) list from which we will find min/max.
    vector<double> rv; // radius list (in global frame of reference) from which we will find min/max.
    vector<double> ratzminv; // radius list (in global frame of reference) at zmin from which we will find min/max.
    vector<double> ratzmaxv; // radius list (in global frame of reference) at zmax from which we will find min/max.

    // mx: (v2+v3)/2 - center, my: (v1+v2)/2 - center
    XYZVector mx = 0.5*( module.basePoly().getVertex(2) + module.basePoly().getVertex(3) ) - center ;
    XYZVector my = 0.5*( module.basePoly().getVertex(1) + module.basePoly().getVertex(2) ) - center ;

    // new vertexes after expansion due to hybrid volumes
    const int npoints = 5; // v0,v1,v2,v3,v4(=v0)
    XYZVector v[npoints-1];
    v[0] = module.center() - (expandedModWidth/modWidth)*mx - (expandedModLength/modLength)*my;
    v[1] = module.center() - (expandedModWidth/modWidth)*mx + (expandedModLength/modLength)*my;
    v[2] = module.center() + (expandedModWidth/modWidth)*mx + (expandedModLength/modLength)*my;
    v[3] = module.center() + (expandedModWidth/modWidth)*mx - (expandedModLength/modLength)*my;

    // Calculate all vertex candidates (8 points)
    XYZVector v_top[npoints];    // module's top surface
    XYZVector v_bottom[npoints]; // module's bottom surface

    for (int ip = 0; ip < npoints-1; ip++) {
      v_top[ip]    = v[ip] + 0.5*expandedModThickness*normal;
      v_bottom[ip] = v[ip] - 0.5*expandedModThickness*normal;

      // for debuging
      vertex.push_back(v_top[ip]);
      vertex.push_back(v_bottom[ip]);

      // Calculate xmin, xmax, ymin, ymax, zmin, zmax
      xv.push_back(v_top[ip].X());
      xv.push_back(v_bottom[ip].X());
      yv.push_back(v_top[ip].Y());
      yv.push_back(v_bottom[ip].Y());
      zv.push_back(v_top[ip].Z());
      zv.push_back(v_bottom[ip].Z());
    }
    // Find min and max
    xmin = *std::min_element(xv.begin(), xv.end());
    xmax = *std::max_element(xv.begin(), xv.end());
    ymin = *std::min_element(yv.begin(), yv.end());
    ymax = *std::max_element(yv.begin(), yv.end());
    zmin = *std::min_element(zv.begin(), zv.end());
    zmax = *std::max_element(zv.begin(), zv.end());


    // Calculate module's mid-points (8 points)
    XYZVector v_mid_top[npoints]; // module's top surface mid-points
    XYZVector v_mid_bottom[npoints]; // module's bottom surface mid-points

    v_top[npoints-1] = v_top[0]; // copy v0 as v4 for convenience
    v_bottom[npoints-1] = v_bottom[0]; // copy v0 as v4 for convenience

    for (int ip = 0; ip < npoints-1; ip++) {
      v_mid_top[ip] = (v_top[ip] + v_top[ip+1]) / 2.0;
      v_mid_bottom[ip] = (v_bottom[ip] + v_bottom[ip+1]) / 2.0;
    }

    // Calculate rmin, rmax, rminatzmin, rmaxatzmax...
    for (int ip = 0; ip < npoints-1; ip++) {

      // module's bottom surface
      if (fabs(v_bottom[ip].Z() - zmin) < 0.001) {
	v_bottom[ip].SetZ(0.); // projection to xy plan.
	ratzminv.push_back(v_bottom[ip].R());
      }
      if (fabs(v_bottom[ip].Z() - zmax) < 0.001) {
	v_bottom[ip].SetZ(0.); // projection to xy plan.
	ratzmaxv.push_back(v_bottom[ip].R());
      }
      v_bottom[ip].SetZ(0.); // projection to xy plan.
      rv.push_back(v_bottom[ip].R());

      // module's top surface
      if (fabs(v_top[ip].Z() - zmin) < 0.001) {
	v_top[ip].SetZ(0.); // projection to xy plan.
	ratzminv.push_back(v_top[ip].R());
      }
      if (fabs(v_top[ip].Z() - zmax) < 0.001) {
	v_top[ip].SetZ(0.); // projection to xy plan.
	ratzmaxv.push_back(v_top[ip].R());
      }
      v_top[ip].SetZ(0.); // projection to xy plan.
      rv.push_back(v_top[ip].R());
    
      // module's bottom surface mid-points
      if (fabs(v_mid_bottom[ip].Z() - zmin) < 0.001) {
	v_mid_bottom[ip].SetZ(0.); // projection to xy plan.
	ratzminv.push_back(v_mid_bottom[ip].R());
      }
      v_mid_bottom[ip].SetZ(0.); // projection to xy plan.
      rv.push_back(v_mid_bottom[ip].R());
    
      // module's top surface mid-points
      if (fabs(v_mid_top[ip].Z() - zmin) < 0.001) {
	v_mid_top[ip].SetZ(0.); // projection to xy plan.
	ratzminv.push_back(v_mid_top[ip].R());
      }
      v_mid_top[ip].SetZ(0.); // projection to xy plan.
      rv.push_back(v_mid_top[ip].R());
    }
    // Find min and max
    rmin = *std::min_element(rv.begin(), rv.end());
    rmax = *std::max_element(rv.begin(), rv.end());
    rminatzmin = *std::min_element(ratzminv.begin(), ratzminv.end());
    rmaxatzmax = *std::max_element(ratzmaxv.begin(), ratzmaxv.end());


    // Material assignment

    ElementsVector matElements = module.getLocalElements();
    ElementsVector::const_iterator meit;
    for (meit = matElements.begin(); meit != matElements.end(); meit++) {

      const MaterialObject::Element* el = *meit;

      // We skip in the case of ...
      // Filter sensor material and unexpected targetVolumes.
      if ( el->componentName() == "Sensor"      ||
	   el->componentName() == "Sensors"     ||
	   el->componentName() == "PS Sensor"   ||
	   el->componentName() == "PS Sensors"  ||
	   el->componentName() == "2S Sensor"   ||
	   el->componentName() == "2S Sensors"     ) {
	continue; // We will not handle sensors in this class.
	          // Indeed, only Sensi is assigned to the Sensor volume. There is no addition of different materials.
	          // Sensi is directly assigned to the sensor volume within Extractor::analyseLayers and Extractor::analyseDiscs.
      } 

      // OUTER TRACKER MODULE
      if (!module.isPixelModule()) {
	if ( el->targetVolume() == InnerSensor  ||
	     el->targetVolume() == OuterSensor     ) {
	  continue; // Still to skip sensors, in case not detected by the component name.
	} else if ( el->targetVolume() >= nTypes   &&
		    el->targetVolume() != HybridFB &&
		    el->targetVolume() != HybridLR &&
		    el->targetVolume() != HybridFBLR_3456  ) {
	  std::cerr << "!!!! ERROR !!!! : Found unexpected targetVolume." << std::endl;
	  std::cerr << "targetVolume " << el->targetVolume() << " is not supported for Outer Barrel modules. Exit." << std::endl;
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

      // PIXEL MODULE
      else {
	if ( el->targetVolume() == PixelModuleSensor ) {
	  continue; // Still to skip sensors, in case not detected by the component name.
	}
	else if ( el->targetVolume() != PixelModuleHybrid &&
		  el->targetVolume() != PixelModuleChip ) {
	  throw PathfulException("!!!! ERROR !!!! : Found unexpected targetVolume, not supported for Pixel Barrel modules.");
	}
	moduleMassWithoutSensors_expected += el->quantityInGrams(module);

	if ( el->targetVolume() == PixelModuleHybrid   ||
	     el->targetVolume() == PixelModuleChip        ) {
          vol[el->targetVolume()]->addMaterial(el->elementName(),el->quantityInGrams(module));
          vol[el->targetVolume()]->addMass(el->quantityInGrams(module));
	}
      }

    }

    // OUTER TRACKER MODULE
    if (!module.isPixelModule()) {
      volumes.push_back(vol[HybridFront]);
      volumes.push_back(vol[HybridBack]);
      volumes.push_back(vol[HybridLeft]);
      volumes.push_back(vol[HybridRight]);
      volumes.push_back(vol[HybridBetween]);
      volumes.push_back(vol[SupportPlate]);
    }
    // PIXEL MODULE
    else {
      volumes.push_back(vol[PixelModuleHybrid]);
      volumes.push_back(vol[PixelModuleChip]);
    }

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
    for (const auto& vit : volumes) {
      if (!(vit->getDensity() > 0.)) continue; 
      Composite comp;
      comp.name    = prefix_material + vit->getName();
      comp.density = vit->getDensity();
      comp.method  = wt;

      double m = 0.0;
      for (const auto& it : vit->getMaterialList()) {
	comp.elements.insert(it);
	m += it.second;
      }
   
      for (auto& elem : comp.elements) {
        elem.second /= m;
      }
      
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
}

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
  void Extractor::analyse(MaterialTable& mt, MaterialBudget& mb, XmlTags& trackerXmlTags, CMSSWBundle& trackerData, CMSSWBundle& otstData, bool wt) {

    std::cout << "Starting analysis..." << std::endl;

    Tracker& tr = mb.getTracker();
    InactiveSurfaces& is = mb.getInactiveSurfaces();
    bool isPixelTracker = tr.isPixelTracker();
    bool hasSubDisks = tr.hasSubDisks();

    //std::vector<std::vector<ModuleCap> >& bc = mb.getBarrelModuleCaps();
    std::vector<std::vector<ModuleCap> >& ec = mb.getEndcapModuleCaps();

    // trackerData
    std::vector<Element>& e = trackerData.elements;
    std::vector<Composite>& c = trackerData.composites;
    std::vector<LogicalInfo>& l = trackerData.logic;
    std::vector<ShapeInfo>& s = trackerData.shapes;
    std::vector<ShapeOperationInfo>& so = trackerData.shapeOps;
    std::vector<PosInfo>& p = trackerData.positions;
    std::vector<AlgoInfo>& a = trackerData.algos;
    std::map<std::string,Rotation>& r = trackerData.rots;
    std::vector<SpecParInfo>& t = trackerData.specs;
    std::vector<RILengthInfo>& ri = trackerData.lrilength;

    // otstData
    std::vector<Composite>& otstComposites = otstData.composites;
    std::vector<LogicalInfo>& otstLogic = otstData.logic;
    std::vector<ShapeInfo>& otstShapes = otstData.shapes;
    std::vector<PosInfo>& otstPositions = otstData.positions;

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
  
    Rotation rot;
    // PLACING MODULE WITHIN A ROD
    // Outer Tracker : local X axis, at e.g. Phi = 0°, is parallel to CMS_Y.
    if (!isPixelTracker) {
      // This rotation places an unflipped module within a rod
      rot.name = xml_OT_places_unflipped_mod_in_rod;
      rot.thetax = 90.0;
      rot.phix = 90.0;
      rot.thetay = 0.0;
      rot.phiy = 0.0;
      rot.thetaz = 90.0;
      rot.phiz = 0.0;
      r.insert(std::pair<const std::string,Rotation>(rot.name,rot));

      // This rotation places a flipped module within a rod.
      // Flip is 180° rotation around local X axis.
      rot.name = xml_OT_places_flipped_mod_in_rod;
      rot.thetax = 90.0;
      rot.phix = 270.0;
      rot.thetay = 0.0;
      rot.phiy = 0.0;
      rot.thetaz = 90.0;
      rot.phiz = 180.0;
      r.insert(std::pair<const std::string,Rotation>(rot.name,rot));
    }
    // Inner Tracker : local X axis, at e.g. Phi = 0°, is antiparallel to CMS_Y.
    else {
      // This rotation places an unflipped module within a rod
      rot.name = xml_PX_places_unflipped_mod_in_rod;
      rot.thetax = 90.0;
      rot.phix = 270.0;
      rot.thetay = 180.0;
      rot.phiy = 0.0;
      rot.thetaz = 90.0;
      rot.phiz = 0.0;
      r.insert(std::pair<const std::string,Rotation>(rot.name,rot));

      // This rotation places a flipped module within a rod.
      // Flip is 180° rotation around local X axis.
      rot.name = xml_PX_places_flipped_mod_in_rod;
      rot.thetax = 90.0;
      rot.phix = 90.0;
      rot.thetay = 180.0;
      rot.phiy = 0.0;
      rot.thetaz = 90.0;
      rot.phiz = 180.0;
      r.insert(std::pair<const std::string,Rotation>(rot.name,rot));
    }

    // Y180 : Rotation of 180° around CMS_Y axis.
    rot.name = xml_Y180;
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
    analyseElements(e, c);
    std::cout << "Elementary materials done." << std::endl;
    // Analyse barrel
    analyseLayers(mt, tr, trackerXmlTags, c, l, s, so, p, a, r, t, ri, wt);
    std::cout << "Barrel layers done." << std::endl;
    // Analyse endcaps
    if(!hasSubDisks){
      analyseDiscs(mt, ec, tr, trackerXmlTags, c, l, s, so, p, a, r, t, ri, wt);
    } else {
      analyseDiscsAndSubDiscs(mt, ec, tr, trackerXmlTags, c, l, s, so, p, a, r, t, ri, wt);
    }
    std::cout << "Endcap discs done." << std::endl;
    // Analyse services
    analyseServices(is, isPixelTracker, trackerXmlTags, c, l, s, p, t);
    std::cout << "Services done." << std::endl;
    // Analyse supports
    analyseSupports(is, isPixelTracker, trackerXmlTags, c, l, s, p, otstComposites, otstLogic, otstShapes, otstPositions, t);
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
  void Extractor::analyseElements(std::vector<Element>& elems, std::vector<Composite>& allComposites) {
    
    /*
    // PREVIOUS CODE
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
    */

    const MaterialsTable& myTable = MaterialsTable::instance();

    // FROM TABLE: CHEMICAL ELEMENTS
    const ChemicalElementMap& allChemicalElements = myTable.getAllChemicalElements();

    for (const auto& elemIt : allChemicalElements) {
      const std::string elementName = elemIt.first;
      const ChemicalElement& elem = elemIt.second;

      Element e;
      e.tag = elementName;
      e.density = elem.getDensity() * 1000.;  // g/cm3
      e.atomic_number = elem.getAtomicNumber();
      e.atomic_weight = elem.getAtomicWeight();
      elems.push_back(e);
    }


    // FROM TABLE: MIXTURES (INCLUDE COMPOUNDS)
    const ChemicalMixtureMap& allChemicalMixtures = myTable.getAllChemicalMixtures();
    const int maxMaterialsTreeHierarchyLevel = myTable.getMaxMaterialsTreeHierarchyLevel();

    // This is to print the mixtures in a specific order: from lowest to highest materialsTreeHierarchyLevel.
    // This is because the materials a mixture is made of, need to be printed first (on top of the XML file).
    for (int materialsTreeHierarchyLevel = 1; materialsTreeHierarchyLevel <= maxMaterialsTreeHierarchyLevel; ++materialsTreeHierarchyLevel) {

      for (const auto& mixIt : allChemicalMixtures) {	
	const ChemicalMixture& mix = mixIt.second;

	// Select the mixtures which have hierarchy level of interest.
	if (mix.getMaterialsTreeHierarchyLevel() == materialsTreeHierarchyLevel) {
	  const std::string mixtureName = mixIt.first;

	  Composite comp;
	  comp.name = xml_tkLayout_material + mixtureName;
	  comp.density = mix.getDensity() * 1000.;  // g/cm3
	  comp.method = wt;  // to do: USE ATOMIC FORMULA METHOD FOR COMPOUNDS ?

	  const MassComposition& fractions = mix.getMassComposition();
	  for (const auto& fractionIt : fractions) {
	    comp.elements.insert(std::make_pair(fractionIt.first, fractionIt.second));
	  }

	  comp.isMixture = true;
	  allComposites.push_back(comp); // Add composite to the queue to be printed
	}
      }

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
	  lname << trackerXmlTags.tracker << xml_layer << layer; // e.g. OTLayer1
	  // module name
	  std::ostringstream mname;
	  mname << trackerXmlTags.tracker << lname.str() << xml_R << modRing << xml_barrel_module; // .e.g. OTLayer1R1Bmodule
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

    int first = -999; // Invalid value: it should never be used anyways
    bool hasfirst = false;
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
	  dname << trackerXmlTags.tracker << xml_disc << layer; // e.g. OTDisc6
	  // module name
	  std::ostringstream mname;
	  mname << dname.str() << xml_R << modRing << xml_endcap_module; // e.g. OTDisc6R1EModule
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

      if (!hasfirst)
	logERROR("I am about to use 'first' variable, but I never set it ('hasfirst' is false). This should NOT happen. "  
          "first = " + any2str(first));
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
    SpecParInfo rocdims, lspec, rspec, srspec, trspec, sspec, otcspec, mspec;    
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
    // Module 'stack' when there are two sensors in 3D module
    otcspec.name = trackerXmlTags.topo_bmodulecomb_name + xml_par_tail;
    otcspec.parameter.first = xml_tkddd_structure;
    otcspec.parameter.second = trackerXmlTags.topo_bmodulecomb_value;
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
      
      // is the layer tilted?
      bool isTilted = lagg.getBarrelLayers()->at(layer - 1)->isTilted();

      // is the layer skewed?
      const bool isSkewedLayer = lagg.getBarrelLayers()->at(layer - 1)->isSkewedForInstallation();

      // check that the layer is not both tilted and skewed
      if (isTilted && isSkewedLayer) { logERROR("Layer is both tilted and skewed: this is not supported."); }

      bool isTimingLayer = lagg.getBarrelLayers()->at(layer - 1)->isTiming();
      // TO DO : NOT THAT EXTREMA OF MODULES WITH HYBRIDS HAVE BEEN INTRODUCED INTO DETECTORMODULE (USED TO BE IN MODULE COMPLEX CLASS ONLY)
      // ALL THIS SHOULD BE DELETED AND REPLACED BY module->minRwithHybrids() GETTERS.

      // Calculate geometrical extrema of rod (straight layer), or of rod part + tilted ring (tilted layer)
      // straight layer : x and y extrema of rod
      double xmin = std::numeric_limits<double>::max();
      double xmax = 0;
      double ymin = std::numeric_limits<double>::max();
      double ymax = 0;
      // straight or tilted layer : z and r extrema of rod (straight layer), or of {rod part + tilted ring} (tilted layer)
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
      double firstPhiRodRadius = 0;   // first rod encountered in phi by the visitor.
      double nextPhiRodRadius = 0;  // next rod encountered in phi by the visitor.

      // SKEWED INSTALLATION MODE: COLLECT INFO
      // (-X) side:
      // Non-skewed ladders
      double unskewedLaddersAtMinusXSideCentersMinPhi = 0.;
      double unskewedLaddersAtMinusXSideCentersMaxPhi = 0.;
      bool isPhiComputationAtMinusXSideInitialized = false;
      int countUnskewedLaddersAtMinusXSide = 0;
      // Skewed ladder
      double skewedLadderAtMinusXSideCenterRadius = 0.;
      double skewedLadderAtMinusXSideCenterPhi = 0.;
      double skewedLadderAtMinusXSideSkewAngle = 0.;

      // (+X) side:
      // Non-skewed ladders
      double unskewedLaddersAtPlusXSideCentersMinPhi = 0.;
      double unskewedLaddersAtPlusXSideCentersMaxPhi = 0.;
      bool isPhiComputationAtPlusXSideInitialized = false;
      int countUnskewedLaddersAtPlusXSide = 0;
      // Skewed ladder
      double skewedLadderAtPlusXSideCenterRadius = 0.;
      double skewedLadderAtPlusXSideCenterPhi = 0.;
      double skewedLadderAtPlusXSideSkewAngle = 0.;
      
      // loop on module caps
      for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	const int modRing = iiter->getModule().uniRef().ring;
	// only positive side	
	if (iiter->getModule().uniRef().side > 0 
	    && (
		// Skewed layer: look at all ladders to find radial envelope
		isSkewedLayer     
		// Non-skewed layer: take only modules with uniref phi == 1 or 2
		// WARING: This assumes that there is a 2 rod - periodicity.
		|| (!isSkewedLayer && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2))
		)) {

	  if (isSkewedLayer && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2) && iiter->getModule().isSkewed()) {
	    logERROR("Nothing gonna work, old code assumes that modules placed with uniRef().phi == 1 or 2 are non-skewed.");
	  }
	  
	  // layer name
	  std::ostringstream lname;
	  lname << trackerXmlTags.tracker << xml_layer << layer; // e.g. OTLayer1
	  // module name
	  std::ostringstream mname;
	  mname << lname.str() << xml_R << modRing << xml_barrel_module; //.e.g. OTLayer1R1BModule
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

	  // WARNING: THIS ONLY LOOKS, ON A GIVEN LAYER, AT THE FIRST 2 MODULES. 
	  // IT ASSUMES ALL OTHER PLACEMENTS ARE SIMILAR.
	  // first rod encountered in phi
	  if (iiter->getModule().uniRef().phi == 1 && (modRing == 1 || modRing == 2)) { 
	    firstPhiRodRadius = firstPhiRodRadius + iiter->getModule().center().Rho() / 2; 
	  }
	  // next rod encountered in phi
	  if (iiter->getModule().uniRef().phi == 2 && (modRing == 1 || modRing == 2)) { 
	    nextPhiRodRadius = nextPhiRodRadius + iiter->getModule().center().Rho() / 2; 
	  }

	  // Skewed layer: take relevant info
	  if (isSkewedLayer && modRing == 1) {
	    // Assumes all rings are identical
	    const bool isAtPositiveXSide = iiter->getModule().isAtPlusXSide();

	    // (-X) side
	    if (!isAtPositiveXSide) { 
	      // Non-skewed ladders
	      if (!iiter->getModule().isSkewed()) {
		// initialize
		if (!isPhiComputationAtMinusXSideInitialized) {
		  unskewedLaddersAtMinusXSideCentersMinPhi = iiter->getModule().center().Phi();
		  unskewedLaddersAtMinusXSideCentersMaxPhi = iiter->getModule().center().Phi();
		  isPhiComputationAtMinusXSideInitialized = true;
		}
		else {
		  unskewedLaddersAtMinusXSideCentersMinPhi = moduloMin(iiter->getModule().center().Phi(), unskewedLaddersAtMinusXSideCentersMinPhi, 2.*M_PI);
		  unskewedLaddersAtMinusXSideCentersMaxPhi = moduloMax(iiter->getModule().center().Phi(), unskewedLaddersAtMinusXSideCentersMaxPhi, 2.*M_PI);
		}
		countUnskewedLaddersAtMinusXSide++;
	      }
	      // Skewed ladder
	      else {
		skewedLadderAtMinusXSideCenterRadius = iiter->getModule().center().Rho();
		skewedLadderAtMinusXSideCenterPhi = iiter->getModule().center().Phi();
		skewedLadderAtMinusXSideSkewAngle = iiter->getModule().skewAngle();
	      }
	    }

	    // (+X) side
	    else {
	      // Non-skewed ladders
	      if (!iiter->getModule().isSkewed()) {
		// initialize
		if (!isPhiComputationAtPlusXSideInitialized) {
		  unskewedLaddersAtPlusXSideCentersMinPhi = iiter->getModule().center().Phi();
		  unskewedLaddersAtPlusXSideCentersMaxPhi = iiter->getModule().center().Phi();
		  isPhiComputationAtPlusXSideInitialized = true;
		}
		else {
		  unskewedLaddersAtPlusXSideCentersMinPhi = moduloMin(iiter->getModule().center().Phi(), unskewedLaddersAtPlusXSideCentersMinPhi, 2.*M_PI);
		  unskewedLaddersAtPlusXSideCentersMaxPhi = moduloMax(iiter->getModule().center().Phi(), unskewedLaddersAtPlusXSideCentersMaxPhi, 2.*M_PI);
		}
		countUnskewedLaddersAtPlusXSide++;
	      }
	      // Skewed ladder
	      else {
		skewedLadderAtPlusXSideCenterRadius = iiter->getModule().center().Rho();
		skewedLadderAtPlusXSideCenterPhi = iiter->getModule().center().Phi();
		skewedLadderAtPlusXSideSkewAngle = iiter->getModule().skewAngle();
	      }
	    }

	  } // end: skewed layer info 

	} // end: select a few modules only
      } // end: loop on modules

      if ((rmax - rmin) == 0.0) continue;


      shape.type = bx; // box
      shape.rmin = 0.0;
      shape.rmax = 0.0;


      // for material properties
      double rtotal = 0.0, itotal = 0.0;     
      int count = 0;
      ril.index = layer;
      

      std::ostringstream lname, ladderName, unflippedLadderName, pconverter;
      lname << trackerXmlTags.tracker << xml_layer << layer; // e.g. OTLayer1
      if (!isPixelTracker) ladderName << trackerXmlTags.tracker << xml_layer << layer << xml_rod; // e.g. OTLayer1Rod
      else {
	ladderName << trackerXmlTags.tracker << xml_layer << layer << xml_rod << xml_flipped; // e.g. OTLayer1RodFlipped
	unflippedLadderName << trackerXmlTags.tracker << xml_layer << layer << xml_rod << xml_unflipped; // e.g. OTLayer1RodUnflipped
      }

      std::string places_unflipped_mod_in_rod = (!isPixelTracker ? xml_OT_places_unflipped_mod_in_rod : xml_PX_places_unflipped_mod_in_rod);
      std::string places_flipped_mod_in_rod = (!isPixelTracker ? xml_OT_places_flipped_mod_in_rod : xml_PX_places_flipped_mod_in_rod);

      double firstPhiRodMeanPhi = 0;
      double nextPhiRodMeanPhi = 0;

      std::map<std::tuple<int, int, int, int >, std::string > timingModuleNames;
      bool newTimingModuleType = true;
      int timingModuleCopyNumber = 0;


      bool hasPhiForbiddenRanges = lagg.getBarrelLayers()->at(layer - 1)->phiForbiddenRanges.state();
      std::map<int, double> phiForbiddenRanges;

      // information on tilted rings, indexed by ring number
      std::map<int, BTiltedRingInfo> rinfoplus; // positive-z side
      std::map<int, BTiltedRingInfo> rinfominus; // negative-z side


      // LOOP ON MODULE CAPS 
      for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {

	if (hasPhiForbiddenRanges) {
	  int numRods = lagg.getBarrelLayers()->at(layer - 1)->numRods();
	  int forbiddenPhiUpperAIndex = numRods / 2;
	  int forbiddenPhiLowerBIndex = forbiddenPhiUpperAIndex + 1;

	  std::vector<int> phiForbiddenRangesIndexes;
	  phiForbiddenRangesIndexes.push_back(1);
	  phiForbiddenRangesIndexes.push_back(forbiddenPhiUpperAIndex);
	  phiForbiddenRangesIndexes.push_back(forbiddenPhiLowerBIndex);
	  phiForbiddenRangesIndexes.push_back(numRods);

	  for (const auto& i : phiForbiddenRangesIndexes) {
	    if (iiter->getModule().uniRef().phi == i) {
	      auto found = phiForbiddenRanges.find(i);
	      if (found == phiForbiddenRanges.end()) {
		phiForbiddenRanges.insert(std::make_pair(i, iiter->getModule().center().Phi()));
	      }
	    }
	  }
	}


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
	  std::ostringstream mnameBase;
	  std::ostringstream mnameNeg;
	  if (!iiter->getModule().isTimingModule()) { mname << lname.str() << xml_R << modRing << xml_barrel_module; } // e.g. OTLayer1R1Bmodule
	  // Timing layer : insert both +Z and -Z sides
	  else {
	    mnameBase << lname.str() << xml_R << modRing << xml_barrel_module;
	    mname << mnameBase.str() << xml_positive_z;
	    mnameNeg << mnameBase.str() << xml_negative_z;
	  }


	  // build module volumes, with hybrids taken into account
	  std::string parentName = mname.str();
	  ModuleComplex modcomplex(mname.str(),parentName,*iiter);
	  modcomplex.buildSubVolumes();
#ifdef __DEBUGPRINT__
	  modcomplex.print();
#endif

	  std::vector<ModuleCap>::iterator partner = findPartnerModule(iiter, oiter->end(), modRing);
	  

	  // ROD 1 (STRAIGHT LAYER), OR ROD 1 + MODULES WITH UNIREF PHI == 1 OF THE TILTED RINGS (TILTED LAYER)
	  if (iiter->getModule().uniRef().phi == 1) {

	    firstPhiRodMeanPhi = iiter->getModule().center().Phi();

            std::ostringstream ringname;
	    ringname << lname.str() << xml_ring << modRing;


	    // MODULE

	    // Special case : timing module.
	    // Only define a new module if no module of the same type already defined.
	    // This is for shortening the size of the XML files.
	    std::tuple<int, int, int, int> crystalProperties;
	    std::string childName;
	    if (iiter->getModule().isTimingModule()) {
	      if (iiter->getModule().numSensors() > 0) {
		double crystalWidth = iiter->getModule().sensors().front().crystalWidth();
		double crystalLength = iiter->getModule().sensors().front().crystalLength();
		double crystalThickness = iiter->getModule().sensors().front().crystalThickness();
		double crystalTiltAngle = iiter->getModule().sensors().front().crystalTiltAngle();
		crystalProperties = std::make_tuple(round(crystalWidth * 10.), round(crystalLength * 10.), round(crystalThickness * 10.), round(crystalTiltAngle * 10.));

		// Define new timing module type.
		if (timingModuleNames.count(crystalProperties) == 0) {
		  newTimingModuleType = true;
		  timingModuleNames.insert(std::make_pair(crystalProperties, mnameBase.str()));
		}
		// Uses an already defined timing module with same properties.
		else { newTimingModuleType = false; }
	      }
	      else { std::cout << "ERROR !! Barrel timing module with no sensor." << std::endl; }
	    }

	    // Standard case
	    if (!iiter->getModule().isTimingModule() || newTimingModuleType) {
	      // For SolidSection in tracker.xml : module's box shape
	      shape.name_tag = mname.str();
	      shape.dx = modcomplex.getExpandedModuleWidth()/2.0;
              //if(iiter->getModule().isPixelModule() && iiter->getModule().numSensors()==2){
              //    shape.dy = modcomplex.getExpandedModuleLengthPixelDoubleSens()/2.0;
              //} else {
                  shape.dy = modcomplex.getExpandedModuleLength()/2.0;
              //}
	      shape.dz = modcomplex.getExpandedModuleThickness()/2.0;
	      s.push_back(shape);
	    
	      // For LogicalPartSection in tracker.xml : module's material
	      logic.material_tag = xml_material_air;
	      logic.name_tag = mname.str();
	      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
	      l.push_back(logic);

	      // Timing layer : insert both +Z and -Z sides
	      if (iiter->getModule().isTimingModule()) {
		shape.name_tag = mnameNeg.str();
		s.push_back(shape);
		logic.name_tag = mnameNeg.str();
		logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
		l.push_back(logic);
	      }	    
	    }
            	  

	    // For PosPart section in tracker.xml : module's positions in rod (straight layer) or rod part (tilted layer)
            if (!isTilted || (isTilted && (tiltAngle == 0))) {
	      pos.parent_tag = trackerXmlTags.nspace + ":" + ladderName.str();
	      // DEFINE CHILD : MODULE TO BE PLACED IN A ROD.
	      // Standard case
	      if (!iiter->getModule().isTimingModule()) {
		pos.child_tag = trackerXmlTags.nspace + ":" + mname.str();
	      }
	      // For timing barrel layer : Child, to place in a rod, can be a copy of an existing module.	      
	      else {
		// Define new timing module child.
		if (newTimingModuleType) {
		  childName = trackerXmlTags.nspace + ":" + mnameBase.str();
		  timingModuleCopyNumber = 1;
		}
		// Uses an already defined timing module with same properties.
		else {
		  childName = trackerXmlTags.nspace + ":" + timingModuleNames.at(crystalProperties);
		  timingModuleCopyNumber += 1;
		}
	      }	// end of timing layer special case      
	      pos.trans.dx = iiter->getModule().center().Rho() - firstPhiRodRadius;
	      pos.trans.dz = iiter->getModule().center().Z();

	      if (!iiter->getModule().flipped()) { pos.rotref = trackerXmlTags.nspace + ":" + places_unflipped_mod_in_rod; }
	      else { pos.rotref = trackerXmlTags.nspace + ":" + places_flipped_mod_in_rod; }
	      pos.copy = (!iiter->getModule().isTimingModule() ? 1 : timingModuleCopyNumber);
	      if (iiter->getModule().isTimingModule()) pos.child_tag = childName + xml_positive_z;
	      p.push_back(pos);
	      
	      // This is a copy of the BModule on -Z side
	      if (partner != oiter->end()) {
		pos.trans.dx = partner->getModule().center().Rho() - firstPhiRodRadius;
		pos.trans.dz = partner->getModule().center().Z();

		if (!partner->getModule().flipped()) { pos.rotref = trackerXmlTags.nspace + ":" + places_unflipped_mod_in_rod; }
		else { pos.rotref = trackerXmlTags.nspace + ":" + places_flipped_mod_in_rod; }
		pos.copy = (!iiter->getModule().isTimingModule() ? 2 : timingModuleCopyNumber);
		if (iiter->getModule().isTimingModule()) pos.child_tag = childName + xml_negative_z;
		p.push_back(pos);
	      }
	      // reset
	      pos.copy = 1;
	      pos.rotref = "";
	    }
	  }


	  if (iiter->getModule().uniRef().phi == 2) {
	    if (!isTilted || (isTilted && (tiltAngle == 0))) {
	      nextPhiRodMeanPhi = iiter->getModule().center().Phi();
	      if (isPixelTracker) {
		pos.parent_tag = trackerXmlTags.nspace + ":" + unflippedLadderName.str();
		pos.child_tag = trackerXmlTags.nspace + ":" + mname.str();
		pos.trans.dx = iiter->getModule().center().Rho() - nextPhiRodRadius;
		pos.trans.dz = iiter->getModule().center().Z();
		if (!iiter->getModule().flipped()) { pos.rotref = trackerXmlTags.nspace + ":" + places_unflipped_mod_in_rod; }
		else { pos.rotref = trackerXmlTags.nspace + ":" + places_flipped_mod_in_rod; }
		p.push_back(pos);
	      
		// This is a copy of the BModule on -Z side
		if (partner != oiter->end()) {
		  pos.trans.dx = partner->getModule().center().Rho() - nextPhiRodRadius;
		  pos.trans.dz = partner->getModule().center().Z();
		  if (!partner->getModule().flipped()) { pos.rotref = trackerXmlTags.nspace + ":" + places_unflipped_mod_in_rod; }
		  else { pos.rotref = trackerXmlTags.nspace + ":" + places_flipped_mod_in_rod; }
		  pos.copy = 2; 
		  p.push_back(pos);
		  pos.copy = 1;
		}
		pos.rotref = "";
	      }
	    }
	  }


	  if (iiter->getModule().uniRef().phi == 1) {
	    if (!iiter->getModule().isTimingModule() || newTimingModuleType) {

	      std::ostringstream ringname;
	      ringname << lname.str() << xml_ring << modRing;

	      // Topology
	      sspec.partselectors.push_back(mname.str());
	      sspec.moduletypes.push_back(minfo_zero);
	      if (iiter->getModule().isTimingModule()) { 
		sspec.partselectors.push_back(mnameNeg.str());
		sspec.moduletypes.push_back(minfo_zero);
	      }


             if(iiter->getModule().numSensors()==2){
               otcspec.partselectors.push_back(mname.str());
               otcspec.moduletypes.push_back(minfo_zero);
             }
	      


	      // WAFER
	      string xml_base_lowerupper = "";
	      if (iiter->getModule().numSensors() == 2 && !iiter->getModule().isPixelModule()) {
		xml_base_lowerupper = xml_base_lower;
		shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_waf;
	      }	    
              else if (iiter->getModule().numSensors() == 2 && iiter->getModule().isPixelModule()){
                xml_base_lowerupper = xml_base_one;
                shape.name_tag = mname.str()+xml_InnerPixel+xml_base_lowerupper+xml_base_waf;
              }
	      else {
		if (iiter->getModule().isPixelModule()) shape.name_tag = mname.str() + xml_InnerPixel + xml_base_waf;
		else if (iiter->getModule().isTimingModule()) shape.name_tag = mname.str() + xml_timing + xml_base_waf;
		else { std::cerr << "Wafer : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	      }

	      // SolidSection
	      shape.name_tag = shape.name_tag;
	      shape.dx = iiter->getModule().area() / iiter->getModule().length() / 2.0;
	      shape.dy = iiter->getModule().length() / 2.0;
             //When we have two sensors, need an additional factor 2 to account for the fact that there are two sensors. NB centralDeadAreaLength is the central dead area + the air gap
              if(iiter->getModule().isPixelModule() && iiter->getModule().numSensors()==2){
                    shape.dy = (iiter->getModule().length() - iiter->getModule().centralDeadAreaLength()+iiter->getModule().outerSensorExtraLength())/4.0;
              }

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
	      pos.trans.dx = 0 ;
	      pos.trans.dy = -iiter->getModule().offsetForSensors();
	      pos.trans.dz = /*shape.dz*/ - iiter->getModule().dsDistance() / 2.0; 
	      p.push_back(pos);

	      // Timing layer : insert both +Z and -Z sides
	      if (iiter->getModule().isTimingModule()) {
		shape.name_tag = mnameNeg.str() + xml_timing + xml_base_waf;
		s.push_back(shape);
		logic.name_tag = shape.name_tag;
		logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
		l.push_back(logic);
		pos.parent_tag = trackerXmlTags.nspace + ":" + mnameNeg.str();
		pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
		p.push_back(pos);
	      }

	      if (iiter->getModule().numSensors() == 2 && !iiter->getModule().isPixelModule()) {

		xml_base_lowerupper = xml_base_upper;

		// SolidSection
		shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_waf;
		shape.dy = (iiter->getModule().length() + iiter->getModule().outerSensorExtraLength()) / 2.0;

		s.push_back(shape);

		// LogicalPartSection
		logic.name_tag = shape.name_tag;
		logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
		l.push_back(logic);

		// PosPart section
		pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
	        pos.trans.dy = iiter->getModule().offsetForSensors();
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
	      } else if (iiter->getModule().numSensors()==2 && iiter->getModule().isPixelModule()){
		xml_base_lowerupper = xml_base_two;



		// SolidSection
                shape.name_tag = mname.str()+xml_InnerPixel+xml_base_lowerupper+xml_base_waf;
                //When we have two sensors, need an additional factor 2 to account for the fact that there are two sensors. NB centralDeadAreaLength is the central dead area + the air gap
		shape.dy = (iiter->getModule().length() - iiter->getModule().centralDeadAreaLength() + iiter->getModule().outerSensorExtraLength()) / 4.0; 
		s.push_back(shape);

		// LogicalPartSection
		logic.name_tag = shape.name_tag;
		logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
		l.push_back(logic);

		// PosPart section
		pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
	        pos.trans.dy = iiter->getModule().offsetForSensors();
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
	      if ( !iiter->getModule().isTimingModule()) {
		if (iiter->getModule().numSensors() == 2) xml_base_lowerupper = iiter->getModule().isPixelModule() ? xml_base_one :  xml_base_lower;
		if (iiter->getModule().moduleType() == "ptPS") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_ps + xml_base_pixel + xml_base_act;
		else if (iiter->getModule().moduleType() == "pt2S") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_2s+ xml_base_act;
		else if (iiter->getModule().isTimingModule()) shape.name_tag = mname.str() + xml_timing + xml_base_act;
		else if (iiter->getModule().isPixelModule()) {
		  if (!iiter->getModule().is3DPixelModule()) { shape.name_tag = mname.str() + xml_InnerPixel + xml_base_lowerupper + xml_base_Act; }
		  else { shape.name_tag = mname.str() + xml_InnerPixel + xml_3D + xml_base_Act + xml_base_lowerupper; }
		}
		else { std::cerr << "Active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	    
		// SolidSection
		shape.dx = iiter->getModule().area() / iiter->getModule().length() / 2.0;
		shape.dy = iiter->getModule().length() / 2.0;
                //When we have two sensors, need an additional factor 2 to account for the fact that there are two sensors. NB centralDeadAreaLength is the central dead area + the air gap
                if(iiter->getModule().isPixelModule() && iiter->getModule().numSensors()==2){
                    shape.dy = (iiter->getModule().length() - iiter->getModule().centralDeadAreaLength() + iiter->getModule().outerSensorExtraLength())/4.0;
                }
		shape.dz = iiter->getModule().sensorThickness() / 2.0;
		s.push_back(shape);   

		// LogicalPartSection
		logic.name_tag = shape.name_tag;
		logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
		logic.material_tag = xml_fileident + ":" + xml_tkLayout_material + (iiter->getModule().isTimingModule() && iiter->getModule().subdet() == BARREL ? xml_sensor_LYSO : xml_sensor_silicon);
		l.push_back(logic);

		// PosPart section
		if (iiter->getModule().numSensors() == 2 && !iiter->getModule().isPixelModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_base_lowerupper + xml_base_waf;
                else if (iiter->getModule().numSensors() == 2 && iiter->getModule().isPixelModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() +xml_InnerPixel+ xml_base_lowerupper + xml_base_waf;
		else {
		  if (iiter->getModule().isTimingModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_timing + xml_base_waf;
		  else if (iiter->getModule().isPixelModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_InnerPixel + xml_base_waf;
		  else { std::cerr << "Positioning active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
		}
		pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
                pos.trans.dy = 0.0;
		pos.trans.dz = 0.0;
#ifdef __FLIPSENSORS_IN__ // Flip INNER sensors
		pos.rotref = trackerXmlTags.nspace + ":" + rot_sensor_tag;
#endif
		p.push_back(pos);

		// Topology
		minfo.name		= iiter->getModule().moduleType();
		mspec.partselectors.push_back(shape.name_tag);
                minfo.bricked   = iiter->getModule().innerSensor().isBricked();
		minfo.rocrows	= any2str<int>(iiter->getModule().innerSensor().numROCRows());  // in case of single sensor module innerSensor() and outerSensor() point to the same sensor
		minfo.roccols	= any2str<int>(iiter->getModule().innerSensor().numROCCols());
                if(iiter->getModule().isPixelModule() && iiter->getModule().numSensors()==2){
                    minfo.roccols = any2str<int>(iiter->getModule().innerSensor().numROCCols()/2-2); //Hack for double-sensor modules (pix 3D). This could be made more general
                }
		minfo.rocx		= any2str<int>(iiter->getModule().innerSensor().numROCX());
		minfo.rocy		= any2str<int>(iiter->getModule().innerSensor().numROCY());
		mspec.moduletypes.push_back(minfo);

		if (iiter->getModule().numSensors() == 2) { 

		  xml_base_lowerupper = iiter->getModule().isPixelModule()? xml_base_two : xml_base_upper;

		  // SolidSection
		  if (iiter->getModule().moduleType() == "ptPS") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_ps + xml_base_strip + xml_base_act;
		  else if (iiter->getModule().moduleType() == "pt2S") shape.name_tag = mname.str() + xml_base_lowerupper + xml_base_2s+ xml_base_act;
                  else if (iiter->getModule().isPixelModule()){
                      if (!iiter->getModule().is3DPixelModule()) { shape.name_tag = mname.str() + xml_InnerPixel + xml_base_lowerupper + xml_base_Act; }
                      else { shape.name_tag = mname.str() + xml_InnerPixel + xml_3D + xml_base_Act + xml_base_lowerupper; }
                  }
		  else { std::cerr << "Active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
		  shape.dy = (iiter->getModule().length() + iiter->getModule().outerSensorExtraLength()) / 2.0;
                  //When we have two sensors, need an additional factor 2 to account for the fact that there are two sensors. NB centralDeadAreaLength is the central dead area + the air gap
                  if(iiter->getModule().isPixelModule() && iiter->getModule().numSensors()==2){
                      shape.dy = (iiter->getModule().length() - iiter->getModule().centralDeadAreaLength() + iiter->getModule().outerSensorExtraLength())/4.0;
                  }
		  s.push_back(shape);

		  // LogicalPartSection
		  logic.name_tag = shape.name_tag;
		  logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
		  l.push_back(logic);

		  // PosPart section
                  if (iiter->getModule().isPixelModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() +xml_InnerPixel+ xml_base_lowerupper + xml_base_waf;
                  else pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_base_lowerupper + xml_base_waf;
		  pos.child_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
#ifdef __FLIPSENSORS_OUT__ // Flip OUTER sensors
		  pos.rotref = trackerXmlTags.nspace + ":" + rot_sensor_tag;
#endif
		  p.push_back(pos);

		  // Topology
		  mspec.partselectors.push_back(shape.name_tag);
                  minfo.bricked = iiter->getModule().outerSensor().isBricked();
		  minfo.rocrows	= any2str<int>(iiter->getModule().outerSensor().numROCRows());
		  minfo.roccols	= any2str<int>(iiter->getModule().outerSensor().numROCCols());
                  if(iiter->getModule().isPixelModule() && iiter->getModule().numSensors()==2){
                      minfo.roccols = any2str<int>(iiter->getModule().outerSensor().numROCCols()/2-2); //Hack for double-sensor pixel modules (3D). This could be made more general
                  }
		  minfo.rocx	= any2str<int>(iiter->getModule().outerSensor().numROCX());
		  minfo.rocy	= any2str<int>(iiter->getModule().outerSensor().numROCY());
		  mspec.moduletypes.push_back(minfo);
		} // End of replica for Pt-modules
	      }

	      // Barrel timing modules : add crystals volumes within the wafer volume
	      else {
		if (iiter->getModule().numSensors() > 0) {

		  int numCrystalsX = iiter->getModule().sensors().front().numCrystalsX();
		  int numCrystalsY = iiter->getModule().sensors().front().numCrystalsY();

		  double alveolaWidth = iiter->getModule().sensors().front().alveolaWidth();
		  double alveolaLength = iiter->getModule().sensors().front().alveolaLength();
		  double alveolaShift = iiter->getModule().sensors().front().alveolaShift();

		  std::string crystalName = mnameBase.str() + xml_timing + xml_base_act;

		  double crystalWidth = iiter->getModule().sensors().front().crystalWidth();
		  double crystalLength = iiter->getModule().sensors().front().crystalLength();
		  double crystalThickness = iiter->getModule().sensors().front().crystalThickness();
		  double crystalTiltAngle = iiter->getModule().sensors().front().crystalTiltAngle();
		  		 		 
		  // SolidSection		  
		  shape.name_tag = crystalName;
		  shape.dx = crystalWidth / 2.0;
		  shape.dy = crystalLength / 2.0;
		  shape.dz = crystalThickness / 2.0;
		  s.push_back(shape);

		  // LogicalPartSection
		  logic.name_tag = shape.name_tag;
		  logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
		  logic.material_tag = xml_fileident + ":" + xml_tkLayout_material + xml_sensor_LYSO;
		  l.push_back(logic);

		  // LOOP ON ALL CRYSTALS IN THE MODULE
		  for (int j = 0; j < numCrystalsY; j++) {
		    for (int i = 0; i < numCrystalsX; i++) {		  	         
	    
		      // PosPart section		 
		      pos.child_tag = trackerXmlTags.nspace + ":" + crystalName;
		      int crystalIndex = j * numCrystalsX + i + 1;
		      pos.copy = crystalIndex;
		      double midX = numCrystalsX / 2 - 0.5;
		      pos.trans.dx = (i - midX) * alveolaWidth;
		      double midY = numCrystalsY / 2 - 0.5;
		      pos.trans.dy = (j - midY) * alveolaLength;
		      pos.trans.dz = pow(-1., i + j) * alveolaShift;

		      addTiltedModuleRot(r, crystalTiltAngle);
		      std::ostringstream tilt;
		      tilt << round(crystalTiltAngle);

		      pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_timing + xml_base_waf;
		      pos.rotref = trackerXmlTags.nspace + ":" + xml_positive_z_tilted_mod_rot + tilt.str();
		      p.push_back(pos);

		      // -Z side
		      pos.parent_tag = trackerXmlTags.nspace + ":" + mnameNeg.str() + xml_timing + xml_base_waf;
		      pos.rotref = trackerXmlTags.nspace + ":" + xml_negative_z_tilted_mod_rot + tilt.str();
		      p.push_back(pos);

		      // reset
		      pos.copy = 1;
		      pos.trans.dx = 0.;
		      pos.trans.dy = 0.;
		      pos.trans.dz = 0.;
		      pos.rotref = "";
		    }
		  } // loop on crystals

		  // Topology
		  if (std::find(mspec.partselectors.begin(), mspec.partselectors.end(), crystalName) == mspec.partselectors.end()) {
		    minfo.name		= iiter->getModule().moduleType();
		    mspec.partselectors.push_back(crystalName);
		    minfo.bricked	= iiter->getModule().innerSensor().isBricked();
		    minfo.rocrows	= any2str<int>(iiter->getModule().innerSensor().numROCRows());
		    minfo.roccols	= any2str<int>(iiter->getModule().innerSensor().numROCCols());
		    minfo.rocx		= any2str<int>(iiter->getModule().innerSensor().numROCX());
		    minfo.rocy		= any2str<int>(iiter->getModule().innerSensor().numROCY());
		    mspec.moduletypes.push_back(minfo);
		  }
		}
		else { std::cout << "ERROR !! Barrel timing module with no sensor." << std::endl; }
	      }

	      modcomplex.addMaterialInfo(c);
	      modcomplex.addShapeInfo(s);
	      modcomplex.addLogicInfo(l);
	      modcomplex.addPositionInfo(p);
	      // -Z side
	      if ( iiter->getModule().isTimingModule()) {
		std::string parentNameNeg = mnameNeg.str();
		ModuleComplex modcomplexNeg(mnameNeg.str(),parentNameNeg,*partner);
		modcomplexNeg.buildSubVolumes();
		modcomplexNeg.addMaterialInfo(c);
		modcomplexNeg.addShapeInfo(s);
		modcomplexNeg.addLogicInfo(l);
		modcomplexNeg.addPositionInfo(p);
	      }
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
		rinf.rmax = modcomplex.getRmax();
		rinf.zmax = modcomplex.getZmax();		
		rinf.rminatzmin = modcomplex.getRminatZmin();
		rinf.rmaxatzmax = modcomplex.getRmaxatZmax();      
		rinfoplus.insert(std::pair<int, BTiltedRingInfo>(modRing, rinf));

		// same ring on negative-z side
		rinf.name = ringname.str() + xml_minus;
		rinf.isZPlus = 0;
		rinf.z1 = - iiter->getModule().center().Z(); // WARNING: this assumes symmetry through (XY) plane!
		rinfominus.insert(std::pair<int, BTiltedRingInfo>(modRing, rinf));
	      }
	    }


	    // material properties
	    rtotal = rtotal + iiter->getRadiationLength();
	    itotal = itotal + iiter->getInteractionLength();
	    count++;
	    //double dt = modcomplex.getExpandedModuleThickness();
	  }	  


	  // ONLY MODULES WITH UNIREF PHI == 2 OF THE TILTED RINGS (TILTED LAYER)
	  if (isTilted && (iiter->getModule().uniRef().phi == 2)) {
	    // Fill the info of the z-positive ring with matching ring number
	    const auto& foundPlusZTiltedRingInfo = rinfoplus.find(modRing);
	    if (foundPlusZTiltedRingInfo != rinfoplus.end()) {
	      BTiltedRingInfo& myInfo = foundPlusZTiltedRingInfo->second;
	      myInfo.fw_flipped = iiter->getModule().flipped();
	      myInfo.r2 = iiter->getModule().center().Rho();
	      myInfo.z2 = iiter->getModule().center().Z();
	      myInfo.startPhiAngle2 = iiter->getModule().center().Phi();
	      // Update tilted ring extrema
	      myInfo.rmin = MIN(myInfo.rmin, modcomplex.getRmin());
	      myInfo.rmax = MAX(myInfo.rmax, modcomplex.getRmax());
	      myInfo.zmin = MIN(myInfo.zmin, modcomplex.getZmin());
	      myInfo.zmax = MAX(myInfo.zmax, modcomplex.getZmax());
	      myInfo.rminatzmin = MIN(myInfo.rminatzmin, modcomplex.getRminatZmin()); 
	      myInfo.rmaxatzmax = MAX(myInfo.rmaxatzmax, modcomplex.getRmaxatZmax());
	    }
	    // Fill the info of the z-negative ring with matching ring number
	    const auto& foundMinusZTiltedRingInfo = rinfominus.find(modRing);
	    if (foundMinusZTiltedRingInfo != rinfominus.end()) {
	      BTiltedRingInfo& myInfo = foundMinusZTiltedRingInfo->second;
	      myInfo.fw_flipped = iiter->getModule().flipped();
	      myInfo.r2 = iiter->getModule().center().Rho();
	      myInfo.z2 = - iiter->getModule().center().Z(); // WARNING: this assumes symmetry through (XY) plane!
	      myInfo.startPhiAngle2 = iiter->getModule().center().Phi();
	      // Update tilted ring extrema
	      myInfo.rmin = MIN(myInfo.rmin, modcomplex.getRmin());
	      myInfo.rmax = MAX(myInfo.rmax, modcomplex.getRmax());
	      myInfo.zmin = MIN(myInfo.zmin, modcomplex.getZmin());
	      myInfo.zmax = MAX(myInfo.zmax, modcomplex.getZmax());
	      myInfo.rminatzmin = MIN(myInfo.rminatzmin, modcomplex.getRminatZmin()); 
	      myInfo.rmaxatzmax = MAX(myInfo.rmaxatzmax, modcomplex.getRmaxatZmax());
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
      shape.name_tag = ladderName.str();
      shape.dx = (ymax - ymin) / 2 + xml_epsilon;
      if (isTilted && !isPixelTracker) shape.name_tag = ladderName.str() + "Full";
      if (isTilted) shape.dx = (flatPartMaxY - flatPartMinY) / 2 + xml_epsilon;
      if (isPixelTracker || isTimingLayer) shape.dx = rodThickness.at(layer) + xml_epsilon;
      shape.dy = (xmax - xmin) / 2 + xml_epsilon;
      if (isTilted) shape.dy = (flatPartMaxX - flatPartMinX) / 2 + xml_epsilon;
      if (isPixelTracker || isTimingLayer) shape.dy = rodWidth.at(layer) + xml_epsilon;
      shape.dz = zmax + xml_epsilon;
      if (isTilted) shape.dz = flatPartMaxZ + xml_epsilon;
      s.push_back(shape);

      // Subtraction of an air volume from the flat part rod container volume, to avoid collision with first tilted ring
      // This trick is only used for the Outer Tracker.
      // For the Inner Tracker, this trick doesn't make sense, since in priciple smallDelta = 0.
      if (isTilted && !isPixelTracker && flatPartNumModules >= 2) {
	shape.name_tag = ladderName.str() + "Air";
	shape.dx = shape.dx / 2.0;
	shape.dy = shape.dy + xml_epsilon;
	shape.dz = (shape.dz - flatPartOneBeforeLastModuleMaxZ) / 2.;
	s.push_back(shape);
	
	shapeOp.name_tag = ladderName.str() + "SubtractionIntermediate";
	shapeOp.type = substract;
	shapeOp.rSolid1 = ladderName.str() + "Full";
	shapeOp.rSolid2 = ladderName.str() + "Air";
	shapeOp.trans.dx = shape.dx + xml_epsilon;
	shapeOp.trans.dy = 0.;
	shapeOp.trans.dz = flatPartMaxZ + xml_epsilon - shape.dz + xml_epsilon;
	so.push_back(shapeOp);

	shapeOp.name_tag = ladderName.str();
	shapeOp.type = substract;
	shapeOp.rSolid1 = ladderName.str() + "SubtractionIntermediate";
	shapeOp.rSolid2 = ladderName.str() + "Air";
	shapeOp.trans.dx = shape.dx + xml_epsilon;
	shapeOp.trans.dy = 0.;
	shapeOp.trans.dz = -shapeOp.trans.dz;
	so.push_back(shapeOp);
      }

      logic.name_tag = ladderName.str();
      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
      logic.material_tag = xml_material_air;
      l.push_back(logic);

      if (isPixelTracker) {
	shape.name_tag = unflippedLadderName.str();
	s.push_back(shape);
	logic.name_tag = unflippedLadderName.str();
	logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
	l.push_back(logic);
	rspec.partselectors.push_back(unflippedLadderName.str());
	srspec.partselectors.push_back(unflippedLadderName.str());
      }

      rspec.partselectors.push_back(ladderName.str());
      rspec.moduletypes.push_back(minfo_zero);
      srspec.partselectors.push_back(ladderName.str());
      srspec.moduletypes.push_back(minfo_zero);



      // rods in layer algorithm(s)

      // Get inner / outer ladders radii
      const double innerLadderCenterRadius = MIN(firstPhiRodRadius, nextPhiRodRadius);  // inner radius
      const double outerLadderCenterRadius = MAX(firstPhiRodRadius, nextPhiRodRadius); // outer radius

      // OUTER TRACKER
      if (!isPixelTracker) {
	if (!hasPhiForbiddenRanges) {
	  alg.name = xml_phialt_algo;
	  alg.parent = trackerXmlTags.nspace + ":" + lname.str();
	  pconverter <<  trackerXmlTags.nspace + ":" + ladderName.str();
	  alg.parameters.push_back(stringParam(xml_childparam, pconverter.str()));
	  alg.parameters.push_back(numericParam(xml_tilt, "90*deg")); // This "tilt" here has nothing to do with the tilt angle of a tilted TBPS.
	  // It is an angle used internally by PhiAltAlgo to shift in Phi the startAngle ( in (X,Y) plane).
	  // 90 deg corresponds to no shift. SHOULD NOT BE MODIFIED!!	  
	  // Is firstPhiRod placed at inner radius?
	  const bool isFirstPhiRodAtInnerRadius = (fabs(firstPhiRodRadius - innerLadderCenterRadius) < xml_epsilon);
	  // The algo (as implemeted in CMSSW) starts by placing the inner radius rod, no matter what!
	  // So need to start placing rods from the inner rod mean phi.
	  const double algoStartPhi = (isFirstPhiRodAtInnerRadius ? firstPhiRodMeanPhi : nextPhiRodMeanPhi); 
	  pconverter.str("");
	  pconverter << algoStartPhi * 180. / M_PI << "*deg"; 
	  alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
	  pconverter.str("");
	  alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
	  pconverter << innerLadderCenterRadius << "*mm";
	  alg.parameters.push_back(numericParam(xml_radiusin, pconverter.str()));
	  pconverter.str("");
	  pconverter << outerLadderCenterRadius << "*mm";
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
	  int numRods = lagg.getBarrelLayers()->at(layer - 1)->numRods();
	  int forbiddenPhiUpperAIndex = numRods / 2;
	  int forbiddenPhiLowerBIndex = forbiddenPhiUpperAIndex + 1;

	  alg.name = xml_phialt_algo;
	  alg.parent = trackerXmlTags.nspace + ":" + lname.str();
	  pconverter <<  trackerXmlTags.nspace + ":" + ladderName.str();
	  alg.parameters.push_back(stringParam(xml_childparam, pconverter.str()));
	  alg.parameters.push_back(numericParam(xml_tilt, "90*deg")); // This "tilt" here has nothing to do with the tilt angle of a tilted TBPS.
	  // It is an angle used internally by PhiAltAlgo to shift in Phi the startAngle ( in (X,Y) plane).
	  // 90 deg corresponds to no shift. SHOULD NOT BE MODIFIED!!
	  pconverter.str("");
	  pconverter << phiForbiddenRanges.at(1) * 180. / M_PI << "*deg";
	  alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
	  pconverter.str("");
	  pconverter << (phiForbiddenRanges.at(forbiddenPhiUpperAIndex) - phiForbiddenRanges.at(1)) * 180. / M_PI << "*deg";
	  alg.parameters.push_back(numericParam(xml_rangeangle, pconverter.str()));
	  // WARNING: Set RadisuIn parameter to the radius of the firstPhiRod.
	  // The algorithm will indeed start by placing (in phi) that firstPhiRod, at radius whatever is assigned to radisuIn.
	  // RadiusIn is just a name, and does not imply that the radius is low !!!!
	  // Look at PhiAltAlgo implementation.
	  pconverter.str("");
	  pconverter << firstPhiRodRadius << "*mm";
	  alg.parameters.push_back(numericParam(xml_radiusin, pconverter.str()));
	  pconverter.str("");
	  pconverter << nextPhiRodRadius << "*mm";
	  alg.parameters.push_back(numericParam(xml_radiusout, pconverter.str()));
	  pconverter.str("");
	  alg.parameters.push_back(numericParam(xml_zposition, "0.0*mm"));
	  pconverter << numRods / 2;
	  alg.parameters.push_back(numericParam(xml_number, pconverter.str()));
	  alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
	  alg.parameters.push_back(numericParam(xml_incrcopyno, "1"));
	  a.push_back(alg);
	  alg.parameters.clear();

	  alg.name = xml_phialt_algo;
	  alg.parent = trackerXmlTags.nspace + ":" + lname.str();
	  pconverter.str("");
	  pconverter <<  trackerXmlTags.nspace + ":" + ladderName.str();
	  alg.parameters.push_back(stringParam(xml_childparam, pconverter.str()));
	  alg.parameters.push_back(numericParam(xml_tilt, "90*deg")); // This "tilt" here has nothing to do with the tilt angle of a tilted TBPS.
	  // It is an angle used internally by PhiAltAlgo to shift in Phi the startAngle ( in (X,Y) plane).
	  // 90 deg corresponds to no shift. SHOULD NOT BE MODIFIED!!
	  pconverter.str("");
	  pconverter << phiForbiddenRanges.at(forbiddenPhiLowerBIndex) * 180. / M_PI << "*deg";
	  alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
	  pconverter.str("");
	  pconverter << (phiForbiddenRanges.at(numRods) - phiForbiddenRanges.at(forbiddenPhiLowerBIndex)) * 180. / M_PI << "*deg";
	  alg.parameters.push_back(numericParam(xml_rangeangle, pconverter.str()));
	  pconverter.str("");
	  pconverter << firstPhiRodRadius << "*mm";
	  alg.parameters.push_back(numericParam(xml_radiusin, pconverter.str()));
	  pconverter.str("");
	  pconverter << nextPhiRodRadius << "*mm";
	  alg.parameters.push_back(numericParam(xml_radiusout, pconverter.str()));
	  pconverter.str("");
	  alg.parameters.push_back(numericParam(xml_zposition, "0.0*mm"));
	  pconverter << numRods / 2;
	  alg.parameters.push_back(numericParam(xml_number, pconverter.str()));
	  pconverter.str("");
	  pconverter << forbiddenPhiLowerBIndex;
	  alg.parameters.push_back(numericParam(xml_startcopyno, pconverter.str()));
	  alg.parameters.push_back(numericParam(xml_incrcopyno, "1"));
	  a.push_back(alg);
	  pconverter.str("");
	  alg.parameters.clear();
	}
      }

      // INNER TRACKER
      else {

	// COMMON ALGORITHM PARAMETERS
	const std::string nameSpace = trackerXmlTags.nspace;
	const std::string parentName = lname.str();
	const XYZVector center = XYZVector(0., 0., 0.);
	const int copyNumberIncrement = 2;


	// NON-SKEWED LAYER
	if (!isSkewedLayer) {
	  
	  // COMMON ALGORITHM PARAMETERS
	  const double rangeAngle = 2. * M_PI;	  
	  const int numberLadders = lagg.getBarrelLayers()->at(layer - 1)->numRods() / 2;	 

	  // FIRST PHI LADDERS ALGORITHM
	  const std::string firstPhiLadderName = ladderName.str();
	  const double firstPhiLadderCenterPhi = firstPhiRodMeanPhi;
	  const double firstPhiLadderRadius = firstPhiRodRadius;
	  const int firstPhiLadderStartCopyNumber = 1;

	  createAndStoreDDTrackerAngularAlgorithmBlock(a,
						       nameSpace, 
						       parentName,
						       firstPhiLadderName,
						       firstPhiLadderCenterPhi,
						       rangeAngle,
						       firstPhiLadderRadius,
						       center,
						       numberLadders,
						       firstPhiLadderStartCopyNumber,
						       copyNumberIncrement);

	  // NEXT PHI LADDERS ALGORITHM
	  const std::string nextPhiLadderName = unflippedLadderName.str();
	  const double nextPhiLadderCenterPhi = nextPhiRodMeanPhi;
	  const double nextPhiLadderRadius = nextPhiRodRadius;
	  const int nextPhiLadderStartCopyNumber = 2;

	  createAndStoreDDTrackerAngularAlgorithmBlock(a,
						       nameSpace, 
						       parentName,
						       nextPhiLadderName,
						       nextPhiLadderCenterPhi,						    
						       rangeAngle,
						       nextPhiLadderRadius,
						       center,
						       numberLadders,
						       nextPhiLadderStartCopyNumber,
						       copyNumberIncrement);
	}

	// SKEWED LAYER
	// WARNING: This assumes that, per (X) side:
	// - there is only 1 skewed ladder (placed at outer radius).
	// - there is an odd number of non-skewed ladders (>=3).

	// ALGO BLOCKS TO PLACE UNSKEWED LADDERS
	else {

	  // (-X) SIDE
	  if (countUnskewedLaddersAtMinusXSide <= 2 || femod(countUnskewedLaddersAtMinusXSide, 2) == 0) { 
	    logERROR("Skewed layer: It is expected to have an odd number (>=3) of non-skewed ladders per (X) side."); 
	  }
	  else {

	    // COMMON ALGORITHM PARAMETERS
	    const double rangeAngleAtMinusXSide = femod(
							femod(unskewedLaddersAtMinusXSideCentersMaxPhi, 2. * M_PI) 
							- femod(unskewedLaddersAtMinusXSideCentersMinPhi, 2. * M_PI)
							, 2. * M_PI);

	    // FLIPPED LADDERS BLOCK (INNER RADIUS)
	    // flipped ladders
	    const std::string flippedLadderName = ladderName.str();
	    // phi
	    const double flippedLadderCenterPhi = unskewedLaddersAtMinusXSideCentersMinPhi; // WARNING: THIS ASSUMES
	                                                                                    // that the non-skewed ladder placed at min phi 
	                                                                                    // is placed at inner radius and flipped!!!
	    // number of flipped ladders
	    const int countFlippedLaddersAtMinusXSide = std::floor(countUnskewedLaddersAtMinusXSide / 2.) + 1; // WARNING: THIS ASSUMES
	                                                                                                       // that there is an odd number of 
                                                                                                               // non-skewed ladders: 
	                                                                                                       // 1 more ladder at inner radius.
	    // start copy number
	    const int flippedLadderStartCopyNumber = 1;

	    // algo block
	    createAndStoreDDTrackerAngularAlgorithmBlock(a,
							 nameSpace, 
							 parentName,
							 flippedLadderName,
							 flippedLadderCenterPhi,
							 rangeAngleAtMinusXSide,
							 innerLadderCenterRadius,
							 center,
							 countFlippedLaddersAtMinusXSide,
							 flippedLadderStartCopyNumber,
							 copyNumberIncrement);	  

	    // UNFLIPPED LADDERS BLOCK (OUTER RADIUS)
	    // non-flipped ladders
	    const std::string unFlippedLadderName = unflippedLadderName.str();
	    // delta phi between 2 consecutive non-skewed ladders 
	    const double deltaPhiAtMinusXSide = rangeAngleAtMinusXSide / (countUnskewedLaddersAtMinusXSide - 1);
	    // phi
	    const double unFlippedLadderCenterPhi = unskewedLaddersAtMinusXSideCentersMinPhi + deltaPhiAtMinusXSide;
	    // number of non-flipped ladders
	    const int countUnFlippedLaddersAtMinusXSide = std::floor(countUnskewedLaddersAtMinusXSide / 2.); // WARNING: THIS ASSUMES
	                                                                                                     // that there is an odd number of 
                                                                                                             // non-skewed ladders: 
	                                                                                                     // 1 less ladder at outer radius.
	    // start copy number
	    const int unFlippedLadderStartCopyNumber = 2;

	    // algo block
	    createAndStoreDDTrackerAngularAlgorithmBlock(a,
							 nameSpace, 
							 parentName,
							 unFlippedLadderName,
							 unFlippedLadderCenterPhi,
							 rangeAngleAtMinusXSide - 2. * deltaPhiAtMinusXSide,
							 outerLadderCenterRadius,
							 center,
							 countUnFlippedLaddersAtMinusXSide,
							 unFlippedLadderStartCopyNumber,
							 copyNumberIncrement);
	  }
	  

	  // (+X) side
	  if (countUnskewedLaddersAtMinusXSide <= 2 || femod(countUnskewedLaddersAtMinusXSide, 2) == 0) { 
	    logERROR("Skewed layer: It is expected to have an odd number (>=3) of non-skewed ladders per (X) side."); 
	  }
	  else {

	    // COMMON ALGORITHM PARAMETERS 
	    const double rangeAngleAtPlusXSide = femod(
						       femod(unskewedLaddersAtPlusXSideCentersMaxPhi, 2. * M_PI) 
						       - femod(unskewedLaddersAtPlusXSideCentersMinPhi, 2. * M_PI)
						       , 2.*M_PI);

	    // FLIPPED LADDERS BLOCK (INNER RADIUS)
	    // flipped ladders
	    const std::string flippedLadderName = ladderName.str();
	    // phi	  
	    const double flippedLadderCenterPhi = unskewedLaddersAtPlusXSideCentersMinPhi; // WARNING: THIS ASSUMES
	                                                                                    // that the non-skewed ladder placed at min phi 
	                                                                                    // is placed at inner radius and flipped!!!
	    // number of flipped ladders
	    const int countFlippedLaddersAtPlusXSide = std::floor(countUnskewedLaddersAtPlusXSide / 2.) + 1; // WARNING: THIS ASSUMES
	                                                                                                     // that there is an odd number of 
                                                                                                             // non-skewed ladders: 
	                                                                                                     // 1 more ladder at inner radius.
	    // start copy number
	    const int flippedLadderStartCopyNumber = countUnskewedLaddersAtMinusXSide + 1;

	    // algo block
	    createAndStoreDDTrackerAngularAlgorithmBlock(a,
							 nameSpace, 
							 parentName,
							 flippedLadderName,
							 flippedLadderCenterPhi,
							 rangeAngleAtPlusXSide,
							 innerLadderCenterRadius,
							 center,
							 countFlippedLaddersAtPlusXSide,
							 flippedLadderStartCopyNumber,
							 copyNumberIncrement);	  

	    // UNFLIPPED LADDERS BLOCK (OUTER RADIUS)
	    // non-flipped ladders
	    const std::string unFlippedLadderName = unflippedLadderName.str();
	    // delta phi between 2 consecutive non-skewed ladders
	    const double deltaPhiAtPlusXSide = rangeAngleAtPlusXSide / (countUnskewedLaddersAtPlusXSide - 1);
	    // phi
	    const double unFlippedLadderCenterPhi = unskewedLaddersAtPlusXSideCentersMinPhi + deltaPhiAtPlusXSide;
	    // number of non-flipped ladders
	    const int countUnFlippedLaddersAtPlusXSide = std::floor(countUnskewedLaddersAtPlusXSide / 2.); // WARNING: THIS ASSUMES
	                                                                                                   // that there is an odd number of 
                                                                                                           // non-skewed ladders: 
	                                                                                                   // 1 less ladder at outer radius.
	    // start copy number
	    const int unFlippedLadderStartCopyNumber = countUnskewedLaddersAtMinusXSide + 2;

	    // algo block
	    createAndStoreDDTrackerAngularAlgorithmBlock(a,
							 nameSpace, 
							 parentName,
							 unFlippedLadderName,
							 unFlippedLadderCenterPhi,
							 rangeAngleAtPlusXSide - 2. * deltaPhiAtPlusXSide,
							 outerLadderCenterRadius,
							 center,
							 countUnFlippedLaddersAtPlusXSide,
							 unFlippedLadderStartCopyNumber,
							 copyNumberIncrement);
	  }


	  const int countTotalUnskewedLadders = countUnskewedLaddersAtMinusXSide + countUnskewedLaddersAtPlusXSide;

	  // SKEWED LADDER, (-X) side
	  pos.parent_tag = nameSpace + ":" + parentName;
	  pos.child_tag = nameSpace  + ":" + unflippedLadderName.str(); // WARNING: THIS ASSUMES
	  // that the skewed ladder is always non-flipped.

	  pos.trans.dx = skewedLadderAtMinusXSideCenterRadius * cos(skewedLadderAtMinusXSideCenterPhi);
	  pos.trans.dy = skewedLadderAtMinusXSideCenterRadius * sin(skewedLadderAtMinusXSideCenterPhi);
	  pos.trans.dz = 0;

	  const double skewRotationAngleInRadAtMinusXSide = skewedLadderAtMinusXSideCenterPhi + skewedLadderAtMinusXSideSkewAngle;	    
	  const std::string skewRotationNameAtMinusXSide = "Z" + any2str(skewRotationAngleInRadAtMinusXSide * 180. / M_PI, xml_angle_name_precision);
	  addRotationAroundZAxis(r, skewRotationNameAtMinusXSide, skewRotationAngleInRadAtMinusXSide);
	  
	  pos.rotref = nameSpace + ":" + skewRotationNameAtMinusXSide;
	  pos.copy = countTotalUnskewedLadders + 1;
	  p.push_back(pos);


	  // SKEWED LADDER, (+X) side
	  pos.parent_tag = nameSpace + ":" + parentName;
	  pos.child_tag = nameSpace  + ":" + unflippedLadderName.str(); // WARNING: THIS ASSUMES
	  // that the skewed ladder is always non-flipped.

	  pos.trans.dx = skewedLadderAtPlusXSideCenterRadius * cos(skewedLadderAtPlusXSideCenterPhi);
	  pos.trans.dy = skewedLadderAtPlusXSideCenterRadius * sin(skewedLadderAtPlusXSideCenterPhi);
	  pos.trans.dz = 0;

	  const double skewRotationAngleInRadAtPlusXSide = skewedLadderAtPlusXSideCenterPhi + skewedLadderAtPlusXSideSkewAngle;
	  const std::string skewRotationNameAtPlusXSide = "Z" + any2str(skewRotationAngleInRadAtPlusXSide * 180. / M_PI, xml_angle_name_precision);
	  addRotationAroundZAxis(r, skewRotationNameAtPlusXSide, skewRotationAngleInRadAtPlusXSide);
	  
	  pos.rotref = nameSpace + ":" + skewRotationNameAtPlusXSide;
	  pos.copy = countTotalUnskewedLadders + 2;
	  p.push_back(pos);

	  // reset
	  pos.trans.dx = 0;
	  pos.trans.dy = 0;
	  pos.trans.dz = 0;
	  pos.rotref = "";
	  pos.copy = 1;
	} // end of skewed layer

      }

      // reset
      shape.dx = 0.0;
      shape.dy = 0.0;
      shape.dyy = 0.0;	  
      pos.trans.dx = 0;
      pos.trans.dy = 0;
      pos.trans.dz = 0;
      pos.rotref = "";
      pos.copy = 1;

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
	      
	      // Tilted ring: first part to be stored
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
	      pconverter << rinfo.r1 << "*mm";
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
	      
	      // Tilted ring: second part to be stored
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
	      pconverter << rinfo.r2 << "*mm";
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
      lspec.moduletypes.push_back(minfo_zero);

      layer++;
    } // end: loop on layers
    if (!lspec.partselectors.empty()) t.push_back(lspec);
    if (!rspec.partselectors.empty()) t.push_back(rspec); 
    if (!srspec.partselectors.empty()) t.push_back(srspec);
    if (!trspec.partselectors.empty()) t.push_back(trspec);
    if (!sspec.partselectors.empty()) t.push_back(sspec);
    if (!otcspec.partselectors.empty()) t.push_back(otcspec);
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
	std::map<int,double> ringSmallAbsZModulesZMax, ringBigAbsZModulesZMin;
	for (const auto& i : ringsIndexes) {
	  ringzmin.insert( {i, std::numeric_limits<double>::max()} );
	  ringzmax.insert( {i, 0.} );
	  ringSmallAbsZModulesZMax.insert( {i, 0.} );
	  ringBigAbsZModulesZMin.insert( {i, std::numeric_limits<double>::max()} );
	}

	// loop on module caps
	for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	  if (iiter->getModule().uniRef().side > 0 && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)) {
	    int modRing = iiter->getModule().uniRef().ring;
	    //disk name
	    std::ostringstream dname;
	    dname << trackerXmlTags.tracker << xml_disc << discNumber; // e.g. OTDisc6
	    // module name
	    std::ostringstream mname;
	    mname << dname.str() << xml_R << modRing << xml_endcap_module; // e.g. OTDisc6R1EModule
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
	    if (iiter->getModule().isSmallerAbsZModuleInRing()) { 
	      ringSmallAbsZModulesZMax.at(modRing) = MAX(ringSmallAbsZModulesZMax.at(modRing), modcomplex.getZmax());
	    }
	    else {
	      ringBigAbsZModulesZMin.at(modRing) = MIN(ringBigAbsZModulesZMin.at(modRing), modcomplex.getZmin());
	    }
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


	//if (zmin > 0) 	
        std::ostringstream dname, pconverter;
	//disk name
        dname << trackerXmlTags.tracker << xml_disc << discNumber; // e.g. OTDisc6

        std::map<int, ERingInfo> rinfo;
	std::set<int> ridx;

      
        std::map<int, std::vector<double>> phi_one;
        std::map<int, std::vector<double>> phi_two;
        std::map<int, std::vector<double>> radius_one;
        std::map<int, std::vector<double>> radius_two;
        std::map<int, std::vector<double>> yaw_one;
        std::map<int, std::vector<double>> yaw_two;

        for (iiter = oiter->begin(); iiter != oiter->end(); iiter++){
          if(!iiter->getModule().inRegularRing()){ //Don't need to check this for every endcap ring, only modules that are NOT in a regular ring
            if(phi_one.count(iiter->getModule().uniRef().ring)==0){
              phi_one[iiter->getModule().uniRef().ring] = std::vector<double>();
            } 
            if(radius_one.count(iiter->getModule().uniRef().ring)==0){
              radius_one[iiter->getModule().uniRef().ring] = std::vector<double>();
            } 
            if(yaw_one.count(iiter->getModule().uniRef().ring)==0){
              yaw_one[iiter->getModule().uniRef().ring] = std::vector<double>();
            } 
            if(phi_two.count(iiter->getModule().uniRef().ring)==0){
              phi_two[iiter->getModule().uniRef().ring] = std::vector<double>();
            }
            if(radius_two.count(iiter->getModule().uniRef().ring)==0){
              radius_two[iiter->getModule().uniRef().ring] = std::vector<double>();
            }
            if(yaw_two.count(iiter->getModule().uniRef().ring)==0){
              yaw_two[iiter->getModule().uniRef().ring] = std::vector<double>();
            }
            if(iiter->getModule().uniRef().phi%2 == 1){
              phi_one[iiter->getModule().uniRef().ring].push_back(iiter->getModule().center().Phi()*180./M_PI);
              radius_one[iiter->getModule().uniRef().ring].push_back(iiter->getModule().center().Rho());
              yaw_one[iiter->getModule().uniRef().ring].push_back(iiter->getModule().yawAngle()*180./M_PI);
            } else {
              phi_two[iiter->getModule().uniRef().ring].push_back(iiter->getModule().center().Phi()*180./M_PI);
              radius_two[iiter->getModule().uniRef().ring].push_back(iiter->getModule().center().Rho());
              yaw_two[iiter->getModule().uniRef().ring].push_back(iiter->getModule().yawAngle()*180./M_PI);
            }
          }
        }

        // LOOP ON MODULE CAPS
        for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	  if (iiter->getModule().uniRef().side > 0 &&  (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)) {
	    // ring number
	    int modRing = iiter->getModule().uniRef().ring;

	    //std::cout << "iiter->getModule().uniRef().phi = " << iiter->getModule().uniRef().phi << " iiter->getModule().uniRef().ring = " << iiter->getModule().uniRef().ring << "iiter->getModule().center().Phi() = " << iiter->getModule().center().Phi() * 180. / M_PI << " iiter->getModule().center().Rho() = " << iiter->getModule().center().Rho() << " iiter->getModule().center().X() = " << iiter->getModule().center().X() << " iiter->getModule().center().Y() = " << iiter->getModule().center().Y() << " iiter->getModule().center().Z() = " << iiter->getModule().center().Z() << " iiter->getModule().flipped() = " << iiter->getModule().flipped() << " iiter->getModule().moduleType() = " << iiter->getModule().moduleType() << std::endl;


	    if (iiter->getModule().uniRef().phi == 1) {
	      // new ring
	      //if (ridx.find(modRing) == ridx.end()) 
	      ridx.insert(modRing);

	      std::ostringstream matname, rname, mname, specname;
	      // ring name
	      rname << dname.str() << xml_ring << modRing; // e.g. OTDisc6Ring1
	      // module name
	      mname << dname.str() << xml_R << modRing << xml_endcap_module; // e.g. OTDisc6R1EModule
 
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
	      if (iiter->getModule().isPixelModule()) shape.name_tag = mname.str() + xml_InnerPixel + xml_base_waf;
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
		shape.dy = (iiter->getModule().length() + iiter->getModule().outerSensorExtraLength()) / 2.0;
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
	      else if (iiter->getModule().isPixelModule()) {
		if (!iiter->getModule().is3DPixelModule()) { shape.name_tag = mname.str() + xml_InnerPixel + xml_base_Act; }
		else { shape.name_tag = mname.str() + xml_InnerPixel + xml_3D + xml_base_Act + xml_base_lowerupper; }
	      }
	      else { std::cerr << "Active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	      shape.dy = iiter->getModule().length() / 2.0;
	      s.push_back(shape);

	      logic.name_tag = shape.name_tag;
	      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
	      logic.material_tag = xml_fileident + ":" + xml_tkLayout_material + xml_sensor_silicon;
	      l.push_back(logic);

	      if (iiter->getModule().numSensors() == 2) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_base_lowerupper + xml_base_waf;
	      else {
		if (iiter->getModule().isTimingModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_timing + xml_base_waf;
		else if (iiter->getModule().isPixelModule())  pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_InnerPixel + xml_base_waf;
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
	      minfo.bricked	= iiter->getModule().innerSensor().isBricked();
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
		shape.dy = (iiter->getModule().length() + iiter->getModule().outerSensorExtraLength()) / 2.0;
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
		minfo.bricked	= iiter->getModule().outerSensor().isBricked();
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
	      
	      // Collect ring info
	      const Ring* myRing = lagg.getEndcapLayers()->at(layer - 1)->ringsMap().at(modRing);
	      ERingInfo myRingInfo;
	      myRingInfo.name = rname.str();
	      myRingInfo.childname = mname.str();
	      myRingInfo.isDiskAtPlusZEnd = iiter->getModule().uniRef().side;
	      myRingInfo.numModules = myRing->numModules();
	      myRingInfo.moduleThickness = modcomplex.getExpandedModuleThickness();
	      myRingInfo.radiusMin  = modcomplex.getRmin();
	      myRingInfo.radiusMid = iiter->getModule().center().Rho();
	      myRingInfo.radiusMax = modcomplex.getRmax();
	      myRingInfo.zMin = ringzmin.at(modRing);
	      myRingInfo.smallAbsZSurfaceZMax = ringSmallAbsZModulesZMax.at(modRing);
	      myRingInfo.bigAbsZSurfaceZMin = ringBigAbsZModulesZMin.at(modRing);
	      myRingInfo.zMax = ringzmax.at(modRing);
	      myRingInfo.zMid = (myRingInfo.zMin + myRingInfo.zMax) / 2.;
	      myRingInfo.isRingOn4Dees = myRing->isRingOn4Dees();
              myRingInfo.isRegularRing = iiter->getModule().inRegularRing();

	      // surface 1 is whatever surface the (phi == 1) module belongs to.
	      myRingInfo.surface1ZMid = iiter->getModule().center().Z();
	      myRingInfo.surface1StartPhi = iiter->getModule().center().Phi();
	      myRingInfo.surface1IsFlipped = iiter->getModule().flipped();
	      rinfo.insert(std::pair<int, ERingInfo>(modRing, myRingInfo));

	      // material properties
	      rtotal = rtotal + iiter->getRadiationLength();
	      itotal = itotal + iiter->getInteractionLength();
	      count++;
	    }

	    if (iiter->getModule().uniRef().phi == 2) {
	      std::map<int,ERingInfo>::iterator it;
	      // Fill the info of the other surface of the same ring.
	      // This assumes that the (phi == 1) module is on another surface.
	      // Hence that the surfaces the modules are assigned to are alternated, as one goes along Phi.
	      it = rinfo.find(modRing);
	      if (it != rinfo.end()) {
		it->second.surface2ZMid = iiter->getModule().center().Z();
		it->second.surface2StartPhi = iiter->getModule().center().Phi();
		it->second.surface2IsFlipped = iiter->getModule().flipped();
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

        for (const auto& ringIndex : ridx) {

	  const auto& found = rinfo.find(ringIndex);
	  if (found != rinfo.end()) {
	    const ERingInfo& myRingInfo = found->second;
	    if (myRingInfo.numModules > 0) {

	      // CASE WHERE SMALLDELTA > BIGDELTA: MODULES OF A GIVEN RING ARE SPREAD ALONG 4 DEES.
	      if (myRingInfo.isRingOn4Dees) {

		const auto& lookForInnerRing = rinfo.find(ringIndex - 1);
		const bool hasAnInnerRing = lookForInnerRing != rinfo.end();
		const auto& lookForOuterRing = rinfo.find(ringIndex + 1);
		const bool hasAnOuterRing = lookForOuterRing != rinfo.end();

		// Full ring cylinder
		shape.name_tag = ((hasAnInnerRing || hasAnOuterRing) ? (myRingInfo.name + "Full") : myRingInfo.name);
		shape.rmin = myRingInfo.radiusMin - xml_epsilon;
		shape.rmax = myRingInfo.radiusMax + xml_epsilon;
		shape.dz = (myRingInfo.zMax - myRingInfo.zMin) / 2.0 + xml_epsilon;
		s.push_back(shape);

		// Is there another ring with smaller radius?
		// If so, a section needs to be removed to avoid clashes.
		if (hasAnInnerRing) {
		  ERingInfo& myInnerRingInfo = lookForInnerRing->second;
		  const double myInnerRingRMax = myInnerRingInfo.radiusMax;

		  shape.name_tag = myRingInfo.name + "InnerCut";
		  shape.rmin = myRingInfo.radiusMin - 2. * xml_epsilon;
		  shape.rmax = myInnerRingRMax + 2. * xml_epsilon;
		  shape.dz = (myRingInfo.bigAbsZSurfaceZMin - myRingInfo.smallAbsZSurfaceZMax) / 2.0 - xml_epsilon;
		  s.push_back(shape);

		  shapeOp.name_tag = (hasAnOuterRing ? (myRingInfo.name + "SubtractionIntermediate") : myRingInfo.name);
		  shapeOp.type = substract;
		  shapeOp.rSolid1 = myRingInfo.name + "Full";
		  shapeOp.rSolid2 = myRingInfo.name + "InnerCut";
		  shapeOp.trans.dx = 0.;
		  shapeOp.trans.dy = 0.;
		  shapeOp.trans.dz = 0.;
		  so.push_back(shapeOp);
		}
		
		// Is there another ring with bigger radius?
		// If so, a section needs to be removed to avoid clashes.
		if (hasAnOuterRing) {
		  ERingInfo& myOuterRingInfo = lookForOuterRing->second;
		   const double myOuterRingRMin = myOuterRingInfo.radiusMin;

		  shape.name_tag = myRingInfo.name + "OuterCut";
		  shape.rmin = myOuterRingRMin - 2. * xml_epsilon;
		  shape.rmax = myRingInfo.radiusMax + 2. * xml_epsilon;
		  shape.dz = (myRingInfo.bigAbsZSurfaceZMin - myRingInfo.smallAbsZSurfaceZMax) / 2.0 - xml_epsilon;
		  s.push_back(shape);

		  shapeOp.name_tag = myRingInfo.name;
		  shapeOp.type = substract;
		  shapeOp.rSolid1 = (hasAnInnerRing ? (myRingInfo.name + "SubtractionIntermediate") : (myRingInfo.name + "Full"));
		  shapeOp.rSolid2 = myRingInfo.name + "OuterCut";
		  shapeOp.trans.dx = 0.;
		  shapeOp.trans.dy = 0.;
		  shapeOp.trans.dz = 0.;
		  so.push_back(shapeOp);
		}	  
	      }

	      // CASE WHERE SMALLDELTA < BIGDELTA: MODULES OF A GIVEN RING ARE PLACED ON BOTH SIDES OF THE SAME DEE.
	      // There are only 2 dees per ring: one per (X) side.
	      else {
		shape.name_tag = myRingInfo.name;
		shape.rmin = myRingInfo.radiusMin - xml_epsilon;
		shape.rmax = myRingInfo.radiusMax + xml_epsilon;
		shape.dz = (myRingInfo.zMax - myRingInfo.zMin) / 2.0 + xml_epsilon;
		s.push_back(shape);
	      }
	    }
	      

	    logic.name_tag = myRingInfo.name;
            logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
            logic.material_tag = xml_material_air;
            l.push_back(logic);

            pos.parent_tag = trackerXmlTags.nspace + ":" + dname.str(); // CUIDADO ended with: + xml_plus;
            pos.child_tag = logic.shape_tag;

	    pos.trans.dz = myRingInfo.zMid - diskZ;
            p.push_back(pos);

            rspec.partselectors.push_back(logic.name_tag);
            rspec.moduletypes.push_back(minfo_zero);

	    // Ring Surface 1 (half the modules)
	    alg.name = xml_trackerring_algo;
            if(!myRingInfo.isRegularRing) alg.name=xml_trackerring_irregular_algo;
            alg.parent = logic.shape_tag;
            alg.parameters.push_back(stringParam(xml_childparam, trackerXmlTags.nspace + ":" + myRingInfo.childname));
            pconverter << (myRingInfo.numModules / 2);
            alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
            pconverter.str("");
            alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
            alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
            alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
            pconverter << myRingInfo.surface1StartPhi * 180. / M_PI << "*deg";
            alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
            pconverter.str("");
            pconverter << myRingInfo.radiusMid << "*mm";
            alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
            pconverter.str("");
            if(!myRingInfo.isRegularRing && (phi_one[ringIndex]).size()>0){
              alg.parameters.push_back(arbitraryLengthVector("phiAngleValues",phi_one[ringIndex]));
              pconverter.str("");
              alg.parameters.push_back(arbitraryLengthVector("yawAngleValues",yaw_one[ringIndex]));
              pconverter.str("");
              alg.parameters.push_back(arbitraryLengthVector("radiusValues",radius_one[ringIndex]));
              pconverter.str("");
            } 
	    alg.parameters.push_back(vectorParam(0, 0, myRingInfo.surface1ZMid - myRingInfo.zMid));
	    pconverter << myRingInfo.isDiskAtPlusZEnd;
	    alg.parameters.push_back(numericParam(xml_iszplus, pconverter.str()));
	    pconverter.str("");
	    alg.parameters.push_back(numericParam(xml_tiltangle, "90*deg"));
	    pconverter << myRingInfo.surface1IsFlipped;
	    alg.parameters.push_back(numericParam(xml_isflipped, pconverter.str()));
	    pconverter.str("");
            a.push_back(alg);
            alg.parameters.clear();

	    // Ring Surface 2 (half the modules)
	    alg.name = xml_trackerring_algo;
            if(!myRingInfo.isRegularRing) alg.name=xml_trackerring_irregular_algo;
            alg.parameters.push_back(stringParam(xml_childparam, trackerXmlTags.nspace + ":" + myRingInfo.childname));
            pconverter << (myRingInfo.numModules / 2);
            alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
            pconverter.str("");
            alg.parameters.push_back(numericParam(xml_startcopyno, "2"));
            alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
            alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
            pconverter << myRingInfo.surface2StartPhi * 180. / M_PI << "*deg";
            alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
            pconverter.str("");
            pconverter << myRingInfo.radiusMid << "*mm";
            alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
            pconverter.str("");
            if(!myRingInfo.isRegularRing && (phi_two[ringIndex]).size()>0 ){
              alg.parameters.push_back(arbitraryLengthVector("phiAngleValues",phi_two[ringIndex]));
              pconverter.str("");
              alg.parameters.push_back(arbitraryLengthVector("yawAngleValues",yaw_two[ringIndex]));
              pconverter.str("");
              alg.parameters.push_back(arbitraryLengthVector("radiusValues",radius_two[ringIndex]));
              pconverter.str("");
            } 
	    alg.parameters.push_back(vectorParam(0, 0, myRingInfo.surface2ZMid - myRingInfo.zMid));
	    pconverter << myRingInfo.isDiskAtPlusZEnd;
	    alg.parameters.push_back(numericParam(xml_iszplus, pconverter.str()));
	    pconverter.str("");
	    alg.parameters.push_back(numericParam(xml_tiltangle, "90*deg"));
	    pconverter << myRingInfo.surface2IsFlipped;
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
   * This is a modification of analyseDiscs, for the case where we have subdisks. It examines the endcap discs in z+
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
  void Extractor::analyseDiscsAndSubDiscs(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& ec, Tracker& tr, XmlTags& trackerXmlTags,
                               std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<ShapeOperationInfo>& so,
			       std::vector<PosInfo>& p, std::vector<AlgoInfo>& a, std::map<std::string,Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, bool wt) {

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
    SpecParInfo rocdims, dspec, sdspec, rspec, sspec, mspec;
    // Disc
    dspec.name = trackerXmlTags.topo_disc_name + xml_par_tail;
    dspec.parameter.first = xml_tkddd_structure;
    dspec.parameter.second = trackerXmlTags.topo_disc_value;
    // SubDisc
    sdspec.name = trackerXmlTags.topo_subdisc_name + xml_par_tail;
    sdspec.parameter.first = xml_tkddd_structure;
    sdspec.parameter.second = trackerXmlTags.topo_subdisc_value;
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
	std::map<int,double> ringSmallAbsZModulesZMax, ringBigAbsZModulesZMin;
	for (const auto& i : ringsIndexes) {
	  ringzmin.insert( {i, std::numeric_limits<double>::max()} );
	  ringzmax.insert( {i, 0.} );
	  ringSmallAbsZModulesZMax.insert( {i, 0.} );
	  ringBigAbsZModulesZMin.insert( {i, std::numeric_limits<double>::max()} );
	}

	// loop on module caps
	for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	  if (iiter->getModule().uniRef().side > 0 && (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)) {
	    int modRing = iiter->getModule().uniRef().ring;
	    //disk name
	    std::ostringstream dname;
	    dname << trackerXmlTags.tracker << xml_disc << discNumber; // e.g. OTDisc6
	    //subdisk name
	    std::ostringstream sdname;
            if(iiter->getModule().uniRef().phi==1){
              sdname << dname.str() << xml_subdisc << "1";
            } else {
              sdname << dname.str() << xml_subdisc << "2";
            }
	    // module name
	    std::ostringstream mname;
            if(iiter->getModule().uniRef().phi==1){
              mname << dname.str() << xml_SD << "1"<< xml_R << modRing << xml_endcap_module; // e.g. OTDisc6SD1R1EModule
            } else {
              mname << dname.str() << xml_SD << "2"<< xml_R << modRing << xml_endcap_module; // e.g. OTDisc6SD2R1EModule
            }
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
	    if (iiter->getModule().isSmallerAbsZModuleInRing()) { 
	      ringSmallAbsZModulesZMax.at(modRing) = MAX(ringSmallAbsZModulesZMax.at(modRing), modcomplex.getZmax());
	    }
	    else {
	      ringBigAbsZModulesZMin.at(modRing) = MIN(ringBigAbsZModulesZMin.at(modRing), modcomplex.getZmin());
	    }
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


	//if (zmin > 0) 	
        std::ostringstream dname, pconverter, sdname;
	//disk name
        dname << trackerXmlTags.tracker << xml_disc << discNumber; // e.g. OTDisc6

        std::map<int, ERingInfo> rinfo;
	std::set<int> ridx;
        int nRings = ringsIndexes.size();

      
        std::map<int, std::vector<double>> phi_one;
        std::map<int, std::vector<double>> phi_two;
        std::map<int, std::vector<double>> radius_one;
        std::map<int, std::vector<double>> radius_two;
        std::map<int, std::vector<double>> yaw_one;
        std::map<int, std::vector<double>> yaw_two;

        for (iiter = oiter->begin(); iiter != oiter->end(); iiter++){
          if(!iiter->getModule().inRegularRing()){ //Don't need to check this for every endcap ring, only modules that are NOT in a regular ring
            if(phi_one.count(iiter->getModule().uniRef().ring)==0){
              phi_one[iiter->getModule().uniRef().ring] = std::vector<double>();
            } 
            if(radius_one.count(iiter->getModule().uniRef().ring)==0){
              radius_one[iiter->getModule().uniRef().ring] = std::vector<double>();
            } 
            if(yaw_one.count(iiter->getModule().uniRef().ring)==0){
              yaw_one[iiter->getModule().uniRef().ring] = std::vector<double>();
            } 
            if(phi_two.count(iiter->getModule().uniRef().ring)==0){
              phi_two[iiter->getModule().uniRef().ring] = std::vector<double>();
            }
            if(radius_two.count(iiter->getModule().uniRef().ring)==0){
              radius_two[iiter->getModule().uniRef().ring] = std::vector<double>();
            }
            if(yaw_two.count(iiter->getModule().uniRef().ring)==0){
              yaw_two[iiter->getModule().uniRef().ring] = std::vector<double>();
            }
            if(iiter->getModule().uniRef().phi%2 == 1){
              phi_one[iiter->getModule().uniRef().ring].push_back(iiter->getModule().center().Phi()*180./M_PI);
              radius_one[iiter->getModule().uniRef().ring].push_back(iiter->getModule().center().Rho());
              yaw_one[iiter->getModule().uniRef().ring].push_back(iiter->getModule().yawAngle()*180./M_PI);
            } else {
              phi_two[iiter->getModule().uniRef().ring].push_back(iiter->getModule().center().Phi()*180./M_PI);
              radius_two[iiter->getModule().uniRef().ring].push_back(iiter->getModule().center().Rho());
              yaw_two[iiter->getModule().uniRef().ring].push_back(iiter->getModule().yawAngle()*180./M_PI);
            }
          }
        }

        // LOOP ON MODULE CAPS
        for (iiter = oiter->begin(); iiter != oiter->end(); iiter++) {
	  if (iiter->getModule().uniRef().side > 0 &&  (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi == 2)) {
	    // ring number
	    int modRing = iiter->getModule().uniRef().ring;

	    if (iiter->getModule().uniRef().phi == 1 || iiter->getModule().uniRef().phi==2) {
	      // new ring  - in SD1 or 2 depending on whether it's phi 1/2
	      int phiuniref = iiter->getModule().uniRef().phi;
	      ridx.insert(modRing); //Just for the index, we'll loop over the subdisk indices as well later

	      std::ostringstream matname, sdname, rname, mname, specname;
              if(phiuniref==1){
                // subdisc name
                sdname << dname.str() << xml_subdisc << "1"; // e.g. OTDisc6SubDisc1
                // ring name
                rname << dname.str() << xml_subdisc << "1" << xml_ring << modRing; // e.g. OTDisc6SubDisc1Ring1
                // module name
                mname << dname.str() << xml_SD << "1" << xml_R << modRing << xml_endcap_module; // e.g. OTDisc6SD1R1EModule
              } else {
                // subdisc name
                sdname << dname.str() << xml_subdisc << "2"; // e.g. OTDisc6SubDisc2
                // ring name
                rname << dname.str() << xml_subdisc << "2" << xml_ring << modRing; // e.g. OTDisc6SubDisc2Ring1
                // module name
                mname << dname.str() << xml_SD << "2" << xml_R << modRing << xml_endcap_module; // e.g. OTDisc6SD2R1EModule
              }
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

	      logic.material_tag = xml_material_air;
	      l.push_back(logic);

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
	      if (iiter->getModule().isPixelModule()) shape.name_tag = mname.str() + xml_InnerPixel + xml_base_waf;
	      else if (iiter->getModule().isTimingModule()) shape.name_tag = mname.str() + xml_timing + xml_base_waf;
	      else { std::cerr << "Wafer : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	    }     

	      shape.dz = iiter->getModule().sensorThickness() / 2.0; // CUIDADO WAS calculateSensorThickness(*iiter, mt) / 2.0;
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
		shape.dy = (iiter->getModule().length() + iiter->getModule().outerSensorExtraLength()) / 2.0;
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
	      else if (iiter->getModule().isPixelModule()) {
		if (!iiter->getModule().is3DPixelModule()) { shape.name_tag = mname.str() + xml_InnerPixel + xml_base_Act; }
		else { shape.name_tag = mname.str() + xml_InnerPixel + xml_3D + xml_base_Act + xml_base_lowerupper; }
	      }
	      else { std::cerr << "Active surface : Unknown module type : " << iiter->getModule().moduleType() << "." << std::endl; }
	      shape.dy = iiter->getModule().length() / 2.0;
	      s.push_back(shape);

	      logic.name_tag = shape.name_tag;
	      logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
	      logic.material_tag = xml_fileident + ":" + xml_tkLayout_material + xml_sensor_silicon;
	      l.push_back(logic);

	      if (iiter->getModule().numSensors() == 2) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_base_lowerupper + xml_base_waf;
	      else {
		if (iiter->getModule().isTimingModule()) pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_timing + xml_base_waf;
		else if (iiter->getModule().isPixelModule())  pos.parent_tag = trackerXmlTags.nspace + ":" + mname.str() + xml_InnerPixel + xml_base_waf;
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
		shape.dy = (iiter->getModule().length() + iiter->getModule().outerSensorExtraLength()) / 2.0;
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
	      
	      // Collect ring info
	      const Ring* myRing = lagg.getEndcapLayers()->at(layer - 1)->ringsMap().at(modRing);
	      ERingInfo myRingInfo;
	      myRingInfo.name = rname.str();
	      myRingInfo.childname = mname.str();
	      myRingInfo.isDiskAtPlusZEnd = iiter->getModule().uniRef().side;
	      myRingInfo.numModules = myRing->numModules();
	      myRingInfo.moduleThickness = modcomplex.getExpandedModuleThickness();
	      myRingInfo.radiusMin  = modcomplex.getRmin();
	      myRingInfo.radiusMid = iiter->getModule().center().Rho();
	      myRingInfo.radiusMax = modcomplex.getRmax();
	      myRingInfo.zMin = ringzmin.at(modRing);
	      myRingInfo.smallAbsZSurfaceZMax = ringSmallAbsZModulesZMax.at(modRing);
	      myRingInfo.bigAbsZSurfaceZMin = ringBigAbsZModulesZMin.at(modRing);
	      myRingInfo.zMax = ringzmax.at(modRing);
	      myRingInfo.zMid = (myRingInfo.zMin + myRingInfo.zMax) / 2.;
	      myRingInfo.isRingOn4Dees = myRing->isRingOn4Dees();
              myRingInfo.isRegularRing = iiter->getModule().inRegularRing();

	      // surface 1 is whatever surface the (phi == 1) module belongs to.
	      myRingInfo.surface1ZMid = iiter->getModule().center().Z();
	      myRingInfo.surface1StartPhi = iiter->getModule().center().Phi();
	      myRingInfo.surface1IsFlipped = iiter->getModule().flipped();
	      rinfo.insert(std::pair<int, ERingInfo>(nRings*(phiuniref-1)+modRing, myRingInfo));

	      // material properties
	      rtotal = rtotal + iiter->getRadiationLength();
	      itotal = itotal + iiter->getInteractionLength();
	      count++;
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
        
       

        std::set<int> sdidx; //Double disc, so if we have subdisks there must be 2
        sdidx.insert(1);
        sdidx.insert(2);
        
        for (const auto& sdIndex : sdidx) {
            std::ostringstream sdname;
            sdname << dname.str() << xml_subdisc << sdIndex; // e.g. OTDisc6SubDisc2
            //First determine subdisk position, based on the first 2 rings
            float subDiskZPosition=0;
            float ringZ=0;
            for (const auto& ringIndex : ridx) {
                if(ringIndex<3){
                  const auto& found = rinfo.find(nRings*(sdIndex-1)+ringIndex);
                  if (found!=rinfo.end()) {
                      const ERingInfo& myRingInfo = found->second;
                      subDiskZPosition+=myRingInfo.surface1ZMid;
                      ringZ+=myRingInfo.zMid;
                  }   
                }
              }
              subDiskZPosition=subDiskZPosition/2.;
              ringZ=ringZ/2.;
                
            for (const auto& ringIndex : ridx) {
            const auto& found = rinfo.find(nRings*(sdIndex-1)+ringIndex);
            if (found != rinfo.end()) {
	      const ERingInfo& myRingInfo = found->second;
              if (myRingInfo.numModules > 0) {

                //When we have subdisks all rings are "flat" and there are no sections to be removed
                shape.name_tag = myRingInfo.name;
                shape.rmin = myRingInfo.radiusMin - xml_epsilon;
                shape.rmax = myRingInfo.radiusMax + xml_epsilon;
                shape.dz = diskThickness/(2*sdidx.size()*ridx.size()); //Bit of a hack, but works for now.
                s.push_back(shape);
	      }
	      

	      logic.name_tag = myRingInfo.name;
              logic.shape_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
              logic.material_tag = xml_material_air;
              l.push_back(logic);



              pos.parent_tag = trackerXmlTags.nspace + ":" + sdname.str(); // CUIDADO ended with: + xml_plus;
              pos.child_tag = logic.shape_tag;

              pos.trans.dz = myRingInfo.surface1ZMid-subDiskZPosition;
              p.push_back(pos);

              rspec.partselectors.push_back(logic.name_tag);
              rspec.moduletypes.push_back(minfo_zero);

	      // Ring Surface 1 (half the modules - subdisk 1)
	      if(sdIndex==1){
                alg.name = xml_trackerring_algo;
                if(!myRingInfo.isRegularRing) alg.name=xml_trackerring_irregular_algo;
                alg.parent = logic.shape_tag;
                alg.parameters.push_back(stringParam(xml_childparam, trackerXmlTags.nspace + ":" + myRingInfo.childname));
                pconverter << (myRingInfo.numModules / 2);
                alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
                pconverter.str("");
                alg.parameters.push_back(numericParam(xml_startcopyno, "1"));
                alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
                alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
                pconverter << myRingInfo.surface1StartPhi * 180. / M_PI << "*deg";
                alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
                pconverter.str("");
                pconverter << myRingInfo.radiusMid << "*mm";
                alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
                pconverter.str("");
                if(!myRingInfo.isRegularRing && (phi_one[ringIndex]).size()>0){
                  alg.parameters.push_back(arbitraryLengthVector("phiAngleValues",phi_one[ringIndex]));
                  pconverter.str("");
                  alg.parameters.push_back(arbitraryLengthVector("yawAngleValues",yaw_one[ringIndex]));
                  pconverter.str("");
                  alg.parameters.push_back(arbitraryLengthVector("radiusValues",radius_one[ringIndex]));
                  pconverter.str("");
                } 
	        alg.parameters.push_back(vectorParam(0, 0, 0));
	        pconverter << myRingInfo.isDiskAtPlusZEnd;
	        alg.parameters.push_back(numericParam(xml_iszplus, pconverter.str()));
	        pconverter.str("");
	        alg.parameters.push_back(numericParam(xml_tiltangle, "90*deg"));
	        pconverter << myRingInfo.surface1IsFlipped;
                alg.parameters.push_back(numericParam(xml_isflipped, pconverter.str()));
                pconverter.str("");
                a.push_back(alg);
                alg.parameters.clear();
              }

	      // Ring Surface 2 (half the modules - subdisk 2)
	      if(sdIndex==2){
                alg.name = xml_trackerring_algo;
                if(!myRingInfo.isRegularRing) alg.name=xml_trackerring_irregular_algo;
                alg.parent = logic.shape_tag;
                alg.parameters.push_back(stringParam(xml_childparam, trackerXmlTags.nspace + ":" + myRingInfo.childname));
                pconverter << (myRingInfo.numModules / 2);
                alg.parameters.push_back(numericParam(xml_nmods, pconverter.str()));
                pconverter.str("");
                alg.parameters.push_back(numericParam(xml_startcopyno, "2"));
                alg.parameters.push_back(numericParam(xml_incrcopyno, "2"));
                alg.parameters.push_back(numericParam(xml_rangeangle, "360*deg"));
                pconverter << myRingInfo.surface1StartPhi * 180. / M_PI << "*deg";
                alg.parameters.push_back(numericParam(xml_startangle, pconverter.str()));
                pconverter.str("");
                pconverter << myRingInfo.radiusMid << "*mm";
                alg.parameters.push_back(numericParam(xml_radius, pconverter.str()));
                pconverter.str("");
                if(!myRingInfo.isRegularRing && (phi_two[ringIndex]).size()>0 ){
                  alg.parameters.push_back(arbitraryLengthVector("phiAngleValues",phi_two[ringIndex]));
                  pconverter.str("");
                  alg.parameters.push_back(arbitraryLengthVector("yawAngleValues",yaw_two[ringIndex]));
                  pconverter.str("");
                  alg.parameters.push_back(arbitraryLengthVector("radiusValues",radius_two[ringIndex]));
                  pconverter.str("");
                } 
	        alg.parameters.push_back(vectorParam(0, 0, 0));
	        pconverter << myRingInfo.isDiskAtPlusZEnd;
	        alg.parameters.push_back(numericParam(xml_iszplus, pconverter.str()));
	        pconverter.str("");
	        alg.parameters.push_back(numericParam(xml_tiltangle, "90*deg"));
	        pconverter << myRingInfo.surface1IsFlipped;
	        alg.parameters.push_back(numericParam(xml_isflipped, pconverter.str()));
	        pconverter.str("");
                a.push_back(alg);
                alg.parameters.clear();
              }
            }
          }
          //subdisc
          shape.name_tag = sdname.str();
          shape.rmin = rmin - 2 * xml_epsilon;
          shape.rmax = rmax + 2 * xml_epsilon;
          shape.dz = (diskThickness / (2.0*sdidx.size())) + 2 * xml_epsilon; //(zmax - zmin) / 2.0; This is a hack, but it works for now
          s.push_back(shape);

          shape.name_tag = sdname.str();
          logic.name_tag = shape.name_tag; // CUIDADO ended with + xml_plus;
          logic.shape_tag = trackerXmlTags.nspace + ":" + shape.name_tag;
          logic.material_tag = xml_material_air;
          l.push_back(logic);

          pos.parent_tag = trackerXmlTags.nspace + ":" + dname.str();
          pos.child_tag = trackerXmlTags.nspace + ":" + logic.name_tag;
          pos.trans.dz = subDiskZPosition - ringZ;
          p.push_back(pos);

          sdspec.partselectors.push_back(logic.name_tag);
          sdspec.moduletypes.push_back(minfo_zero);
          sdspec.partextras.push_back(logic.extra);

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
	pos.trans.dz = diskZ - xml_z_pixfwd;
        p.push_back(pos);

        dspec.partselectors.push_back(logic.name_tag);
        dspec.moduletypes.push_back(minfo_zero);
        dspec.partextras.push_back(logic.extra);
	discNumber++;
      }
      layer++;
    }
    if (!dspec.partselectors.empty()) t.push_back(dspec);
    if (!sdspec.partselectors.empty()) t.push_back(sdspec);
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
	  if (!isPixelTracker) startEndcaps = xml_outerTrackerEndcapsMinZ;
	  else startEndcaps = xml_innerTrackerEndcapsMinZ;  // Tilted Inner Tracker : startEndcaps = xml_innerTiltedTrackerEndcapsMinZ;
          
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
	    pos.rotref = trackerXmlTags.nspace + ":" + xml_Y180;
	    p.push_back(pos);
	  }

	  // ENDCAPS services
	  else {
	    // VERY IMPORTANT : Barrel timing Layer : Services from Tracker endcaps removed. 
	    // TEMPORARY !! This is because the Timing Layer is put in the tracker, and hence tracker and timing services are grouped togetehr.
	    //if (!isTimingLayout || (isTimingLayout && (iter->getInnerRadius() + iter->getRWidth()) <= 1160.)) { // BTL

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
		pos.rotref = trackerXmlTags.nspace + ":" + xml_Y180;
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
	      //}
	      //else { std::cout << "VERY IMPORTANT : Barrel timing Layer : Services from Tracker endcaps removed. TEMPORARY !! This is because the Timing Layer is put in the tracker, and hence tracker and timing services are grouped togetehr. TO DO : Place timing layer in an independant container." << std::endl; }
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
				  std::vector<ShapeInfo>& s, std::vector<PosInfo>& p,
				  std::vector<Composite>& otstComposites, std::vector<LogicalInfo>& otstLogic, 
				  std::vector<ShapeInfo>& otstShapes, std::vector<PosInfo>& otstPositions,
				  std::vector<SpecParInfo>& t, 
				  bool wt) {
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

   
      // Test whether the support is the OTST.
      const std::map<LocalElement, double, ElementNameCompare>& allMasses = iter->getLocalElementsDetails();
      bool isOTST = false;
      for (const auto& massIt : allMasses) {
	const LocalElement& myElement = massIt.first;
	const std::string subdetectorName = myElement.subdetectorName();
	if (subdetectorName == xml_OTST) { isOTST = true; }
	if (isOTST && subdetectorName != xml_OTST) {
	  std::cout << "VERY IMPORTANT!! NOT SUPPORTED: SAME volume contains OTST element(s), but also other element(s) assigned to different subdetector(s)." << std::endl;
	}
      }


      if ((iter->getZOffset() + iter->getZLength()) > 0) { // && shapename.str() != "supportR1191Z1325")
        if ( iter->getLocalMasses().size() ) {

	  // Create composite
          if (!isOTST) { c.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter)); }
	  else { otstComposites.push_back(createComposite(matname.str(), compositeDensity(*iter), *iter)); }

	  // TO DO : CALCULATION OF OUTERMOST SHAPES BOUNDARIES
	  double startEndcaps;
	  if (!isPixelTracker) startEndcaps = xml_outerTrackerEndcapsMinZ;
	  else startEndcaps = xml_innerTrackerEndcapsMinZ;  // Tilted Inner Tracker : startEndcaps = xml_innerTiltedTrackerEndcapsMinZ;


	  // BARREL supports
	  if ( (!isPixelTracker && iter->getZOffset() < startEndcaps)
	       || (isPixelTracker && ((iter->getZOffset() + iter->getZLength() / 2.0) < startEndcaps))
	       || isOTST
	       ) {

	    shape.name_tag = shapename.str();
	    shape.dz = iter->getZLength() / 2.0;
	    shape.rmin = iter->getInnerRadius();
	    shape.rmax = shape.rmin + iter->getRWidth();
	    if (!isOTST) { s.push_back(shape); }
	    else { otstShapes.push_back(shape); }
	    if (shape.rmin < supportBarrelRMin) supportBarrelRMin = shape.rmin;
	    if (shape.rmax > supportBarrelRMax) supportBarrelRMax = shape.rmax;

	    logic.name_tag = shapename.str();
	    logic.shape_tag = (!isOTST ? trackerXmlTags.nspace : xml_otstident) + ":" + shapename.str();
	    logic.material_tag = (!isOTST ? trackerXmlTags.nspace : xml_otstident) + ":" + matname.str();
	    if (!isOTST) { l.push_back(logic); }
	    else { otstLogic.push_back(logic); }

	    if (!isOTST) { pos.parent_tag = xml_pixbarident + ":" + trackerXmlTags.bar; }
	    else { pos.parent_tag = trackerXmlTags.nspace + ":" + xml_tracker; }
	    pos.child_tag = logic.shape_tag;
	    pos.trans.dz = iter->getZOffset() + shape.dz;
	    if (!isOTST) { p.push_back(pos); }
	    else { otstPositions.push_back(pos); }
	    pos.copy = 2;
	    pos.trans.dz = -pos.trans.dz;
	    pos.rotref = trackerXmlTags.nspace + ":" + xml_Y180;
	    if (!isOTST) { p.push_back(pos);  }
	    else { otstPositions.push_back(pos); }
	  }

	  // ENDCAPS supports
	  else {

	    // cut in 2 the services that belong to both Barrel and Endcaps mother volumes
	    if (isPixelTracker && iter->getZOffset() < startEndcaps) {
	      std::ostringstream shapenameBarrel, shapenameEndcaps;
	      shapenameBarrel << xml_base_lazy << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(iter->getZLength() / 2.0 + iter->getZOffset()) << "BarrelPart";
	      shapenameEndcaps << xml_base_lazy << "R" << (int)(iter->getInnerRadius()) << "Z" << (int)(iter->getZLength() / 2.0 + iter->getZOffset()) << "EndcapsPart";

	      // Barrel part
	      shape.name_tag = shapenameBarrel.str();
	      shape.dz = (startEndcaps - iter->getZOffset()) / 2.0 - xml_epsilon;
	      shape.rmin = iter->getInnerRadius();
	      shape.rmax = shape.rmin + iter->getRWidth();
	      s.push_back(shape);
	      if (shape.rmin < supportBarrelRMin) supportBarrelRMin = shape.rmin;
	      if (shape.rmax > supportBarrelRMax) supportBarrelRMax = shape.rmax;

	      logic.name_tag = shapenameBarrel.str();
	      logic.shape_tag = trackerXmlTags.nspace + ":" + shapenameBarrel.str();
	      logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
	      l.push_back(logic);

	      pos.parent_tag = xml_pixbarident + ":" + trackerXmlTags.bar;
	      pos.child_tag = logic.shape_tag;
	      pos.trans.dz = iter->getZOffset() + shape.dz + xml_epsilon;
	      p.push_back(pos);
	      pos.copy = 2;
	      pos.trans.dz = -pos.trans.dz;
	      pos.rotref = trackerXmlTags.nspace + ":" + xml_Y180;
	      p.push_back(pos);

	      pos.copy = 1;
	      pos.rotref.clear();

	      // Endcaps part
	      shape.name_tag = shapenameEndcaps.str();
	      shape.dz = (iter->getZOffset() + iter->getZLength() - startEndcaps) / 2.0 - xml_epsilon;
	      shape.rmin = iter->getInnerRadius();
	      shape.rmax = shape.rmin + iter->getRWidth();
	      s.push_back(shape);    
	      if (shape.rmin < supportEndcapsRMin) supportEndcapsRMin = shape.rmin;
	      if (shape.rmax > supportEndcapsRMax) supportEndcapsRMax = shape.rmax;

	      logic.name_tag = shapenameEndcaps.str();
	      logic.shape_tag = trackerXmlTags.nspace + ":" + shapenameEndcaps.str();
	      logic.material_tag = trackerXmlTags.nspace + ":" + matname.str();
	      l.push_back(logic);

	      pos.parent_tag = xml_pixfwdident + ":" + trackerXmlTags.fwd;
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



  /* 
   * PRIVATE
   */


  // Add a tilted module rotation to the map of rotations
  // This rotation is defined as the following : axial rotation, with axis X (local X of the module), and angle tiltAngle.
  void Extractor::addTiltedModuleRot( std::map<std::string,Rotation>& rotations, double tiltAngle) {
    std::ostringstream tilt;
    tilt << round(tiltAngle);

    std::string positiveZName = xml_positive_z_tilted_mod_rot + tilt.str();
    if (rotations.count(positiveZName) == 0) {
      Rotation rot;
      rot.name = positiveZName;
      rot.thetax = 90.0;
      rot.phix = 0.0;
      rot.thetay = 90.0 + tiltAngle;
      rot.phiy = 90.0;
      rot.thetaz = tiltAngle;
      rot.phiz = 90.0;
      rotations.insert(std::pair<const std::string,Rotation>(rot.name,rot));
    }

    std::string negativeZName = xml_negative_z_tilted_mod_rot + tilt.str();
    if (rotations.count(negativeZName) == 0) {
      Rotation rot;
      rot.name = negativeZName;
      rot.thetax = 90.0;
      rot.phix = 0.0;
      rot.thetay = 90.0 - tiltAngle;
      rot.phiy = 90.0;
      rot.thetaz = tiltAngle;
      rot.phiz = 270.0;
      rotations.insert(std::pair<const std::string,Rotation>(rot.name,rot));
    }
  }


  /* Add rotation around local Y axis.
     For example, this is used to skew modules: modules are rotated around their local Y axis.
  */
  void Extractor::addRotationAroundZAxis(std::map<std::string,Rotation>& storedRotations, 
					 const std::string rotationName,
					 const double rotationAngleInRad) const {

    if (storedRotations.find(rotationName) == storedRotations.end()) {
      Rotation rot;
      rot.name = rotationName;
      rot.thetax = 90.;
      rot.phix = rotationAngleInRad * 180. / M_PI;
      rot.thetay = 90.;
      rot.phiy = rot.phix + 90.;
      rot.thetaz = 0.;
      rot.phiz = 0.;
      storedRotations.insert(std::pair<std::string, Rotation>(rotationName, rot));
    }
  } 



  void Extractor::createAndStoreDDTrackerAngularAlgorithmBlock(std::vector<AlgoInfo>& storedAlgorithmBlocks,
							       const std::string nameSpace, 
							       const std::string parentName,
							       const std::string childName,
							       const double startAngleInRad,
							       const double rangeAngleInRad,
							       const double radius,
							       const XYZVector& center,
							       const int numCopies,
							       const int startCopyNumber,
							       const int copyNumberIncrement) {
    AlgoInfo myAlgo; // would obviously be better to build the algo directly here instead of having a constructor with no argument!

    // Algo name
    myAlgo.name = xml_angular_algo;

    // Parent volume name
    myAlgo.parent = nameSpace + ":" + parentName;

    // Child volume name
    myAlgo.parameters.push_back(stringParam(xml_childparam, nameSpace + ":" + childName));
	 
    // Volume with copy number 1, center phi angle
    myAlgo.parameters.push_back(numericParam(xml_startangle, any2str(startAngleInRad * 180. / M_PI, xml_angle_precision) + "*deg"));
	 
    // Difference between first and last volumes' centers phi angles
    myAlgo.parameters.push_back(numericParam(xml_rangeangle, any2str(rangeAngleInRad * 180. / M_PI, xml_angle_precision) + "*deg"));

    // All volumes are assumed to be placed at same radius in that algorithm
    myAlgo.parameters.push_back(numericParam(xml_radius, any2str(radius) + "*mm"));

    // Shift of children with respect to parent, AFTER rotation is done
    myAlgo.parameters.push_back(vectorParam(center.X(), center.Y(), center.Z()));

    // Number of children copies
    myAlgo.parameters.push_back(numericParam(xml_nmods, any2str(numCopies)));
	  
    // Start copy number
    myAlgo.parameters.push_back(numericParam(xml_startcopyno, any2str(startCopyNumber)));

    // Copy number increment
    myAlgo.parameters.push_back(numericParam(xml_incrcopyno, any2str(copyNumberIncrement)));

    // Store algo
    storedAlgorithmBlocks.push_back(myAlgo);
  }


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
      if (!nosensors || ( (it.first.compare(xml_sensor_silicon) != 0) && (it.first.compare(xml_sensor_LYSO) != 0)) ) { // SenSi element is not treated here. 
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

    comp.normalizedRIRatioPerMechanicalCategory = mp.getNormalizedRIRatioPerMechanicalCategory();

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
    res << xml_algorithm_vector_open << x << "*mm, " << y << "*mm, " << z << "*mm" << xml_algorithm_vector_close;
    return res.str();
  }

 std::string Extractor::arbitraryLengthVector(std::string name, std::vector<double> invec){
    std::ostringstream res;
    std::string vector_opening = "<Vector name=\""+name+"\" type=\"numeric\" nEntries=\""+std::to_string(invec.size())+"\">";
    if(invec.size() > 0){
      res << vector_opening << invec.at(0);
      for(unsigned int i=1; i<invec.size();i++){
        res << "," << invec.at(i);
      }
      res << xml_algorithm_vector_close;
    }
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
        if (it->first.compare(xml_sensor_silicon) != 0 && it->first.compare(xml_sensor_LYSO) != 0) m += it->second;
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
    //radEstimate = compute_X0(Z, A_estimate);  
    //double intEstimate = pow(A_estimate, 1./3.)*A_a+A_b;

    // On CMSSW side, radiation lengths and nuclear interaction lengths will be recomputed from Z, A and density.
    // Now, a question is whether the computed radiation and interaction lengths on CMSSW side (similar to radEstimate and intEstimate) are closed to the values we initially had at hand (double radiationLength and interactionLength) !
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

  

  const double ModuleComplexHelpers::computeExpandedModWidth(const double moduleWidth, 
							     const double serviceHybridWidth, 
							     const double deadAreaExtraWidth,
							     const double chipNegativeXExtraWidth,
							     const double chipPositiveXExtraWidth) {
    const double totalServiceHybridWidth = 2. * serviceHybridWidth; // OT case
    const double totalDeadAreaExtraWidth = 2. * deadAreaExtraWidth; // IT case: around sensor
    const double totalChipExtraWidth = 2. * MAX(chipNegativeXExtraWidth, chipPositiveXExtraWidth); // IT case: around chip
    const double expandedModWidth = moduleWidth + totalServiceHybridWidth + MAX(totalDeadAreaExtraWidth, totalChipExtraWidth);
    return expandedModWidth;
  }

  const double ModuleComplexHelpers::computeExpandedModLength(const double moduleLength, 
							      const double frontEndHybridWidth, 
							      const double deadAreaExtraLength) {
    const double totalFrontEndHybridWidth = 2. * frontEndHybridWidth; // OT case
    const double totalDeadAreaExtraLength = 2.* deadAreaExtraLength;  // IT case
    const double expandedModLength = moduleLength + totalFrontEndHybridWidth + totalDeadAreaExtraLength;
    return expandedModLength;
  }





  ModuleComplex::ModuleComplex(std::string moduleName,
                               std::string parentName,
                               ModuleCap&  modcap        ) : modulecap(modcap),
							     module(modcap.getModule()),
							     moduleId(moduleName),
                                                             parentId(parentName),
							     modThickness(module.thickness()),
                                                             sensorThickness(module.sensorThickness()),
                                                             sensorDistance(module.dsDistance()),
							     modWidth(module.area()/module.length()),
                                                             modLength(module.length() + module.outerSensorExtraLength()),
                                                             modLengthDoubleSens(2*module.length() + module.outerSensorExtraLength()+module.centralDeadAreaLength()),
                                                             frontEndHybridWidth(module.frontEndHybridWidth()),
                                                             serviceHybridWidth(module.serviceHybridWidth()),
                                                             hybridThickness(module.hybridThickness()),
                                                             supportPlateThickness(module.supportPlateThickness()),
                                                             chipThickness(module.chipThickness()),
							     deadAreaExtraLength(module.deadAreaExtraLength()),
							     deadAreaExtraWidth(module.deadAreaExtraWidth()),
                                                             centralDeadAreaLength(module.centralDeadAreaLength()),
							     chipNegativeXExtraWidth(module.chipNegativeXExtraWidth()),
							     chipPositiveXExtraWidth(module.chipPositiveXExtraWidth()),
                                                             hybridTotalMass(0.),
                                                             hybridTotalVolume_mm3(-1.),
                                                             hybridFrontAndBackVolume_mm3(-1.),
                                                             hybridLeftAndRightVolume_mm3(-1.),
							     deadAreaTotalVolume_mm3(-1.),
                                                             moduleMassWithoutSensors_expected(0.),
                                                             expandedModWidth(ModuleComplexHelpers::computeExpandedModWidth(modWidth, 
															    serviceHybridWidth, 
															    deadAreaExtraWidth,
															    chipNegativeXExtraWidth,
															    chipPositiveXExtraWidth)
									      ),
		      
                                                             expandedModLength(ModuleComplexHelpers::computeExpandedModLength(modLength,
															      frontEndHybridWidth,
															      deadAreaExtraLength)
									       ),
                                                             expandedModLengthPixDoubleSens(ModuleComplexHelpers::computeExpandedModLength(modLengthDoubleSens,
															      frontEndHybridWidth,
															      deadAreaExtraLength)
									       ),
							     center(module.center()),
                                                             normal(module.normal()),
                                                             prefix_material(xml_hybrid_comp) {
    if (!module.isPixelModule()) {
      if (!module.isTimingModule()) expandedModThickness = sensorDistance + 2.0 * (supportPlateThickness + sensorThickness);
      else expandedModThickness = sensorThickness + 2.0 * MAX(supportPlateThickness, hybridThickness);
      prefix_xmlfile = xml_fileident + ":";
      nTypes = 9;
    }
    else {
      expandedModThickness = sensorThickness + 2.0 * MAX(chipThickness, hybridThickness);
      prefix_xmlfile = xml_PX_fileident + ":";
      nTypes = 11;
    }
  }

  ModuleComplex::~ModuleComplex() {
    std::vector<Volume*>::const_iterator vit;
    for ( vit = volumes.begin(); vit != volumes.end(); vit++ ) delete *vit;
  }

  void ModuleComplex::buildSubVolumes() {
    Volume* vol[nTypes];
    if (!module.isPixelModule()) {

      if (!module.isTimingModule()) {
	//                                                   (FOR XML EXPORT ONLY)
	//                                           MATERIAL ASSIGNMENT IN TARGET VOLUMES
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
	//  SupportPlate(8) thickness is 0 for 2S modules
    
    
	// Combinations
	// Volume (34) is volume (3) + volume (4)
	// Volume (56) is volume (5) + volume (6)
	// Volume (0) is volume (3) + volume (4) + volume (5) + volume (6)
    

	//Unused pointers
	vol[xml_HybridFBLR_0] = 0;
	vol[xml_InnerSensor]  = 0;
	vol[xml_OuterSensor]  = 0;

	double dx = serviceHybridWidth;              
	double dy = modLength; 
	double dz = hybridThickness;  
	double posx = (modWidth+serviceHybridWidth)/2.;
	double posy = 0.;
	double posz = 0.;
	// Hybrid FrontSide Volume
	vol[xml_HybridFront] = new Volume(moduleId+"FSide",xml_HybridFront,parentId,dx,dy,dz,posx,posy,posz);

	posx = -(modWidth+serviceHybridWidth)/2.;
	posy = 0.;
	posz = 0.;
	// Hybrid BackSide Volume
	vol[xml_HybridBack] = new Volume(moduleId+"BSide",xml_HybridBack,parentId,dx,dy,dz,posx,posy,posz);

	dx = modWidth+2*serviceHybridWidth;  
	dy = frontEndHybridWidth;
	posx = 0.;
	posy = (modLength+frontEndHybridWidth)/2.;
	posz = 0.;
	// Hybrid LeftSide Volume
	vol[xml_HybridLeft] = new Volume(moduleId+"LSide",xml_HybridLeft,parentId,dx,dy,dz,posx,posy,posz);

	posx = 0.;
	posy = -(modLength+frontEndHybridWidth)/2.;
	posz = 0.;
	// Hybrid RightSide Volume
	vol[xml_HybridRight] = new Volume(moduleId+"RSide",xml_HybridRight,parentId,dx,dy,dz,posx,posy,posz);

	dx = modWidth; 
	dy = modLength; 
	posx = 0.;
	posy = 0.;
	posz = 0.;
	// Hybrid Between Volume
	vol[xml_HybridBetween] = new Volume(moduleId+"Between",xml_HybridBetween,parentId,dx,dy,dz,posx,posy,posz);

	dx = expandedModWidth;  
	dy = expandedModLength; 
	dz = supportPlateThickness;
	posx = 0.;
	posy = 0.;
	posz = - ( ( sensorDistance + supportPlateThickness )/2. + sensorThickness ); 
	// SupportPlate
	vol[xml_SupportPlate] = new Volume(moduleId+"SupportPlate",xml_SupportPlate,parentId,dx,dy,dz,posx,posy,posz);
      }

      else {
	//                                                   TIMING MODULE
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
	//          
	//  ====== ================ ====== Hybrids    +----> x
	//         ---------------- Sensor
	//  ============================== 
	//          SupportPlate(8)                      
	//
	//  R(6) and L(5) are Front-End Hybrids.
	//  B(4) and F(3) are Service Hybdrids.
	//
  
    
	//Unused pointers
	vol[xml_HybridFBLR_0] = 0;
	vol[xml_InnerSensor]  = 0;
	vol[xml_OuterSensor]  = 0;

	double dx = serviceHybridWidth;              
	double dy = modLength; 
	double dz = hybridThickness;  
	double posx = (modWidth+serviceHybridWidth)/2.;
	double posy = 0.;
	double posz = sensorThickness / 2. + hybridThickness / 2.;
	// Hybrid FrontSide Volume
	vol[xml_HybridFront] = new Volume(moduleId+"FSide",xml_HybridFront,parentId,dx,dy,dz,posx,posy,posz);

	posx = -(modWidth+serviceHybridWidth)/2.;
	posy = 0.;
	posz = sensorThickness / 2. + hybridThickness / 2.; 
	// Hybrid BackSide Volume
	vol[xml_HybridBack] = new Volume(moduleId+"BSide",xml_HybridBack,parentId,dx,dy,dz,posx,posy,posz);

	dx = modWidth+2*serviceHybridWidth;  
	dy = frontEndHybridWidth;
	posx = 0.;
	posy = (modLength+frontEndHybridWidth)/2.;
	posz = sensorThickness / 2. + hybridThickness / 2.; 
	// Hybrid LeftSide Volume
	vol[xml_HybridLeft] = new Volume(moduleId+"LSide",xml_HybridLeft,parentId,dx,dy,dz,posx,posy,posz);

	posx = 0.;
	posy = -(modLength+frontEndHybridWidth)/2.;
	posz = sensorThickness / 2. + hybridThickness / 2.; 
	// Hybrid RightSide Volume
	vol[xml_HybridRight] = new Volume(moduleId+"RSide",xml_HybridRight,parentId,dx,dy,dz,posx,posy,posz);

	dx = modWidth; 
	dy = modLength; 
	posx = 0.;
	posy = 0.;
	posz = sensorThickness / 2. + hybridThickness / 2.; 
	// Hybrid Between Volume
	vol[xml_HybridBetween] = new Volume(moduleId+"Between",xml_HybridBetween,parentId,dx,dy,dz,posx,posy,posz);

	dx = expandedModWidth;  
	dy = expandedModLength; 
	dz = supportPlateThickness;
	posx = 0.;
	posy = 0.;
	posz = - sensorThickness / 2. - supportPlateThickness / 2.;
	// SupportPlate
	vol[xml_SupportPlate] = new Volume(moduleId+"SupportPlate",xml_SupportPlate,parentId,dx,dy,dz,posx,posy,posz);
      }
    }

    else {
      //                                                   (FOR XML EXPORT ONLY)
      //                                           MATERIAL ASSIGNMENT IN TARGET VOLUMES
      //                                                      PIXEL MODULE
      //
      //  Top View 
      //
      //    --------------------       
      //    |        (7)       |
      //    --------------------        y
      //    |   |          |   |        ^
      //    |   |          |   |        |
      //    |   |          |   |        |
      //    |   |          |   |        |
      //    |   |          |   |        |
      //    |   |  Sensor  |   |        |
      //    |(6)|    (2)   |(5)|        |
      //    |   |          |   |        |
      //    |   |          |   |        |
      //    |   |          |   |        +----> x
      //    --------------------    
      //    |        (8)       |                        
      //    --------------------                        z
      //                                                ^
      //  Side View                                     |
      //         ================          Hybrid  (1)  +----> x
      //     (6) ---------------- (5)      Sensor  (2)
      //     =========================     Chip    (3)
      //
      // Chip(3) volume can contain Bumps and any other material (supports, etc) for simplification.
      //
      // Volumes (5), (6), (7), and (8) are volumes of same thickness as the sensor, located around the sensor.
      // They are used for Si inactive areas.
      
      // Combinations
      // Volume (4) is volume (5) + volume (6) + volume (7) + volume (8).
      // Volume (0) does not target anything.


      //Unused pointers
      vol[xml_PixelModuleNull] = 0;
      vol[xml_PixelModuleDeadArea] = 0;
      vol[xml_PixelModuleDeadAreaFrontOfCentre] = 0;
      vol[xml_PixelModuleDeadAreaBackOfCentre] = 0;

      // Hybrid Volume (Top Inactive)
      const double myHybridWidth = modWidth;
      const double myHybridLength = modLength;
      /*if(module.numSensors()==2){
          myHybridLength = modLengthDoubleSens;
      }*/
      const double myHybridThickness = hybridThickness; 
      const double myHybridPosX = 0.;
      const double myHybridPosY = 0.;
      const double myHybridPosZ = sensorThickness / 2. + hybridThickness / 2.;
      vol[xml_PixelModuleHybrid] = new Volume(moduleId + "Hybrid", xml_PixelModuleHybrid, parentId, 
					  myHybridWidth, myHybridLength, myHybridThickness, 
					  myHybridPosX, myHybridPosY, myHybridPosZ);

      // Dead area Right (Inactive silicon around sensor)
      const double myDeadAreaRightWidth = deadAreaExtraWidth;
      const double myDeadAreaRightLength = modLength;
      /*if(module.numSensors()==2){
          myDeadAreaRightLength = modLengthDoubleSens;
      }*/
      const double myDeadAreaRightThickness = sensorThickness; 
      const double myDeadAreaRightPosX = (modWidth + deadAreaExtraWidth) / 2.;
      const double myDeadAreaRightPosY = 0.;
      const double myDeadAreaRightPosZ = 0.;
      vol[xml_PixelModuleDeadAreaRight] = new Volume(moduleId + "DeadAreaRight", xml_PixelModuleDeadAreaRight, parentId, 
					  myDeadAreaRightWidth, myDeadAreaRightLength, myDeadAreaRightThickness, 
					  myDeadAreaRightPosX, myDeadAreaRightPosY, myDeadAreaRightPosZ);



      // Dead area Left (Inactive silicon around sensor)
      const double myDeadAreaLeftWidth = deadAreaExtraWidth;
      const double myDeadAreaLeftLength = modLength;
      const double myDeadAreaLeftThickness = sensorThickness; 
      const double myDeadAreaLeftPosX = -(modWidth + deadAreaExtraWidth) / 2.;
      const double myDeadAreaLeftPosY = 0.;
      const double myDeadAreaLeftPosZ = 0.;
      vol[xml_PixelModuleDeadAreaLeft] = new Volume(moduleId + "DeadAreaLeft", xml_PixelModuleDeadAreaLeft, parentId, 
					  myDeadAreaLeftWidth, myDeadAreaLeftLength, myDeadAreaLeftThickness, 
					  myDeadAreaLeftPosX, myDeadAreaLeftPosY, myDeadAreaLeftPosZ);


      
      // Dead area Front (Inactive silicon around sensor)
      const double myDeadAreaFrontWidth = modWidth + 2. * deadAreaExtraWidth;
      const double myDeadAreaFrontLength = deadAreaExtraLength;
      const double myDeadAreaFrontThickness = sensorThickness; 
      const double myDeadAreaFrontPosX = 0.;
      const double myDeadAreaFrontPosY = (modLength + deadAreaExtraLength) / 2.;
      const double myDeadAreaFrontPosZ = 0.;
      vol[xml_PixelModuleDeadAreaFront] = new Volume(moduleId + "DeadAreaFront", xml_PixelModuleDeadAreaFront, parentId, 
					  myDeadAreaFrontWidth, myDeadAreaFrontLength, myDeadAreaFrontThickness, 
					  myDeadAreaFrontPosX, myDeadAreaFrontPosY, myDeadAreaFrontPosZ);

      // Dead area Front of centre (Inactive silicon around sensor)
      if(module.numSensors()==2){
        const double myDeadAreaFrontOfCentreWidth = modWidth;
        const double myDeadAreaFrontOfCentreLength = deadAreaExtraLength;//Assume the same silicon dead area as for the front - this could be made more general
        const double myDeadAreaFrontOfCentreThickness = sensorThickness; 
        const double myDeadAreaFrontOfCentrePosX = 0.;
        const double myDeadAreaFrontOfCentrePosY = (centralDeadAreaLength - deadAreaExtraLength) / 2.; //NB centralDeadAreaLength is total central dead area + air gap.
        const double myDeadAreaFrontOfCentrePosZ = 0.;
        vol[xml_PixelModuleDeadAreaFrontOfCentre] = new Volume(moduleId + "DeadAreaFrontOfCentre", xml_PixelModuleDeadAreaFrontOfCentre, parentId, 
          					  myDeadAreaFrontOfCentreWidth, myDeadAreaFrontOfCentreLength, myDeadAreaFrontOfCentreThickness, 
	        				  myDeadAreaFrontOfCentrePosX, myDeadAreaFrontOfCentrePosY, myDeadAreaFrontOfCentrePosZ);
      }


      // Dead area Back (Inactive silicon around sensor)
      const double myDeadAreaBackWidth = modWidth + 2. * deadAreaExtraWidth;
      const double myDeadAreaBackLength = deadAreaExtraLength;
      const double myDeadAreaBackThickness = sensorThickness; 
      const double myDeadAreaBackPosX = 0.;
      const double myDeadAreaBackPosY = -(modLength + deadAreaExtraLength) / 2.;
      const double myDeadAreaBackPosZ = 0.;
      vol[xml_PixelModuleDeadAreaBack] = new Volume(moduleId + "DeadAreaBack", xml_PixelModuleDeadAreaBack, parentId, 
					  myDeadAreaBackWidth, myDeadAreaBackLength, myDeadAreaBackThickness, 
					  myDeadAreaBackPosX, myDeadAreaBackPosY, myDeadAreaBackPosZ);


      // Dead area Back of centre (Inactive silicon around sensor)
      if(module.numSensors()==2){
        const double myDeadAreaBackOfCentreWidth = modWidth;
        const double myDeadAreaBackOfCentreLength = deadAreaExtraLength;//Assume the same silicon dead area as for the front - this could be made more general
        const double myDeadAreaBackOfCentreThickness = sensorThickness; 
        const double myDeadAreaBackOfCentrePosX = 0.;
        const double myDeadAreaBackOfCentrePosY = -(centralDeadAreaLength - deadAreaExtraLength) / 2.; //NB centralDeadAreaLength is total central dead area + air gap
        const double myDeadAreaBackOfCentrePosZ = 0.;
        vol[xml_PixelModuleDeadAreaBackOfCentre] = new Volume(moduleId + "DeadAreaBackOfCentre", xml_PixelModuleDeadAreaBackOfCentre, parentId, 
  		         			  myDeadAreaBackOfCentreWidth, myDeadAreaBackOfCentreLength, myDeadAreaBackOfCentreThickness, 
	         				  myDeadAreaBackOfCentrePosX, myDeadAreaBackOfCentrePosY, myDeadAreaBackOfCentrePosZ);

      }


      // Chip Volume (Bottom Inactive)
      const double myChipWidth = modWidth + chipNegativeXExtraWidth + chipPositiveXExtraWidth;
      const double myChipLength = modLength;
      const double myChipThickness = chipThickness; 
      const double myChipPosX = (chipPositiveXExtraWidth - chipNegativeXExtraWidth) / 2.;
      const double myChipPosY = 0.;
      const double myChipPosZ = - sensorThickness / 2. - chipThickness / 2.;
      vol[xml_PixelModuleChip] = new Volume(moduleId + "Chip", xml_PixelModuleChip, parentId,  
					myChipWidth, myChipLength, myChipThickness, 
					myChipPosX, myChipPosY, myChipPosZ);
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
    XYZVector mx = (0.5*( module.basePoly().getVertex(2) + module.basePoly().getVertex(3) ) - center).Unit() ;
    XYZVector my = (0.5*( module.basePoly().getVertex(1) + module.basePoly().getVertex(2) ) - center).Unit() ;

    // new vertexes after expansion due to hybrid volumes
    const int npoints = 5; // v0,v1,v2,v3,v4(=v0)
    XYZVector v[npoints-1];
    v[0] = module.center() - expandedModWidth/2. * mx - expandedModLength/2. * my;
    v[1] = module.center() - expandedModWidth/2. * mx + expandedModLength/2. * my;
    v[2] = module.center() + expandedModWidth/2. * mx + expandedModLength/2. * my;
    v[3] = module.center() + expandedModWidth/2. * mx - expandedModLength/2. * my;
    /*if(module.numSensors()==2){
        v[0] = module.center() - expandedModWidth/2. * mx - expandedModLengthPixDoubleSens/2. * my;
        v[1] = module.center() - expandedModWidth/2. * mx + expandedModLengthPixDoubleSens/2. * my;
        v[2] = module.center() + expandedModWidth/2. * mx + expandedModLengthPixDoubleSens/2. * my;
        v[3] = module.center() + expandedModWidth/2. * mx - expandedModLengthPixDoubleSens/2. * my;
    }*/

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
	if ( el->targetVolume() == xml_InnerSensor  ||
	     el->targetVolume() == xml_OuterSensor     ) {
	  continue; // Still to skip sensors, in case not detected by the component name.
	} else if ( el->targetVolume() >= nTypes   &&
		    el->targetVolume() != xml_HybridFB &&
		    el->targetVolume() != xml_HybridLR &&
		    el->targetVolume() != xml_HybridFBLR_3456  ) {
	  std::cerr << "!!!! ERROR !!!! : Found unexpected targetVolume." << std::endl;
	  std::cerr << "targetVolume " << el->targetVolume() << " is not supported for Outer Barrel modules. Exit." << std::endl;
	  std::exit(1);
	}

	moduleMassWithoutSensors_expected += el->quantityInGrams(module);

	if ( el->targetVolume() == xml_HybridFront   ||
	     el->targetVolume() == xml_HybridBack    ||
	     el->targetVolume() == xml_HybridLeft    ||
	     el->targetVolume() == xml_HybridRight   ||
	     el->targetVolume() == xml_HybridBetween ||
	     el->targetVolume() == xml_SupportPlate     ) {
          vol[el->targetVolume()]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module));
          vol[el->targetVolume()]->addMass(el->quantityInGrams(module));
	} 
	else if ( el->targetVolume() == xml_HybridFB ) { 
          if (hybridFrontAndBackVolume_mm3 < 0) { // Need only once
            hybridFrontAndBackVolume_mm3 = vol[xml_HybridFront]->getVolume()
	      + vol[xml_HybridBack]->getVolume();
          }
          vol[xml_HybridFront]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_HybridFront]->getVolume()/hybridFrontAndBackVolume_mm3);
          vol[xml_HybridBack]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_HybridBack]->getVolume()/hybridFrontAndBackVolume_mm3);
          vol[xml_HybridFront]->addMass(el->quantityInGrams(module)*vol[xml_HybridFront]->getVolume()/hybridFrontAndBackVolume_mm3);
          vol[xml_HybridBack]->addMass(el->quantityInGrams(module)*vol[xml_HybridBack]->getVolume()/hybridFrontAndBackVolume_mm3);
	} 
	else if ( el->targetVolume() == xml_HybridLR ) {
          if (hybridLeftAndRightVolume_mm3 < 0) { // Need only once
            hybridLeftAndRightVolume_mm3 = vol[xml_HybridLeft]->getVolume()
	      + vol[xml_HybridRight]->getVolume();
          }
          vol[xml_HybridLeft]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_HybridLeft]->getVolume()/hybridLeftAndRightVolume_mm3);
          vol[xml_HybridRight]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_HybridRight]->getVolume()/hybridLeftAndRightVolume_mm3);
          vol[xml_HybridLeft]->addMass(el->quantityInGrams(module)*vol[xml_HybridLeft]->getVolume()/hybridLeftAndRightVolume_mm3);
          vol[xml_HybridRight]->addMass(el->quantityInGrams(module)*vol[xml_HybridRight]->getVolume()/hybridLeftAndRightVolume_mm3);
	} 
	else if ( el->targetVolume() == xml_HybridFBLR_0 || el->targetVolume() == xml_HybridFBLR_3456 ) { // Uniformly Distribute
	  if (hybridTotalVolume_mm3 < 0) { // Need only once
            hybridTotalVolume_mm3 = vol[xml_HybridFront]->getVolume()
	      + vol[xml_HybridBack]->getVolume()
	      + vol[xml_HybridLeft]->getVolume()
	      + vol[xml_HybridRight]->getVolume();
          }
          vol[xml_HybridFront]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_HybridFront]->getVolume()/hybridTotalVolume_mm3);
          vol[xml_HybridBack]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_HybridBack]->getVolume()/hybridTotalVolume_mm3);
          vol[xml_HybridLeft]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_HybridLeft]->getVolume()/hybridTotalVolume_mm3);
          vol[xml_HybridRight]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_HybridRight]->getVolume()/hybridTotalVolume_mm3);    
          // Uniform density distribution and consistent with total mass
          vol[xml_HybridFront]->addMass(el->quantityInGrams(module)*vol[xml_HybridFront]->getVolume()/hybridTotalVolume_mm3); 
          vol[xml_HybridBack]->addMass(el->quantityInGrams(module)*vol[xml_HybridBack]->getVolume()/hybridTotalVolume_mm3);   
          vol[xml_HybridLeft]->addMass(el->quantityInGrams(module)*vol[xml_HybridLeft]->getVolume()/hybridTotalVolume_mm3);   
          vol[xml_HybridRight]->addMass(el->quantityInGrams(module)*vol[xml_HybridRight]->getVolume()/hybridTotalVolume_mm3);
  	}
      }

      // PIXEL MODULE
      else {
	if ( el->targetVolume() == xml_PixelModuleSensor ) {
	  continue; // Still to skip sensors, in case not detected by the component name.
	}
	else if ( el->targetVolume() != xml_PixelModuleHybrid &&
		  el->targetVolume() != xml_PixelModuleChip &&
		  el->targetVolume() != xml_PixelModuleDeadAreaRight && 
		  el->targetVolume() != xml_PixelModuleDeadAreaFrontOfCentre && 
		  el->targetVolume() != xml_PixelModuleDeadAreaLeft && 
		  el->targetVolume() != xml_PixelModuleDeadAreaBackOfCentre && 
		  el->targetVolume() != xml_PixelModuleDeadAreaFront && 
		  el->targetVolume() != xml_PixelModuleDeadAreaBack && 
		  el->targetVolume() != xml_PixelModuleDeadArea
		  ) {
	  throw PathfulException("!!!! ERROR !!!! : Found unexpected targetVolume, not supported for Pixel Barrel modules.");
	}
	moduleMassWithoutSensors_expected += el->quantityInGrams(module);

	if ( el->targetVolume() == xml_PixelModuleHybrid   ||
	     el->targetVolume() == xml_PixelModuleChip ||
	     el->targetVolume() == xml_PixelModuleDeadAreaRight ||     
	     el->targetVolume() == xml_PixelModuleDeadAreaFrontOfCentre ||     
	     el->targetVolume() == xml_PixelModuleDeadAreaLeft ||
	     el->targetVolume() == xml_PixelModuleDeadAreaBackOfCentre ||
	     el->targetVolume() == xml_PixelModuleDeadAreaFront ||
	     el->targetVolume() == xml_PixelModuleDeadAreaBack
	     ) {
          vol[el->targetVolume()]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module));
          vol[el->targetVolume()]->addMass(el->quantityInGrams(module));
	}

	else if ( el->targetVolume() == xml_PixelModuleDeadArea) { // Uniformly Distribute

	  if (deadAreaTotalVolume_mm3 < 0) { // Need only once
            deadAreaTotalVolume_mm3 = vol[xml_PixelModuleDeadAreaRight]->getVolume()
	      + vol[xml_PixelModuleDeadAreaLeft]->getVolume()
	      + vol[xml_PixelModuleDeadAreaFront]->getVolume()
	      + vol[xml_PixelModuleDeadAreaBack]->getVolume();
            if (module.numSensors()==2){
              deadAreaTotalVolume_mm3 += vol[xml_PixelModuleDeadAreaFrontOfCentre]->getVolume() 
	        + vol[xml_PixelModuleDeadAreaBackOfCentre]->getVolume();
            }
          }

          vol[xml_PixelModuleDeadAreaRight]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaRight]->getVolume()/deadAreaTotalVolume_mm3);
          vol[xml_PixelModuleDeadAreaLeft]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaLeft]->getVolume()/deadAreaTotalVolume_mm3);
          if(module.numSensors()==2){
            vol[xml_PixelModuleDeadAreaFrontOfCentre]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaFrontOfCentre]->getVolume()/deadAreaTotalVolume_mm3);
            vol[xml_PixelModuleDeadAreaBackOfCentre]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaBackOfCentre]->getVolume()/deadAreaTotalVolume_mm3);
          }


          vol[xml_PixelModuleDeadAreaFront]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaFront]->getVolume()/deadAreaTotalVolume_mm3);
          vol[xml_PixelModuleDeadAreaBack]->addMaterial(el->elementName(), el->componentName(), el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaBack]->getVolume()/deadAreaTotalVolume_mm3);

          // Uniform density distribution and consistent with total mass
          vol[xml_PixelModuleDeadAreaRight]->addMass(el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaRight]->getVolume()/deadAreaTotalVolume_mm3); 
          vol[xml_PixelModuleDeadAreaLeft]->addMass(el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaLeft]->getVolume()/deadAreaTotalVolume_mm3);   
          if(module.numSensors()==2){
            vol[xml_PixelModuleDeadAreaBackOfCentre]->addMass(el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaBackOfCentre]->getVolume()/deadAreaTotalVolume_mm3); 
            vol[xml_PixelModuleDeadAreaFrontOfCentre]->addMass(el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaFrontOfCentre]->getVolume()/deadAreaTotalVolume_mm3);   
          }
          vol[xml_PixelModuleDeadAreaFront]->addMass(el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaFront]->getVolume()/deadAreaTotalVolume_mm3);   
          vol[xml_PixelModuleDeadAreaBack]->addMass(el->quantityInGrams(module)*vol[xml_PixelModuleDeadAreaBack]->getVolume()/deadAreaTotalVolume_mm3);
  	}

      }

    }

    // OUTER TRACKER MODULE
    if (!module.isPixelModule()) {
      volumes.push_back(vol[xml_HybridFront]);
      volumes.push_back(vol[xml_HybridBack]);
      volumes.push_back(vol[xml_HybridLeft]);
      volumes.push_back(vol[xml_HybridRight]);
      volumes.push_back(vol[xml_HybridBetween]);
      volumes.push_back(vol[xml_SupportPlate]);
    }
    // PIXEL MODULE
    else {
      volumes.push_back(vol[xml_PixelModuleHybrid]);
      volumes.push_back(vol[xml_PixelModuleChip]);
      volumes.push_back(vol[xml_PixelModuleDeadAreaRight]);
     if(module.numSensors()==2){
        volumes.push_back(vol[xml_PixelModuleDeadAreaFrontOfCentre]);
        volumes.push_back(vol[xml_PixelModuleDeadAreaBackOfCentre]);
     }
      volumes.push_back(vol[xml_PixelModuleDeadAreaLeft]);
      volumes.push_back(vol[xml_PixelModuleDeadAreaFront]);
      volumes.push_back(vol[xml_PixelModuleDeadAreaBack]);
    }

    checkSubVolumes();
  }

  void ModuleComplex::checkSubVolumes() {
    for (auto& vit : volumes) {
      vit->check();
    }
  }

  void ModuleComplex::addShapeInfo(std::vector<ShapeInfo>& vec) {
    ShapeInfo ele;
    ele.type = bx; // Box 
    for (const auto& myVolume : volumes) {
      if (!myVolume->isValid()) continue;
      ele.name_tag = myVolume->getName();
      ele.dx = myVolume->getDx()/2.; // half length
      ele.dy = myVolume->getDy()/2.; // half length
      ele.dz = myVolume->getDz()/2.; // half length
      vec.push_back(ele);
    }
  }

  void ModuleComplex::addLogicInfo(std::vector<LogicalInfo>& vec) {
    LogicalInfo ele;
    for (const auto& myVolume : volumes) {
      if (!myVolume->isValid()) continue;
      ele.name_tag     = myVolume->getName();
      ele.shape_tag    = prefix_xmlfile + myVolume->getName(); 
      ele.material_tag = prefix_xmlfile + prefix_material + myVolume->getName();
      vec.push_back(ele);
    }
  }

  void ModuleComplex::addPositionInfo(std::vector<PosInfo>& vec) {
    PosInfo ele;
    ele.copy = 1;
    for (const auto& myVolume : volumes) {
      if (!myVolume->isValid()) continue;
      ele.parent_tag = prefix_xmlfile + myVolume->getParentName();
      ele.child_tag  = prefix_xmlfile + myVolume->getName();
      ele.trans.dx   = myVolume->getX();
      ele.trans.dy   = myVolume->getY();
      ele.trans.dz   = myVolume->getZ();
      vec.push_back(ele);
    }
  }

  void ModuleComplex::addMaterialInfo(std::vector<Composite>& vec) {
    for (const auto& myVolume : volumes) {
      if (!myVolume->isValid()) continue;

      Composite comp;
      comp.name    = prefix_material + myVolume->getName();
      comp.density = myVolume->getDensity();
      comp.method  = wt;

      double m = 0.0;
      for (const auto& it : myVolume->getMaterialList()) {
	comp.elements.insert(it);
	m += it.second;
      }
   
      for (auto& elem : comp.elements) {
        elem.second /= m;
      }

      comp.normalizedRIRatioPerMechanicalCategory = myVolume->getNormalizedRIRatioPerMechanicalCategory();
      
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
    double moduleTotalMass = 0.;
    for (const auto& myVolume : volumes) {
      myVolume->print();
      moduleTotalMass += myVolume->getMass();
    }
    std::cerr << "  Module Total Mass = " << moduleTotalMass 
              << " (" << moduleMassWithoutSensors_expected << " is expected.)" << std::endl; 
  }

}

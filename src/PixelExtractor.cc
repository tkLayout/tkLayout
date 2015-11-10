#include<PixelExtractor.h>
#include <ModuleCap.h>
#include<iomanip>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <string>
using boost::property_tree::ptree;
using boost::property_tree::write_xml;
using boost::property_tree::xml_writer_settings;
using boost::property_tree::xml_writer_make_settings;

namespace insur {

  void PixelExtractor::analyse( MaterialTable& mt,MaterialBudget& mb ) {

    std::string xml_bpix ="BPix_";
    std::string xml_lay ="Layer_";
    std::string xml_rod ="Rod_";

    Tracker& pix = mb.getTracker();
    InactiveSurfaces& inatvser = mb.getInactiveSurfaces();

    createPixelBarrelActive(pix);

    createPixelBarrelServices(inatvser); //ok

    createPixelEndcapActive(pix);

    createPixelEndcapServices(inatvser);//

    analyseElements(mt);

  }

  ////////////////////////////////////////////////////
  std::string PixelExtractor::whichShape( ShapeType s ) {
    if( s == ShapeType::bx )
      return "Box";
    else if( s == ShapeType::tb )
      return "Tubs";
    else if( s == ShapeType::tp )
      return "Trapezoid";
    else if( s == ShapeType::pc )
      return "Polycone";

  }
  //////////////////////////////////////////////////////
  void PixelExtractor::createPixelBarrelActive(Tracker& pix){
    LayerAggregator lagg;
    pix.accept(lagg);
    std::vector<Layer*>* bLayers = lagg.getBarrelLayers();
    numBarrelLayers = bLayers->size();
    std::vector<ShapeInfo> shapevec;

    std::vector<std::string> partselectors;

    std::vector<double> layerRadii;

    ModuleROCInfo minfo_zero={}; 
    SpecParInfo layer_spec,rod_spec;
    //Layer
    layer_spec.name = xml_subdet_layer + xml_par_tail;
    layer_spec.parameter.first = xml_tkddd_structure;
    layer_spec.parameter.second = xml_phaseII_pixbar + xml_layer;

    //Rod
    rod_spec.name = xml_phaseII_pixbar + xml_rod + xml_par_tail;
    rod_spec.parameter.first = xml_tkddd_structure;
    rod_spec.parameter.second = xml_phaseII_pixbar + xml_rod;

    //have to implement pixel mother volume!!!
    for( int i = 0; i<bLayers->size(); i++ ) {
      //Layer
      std::stringstream stemp;
      stemp << "Layer" << i+1;
      ShapeInfo layer_shape;
      layer_shape.type = tb;
      layer_shape.name_tag = stemp.str();
      layer_shape.rmin=bLayers->at(i)->minR();
      layer_shape.rmax=bLayers->at(i)->maxR();
      layer_shape.dz=bLayers->at(i)->maxZ();
      layer_shape.dx=0.;
      layer_shape.dy=0.;
      layer_shape.dxx=0.;
      layer_shape.dyy=0.;
      std::vector<std::pair<double, double> > temp;
      layer_shape.rzup=temp;
      layer_shape.rzdown=temp;
      cmsswXmlInfo.shapes.push_back(layer_shape);

      layerRadii.push_back(layer_shape.rmin);
      layerRadii.push_back(layer_shape.rmax);

      LogicalInfo layer_logic;
      layer_logic.name_tag = layer_shape.name_tag;
      layer_logic.shape_tag =  xml_phaseII_pixbar + ":" + layer_shape.name_tag;
      layer_logic.material_tag = "materials:Air";
      cmsswXmlInfo.logic.push_back(layer_logic);

      layer_spec.partselectors.push_back(layer_logic.name_tag);
      layer_spec.moduletypes.push_back(minfo_zero);

      PosInfo layer_pos;
      layer_pos.copy = 1;
      layer_pos.trans.dx = 0.0;
      layer_pos.trans.dy = 0.0;
      layer_pos.trans.dz = 0.0;
      layer_pos.parent_tag = xml_phaseII_pixbar;
      layer_pos.child_tag = layer_logic.shape_tag;
      cmsswXmlInfo.positions.push_back(layer_pos);

      //Rod
      const RodPair& r = bLayers->at(i)->rods().at(0);
      std::pair<const BmoduleContainer&,const BmoduleContainer&> mods = r.modules();
      auto& zplusmod = mods.first.at(0);
      double mod0Thick = zplusmod.sensorThickness() + zplusmod.hybridThickness();

      std::stringstream stemp_Rod;
      stemp_Rod << "Rod" << i+1;
      ShapeInfo rod_shape;
      rod_shape.type = bx;
      rod_shape.name_tag = stemp_Rod.str();
      //explicitly added module thickness twice to remind that chip thickness = module thickness
      rod_shape.dx = mods.first.at(0).thickness()/2. + mods.first.at(0).hybridThickness()/2. + + mods.first.at(0).thickness()/2.;
      rod_shape.dy = zplusmod.area() /zplusmod.length() / 2.0;//mod0Thick;
      rod_shape.dz = r.maxZ();
      rod_shape.rmin = 0.;
      rod_shape.rmax = 0.;
      cmsswXmlInfo.shapes.push_back(rod_shape);

      LogicalInfo rod_logic;
      rod_logic.name_tag = rod_shape.name_tag;
      rod_logic.shape_tag =  xml_phaseII_pixbar + ":" + rod_shape.name_tag;
      rod_logic.material_tag = "materials:Air";
      cmsswXmlInfo.logic.push_back(rod_logic);

      rod_spec.partselectors.push_back(rod_logic.name_tag);
      rod_spec.moduletypes.push_back(minfo_zero);
      //Rods will be placed with algorithm
      AlgoInfo rod_alg;

      rod_alg.name = xml_phialt_algo;
      rod_alg.parent = layer_shape.name_tag;
      std::stringstream par,parmap;

      rod_alg.parameter_map["ChildName"] = {rod_shape.name_tag,AlgoPartype::st};

      rod_alg.parameter_map["Tilt"] = {"90*deg",AlgoPartype::num};

      parmap << bLayers->at(i)->startAngle();
      rod_alg.parameter_map["StartAngle"] = {parmap.str(),AlgoPartype::num};
      parmap.str("");

      rod_alg.parameter_map["RangeAngle"] = {"360*deg",AlgoPartype::num};

      //not correct !chk extractor.cc//use rod minR,maxR

      parmap << layer_shape.rmin;
      rod_alg.parameter_map["RadiusIn"] = {parmap.str(),AlgoPartype::num};
      parmap.str("");

      parmap << layer_shape.rmax;
      rod_alg.parameter_map["RadiusOut"] = {parmap.str(),AlgoPartype::num};
      parmap.str("");

      rod_alg.parameter_map["ZPosition"] = {"0.0*mm",AlgoPartype::num};

      parmap << bLayers->at(i)->numRods();
      rod_alg.parameter_map["Number"] = {parmap.str(),AlgoPartype::num};
      parmap.str("");

      rod_alg.parameter_map["StartCopyNo"] = {"1",AlgoPartype::num};

      rod_alg.parameter_map["IncrCopyNo"] = {"1",AlgoPartype::num};

      rod_alg.parameters.push_back(par.str());
      cmsswXmlInfo.algos.push_back(rod_alg);
 
     double rod_rmin = bLayers->at(i)->minR();
      //double rod_dR = r.thickness();//
      //explicitly added module thickness twice to remind that chip thickness = module thickness
      //rod_dR = mod thickness + hybrid thickness + chip thickness(= mod thickness)
      double rod_dR = mods.first.at(0).thickness() + mods.first.at(0).hybridThickness() + mods.first.at(0).thickness();
      //creating a vector of BModules in sorted Z position(- to +)
      std::vector<BarrelModule> rod_modules;
      for( int j = mods.second.size() - 1 ; j >= 0 ; j--) {
        auto& mod = mods.second.at(j);
        rod_modules.push_back( mod );
      }

      for( unsigned int j =0; j < mods.first.size(); j++ ) {
        auto& mod = mods.first.at(j);
        rod_modules.push_back( mod );
      }

      if( r.myid() == 1) {
          std::stringstream barrelRmatpathCommon;
          barrelRmatpathCommon << "//" << xml_phaseII_pixbar << "/"
                               << layer_shape.name_tag << "/"
                               << rod_shape.name_tag << "/";
        getPixelBarrelModuleInfo( rod_modules ,i+1,rod_shape.name_tag,rod_rmin, rod_dR, barrelRmatpathCommon.str());
      }
    }

    //layer
    cmsswXmlInfo.specs.push_back(layer_spec);
    //rod
    cmsswXmlInfo.specs.push_back(rod_spec);
  }

  void PixelExtractor::getPixelBarrelModuleInfo( std::vector<BarrelModule> moduleVec ,int LayerNum,std::string rodName,
                                                 double rod_Rmin, double rod_dR, std::string barrelRmatpathCommon){

    std::vector<std::string> rmatPathvec;
    ShapeInfo modbox_shape,mod_shape,modwafer_shape,modactive_shape,hybrid_shape,chip_shape;
    modbox_shape.type = bx;
    modbox_shape.rmin = 0.;
    modbox_shape.rmax = 0.;

    mod_shape.type = bx;
    mod_shape.rmin = 0.;
    mod_shape.rmax = 0.;

    modwafer_shape.type = bx;
    modwafer_shape.rmin = 0.;
    modwafer_shape.rmax = 0.;


    modactive_shape.type = bx;
    modactive_shape.rmin = 0.;
    modactive_shape.rmax = 0.;

    hybrid_shape.type = bx;
    hybrid_shape.rmin = 0.;
    hybrid_shape.rmax = 0.;

    chip_shape.type = bx;
    chip_shape.rmin = 0.;
    chip_shape.rmax = 0.;

    LogicalInfo modbox_logic,mod_logic,modwafer_logic,modactive_logic,hybrid_logic,chip_logic;

    PosInfo modbox_pos,mod_pos,modwafer_pos,modactive_pos,hybrid_pos,chip_pos;

    modbox_pos.parent_tag = rodName;
    //mod_pos.parent_tag = rodName;
    //hybrid_pos.parent_tag = rodName;


    ModuleROCInfo minfo;
    ModuleROCInfo minfo_zero={}; 
    //Module
    SpecParInfo module_spec;
    module_spec.parameter.first = xml_tkddd_structure;
    module_spec.parameter.second = xml_phaseII_pixbardet;

    for( unsigned int i = 0; i<moduleVec.size(); i++ ) {
      auto& module = moduleVec.at(i);

      std::stringstream stemp_mod;

      module_spec.partselectors.clear();
      module_spec.moduletypes.clear();

      //Module part
      stemp_mod << xml_bmodbox << i+1 << xml_layer << LayerNum ;

      modbox_shape.name_tag = stemp_mod.str();
      modbox_shape.dx = module.area() /module.length() / 2.0;
      modbox_shape.dy = module.length() / 2.0;//mod0Thick;
      //explicitly added module thickness twice to remind that chip thickness = module thickness
      modbox_shape.dz = module.thickness()/ 2.0 + module.hybridThickness()/2. + module.thickness()/ 2.0;


      stemp_mod.str("");
      stemp_mod << xml_barrel_module << i+1 << xml_layer << LayerNum << xml_phaseII_pixeldetTag;
      mod_shape.name_tag = stemp_mod.str();
      mod_shape.dx = module.area() /module.length() / 2.0;
      mod_shape.dy = module.length() / 2.0;//mod0Thick;
      mod_shape.dz = module.thickness()/ 2.0;


      hybrid_shape.name_tag = mod_shape.name_tag + xml_phaseII_pixelHybridTag;
      hybrid_shape.dx = mod_shape.dx;
      hybrid_shape.dy = mod_shape.dy;//mod0Thick;
      hybrid_shape.dz = module.hybridThickness()/2.;

      chip_shape.name_tag = mod_shape.name_tag + xml_phaseII_pixelChipTag;
      chip_shape.dx = module.area() /module.length() / 2.0;
      chip_shape.dy = module.length() / 2.0;//mod0Thick;
      chip_shape.dz = module.thickness()/ 2.0;


      modbox_logic.name_tag = modbox_shape.name_tag;
      modbox_logic.shape_tag =  xml_phaseII_pixbar + ":" + modbox_shape.name_tag;
      modbox_logic.material_tag = "materials:Air";

      mod_logic.name_tag = mod_shape.name_tag;
      mod_logic.shape_tag =  xml_phaseII_pixbar + ":" + mod_shape.name_tag;
      mod_logic.material_tag = "materials:Air";


      hybrid_logic.name_tag = hybrid_shape.name_tag;
      hybrid_logic.shape_tag =  xml_phaseII_pixbar + ":" + hybrid_shape.name_tag;
      stringstream hybrid_logic_mat;
      hybrid_logic_mat << "pixel:topInactiveComposite" << xml_barrel_module << i+1 << xml_layer << LayerNum ;
      hybrid_logic.material_tag = hybrid_logic_mat.str();

      chip_logic.name_tag = chip_shape.name_tag;
      chip_logic.shape_tag =  xml_phaseII_pixbar + ":" + chip_shape.name_tag;
      stringstream chip_logic_mat;
      chip_logic_mat << "pixel:bottomInactiveComposite" << xml_barrel_module << i+1 << xml_layer << LayerNum ;
      chip_logic.material_tag = chip_logic_mat.str();


      modbox_pos.copy = 1;
      modbox_pos.trans.dx = 0.;
      modbox_pos.trans.dy = 0.;
      modbox_pos.trans.dz = module.maxZ() - mod_shape.dy;
      modbox_pos.rotref = xml_places_unflipped_mod_in_rod;
      modbox_pos.child_tag = modbox_logic.shape_tag;

//Volume positioning inside modulebox
//     Hybrid --> TopInactive
//     Module --> Active
//     Chip   --> BottomInactive
/////////////////////////////////////
     
      chip_pos.copy = 1;
      chip_pos.trans.dx = 0.;
      chip_pos.trans.dy = 0.;
      chip_pos.trans.dz = -modbox_shape.dz + chip_shape.dz;
      chip_pos.child_tag = chip_logic.shape_tag;
      chip_pos.parent_tag = modbox_shape.name_tag;

      mod_pos.copy = 1;
      //  ss << "ModRho=" << zplusmod.center().Rho() << " rod_rmin + rod_dR/2.=" << rod_rmin + rod_dR/2.  << std::endl;
      //if( zplusmod.center().Rho() > (rod_rmin + rod_dR/2.) )
      //  mod_pos.trans.dx = rod_dR / 2.0 - mod_shape.dz;
      //else
      //mod_pos.trans.dx = mod_shape.dz - rod_dR / 2.0;
      mod_pos.trans.dx = 0.;
      mod_pos.trans.dy = 0.0;
      mod_pos.trans.dz = chip_pos.trans.dz + chip_shape.dz + mod_shape.dz;
      mod_pos.child_tag = mod_logic.shape_tag;
      mod_pos.parent_tag = modbox_shape.name_tag;
 
      hybrid_pos.copy = 1;
      hybrid_pos.trans.dx = 0.;
      hybrid_pos.trans.dy = 0.;
      hybrid_pos.trans.dz = mod_pos.trans.dz + mod_shape.dz + hybrid_shape.dz;
      hybrid_pos.child_tag = hybrid_logic.shape_tag;
      hybrid_pos.parent_tag = modbox_shape.name_tag;

      if( module.flipped() ) modbox_pos.rotref = xml_places_flipped_mod_in_rod;

      cmsswXmlInfo.shapes.push_back(modbox_shape);
      cmsswXmlInfo.logic.push_back(modbox_logic);
      cmsswXmlInfo.positions.push_back(modbox_pos);

      cmsswXmlInfo.shapes.push_back(mod_shape);
      cmsswXmlInfo.logic.push_back(mod_logic);
      cmsswXmlInfo.positions.push_back(mod_pos);

      cmsswXmlInfo.shapes.push_back(hybrid_shape);
      cmsswXmlInfo.logic.push_back(hybrid_logic);
      cmsswXmlInfo.positions.push_back(hybrid_pos);

      cmsswXmlInfo.shapes.push_back(chip_shape);
      cmsswXmlInfo.logic.push_back(chip_logic);
      cmsswXmlInfo.positions.push_back(chip_pos);

      //////////////////active part//////////////////////////

      modwafer_shape.name_tag = mod_shape.name_tag + "Wafer";
      modwafer_shape.dx = mod_shape.dx;
      modwafer_shape.dy = mod_shape.dy;
      modwafer_shape.dz = module.sensors().front().sensorThickness()/2.;


      modwafer_logic.name_tag = modwafer_shape.name_tag;
      modwafer_logic.shape_tag =  xml_phaseII_pixbar + ":" + modwafer_shape.name_tag;
      modwafer_logic.material_tag = "pixel:SenSi";//??or air

      modwafer_pos.copy = 1;
      modwafer_pos.trans.dx = 0.;
      modwafer_pos.trans.dy = 0.;
      modwafer_pos.trans.dz = 0.;
      modwafer_pos.parent_tag = mod_pos.child_tag;
      modwafer_pos.child_tag = modwafer_logic.shape_tag;

      //inactive space around sensor boundary
      double cut_dx = 0., cut_dy = 0., cut_dz = 0.;

      modactive_shape.name_tag = mod_shape.name_tag + "Active";
      modactive_shape.dx = modwafer_shape.dx - cut_dx;
      modactive_shape.dy = modwafer_shape.dy - cut_dy;
      modactive_shape.dz = modwafer_shape.dz - cut_dz;


      modactive_logic.name_tag = modactive_shape.name_tag;
      modactive_logic.shape_tag =  xml_phaseII_pixbar + ":" + modactive_shape.name_tag;
      modactive_logic.material_tag = "pixel:SenSi";

      modactive_pos.copy = 1;
      modactive_pos.trans.dx = 0.;
      modactive_pos.trans.dy = 0.;
      modactive_pos.trans.dz = 0.;
      modactive_pos.parent_tag = modwafer_pos.child_tag;
      modactive_pos.child_tag = modactive_logic.shape_tag;

      cmsswXmlInfo.shapes.push_back(modwafer_shape);
      cmsswXmlInfo.logic.push_back(modwafer_logic);        
      cmsswXmlInfo.positions.push_back(modwafer_pos);

      cmsswXmlInfo.shapes.push_back(modactive_shape);
      cmsswXmlInfo.logic.push_back(modactive_logic);

      module_spec.name=modactive_logic.name_tag+"Par";
      module_spec.partselectors.push_back(modactive_logic.name_tag);
      minfo.name		= module.moduleType();
      minfo.rocrows	= any2str<int>(module.innerSensor().numROCRows());  // in case of single sensor module innerSensor() and outerSensor() point to the same sensor
      minfo.roccols	= any2str<int>(module.innerSensor().numROCCols());
      minfo.rocx		= any2str<int>(module.innerSensor().numROCX());
      minfo.rocy		= any2str<int>(module.innerSensor().numROCY());
      module_spec.moduletypes.push_back(minfo);

      cmsswXmlInfo.positions.push_back(modactive_pos);
        //analyseCompositeElements( hybrid_logic.material_tag, compositeDensity(*module.getModuleCap(), true),
        //  *module.getModuleCap(), true, ss);
      addMaterialToVolume( hybrid_logic.material_tag, module, volumeType::hybrid,
                           (8.*hybrid_shape.dx*hybrid_shape.dy*hybrid_shape.dz/1000.));
      addMaterialToVolume( chip_logic.material_tag, module, volumeType::chip,
                           (8.*chip_shape.dx*chip_shape.dy*chip_shape.dz/1000.));
      cmsswXmlInfo.specs.push_back(module_spec);
              
      barrelRmatpath.push_back( barrelRmatpathCommon + modbox_shape.name_tag +"/" 
                                + mod_shape.name_tag + "/" + modwafer_shape.name_tag + "/" + 
                                modactive_shape.name_tag );
    }
  }
  void PixelExtractor::createPixelBarrelServices(InactiveSurfaces& inatcvser) {
    //Inactive elements part
    //InactiveSurfaces& inatvser = pix.getInactiveSurfaces();
    std::vector<InactiveElement>& bser = inatcvser.getBarrelServices();
    ShapeInfo bser_shape;
    bser_shape.type = tb;
    bser_shape.dx = 0.0;
    bser_shape.dy = 0.0;
    bser_shape.dyy = 0.0;

    LogicalInfo bser_logic;

    PosInfo bser_pos;
    bser_pos.copy = 1;
    bser_pos.trans.dx = 0.0;
    bser_pos.trans.dy = 0.0;

    Composite comp;

    comp.method = wt;

    for( auto& bs: bser ) {

      std::stringstream sname,matname;
      sname << xml_phaseII_pixbar << xml_base_serf << "R" << (int)(bs.getInnerRadius()) << "Z" << (int)(bs.getZOffset());
      bser_shape.name_tag = sname.str();
      bser_shape.rmin = bs.getInnerRadius();
      bser_shape.rmax = bs.getInnerRadius() + bs.getRWidth();
      bser_shape.dz = bs.getZLength() / 2.0;

      cmsswXmlInfo.shapes.push_back(bser_shape);

      matname << xml_base_serfcomp
        << "R" << (int)(bs.getInnerRadius())
        << "Z" << (int)(bs.getZLength());
      bser_logic.name_tag = sname.str();
      bser_logic.shape_tag =  xml_phaseII_pixbar + ":" + sname.str();
      bser_logic.material_tag = xml_phaseII_pixbar + ":" + matname.str();

      cmsswXmlInfo.logic.push_back(bser_logic);

      bser_pos.parent_tag = xml_phaseII_pixbar;
      bser_pos.child_tag = bser_logic.shape_tag;
      bser_pos.trans.dz = bs.getZOffset() + bser_shape.dz;

      cmsswXmlInfo.positions.push_back(bser_pos);

      std::stringstream bsermod;

      comp.name = bser_logic.material_tag;

      comp.density = compositeDensity(bs);
      analyseCompositeElements( bser_logic.material_tag, comp.density,bs, false);
    }
 }
 
 //Functions for endcap
  void PixelExtractor::createPixelEndcapActive(Tracker& pix) {
    LayerAggregator lagg;
    pix.accept(lagg);
    std::vector<Disk*>* eDisks = lagg.getEndcapDisks();
    
    std::vector<std::string> partselectors;

    ModuleROCInfo minfo_zero={}; 
    SpecParInfo disc_spec, ring_spec, module_spec;
    // Disk
    disc_spec.name = xml_subdet_wheel + xml_par_tail;
    disc_spec.parameter.first = xml_tkddd_structure;
    disc_spec.parameter.second = xml_phaseII_pixecap + xml_disc;
    // Ring
    ring_spec.name = xml_subdet_ring + xml_par_tail;
    ring_spec.parameter.first = xml_tkddd_structure;
    ring_spec.parameter.second = xml_phaseII_pixecap + "Panel";

    //have to implement pixel mother volume!!!
    for( int i = 0; i<eDisks->size(); i++ ) {
      if( ( eDisks->at(i)->maxZ() + eDisks->at(i)->minZ() ) / 2.0 < 0.) continue;
      //Layer
      std::stringstream stemp;
      stemp << "Disc" << i+1;
      ShapeInfo disc_shape;

      disc_shape.type = tb;
      disc_shape.name_tag = stemp.str();
      disc_shape.rmin = eDisks->at(i)->minR();
      disc_shape.rmax = eDisks->at(i)->maxR();
      disc_shape.dz = eDisks->at(i)->thickness()/2.;


      LogicalInfo disc_logic; 
      disc_logic.name_tag = disc_shape.name_tag;
      disc_logic.shape_tag =  xml_phaseII_pixecap + ":" + disc_shape.name_tag;
      disc_logic.material_tag = "materials:Air";

      PosInfo disc_pos;
      disc_pos.copy = 1;
      disc_pos.trans.dx = 0.0;
      disc_pos.trans.dy = 0.0;
      //will the mod+hybrid fit inside the disc thickness?check!!
      disc_pos.trans.dz = ( eDisks->at(i)->maxZ() + eDisks->at(i)->minZ() ) / 2.0;// - xml_z_pixfwd;
      disc_pos.parent_tag = xml_phaseII_pixecap;
      disc_pos.child_tag = disc_logic.shape_tag;

      cmsswXmlInfo.shapes.push_back(disc_shape);
      cmsswXmlInfo.logic.push_back(disc_logic);

      disc_spec.partselectors.push_back(disc_logic.name_tag);
      disc_spec.moduletypes.push_back(minfo_zero);

      cmsswXmlInfo.positions.push_back(disc_pos);

      //Rings
      const EndcapRings ringVec = eDisks->at(i)->rings();

      ShapeInfo ring_shape;
      ring_shape.type = tb;

      ring_shape.dx = 0.0;
      ring_shape.dy = 0.0;
      ring_shape.dyy = 0.0;


      LogicalInfo ring_logic;
      ring_logic.material_tag = "materials:Air";
      discRingpair.push_back({i+1,ringVec.size()});
      for( unsigned int ri = 0; ri<ringVec.size(); ri++) {
        const EmoduleContainer& emodules = ringVec.at(ri).modules();

        std::stringstream rname;
        rname << "Ring" << ri+1 << disc_shape.name_tag;
        ring_shape.name_tag = rname.str();
        //explicitly added mod thickness twice to remind that chip thickness= mod thickness
        ring_shape.dz = eDisks->at(i)->maxRingThickness() / 2.0 + emodules.at(0).hybridThickness()/2. + eDisks->at(i)->maxRingThickness() / 2.0;
        ring_shape.rmin = emodules.at(0).minR();//minR of ring = minR of mod
        ring_shape.rmax = emodules.at(0).maxR();//maxR for mod > maxR for ring. But ring has to house the mod.So!!


        cmsswXmlInfo.shapes.push_back(ring_shape);

        ring_logic.name_tag = ring_shape.name_tag;
        ring_logic.shape_tag = xml_phaseII_pixecap + ":" + ring_shape.name_tag;
        cmsswXmlInfo.logic.push_back(ring_logic); 

        ring_spec.partselectors.push_back(ring_logic.name_tag);
        ring_spec.moduletypes.push_back(minfo_zero);

        PosInfo ring_pos;
        ring_pos.copy = 1;
        ring_pos.parent_tag = xml_phaseII_pixecap + ":" + disc_shape.name_tag;
        ring_pos.child_tag = xml_phaseII_pixecap + ":" + ring_logic.name_tag;
        ring_pos.trans.dx = 0.;
        ring_pos.trans.dy = 0.;
        if( emodules.at(0).center().Z() < (eDisks->at(i)->maxZ() + eDisks->at(i)->minZ()) / 2.0 )
          ring_pos.trans.dz = (eDisks->at(i)->minZ() - eDisks->at(i)->maxZ()) / 2.0 + ring_shape.dz;
        else
          ring_pos.trans.dz = (eDisks->at(i)->maxZ() - eDisks->at(i)->minZ()) / 2.0 - ring_shape.dz;
        cmsswXmlInfo.positions.push_back(ring_pos);

        //Endcap Module Info
        std::stringstream ecapRmatpathCommon;
          ecapRmatpathCommon << "//" << xml_phaseII_pixecap << "/"
          << disc_shape.name_tag << "/"
          << ring_shape.name_tag << "/";
        getPixelEndcapModuleInfo( emodules.at(0) ,ri+1, i+1,ring_shape.name_tag, ring_shape.dz, 0.,ecapRmatpathCommon.str());

        AlgoInfo ring_algo;
        ring_algo.name = xml_trackerring_algo;//track:DDTrackerRingAlgo
        ring_algo.parent = ring_logic.shape_tag;
        std::stringstream emodname,algopar;
        emodname << xml_emodbox << ri+1 << "Disc" << i+1;

        algopar << xml_phaseII_pixecap + ":" + emodname.str();
        ring_algo.parameter_map[xml_childparam] = {algopar.str(),AlgoPartype::st};

        algopar.str("");
        algopar << ringVec.at(ri).nModules()/2;
        ring_algo.parameter_map[xml_nmods]={algopar.str(),AlgoPartype::num};

        ring_algo.parameter_map[xml_startcopyno]={"1",AlgoPartype::num};

        ring_algo.parameter_map[xml_incrcopyno]={"2",AlgoPartype::num};

        ring_algo.parameter_map[xml_rangeangle]={"360*deg",AlgoPartype::num};        

        algopar.str("");
        algopar << emodules.at(0).center().Phi();
        ring_algo.parameter_map[xml_startangle]={algopar.str(),AlgoPartype::num};

        algopar.str("");
        algopar << emodules.at(0).center().Rho();
        ring_algo.parameter_map[xml_radius]={algopar.str(),AlgoPartype::num};

        algopar.str("");
        ring_algo.vecpar.name="Center";
        ring_algo.vecpar.type="Numeric";
        ring_algo.vecpar.nEntries="3";
        ring_algo.vecpar.values.push_back(0);
        ring_algo.vecpar.values.push_back(0);
        ring_algo.vecpar.values.push_back(ring_shape.dz - emodules.at(0).thickness() / 2.0); 

        ring_algo.parameter_map[xml_iszplus]={"1",AlgoPartype::num};
        ring_algo.parameter_map[xml_tiltangle]={"90*deg",AlgoPartype::num};
        ring_algo.parameter_map[xml_isflipped]={"0",AlgoPartype::num};
        //push the first child module into algorithm
        //this module is not flipped
        //the copy numbers for these modules will be odd
        cmsswXmlInfo.algos.push_back(ring_algo);

        ring_algo.parameters.clear();

        algopar<< xml_phaseII_pixecap + ":" + emodname.str();
        ring_algo.parameter_map[xml_childparam] = {algopar.str(),AlgoPartype::st};     

        algopar.str(""); 
        algopar << ringVec.at(ri).nModules()/2;
        ring_algo.parameter_map[xml_nmods]={algopar.str(),AlgoPartype::num};

        ring_algo.parameter_map[xml_startcopyno]={"2",AlgoPartype::num};

        ring_algo.parameter_map[xml_incrcopyno]={"2",AlgoPartype::num};

        ring_algo.parameter_map[xml_rangeangle]={"360*deg",AlgoPartype::num};

        algopar.str("");
        algopar << (emodules.at(0).center().Phi() + 2 * M_PI / (double)(ringVec.at(ri).nModules()));

        ring_algo.parameter_map[xml_startangle]={algopar.str(),AlgoPartype::num};

        algopar.str("");
        algopar << emodules.at(0).center().Rho();

        ring_algo.parameter_map[xml_radius]={algopar.str(),AlgoPartype::num};

        algopar.str("");
        ring_algo.vecpar.values.clear();
        ring_algo.vecpar.values.push_back(0);
        ring_algo.vecpar.values.push_back(0);
        ring_algo.vecpar.values.push_back(emodules.at(0).thickness() / 2.0 - ring_shape.dz); 
        //push the flipped child module
        //copy numbers are even
        ring_algo.parameter_map[xml_iszplus]={"1",AlgoPartype::num};
        ring_algo.parameter_map[xml_tiltangle]={"90*deg",AlgoPartype::num};
        ring_algo.parameter_map[xml_isflipped]={"1",AlgoPartype::num};
  
        cmsswXmlInfo.algos.push_back(ring_algo);
      }
      
    }
    //disc
    cmsswXmlInfo.specs.push_back(disc_spec);

    //ring
    cmsswXmlInfo.specs.push_back(ring_spec);
  }

    void PixelExtractor::getPixelEndcapModuleInfo(EndcapModule emodule ,int ringNo, int discno,
                                                  std::string ringName, double ring_thickness, double rod_dR,
                                                  std::string ecapRmatpathCommon) {

    ShapeInfo emodbox_shape,emod_shape,emodwafer_shape,emodactive_shape,ehybrid_shape,chip_shape;// , flip_shape;

    LogicalInfo emodbox_logic,emod_logic,emodwafer_logic,emodactive_logic,ehybrid_logic,chip_logic; //,flip_logic;

    PosInfo emodbox_pos,emod_pos,emodwafer_pos,emodactive_pos,ehybrid_pos,chip_pos; //,flip_pos;

    ModuleROCInfo minfo;
    //Module
    SpecParInfo module_spec;//, flip_spec;
    module_spec.parameter.first = xml_tkddd_structure;
    module_spec.parameter.second = xml_phaseII_pixecapdet;

    std::stringstream emodbox_name;
    emodbox_name << xml_emodbox << ringNo << "Disc" << discno;
    emodbox_shape.name_tag = emodbox_name.str();
    emodbox_shape.type = emodule.shape() == RECTANGULAR ? bx : tp;
    emodbox_shape.dx = emodule.minWidth() / 2.0;
    emodbox_shape.dxx = emodule.maxWidth() / 2.0;
    emodbox_shape.dy = emodule.length() / 2.0;
    emodbox_shape.dyy = emodule.length() / 2.0;
    //chip thickness = module thickness
    emodbox_shape.dz = emodule.thickness() / 2.0 + emodule.hybridThickness()/2. + emodule.thickness() / 2.0;
    emodbox_shape.rmin = 0.;
    emodbox_shape.rmax = 0.;
    cmsswXmlInfo.shapes.push_back(emodbox_shape);//no flip module
    //flip_shape  = emodbox_shape;
    //flip_shape.name_tag = emodbox_shape.name_tag + "FLIPPED";
    //cmsswXmlInfo.shapes.push_back(flip_shape);//flip module
    

    std::stringstream e_name;
    e_name << xml_endcap_module << ringNo << "Disc" << discno << xml_phaseII_pixeldetTag;
    emod_shape.name_tag = e_name.str();
    emod_shape.type = emodbox_shape.type;
    emod_shape.dx = emodbox_shape.dx;
    emod_shape.dxx = emodbox_shape.dxx;
    emod_shape.dy = emodbox_shape.dy;
    emod_shape.dyy = emodbox_shape.dyy;
    emod_shape.dz = emodule.thickness() / 2.0;
    emod_shape.rmin = 0.;
    emod_shape.rmax = 0.;
    cmsswXmlInfo.shapes.push_back(emod_shape);
    //flip_shape  = emod_shape;
    //flip_shape.name_tag = emod_shape.name_tag + "FLIPPED";
    //cmsswXmlInfo.shapes.push_back(flip_shape);//flip module
    

    ehybrid_shape.name_tag = emod_shape.name_tag + xml_phaseII_pixelHybridTag;
    ehybrid_shape.type = emod_shape.type;
    ehybrid_shape.dx = emod_shape.dx;
    ehybrid_shape.dxx = emod_shape.dxx;
    ehybrid_shape.dy = emod_shape.dy;
    ehybrid_shape.dyy = emod_shape.dyy;
    ehybrid_shape.dz = emodule.hybridThickness()/2.;
    cmsswXmlInfo.shapes.push_back(ehybrid_shape);
    //flip_shape  = ehybrid_shape;
    //flip_shape.name_tag = ehybrid_shape.name_tag + "FLIPPED";
    //cmsswXmlInfo.shapes.push_back(flip_shape);//flip module

    chip_shape.name_tag = emod_shape.name_tag + xml_phaseII_pixelChipTag;
    chip_shape.type = emodbox_shape.type;
    chip_shape.dx = emodbox_shape.dx;
    chip_shape.dxx = emodbox_shape.dxx;
    chip_shape.dy = emodbox_shape.dy;
    chip_shape.dyy = emodbox_shape.dyy;
    chip_shape.dz = emodule.thickness() / 2.0;
    chip_shape.rmin = 0.;
    chip_shape.rmax = 0.;
    cmsswXmlInfo.shapes.push_back(chip_shape);
    //flip_shape  = chip_shape;
    //flip_shape.name_tag = chip_shape.name_tag + "FLIPPED";
    //cmsswXmlInfo.shapes.push_back(flip_shape);//flip module


    //*************Logical Volumes**********//
    emodbox_logic.name_tag = emodbox_shape.name_tag;
    emodbox_logic.shape_tag = xml_phaseII_pixecap + ":" + emodbox_shape.name_tag;
    emodbox_logic.material_tag = "materials:Air";
    cmsswXmlInfo.logic.push_back(emodbox_logic);
    //flip_logic = emodbox_logic;
    //flip_logic.name_tag = emodbox_logic.name_tag + "FLIPPED";
    //flip_logic.shape_tag = emodbox_logic.shape_tag + "FLIPPED";
    //cmsswXmlInfo.logic.push_back(flip_logic);
     
    emod_logic.name_tag = emod_shape.name_tag;
    emod_logic.shape_tag =  xml_phaseII_pixecap + ":" + emod_shape.name_tag;
    emod_logic.material_tag = "materials:Air";
    cmsswXmlInfo.logic.push_back(emod_logic);
    //flip_logic = emod_logic;
    //flip_logic.name_tag = emod_logic.name_tag + "FLIPPED";
    //flip_logic.shape_tag = emod_logic.shape_tag + "FLIPPED";
    //cmsswXmlInfo.logic.push_back(flip_logic);

    ehybrid_logic.name_tag = ehybrid_shape.name_tag;
    ehybrid_logic.shape_tag =  xml_phaseII_pixecap + ":" + ehybrid_shape.name_tag;
    stringstream ehybrid_logic_mat;
    ehybrid_logic_mat << "pixel:topInactiveCompositeEModule" << ringNo << "Disc" << discno;
    ehybrid_logic.material_tag = ehybrid_logic_mat.str();
    cmsswXmlInfo.logic.push_back(ehybrid_logic);
    //flip_logic = ehybrid_logic;
    //flip_logic.name_tag = ehybrid_logic.name_tag + "FLIPPED";
    //flip_logic.shape_tag = ehybrid_logic.shape_tag + "FLIPPED";
    //cmsswXmlInfo.logic.push_back(flip_logic);

    chip_logic.name_tag = chip_shape.name_tag;
    chip_logic.shape_tag =  xml_phaseII_pixecap + ":" + chip_shape.name_tag;
    stringstream chip_logic_mat;
    chip_logic_mat << "pixe:bottomInactiveCompositeEModule" << ringNo << "Disc" << discno;
    chip_logic.material_tag = chip_logic_mat.str();
    cmsswXmlInfo.logic.push_back(chip_logic);
    //flip_logic = chip_logic;
    //flip_logic.name_tag = chip_logic.name_tag + "FLIPPED";
    //flip_logic.shape_tag = chip_logic.shape_tag + "FLIPPED";
    //cmsswXmlInfo.logic.push_back(flip_logic);


    chip_pos.copy = 1;
    chip_pos.parent_tag = emodbox_shape.name_tag;
    chip_pos.trans.dx = 0.;
    chip_pos.trans.dy = 0.;
    chip_pos.trans.dz = -emodbox_shape.dz + chip_shape.dz;
    chip_pos.child_tag = chip_logic.shape_tag;
    cmsswXmlInfo.positions.push_back(chip_pos);
    //flip_pos = chip_pos;
    //flip_pos.parent_tag += "FLIPPED";
    //flip_pos.child_tag += "FLIPPED";
    //flip_pos.trans.dz *= -1.;
    //cmsswXmlInfo.positions.push_back(flip_pos);

    emod_pos.copy = 1;
    emod_pos.parent_tag = emodbox_shape.name_tag;
    emod_pos.trans.dx = 0.;
    emod_pos.trans.dy = 0.;
    emod_pos.trans.dz = chip_pos.trans.dz + chip_shape.dz + emod_shape.dz;
    emod_pos.child_tag = emod_logic.shape_tag;
    cmsswXmlInfo.positions.push_back(emod_pos);
    //flip_pos = emod_pos;
    //flip_pos.parent_tag += "FLIPPED";
    //flip_pos.child_tag += "FLIPPED";
    //flip_pos.trans.dz *= -1.;
    //cmsswXmlInfo.positions.push_back(flip_pos);

    ehybrid_pos.copy = 1;
    ehybrid_pos.parent_tag = emodbox_shape.name_tag;
    ehybrid_pos.trans.dx = 0.;
    ehybrid_pos.trans.dy = 0.;
    ehybrid_pos.trans.dz = emod_pos.trans.dz + emod_shape.dz + ehybrid_shape.dz;
    ehybrid_pos.child_tag = ehybrid_logic.shape_tag;
    cmsswXmlInfo.positions.push_back(ehybrid_pos);
    //flip_pos = ehybrid_pos;
    //flip_pos.parent_tag += "FLIPPED";
    //flip_pos.child_tag += "FLIPPED";
    //flip_pos.trans.dz *= -1.;
    //cmsswXmlInfo.positions.push_back(flip_pos);


    //////////////////active part//////////////////////////
    //inactive space around sensor boundary
    double cut_dx = 0., cut_dy = 0., cut_dz = 0.;

    emodwafer_shape.type = emod_shape.type;
    emodwafer_shape.rmin = 0.;
    emodwafer_shape.rmax = 0.;
    emodwafer_shape.name_tag = emod_shape.name_tag + "Wafer";
    emodwafer_shape.dx = emod_shape.dx;
    emodwafer_shape.dx = emod_shape.dxx;
    emodwafer_shape.dy = emod_shape.dy;
    emodwafer_shape.dyy = emod_shape.dyy;
    emodwafer_shape.dz = emodule.sensors().front().sensorThickness()/2.;
    cmsswXmlInfo.shapes.push_back(emodwafer_shape);
    //flip_shape  = emodwafer_shape;
    //flip_shape.name_tag = emodwafer_shape.name_tag + "FLIPPED";
    //cmsswXmlInfo.shapes.push_back(flip_shape);//flip module

    emodactive_shape.type = emod_shape.type;
    emodactive_shape.rmin = 0.;
    emodactive_shape.rmax = 0.;
    emodactive_shape.name_tag = emod_shape.name_tag + "Active";
    emodactive_shape.dx = emodwafer_shape.dx - cut_dx;
    emodactive_shape.dy = emodwafer_shape.dy - cut_dy;
    emodactive_shape.dz = emodwafer_shape.dz - cut_dz;
    cmsswXmlInfo.shapes.push_back(emodwafer_shape);
    //flip_shape  = emodactive_shape;
    //flip_shape.name_tag = emodactive_shape.name_tag + "FLIPPED";
    //cmsswXmlInfo.shapes.push_back(flip_shape);//flip module

    emodwafer_logic.name_tag = emodwafer_shape.name_tag;
    emodwafer_logic.shape_tag =  xml_phaseII_pixecap + ":" + emodwafer_shape.name_tag;
    emodwafer_logic.material_tag = "pixel:SenSi";//??or air
    cmsswXmlInfo.logic.push_back(emodwafer_logic);
    //flip_logic = emodwafer_logic;
    //flip_logic.name_tag = emodwafer_logic.name_tag + "FLIPPED";
    //flip_logic.shape_tag = emodwafer_logic.shape_tag + "FLIPPED";
    //cmsswXmlInfo.logic.push_back(flip_logic);


    emodactive_logic.name_tag = emodactive_shape.name_tag;
    emodactive_logic.shape_tag =  xml_phaseII_pixecap + ":" + emodactive_shape.name_tag;
    emodactive_logic.material_tag = "pixel:SenSi";
    cmsswXmlInfo.logic.push_back(emodactive_logic);
    //flip_logic = emodactive_logic;
    //flip_logic.name_tag = emodactive_logic.name_tag + "FLIPPED";
    //flip_logic.shape_tag = emodactive_logic.shape_tag + "FLIPPED";
    //cmsswXmlInfo.logic.push_back(flip_logic);

    emodwafer_pos.copy = 1;
    emodwafer_pos.trans.dx = 0.;
    emodwafer_pos.trans.dy = 0.;
    emodwafer_pos.trans.dz = 0.;
    emodwafer_pos.parent_tag = emod_pos.child_tag;
    emodwafer_pos.child_tag = emodwafer_logic.shape_tag;
    cmsswXmlInfo.positions.push_back(emod_pos);
    //flip_pos = emodwafer_pos;
    //flip_pos.parent_tag += "FLIPPED";
    //flip_pos.child_tag += "FLIPPED";
    //flip_pos.trans.dz *= -1.;
    //cmsswXmlInfo.positions.push_back(flip_pos);


    emodactive_pos.copy = 1;
    emodactive_pos.trans.dx = 0.;
    emodactive_pos.trans.dy = 0.;
    emodactive_pos.trans.dz = 0.;
    emodactive_pos.parent_tag = emodwafer_pos.child_tag;
    emodactive_pos.child_tag = emodactive_logic.shape_tag;
    cmsswXmlInfo.positions.push_back(emodactive_pos);
    //flip_pos = emodactive_pos;
    //flip_pos.parent_tag += "FLIPPED";
    //flip_pos.child_tag += "FLIPPED";
    //flip_pos.trans.dz *= -1.;
    //cmsswXmlInfo.positions.push_back(flip_pos);

    module_spec.name = emodactive_logic.name_tag + "Par";
    module_spec.partselectors.push_back(emodactive_logic.name_tag);
    minfo.name		= emodule.moduleType();
    minfo.rocrows	= any2str<int>(emodule.innerSensor().numROCRows());  // in case of single sensor module innerSensor() and outerSensor() point to the same sensor
    minfo.roccols	= any2str<int>(emodule.innerSensor().numROCCols());
    minfo.rocx		= any2str<int>(emodule.innerSensor().numROCX());
    minfo.rocy		= any2str<int>(emodule.innerSensor().numROCY());
    module_spec.moduletypes.push_back(minfo);
    //flip_spec = module_spec;
    //flip_spec.partselectors.clear();
    //flip_spec.name = emodactive_logic.name_tag + "FLIPPEDPar";
    //flip_spec.partselectors.push_back( emodactive_logic.name_tag + "FLIPPED" );
   
    cmsswXmlInfo.specs.push_back(module_spec);
    //cmsswXmlInfo.specs.push_back(flip_spec);

    ecapRmatpath.push_back( ecapRmatpathCommon + "/" +
                            emodbox_shape.name_tag + "/" +
                            emod_shape.name_tag + "/" +
                            emodwafer_shape.name_tag + "/" +
                            emodactive_shape.name_tag );

    /*ecapRmatpath.push_back( ecapRmatpathCommon + "/" +
                            emodbox_shape.name_tag + "FLIPPED/" +
                            emod_shape.name_tag + "FLIPPED/" +
                            emodwafer_shape.name_tag + "FLIPPED/" +
                            emodactive_shape.name_tag + "FLIPPED" );
     */
        
    //analyseCompositeElements( ehybrid_logic.material_tag, compositeDensity(*emodule.getModuleCap(), true),
    //    *emodule.getModuleCap(), true);
    addMaterialToVolume(ehybrid_logic.material_tag, emodule, volumeType::hybrid,
                            (8.*ehybrid_shape.dx*ehybrid_shape.dy*ehybrid_shape.dz/1000.));
    addMaterialToVolume(chip_logic.material_tag, emodule, volumeType::chip,
                            (8.*chip_shape.dx*chip_shape.dy*chip_shape.dz/1000.));

  }
  void PixelExtractor::createPixelEndcapServices(InactiveSurfaces& inatcvser) {
    //Inactive elements part
    //InactiveSurfaces& inatvser = pix.getInactiveSurfaces();
    std::vector<InactiveElement>& bser = inatcvser.getEndcapServices();
    ShapeInfo bser_shape;
    bser_shape.type = tb;
    bser_shape.dx = 0.0;
    bser_shape.dy = 0.0;
    bser_shape.dyy = 0.0;


    LogicalInfo bser_logic;

    PosInfo bser_pos;
    bser_pos.copy = 1;
    bser_pos.trans.dx = 0.0;
    bser_pos.trans.dy = 0.0;

    Composite comp;

    comp.method = wt;

    for( auto& bs: bser ) {

      std::stringstream sname,matname;
      sname << xml_phaseII_pixecap  << xml_base_serf << "R" << (int)(bs.getInnerRadius())
        << "Z" << (int)( bs.getZOffset() + bs.getZLength() / 2.0 );
      bser_shape.name_tag = sname.str();
      bser_shape.rmin = bs.getInnerRadius();
      bser_shape.rmax = bs.getInnerRadius() + bs.getRWidth();
      bser_shape.dz = bs.getZLength() / 2.0;

      cmsswXmlInfo.shapes.push_back(bser_shape);

      matname << xml_base_serfcomp
        << "R" << (int)(bs.getInnerRadius())
        << "Z" << (int)(bs.getZLength());

      bser_logic.name_tag = sname.str();
      bser_logic.shape_tag =  xml_phaseII_pixecap + ":" + sname.str();
      bser_logic.material_tag = xml_phaseII_pixecap + ":" + matname.str();

      cmsswXmlInfo.logic.push_back(bser_logic);

      bser_pos.parent_tag = xml_phaseII_pixecap;
      bser_pos.child_tag = bser_logic.shape_tag;
      bser_pos.trans.dz = bs.getZOffset() + bser_shape.dz;

      cmsswXmlInfo.positions.push_back(bser_pos);

      std::stringstream bsermod;

      comp.name = bser_logic.material_tag;
      comp.density = compositeDensity(bs);;

      analyseCompositeElements( bser_logic.material_tag, comp.density,bs, false );
    }


  }

  int PixelExtractor::getAtomicNumber(double x0, double A) {
    // TODO: understand this: why do we need to
    // get an integer value as output?
    double d = 4 - 4 * (1.0 - 181.0 * A / x0);
    if (d > 0) return int(floor((sqrt(d) - 2.0) / 2.0 + 0.5));
    else return -1;
  }

  void PixelExtractor::analyseElements(MaterialTable& mattab) {
    for (unsigned int i = 0; i < mattab.rowCount(); i++) {
      Element e;
      MaterialRow& r = mattab.getMaterial(i);
      e.tag = r.tag;
      e.density = r.density;
      e.atomic_weight = pow((r.ilength / 35.), 3); // magic!
      e.atomic_number = getAtomicNumber(r.rlength, e.atomic_weight);
      //elems.push_back(e);
      cmsswXmlInfo.elements.push_back(e);
    }
  }

  void PixelExtractor::analyseCompositeElements( std::string name, double density, MaterialProperties& mp, bool nosensors ) {
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
    //////
    cmsswXmlInfo.composites.push_back(comp);
  }

//function to create composite material for different target volumes
  void PixelExtractor::addMaterialToVolume( std::string name, Module& module,int targetVoltag,double volume) {
     material::ElementsVector& modMaterialVec = module.getLocalElements();
    Composite comp;
    comp.name = name;
    comp.method = wt;
    double m = 0.0;
    for( auto& e : modMaterialVec) {
      if( e->targetVolume() == targetVoltag ) {
        std::pair<std::string, double> ele;
        ele.first = e->elementName();
        ele.second = e->quantityInGrams(module);
        m += e->quantityInGrams(module);
        bool flag = false;
        for( auto& ev : comp.elements ) {
            if(ev.first == ele.first) {
                ev.second += ele.second;
                flag = true;
                break;
            }
        }
          
        if(!flag)    comp.elements.push_back(ele);
      }
    }
    comp.density = m/volume;
    for( auto& ev : comp.elements )
        ev.second /= m;
    cmsswXmlInfo.composites.push_back(comp);

  }

  double PixelExtractor::compositeDensity(ModuleCap& mc, bool nosensors) {
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

  double PixelExtractor::compositeDensity(InactiveElement& ie) {
    double d = ie.getRWidth() + ie.getInnerRadius();
    d = d * d - ie.getInnerRadius() * ie.getInnerRadius();
    d = 1000 * ie.getTotalMass() / (M_PI * ie.getZLength() * d);
    return d;
  }

  void PixelExtractor::printXml(mainConfigHandler& mainConfiguration, std::string outsubdir) {
   
    std::string xmlpath = mainConfiguration.getXmlDirectory() + "/" + outsubdir + "/";
    std::cout<< "Xmls to be produced here=" << xmlpath<< std::endl;

    std::vector<ShapeInfo>& shapes = cmsswXmlInfo.shapes;
    std::vector<LogicalInfo>& logic = cmsswXmlInfo.logic;
    std::vector<PosInfo>& positions = cmsswXmlInfo.positions;
    std::vector<AlgoInfo>& algos = cmsswXmlInfo.algos;
    std::vector<Composite>& composites = cmsswXmlInfo.composites;
    std::vector<Element>& elements = cmsswXmlInfo.elements;

    ptree tree;
    tree.add("DDDefinition.<xmlattr>.xmlns", "http://www.cern.ch/cms/DDL");
    tree.add("DDDefinition.<xmlattr>.xmlns:xsi", "http://www.cern.ch/www.w3.org/2001/XMLSchema-instance");
    tree.add("DDDefinition.<xmlattr>.xsi:schemaLocation",
        "http://www.cern.ch/cms/DDL ../../../DetectorDescription/Schema/DDLSchema.xsd");

    ptree& matSec = tree.add("DDDefinition.MaterialSection", "");
    matSec.add("<xmlattr>.label", "pixel.xml");

    for( auto& e: elements ) {
      ptree& elem = matSec.add("ElementaryMaterial","");
      elem.add("<xmlattr>.name",e.tag);
      elem.add("<xmlattr>.symbol",e.tag);
      elem.add("<xmlattr>.atomicNumber",e.atomic_number);
      std::stringstream ss;
      ss << e.atomic_weight << "*g/mole";
      elem.add("<xmlattr>.atomicWeight",ss.str());
      ss.str("");
      ss << e.density << "*g/cm3";
      elem.add("<xmlattr>.density",ss.str());
    }

    for( auto& c: composites ) {
      ptree& comp = matSec.add("CompositeMaterial", "");
      comp.add("<xmlattr>.name",c.name);
      std::stringstream ss;
      ss << std::setprecision(3) << c.density << "*g/cm3";
      comp.add("<xmlattr>.density",ss.str());
      ss.str("");
      comp.add("<xmlattr>.method","mixture by weight");

          for( auto& e: c.elements ) {
          ptree& elem = comp.add("MaterialFraction", "");
          elem.add("<xmlattr>.fraction",e.second);
          ptree& m = elem.add("rMaterial", "");
          m.add("<xmlattr>.name", "pixel:" + e.first);
    
          }
    }

    ptree& solidSec = tree.add("DDDefinition.SolidSection", "");
    solidSec.add("<xmlattr>.label", "pixel.xml");
     
    for( auto& s: shapes) {
        if( s.type == ShapeType::bx ) {
        ptree& solid = solidSec.add("Box", "");
        solid.add("<xmlattr>.name",s.name_tag);
        std::stringstream ss;
        ss << std::setprecision(3) << s.dx << "*mm";
        solid.add("<xmlattr>.dx",ss.str());
        ss.str("");
        ss << std::setprecision(3) << s.dy << "*mm";
        solid.add("<xmlattr>.dy",ss.str());
        ss.str("");
        ss << std::setprecision(3) << s.dz << "*mm";
        solid.add("<xmlattr>.dz",ss.str());
      } else  if( s.type == ShapeType::tb ) {
        ptree& solid = solidSec.add("Tubs", "");
        solid.add("<xmlattr>.name",s.name_tag);
        std::stringstream ss;
        ss << std::setprecision(4) << s.rmin << "*mm";
        solid.add("<xmlattr>.rMin",ss.str());
        ss.str("");
        ss << std::setprecision(4) << s.rmax << "*mm";
        solid.add("<xmlattr>.rMax",ss.str());
        ss.str("");
        ss << std::setprecision(4) << s.dz << "*mm";
        solid.add("<xmlattr>.dz",ss.str());
        solid.add("<xmlattr>.startPhi","0*deg");
        solid.add("<xmlattr>.deltaPhi","360*deg");
      }
    }

    ptree& rotSec = tree.add("DDDefinition.RotationSection", "");
    rotSec.add("<xmlattr>.label", "pixel.xml");
    
    ptree& rotation1 = rotSec.add("Rotation","");
    rotation1.add("<xmlattr>.name","HCZ2YX");
    rotation1.add("<xmlattr>.thetaX","90*deg");
    rotation1.add("<xmlattr>.phiX","270*deg");
    rotation1.add("<xmlattr>.thetaY","180*deg");
    rotation1.add("<xmlattr>.phiY","0*deg");
    rotation1.add("<xmlattr>.thetaZ","90*deg");
    rotation1.add("<xmlattr>.phiZ","0*deg");
 
    ptree& rotation2 = rotSec.add("Rotation","");
    rotation2.add("<xmlattr>.name","FlippedHCZ2YX");
    rotation2.add("<xmlattr>.thetaX","90*deg");
    rotation2.add("<xmlattr>.phiX","270*deg");
    rotation2.add("<xmlattr>.thetaY","0*deg");
    rotation2.add("<xmlattr>.phiY","0*deg");
    rotation2.add("<xmlattr>.thetaZ","90*deg");
    rotation2.add("<xmlattr>.phiZ","180*deg");

    ptree& rotation3 = rotSec.add("Rotation","");
    rotation3.add("<xmlattr>.name","FLIP");
    rotation3.add("<xmlattr>.thetaX","90*deg");
    rotation3.add("<xmlattr>.phiX","180*deg");
    rotation3.add("<xmlattr>.thetaY","90*deg");
    rotation3.add("<xmlattr>.phiY","90*deg");
    rotation3.add("<xmlattr>.thetaZ","`180*deg");
    rotation3.add("<xmlattr>.phiZ","0*deg");
      
    ptree& logicSec = tree.add("DDDefinition.LogicalPartSection", "");
    logicSec.add("<xmlattr>.label", "pixel.xml");

    for( auto& l: logic) {
      ptree& logical = logicSec.add("LogicalPart","");
      logical.add("<xmlattr>.name",l.name_tag);
      logical.add("<xmlattr>.category", "unspecified");

      ptree& rsolid = logical.add("rSolid","");
      rsolid.add("<xmlattr>.name",l.shape_tag); 

      ptree& rmat = logical.add("rMaterial","");
      rmat.add("<xmlattr>.name",l.material_tag); 

    }

    ptree& posSec = tree.add("DDDefinition.PosPartSection", "");
    posSec.add("<xmlattr>.label", "pixel.xml");

    for( auto& p: positions) {

      ptree& position = posSec.add("PosPart","");
      position.add("<xmlattr>.copyNumber",p.copy);

      ptree& parent = position.add("rParent","");
      parent.add("<xmlattr>.name",p.parent_tag);

      ptree& child = position.add("rChild","");
      child.add("<xmlattr>.name",p.child_tag);

      if( p.rotref != "" ) {
        ptree& rot = position.add("rRotation","");
        rot.add("<xmlattr>.name",p.rotref);
      }

      if( p.trans.dx != 0 || p.trans.dy != 0 || p.trans.dz != 0) {
        ptree& translation = position.add("Translation","");
        std::stringstream ss;
        ss << std::setprecision(3) << p.trans.dx << "*mm";
        translation.add("<xmlattr>.dx",ss.str());
        ss.str("");
        ss << std::setprecision(3) << p.trans.dy << "*mm";
        translation.add("<xmlattr>.dy",ss.str());
        ss.str("");
        ss << std::setprecision(3) << p.trans.dz << "*mm";
        translation.add("<xmlattr>.dz",ss.str());
      }
    }

    for( auto& a: algos ) {
      ptree& algo = posSec.add("Algorithm","");
      algo.add("<xmlattr>.name",a.name);

      ptree& parent = algo.add("rParent","");
      parent.add("<xmlattr>.name",a.parent);
      if( !a.parameter_map.empty() ) {
        for( auto& p: a.parameter_map ) { 
          std::string ptype; 
          if( p.second.second == AlgoPartype::st ) 
            ptype = "String";
          else if( p.second.second == AlgoPartype::num )
            ptype = "Numeric";
          ptree& algoPar = algo.add(ptype,"");
          algoPar.add("<xmlattr>.name",p.first);
          algoPar.add("<xmlattr>.value",p.second.first);
        } 
      }
      if( a.vecpar.name != "" ) {
        ptree& algovPar = algo.add("Vector","");
        algovPar.add("<xmlattr>.name",a.vecpar.name);
        algovPar.add("<xmlattr>.type",a.vecpar.type);
        algovPar.add("<xmlattr>.nEntries",a.vecpar.nEntries);
        std::stringstream ss;
        for( unsigned int i = 0; i<a.vecpar.values.size() - 1; i++ )
          ss << a.vecpar.values[i] << ",";
        ss << a.vecpar.values[a.vecpar.values.size() - 1];
        algovPar.add("<xmltext>",ss.str());

      }
    }
    //xml_writer_settings<std::string> settings(' ', 1);//for new boost version
    xml_writer_settings<char> settings(' ', 1);
    write_xml(xmlpath+"pixel_test.xml", tree, std::locale(), settings);
 
    ///////////////writing pixel structure Topology////////////////////
    std::vector<SpecParInfo>& specs = cmsswXmlInfo.specs;
    ptree tree_topo;
    tree_topo.add("DDDefinition.<xmlattr>.xmlns", "http://www.cern.ch/cms/DDL");
    tree_topo.add("DDDefinition.<xmlattr>.xmlns:xsi", "http://www.cern.ch/www.w3.org/2001/XMLSchema-instance");
    tree_topo.add("DDDefinition.<xmlattr>.xsi:schemaLocation",
        "http://www.cern.ch/cms/DDL ../../../DetectorDescription/Schema/DDLSchema.xsd");
    ptree& specParSec = tree_topo.add("DDDefinition.SpecParSection", "");
    specParSec.add("<xmlattr>.label", xml_specpars_label);//xml_specpars_label=spec-p```ars2.xml

    ptree& spec = specParSec.add("SpecPar","");
    spec.add("<xmlattr>.name","FullTrackerPar");

    ptree& partSel = spec.add("PartSelector","");
    partSel.add("<xmlattr>.namepath","//Tracker");

    ptree& param = spec.add("Parameter","");
    param.add("<xmlattr>.name","TkDDDStructure");
    param.add("<xmlattr>.value","FullTracker");

    for( auto& s: specs ) {
      if( s.name.find("Module") != std::string::npos )      continue;
      // if( s.name.find("ROUHits") != std::string::npos)      continue;
      ptree& spec = specParSec.add("SpecPar","");
      spec.add("<xmlattr>.name",s.name);  
      for( auto& p: s.partselectors ) {
        ptree& partSel = spec.add("PartSelector","");
        partSel.add("<xmlattr>.path","//"+ p);
      }
      ptree& param = spec.add("Parameter","");
      param.add("<xmlattr>.name",s.parameter.first);
      param.add("<xmlattr>.value",s.parameter.second);
    }

    std::vector<std::string> barrel_partselectors;
    std::vector<std::string> endcap_partselectors;

    for( auto& s: specs ) {
      if( s.name.find("Module") == std::string::npos )                 continue;
      else if( s.name.find("BModule") != std::string::npos )      barrel_partselectors.push_back(s.name);
      else if( s.name.find("EModule") != std::string::npos )      endcap_partselectors.push_back(s.name);
      ptree& spec = specParSec.add("SpecPar","");
      spec.add("<xmlattr>.name",s.name);  
      for( auto& p: s.partselectors ) {
        ptree& partSel = spec.add("PartSelector","");
        partSel.add("<xmlattr>.path","//"+ p);
      }
      ptree& param = spec.add("Parameter","");
      param.add("<xmlattr>.name",s.parameter.first);
      param.add("<xmlattr>.value",s.parameter.second);
      for(auto m: s.moduletypes){
        ptree& param1 = spec.add("Parameter","");
        param1.add("<xmlattr>.PixelROCRows", m.rocrows);
        ptree& param2 = spec.add("Parameter","");
        param2.add("<xmlattr>.PixelROCCols", m.roccols);
        ptree& param3 = spec.add("Parameter","");
        param3.add("<xmlattr>.PixelROC_X", m.rocx);
        ptree& param4 = spec.add("Parameter","");
        param4.add("<xmlattr>.PixelROC_Y", m.rocy);
      } 
    }
    write_xml(xmlpath+"pixelStructureTopology_test.xml", tree_topo, std::locale(), settings);

    //sensor portion
    ptree tree_sense;
    tree_sense.add("DDDefinition.<xmlattr>.xmlns", "http://www.cern.ch/cms/DDL");
    tree_sense.add("DDDefinition.<xmlattr>.xmlns:xsi", "http://www.cern.ch/www.w3.org/2001/XMLSchema-instance");
    tree_sense.add("DDDefinition.<xmlattr>.xsi:schemaLocation",
        "http://www.cern.ch/cms/DDL ../../../DetectorDescription/Schema/DDLSchema.xsd");
    ptree& specParSensorSec = tree_sense.add("DDDefinition.SpecParSection", "");
    specParSensorSec.add("<xmlattr>.label", xml_specpars_label);//xml_specpars_label=spec-pars2.xml
    //Barrel Part 
    ptree& barrel_specSensor = specParSensorSec.add("SpecPar","");
    barrel_specSensor.add("<xmlattr>.name","ROUHitsTracker" + xml_phaseII_pixbar);  
    for( auto& p: barrel_partselectors ) {
      ptree& partSel = barrel_specSensor.add("PartSelector","");
      p.erase(p.end()-3,p.end());
      partSel.add("<xmlattr>.path","//"+ p);
    }
    ptree& barrel_param1 = barrel_specSensor.add("Parameter","");
    barrel_param1.add("<xmlattr>.name","SensitiveDetector");
    barrel_param1.add("<xmlattr>.value","TkAccumulatingSensitiveDetector");

    ptree& barrel_param2 = barrel_specSensor.add("Parameter","");
    barrel_param2.add("<xmlattr>.name","ReadOutName");
    barrel_param2.add("<xmlattr>.value","TrackerHits" + xml_phaseII_pixbar);
   
    //Endcap Part 
    ptree& endcap_specSensor = specParSensorSec.add("SpecPar","");
    endcap_specSensor.add("<xmlattr>.name","ROUHitsTracker" + xml_phaseII_pixecap);  
    for( auto& p: endcap_partselectors ) {
      ptree& partSel = endcap_specSensor.add("PartSelector","");
      p.erase(p.end()-3,p.end());
      partSel.add("<xmlattr>.path","//"+ p);
    }
    ptree& endcap_param1 = endcap_specSensor.add("Parameter","");
    endcap_param1.add("<xmlattr>.name","SensitiveDetector");
    endcap_param1.add("<xmlattr>.value","TkAccumulatingSensitiveDetector");

    ptree& endcap_param2 = endcap_specSensor.add("Parameter","");
    endcap_param2.add("<xmlattr>.name","ReadOutName");
    endcap_param2.add("<xmlattr>.value","TrackerHits" +xml_phaseII_pixecap);

    write_xml(xmlpath+"pixelsens_test.xml", tree_sense, std::locale(), settings);
    
    //Prodcut portion
    ptree tree_prodCut;
    tree_prodCut.add("DDDefinition.<xmlattr>.xmlns", "http://www.cern.ch/cms/DDL");
    tree_prodCut.add("DDDefinition.<xmlattr>.xmlns:xsi", "http://www.cern.ch/www.w3.org/2001/XMLSchema-instance");
    tree_prodCut.add("DDDefinition.<xmlattr>.xsi:schemaLocation",
        "http://www.cern.ch/cms/DDL ../../../DetectorDescription/Schema/DDLSchema.xsd");
    ptree& specParProdSec = tree_prodCut.add("DDDefinition.SpecParSection", "");
    specParProdSec.add("<xmlattr>.label", "trackerProdCuts.xml");//xml_specpars_label=spec-pars2.xml
    specParProdSec.add("<xmlattr>.eval","true");

    ptree& dead_spec = specParProdSec.add("SpecPar","");
    dead_spec.add("<xmlattr>.name","tracker-dead-pixel");

    ptree& barrel_partSel = dead_spec.add("PartSelctor","");
    barrel_partSel.add("<xmlattr>.path","//" + xml_phaseII_pixbar );

    ptree& endcap_partSel = dead_spec.add("PartSelector","");
    endcap_partSel.add("<xmlattr>.path","//" + xml_phaseII_pixecap);

    ptree& dead_param1 = dead_spec.add("Parameter","");
    dead_param1.add("<xmlattr>.name","CMSCutsRegion");
    dead_param1.add("<xmlattr>.value","TrackerPixelDeadRegion");
    dead_param1.add("<xmlattr>.eval","false");

    ptree& dead_param2 = dead_spec.add("Parameter","");
    dead_param2.add("<xmlattr>.name","ProdCutsForElectrons");
    dead_param2.add("<xmlattr>.value","1*mm");

    ptree& dead_param3 = dead_spec.add("Parameter","");
    dead_param3.add("<xmlattr>.name","ProdCutsForPositrons");
    dead_param3.add("<xmlattr>.value","1*mm");

    ptree& dead_param4 = dead_spec.add("Parameter","");
    dead_param4.add("<xmlattr>.name","ProdCutsForGamma");
    dead_param4.add("<xmlattr>.value","1*mm");


    ptree& sens_spec = specParProdSec.add("SpecPar","");
    sens_spec.add("<xmlattr>.name","tracker-sens-pixel");
    for( auto& p: barrel_partselectors ) {
      ptree& partSel = sens_spec.add("PartSelector","");
      partSel.add("<xmlattr>.path","//"+ p);
    }
    for( auto& p: endcap_partselectors ) {
      ptree& partSel = sens_spec.add("PartSelector","");
      partSel.add("<xmlattr>.path","//"+ p);
    }
    ptree& sens_param1 = sens_spec.add("Parameter","");
    sens_param1.add("<xmlattr>.name","CMSCutsRegion");
    sens_param1.add("<xmlattr>.value","TrackerPixelSensRegion");
    sens_param1.add("<xmlattr>.eval","false");

    ptree& sens_param2 = sens_spec.add("Parameter","");
    sens_param2.add("<xmlattr>.name","ProdCutsForElectrons");
    sens_param2.add("<xmlattr>.value","0.01*mm");

    ptree& sens_param3 = sens_spec.add("Parameter","");
    sens_param3.add("<xmlattr>.name","ProdCutsForPositrons");
    sens_param3.add("<xmlattr>.value","0.01*mm");

    ptree& sens_param4 = sens_spec.add("Parameter","");
    sens_param4.add("<xmlattr>.name","ProdCutsForGamma");
    sens_param4.add("<xmlattr>.value","0.01*mm");

    write_xml(xmlpath+"pixelProdCuts_test.xml", tree_prodCut, std::locale(), settings);
 
   //Reco Material
    ptree tree_recoMat;
    tree_recoMat.add("DDDefinition.<xmlattr>.xmlns", "http://www.cern.ch/cms/DDL");
    tree_recoMat.add("DDDefinition.<xmlattr>.xmlns:xsi", "http://www.cern.ch/www.w3.org/2001/XMLSchema-instance");
    tree_recoMat.add("DDDefinition.<xmlattr>.xsi:schemaLocation",
        "http://www.cern.ch/cms/DDL ../../../DetectorDescription/Schema/DDLSchema.xsd");
    ptree& specParRecoSec = tree_recoMat.add("DDDefinition.SpecParSection", "");
    specParRecoSec.add("<xmlattr>.label", "spec-pars2.xml");//xml_specpars_label=spec-pars2.xml
   
      std::string pixbarRecoSpeccommon = "TrackerRecMaterial" + xml_phaseII_pixbar + xml_layer;
      for( unsigned int l = 1; l<=numBarrelLayers; l++){
          std::stringstream stemp;
          stemp << pixbarRecoSpeccommon << l;
          ptree& specReco = specParRecoSec.add("SpecPar","");
          specReco.add("<xmlattr>.name",stemp.str());
          specReco.add("<xmlattr>.eval","true");
          for(auto& e:barrelRmatpath) {
            std::stringstream sl;
              sl << "Layer" << l;
              if( e.find(sl.str()) != std::string::npos ){
                 ptree& part = specReco.add("PartSelector","");
                 part.add("<xmlattr>.path",e);
              }
          }
      }
      
      std::string pixfwdRecoSpeccommon = "TrackerRecMaterial" + xml_phaseII_pixecap + "Disk";
      for( unsigned int d = 0; d<discRingpair.size(); d++){
          std::stringstream stemp;
          stemp << pixfwdRecoSpeccommon << discRingpair[d].first << "Fw";
          ptree& specReco = specParRecoSec.add("SpecPar","");
          specReco.add("<xmlattr>.name",stemp.str());
          specReco.add("<xmlattr>.eval","true");
          for( unsigned int r = 0; r<discRingpair[d].second; r++){
            for(auto& e:ecapRmatpath) {
                std::stringstream sl;
                sl << "Ring" << r+1 << "Disc" << discRingpair[d].first;
                if( e.find(sl.str()) != std::string::npos ){
                    ptree& part = specReco.add("PartSelector","");
                    part.add("<xmlattr>.path",e);
                }
            }
         }
      }
      write_xml(xmlpath+"pixelRecoMaterial_test.xml", tree_recoMat, std::locale(), settings);
  }
}//for namespace

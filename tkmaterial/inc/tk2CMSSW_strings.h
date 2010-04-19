#ifndef _TK2CMSSW_STRINGS_H
#define _TK2CMSSW_STRINGS_H

/**
 * @file tk2CMSSW_strings.h
 * @brief This header file lists XML text constants that are combined and concatenated with various parameter values during translation of a geometry to CMSSW XML.
 */

#include <string>

namespace insur {
    /**
     * Numeric constants
     */
    static const int xml_prec = 3;
    static const int xml_roc_rows = 128;
    static const int xml_roc_cols = 1;
    /**
     * XML tags and attributes
     */
    static const std::string xml_preamble = "<?xml version=\"1.0\"?>\n<DDDefinition xmlns=\"http://www.cern.ch/cms/DDL\" xmlns:xsi=\"http://www.cern.ch/www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.cern.ch/cms/DDL ../../../DetectorDescription/Schema/DDLSchema.xsd\">\n";
    static const std::string xml_defclose = "</DDDefinition>\n";
    static const std::string xml_general_inter = "\">\n";
    static const std::string xml_general_endline = "\"/>\n";
    static const std::string xml_const_section = "<ConstantsSection label=\"tobrodpar.xml\" eval=\"true\">\n<Constant name=\"BackPlaneDz\" value=\"0.015*mm\"/>\n</ConstantsSection>\n";
    static const std::string xml_recomat_parameters = "<Parameter name=\"TrackerRadLength\" value=\"0.01\"/>\n<Parameter name=\"TrackerXi\" value=\"0.0001";
    static const std::string xml_recomat_radlength = "TrackerRadLength";
    static const std::string xml_recomat_xi = "TrackerXi";
    static const std::string xml_eval_true = "\" eval=\"true\">\n";
    static const std::string xml_material_section_open = "<MaterialSection label=\"";
    static const std::string xml_material_section_close = "</MaterialSection>\n";
    static const std::string xml_rotation_section_open = "<RotationSection label=\"";
    static const std::string xml_rotation_section_close = "</RotationSection>\n";
    static const std::string xml_logical_part_section_open = "<LogicalPartSection label=\"";
    static const std::string xml_logical_part_section_close = "</LogicalPartSection>\n";
    static const std::string xml_solid_section_open = "<SolidSection label=\"";
    static const std::string xml_solid_section_close = "</SolidSection>\n";
    static const std::string xml_pos_part_section_open = "<PosPartSection label=\"";
    static const std::string xml_pos_part_section_close = "</PosPartSection>\n";
    static const std::string xml_spec_par_section_open = "<SpecParSection label=\"";
    static const std::string xml_spec_par_section_close = "</SpecParSection>\n";
    static const std::string xml_algorithm_open = "<Algorithm name=\"";
    static const std::string xml_algorithm_parent = "\">\n<rParent name=\"";
    static const std::string xml_algorithm_string = "<String name=\"";
    static const std::string xml_algorithm_numeric = "<Numeric name=\"";
    static const std::string xml_algorithm_vector_open = "<Vector name=\"Center\" type=\"numeric\" nEntries=\"3\">";
    static const std::string xml_algorithm_vector_close = "</Vector>\n";
    static const std::string xml_algorithm_value = "\" value=\"";
    static const std::string xml_algorithm_close = "</Algorithm>\n";
    static const std::string xml_elementary_material_open = "<ElementaryMaterial name=\"";
    static const std::string xml_elementary_material_first_inter = "\" symbol=\"";
    static const std::string xml_elementary_material_second_inter = "\" atomicNumber=\"";
    static const std::string xml_elementary_material_third_inter = "\" atomicWeight=\"";
    static const std::string xml_elementary_material_fourth_inter = "*g/mole\" density=\"";
    static const std::string xml_elementary_material_close = "*g/cm3\"/>\n";
    static const std::string xml_composite_material_open = "<CompositeMaterial name=\"";
    static const std::string xml_composite_material_first_inter = "\" density=\"";
    static const std::string xml_composite_material_second_inter = "*g/cm3\" method=\"";
    static const std::string xml_composite_material_close = "</CompositeMaterial>\n";
    static const std::string xml_material_fraction_open = "<MaterialFraction fraction=\"";
    static const std::string xml_material_fraction_inter = "\">\n<rMaterial name=\"";
    static const std::string xml_material_fraction_close = "\"/>\n</MaterialFraction>\n";
    static const std::string xml_logical_part_open = "<LogicalPart name=\"";
    static const std::string xml_logical_part_first_inter = "\" category=\"unspecified\">\n<rSolid name=\"";
    static const std::string xml_logical_part_second_inter = "\"/>\n<rMaterial name=\"";
    static const std::string xml_logical_part_close = "\"/>\n</LogicalPart>\n";
    static const std::string xml_box_open = "<Box name=\"";
    static const std::string xml_box_first_inter = "\" dx=\"";
    static const std::string xml_box_second_inter = "*mm\" dy=\"";
    static const std::string xml_box_third_inter = "*mm\" dz=\"";
    static const std::string xml_box_close = "*mm\"/>\n";
    static const std::string xml_trapezoid_open = "<Trd1 name=\"";
    static const std::string xml_trapezoid_first_inter = "\" dx1=\"";
    static const std::string xml_trapezoid_second_inter = "*mm\" dx2=\"";
    static const std::string xml_trapezoid_third_inter = "*mm\" dy1=\"";
    static const std::string xml_trapezoid_fourth_inter = "*mm\" dy2=\"";
    static const std::string xml_trapezoid_fifth_inter = "*mm\" dz=\"";
    static const std::string xml_trapezoid_close = "*mm\"/>\n";
    static const std::string xml_tubs_open = "<Tubs name=\"";
    static const std::string xml_tubs_first_inter = "\" rMin=\"";
    static const std::string xml_tubs_second_inter = "*mm\" rMax=\"";
    static const std::string xml_tubs_third_inter = "*mm\" dz=\"";
    static const std::string xml_tubs_close = "*mm\" startPhi=\"0*deg\" deltaPhi=\"360*deg\"/>\n";
    static const std::string xml_polycone_open = "<Polycone name=\"";
    static const std::string xml_polycone_inter = "\" startPhi=\"0*deg\" deltaPhi=\"360*deg\">\n";
    static const std::string xml_polycone_close = "</Polycone>\n";
    static const std::string xml_rzpoint_open = "<RZPoint r=\"";
    static const std::string xml_rzpoint_inter = "*mm\" z=\"";
    static const std::string xml_rzpoint_close = "*mm\"/>\n";
    static const std::string xml_pos_part_open = "<PosPart copyNumber=\"";
    static const std::string xml_pos_part_first_inter = "\">\n<rParent name=\"";
    static const std::string xml_pos_part_second_inter = "\"/>\n<rChild name=\"";
    static const std::string xml_pos_part_third_inter = "<rRotation name=\"";
    static const std::string xml_pos_part_close = "</PosPart>\n";
    static const std::string xml_rotation_open = "<Rotation name=\"";
    static const std::string xml_rotation_first_inter = "\" thetaX=\"";
    static const std::string xml_rotation_second_inter = "*deg\" phiX=\"";
    static const std::string xml_rotation_third_inter = "*deg\" thetaY=\"";
    static const std::string xml_rotation_fourth_inter = "*deg\" phiY=\"";
    static const std::string xml_rotation_fifth_inter = "*deg\" thetaZ=\"";
    static const std::string xml_rotation_sixth_inter = "*deg\" phiZ=\"";
    static const std::string xml_rotation_close = "*deg\"/>\n";
    static const std::string xml_translation_open = "<Translation x=\"";
    static const std::string xml_translation_first_inter = "*mm\" y=\"";
    static const std::string xml_translation_second_inter = "*mm\" z=\"";
    static const std::string xml_translation_close = "*mm\"/>\n";
    static const std::string xml_spec_par_open = "<SpecPar name=\"";
    static const std::string xml_spec_par_selector = "<PartSelector path=\"//";
    static const std::string xml_spec_par_parameter_first = "<Parameter name=\"";
    static const std::string xml_spec_par_parameter_second = "\" value=\"";
    static const std::string xml_spec_par_close = "\"/>\n</SpecPar>\n";
    /**
     * Output filenames
     */
    static const std::string xml_trackerfile = "tracker.xml";
    static const std::string xml_pixbarfile = "pixbar.xml";
    static const std::string xml_pixfwdfile = "pixfwd.xml";
    static const std::string xml_topologyfile = "trackerStructureTopology.xml";
    static const std::string xml_prodcutsfile = "trackerProdCuts.xml";
    static const std::string xml_trackersensfile = "trackersens.xml";
    static const std::string xml_recomatfile = "trackerRecoMaterial.xml";
    static const std::string xml_tmppath = "tmp";
    /**
     * Naming conventions and variable names
     */
    static const std::string xml_insert_marker = "<!--mid point marker-->";
    static const std::string xml_specpars_label = "spec-pars2.xml";
    static const std::string xml_base_act = "active";
    static const std::string xml_base_waf = "wafer";
    static const std::string xml_base_serf = "service";
    static const std::string xml_base_lazy = "support";
    static const std::string xml_layer = "Layer";
    static const std::string xml_disc = "Disc";
    static const std::string xml_rod = "Rod";
    static const std::string xml_ring = "Ring";
    static const std::string xml_plus = "Plus";
    static const std::string xml_minus = "Minus";
    static const std::string xml_barrel_module = "BModule";
    static const std::string xml_endcap_module = "EModule";
    static const std::string xml_base_actcomp = "modulecomposite";
    static const std::string xml_base_serfcomp = "servicecomposite";
    static const std::string xml_base_lazycomp = "supportcomposite";
    static const std::string xml_material_air = "materials:Air";
    static const std::string xml_sensor_silicon = "SenSi";
    static const std::string xml_pixbarident = "pixbar";
    static const std::string xml_pixfwdident = "pixfwd";
    static const std::string xml_fileident = "tracker";
    static const std::string xml_pixbar = "PixelBarrel";
    //TODO: define pixel forward root name
    static const std::string xml_tracker = "Tracker";
    static const std::string xml_tob = "TOB";
    static const std::string xml_tid = "TID";
    static const std::string xml_tobalgo = "track:DDTrackerPhiAltAlgo";
    static const std::string xml_ecalgo = "track:DDTrackerAngular";
    static const std::string xml_param_string = "String";
    static const std::string xml_param_numeric = "Numeric";
    static const std::string xml_childparam = "ChildName";
    static const std::string xml_tilt = "Tilt";
    static const std::string xml_startangle = "StartAngle";
    static const std::string xml_rangeangle = "RangeAngle";
    static const std::string xml_radiusin = "RadiusIn";
    static const std::string xml_radiusout = "RadiusOut";
    static const std::string xml_zposition = "ZPosition";
    static const std::string xml_number = "Number";
    static const std::string xml_startcopyno = "StartCopyNo";
    static const std::string xml_incrcopyno = "IncrCopyNo";
    static const std::string xml_radius = "Radius";
    static const std::string xml_nmods = "N";
    static const std::string xml_tkddd_structure = "TkDDDStructure";
    static const std::string xml_full_tracker = "FullTracker";//maz be obsolete
    static const std::string xml_det_layer = "TOBLayer";
    static const std::string xml_det_rod = "TOBRod";
    static const std::string xml_det_tobdet = "TOBDet";
    static const std::string xml_tob_subdet = "TOBSubDet";
    static const std::string xml_subdet_layer = "PixelBarrelLayer";
    static const std::string xml_subdet_rod = "PixelBarrelLadder";
    static const std::string xml_subdet_tobdet = "PixelBarrelModule";
    static const std::string xml_det_wheel = "TIDWheel";
    static const std::string xml_det_ring = "TIDRing";
    static const std::string xml_det_tiddet = "TIDDet";
    static const std::string xml_tid_subdet = "TIDSubDet";
    static const std::string xml_subdet_wheel = "PixelEndcapDisk";
    static const std::string xml_subdet_ring = "PixelEndcapPanel";
    static const std::string xml_subdet_tiddet = "PixelEndcapDet";
    static const std::string xml_apv_head = "TrackerAPVNumber";
    static const std::string xml_roc_x = "PixelROC_X";
    static const std::string xml_roc_y = "PixelROC_Y";
    static const std::string xml_par_tail = "Par";
    static const std::string xml_tob_prefix = "TrackerRecMaterialTOB";
    static const std::string xml_tid_prefix = "TrackerRecMaterialTIDDisk";
    static const std::string xml_forward = "Fw";
    static const std::string xml_backward = "Bw";
    static const std::string xml_barrel_tilt = "HCZ2YX";
    /**
     * CMSSW constants
     */
    static const std::string xml_zv3 = "[Zv3]";
}
#endif /* _TK2CMSSW_STRINGS_H */

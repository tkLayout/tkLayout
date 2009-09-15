#ifndef _TK2CMSSW_STRINGS_H
#define _TK2CMSSW_STRINGS_H

/**
 * @file tk2CMSSW_strings.h
 * @brief This header file lists XML text constants that are combined and concatenated with various parameter values during translation of a geometry to CMSSW XML.
 */

#include <string>

namespace insur {
    static const std::string xml_preamble = "<?xml version=1.0?><DDDefinition xmlns=\"http://www.cern.ch/cms/DDL\" xmlns:xsi=\"http://www.cern.ch/www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.cern.ch/cms/DDL ../../../DetectorDescription/Schema/DDLSchema.xsd\">";
    static const std::string xml_defclose = "</DDDefinition>";
    static const std::string xml_material_section_open = "<MaterialSection label=\"";
    static const std::string xml_material_section_inter = "\">\n";
    static const std::string xml_material_section_close = "</MaterialSection>\n";
    static const std::string xml_logical_part_section_open = "<LogicalPartSection label=\"";
    static const std::string xml_logical_part_section_inter = "\">\n";
    static const std::string xml_logical_part_section_close = "</LogicalPartSection>\n";
    static const std::string xml_solid_section_open = "<SolidSection label=\"";
    static const std::string xml_solid_section_inter = "\">\n";
    static const std::string xml_solid_section_close = "</SolidSection>\n";
    static const std::string xml_pos_part_section_open = "<PosPartSection label=\"";
    static const std::string xml_pos_part_section_inter = "\">\n";
    static const std::string xml_pos_part_section_close = "</PosPartSection>\n";
    static const std::string xml_elementary_material_open = "<ElementaryMaterial name=\"";
    static const std::string xml_elementary_material_first_inter = "\" symbol=\"";
    static const std::string xml_elementary_material_second_inter = "\" atomicNumber=\"";
    static const std::string xml_elementary_material_third_inter = "\" atomicWeight=\"";
    static const std::string xml_elementary_material_fourth_inter = "*g/mole\" density=\"";
    static const std::string xml_elementary_material_close = "*mg/cm3\"/>\n";
    static const std::string xml_composite_material_open = "<CompositeMaterial name=\"";
    static const std::string xml_composite_material_first_inter = "\" density=\"";
    static const std::string xml_composite_material_second_inter = "*g/cm3\" method=\"";
    static const std::string xml_composite_material_third_inter = "\">\n";
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
    static const std::string xml_pos_part_open = "<PosPart copyNumber=\"";
    static const std::string xml_pos_part_first_inter = "\">\n<rParent name=\"";
    static const std::string xml_pos_part_second_inter = "\"/>\n<rChild name=\"";
    static const std::string xml_pos_part_third_inter = "\"/>\n";
    static const std::string xml_pos_part_close = "</PosPart>\n";
    static const std::string xml_rotation_open = "<Rotation name=\"";
    static const std::string xml_rotation_first_inter = "\" phiX=\"";
    static const std::string xml_rotation_second_inter = "*rad\" phiY=\"";
    static const std::string xml_rotation_third_inter = "*rad\" phiZ=\"";
    static const std::string xml_rotation_fourth_inter = "*rad\" thetaX=\"";
    static const std::string xml_rotation_fifth_inter = "*rad\" thetaY=\"";
    static const std::string xml_rotation_sixth_inter = "*rad\" thetaZ=\"";
    static const std::string xml_rotation_close = "*rad\"/>\n";
    static const std::string xml_translation_open = "<Translation x=\"";
    static const std::string xml_translation_first_inter = "*mm\" y=\"";
    static const std::string xml_translation_second_inter = "*mm\" z=\"";
    static const std::string xml_translation_close = "*mm\"/>\n";
}
#endif /* _TK2CMSSW_STRINGS_H */

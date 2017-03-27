#ifndef _TK2CMSSW_DATATYPES_H
#define _TK2CMSSW_DATATYPES_H

/**
 * @file tk2CMSSW_datatypes.h
 * @brief This header file lists the internal structs and enumerations needed to translate a full tracker and material budget to CMSSW XML
 */

#include <string>
#include <vector>
#include <map>
#include <math.h> 
#include <tk2CMSSW_strings.hh>

namespace insur {
    /**
     * @enum CompType A list of the options for composite materials in CMSSW: mixture by weight, by volume or by atomic proportion
     */
    enum CompType { wt, vl, ap };
    /**
     * @enum ShapeType A list of possible shape types: box, tubs, cones, tdr1 (trapezoid) and polycone
     */
    enum ShapeType { bx, tb, co, tp, pc };

    enum ShapeOperationType { uni, intersec, substract };

    enum AlgoPartype { st,num,vec};

    /**
     * @struct Rotation
     * @brief This struct collects the parameters that describe a rotation in 3D.
     * @param name The name of the rotation, if it has one
     * @param phix (desc here) TODO
     * @param phiy (desc here) TODO
     * @param phiz (desc here) TODO
     * @param thetax (desc here) TODO
     * @param thetay (desc here) TODO
     * @param thetaz (desc here) TODO
     */
    struct Rotation {
        std::string name;
        double phix;
        double phiy;
        double phiz;
        double thetax;
        double thetay;
        double thetaz;
    };
    /**
     * @struct Translation
     * @brief This struct collects the parameters (in mm) that describe a translation in 3D.
     * @param dx The shift along the x-axis
     * @param dy The shift along the y-axis
     * @param dz The shift along the z-axis
     */
    struct Translation {
        double dx;
        double dy;
        double dz;
    };
    /**
     * @struct Element
     * @brief This struct collects some of the parameters of a chemical element.
     * @param tag The name of the element
     * @param density The average density of the element (in g/cm3)
     * @param atomic_number The atomic number of the element
     * @param atomic_weight The average atomic weight of the element (in g/mole)
     */
    struct Element {
        std::string tag;
        double density;
        int atomic_number;
        double atomic_weight;
    };
    /**
     * @struct Composite
     * @brief This struct collects some global properties and the chemical elements that make up a composite material.
     * @param name The name of the composite material
     * @param density The overall density of the composite material (in g/cm3)
     * @param method The type of mixture: by weight, by volume or by atomic proportion
     * @param elements A list of elements in the compound, key = name of the element, value = massic fraction of the whole
     */
    struct Composite {
        std::string name;
        double density;
        CompType method;
        std::map<std::string, double> elements;

        // This is to avoid the duplicated descriptions of composite materials in the XMLs (very significant effect on total XML size)
        // 2 components are said equal if they have same total density, same mixture method, and exactly same composing elements.
        // No matter the name of the component !
        bool operator==(const Composite& otherComp) const {
	  if ((fabs(density - otherComp.density) > xml_composite_density_tolerance)  // not same density ?
	      || (method != otherComp.method)  // not same mixture method ?
	      || (elements.size() != otherComp.elements.size())) {  // not same number of elements ?
	    return false;
  	  }
	  for (const auto& elem : otherComp.elements) {  // for a given element :
	    if ((elements.find(elem.first) == elements.end())  // not found in the composite ?
	        || (fabs(elements.at(elem.first) - elem.second) > xml_composite_ratio_tolerance)) { // not same massic ratio ?
	      return false;
	  }
	}
	return true;
      }
    };
    /**
     * @struct LogicalInfo
     * @brief This struct corresponds to a <i>LogicalPart</i> block in CMSSW XML.
     * @param name_tag The name of the logical part
     * @param shape_tag A name referencing the associated shape
     *@param material_tag A name referencing the associated material
     */
    struct LogicalInfo {
        std::string name_tag;
        std::string shape_tag;
        std::string material_tag;
        std::string extra; // CUIDADO this was introduced to add the Plus / Minus attribute to discs, not really useful for anything else at least for now
    };
    /**
     * @struct ShapeInfo
     * @brief This struct corresponds to a geometrical shape in a <i>SolidSection</i> block in CMSSW XML (all measurements in mm).
     * @param type The shape type: one of box, tube section, cone section, trapezoid or polycone
     * @param name_tag The name of the shape volume
     * @param dx Half the width along the x-axis of the bottom side (boxes and trapezoids)
     * @param dxx Half the width along the x-axis at the top side (trapezoids only)
     * @param dy Half the height along the y-axis of the left side (boxes and trapezoids)
     * @param dyy Half the height along the y-axis of the right side (trapezoids only)
     * @param dz Half the depth along the z-axis (boxes, trapezoids, tub and cone sections)
     * @param rmin The minimal radius as measured from the z-axis (tube sections only) 
     * @param rmax The maximum radius as measured from the z-axis (tube sections only)
     * @param rmin1 The minimal radius as measured from the z-axis, for the smallest-z section (cone sections only) 
     * @param rmax1 The maximum radius as measured from the z-axis, for the smallest-z section (cone sections only)
     * @param rmin2 The minimal radius as measured from the z-axis, for the biggest-z section (cone sections only) 
     * @param rmax2 The maximum radius as measured from the z-axis, for the biggest-z section (cone sections only)
     * @param rzup A vector of pairs collecting the first half of an ordered sequence of points in <i>(r, z)</i> (polycone only)
     * @param rzdown A vector of pairs collecting the second half of an ordered sequence of points in <i>(r, z)</i> (polycone only)
     */
    struct ShapeInfo {
        ShapeType type;
        std::string name_tag;
        double dx;
        double dxx;
        double dy;
        double dyy;
        double dz;
        double rmin;
        double rmax;
        double rmin1;
        double rmax1;
        double rmin2;
        double rmax2;
      std::vector<std::pair<double, double> > rzup;
      std::vector<std::pair<double, double> > rzdown;
    };
    /**
     * @struct ShapeOperationInfo
     * @brief This struct corresponds to operations on geometrical shapes in a <i>SolidSection</i> block in CMSSW XML.
     * @param type The type of operation on volumes
     * @param name_tag The name of the result volume of the operation
     * @param rSolid1 The name of one of the volume the operation is made on
     * @param rSolid2 The name of a second volume the operation is made on
     * @param trans A translation applied to the child volume within the parent volume's coordinate system
     */
    struct ShapeOperationInfo {
        ShapeOperationType type;
        std::string name_tag;
        std::string rSolid1;
        std::string rSolid2;
        Translation trans;
    };
    /**
     * @struct PosInfo
     * @brief This struct corresponds to a <i>PosPart</i> block in CMSSW XML, positioning a volume in 3D space and within the hierarchy of <i>LogicalPart</i>s.
     * @param parent_tag A name referencing the parent logical part
     * @param child_tag A name referencing the child logical part
     * @param copy The number of the duplicated original shape
     * @param rot A rotation applied to the child volume within the parent volume's coordinate system
     * @param trans A translation applied to the child volume within the parent volume's coordinate system
     */
    struct PosInfo {
        std::string parent_tag;
        std::string child_tag;
        int copy;
        std::string rotref;
        Translation trans;
    };

    struct VecInfo {
      std::string name;
      std::string type;
      std::string nEntries;
      std::vector<double> values;
    };

    /**
     * @struct AlgoInfo
     * @brief This struct collects the parameters that go into an <i>Algorithm</i> block in CMSSW XML.
     * @param name A name referencing the volume placement algorithm to be used
     * @param parent A logical part name referencing the parent (container) volume and coordinate system
     * @param parameters A series of pre-formatted parameter strings that the algorithm receives as arguments
     */
    struct AlgoInfo {
        std::string name;
        std::string parent;
        std::vector<std::string> parameters;
        std::map<std::string,std::pair<std::string,AlgoPartype> > parameter_map;
        VecInfo vecpar;
    };
    /**
     * @struct ERingInfo
     * @brief This is a struct to collect temporary information about an endcap ring and the modules within it.
     * @param name The logical part name that identifies the ring
     * @param childname The logical part name that identifies the modules contained in the ring
     * @param fw Is it the forward ring of the disk ?
     * @param isZPlus Is the ring (and disk) in the positive-z side ?
     * @param fw_flipped Are modules in the forward part (big |z|) of the ring flipped ?
     * @param phi The angle <i>phi</i> in the x/y-plane of the first module on the ring
     * @param modules The number of modules within the ring
     * @param mthk Thickness of one of the ring's modules (hybrids included)
     * @param rmin The minimum radius of the ring, as measured from the z-axis 
     * @param rmid The radius of the module mean point, as measured from the z-axis
     * @param rmax The maximum radius of the ring, as measured from the z-axis
     */
    struct ERingInfo {
        std::string name;
        std::string childname;
        bool fw;
        bool isZPlus;
        bool fw_flipped;
        int modules;
        double mthk;
        double rmin;
        double rmid;
        double rmax;
        double zmin;
        double zmax;
        double zfw;
        double startPhiAnglefw;  // in RAD
        double zbw;
        double startPhiAnglebw;  // in RAD
    };
    /**
     * @struct BTiltedRingInfo
     * @brief This is a struct to collect temporary information about a barrel tilted ring and the modules within it.
     * @param name The logical part name that identifies the ring
     * @param childname The logical part name that identifies the modules contained in the ring
     * @param isZPlus Is the ring in the positive-z side ?
     * @param tiltAngle Module's tilt angle
     * @param bw_flipped Are modules in the backward part (small |z|) of the ring flipped ?
     * @param fw_flipped Are modules in the forward part (big |z|) of the ring flipped ?
     * @param phi The angle <i>phi</i> in the x/y-plane of the first module on the ring
     * @param modules The number of modules within the ring
     * @param r1 The radius of the module mean point, as measured from the z-axis, for the inner part (small |z|) of the ring
     * @param r2 The radius of the module mean point, as measured from the z-axis, for the outer part (big |z|) of the ring
     * @param z1 z of the module's mean point, for the inner part (small |z|) of the ring
     * @param z2 z of the module's mean point, for the outer part (big |z|) of the ring
     * @param rmin The global minimum radius of the ring, as measured from the z-axis
     * @param rmax The global maximum radius of the ring, as measured from the z-axis
     * @param zmin The minimum z of the ring
     * @param zmax The maximum z of the ring
     * @param rminatzmin The minimum radius of the ring at z=zmin, as measured from the z-axis
     * @param rmaxatzmax The maximum radius of the ring at z=zmax, as measured from the z-axis
     */
    struct BTiltedRingInfo {
        std::string name;
        std::string childname;
        bool isZPlus;
        double tiltAngle;      // in DEG
        bool bw_flipped;
        bool fw_flipped;  
        int modules;
        double r1;
        double z1;
        double startPhiAngle1; // in RAD
        double r2;
        double z2;
        double startPhiAngle2; // in RAD
        double rmin;
        double rmax;
        double zmin;
        double zmax;
        double rminatzmin;
        double rmaxatzmax;     
    };
     /*
      * * @struct ModuleROCInfo
      * * @brief The information in this struct is a parameter in SpecParInfo.
      * * @param name The module type
      * * @param rocrows The number of ROCRows
      * * @param roccols The number of ROCCols
      * * @param rocx The number of ROC_X
      * * @param rocy The number of ROC_Y
      * */
    struct ModuleROCInfo {
      std::string name;
      std::string rocrows;
      std::string roccols;
      std::string rocx;
      std::string rocy;
    };

    /**
     * @struct SpecParInfo
     * @brief The information in this struct translates to a <i>SpecPar</i> block in CMSSW XML.
     * @param name The name of the <i>SpecPar</i> block
     * @param parameter The property within CMSSW that the block refers to
     * @param partselectors A list of logical part names that have the property named in <i>parameter</i>
     */
    struct SpecParInfo {
        std::string name;
        std::pair<std::string, std::string> parameter;
        std::vector<std::string> partselectors;
        std::vector<std::string> partextras; // CUIDADO this was introduced to add the Plus / Minus attribute to discs, not really useful for anything else at least for now
        //std::vector<std::string> moduletypes;
        std::vector<ModuleROCInfo> moduletypes;
    };

    /**
     * @struct RILengthInfo
     * @brief 
     * @param barrel 
     * @param index 
     * @param rlength 
     * @param ilength 
     */
    struct RILengthInfo {
        bool barrel;
        int index;
        double rlength;
        double ilength;
    };
    /**
     * @struct PathInfo 
     * @brief 
     * @param 
     */
    struct PathInfo {
        std::string block_name;
        int layer;
        bool barrel;
        std::vector<std::string> paths;
    };
    /**
     * @struct CMSSWBundle
     * @brief 
     * @param elements 
     * @param composites 
     * @param logic 
     * @param shapes 
     * @param positions 
     * @param algos 
     * @param rots 
     * @param specs 
     */
    struct CMSSWBundle {
      std::vector<Element> elements;
      std::vector<Composite> composites;
      std::vector<LogicalInfo> logic;
      std::vector<ShapeInfo> shapes;
      std::vector<ShapeOperationInfo> shapeOps;
      std::vector<PosInfo> positions;
      std::vector<AlgoInfo> algos;
      std::map<std::string,Rotation> rots;
      std::vector<SpecParInfo> specs;
      std::vector<RILengthInfo> lrilength;
    };
}
#endif /* _TK2CMSSW_DATATYPES_H */

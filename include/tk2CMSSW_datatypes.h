#ifndef _TK2CMSSW_DATATYPES_H
#define _TK2CMSSW_DATATYPES_H

/**
 * @file tk2CMSSW_datatypes.h
 * @brief This header file lists the internal structs and enumerations needed to translate a full tracker and material budget to CMSSW XML
 */

#include <string>
#include <vector>

namespace insur {
    /**
     * @enum CompType A list of the options for composite materials in CMSSW: mixture by weight, by volume or by atomic proportion
     */
    enum CompType { wt, vl, ap };
    /**
     * @enum ShapeType A list of possible shape types: box, tubs, tdr1 (trapezoid) and polycone
     */
    enum ShapeType { bx, tb, tp, pc };
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
     * @param elements A list of pairs indicating the elements in the compound by name and fraction of the whole
     */
    struct Composite {
        std::string name;
        double density;
        CompType method;
        std::vector<std::pair<std::string, double> > elements;
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
     * @param type The shape type: one of box, tube section, trapezoid or polycone
     * @param name_tag The name of the shape volume
     * @param dx Half the width along the x-axis (boxes and trapezoids)
     * @param dy Half the height along the y-axis at the base (boxes and trapezoids)
     * @param dyy Half the height along the y-axis at the top (trapezoids only)
     * @param dz Half the depth along the z-axis (boxes, trapezoids and tube sections)
     * @param rmin The minimal radius as measured from the z-axis (tube sections only) 
     * @param rmax The maximum radius as measured from the z-axis (tube sections only)
     * @param rzup A vector of pairs collecting the first half of an ordered sequence of points in <i>(r, z)</i> (polycone only)
     * @param rzdown A vector of pairs collecting the second half of an ordered sequence of points in <i>(r, z)</i> (polycone only)
     */
    struct ShapeInfo {
        ShapeType type;
        std::string name_tag;
        double dx;
        double dy;
        double dyy;
        double dz;
        double rmin;
        double rmax;
        std::vector<std::pair<double, double> > rzup;
        std::vector<std::pair<double, double> > rzdown;
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
    };
    /**
     * @struct RingInfo
     * @brief This is a struct to collect temporary information about an endcap ring and the modules within it.
     * @param name The logical part name that identifies the ring
     * @param childname The logical part name that identifies the modules contained in the ring
     * @param fw A flag indicating where in z the ring is: true if it is in z+, false if it is in z-
     * @param modules The number of modules within the ring
     * @param rin The minimum radius of the ring as measured from the z-axis
     * @param rout The maximum radius of the ring as measured from the z-axis
     * @param rmid The radius of the module mean point as measured from the z-axis
     * @param mthk The total thickness of one of the ring modules
     * @param phi The angle <i>phi</i> in the x/y-plane of the first module on the ring
     */
    struct RingInfo {
        std::string name;
        std::string childname;
        bool fw;
        int modules;
        double rin;
        double rout;
        double rmid;
        double mthk;
        double phi;
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
        std::vector<std::string> moduletypes;
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
        std::vector<PosInfo> positions;
        std::vector<AlgoInfo> algos;
        std::vector<Rotation> rots;
        std::vector<SpecParInfo> specs;
        std::vector<RILengthInfo> lrilength;
    };
}
#endif /* _TK2CMSSW_DATATYPES_H */

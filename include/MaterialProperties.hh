//
// File:   MaterialProperties.h
// Author: ndemaio
//
// Created on January 30, 2009, 3:50 PM
//

/**
 * @file MaterialProperties.h
 * @brief This is the header of a base class for elements that have an influence on the material budget of a tracker layout
 */

#ifndef _MATERIALPROPERTIES_H
#define	_MATERIALPROPERTIES_H

#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <MaterialTable.hh>

class RILength {
public:
  RILength() {radiation=0; interaction=0;}
  double radiation;
  double interaction;
  RILength& operator+=(const RILength &a);
  const RILength operator+(const RILength &other) const;
};

typedef RILength Material;


namespace insur {
    /**
     * Errors and messages that may be reported during operations on member variables
     */
    static const std::string err_local_mass = "Local mass not found";
    static const std::string msg_mattab_except_local = "Exception other than runtime_error occurred accessing material table for local masses: ";
    /**
     * @class MaterialProperties
     * @brief This is the base class for collections of properties related to the material budget.
     *
     * It encapsulates the main parameters of interest, namely the overall density, radiation length
     * and interaction length of a tracker element, as well as the influence of the materials that make
     * up this building block. Access functions are provided where appropriate. But unless the object
     * is cloned the overall parameters should typically be calculated from a list of materials and their
     *  properties rather than set explicitly. Some of the access functions for individual materials may
     * throw exceptions if the requested material does not appear on the list.
     */
    class MaterialProperties {
    public:
        /**
         * @enum Category A list of logical categories within the detector geometry; a single element belongs to exactly one of them
         */

        enum Category {no_cat, b_mod, e_mod, b_ser, e_ser, b_sup, e_sup, o_sup, t_sup, u_sup};
        MaterialProperties();
        //virtual ~MaterialProperties() {} 
        // bureaucracy
        Category getCategory();
        void setCategory(Category c);
        // to be used by the subclasses
        virtual double getSurface() const;
        virtual double getLength() const;
        // material mass handling
        const std::map<std::string, double>& getLocalMasses() const;
        const std::map<std::string, double>& getLocalMassesComp() const;
        double getLocalMass(std::string tag); // throws exception
        double getLocalMassComp(std::string tag); // throws exception
        void addLocalMass(std::string tag, std::string comp, double ms, int minZ = -777);
        void addLocalMass(std::string tag, double ms);
        unsigned int localMassCount();
        unsigned int localMassCompCount();
        void clearMassVectors();
        void copyMassVectors(MaterialProperties& mp);
        // calculated output values
        double getTotalMass() const;
        double getLocalMass();
        double getRadiationLength();
        double getInteractionLength();
        RILength getMaterialLengths();
        const std::map<std::string, RILength>& getComponentsRI() const;
        // output calculations
        void calculateTotalMass(double offset = 0);
        void calculateLocalMass(double offset = 0);
        void calculateRadiationLength(MaterialTable& materials, double offset = 0.0);
        void calculateInteractionLength(MaterialTable& materials,  double offset = 0.0);
        void calculateRadiationLength(double offset = 0.0);
        void calculateInteractionLength(double offset = 0.0);
        // tracking information
        bool track();
        void track(bool tracking_on);
        // base print function
        void print();

    protected:
        // init flags and tracking
        bool msl_set, trck;
        // geometry-dependent parameters
        Category cat;
        std::map<std::string, double> localmasses;
        std::map<std::string, double> localmassesComp;

        std::map<std::string, std::map<std::string, double> > localCompMats; // format here is <component name string, <material name, mass> >

        std::map<std::string, RILength> componentsRI;  // component-by-component radiation and interaction lengths
        // complex parameters (OUTPUT)
        double total_mass, local_mass, r_length, i_length;
        // internal help
        std::string getSuperName(std::string name) const;
        std::string getSubName(std::string name) const;
    };
}
#endif	/* _MATERIALPROPERTIES_H */


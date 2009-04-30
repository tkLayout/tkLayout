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
#include <MaterialTable.h>
namespace insur {
    static const std::string err_local_mass = "Local mass not found";
    static const std::string err_exiting_mass = "Exiting mass not found";
    /**
     * @class MaterialProperties
     * @brief This is the base class for collections of properties related to the material budget.
     *
     * It encapsulates the main parameters of interest, namely the overall density, radiation length
     * and interaction length of a tracker element, as well as the influence of the materials that make
     * up this building block. Access functions are provided where appropriate. But unless the object
     * is cloned the overall parameters should typically be calculated from a list of materials and their
     *  properties rather than set explicitly.
     */
    class MaterialProperties {
    public:
        enum Category {no_cat, b_mod, e_mod, b_ser, e_ser, b_sup, e_sup, t_sup, u_sup};
        MaterialProperties();
        virtual ~MaterialProperties() {}
        // bureaucracy
        Category getCategory();
        void setCategory(Category c);
        // to be used by the subclasses
        virtual double getSurface();
        // material mass handling
        double getLocalMass(std::string tag); // throws exception
        double getLocalMass(int index); // throws exception
        std::string getLocalTag(int index);
        double getExitingMass(std::string tag); // throws exception
        double getExitingMass(int index); // throws exception
        std::string getExitingTag(int index);
        void setLocalMass(std::string tag, double ms);
        void addLocalMass(std::string tag, double ms);
        void setExitingMass(std::string tag, double ms);
        void addExitingMass(std::string tag, double ms);
        uint localMassCount();
        uint exitingMassCount();
        void clearMassVectors();
        void copyMassVectors(MaterialProperties& mp);
        // calculated output values
        double getTotalMass();
        double getLocalMass();
        double getExitingMass();
        double getRadiationLength();
        double getInteractionLength();
        // output calculations
        void calculateTotalMass(double offset = 0);
        void calculateLocalMass(double offset = 0);
        void calculateExitingMass(double offset = 0);
        void calculateRadiationLength(MaterialTable& materials, double offset = 0);
        void calculateInteractionLength(MaterialTable& materials,  double offset = 0);
        // base print function
        void print();
    protected:
        // init flag
        bool msl_set, mse_set;
        // geometry-dependent parameters
        Category cat;
        std::vector<std::pair<std::string, double> > localmasses, exitingmasses;
        // complex parameters (OUTPUT)
        double total_mass, local_mass, exiting_mass, r_length, i_length;
        // internal help
        void setLocalMass(std::pair<std::string, double> ms);
        void addLocalMass(std::pair<std::string, double> ms);
        void setExitingMass(std::pair<std::string, double> ms);
        void addExitingMass(std::pair<std::string, double> ms);
        int findLocalIndex(std::string tag);
        int findExitingIndex(std::string tag);
        bool newLocalMaterial(std::string tag);
        bool newExitingMaterial(std::string tag);
    };
}
#endif	/* _MATERIALPROPERTIES_H */


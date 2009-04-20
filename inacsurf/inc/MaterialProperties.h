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
        MaterialProperties() { ms_set = false; total_mass = -1; local_mass = -1; exiting_mass = -1; r_length = -1; i_length = -1; }
        virtual ~MaterialProperties() {}
        // to be used by the subclasses
        virtual double getSurface();
        // material mass handling
        double getLocalMass(std::string tag); // throws exception
        double getLocalMass(int index); // throws exception
        double getExitingMass(std::string tag); // throws exception
        double getExitingMass(int index); // throws exception
        void setLocalMass(std::string tag, double tk);
        void addLocalMass(std::string tag, double tk);
        void setExitingMass(std::string tag, double tk);
        void addExitingMass(std::string tag, double tk);
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
        bool ms_set;
        // geometry-dependent parameters
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


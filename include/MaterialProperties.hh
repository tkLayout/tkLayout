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

  enum MechanicalCategory { UNKNOWN, MODULE, CABLING, SUPPORT_AND_COOLING };

  static const std::string mechanical_module = "Module:";
  static const std::string mechanical_cabling = "Cabling:";
  static const std::string mechanical_cooling = "Cooling";
  static const std::string mechanical_support = "Supports Mechanics:";

  static const MechanicalCategory computeMechanicalCategory(const std::string componentName) {
    if (componentName.find(mechanical_module) != std::string::npos) return MechanicalCategory::MODULE;
    else if (componentName.find(mechanical_cabling) != std::string::npos) return MechanicalCategory::CABLING;
    else if (componentName.find(mechanical_support) != std::string::npos) return MechanicalCategory::SUPPORT_AND_COOLING;
    else if (componentName.find(mechanical_cooling) != std::string::npos) return MechanicalCategory::SUPPORT_AND_COOLING;
    else return MechanicalCategory::UNKNOWN;
  }

  /* Normalize the mechanical categories' contributions to RL and IL.
   */
  static void normalizeRIRatio(std::map<MechanicalCategory, std::pair<double, double> >& normalizedRIRatioPerMechanicalCategory) {
    double nonNormalizedRadiationLengthSum = 0.;
    double nonNormalizedInteractionLengthSum = 0.;

    // GET THE SUM OF THE MECHANICAL CATEGORIES' RI CONTRIBUTIONS
    for (const auto& mechanicalCategoryIt : normalizedRIRatioPerMechanicalCategory) {
      const std::pair<double, double>& nonNormalizedRI = mechanicalCategoryIt.second;

      // Sum of contributions to RL
      const double nonNormalizedRadiationLength = nonNormalizedRI.first;
      nonNormalizedRadiationLengthSum += nonNormalizedRadiationLength;

      // Sum of contributions to IL
      const double nonNormalizedInteractionLength = nonNormalizedRI.second;
      nonNormalizedInteractionLengthSum += nonNormalizedInteractionLength;
    }

    // NORMALIZE
    for (auto& mechanicalCategoryIt : normalizedRIRatioPerMechanicalCategory) {
      std::pair<double, double>& nonNormalizedRI = mechanicalCategoryIt.second;

      // Normalize each RL contribution
      double& nonNormalizedRadiationLength = nonNormalizedRI.first;
      nonNormalizedRadiationLength /= nonNormalizedRadiationLengthSum;

      // Normalize each IL contribution
      double& nonNormalizedInteractionLength = nonNormalizedRI.second;	    
      nonNormalizedInteractionLength /= nonNormalizedInteractionLengthSum;
    }
  }



  class LocalElement {
  public:
    LocalElement(const std::string subdetectorName, const std::string componentName, const std::string elementName) :
      subdetectorName_(subdetectorName),
      componentName_(componentName),
      elementName_(elementName) 
    { 
      category_ = computeMechanicalCategory(componentName);
    }

    const std::string subdetectorName() const { return subdetectorName_; }
    const MechanicalCategory mechanicalCategory() const { return category_; }
    const std::string componentName() const { return componentName_; }
    const std::string elementName() const { return elementName_; }

 
  protected:
    std::string subdetectorName_;
    MechanicalCategory category_;
    std::string componentName_;
    std::string elementName_;
  };


  struct ElementNameCompare {
    bool operator() (const LocalElement& localA, const LocalElement& localB) const {
      if (localA.subdetectorName() != localB.subdetectorName()) return (localA.subdetectorName() < localB.subdetectorName());
      else {
	if (localA.mechanicalCategory() != localB.mechanicalCategory()) return (localA.mechanicalCategory() < localB.mechanicalCategory());
	else {
	  if (localA.componentName() != localB.componentName()) return (localA.componentName() < localB.componentName());
	  else return (localA.elementName() < localB.elementName());
	}
      }
    }
    // NB: No need of operator==
    // (a == b)  <=>  ( !(a<b) && !(b<a) )
  };

  struct ComponentNameCompare {
    bool operator() (const LocalElement& localA, const LocalElement& localB) const {
      if (localA.subdetectorName() != localB.subdetectorName()) return (localA.subdetectorName() < localB.subdetectorName());
      else {
	if (localA.mechanicalCategory() != localB.mechanicalCategory()) return (localA.mechanicalCategory() < localB.mechanicalCategory());
	else {
	  return (localA.componentName() < localB.componentName());
	}
      }
    }
    // NB: No need of operator==
    // (a == b)  <=>  ( !(a<b) && !(b<a) )
  };

  



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
        virtual ~MaterialProperties() {} 
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

      // weights
      const std::map<LocalElement, double, ElementNameCompare> getLocalElementsDetails() const { return localMassesDetails_; }
      // RI
      const std::map<LocalElement, RILength, ComponentNameCompare>& getComponentsRI() const { return componentsRI_; }
      // Return the mechanical categories' normalized contributions to RI.
      const std::map<MechanicalCategory, std::pair<double, double> >& getNormalizedRIRatioPerMechanicalCategory() {
	// Normalize first.
	insur::normalizeRIRatio(normalizedRIRatioPerMechanicalCategory_);
	return normalizedRIRatioPerMechanicalCategory_; 
      }

        void addLocalMass(const std::string subdetectorName, const std::string tag, const std::string comp, double ms, int minZ = -777);
        void addLocalMass(const std::string subdetectorName, const std::string tag, double ms);
      //void addLocalMass(std::string tag, std::string comp, double ms, int minZ = -777);
      //void addLocalMass(std::string tag, double ms);
        unsigned int localMassCount();
        unsigned int localMassCompCount();
        void clearMassVectors();
        void copyMassVectors(MaterialProperties& mp);
        // calculated output values
        double getTotalMass() const;
        double getLocalMass();
      const double getMechanicalModuleWeight() const;
        double getRadiationLength();
        double getInteractionLength();
        RILength getMaterialLengths();
       
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

        // THIS SHOULD REPLACE localmasses, localmassesComp, and so on. All desired info is accessed from LocalElementDetails:
      std::map<LocalElement, double, ElementNameCompare> localMassesDetails_;
      std::map<LocalElement, RILength, ComponentNameCompare> componentsRI_;  // component-by-component radiation and interaction lengths
      std::map<MechanicalCategory, std::pair<double, double> > normalizedRIRatioPerMechanicalCategory_;

        std::map<std::string, double> localmassesComp;
        std::map<std::string, std::map<std::string, double> > localCompMats; // format here is <component name string, <material name, mass> >

        // complex parameters (OUTPUT)
        double total_mass, local_mass, r_length, i_length;
        // internal help
        std::string getSuperName(std::string name) const;
        std::string getSubName(std::string name) const;
    };
}
#endif	/* _MATERIALPROPERTIES_H */


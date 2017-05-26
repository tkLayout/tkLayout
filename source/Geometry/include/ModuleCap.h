// 
// File:   ModuleCap.h
// Author: ndemaio
//
// Created on January 30, 2009, 2:44 PM
//

/**
 * @file ModuleCap.h
 * @brief This is the add-on class that collects material properties for a module
 */

#ifndef INCLUDE_MODULECAP_H_
#define INCLUDE_MODULECAP_H_

#include "MaterialProperties.h"

// Fwd declaration
class DetectorModule;

/**
* @class ModuleCap
* @brief This is the header of a class associated to an instance of <i>Module</i> and containing its material parameters.
*
* The default constructor without arguments is private because an instance of this class only makes sense
* in connection with an existing active surface. The only way to instantiate this class is by giving it a reference
* to a previously created module. The combination will then include the geometry information in the <i>Module</i>
* class and the information to calculate the material budget in this one.
*/
class ModuleCap : public MaterialProperties {

public:

  //! Constructer
  ModuleCap(DetectorModule& mod);

  //! Destructor
  virtual ~ModuleCap();

  //! Get module which the material cap relates to
  DetectorModule& getModule();

  //! Get module surface
  virtual double getSurface() const;

  //! Get module lenght
  virtual double getLength() const;

  //! Print module properties
  virtual void print();

protected:
  DetectorModule& module; //!< Reference to a module

}; // Class

#endif	/* INCLUDE_MODULECAP_H_ */


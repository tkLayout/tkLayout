/**
 * @file ModuleCap.cc
 * @brief This class bundles the material properties of an active module.
 */

#include "ModuleCap.h"

#include "DetectorModule.h"
#include <global_constants.h>

// public
/**
 * The constructor links the cap to an existing module.
 */
ModuleCap::ModuleCap(DetectorModule& mod) : module(mod) {
    mod.setModuleCap(this);
}

/**
 * Nothing to do for the destructor...
 */
ModuleCap::~ModuleCap() {}

/**
 * Get the module object that this cap refers to.
 * @return A reference to the registered module
 */
DetectorModule& ModuleCap::getModule() {
    return module;
}

/**
 * Get the module surface.
 * @return The module surface that is relevant for tracking
 */
double ModuleCap::getSurface() const {
    return module.area();
}

/**
 * Get the module surface.
 * @return The module surface that is relevant for tracking
 */
double ModuleCap::getLength() const {
    return module.length();
}


/**
 * This function prints a summary of the object to <i>cout</i>
 */
void ModuleCap::print() {
      MaterialProperties::print();
      std::cout << "Associated module is of type " << module.moduleType() << "." << std::endl;
  }

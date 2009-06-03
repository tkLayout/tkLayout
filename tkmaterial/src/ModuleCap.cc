/**
 * @file ModuleCap.cc
 * @brief This class bundles the material properties of an active module.
 */

#include <ModuleCap.h>
namespace insur {
    // public
    ModuleCap::ModuleCap(Module& mod) { 
        module = &mod;
    }
    
    ModuleCap::~ModuleCap() {}
    
    Module& ModuleCap::getModule() {
        return *module;
    }
    
    double ModuleCap::getSurface() {
        return module->getArea();
    }
    
    void ModuleCap::print() {
        MaterialProperties::print();
        // TODO: finish
    }
}

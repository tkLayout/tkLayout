/**
 * @file ModuleCap.cc
 * @brief This class bundles the material properties of an active module.
 */

#include <ModuleCap.hh>

namespace insur {

    // public
    /**
     * The constructor links the cap to an existing module.
     */
    ModuleCap::ModuleCap(Module& mod, int id) {
        m_module = &mod;
        mod.setModuleCap(this);
        m_layerOrDiscID = id;
        m_detName = "Undefined";
    }
    
    /**
     * Nothing to do for the destructor...
     */
    ModuleCap::~ModuleCap() {}
    
    /**
     * Get the module object that this cap refers to.
     * @return A reference to the registered module
     */
    Module& ModuleCap::getModule() {
        return *m_module;
    }
    
    /**
     * Get the module surface.
     * @return The module surface that is relevant for tracking
     */
    double ModuleCap::getSurface() const {
        return m_module->area();
    }

    /**
     * Get the module surface.
     * @return The module surface that is relevant for tracking
     */
    double ModuleCap::getLength() const {
        return m_module->length();
    }

    
    /**
     * This function prints a summary of the object to <i>cout</i>
     */
    void ModuleCap::print() {
        MaterialProperties::print();
        std::cout << "Associated module is of type " << m_module->moduleType() << "." << std::endl;
    }
}

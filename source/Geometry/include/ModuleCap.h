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

#ifndef _MODULECAP_H
#define	_MODULECAP_H

//#include <module.hh>
#include "Module.h"
#include <MaterialProperties.h>
#include <global_constants.h>

namespace insur {
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
        ModuleCap(Module& mod);
        virtual ~ModuleCap();
        Module& getModule();
        virtual double getSurface() const;
        virtual double getLength() const;
        virtual void print();
    protected:
        Module* module;
    };
}
#endif	/* _MODULECAP_H */


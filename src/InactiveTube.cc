/**
 * @file InactiveTube.cc
 * @brief This is the derived class implementation of a single tube-shaped inactive element
 */

#include <InactiveTube.hh>
namespace insur {
    /**
     * The constructors make sure the UP/DOWN flag is set correctly (it is provided by the parent class, but not set).
     */
    InactiveTube::InactiveTube() : InactiveElement() {
        is_vertical = false;
    }
    
    InactiveTube::InactiveTube(InactiveElement& previous) : InactiveElement(previous) {
        is_vertical = false;
    }
    
     /**
     * Nothing to do for the destructor...
     */
    InactiveTube::~InactiveTube() {}
    
    /**
     * This function prints a representation of the element to <i>cout</i>
     */
    void InactiveTube::print() {
        InactiveElement::print();
        if (sanityCheck()) std::cout << "Volume is sane." << std::endl;
        else std::cout << "WARNING: Volume is not sane!" << std::endl;
    }
    
     /**
     * The sanity check virtual function tests if the object has the geometric properties of a tube.
     * @return True if the length is greater than or equal to the width, false otherwise
     */
    bool InactiveTube::sanityCheck() {
        return (z_length >= w_radius) && !is_vertical;
    }
}

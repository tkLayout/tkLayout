/**
 * @file InactiveElement.cc
 * @brief This is the base class implementation of a single inactive element
 */

#include <InactiveElement.h>
namespace insur {
    /*-----public functions-----*/
    
    InactiveElement::InactiveElement() {
        is_final = false;
        feeder_type = no_in;
        feeder_index = -1;
        neighbour_type = no_in;
        neighbour_index = -1;
    }
    
    double InactiveElement::getSurface() {
        if (isVertical()) return ((i_radius + w_radius) * (i_radius + w_radius) - i_radius * i_radius) * PI;
        else return 2 * PI * (i_radius + w_radius / 2.0) * z_length;
    }
    
    /**
     * Query the orientation of the object.
     * @return True if the element points up or down, false if it points sideways
     */
    bool InactiveElement::isVertical() { return is_vertical; }
    
    /**
     * Set the orientation flag of the object.
     * @param vertical The new value for the up/down flag
     */
    void InactiveElement::setVertical(bool vertical) { is_vertical = vertical; }
    
    bool InactiveElement::isFinal() { return is_final; }
    
    void InactiveElement::setFinal(bool final) { is_final = final; }
    
    /**
     * Get the distance of this object's leftmost point to the xy-plane.
     * @return The offset from the origin along the z-axis
     */
    double InactiveElement::getZOffset() { return z_offset; }
    
    /**
     * Set the distance of this object's leftmost point to the xy-plane.
     * @param zoffset The offset from the origin along the z-axis
     */
    void InactiveElement::setZOffset(double zoffset) { z_offset = zoffset; }
    
    /**
     * Get the length of the element.
     * @return The total length of the element along the z-axis
     */
    double InactiveElement::getZLength() { return z_length; }
    
    /**
     * Set the length of the element.
     * @param zlength The total length of the element along the z-axis
     */
    void InactiveElement::setZLength(double zlength) { z_length = zlength; }
    
    /**
     * Get the inner radius of the element.
     * @return The distance from the z-axis to the innermost point of the element
     */
    double InactiveElement::getInnerRadius() { return i_radius; }
    
    /**
     * Set the inner radius of the element.
     * @param iradius The distance from the z-axis to the innermost point of the element
     */
    void InactiveElement::setInnerRadius(double iradius) { i_radius = iradius; }
    
    /**
     * Get the width of the element.
     * @return The distance from the innermost to the outermost point of the element in the xy-plane
     */
    double InactiveElement::getRWidth() { return w_radius; }
    
    /**
     * Set the width of the element.
     * @param rwidth The distance from the innermost to the outermost point of the element in the xy-plane
     */
    void InactiveElement::setRWidth(double rwidth) { w_radius = rwidth; }
    
    int InactiveElement::getFeederIndex() { return feeder_index; }
    
    void InactiveElement::setFeederIndex(int layer) { feeder_index = layer; }
    
    InactiveElement::InType InactiveElement::getFeederType() { return feeder_type; }
    
    void InactiveElement::setFeederType(InType type) { feeder_type = type; }
    
    int InactiveElement::getNeighbourIndex() { return neighbour_index; }
    
    void InactiveElement::setNeighbourIndex(int previous) { neighbour_index = previous; }
    
    InactiveElement::InType InactiveElement::getNeighbourType() { return neighbour_type; }
    
    void InactiveElement::setNeighbourType(InactiveElement::InType type) { neighbour_type = type; }
    
    /**
     * Set the total mass of the inactive element.
     * @param mass The new overall mass
     */
    void InactiveElement::setTotalMass(double mass) { total_mass = mass; }
    
    /**
     * Set the total mass of the inactive element.
     * @param mass The new overall mass
     */
    void InactiveElement::setLocalMass(double mass) { local_mass = mass; }
    
    /**
     * Set the total mass of the inactive element.
     * @param mass The new overall mass
     */
    void InactiveElement::setExitingMass(double mass) { exiting_mass = mass; }
    
    /**
     * Set the radiation length of the inactive element.
     * @param rlength The new overall radiation length, averaged over all the different material that occur in the inactive element
     */
    void InactiveElement::setRadiationLength(double rlength) { r_length = rlength; }
    
    /**
     * Set the interaction length of the inactive element.
     * @param ilength The new overall interaction length, averaged over all the different material that occur in the inactive element
     */
    void InactiveElement::setInteractionLength(double ilength) { i_length = ilength; }
    
    /**
     * Calculate and return the Eta range of the element
     * @return The pair <i>(Eta_min, Eta_max)</i>
     */
    std::pair<double, double> InactiveElement::getEtaMinMax() {
        std::pair<double, double> res;
        double theta0, theta1;
        if ((getZOffset() < 0) && (getZOffset() + getZLength() > 0)) {
            theta0 = atan(getInnerRadius() / (-1 * getZOffset()));
            theta0 = PI - theta0;
            theta1 = atan(getInnerRadius() / (getZOffset() + getZLength()));
        }
        else {
            if (isVertical()) {
                theta0 = atan((getInnerRadius() + getRWidth()) / (getZLength() / 2.0 + getZOffset()));
                theta1 = atan(getInnerRadius() / (getZLength() / 2.0 + getZOffset()));
            }
            else {
                theta0 = atan((getRWidth() / 2.0 + getInnerRadius()) / getZOffset());
                theta1 = atan((getRWidth() / 2.0 + getInnerRadius()) / (getZOffset() + getZLength()));
            }
        }
        res.first = -1 * log(tan(theta0 / 2.0));
        res.second = -1 * log(tan(theta1 / 2.0));
        return res;
    }
    
    /**
     * Print the geometry-specific parameters of the inactive element including the orientation.
     */
    void InactiveElement::print() {
        MaterialProperties::print();
        std::cout << "Inactive element properties (current state)" << std::endl;
        std::cout << "z_offset = " << z_offset << std::endl;
        std::cout << "z_length = " << z_length << std::endl;
        std::cout << "i_radius = " <<i_radius  << std::endl;
        std::cout << "w_radius = " << w_radius << std::endl;
        if (isVertical()) std::cout << "Volume is vertical." << std::endl;
        else std::cout << "Volume is horizontal" << std::endl;
    }
}

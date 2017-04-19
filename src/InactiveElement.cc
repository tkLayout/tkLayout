#include <cmath>
#include <iostream>

/**
 * @file InactiveElement.cc
 * @brief This is the base class implementation of a single inactive element
 */

#include <InactiveElement.hh>
namespace insur {
    /*-----public functions-----*/
    /**
     * The constructor sets some defaults: no neighbours, intermediate element.
     */
    InactiveElement::InactiveElement() {
        is_final = false;
        feeder_type = no_in;
        feeder_index = -1;
        neighbour_type = no_in;
        neighbour_index = -1;
    }
    
    /**
     * Get the surface of the element that is relevant for material budget calculation.
     * @return The average cylinder surface for a horizontal tube, the disc surface for a vertical disc
     */
    double InactiveElement::getSurface() const {
        if (isVertical()) return ((i_radius + w_radius) * (i_radius + w_radius) - i_radius * i_radius) * M_PI;
        else return 2 * M_PI * (i_radius + w_radius / 2.0) * z_length;
    }
    
    /**
     * Get the orientation of the object.
     * @return True if the element points up or down, false if it points sideways
     */
    bool InactiveElement::isVertical() const { return is_vertical; }
    
    /**
     * Set the orientation flag of the object.
     * @param vertical The new value for the up/down flag
     */
    void InactiveElement::setVertical(bool vertical) { is_vertical = vertical; }
    
    /**
     * Check if the content of this element travels out of the tracking volume after this.
     * @return True if the element is not a neighbour to anything, false otherwise
     */
    bool InactiveElement::isFinal() { return is_final; }
    
    /**
     * Set if the content of this element travels out of the tracking volume after this.
     * @param final Should be true if the element is not a neighbour to anything and false otherwise
     */
    void InactiveElement::setFinal(bool final) { is_final = final; }
    
    /**
     * Get the distance of this object's leftmost point to the xy-plane.
     * @return The offset from the origin along the z-axis
     */
    double InactiveElement::getZOffset() const { return z_offset; }
    
    /**
     * Set the distance of this object's leftmost point to the xy-plane.
     * @param zoffset The offset from the origin along the z-axis
     */
    void InactiveElement::setZOffset(double zoffset) { z_offset = zoffset; }
    
    /**
     * Get the length of the element.
     * @return The total length of the element along the z-axis
     */
    double InactiveElement::getZLength() const { return z_length; }
    
    /**
     * Set the length of the element.
     * @param zlength The total length of the element along the z-axis
     */
    void InactiveElement::setZLength(double zlength) { z_length = zlength; }
    
    /**
     * Get the inner radius of the element.
     * @return The distance from the z-axis to the innermost point of the element
     */
    double InactiveElement::getInnerRadius() const { return i_radius; }
    
    /**
     * Set the inner radius of the element.
     * @param iradius The distance from the z-axis to the innermost point of the element
     */
    void InactiveElement::setInnerRadius(double iradius) { i_radius = iradius; }
    
    /**
     * Get the width of the element.
     * @return The distance from the innermost to the outermost point of the element in the xy-plane
     */
    double InactiveElement::getRWidth() const { return w_radius; }
    
    /**
     * Set the width of the element.
     * @param rwidth The distance from the innermost to the outermost point of the element in the xy-plane
     */
    void InactiveElement::setRWidth(double rwidth) { w_radius = rwidth; }
    
    double InactiveElement::getLength() const {
      if(is_vertical) {
        return getRWidth();
      } else {
        return getZLength();
      }
    }

    double InactiveElement::getVolume() const {
      return M_PI * w_radius * (w_radius + 2*i_radius) * z_length;
    }

    /**
     * Get the index of the element's feeder volume.
     * @return The index within the tracker object's layer or disc vector, or of the service volume; -1 if there is none
     */
    int InactiveElement::getFeederIndex() { return feeder_index; }
    
    /**
     * Set the index of the element's feeder volume.
     * @param layer The index within the tracker object's layer or disc vector, or of the service volume
     */
    void InactiveElement::setFeederIndex(int layer) { feeder_index = layer; }
    
    /**
     * Get the type of the element's feeder volume.
     * @return The type of feeder as listed in the enumeration <i>InType</i>
     */
    InactiveElement::InType InactiveElement::getFeederType() { return feeder_type; }
    
    /**
     * Set the type of the element's feeder volume.
     * @param type The type of feeder as listed in the enumeration <i>InType</i>
     */
    void InactiveElement::setFeederType(InType type) { feeder_type = type; }
    
    /**
     * Get the index of the element's neighbour volume.
     * @return The index of the previous service volume
     */
    int InactiveElement::getNeighbourIndex() { return neighbour_index; }
    
    /**
     * Set the index of the element's neighbour volume.
     * @param previous The index of the previous service volume
     */
    void InactiveElement::setNeighbourIndex(int previous) { neighbour_index = previous; }
    
    /**
     * Get the type of the element's neighbour volume.
     * @return The type of neighbour as listed in the enumeration <i>InType</i>
     */
    InactiveElement::InType InactiveElement::getNeighbourType() { return neighbour_type; }
    
    /**
     * Set the type of the element-s neighbour volume.
     * @param type The type of neighbour as listed in the enumeration <i>InType</i>
     */
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
        // volumes crossing z=0
        if ((getZOffset() < 0) && (getZOffset() + getZLength() > 0)) {
            // lower left of tube wall above z-axis
            theta0 = atan(getInnerRadius() / (-1 * getZOffset()));
            theta0 = M_PI - theta0;
            // lower right of tube wall above z-axis
            theta1 = atan(getInnerRadius() / (getZOffset() + getZLength()));
        }
        // volumes on either side of z=0
        else {
            // rings
            if (isVertical()) {
                // upper centre of tube wall above z-axis
                theta0 = atan((getInnerRadius() + getRWidth()) / (getZLength() / 2.0 + getZOffset()));
                // lower centre of tube wall above z-axis
                theta1 = atan(getInnerRadius() / (getZLength() / 2.0 + getZOffset()));
            }
            // tubes
            else {
                // centre left of tube wall above z-axis
                theta0 = atan((getRWidth() / 2.0 + getInnerRadius()) / getZOffset());
                // centre right of tube wall above z-axis
                theta1 = atan((getRWidth() / 2.0 + getInnerRadius()) / (getZOffset() + getZLength()));
            }
        }
        // convert angle theta to pseudorapidity eta
        res.first = -1 * log(tan(theta0 / 2.0));    // TODO change to the default converters
        res.second = -1 * log(tan(theta1 / 2.0));   // TODO change to the default converters
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

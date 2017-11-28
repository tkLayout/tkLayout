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
     * @return Max Z.
     */
    const double InactiveElement::getZMax() const { return z_offset + z_length; }
    
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

    /**
     * @return Max Rho.
     */
    const double InactiveElement::getOuterRadius() const { return i_radius + w_radius; }
    
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
   * Calculate and return the Eta range of the element, with respect to the origin.
   * @return The pair <i>(Eta_min, Eta_max)</i>
   */
    std::pair<double, double> InactiveElement::getEtaMinMax() const {
      std::pair<double, double> res;
      double theta0, theta1;

      // rings
      if (isVertical()) {
	// up point
	theta0 = atan((getInnerRadius() + getRWidth()) / (getZLength() / 2.0 + getZOffset()));
	// down point
	theta1 = atan(getInnerRadius() / (getZLength() / 2.0 + getZOffset()));
      }

      // tubes
      else {
	// left point
	theta0 = atan((getRWidth() / 2.0 + getInnerRadius()) / getZOffset());
	// right point
	theta1 = atan((getRWidth() / 2.0 + getInnerRadius()) / (getZOffset() + getZLength()));
      }

      // convert angle theta to pseudorapidity eta
      res.first = -1 * log(tan(theta0 / 2.0));    // TODO change to the default converters
      res.second = -1 * log(tan(theta1 / 2.0));   // TODO change to the default converters
      return res;
    }

  

  /**
   * Checks whether track hits the inactive element -> if yes, return true with passed hit position vector & material. 
   * Only implemeted for eta >= 0. (should work out of the box for eta < 0. but need to be tested in that case).
   * Assumes tracks are straight lines!
   * IP should be on (Z) axis, and track not parallel to (Z) axis.
   */
  const bool InactiveElement::checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir, XYZVector& hitPos, Material& hitMaterial) {

    // Initialize
    hitPos.SetX(0.);
    hitPos.SetY(0.);
    hitPos.SetZ(0.);
    hitMaterial.radiation   = 0.;
    hitMaterial.interaction = 0.;

    // Calculate hit only if inactive element non-transparent, otherwise no effect in tracking or material budget etc.
    if (getRadiationLength() != 0. || getInteractionLength() != 0.) {

      // eta >= 0.
      if (trackDir.Eta() >= 0.) {
	// IP should be on (Z) axis, and track not parallel to (Z) axis.
	if (trackOrig.Rho() == 0. && trackDir.Rho() > geom_zero) {	

	  // IP position + track direction parameters of interest for the calculation
	  const double trackOrigZ = trackOrig.Z();
	  const double trackDirZ = trackDir.Z();
	  const double trackDirRho = trackDir.Rho();

	  // Geometric extrema of the inactive volume in (RZ) plane
	  const double innerZ = getZOffset();
	  const double outerZ = getZMax();
	  const double innerRho = getInnerRadius();
	  const double outerRho = getOuterRadius();

	  // PART A
	  // COMPUTE THE PROJECTION OF THE TRACK IN THE 4 EXTREMA SURFACES OF THE INACTIVE VOLUME.
	  // Hit on inner (Z) plane
	  const double kInnerZPlane = (fabs(trackDirZ) > geom_zero ? (innerZ - trackOrigZ) / trackDirZ : -1.);
	  XYZVector hitOnInnerZPlane;
	  if (kInnerZPlane > 0.) hitOnInnerZPlane = trackOrig + kInnerZPlane * trackDir;

	  // Hit on outer (Z) plane
	  const double kOuterZPlane = (fabs(trackDirZ) > geom_zero ? (outerZ - trackOrigZ) / trackDirZ : -1.);
	  XYZVector hitOnOuterZPlane;
	  if (kOuterZPlane > 0.) hitOnOuterZPlane = trackOrig + kOuterZPlane * trackDir;

	  // Hit on inner radius cylinder
	  const double kInnerRhoCylinder = innerRho / trackDirRho;
	  const XYZVector hitOnInnerRhoCylinder = trackOrig + kInnerRhoCylinder * trackDir;

	  // Hit on outer radius cylinder
	  const double kOuterRhoCylinder = outerRho / trackDirRho;
	  const XYZVector hitOnOuterRhoCylinder = trackOrig + kOuterRhoCylinder * trackDir;
	

	  // PART B
	  // FIND OUT WHETHER THE TRACK CROSSES THE INSIDE OF THE VOLUME.
	  // IF SO, COMPUTES HIT POSITION AND PATH LENGTH.
	  bool hitFound = false;
	  double hitPathLength = 0.;  // length of the track path crossing the volume.

	  // Track enters and exits the volume between extrema radii cylinders.
	  if ( kInnerZPlane > 0. && kOuterZPlane > 0.
	       && (hitOnInnerZPlane.Rho() >= innerRho && hitOnInnerZPlane.Rho() <= outerRho)
	       && (hitOnOuterZPlane.Rho() >= innerRho && hitOnOuterZPlane.Rho() <= outerRho)
	       ) {
	    hitFound = true;
	    hitPos = (hitOnInnerZPlane + hitOnOuterZPlane) / 2.;
	    hitPathLength = (hitOnOuterZPlane - hitOnInnerZPlane).R();
	  }
	  // Top left corner.
	  else if (kInnerZPlane > 0. 
		   && hitOnInnerZPlane.Rho() >= innerRho && hitOnInnerZPlane.Rho() <= outerRho) {
	    const double exitFactor = outerRho / hitOnInnerZPlane.Rho();
	    const XYZVector hitExit = trackOrig + (kInnerZPlane * exitFactor) * trackDir;
	    hitFound = true;
	    hitPos = (hitOnInnerZPlane + hitExit) / 2.;
	    hitPathLength = (hitExit - hitOnInnerZPlane).R();
	  }
	  // Bottom right corner.
	  else if (kOuterZPlane > 0. 
		   && hitOnOuterZPlane.Rho() >= innerRho && hitOnOuterZPlane.Rho() <= outerRho) {
	    const double exitFactor = innerRho / hitOnOuterZPlane.Rho();
	    const XYZVector hitExit =  trackOrig + (kOuterZPlane * exitFactor) * trackDir;
	    hitFound = true;
	    hitPos = (hitExit + hitOnOuterZPlane) / 2.;
	    hitPathLength = (hitOnOuterZPlane - hitExit).R();
	  }
	  // Track enters and exists the volume between extrema Z planes.
	  else if ( (hitOnInnerRhoCylinder.Z() >= innerZ && hitOnInnerRhoCylinder.Z() <= outerZ)
		    && (hitOnOuterRhoCylinder.Z() >= innerZ && hitOnOuterRhoCylinder.Z() <= outerZ)
		    ) {
	    hitFound = true;
	    hitPos = (hitOnInnerRhoCylinder + hitOnOuterRhoCylinder) / 2.;
	    hitPathLength = (hitOnOuterRhoCylinder - hitOnInnerRhoCylinder).R();
	  }

	  // If hit was found, the computed path length allows to compute Material estimate.
	  if (hitFound) {
	    // Cross-check that hit is inside volume
	    if (hitPos.Rho() >= innerRho && hitPos.Rho() <= outerRho && hitPos.Z() >= innerZ && hitPos.Z() <= outerZ) {
	      // Cross-check that path length is > 0.
	      if (hitPathLength >= 0.) {
		// This normalizationFactor comes from the fact that MB is by default multiplied by length or width.
		// TO DO: Is that actually true?
		const double normalizationFactor = 1. / (isVertical() ? getZLength() : getRWidth());
		// Compute MB
		hitMaterial.radiation   = getRadiationLength() * normalizationFactor * hitPathLength;
		hitMaterial.interaction = getInteractionLength() * normalizationFactor * hitPathLength;
	      }
	      else { logERROR("InactiveElement::checkTrackHits : Computed a hitPathLength < 0."); }
	    }
	    else { logERROR("InactiveElement::checkTrackHits : Created a hit which is not inside Inactive Volume."); }
	  }

	}
	else { 
	  logERROR("InactiveElement::checkTrackHits : trackOrig.Rho() = " 
		   + any2str(trackOrig.Rho()) 
		   + " trackDir.Rho() = " 
		   + any2str(trackDir.Rho())
		   + ". IP not on (Z) axis or track parallel to (Z) axis not supported."
		   ); 
	}   
      }
      else { logERROR("InactiveElement::checkTrackHits : Try to compute inactive MB on eta < 0., which is not supported"); }

    }
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

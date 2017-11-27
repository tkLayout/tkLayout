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


  const bool InactiveElement::checkTrackHits(const XYZVector& trackOrig, const double& trackEta) const {

    double etaMin, etaMax;

    // rings
    if (isVertical()) {
      // outer radius point
      const double outerRadius = getOuterRadius();
      const double averageZ = getZOffset() + getZLength() / 2.;
      const XYZVector outerRadiusPoint(outerRadius, 0., averageZ);
      const XYZVector& origToOuterRadiusPoint = outerRadiusPoint - trackOrig;
      etaMin = origToOuterRadiusPoint.Eta();

      // inner radius point
      const double innerRadius = getInnerRadius();
      const XYZVector innerRadiusPoint(innerRadius, 0., averageZ);
      const XYZVector& origToInnerRadiusPoint = innerRadiusPoint - trackOrig;
      etaMax = origToInnerRadiusPoint.Eta();
    }

    // tubs
    else {
      // inner Z point
      const double innerZ = getZOffset();
      const double averageRadius = getInnerRadius() + getRWidth() / 2.;
      const XYZVector innerZPoint(averageRadius, 0., innerZ);
      const XYZVector& origToInnerZPoint = innerZPoint - trackOrig;
      etaMin = origToInnerZPoint.Eta();

      // outer Z point
      const double outerZ = getZMax();   
      const XYZVector outerZPoint(averageRadius, 0., outerZ);
      const XYZVector& origToOuterZPoint = outerZPoint - trackOrig;
      etaMax = origToOuterZPoint.Eta();
    }

    return (trackEta > etaMin && trackEta < etaMax);
  }

  

  //
  // Check if track hit the inactive element -> if yes, return true with passed material & hit position vector
  // Assume (0.,0., Z) origin and straigh tracks
  //
  const bool InactiveElement::checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir, XYZVector& hitPos, Material& hitMaterial) {
    // check trackOrigRho = 0
    // check eta > 0
    // trackDirRho != 0.0001 et > 0

    // Initialize: hit was found, material, hitPos & relative hit path length wrt perpendicular passage
    hitPos.SetX(0.);
    hitPos.SetY(0.);
    hitPos.SetZ(0.);
    hitMaterial.radiation   = 0.;
    hitMaterial.interaction = 0.;

    const double trackOrigZ = trackOrig.Z();
    const double trackDirZ = trackDir.Z();
    const double trackDirRho = trackDir.Rho();

    const double innerZ = getZOffset();
    const double outerZ = getZMax();
    const double innerRho = getInnerRadius();
    const double outerRho = getOuterRadius();


    // Hit on inner Z plane
    double trackDirFactorA;
    if (fabs(trackDirZ) != 0.0001) trackDirFactorA = (innerZ - trackOrigZ) / trackDirZ;
    else trackDirFactorA = innerRho / trackDirRho;
    const XYZVector hitOnInnerZPlane = trackOrig + trackDirFactorA * trackDir;

    // Hit on outer Z plane
    double trackDirFactorB;
    if (fabs(trackDirZ) != 0.0001) trackDirFactorB = (outerZ - trackOrigZ) / trackDirZ;
    else trackDirFactorB = outerRho / trackDirRho;
    const XYZVector hitOnOuterZPlane = trackOrig + trackDirFactorB * trackDir;

    // Hit on inner radius cylinder
    const double trackDirFactorC = innerRho / trackDirRho;
    const XYZVector hitOnInnerRhoCylinder = trackOrig + trackDirFactorC * trackDir;

    // Hit on outer radius cylinder
    const double trackDirFactorD = outerRho / trackDirRho;
    const XYZVector hitOnInnerRhoCylinder = trackOrig + trackDirFactorD * trackDir;


    bool hitFound = false;
    double hitPathLength = 0.;

    if ( (hitOnInnerZPlane.Rho() >= innerRadius && hitOnInnerZPlane.Rho() <= outerRadius)
	 && (hitOnOuterZPlane.Rho() >= innerRadius && hitOnOuterZPlane.Rho() <= outerRadius)
	 ) {
      hitFound = true;
      hitPos = (hitOnInnerZPlane + hitOnOuterZPlane) / 2.;
      hitPathLength = hitOnOuterZPlane.R() - hitOnInnerZPlane.R();
    }
    else if (hitOnInnerZPlane.Rho() >= innerRadius && hitOnInnerZPlane.Rho() <= outerRadius) {
      const XYZVector hitExit =  hitOnInnerZPlane * outerRho / hitOnInnerZPlane.Rho();
      hitFound = true;
      hitPos = (hitOnInnerZPlane + hitExit) / 2.;
      hitPathLength = hitExit.R() - hitOnInnerZPlane.R();
    }
    else if (hitOnOuterZPlane.Rho() >= innerRadius && hitOnOuterZPlane.Rho() <= outerRadius) {
      const XYZVector hitExit =  hitOnOuterZPlane * innerRho / hitOnOuterZPlane.Rho;
      hitFound = true;
      hitPos = (hitExit + hitOnOuterZPlane) / 2.;
      hitPathLength = hitOnOuterZPlane.R() - hitExit.R();
    }
    else if ( (hitOnInnerRhoCylinder.Z() >= innerZ && hitOnInnerRhoCylinder.Z() <= outerZ)
	      && (hitOnOuterRhoCylinder.Z() >= innerZ && hitOnOuterRhoCylinder.Z() <= outerZ)
	      ) {
      hitFound = true;
      hitPos = (hitOnInnerRhoCylinder + hitOnOuterRhoCylinder) / 2.;
      hitPathLength = hitOnOuterRhoCylinder.R() - hitOnInnerRhoCylinder.R();
    }



    if (hitFound) {
      const double normalizationFactor = 1. / (isVertical() ? getZLength() : getRWidth());
      hitMaterial.radiation   = getRadiationLength() * normalizationFactor * hitPathLength;
      hitMaterial.interaction = getInteractionLength() * normalizationFactor * hitPathLength;
    }

	// check hit is in volume
	//check pathLength is > 0

  }












  //
  // Check if track hit the inactive element -> if yes, return true with passed material & hit position vector
  //
  const bool InactiveElement::checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir, Material& hitMaterial, XYZVector& hitPos, double& factor) {
    // Initialize: hit was found, material, hitPos & relative hit path length wrt perpendicular passage
    hitMaterial.radiation   = 0.;
    hitMaterial.interaction = 0.;
    hitPos.SetX(0.);
    hitPos.SetY(0.);
    hitPos.SetZ(0.);

    bool      hitFound         = false;
    double    hitRelPathLength = 0;

    if (getRadiationLength() == 0 && getInteractionLength() == 0) {
      std::cout << "(getRadiationLength() == 0 && getInteractionLength() == 0)" << std::endl;
    }

    // Calculate hit only if inactive element non-transparent, otherwise no effect in tracking or material budget etc.
    if (getRadiationLength()!=0 || getInteractionLength()!=0) {

      // Disc
      if (isVertical()) {

	// Find number k as: vec_orig + k*vec_dir = intersection, i.e. z position of (vec_orig + k*vec_dir) equals to disc z position
	double innerZPos = getZOffset();
	double outerZPos = getZOffset() + getZLength();
	double kInner    = (innerZPos - trackOrig.z())/trackDir.z();
	double kOuter    = (outerZPos - trackOrig.z())/trackDir.z();

	// Assume origin to be "before" the disc (reasonable assumption for tkLayout)
	if (fabs(trackOrig.z())>innerZPos) {

	  std::cout << "AHH innerZPos < fabs(trackOrig.z())" << std::endl;
	  //logWARNING("InactiveElement::checkTrackHits - track origin inside tube: zInner= "+any2str(innerZPos,1)+", zOuter= "+any2str(outerZPos,1)+", check!");
	  //return false;
	}

	// Take only positive solution, so in the given direction
	if (kInner>0 && kOuter>0) {

	  XYZVector vecInner   = trackOrig + kInner*trackDir;
	  XYZVector vecOuter   = trackOrig + kOuter*trackDir;

	  // Track passes through "central" part of the disc
	  if (vecInner.rho()>=getInnerRadius() && vecInner.rho()<=getOuterRadius() && vecOuter.rho()>=getInnerRadius() && vecOuter.rho()<=getOuterRadius()) {

	    std::cout << "Disc: track IN BETWEEN" << std::endl;
	    std::cout << "kInner = " << kInner << "vecInner.rho() = " << vecInner.rho() << "vecInner.z() = " << vecInner.z() << std::endl;
	    std::cout << "kOuter = " << kOuter << "vecOuter.rho() = " << vecOuter.rho() << "vecOuter.z() = " << vecOuter.z() << std::endl;


	    hitFound = true;
	  }
	  // Track passes the inner z, but not the outer z
	  else if (vecInner.rho()>=getInnerRadius() && vecInner.rho()<=getOuterRadius()) {

	    // kOuter fixed then by getInnerRadius() or getOuterRadius()
	    /*double a  = trackDir.rho()*trackDir.rho();
	    double b  = 2*(trackOrig.x()*trackDir.x() + trackOrig.y()*trackDir.y());
	    double ci = trackOrig.rho()*trackOrig.rho() - getInnerRadius()*getInnerRadius();
	    double Di = b*b - 4*a*ci;
	    double co = trackOrig.rho()*trackOrig.rho() - getOuterRadius()*getOuterRadius();
	    double Do = b*b - 4*a*co;
	    if      (Di>0 && (-b +sqrt(Di)/2./a)>0) kOuter = (-b + sqrt(Di))/2./a; // Take only positive solution, so in the given direction
	    else if (Do>0 && (-b +sqrt(Do)/2./a)>0) kOuter = (-b + sqrt(Do))/2./a; // Take only positive solution, so in the given direction))
	    else {

	      logWARNING("InactiveElement::checkTrackHits - track passes dic inner z, but not outer z and can't find the z-pos -> seems as a bug!");
	      return false;
	      }*/

	    // SIMPLE!!
	    kOuter = kInner * (getOuterRadius() / vecInner.Rho()); 

	    hitFound = true;
	    vecOuter = trackOrig + kOuter*trackDir;

	    std::cout << "Disc: track Above cut" << std::endl;
	    std::cout << "kInner = " << kInner << "vecInner.rho() = " << vecInner.rho() << "vecInner.z() = " << vecInner.z() << std::endl;
	    std::cout << "kOuter = " << kOuter << "vecOuter.rho() = " << vecOuter.rho() << "vecOuter.z() = " << vecOuter.z() << std::endl;
	  }
	  // Track passes the outer z, but not the inner z
	  else if (vecOuter.rho()>=getInnerRadius() && vecOuter.rho()<=getOuterRadius()) {

	    // kInner fixed then by getInnerRadius() or getOuterRadius()
	    /*double a  = trackDir.rho()*trackDir.rho();
	    double b  = 2*(trackOrig.x()*trackDir.x() + trackOrig.y()*trackDir.y());
	    double ci = trackOrig.rho()*trackOrig.rho() - getInnerRadius()*getInnerRadius();
	    double Di = b*b - 4*a*ci;
	    double co = trackOrig.rho()*trackOrig.rho() - getOuterRadius()*getOuterRadius();
	    double Do = b*b - 4*a*co;
	    if      (Di>0 && (-b +sqrt(Di)/2./a)>0) kInner = (-b + sqrt(Di))/2./a; // Take only positive solution, so in the given direction
	    else if (Do>0 && (-b +sqrt(Do)/2./a)>0) kInner = (-b + sqrt(Do))/2./a; // Take only positive solution, so in the given direction))
	    else {

	      logWARNING("InactiveElement::checkTrackHits - track passes dic outer z, but not inner z and can't find the z-pos -> seems as a bug!");
	      return false;
	      }*/

	    // SIMPLE!!
	    kInner = kOuter * (getInnerRadius() / vecOuter.Rho()); 

	    hitFound = true;
	    vecInner = trackOrig + kInner*trackDir;

	    std::cout << "Disc: track below cut" << std::endl;
	    std::cout << "kInner = " << kInner << "vecInner.rho() = " << vecInner.rho() << "vecInner.z() = " << vecInner.z() << std::endl;
	    std::cout << "kOuter = " << kOuter << "vecOuter.rho() = " << vecOuter.rho() << "vecOuter.z() = " << vecOuter.z() << std::endl;
	  }

	  else {
	    std::cout << "Disc: very steep angle, nottaken into account!!!!" << std::endl;
	  }

	  if (hitFound) {

	    hitFound                = true;
	    hitRelPathLength        = fabs(kOuter-kInner)*trackDir.r()/getZLength();
	    hitMaterial.radiation   = getRadiationLength()*hitRelPathLength;
	    hitMaterial.interaction = getInteractionLength()*hitRelPathLength;
	    hitPos                  = (vecInner + vecOuter)/2.;
	  }
	}
      }
      // Tube
      else {

	// Find number k as: vec_orig + k*vec_dir = intersection, i.e. radial position of (vec_orig + k*vec_dir) equals radial position
	double innerRPos = getInnerRadius();
	double outerRPos = getOuterRadius();
	double kInner    = -1;
	double kOuter    = -1;

	// Assume origin to be inside the tube (reasonable assumption for tkLayout)
	if (trackOrig.rho()>innerRPos) {
	  logWARNING("InactiveElement::checkTrackHits - track origin outside of tube inner radius= "+any2str(innerRPos,1)+", check!");
	  return false;
	}

	// Calculate kInner & kOuter: (vec_orig_X + k*vec_dir_X)^2 + (vec_orig_Y + k*vec_dir_Y)^2 = r^2
	/*double a = trackDir.rho()*trackDir.rho();
	double b = 2*(trackOrig.x()*trackDir.x() + trackOrig.y()*trackDir.y());
	double c = trackOrig.rho()*trackOrig.rho() - innerRPos*innerRPos;
	double D = b*b - 4*a*c;
	if (D>0) kInner = (-b + sqrt(D))/2/a; // Take only positive solution, so in the given direction

	c = trackOrig.rho()*trackOrig.rho() - outerRPos*outerRPos;
	D = b*b - 4*a*c;
	if (D>0) kOuter = (-b + sqrt(D))/2/a; // Take only positive solution, so in the given direction*/

	// SIMPLE!!
	double innerRadius = getInnerRadius();
	double outerRadius = getOuterRadius();
	kInner    = innerRadius/trackDir.rho();
	kOuter    = outerRadius/trackDir.rho();


	// Due to condition on trackOrig, both solutions exist
	XYZVector vecInner   = trackOrig + kInner*trackDir;
	XYZVector vecOuter   = trackOrig + kOuter*trackDir;

	// Track passes through "central" part of the tube
	if (vecInner.z()>=getZOffset() && vecInner.z()<=(getZOffset()+getZLength()) &&
	    vecOuter.z()>=getZOffset() && vecOuter.z()<=(getZOffset()+getZLength())) {

	  hitFound = true;

	  std::cout << "Tub: track IN BETWEEN" << std::endl;
	  std::cout << "kInner = " << kInner << "vecInner.rho() = " << vecInner.rho() << "vecInner.z() = " << vecInner.z() << std::endl;
	  std::cout << "kOuter = " << kOuter << "vecOuter.rho() = " << vecOuter.rho() << "vecOuter.z() = " << vecOuter.z() << std::endl;

	}
	// Track passes the inner radius, but not the outer radius
	else if (vecInner.z()>=getZOffset() && vecInner.z()<=(getZOffset()+getZLength())) {

	  // kOuter fixed then by getZOffset() + getZLength() or getZOffset()
	  /*if     (((getZOffset() + getZLength()) - trackOrig.z())/trackDir.z()>0 ) kOuter = ((getZOffset() + getZLength()) - trackOrig.z())/trackDir.z();
	  else if((getZOffset() - trackOrig.z())/trackDir.z()>0 )               kOuter = (getZOffset() - trackOrig.z())/trackDir.z();
	  else {
	    logWARNING("InactiveElement::checkTrackHits - track passes tube inner radius, but not outer radius and can't find the radius -> seems as a bug!");
	    return false;
	    }*/

	  // SIMPLE!!
	  kOuter = kInner * (getZMax() - trackOrig.z()) / (vecInner.Z() - trackOrig.z()); 

	  hitFound = true;
	  vecOuter = trackOrig + kOuter*trackDir;

	  std::cout << "Tub: track right corner" << std::endl;
	  std::cout << "kInner = " << kInner << "vecInner.rho() = " << vecInner.rho() << "vecInner.z() = " << vecInner.z() << std::endl;
	  std::cout << "kOuter = " << kOuter << "vecOuter.rho() = " << vecOuter.rho() << "vecOuter.z() = " << vecOuter.z() << std::endl;

	}
	// Track passes the Z leftmost corner
	else if (vecOuter.z()>=getZOffset() && vecOuter.z()<=(getZOffset()+getZLength())) {

	  // kInner fixed then by getZOffset() + getZLength() or getZOffset()
	  /* if     (((getZOffset() + getZLength()) - trackOrig.z())/trackDir.z()>0 ) kInner = ((getZOffset() + getZLength()) - trackOrig.z())/trackDir.z();
	  else if((getZOffset() - trackOrig.z())/trackDir.z()>0 )               kInner = (getZOffset() - trackOrig.z())/trackDir.z();
	  else {
	    logWARNING("InactiveElement::checkTrackHits - track passes tube outer radius, but not inner radius and can't find the radius -> seems as a bug!");
	    return false;
	    }*/

	  // SIMPLE!!
	  kInner = kOuter * (getZOffset() - trackOrig.z()) / (vecOuter.Z() - trackOrig.z()); 

	  hitFound = true;
	  vecInner = trackOrig + kInner*trackDir;

	  std::cout << "Tub: track left corner" << std::endl;
	  std::cout << "kInner = " << kInner << "vecInner.rho() = " << vecInner.rho() << "vecInner.z() = " << vecInner.z() << std::endl;
	  std::cout << "kOuter = " << kOuter << "vecOuter.rho() = " << vecOuter.rho() << "vecOuter.z() = " << vecOuter.z() << std::endl;
	}

	else {
	    std::cout << "Tube: very steep angle, nottaken into account!!!!" << std::endl;
	  }

	if (hitFound) {

	  hitFound                = true;
	  hitRelPathLength        = fabs(kOuter-kInner)*trackDir.r()/getRWidth();
	  hitMaterial.radiation   = getRadiationLength()*hitRelPathLength;
	  hitMaterial.interaction = getInteractionLength()*hitRelPathLength;
	  hitPos                  = (vecInner + vecOuter)/2.;
	}

      } // Tube
    }

    factor = hitRelPathLength;
    std::cout << "hitRelPathLength = " << hitRelPathLength << std::endl;

    return hitFound;
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

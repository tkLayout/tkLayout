/**
 * @file InactiveElement.cc
 * @brief This is the base class implementation of a single inactive element (tubes & discs).
 */
#include "InactiveElement.h"

#include <cmath>
#include <iostream>

#include <Math/GenVector/Transform3D.h>
#include "MessageLogger.h"
#include "Units.h"

/**
 * The constructor sets some defaults: no neighbours, intermediate element.
 */
InactiveElement::InactiveElement(double zOffset, double zLength, double rInner, double rWidth, bool isVertical) :
 m_zOffset(zOffset),
 m_zLength(zLength),
 m_rInner(rInner),
 m_rWidth(rWidth),
 m_rOuter(rInner+rWidth),
 m_isVertical(isVertical)
{

  m_feederType     = no_in;
  m_feederIndex    = -1;
  m_neighbourType  = no_in;
  m_neighbourIndex = -1;

  if (m_isVertical) m_surface = M_PI * (m_rOuter*m_rOuter - m_rInner*m_rInner);
  else              m_surface = 2*M_PI * (m_rInner+m_rOuter)/2. * m_zLength;

  m_volume = M_PI * m_rWidth * (m_rWidth + 2*m_rInner) * m_zLength;
}

/**
 * The copy constructor
 */
InactiveElement::InactiveElement(InactiveElement& element)
{
  m_zOffset   = element.m_zOffset;
  m_zLength   = element.m_zLength;
  m_rInner    = element.m_rInner;
  m_rWidth    = element.m_rWidth;
  m_rOuter    = element.m_rOuter;
  m_isVertical= element.m_isVertical;
  m_surface   = element.m_surface;
  m_volume    = element.m_volume;

  m_feederType     = no_in;
  m_feederIndex    = -1;
  m_neighbourType  = no_in;
  m_neighbourIndex = -1;
}
 

/**
 * Checks whether track hits the inactive element -> if yes, return true with passed hit position vector & material. 
 * Only implemeted for eta >= 0. (should work out of the box for eta < 0. but need to be tested in that case).
 * Assumes tracks are straight lines!
 * IP should be on (Z) axis, and track not parallel to (Z) axis.
 */

bool InactiveElement::checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir, Material& hitMaterial, XYZVector& hitPos) const {

  // Initialize
  bool hitFound = false;
  hitPos.SetX(0.);
  hitPos.SetY(0.);
  hitPos.SetZ(0.);
  hitMaterial.radiation   = 0.;
  hitMaterial.interaction = 0.;

  const double geom_zero = 1E-6; // devLite: add in global constants

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
	const double outerZ = m_zOffset + m_zLength; // devLite: add method
	const double innerRho = getInnerRadius();
	const double outerRho = m_rInner + m_rWidth; // devLite: add method

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
  return hitFound;
}


/**
 * Calculate and return the Eta range of the element
 * @return The pair <i>(Eta_min, Eta_max)</i>
 */
std::pair<double, double> InactiveElement::getEtaMinMax() const {

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
  std::cout << "zOffset = " << m_zOffset << std::endl;
  std::cout << "zLength = " << m_zLength << std::endl;
  std::cout << "rInner  = " << m_rInner  << std::endl;
  std::cout << "rWidth  = " << m_rWidth  << std::endl;
  if (isVertical()) std::cout << "Volume is vertical." << std::endl;
  else              std::cout << "Volume is horizontal" << std::endl;
}

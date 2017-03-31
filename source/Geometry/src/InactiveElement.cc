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
    
//
// Check if track hit the inactive element -> if yes, return true with passed material & hit position vector
//
bool InactiveElement::checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir, Material& hitMaterial, XYZVector& hitPos) const
{
  // Initialize: hit was found, material, hitPos & relative hit path length wrt perpendicular passage
  hitMaterial.radiation   = 0.;
  hitMaterial.interaction = 0.;
  hitPos.SetX(0.);
  hitPos.SetY(0.);
  hitPos.SetZ(0.);

  bool      hitFound         = false;
  double    hitRelPathLength = 0;

  // Calculate hit only if inactive element non-transparent, otherwise no effect in tracking or material budget etc.
  if (getRadiationLength()!=0 || getInteractionLength()!=0) {

    // Disc
    if (m_isVertical) {

      // Find number k as: vec_orig + k*vec_dir = intersection, i.e. z position of (vec_orig + k*vec_dir) equals to disc z position
      double innerZPos = m_zOffset;
      double outerZPos = m_zOffset+m_zLength;
      double kInner    = (innerZPos - trackOrig.z())/trackDir.z();
      double kOuter    = (outerZPos - trackOrig.z())/trackDir.z();

      // Assume origin to be "before" the disc (reasonable assumption for tkLayout)
      if (fabs(trackOrig.z())>innerZPos) {

        logWARNING("InactiveElement::checkTrackHits - track origin inside tube: zInner= "+any2str(innerZPos/Units::mm,1)+", zOuter= "+any2str(outerZPos/Units::mm,1)+", check!");
        return false;
      }

      // Take only positive solution, so in the given direction
      if (kInner>0 && kOuter>0) {

        XYZVector vecInner   = trackOrig + kInner*trackDir;
        XYZVector vecOuter   = trackOrig + kOuter*trackDir;

        // Track passes through "central" part of the disc
        if (vecInner.rho()>=m_rInner && vecInner.rho()<=m_rOuter && vecOuter.rho()>=m_rInner && vecOuter.rho()<=m_rOuter) {

          hitFound = true;
        }
        // Track passes the inner z, but not the outer z
        else if (vecInner.rho()>=m_rInner && vecInner.rho()<=m_rOuter) {

          // kOuter fixed then by m_rInner or m_rOuter
          double a  = trackDir.rho()*trackDir.rho();
          double b  = 2*(trackOrig.x()*trackDir.x() + trackOrig.y()*trackDir.y());
          double ci = trackOrig.rho()*trackOrig.rho() - m_rInner*m_rInner;
          double Di = b*b - 4*a*ci;
          double co = trackOrig.rho()*trackOrig.rho() - m_rOuter*m_rOuter;
          double Do = b*b - 4*a*co;
          if      (Di>0 && (-b +sqrt(Di)/2./a)>0) kOuter = (-b + sqrt(Di))/2./a; // Take only positive solution, so in the given direction
          else if (Do>0 && (-b +sqrt(Do)/2./a)>0) kOuter = (-b + sqrt(Do))/2./a; // Take only positive solution, so in the given direction))
          else {

            logWARNING("InactiveElement::checkTrackHits - track passes dic inner z, but not outer z and can't find the z-pos -> seems as a bug!");
            return false;
          }

          hitFound = true;
          vecOuter = trackOrig + kOuter*trackDir;
        }
        // Track passes the outer z, but not the inner z
        else if (vecOuter.rho()>=m_rInner && vecOuter.rho()<=m_rOuter) {

          // kInner fixed then by m_rInner or m_rOuter
          double a  = trackDir.rho()*trackDir.rho();
          double b  = 2*(trackOrig.x()*trackDir.x() + trackOrig.y()*trackDir.y());
          double ci = trackOrig.rho()*trackOrig.rho() - m_rInner*m_rInner;
          double Di = b*b - 4*a*ci;
          double co = trackOrig.rho()*trackOrig.rho() - m_rOuter*m_rOuter;
          double Do = b*b - 4*a*co;
          if      (Di>0 && (-b +sqrt(Di)/2./a)>0) kInner = (-b + sqrt(Di))/2./a; // Take only positive solution, so in the given direction
          else if (Do>0 && (-b +sqrt(Do)/2./a)>0) kInner = (-b + sqrt(Do))/2./a; // Take only positive solution, so in the given direction))
          else {

            logWARNING("InactiveElement::checkTrackHits - track passes dic outer z, but not inner z and can't find the z-pos -> seems as a bug!");
            return false;
          }

          hitFound = true;
          vecInner = trackOrig + kInner*trackDir;
        }

        if (hitFound) {

          hitFound                = true;
          hitRelPathLength        = fabs(kOuter-kInner)*trackDir.r()/m_zLength;
          hitMaterial.radiation   = getRadiationLength()*hitRelPathLength;
          hitMaterial.interaction = getInteractionLength()*hitRelPathLength;
          hitPos                  = (vecInner + vecOuter)/2.;
        }
      }
    }
    // Tube
    else {

      // Find number k as: vec_orig + k*vec_dir = intersection, i.e. radial position of (vec_orig + k*vec_dir) equals radial position
      double innerRPos = m_rInner;
      double outerRPos = m_rOuter;
      double kInner    = -1;
      double kOuter    = -1;

      // Assume origin to be inside the tube (reasonable assumption for tkLayout)
      if (trackOrig.rho()>innerRPos) {
        logWARNING("InactiveElement::checkTrackHits - track origin outside of tube inner radius= "+any2str(innerRPos/Units::mm,1)+", check!");
        return false;
      }

      // Calculate kInner & kOuter: (vec_orig_X + k*vec_dir_X)^2 + (vec_orig_Y + k*vec_dir_Y)^2 = r^2
      double a = trackDir.rho()*trackDir.rho();
      double b = 2*(trackOrig.x()*trackDir.x() + trackOrig.y()*trackDir.y());
      double c = trackOrig.rho()*trackOrig.rho() - innerRPos*innerRPos;
      double D = b*b - 4*a*c;
      if (D>0) kInner = (-b + sqrt(D))/2/a; // Take only positive solution, so in the given direction

      c = trackOrig.rho()*trackOrig.rho() - outerRPos*outerRPos;
      D = b*b - 4*a*c;
      if (D>0) kOuter = (-b + sqrt(D))/2/a; // Take only positive solution, so in the given direction

      // Due to condition on trackOrig, both solutions exist
      XYZVector vecInner   = trackOrig + kInner*trackDir;
      XYZVector vecOuter   = trackOrig + kOuter*trackDir;

      // Track passes through "central" part of the tube
      if (vecInner.z()>=m_zOffset && vecInner.z()<=(m_zOffset+m_zLength) &&
          vecOuter.z()>=m_zOffset && vecOuter.z()<=(m_zOffset+m_zLength)) {

        hitFound = true;
      }
      // Track passes the inner radius, but not the outer radius
      else if (vecInner.z()>=m_zOffset && vecInner.z()<=(m_zOffset+m_zLength)) {

        // kOuter fixed then by m_zOffset + m_zLength or m_zOffset
        if     (((m_zOffset + m_zLength) - trackOrig.z())/trackDir.z()>0 ) kOuter = ((m_zOffset + m_zLength) - trackOrig.z())/trackDir.z();
        else if((m_zOffset - trackOrig.z())/trackDir.z()>0 )               kOuter = (m_zOffset - trackOrig.z())/trackDir.z();
        else {
          logWARNING("InactiveElement::checkTrackHits - track passes tube inner radius, but not outer radius and can't find the radius -> seems as a bug!");
          return false;
        }

        hitFound = true;
        vecOuter = trackOrig + kOuter*trackDir;
      }
      // Track passes the Z leftmost corner
      else if (vecOuter.z()>=m_zOffset && vecOuter.z()<=(m_zOffset+m_zLength)) {

        // kInner fixed then by m_zOffset + m_zLength or m_zOffset
        if     (((m_zOffset + m_zLength) - trackOrig.z())/trackDir.z()>0 ) kInner = ((m_zOffset + m_zLength) - trackOrig.z())/trackDir.z();
        else if((m_zOffset - trackOrig.z())/trackDir.z()>0 )               kInner = (m_zOffset - trackOrig.z())/trackDir.z();
        else {
          logWARNING("InactiveElement::checkTrackHits - track passes tube outer radius, but not inner radius and can't find the radius -> seems as a bug!");
          return false;
        }

        hitFound = true;
        vecInner = trackOrig + kInner*trackDir;
      }

      if (hitFound) {

        hitFound                = true;
        hitRelPathLength        = fabs(kOuter-kInner)*trackDir.r()/m_rWidth;
        hitMaterial.radiation   = getRadiationLength()*hitRelPathLength;
        hitMaterial.interaction = getInteractionLength()*hitRelPathLength;
        hitPos                  = (vecInner + vecOuter)/2.;
      }

    } // Tube
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

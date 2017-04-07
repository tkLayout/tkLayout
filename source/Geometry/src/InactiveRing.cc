/**
 * @file InactiveRing.cc
 * @brief This is the derived class implementation of a single ring-shaped inactive element
 */

#include "InactiveRing.h"

#include "MessageLogger.h"

/**
 * This function prints a representation of the element to <i>cout</i>
 */
void InactiveRing::print() {

  InactiveElement::print();
  if (sanityCheck()) logINFO("InactiveRing volume is sane");
  else               logWARNING("InactiveRing volume is not sane!");
}
    
/**
 * The sanity check virtual function tests if the object has the geometric properties of a ring.
 * @return True if the length is less than or equal to the width, false otherwise
 */
bool InactiveRing::sanityCheck() {
  return (m_zLength <= m_rWidth) && m_isVertical;
}

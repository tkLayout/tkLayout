/**
 * @file InactiveTube.cc
 * @brief This is the derived class implementation of a single tube-shaped inactive element
 */

#include "InactiveTube.h"
#include "MessageLogger.h"

/**
 * This function prints a representation of the element to <i>cout</i>
 */
void InactiveTube::print() {
  InactiveElement::print();
  if (sanityCheck()) logINFO("InactiveTube volume is sane");
  else               logWARNING("InactiveTube volume is not sane!");
}
    
/**
 * The sanity check virtual function tests if the object has the geometric properties of a tube.
 * @return True if the length is greater than or equal to the width, false otherwise
 */
bool InactiveTube::sanityCheck() {
  return (m_zLength >= m_rWidth) && !m_isVertical;
}

//
// File:   InactiveRing.h
// Author: ndemaio
//
// Created on October 13, 2008, 2:53 PM
//
#ifndef INCLUDE_INACTIVERING_H_
#define INCLUDE_INACTIVERING_H_

#include "InactiveElement.h"

/**
 * @class InactiveRing
 * @brief The only thing that this class adds to its parent is a check that it is a ring rather than a tube.
 */
class InactiveRing : public InactiveElement {

 public:
  InactiveRing(double zOffset, double zLength, double rInner, double rWidth) : InactiveElement(zOffset, zLength, rInner, rWidth, true) {};
  InactiveRing(InactiveElement& element) : InactiveElement(element) {};
  virtual ~InactiveRing() {};

  virtual void print();
  virtual bool sanityCheck();

}; // Class

#endif	/* INCLUDE_INACTIVERING_H_ */


//
// File:   InactiveTube.h
// Author: ndemaio
//
// Created on October 13, 2008, 2:53 PM
//
#ifndef INCLUDE_INACTIVETUBE_H_
#define INCLUDE_INACTIVETUBE_H_

#include "InactiveElement.h"

/**
 * @class InactiveTube
 * @brief The only thing that this class adds to its parent is a check that it is a tube rather than a ring.
 */
class InactiveTube : public InactiveElement {

 public:
  InactiveTube(double zOffset, double zLength, double rInner, double rWidth) : InactiveElement(zOffset, zLength, rInner, rWidth, false) {};
  InactiveTube(InactiveElement& element) : InactiveElement(element) {};
  virtual ~InactiveTube() {};

  virtual void print();
  virtual bool sanityCheck();

}; // Class

#endif	/* INCLUDE_INACTIVETUBE_H_ */


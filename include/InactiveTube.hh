//
// File:   InactiveTube.h
// Author: ndemaio
//
// Created on October 13, 2008, 2:53 PM
//

/**
 * @file InactiveTube.h
 * @brief This is the derived class header file for a single tube-shaped inactive element
 */

#ifndef _INACTIVETUBE_H
#define	_INACTIVETUBE_H

#include <InactiveElement.hh>
namespace insur {
  /**
   * @class InactiveTube
   * @brief The only thing that this class adds to its parent is a check that it is a tube rather than a ring.
   */
  class InactiveTube : public InactiveElement {
  public:
    InactiveTube();
    InactiveTube(InactiveElement& previous);
    virtual ~InactiveTube();
    virtual void print();
    virtual bool sanityCheck();
  };
}
#endif	/* _INACTIVETUBE_H */


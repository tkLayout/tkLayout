//
// File:   InactiveRing.h
// Author: ndemaio
//
// Created on October 13, 2008, 2:53 PM
//

/**
 * @file InactiveRing.h
 * @brief This is the derived class header file for a single ring-shaped inactive element
 */

#ifndef _INACTIVERING_H
#define	_INACTIVERING_H

#include <InactiveElement.hh>
namespace insur {
  /**
   * @class InactiveRing
   * @brief The only thing that this class adds to its parent is a check that it is a ring rather than a tube.
   */
  class InactiveRing : public InactiveElement {
  public:
    InactiveRing();
    InactiveRing(InactiveElement& previous);
    virtual ~InactiveRing();
    virtual void print();
    virtual bool sanityCheck();
  };
}
#endif	/* _INACTIVERING_H */


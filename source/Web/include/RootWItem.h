/*
 * RootWItem.h
 *
 *  Created on: 21. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */

#ifndef INCLUDE_ROOTWITEM_H_
#define INCLUDE_ROOTWITEM_H_

#include <ostream>

/*
 * @class RootWItem
 * Pure virtual class, from which all html items (images, texts, ...) must be derived.
 */
class RootWItem {

 public:

  //! Default constructor
  RootWItem();

  //! Default destructor
  virtual ~RootWItem();

  //! Pure virtual dump method to be implemented by derived classes
  virtual std::ostream& dump(std::ostream& output) = 0;

}; // Class

#endif /* INCLUDE_ROOTWITEM_H_ */

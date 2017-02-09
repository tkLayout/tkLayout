/*
 * RootWText.h
 *
 *  Created on: 21. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */
#ifndef SOURCE_WEB_INCLUDE_ROOTWTEXT_H_
#define SOURCE_WEB_INCLUDE_ROOTWTEXT_H_

#include <string>
#include <ostream>
#include <sstream>

#include "RootWItem.h"

/*
 * @class RootWText
 * A simple text - no html formatting
 */
class RootWText: public RootWItem {

 public:

  RootWText() : RootWItem()                    {m_text.clear();};
  RootWText(std::string newText) : RootWItem() {m_text.str(newText); };
  ~RootWText() {};

  //! Add text to the container
  void addText(std::string newText);

  //! Dump the text
  virtual std::ostream& dump(std::ostream& output);

 protected:

  std::stringstream m_text;

}; // Class

#endif /* SOURCE_WEB_INCLUDE_ROOTWTEXT_H_ */

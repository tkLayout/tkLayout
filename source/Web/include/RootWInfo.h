/*
 * RootWInfo.h
 *
 *  Created on: 21. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */
#ifndef INCLUDE_ROOTWINFO_H_
#define INCLUDE_ROOTWINFO_H_

#include <string>
#include <ostream>

#include "RootWItem.h"

/*
 * @class RootWInfo
 * A simple html info represented as a line
 */
class RootWInfo: public RootWItem {

 public:

  RootWInfo() : RootWItem()                                           {m_description="missing_description"; m_value="missing_value";};
  RootWInfo(std::string description) : RootWItem()                    {m_description=description; m_value="missing_value";};
  RootWInfo(std::string description, std::string value) : RootWItem() {m_description=description; m_value=value;};
  ~RootWInfo() {};

  std::string setDescription(std::string newText) { m_description = newText; return newText; };
  std::string setValue(std::string newText);
  std::string setValue(int number);
  std::string setValue(double number, int precision);
  std::string addValueText(std::string newText) {m_value+=newText; return m_value;};

  //! Dump method - printing info in the html format
  virtual std::ostream& dump(std::ostream& output);

 protected:
  std::string m_description;
  std::string m_value;

}; // Class

#endif /* INCLUDE_ROOTWINFO_H_ */

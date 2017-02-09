/*
 * RootWInfo.cc
 *
 *  Created on: 21. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.cc to standalone file
 */
#include "RootWInfo.h"

#include <iomanip>
#include <sstream>

std::string RootWInfo::setValue(std::string newText) {
  m_value = newText;
  return newText;
}

std::string RootWInfo::setValue(int number) {
  std::stringstream myNum;
  myNum.clear();
  myNum << std::dec << number;
  return(setValue(myNum.str()));
}

std::string RootWInfo::setValue(double number, int precision) {
  std::stringstream myNum_;
  myNum_.clear();
  myNum_ << std::dec << std::fixed << std::setprecision(precision) << number;
  return(setValue(myNum_.str()));
}

//
// Dump method - printing info in the html format
//
std::ostream& RootWInfo::dump(std::ostream& output) {

  output << "<b>" << m_description << ":</b> "
         << m_value << "</a></tt><br/>";
  return output;
}




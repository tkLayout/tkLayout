/*
 * RootWText.cc
 *
 *  Created on: 21. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.cc to standalone file
 */
#include "RootWText.h"

//
// Add text to the container
//
void RootWText::addText(std::string newText) {
  m_text << newText;
};

//
// Dump the text
//
std::ostream& RootWText::dump(std::ostream& output) {

  output << m_text.str();
  return output;
};



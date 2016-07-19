/*
 * CsvTextBuilder.cc
 *
 *  Created on: 16. 10. 2015
 *      Author: Drasal (CERN)
 */
#include <CsvTextBuilder.h>
#include <sstream>
#include <iostream>

/*
 * Add an element (element+;) to a csv formatted text (map of strings),
 * id identifies a given section in the text (label etc.)
 */
void CsvTextBuilder::addCsvElement(std::string id, std::string element) {

  auto itMap = m_text.find(id);
  if (itMap != m_text.end()) m_text[id] += element+std::string(m_csvSeparator);
  else                       m_text[id]  = element+std::string(m_csvSeparator);
}

void CsvTextBuilder::addCsvElement(std::string id, double element) {

  std::ostringstream myElement;
  myElement.str("");
  myElement << element;
  addCsvElement(id, myElement.str());
}

/*
 * Remove last csv separator character and add end-of-line character to a csv formatted text
 */
void CsvTextBuilder::addCsvEOL(std::string id) {

  if (m_text[id].size()>0) {

    // Remove last semi-colon if exists
    std::string::size_type iPos = m_text[id].find_last_of(m_csvSeparator);
    if (iPos!=std::string::npos) m_text[id].erase(iPos, m_text[id].length());
  }
  m_text[id].push_back('\n');
}

/*
 * Read text block content
 */
std::string CsvTextBuilder::getCsvText(std::string id) const {

  if (existCsvText(id)) return m_text.at(id);
  else                  return "";
}
/*
 * Clear text block content
 */
void CsvTextBuilder::clearCsVText(std::string id) {

  m_text.erase(id);
}

/*
 *  Test if text block exists
 */
bool CsvTextBuilder::existCsvText(std::string id) const {

  if (m_text.find(id)!=m_text.end()) return true;
  else                               return false;
}

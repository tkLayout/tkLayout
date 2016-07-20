/*
 * CsvTextBuilder.h
 *
 *  Created on: 16. 10. 2015
 *      Author: Drasal (CERN)
 */
#ifndef INCLUDE_CSVTEXTBUILDER_H_
#define INCLUDE_CSVTEXTBUILDER_H_

#include <map>
#include <string>

/*
 * @class CsvTextBuilder
 * @brief Build a named list of csv formatted blocks of text.
 * @details Build a named list of csv formatted blocks of text. Each block is identified with a unique id. For example,
 * the first block may be a "label", the second an "information about inner tracker", ... Use EOL method to add end-of-lign
 * character and start writing text to a new line.
 */
class CsvTextBuilder {

 public:
  //! Constructor - default
  CsvTextBuilder() : m_csvSeparator(";") {}
  //! Constructor - set csv separator
  CsvTextBuilder(std::string csvSeparator) : m_csvSeparator(csvSeparator) {}
  //! Add an element (element+;) to a csv formatted text (map of strings), id identifies a given section in the text (label etc.)
  void addCsvElement(std::string id, std::string element);
  void addCsvElement(std::string id,  double element);
  //! Remove last ';' character and add end-of-line character to a csv formatted text
  void addCsvEOL(std::string id);
  //! Read text block content
  std::string getCsvText(std::string id) const;
  //! Clear text block content
  void clearCsVText(std::string id);
  //! Test if text block exists
  bool existCsvText(std::string id) const;

  // Iterators over ids & individual text blocks
  std::map<std::string, std::string>::iterator       getCsvTextBegin() { return m_text.begin();}
  std::map<std::string, std::string>::const_iterator getCsvTextBegin() const { return m_text.begin();}
  std::map<std::string, std::string>::iterator       getCsvTextEnd()   { return m_text.end();}
  std::map<std::string, std::string>::const_iterator getCsvTextEnd() const { return m_text.end();}

 private:
  std::string                        m_csvSeparator; //!< csv separator character
  std::map<std::string, std::string> m_text;         //!< by ID one identifies the text block


};
#endif /* INCLUDE_CSVTEXTBUILDER_H_ */

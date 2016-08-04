/*
 * RootWTable.h
 *
 *  Created on: 21. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */

#ifndef INCLUDE_ROOTWTABLE_H_
#define INCLUDE_ROOTWTABLE_H_

#include <map>
#include <string>
#include <ostream>

#include "RootWItem.h"

/*
 * @class RootWTable
 * An html table created from string content added directly to a given table cell or automatically
 * added to a new cell based on internal row/column counting.
 */
class RootWTable : public RootWItem {

 public:

  //! Default constructor
  RootWTable();
  //!Default destructor
  ~RootWTable() {};

  //! Set content to given row & column
  void setContent(int row, int column, std::string content);
  void setContent(int row, int column, int number);
  void setContent(int row, int column, double number, int precision);
  //void setContent(const rootWTableContent& newContent) { m_tableContent = newContent; };

  //! Set ROOT defined color to color the text in a given cell
  void setColor(int row, int column, int newColor);

  //! Set content to current row & new column
  std::pair<int, int> addContent(std::string content);
  std::pair<int, int> addContent(int number);
  std::pair<int, int> addContent(double number, int precision);

  //! Set new row to be used in addContent() method
  std::pair<int, int> newLine();

  //! Dump method - printing info in the html format
  virtual std::ostream& dump(std::ostream& output);

 private:

  std::map<std::pair<int,int>, std::string> m_tableContent;
  std::map<std::pair<int,int>, int>         m_tableContentColor;

  int m_serialRow, m_serialCol;

}; // Class

#endif /* INCLUDE_ROOTWTABLE_H_ */

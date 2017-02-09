/*
 * RootWTable.cc
 *
 *  Created on: 21. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.cc to standalone file
 */
#include "RootWTable.h"

#include <iomanip>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TColor.h"
#include "MessageLogger.h"

using std::pair;

//
// Default constructor
//
RootWTable::RootWTable() {
  m_serialRow = 0;
  m_serialCol = 0;
}

//
// Dump method - printing info in the html format
//
std::ostream& RootWTable::dump(std::ostream& output) {
  RootWTable& myRootWTable = (*this);

  int minRow = m_tableContent.begin()->first.first; int maxRow=0;
  int minCol = m_tableContent.begin()->first.second; int maxCol=0;
  bool firstNotFound = true;

  pair<int, int> myIndex;
  auto& myTableContent = myRootWTable.m_tableContent;

  //std::cerr << "Table: "; // debug
  for (auto tableContentIt = m_tableContent.begin(); tableContentIt != m_tableContent.end(); ++tableContentIt) {

    myIndex = (*tableContentIt).first;
    //std::cerr << "(" << myIndex.first << ", " << myIndex.second << ")" << std::endl; // debug
    if (firstNotFound) {
      firstNotFound = false;
      minRow = myIndex.first;
      maxRow = myIndex.first;
      minCol = myIndex.second;
      maxCol = myIndex.second;
    } else {
      minRow = ( myIndex.first < minRow ) ? myIndex.first : minRow;
      maxRow = ( myIndex.first > maxRow ) ? myIndex.first : maxRow;
      minCol = ( myIndex.second < minCol ) ? myIndex.second : minCol;
      maxCol = ( myIndex.second > maxCol ) ? myIndex.second : maxCol;
    }
  }
  //std::cerr << "Size:" << m_tableContent.size() << "; Rows: " << minRow << "-" << maxRow <<"; Cols: " << minCol << "-" << maxCol << endl; // debug
  if (firstNotFound) return output;
  std::string myColorCode;
  std::string myCellCode;
  int myColorIndex;

  output << "<table>";
  for (int iRow = minRow; iRow<=maxRow; ++iRow) {
    output << "<tr>";
    for (int iCol = minCol; iCol<=maxCol; ++iCol) {
      if ((iRow==minRow)&&(iRow==0)) myCellCode = "th";
      else myCellCode = "td";
      output << "<" << myCellCode;
      myColorIndex = m_tableContentColor[std::make_pair(iRow, iCol)];
      if (myColorIndex!=0) {

        TColor* colour = gROOT->GetColor(myColorIndex);
        if (colour!=nullptr) output << " style=\"color:" << colour->AsHexString() << ";\" " << std::endl;
        else {

          std::ostringstream message;
          message << "RootWTable::dump - Problem with color index: " << myColorIndex << std::endl;
          logERROR(message);
        }
      }
      output << ">";
      output << myTableContent[std::make_pair(iRow, iCol)];
      output << "</" << myCellCode << ">" << " ";
    }
    output << "</tr>";
  }
  output << "</table>" << std::endl;

  return output;
}

//
// Set ROOT defined color to color the text in a given cell
//
void RootWTable::setColor(int row, int column, int newColor) {

  m_tableContentColor[std::make_pair(row, column)] = newColor;
}

//
// Set content to given row & column
//
void RootWTable::setContent(int row, int column, std::string content) {
  // std::cerr << "setContent("<<row<<", "<<column<<", "<<content<<")"<<endl; // debug
  m_tableContent[std::make_pair(row, column)] = content;
}
void RootWTable::setContent(int row, int column, int number) {
  // std::cerr << "setContent("<<row<<", "<<column<<", "<<number<<")"<<endl; // debug
  std::stringstream myNum;
  myNum.clear();
  myNum << std::dec << number;
  m_tableContent[std::make_pair(row, column)] = myNum.str();
}
void RootWTable::setContent(int row, int column, double number, int precision) {
  // std::cerr << "setContent("<<row<<", "<<column<<", "<<number<<")"<<endl; // debug
  std::stringstream myNum;
  myNum.clear();
  myNum << std::dec << std::fixed << std::setprecision(precision) << number;
  m_tableContent[std::make_pair(row, column)] = myNum.str();
}

//
// Set content to current row & new column
//
pair<int, int> RootWTable::addContent(std::string myContent) {
  setContent(m_serialRow, m_serialCol++, myContent);
  return std::make_pair(m_serialRow, m_serialCol);
}
pair<int, int> RootWTable::addContent(int number) {
  std::stringstream myNum;
  myNum.clear();
  myNum << std::dec << number;
  return(addContent(myNum.str()));
}
pair<int, int> RootWTable::addContent(double number, int precision) {
  std::stringstream myNum;
  myNum.clear();
  myNum << std::dec << std::fixed << std::setprecision(precision) << number;
  return(addContent(myNum.str()));
}

//
// Set new row to be used in addContent() method
//
pair<int, int> RootWTable::newLine() {
  m_serialRow++;
  m_serialCol=0;
  return std::make_pair(m_serialRow, m_serialCol);
}



#ifndef SUMMARYTABLE_H
#define SUMMARYTABLE_H

#include <map>
#include <string>

#include "global_funcs.hh"


/**
 * @class SummaryTable
 * @brief A generic object to build summary tables
 */

class SummaryTable {
public:
  SummaryTable() : numRows_(0), numColumns_(0), rowOffset_(0), columnOffset_(0), precision_(-1), summaryCellPosition_(0,0), summaryLabelPosition_(0,0) {};
  void setHeader(std::string rowHeader, std::string columnHeader, int rowOffset = 0, int columnOffset = 0) { // has to be called before filling the table with content or the row and column numbering will not be correctly set
    rowOffset_ = rowOffset; columnOffset_ = columnOffset;
    summaryTable[std::make_pair(0,0)] = columnHeader + " &rarr;<br>" + rowHeader + " &darr;";
  }
  void setPrecision(int precision) { precision_ = precision; } // has to be called before filling the table or conversions from floating point won't have the desired precision

  template<typename T> void setCell(int row, int column, const T& content) { setCell(row, column, any2str(content, precision_)); }
  template<typename T, typename BinaryOp> void setCell(int row, int column, const T& content, BinaryOp binop) { setCell(row, column, binop(hasCell(row, column) ? str2any<T>(getCell(row, column)) : T(), content)); }

  template<typename T> void setSummaryCell(std::string label, const T& content) { setSummaryCell(label, any2str(content, precision_)); }

  std::string getCell(int row, int column) { return summaryTable[std::make_pair(row,column)];} // this actually alters the map if the cell's not there = DANGEROUS

  bool hasCell(int row, int column) const { return summaryTable.count(std::make_pair(row,column)); }  // tests whether a cell has already been inserted = SAFE
  bool hasSummaryCell() const { return summaryCellPosition_ > std::make_pair(0, 0); }

  std::map<std::pair<int, int>, std::string>& getContent() { return summaryTable; }

  void clear() { summaryTable.clear(); }
private:
  std::map<std::pair<int, int>, std::string> summaryTable;
  int numRows_, numColumns_;
  int rowOffset_, columnOffset_; // from which number rows and columns headers should start
  int precision_; // precision to convert floating point numbers with
  std::pair<int, int> summaryCellPosition_, summaryLabelPosition_;
};


typedef std::map<std::string, SummaryTable> MultiSummaryTable;

#endif

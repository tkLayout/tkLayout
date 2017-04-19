#include "SummaryTable.hh"

template<> void SummaryTable::setCell<std::string>(const int row, const int column, const std::string& content) {
  if (column > 0 && !hasCell(0, column)) summaryTable[std::make_pair(0, column)] = any2str(column + columnOffset_);
  if (row > 0 && !hasCell(row, 0)) summaryTable[std::make_pair(row, 0)] = any2str(row + rowOffset_);
  summaryTable[std::make_pair(row, column)]=content;
  numRows_ = row+1 > numRows_ ? row+1 : numRows_;
  numColumns_ = column+1 > numColumns_ ? column+1 : numColumns_;
}

template<> void SummaryTable::setSummaryCell<std::string>(std::string label, const std::string& content) {
  if (!hasSummaryCell()) {
    if (numRows_ > 2 && numColumns_ > 2) {
      summaryLabelPosition_ = std::make_pair(numRows_, 0); 
      summaryCellPosition_ = std::make_pair(numRows_++, numColumns_++);
    } else if (numRows_ >= 2 && numColumns_ == 2) {
      summaryLabelPosition_ = std::make_pair(numRows_, 0); 
      summaryCellPosition_ = std::make_pair(numRows_++, 1);
    } else if (numRows_ == 2 && numColumns_ > 2) {
      summaryLabelPosition_ = std::make_pair(0, numColumns_); 
      summaryCellPosition_ = std::make_pair(1, numColumns_++);
    }
  }
  summaryTable[summaryLabelPosition_] = label;
  summaryTable[summaryCellPosition_] = content;
}

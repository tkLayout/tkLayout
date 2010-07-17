#ifndef MESSAGELOGGER_H
#define MESSAGELOGGER_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

class MessageLogger {
 public:
  ~MessageLogger();
  MessageLogger();
  MessageLogger(string newObjectName);
  bool addMessage(string message, int level=UNKNOWN);
  string getCompleteLog(int level);
  string getLatestLog(int level);
  // NumberOfLevels should always be the last here
  enum {UNKNOWN, ERROR, WARNING, INFO, DEBUG, NumberOfLevels};
  static string getLevelName(int level);
  bool hasEmptyLog(int level) { return (!wasModified[level]); };
  string getObjectName() { return objectName; };
 protected:
  ostringstream tempString;
 private:
  string objectName;
  bool wasModified[NumberOfLevels];
  ostringstream logStringStream[NumberOfLevels];
  ostringstream partialLogStringStream[NumberOfLevels];
};

#endif

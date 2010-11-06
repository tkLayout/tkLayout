#ifndef MESSAGELOGGER_H
#define MESSAGELOGGER_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

class LogMessage {
 public:
  LogMessage() {};
  ~LogMessage() {};
  string sender;
  int level;
  string message;
};

class MessageLogger {
 public:
  ~MessageLogger();
  MessageLogger();
  MessageLogger(string newObjectName);
  bool addMessage(string message, int level=UNKNOWN);
  bool addMessage(ostringstream& message, int level=UNKNOWN);
  //string getCompleteLog(int level);
  static string getLatestLog();
  static string getLatestLog(int level);
  // NumberOfLevels should always be the last here
  enum {UNKNOWN, ERROR, WARNING, INFO, DEBUG, NumberOfLevels};
  static std::string shortLevelCode[];
  static string getLevelName(int level);
  static bool hasEmptyLog(int level);
  string getObjectName() { return objectName; };
  string setObjectName(string newObjectName) { return objectName=newObjectName; };
 protected:
  ostringstream tempString;
 private:
  string objectName;
  //static bool wasModified[];
  static std::vector<LogMessage> logMessageV;
  static int countInstances;
  static int messageCounter[];
};

#endif

#ifndef MESSAGELOGGER_H
#define MESSAGELOGGER_H

#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <sstream>

#define logERROR(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::ERROR)
#define logWARNING(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::WARNING)
#define logINFO(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::INFO)
#define logDEBUG(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::DEBUG)

using namespace std;

class LogMessage {
 public:
  LogMessage() {};
  ~LogMessage() {};
  int level;
  string message;
};

class MessageLogger {
 public:
  static MessageLogger* instance();
  bool addMessage(string sourceFunction, string message, int level=UNKNOWN);
  bool addMessage(string sourceFunction, ostringstream& message, int level=UNKNOWN);
  static string getLatestLog();
  static string getLatestLog(int level);
  // NumberOfLevels should always be the last here
  enum {UNKNOWN, ERROR, WARNING, INFO, DEBUG, NumberOfLevels};
  static std::string shortLevelCode[];
  static string getLevelName(int level);
  static bool hasEmptyLog(int level);
  void setScreenLevel(int screenLevel) { screenLevel_ = screenLevel; }
 private:
  ~MessageLogger();
  MessageLogger();
  MessageLogger(MessageLogger const&){};
  static MessageLogger* myInstance_;
  static std::vector<LogMessage> logMessageV;
  static int countInstances;
  static int messageCounter[];
  int screenLevel_;
};

#endif

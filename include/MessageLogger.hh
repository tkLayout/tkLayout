#ifndef MESSAGELOGGER_H
#define MESSAGELOGGER_H

#include <iostream>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>
#include <sstream>

#define logERROR(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::ERROR)
#define logWARNING(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::WARNING)
#define logINFO(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::INFO)
#define logDEBUG(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::DEBUG)

#define logUniqueERROR(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::ERROR, MessageLogger::UNIQUE)
#define logUniqueWARNING(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::WARNING, MessageLogger::UNIQUE)
#define logUniqueINFO(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::INFO, MessageLogger::UNIQUE)
#define logUniqueDEBUG(message) MessageLogger::instance()->addMessage(__func__, message, MessageLogger::DEBUG, MessageLogger::UNIQUE)

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
  bool addMessage(string sourceFunction, string message, int level=UNKNOWN, bool unique=false);
  bool addMessage(string sourceFunction, ostringstream& message, int level=UNKNOWN, bool unique=false);
  static string getLatestLog();
  static string getLatestLog(int level);
  // NumberOfLevels should always be the last here
  enum {UNKNOWN, ERROR, WARNING, INFO, DEBUG, NumberOfLevels};
  static const bool UNIQUE = true;
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
  std::set<std::string> uniqueMessages;
};

#endif

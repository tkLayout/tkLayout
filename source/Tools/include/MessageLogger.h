#ifndef MESSAGELOGGER_H
#define MESSAGELOGGER_H

#include <iostream>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>
#include <sstream>

//! Log message with ERROR verbosity
#define logERROR(message)   MessageLogger::getInstance().addMessage(__func__, message, MessageLogger::ERROR)
//! Log message with WARNING verbosity
#define logWARNING(message) MessageLogger::getInstance().addMessage(__func__, message, MessageLogger::WARNING)
//! Log message with INFO verbosity
#define logINFO(message)    MessageLogger::getInstance().addMessage(__func__, message, MessageLogger::INFO)
//! Log message with DEBUG verbosity
#define logDEBUG(message)   MessageLogger::getInstance().addMessage(__func__, message, MessageLogger::DEBUG)

#define logUniqueERROR(message)   MessageLogger::getInstance().addMessage(__func__, message, MessageLogger::ERROR, MessageLogger::UNIQUE)
#define logUniqueWARNING(message) MessageLogger::getInstance().addMessage(__func__, message, MessageLogger::WARNING, MessageLogger::UNIQUE)
#define logUniqueINFO(message)    MessageLogger::getInstance().addMessage(__func__, message, MessageLogger::INFO, MessageLogger::UNIQUE)
#define logUniqueDEBUG(message)   MessageLogger::getInstance().addMessage(__func__, message, MessageLogger::DEBUG, MessageLogger::UNIQUE)

using namespace std;

/*
 * @class LogMessage
 * @brief Message class - container keeping message info + verbosity level
 */
class LogMessage {

 public:

  //! Constructor setting message + verbosity level
  LogMessage(int level, std::string message) : m_level(level), m_message(message) {};

  //! Destructor
  ~LogMessage() {};

  int         m_level;  //!< Verbosity level
  std::string m_message;//!< Message
};

/*
 * @class MessageLogger
 * @brief A singleton class providing general messaging interface. Individual messages are supposed to be logged through
 * various defined macros: logERROR, logWARNING, logINFO and logDEBUG, called at ERROR, WARNING, INFO and DEBUG verbosity
 * level, respectively. Use summary methods to get the overall log for given verbosity level.
 */
class MessageLogger {

 public:

  //! Message logger access method -> get instance of singleton class Message logger
  static MessageLogger& getInstance();

  //! Destructor
  ~MessageLogger();

  //! Add new message (string) to the logger of given level of verbosity
  bool addMessage(string sourceFunction, string message        , int level=UNKNOWN, bool unique=false);

  //! Add new message (stringstream) to the logger of given level of verbosity
  bool addMessage(string sourceFunction, ostringstream& message, int level=UNKNOWN, bool unique=false);

  //! Get summary log
  static std::string getSummaryLog();

  //! Get summary log of given verbosity level
  static std::string getSummaryLog(int level);

  //! Does interface has any log of given verbosity
  static bool hasEmptyLog(int level);

  // Verbosity levels (NumberOfLevels should always be the last here)
  enum {UNKNOWN, ERROR, WARNING, INFO, DEBUG, NumberOfLevels};

  //! Acronyms for verbosity level
  static std::string shortLevelCode[];

  //! Get corresponding verbosity level name
  static std::string getLevelName(int level);

  //! Set screen verbosity level
  void setScreenLevel(int screenLevel) { m_screenLevel = screenLevel; }

  //! Define unique message constant -> used for unique logging
  static const bool UNIQUE = true;

 private:

  //! Private constructor -> necessary for singleton pattern
  MessageLogger();

  static std::vector<LogMessage> s_logMessages;      //!< Container with all messages
  static int                     s_messageCounter[]; //!< Number of messages of given verbosity saved in the logger interface
  int                            m_screenLevel;      //!< Level to be printed to the output

  std::set<std::string>          m_uniqueMessages;   //!< Container of unique messages
};

#endif

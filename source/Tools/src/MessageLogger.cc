#include <MessageLogger.h>

//
// Acronyms for verbosity levele {UNKNOWN, ERROR, WARNING, INFO, DEBUG, NumberOfLevels};
//
std::string MessageLogger::shortLevelCode[] = { "??", "EE", "WW", "II", "DD" };

//
// Static instance of Container with all messages
//
std::vector<LogMessage> MessageLogger::s_logMessages;
int                     MessageLogger::s_messageCounter[NumberOfLevels];

//
// Message logger access method -> get instance of singleton class Message logger
//
MessageLogger& MessageLogger::getInstance() {

  static MessageLogger s_instance;
  return s_instance;
}

//
// Private constructor -> zero message counter for all levels
//
MessageLogger::MessageLogger() : m_screenLevel(ERROR)
{
  for (unsigned int i=0; i<NumberOfLevels; ++i) s_messageCounter[i]=0;
}

//
// Add new message (string) to the logger of given level of verbosity
//
bool MessageLogger::addMessage(std::string sourceFunction, std::string message, int level /*=UNKNOWN*/, bool unique /*=false*/ ) {

  if(unique) {

    if(m_uniqueMessages.count(message) == 0) {

      m_uniqueMessages.insert(message);
    } else return false;
  }

  // Print to the screen if level below or equal to the current screen verbosity level
  if (level<=m_screenLevel) std::cout << "(" + shortLevelCode[level] + ") " << sourceFunction<<": " << message << std::endl;

  // If DEBUG level, print function name, from which Logger interface called
  if (level==DEBUG) { // TODO: this could be more efficient

       for (int i=sourceFunction.length(); i<20; ++i) message=" "+ message;
       message = "[" + sourceFunction+"]: " + message;
  }

  // Save message in the logger for html output etc.
  if ((level>=0) && (level<NumberOfLevels)) {

    LogMessage newMessage(level, message);
    s_logMessages.push_back(newMessage);
    s_messageCounter[level]++;

    return true;
  }
  else return false;
}

//
// Add new message (stringstream) to the logger of given level of verbosity
//
bool MessageLogger::addMessage(std::string sourceFunction, std::ostringstream& message, int level /*=UNKNOWN*/, bool unique /*=false*/ ) {

  return addMessage(sourceFunction, message.str(), level, unique);
}

//
// Check that if log interface contains any message of given level
//
bool MessageLogger::hasEmptyLog(int level) {

  if ((level>=0) && (level<NumberOfLevels)) {
    return (s_messageCounter[level]==0);
  }
  return true;
}

//
// Get summary log of given verbosity level
//
std::string MessageLogger::getSummaryLog(int level) {

  std::string result="";

  if ((level>=0) && (level<NumberOfLevels)) {
    for (auto itMessage=s_logMessages.begin(); itMessage!=s_logMessages.end(); ) {

      if (itMessage->m_level==level) {

        result += (*itMessage).m_message+"\n";
      	s_messageCounter[itMessage->m_level]--;
        s_logMessages.erase(itMessage);
      }
      else ++itMessage;
    }
  }
  return result;
}

//
// Get summary log
//
std::string MessageLogger::getSummaryLog() {

  std::string result="";

  auto itMessage=s_logMessages.begin();
  while (itMessage!=s_logMessages.end()) {

    result += "(" + shortLevelCode[itMessage->m_level]+ ") " + itMessage->m_message+"\n";
    s_messageCounter[itMessage->m_level]--;
    s_logMessages.erase(itMessage);
    itMessage=s_logMessages.begin();
  }
  return result;
}

//
// Destructor
//
MessageLogger::~MessageLogger() {

    std::cerr << getSummaryLog();
}

//
// Get corresponding verbosity level name
//
std::string MessageLogger::getLevelName(int level) {
  switch (level) {
    case UNKNOWN: return "Unknown";
    case ERROR:   return "Error";
    case WARNING: return "Warning";
    case INFO:    return "Info";
    case DEBUG:   return "Debug";
    default:      return "InvalidCode";
  }
}

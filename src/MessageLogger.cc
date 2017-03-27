#include <MessageLogger.hh>

//bool MessageLogger::wasModified[MessageLogger::NumberOfLevels];
std::vector<LogMessage> MessageLogger::logMessageV;
int MessageLogger::countInstances = 0;

//  enum {UNKNOWN, ERROR, WARNING, INFO, DEBUG, NumberOfLevels};
std::string MessageLogger::shortLevelCode[] = { "??", "EE", "WW", "II", "DD" };
int MessageLogger::messageCounter[NumberOfLevels];

// Global static pointer used to ensure a single instance of the class.
MessageLogger* MessageLogger::myInstance_ = NULL;

// Returns the instance (if already present) or creates one if needed
MessageLogger* MessageLogger::instance() {
  return myInstance_ ? myInstance_ : (myInstance_ = new MessageLogger);
}

MessageLogger::MessageLogger() {
  screenLevel_ = ERROR;
  if (countInstances==0) {
    for (unsigned int i=0; i<NumberOfLevels; ++i)
      messageCounter[i]=0;
  }
  ++countInstances;
}

bool MessageLogger::addMessage(string sourceFunction, string message, int level /*=UNKNOWN*/, bool unique /*=false*/ ) {
  if(unique) {
    if(uniqueMessages.count(message) == 0) {
      uniqueMessages.insert(message);
    } else {
      return false;
    }
  }
    
  if (level<=screenLevel_) {
    std::cout << "(" + shortLevelCode[level]+ ") "
	      << sourceFunction<<": " << message << std::endl;
  }

  if (level==DEBUG) { // TODO: this could be more efficient
     for (int i=sourceFunction.length(); i<20; ++i) message=" "+ message;
     message = "[" + sourceFunction+"]: " + message;
  }

  if ((level>=0)&&(level<NumberOfLevels)) {
    LogMessage newMessage;
    newMessage.level=level;
    newMessage.message=message;
    logMessageV.push_back(newMessage);
    messageCounter[level]++;
    return true;
  } else return false;
}

bool MessageLogger::addMessage(string sourceFunction, ostringstream& message, int level /*=UNKNOWN*/, bool unique /*=false*/ ) {
  string newMessage = message.str();
  return addMessage(sourceFunction, newMessage, level, unique);
}

bool MessageLogger::hasEmptyLog(int level) {
  if ((level>=0)&&(level<NumberOfLevels)) {
    return (messageCounter[level]==0);
  }
  return true;
}

string MessageLogger::getLatestLog(int level) {
  string result="";
  if ((level>=0)&&(level<NumberOfLevels)) {
    std::vector<LogMessage>::iterator itMessage;
    //std::vector<LogMessage>::iterator nextMessage;
    for (itMessage=logMessageV.begin();
	 itMessage!=logMessageV.end(); ) {
      if (itMessage->level==level) {
        result += (*itMessage).message+"\n";
	//nextMessage=itMessage+1;
	messageCounter[itMessage->level]--;
        logMessageV.erase(itMessage);
	//itMessage=nextMessage;
      } else {
	++itMessage;
      }
    }
  }
  return result;
}

string MessageLogger::getLatestLog() {
  string result="";
  std::vector<LogMessage>::iterator itMessage=logMessageV.begin();
  while (itMessage!=logMessageV.end()) {
    result += "(" + shortLevelCode[itMessage->level]+ ") " + itMessage->message+"\n";
    messageCounter[itMessage->level]--;
    logMessageV.erase(itMessage);
    itMessage=logMessageV.begin();
  }
  return result;
}

MessageLogger::~MessageLogger() {
  if (--countInstances==0)
    std::cerr << getLatestLog();
}

string MessageLogger::getLevelName(int level) {
  switch (level) {
  case UNKNOWN: return "Unknown";
  case ERROR: return "Error";
  case WARNING: return "Warning";
  case INFO: return "Info";
  case DEBUG: return "Debug";
  default: return "InvalidCode";
  }
}

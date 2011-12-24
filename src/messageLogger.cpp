#include <messageLogger.h>

//bool MessageLogger::wasModified[MessageLogger::NumberOfLevels];
std::vector<LogMessage> MessageLogger::logMessageV;
int MessageLogger::countInstances = 0;

//  enum {UNKNOWN, ERROR, WARNING, INFO, DEBUG, NumberOfLevels};
std::string MessageLogger::shortLevelCode[] = { "??", "EE", "WW", "II", "DD" };
int MessageLogger::messageCounter[NumberOfLevels];

MessageLogger::MessageLogger() {
  objectName="unknownObject";
  if (countInstances==0) {
    for (unsigned int i=0; i<NumberOfLevels; ++i)
      messageCounter[i]=0;
  }
  ++countInstances;
}

MessageLogger::MessageLogger(string newObjectName) {
  objectName=newObjectName;
  ++countInstances;
}

bool MessageLogger::addMessage(string message, int level /*=UNKNOWN*/ ) {
  if ((level>=0)&&(level<NumberOfLevels)) {
    LogMessage newMessage;
    newMessage.sender=objectName;
    newMessage.level=level;
    newMessage.message=message;
    logMessageV.push_back(newMessage);
    messageCounter[level]++;
    return true;
  } else return false;
}

bool MessageLogger::addMessage(ostringstream& message, int level /*=UNKNOWN*/ ) {
  string newMessage = message.str();
  return addMessage(newMessage, level);
}

bool MessageLogger::hasEmptyLog(int level) {
  if ((level>=0)&&(level<NumberOfLevels)) {
    return (messageCounter[level]==0);
    
    /*
    for(unsigned int i=0; i<logMessageV.size(); ++i) {
      if (logMessageV[i].level==level) return false;
    }
    */

  }
  return true;
}

string MessageLogger::getLatestLog(int level) {
  string result="";
  if ((level>=0)&&(level<NumberOfLevels)) {
    std::vector<LogMessage>::iterator itMessage;
    std::vector<LogMessage>::iterator nextMessage;
    for (itMessage=logMessageV.begin();
	 itMessage!=logMessageV.end(); ) {
      // std::cerr << "==" << (itMessage->message).c_str() << "==" << std::endl;
      if (itMessage->level==level) {
        result += (*itMessage).sender+": "+(*itMessage).message+"\n";
	nextMessage=itMessage+1;
	messageCounter[itMessage->level]--;
        logMessageV.erase(itMessage);
	itMessage=nextMessage;
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
    result += "(" + shortLevelCode[itMessage->level]+ ") " + itMessage->sender+": "+ itMessage->message+"\n";
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

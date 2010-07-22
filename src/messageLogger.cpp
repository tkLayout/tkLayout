#include <messageLogger.h>

//bool MessageLogger::wasModified[MessageLogger::NumberOfLevels];
std::vector<LogMessage> MessageLogger::logMessageV[MessageLogger::NumberOfLevels];

MessageLogger::MessageLogger() {
  objectName="unknownObject";
}

MessageLogger::MessageLogger(string newObjectName) {
  objectName=newObjectName;
}

bool MessageLogger::addMessage(string message, int level /*=UNKNOWN*/ ) {
  if ((level>=0)&&(level<NumberOfLevels)) {
    LogMessage newMessage;
    newMessage.sender=objectName;
    newMessage.level=level;
    newMessage.message=message;
    logMessageV[level].push_back(newMessage);
    return true;
  } else return false;
}

bool MessageLogger::addMessage(ostringstream& message, int level /*=UNKNOWN*/ ) {
  string newMessage = message.str();
  return addMessage(newMessage, level);
}

bool MessageLogger::hasEmptyLog(int level) {
  bool result=false;
  if ((level>=0)&&(level<NumberOfLevels)) {
    if (logMessageV[level].size()==0) result=true;
  }
  return result;
}

string MessageLogger::getLatestLog(int level) {
  string result="";
  if ((level>=0)&&(level<NumberOfLevels)) {
    for (std::vector<LogMessage>::iterator itMessage=logMessageV[level].begin();
	 itMessage != logMessageV[level].end(); ++itMessage) {
      result += (*itMessage).sender+": "+(*itMessage).message+"\n";
      // std::cerr << "getLatestLog("<<level<<"):"<<(*itMessage).sender << ": " << (*itMessage).message << std::endl; //debug
    }
    //wasModified[level]=false;
    logMessageV[level].clear();
  }
  return result;
}

MessageLogger::~MessageLogger() {
  //for (unsigned int level=0; level<NumberOfLevels; ++level) {
  //  if (!hasEmptyLog(level)) {
  //    std::cerr << "Warning: unread messages for object " << objectName << " at level " << level << std::endl;
  //    std::cerr << getLatestLog(level);
  //  }
  //}
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

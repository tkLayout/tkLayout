#include <messageLogger.h>

//bool MessageLogger::wasModified[MessageLogger::NumberOfLevels];
std::vector<LogMessage> MessageLogger::logMessageV[MessageLogger::NumberOfLevels];

MessageLogger::MessageLogger() {
  objectName="unknownObject";
  for (unsigned int level=0; level<NumberOfLevels; ++level) {
    //wasModified[level]=false;
    //logStringStream[level].str("");
    //partialLogStringStream[level].str("");
  }
}

MessageLogger::MessageLogger(string newObjectName) {
  objectName=newObjectName;
  for (unsigned int level=0; level<NumberOfLevels; ++level) {
    //wasModified[level]=false;
    //logStringStream[level].str("");
    //partialLogStringStream[level].str("");
  }
}

bool MessageLogger::addMessage(string message, int level /*=UNKNOWN*/ ) {
  if ((level>=0)&&(level<NumberOfLevels)) {
    //logStringStream[level] << message;
    //logStringStream[level] << std::endl;
    LogMessage newMessage;
    newMessage.sender=objectName;
    newMessage.level=level;
    newMessage.message=message;
    //partialLogStringStream[level] << message;
    //partialLogStringStream[level] << std::endl;
    //wasModified[level]=true;
    logMessageV[level].push_back(newMessage);
    return true;
  } else return false;
}

/*
string MessageLogger::getCompleteLog(int level) {
  if ((level>=0)&&(level<NumberOfLevels)) {
    wasModified[level]=false;
    return logStringStream[level].str();
  } else return string("");  
  }*/

bool MessageLogger::hasEmptyLog(int level) {
  bool result=false;
  if ((level>=0)&&(level<NumberOfLevels)) {
    // std::cerr << "Size of log: " << logMessageV[level].size() << std::endl; // debug
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
  for (unsigned int level=0; level<NumberOfLevels; ++level) {
    if (!hasEmptyLog(level)) {
      std::cerr << "Warning: unread messages for object " << objectName << " at level " << level << std::endl;
      std::cerr << getLatestLog(level);
    }
  }
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

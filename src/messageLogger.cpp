#include <messageLogger.h>

MessageLogger::MessageLogger() {
  objectName="unknownObject";
  for (unsigned int level=0; level<NumberOfLevels; ++level) {
    wasModified[level]=false;
    logStringStream[level].str("");
    partialLogStringStream[level].str("");
  }
}

MessageLogger::MessageLogger(string newObjectName) {
  objectName=newObjectName;
  for (unsigned int level=0; level<NumberOfLevels; ++level) {
    wasModified[level]=false;
    logStringStream[level].str("");
    partialLogStringStream[level].str("");
  }
}

bool MessageLogger::addMessage(string message, int level /*=UNKNOWN*/ ) {
  if ((level>=0)&&(level<NumberOfLevels)) {
    logStringStream[level] << message;
    logStringStream[level] << std::endl;
    partialLogStringStream[level] << message;
    partialLogStringStream[level] << std::endl;
    wasModified[level]=true;
    return true;
  } else return false;
}

string MessageLogger::getCompleteLog(int level) {
  if ((level>=0)&&(level<NumberOfLevels)) {
    wasModified[level]=false;
    return logStringStream[level].str();
  } else return string("");  
}

string MessageLogger::getLatestLog(int level) {
  string result;
  if ((level>=0)&&(level<NumberOfLevels)) {
    wasModified[level]=false;
    result = partialLogStringStream[level].str();
    partialLogStringStream[level].str("");
    return result;
  } else return string("");  
}

MessageLogger::~MessageLogger() {
  for (unsigned int level=0; level<NumberOfLevels; ++level) {
    if (wasModified[level]) {
      std::cerr << "Warning: unread messages for object " << objectName << " at level " << level << std::endl;
      std::cerr << partialLogStringStream[level].str() << std::endl;
    }
  }
}

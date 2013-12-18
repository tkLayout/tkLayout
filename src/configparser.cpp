#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>
#include <cctype>
#include "configparser.hh"
//#include "module.hh"
//#include "layer.hh"
//#include "tracker.hh"
#include <stdlib.h>

using namespace std;
/*
configParser::configParser() {
  myTracker_= NULL;
}

configParser::~configParser() {

}
*/
string configParser::getTill(istream &inStream, char delimiter, bool singleWord, bool allowNothing=false) {
  string result="";

  if (getline(inStream, result, delimiter)) {
    inStream.get();
    if (singleWord) {
      string dummyString;
      istringstream tempStream(result.c_str());
      tempStream >> result;
      if (tempStream >> dummyString) {
        result="";
        cerr << "ERROR: unexpected extra word found: \""<<dummyString<<"\"";
      }
    }
  } else {
    result="";
    if (!allowNothing) {
      cerr << "ERROR: expecting " << delimiter << " and found end-of-file" << endl;
    }
  }

  return result;
}

// Parse a single parameter of the config file, dividing what's before the '=' to what's after
bool configParser::parseParameter(string &parameterName, string &parameterValue, istream& inStream) {
  string str;
  string myLine;

  myLine=getTill(inStream, ';', false, true);

  if (myLine!="") {
    istringstream myLineStr(myLine);
    parameterName=getTill(myLineStr, '=', true);
    if (parameterName!="") {
      if (myLineStr >> parameterValue) {
        string dummy;
        if (myLineStr >> dummy) {
          std::ostringstream tempString;
          tempString.str(""); tempString << "Ignoring extra parameter value " << dummy;
          logWARNING(tempString);
        }
        return true;
      }
    }
  }
  parameterName="";
  parameterValue="";
  return false;
}

bool configParser::breakParameterName(string& parameterName, int& ringIndex, int& diskIndex) {
  bool result = false;
  string tempParamStr;
  string tempParamInd;
  string dummy;
  int nIndexes;
  istringstream parameterStream(parameterName);
  diskIndex=0;

  getline(parameterStream, tempParamStr, '[');
  getline(parameterStream, tempParamInd, ']');
  getline(parameterStream, dummy);

  // We check that we found parameter[index] and nothing else
  if ((tempParamStr!="")&&(tempParamInd!="")&&(dummy=="")&&(parameterStream.eof())) {
    nIndexes=sscanf(tempParamInd.c_str(), "%d,%d", &ringIndex, &diskIndex);
    if (nIndexes==2) {
      result=true;
    } else if (nIndexes==1) {
      result=true;
      diskIndex=0;
    }
    parameterName=tempParamStr;
  }

  return result;
}

bool configParser::breakParameterName(string& parameterName, string& stringIndex) {
  bool result = false;
  string tempParamStr;
  string tempParamInd;
  string dummy;
  istringstream parameterStream(parameterName);
  stringIndex = "";

  getline(parameterStream, tempParamStr, '[');
  getline(parameterStream, tempParamInd, ']');
  getline(parameterStream, dummy);

  // We check that we found parameter[index] and nothing else
  if ((tempParamStr!="")&&(tempParamInd!="")&&(dummy=="")&&(parameterStream.eof())) {
    parameterName=tempParamStr;
    stringIndex = tempParamInd;
    result = true;
  }

  return result;
}


bool configParser::parseSupportParameters(std::istream& inStream, std::list<std::pair<int, double> >& plist) {
  std::string name, value;
  int index, dummy;
  double mid_z;
  while (!inStream.eof()) {
    while (parseParameter(name, value, inStream)) {
      std::pair<int, double> p;
      if (name.compare("midZ") == 0) index = 0;
      else {
        if (name.compare("midZ") > 0) {
          if (!breakParameterName(name, index, dummy)) return false;
        }
        else return false;
      }
      p.first = index;
      mid_z = atof(value.c_str());
      p.second = mid_z;
      if (mid_z >= 0) plist.push_back(p);
    }
  }
  return true;
}


std::list<std::pair<int, double> >* configParser::parseSupportsFromFile(string fileName) {
  std::list<std::pair<int, double> >* result = NULL;
  std::ifstream infilestream;
  infilestream.open(fileName.c_str());
  if (infilestream.is_open()) {
    result = new std::list<std::pair<int, double> >;
    std::stringstream filecontents;
    std::string lineorword;
    std::string::size_type locComm1;
    std::string::size_type locComm2;
    while (std::getline(infilestream, lineorword)) {
      locComm1 = lineorword.find("//", 0);
      locComm2 = lineorword.find("#", 0);
      if ((locComm1 == string::npos) && (locComm2 == string::npos)) {
        filecontents << lineorword << endl;
      }
      else {
        filecontents << lineorword.substr(0, (locComm1 < locComm2 ? locComm1 : locComm2)) << endl;
      }
    }
    try {
      while (filecontents >> lineorword) {
        if (lineorword == "Support") {
          getTill(filecontents, '{', true);
          std::string configparams = getTill(filecontents, '}', false, true);
          if (configparams != "") {
            std::istringstream paramstream(configparams);
            if (!parseSupportParameters(paramstream, *result)) {
              if (result) delete result;
              return NULL;
            }
          }
        }
      }
      return result;
    }
    catch (std::exception& e) {
      cerr << e.what() << endl;
      if (result) delete result;
      return NULL;
    }
  }
  if (result) delete result;
  return NULL;
}


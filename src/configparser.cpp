#include <fstream>
#include <sstream>
#include <iostream>
#include "configparser.hh"

using namespace std;

configParser::configParser() {
}

configParser::~configParser() {
}

string configParser::getTill(istream &inStream, char delimiter, bool singleWord,bool allowNothing=false) {
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
	  cerr << "WARNING: ignoring extra parameter value " << dummy << endl;
	}
	return true;
      }
    }
  } 
  parameterName="";
  parameterValue="";
  return false;
}

bool configParser::parseBarrel(string myName, istream& inStream) {
  string parameterName;
  string parameterValue;

  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      cout << "\t\"" << parameterName << "\"=\"" << parameterValue << "\";" << endl; // debug
    }
  }
}

bool configParser::parseEndcap(string myName, istream &inStream) {
  string parameterName;
  string parameterValue;

  while (parseParameter(parameterName, parameterValue, inStream)) {
    cout << "\t\"" << parameterName << "\"=\"" << parameterValue << "\";" << endl; // debug    
  }  
}

bool configParser::parseType(string myType) {
  string str;
  string typeConfig;

  if (myType=="Barrel") {
    cout << "CREATING BARREL:\t";
    str=getTill(configFile_, '{', true);
    if (str!="") { 
      cout << str << endl;
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
	istringstream typeStream(typeConfig);
	parseBarrel(str, typeStream);
      }
    }
  } else if (myType=="Endcap") {
    cout << "CREATING ENDCAP:\t";
    str=getTill(configFile_, '{', true);
    if (str!="") { 
      cout << str << endl;
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
	istringstream typeStream(typeConfig);
	parseEndcap(str, typeStream);
      }
    }
  } else {
    cerr << "Error: unknown piece of tracker " << myType;
    return false;
  }

  return true;
}

bool configParser::parseFile(string configFileName) {
  string str;

  if (configFile_.is_open()) {
    cerr << "config file is already open" << endl;
    return false;
  }

  configFile_.open(configFileName.c_str());
  if (configFile_.is_open()) {
    while(configFile_ >> str) {
      if (!parseType(str)) break;
    }
    configFile_.close();
  } else {
    cerr << "Error: could not open config file " << configFileName << endl;
    return false;
  }
  
  return true;
}

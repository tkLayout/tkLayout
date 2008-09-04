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
  
  if (myType=="Tracker") {
    cout << "Reading tracker main parameters and name: ";
    str=getTill(configFile_, '{', true);
    if (str!="") { 
      cout << str << endl;
      typeConfig=getTill(configFile_, '}', false);
      //if (typeConfig!="") {
      //	istringstream typeStream(typeConfig);
      //	parseBarrel(str, typeStream);
      //}
    }
  } else if (myType=="Barrel") {
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

  if (rawConfigFile_.is_open()) {
    cerr << "config file is already open" << endl;
    return false;
  }

  rawConfigFile_.open(configFileName.c_str());
  if (rawConfigFile_.is_open()) {

    // Reset the configFile_ object
    configFile_.~stringstream();                      
    new ( (void*) &configFile_) std::stringstream();

    // Skim comments delimited by // or # and put the skimmed file into configFile_
    string::size_type locComm1; // Location of commenting substring 1
    string::size_type locComm2; // Location of commenting substring 2

    while (getline(rawConfigFile_,str)) {
      locComm1=str.find("//", 0);
      locComm2=str.find("#", 0);
      if ((locComm1==string::npos)&&(locComm2==string::npos)) {
	configFile_ << str << endl;
      } else {
	configFile_ << str.substr(0, (locComm1<locComm2 ? locComm1 : locComm2)) << endl;
      }
    }
    
    
    while(configFile_ >> str) {
      if (!parseType(str)) break;
    }
    rawConfigFile_.close();
  } else {
    cerr << "Error: could not open config file " << configFileName << endl;
    return false;
  }
  
  return true;
}

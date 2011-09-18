#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

#include <stdio.h>
#include <stdlib.h>

#include <boost/regex.hpp>
#include <boost/filesystem/operations.hpp>

#include <sys/types.h>
#include <regex.h>

#include <mainConfigHandler.h>

using namespace std;
using namespace boost;

template <class T> bool from_string(T& t, const std::string& s, 
				    std::ios_base& (*f)(std::ios_base&)) {
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

mainConfigHandler::mainConfigHandler() {
  goodConfigurationRead_ = false;
  //styleDirectory_ = "";
  layoutDirectory_ = "";
  standardDirectory_ = "";
}

bool mainConfigHandler::checkDirectory(string dirName) {
  if (! filesystem::exists(dirName)) {
    cout << "Directory '" << dirName << "' does not exist!" << endl;
    return false;
  }
  if (! filesystem::is_directory(dirName) ) {
    cout << "Directory '" << dirName << "' is not a directory!" << endl;
    return false;    
  }

  return true;
}

bool mainConfigHandler::createConfigurationFileFromQuestions(string& configFileName) {
  string tempString;

  // Clear screen
  //cout << "\033[2J"; // Clears the screen
  //cout << "\033[1;1H"; // Places cursor on line 1

  // I have no configuration, so I must create it
  cout << "Could not find the configuration file "  << configFileName 
       << " maybe this is the first time you run with the new system." << endl;
  cout << "Answer to the following questions to have your configuration file automatically created." << endl;
  cout << "You will be later able to edit it manually, or you can just delete it and answer these questions again." << endl;
  cout << endl;

  cout << "*** What is the web server directory where you want to" << endl
       << "    place your output?" << endl
       << "    Example: " << getenv(HOMEDIRECTORY) << "/www/layouts : ";
  cin >> layoutDirectory_;
  if (!checkDirectory(layoutDirectory_)) return false;
  cout << endl;

  cout << "*** What is the standard output directory?" << endl
       << "    xml files and other various output will be put here" << endl
       << "    Example: " << getenv(HOMEDIRECTORY) << "/tkgeometry : ";
  cin >> standardDirectory_;
  cout << endl;

  cout << "*** Specify the list of transverse momenta to be used for the" << endl
       << "    tracking performance test (in GeV/c)" << endl
       << "    Example: 1, 10, 100 : ";
  cin >> tempString;
  string tempString2;
  getline(cin,tempString2);
  tempString+=tempString2;
  momenta_ = parseDoubleList(tempString);

  tempString = "";
  tempString2 = "";
  cout << "*** Specify the list of transverse momenta to be used for the" << endl
       << "    trigger efficiency performance test (in GeV/c)" << endl
       << "    Example: 1, 2, 5, 10 : ";
  cin >> tempString;
  getline(cin,tempString2);
  tempString+=tempString2;
  triggerMomenta_ = parseDoubleList(tempString2);

  ofstream configFile;
  configFile.open(configFileName.c_str(), ifstream::out);
  if (!configFile.good()) {
    cout << "Could not open " << configFile << " for writing. I quit."  << endl;
    configFile.close();
    return false;
  } else {
    configFile << LAYOUTDIRECTORYDEFINITION << "=\"" << layoutDirectory_ << "\"" << endl;
    configFile << STANDARDDIRECTORYDEFINITION << "=\"" << standardDirectory_ << "\"" << endl;

    configFile << MOMENTADEFINITION << "=\"";
    for (std::vector<double>::iterator it = momenta_.begin(); it!=momenta_.end(); ++it) {
      if (it!=momenta_.begin()) configFile << ", ";
      configFile << std::fixed << std::setprecision(2) << (*it);
    }
    configFile << "\"" << std::endl;

    configFile << TRIGGERMOMENTADEFINITION << "=\"";
    for (std::vector<double>::iterator it = triggerMomenta_.begin(); it!=triggerMomenta_.end(); ++it) {
      if (it!=triggerMomenta_.begin()) configFile << ", ";
      configFile << std::fixed << std::setprecision(2) << (*it);
    }
    configFile << "\"" << std::endl;

    configFile.close();
  }

  return true;
}


bool mainConfigHandler::parseLine(const char* codeLine, string& parameter, string& value) {
  cmatch what;
  regex parseLineExpression("^[ \\t]*([a-zA-Z0-9_]*)=\"([^\"]+)\".*");
  regex parseLineEmpty("^[ \\t]*");

  if (regex_match(codeLine, what, parseLineExpression)) {
    // what[1] contains the parameter name
    // what[2] contains the parameter value
    parameter = what[1];
    value = what[2];
    return true;
  } else if (regex_match(codeLine, what, parseLineEmpty)) {
    // Empty line, ok
    parameter="";
    value="";
    return false;
  } else {
    // Wrong formatting
    cerr << "Cannot understand line: '" << codeLine << "' in the configuration file " << CONFIGURATIONFILENAME << endl;
    return false;
  }
}

vector<double> mainConfigHandler::parseDoubleList(string inString) {
  cmatch what;
  double myDouble;
  vector<double> result;
  regex parseExpression("[^0-9\\.]*([0-9\\.]*)(.*)");

  int escapeCounter=100;
  while ((regex_match(inString.c_str(), what, parseExpression)&&(escapeCounter))) {
    escapeCounter--;
    if (from_string<double>(myDouble, string(what[1]), std::dec)) {
      result.push_back(myDouble);
      inString = what[2];
    } else {
      escapeCounter=0;
    }
  }

  return result;
}

bool mainConfigHandler::readConfigurationFile(ifstream& configFile) {
  char myLine[1024];
  string parameter, value;
  //bool styleFound=false;
  bool layoutFound=false;
  bool xmlFound=false;
  bool momentaFound=false;
  bool triggerMomentaFound=false;

  // Parsing all the lines of the configuration file
  while (configFile.good()) {
    configFile.getline(myLine, 1024);
    if (parseLine(myLine, parameter, value)) {
      // Case insensitive: the configuration file should be a shell script too
      //std::transform(parameter.begin(), parameter.end(), parameter.begin(), ::tolower);
      // If he's defining the style directory, then we map it to styleDirectory_
      //if (parameter==STYLEDIRECTORYDEFINITION) {
      //styleDirectory_ = value;
      //styleFound=true;
      //} 
      if (parameter==LAYOUTDIRECTORYDEFINITION) {
	layoutDirectory_ = value;
	layoutFound = true;
      } else if (parameter==STANDARDDIRECTORYDEFINITION) {
	standardDirectory_ = value;
	xmlFound = true;
      } else if (parameter==MOMENTADEFINITION) {
	momenta_ = parseDoubleList(value);
	momentaFound = true;
      } else if (parameter==TRIGGERMOMENTADEFINITION) {
	triggerMomenta_ = parseDoubleList(value);
	triggerMomentaFound = true;
      } else {
	cerr << "Unknown parameter " << parameter << " in the configuration file " << CONFIGURATIONFILENAME << endl;
      }
    }
  }

  //return (styleFound&&layoutFound&&xmlFound&&momentaFound);
  return (layoutFound&&xmlFound&&momentaFound&&triggerMomentaFound);
}

bool mainConfigHandler::getConfiguration(bool checkDirExists /* = true */) {
  return readConfiguration(checkDirExists);
}

bool mainConfigHandler::getConfiguration(string& layoutDirectory) {
  bool result = readConfiguration(true);
  if (result) {
    //styleDirectory = styleDirectory_;
    layoutDirectory = layoutDirectory_;
  }
  return result;
}

string mainConfigHandler::getConfigFileName() {  
  string homeDirectory = string(getenv(HOMEDIRECTORY));
  return homeDirectory+"/"+CONFIGURATIONFILENAME;
}

bool mainConfigHandler::readConfiguration( bool checkDirExists ) {
  if (goodConfigurationRead_) return true;

  ifstream configFile;
  string homeDirectory = string(getenv(HOMEDIRECTORY));
  string configFileName = getConfigFileName();
  bool goodConfig=false;

  configFile.open(configFileName.c_str(), ifstream::in);
  if (!configFile.good()) {
    configFile.close();
    goodConfig = createConfigurationFileFromQuestions(configFileName);
  } else {
    // I will read the configuration out of that
    goodConfig = readConfigurationFile(configFile);
    configFile.close();
    if (goodConfig) {
      if (checkDirExists) {
	// Check the basic configuration directories
	if (!checkDirectory(layoutDirectory_)) {
	  cout << "You probably need to edit or delete the configuration file " << CONFIGURATIONFILENAME << endl;
	  return false;
	}
	if (!checkDirectory(standardDirectory_)) {
	  cout << "You probably need to edit or delete the configuration file " << CONFIGURATIONFILENAME << endl;
	  return false;
	}
	// Check the mandatory subdirectories in the main
	if (!checkDirectory(getXmlDirectory_())) return false;
	if (!checkDirectory(getMattabDirectory_())) return false;
	if (!checkDirectory(getRootfileDirectory_())) return false;
	if (!checkDirectory(getGraphDirectory_())) return false;
	if (!checkDirectory(getSummaryDirectory_())) return false;
      }
    } else { // not good config read
      cout << "Configuration file '" << configFileName << "' not properly formatted. You probably need to edit or delete it" << endl;
      return false;
    }
  }

  goodConfigurationRead_ = goodConfig;
  return goodConfig;
}

vector<double>& mainConfigHandler::getMomenta() {
  getConfiguration();
  return momenta_;
}

vector<double>& mainConfigHandler::getTriggerMomenta() {
  getConfiguration();
  return triggerMomenta_;
}

string mainConfigHandler::getLayoutDirectory() {
  getConfiguration();
  return getLayoutDirectory_();
}

string mainConfigHandler::getStandardDirectory() {
  getConfiguration();
  return getStandardDirectory_();
}

string mainConfigHandler::getStyleDirectory() {
  getConfiguration();
  return getStyleDirectory_();
}

string mainConfigHandler::getXmlDirectory() {
  getConfiguration();
  return getXmlDirectory_();
}

string mainConfigHandler::getMattabDirectory() {
  getConfiguration();
  return getMattabDirectory_();
}

string mainConfigHandler::getRootfileDirectory() {
  getConfiguration();
  return getRootfileDirectory_();
}

string mainConfigHandler::getGraphDirectory() {
  getConfiguration();
  return getGraphDirectory_();
}

string mainConfigHandler::getSummaryDirectory() {
  getConfiguration();
  return getSummaryDirectory_();
}

string mainConfigHandler::getLayoutDirectory_() { return layoutDirectory_; }
string mainConfigHandler::getStandardDirectory_() { return standardDirectory_; }
string mainConfigHandler::getStyleDirectory_() { return layoutDirectory_+"/"+insur::default_styledir; }
string mainConfigHandler::getXmlDirectory_() { return standardDirectory_+"/"+insur::default_xmlpath; } 
string mainConfigHandler::getMattabDirectory_() { return standardDirectory_+"/"+insur::default_mattabdir; }
string mainConfigHandler::getRootfileDirectory_() { return standardDirectory_+"/"+insur::default_rootfiledir; }
string mainConfigHandler::getGraphDirectory_() { return standardDirectory_+"/"+insur::default_graphdir; }
string mainConfigHandler::getSummaryDirectory_() { return standardDirectory_+"/"+insur::default_summarypath; }


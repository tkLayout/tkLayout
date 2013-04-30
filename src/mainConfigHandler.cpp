#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>

#include <boost/filesystem/operations.hpp>

#include <sys/types.h>

#include <mainConfigHandler.h>
#include <global_funcs.h>

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

  tempString = "";
  tempString2 = "";
  cout << "*** Specify the list of trigger efficiency to be used for the" << endl
       << "    pt threshold find test (in percent: write 100 for full efficiency)" << endl
       << "    Example: 1, 50, 90, 95 : ";
  cin >> tempString;
  getline(cin,tempString2);
  tempString+=tempString2;
  thresholdProbabilities_ = parseDoubleList(tempString2);
  for (vector<double>::iterator it = thresholdProbabilities_.begin();
       it!=thresholdProbabilities_.end(); ++it) (*it)/=100;

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

    configFile << THRESHOLDPROBABILITIESDEFINITION << "=\"";
    for (std::vector<double>::iterator it = thresholdProbabilities_.begin(); it!=thresholdProbabilities_.end(); ++it) {
      if (it!=thresholdProbabilities_.begin()) configFile << ", ";
      configFile << std::fixed << std::setprecision(2) << (*it);
    }
    configFile << "\"" << std::endl;

    configFile.close();
  }

  return true;
}


bool mainConfigHandler::parseLine(const char* codeLine, string& parameter, string& value) {
  std::vector<string> tokens = split(codeLine, "=");
  if (tokens.empty()) {
    parameter = "";
    value = "";
    return false;
  } else if (tokens.size() < 2) { 
    cerr << "Cannot understand line: '" << codeLine << "' in the configuration file " << CONFIGURATIONFILENAME << endl;
    return false;
  } else {
    parameter = trim(tokens.at(0), " \"\n\t");
    value = trim(tokens.at(1), " \"\n\t");
    return true;
  }
}

vector<double> mainConfigHandler::parseDoubleList(string inString) {
  std::vector<std::string> tokens = split(inString, ",");
  std::vector<double> result;
  std::transform(tokens.begin(), tokens.end(), std::back_inserter(result), str2any<double>);

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
  bool thresholdProbabilitiesFound=false;

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
      } else if (parameter==THRESHOLDPROBABILITIESDEFINITION) {
	thresholdProbabilities_ = parseDoubleList(value);
	thresholdProbabilitiesFound = true;
      } else {
	cerr << "Unknown parameter " << parameter << " in the configuration file " << CONFIGURATIONFILENAME << endl;
      }
    }
  }

  //return (styleFound&&layoutFound&&xmlFound&&momentaFound);
  return (layoutFound&&xmlFound&&momentaFound&&triggerMomentaFound&&thresholdProbabilitiesFound);
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
  char* specialConfigFile = getenv(CONFIGURATIONFILENAMEDEFINITION);
  if (specialConfigFile) return string(specialConfigFile);
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

vector<double>& mainConfigHandler::getThresholdProbabilities() {
  getConfiguration();
  return thresholdProbabilities_;
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

string mainConfigHandler::getIrradiationDirectory() {
  getConfiguration();
  return getIrradiationDirectory_();
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

string mainConfigHandler::getDefaultMaterialsDirectory() {
  getConfiguration();
  return getDefaultMaterialsDirectory_();
}

string mainConfigHandler::getLayoutDirectory_() { return layoutDirectory_; }
string mainConfigHandler::getStandardDirectory_() { return standardDirectory_; }
string mainConfigHandler::getStyleDirectory_() { return layoutDirectory_+"/"+insur::default_styledir; }
string mainConfigHandler::getXmlDirectory_() { return standardDirectory_+"/"+insur::default_xmlpath; } 
string mainConfigHandler::getMattabDirectory_() { return standardDirectory_+"/"+insur::default_mattabdir; }
string mainConfigHandler::getIrradiationDirectory_() { return standardDirectory_+"/"+insur::default_irradiationdir; }
string mainConfigHandler::getRootfileDirectory_() { return standardDirectory_+"/"+insur::default_rootfiledir; }
string mainConfigHandler::getGraphDirectory_() { return standardDirectory_+"/"+insur::default_graphdir; }
string mainConfigHandler::getSummaryDirectory_() { return standardDirectory_+"/"+insur::default_summarypath; }
string mainConfigHandler::getDefaultMaterialsDirectory_() { return standardDirectory_+"/"+insur::default_materialsdir; }


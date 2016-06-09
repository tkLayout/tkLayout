/**
 * @file MainConfigHandler.cpp
 * @brief implementation of class that take care of creating and reading the configuration file .tkgeometryrc
 */
#include <MainConfigHandler.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

#include <stdio.h>
#include <stdlib.h>

#include <boost/filesystem/operations.hpp>
#include <global_constants.h>
#include <global_funcs.h>

#include <sys/types.h>

using namespace std;
using namespace boost;

template <class T> bool from_string(T& t, const std::string& s, 
                                    std::ios_base& (*f)(std::ios_base&)) {
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

MainConfigHandler::MainConfigHandler() {
  goodConfigurationRead_ = false;
  //styleDirectory_ = "";
  binDirectory_ = "";
  layoutDirectory_ = "";
  standardDirectory_ = "";
}

MainConfigHandler& MainConfigHandler::getInstance() {

  static MainConfigHandler s_instance;
  return s_instance;
}

bool MainConfigHandler::checkDirectory(string dirName) {
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

void MainConfigHandler::askBinDirectory() {
  cout << "*** What is the bin directory where you want to" << endl
      << "    place your executables?" << endl
      << "    Example: " << getenv(HOMEDIRECTORY) << "/bin : ";
  cin >> binDirectory_;
}

void MainConfigHandler::askLayoutDirectory() {
  cout << "*** What is the web server directory where you want to" << endl
    << "    place your output?" << endl
    << "    Example: " << getenv(HOMEDIRECTORY) << "/www/layouts : ";
  cin >> layoutDirectory_;
}

void MainConfigHandler::askStandardDirectory() {
  cout << "*** What is the standard output directory?" << endl
      << "    xml files and other various output will be put here" << endl
      << "    Example: " << getenv(HOMEDIRECTORY) << "/tkgeometry : ";
  cin >> standardDirectory_;
}

void MainConfigHandler::askMomenta() {
  string tempString = "";
  string tempString2 = "";

  cout << "*** Specify the list of transverse momenta to be used for the" << endl
      << "    tracking performance test (in GeV/c)" << endl
      << "    Example: 1, 10, 100 : ";
  cin >> tempString;

  getline(cin,tempString2);
  tempString+=tempString2;
  momenta_ = parseDoubleList(tempString);
}

void MainConfigHandler::askTriggerMomenta() {
  string tempString = "";
  string tempString2 = "";

  cout << "*** Specify the list of transverse momenta to be used for the" << endl
      << "    trigger efficiency performance test (in GeV/c)" << endl
      << "    Example: 1, 2, 5, 10 : ";
  cin >> tempString;
  getline(cin,tempString2);
  tempString+=tempString2;
  triggerMomenta_ = parseDoubleList(tempString);
}

void MainConfigHandler::askThresholdProbabilities() {
  string tempString = "";
  string tempString2 = "";

  cout << "*** Specify the list of trigger efficiency to be used for the" << endl
      << "    pt threshold find test (in percent: write 100 for full efficiency)" << endl
      << "    Example: 1, 50, 90, 95 : ";
  cin >> tempString;
  getline(cin,tempString2);
  tempString+=tempString2;
  thresholdProbabilities_ = parseDoubleList(tempString);
  for (vector<double>::iterator it = thresholdProbabilities_.begin();
      it!=thresholdProbabilities_.end(); ++it) (*it)/=100;
}


bool MainConfigHandler::createConfigurationFileFromQuestions(string& configFileName) {

  // Clear screen
  //cout << "\033[2J"; // Clears the screen
  //cout << "\033[1;1H"; // Places cursor on line 1

  // I have no configuration, so I must create it
  cout << "Could not find the configuration file "  << configFileName
    << " maybe this is the first time you run with the new system." << endl;
  cout << "Answer to the following questions to have your configuration file automatically created." << endl;
  cout << "You will be later able to edit it manually, or you can just delete it and answer these questions again." << endl;
  cout << endl;

  askBinDirectory();
  if (!checkDirectory(binDirectory_)) return false;
  cout << endl;

  askLayoutDirectory();
  if (!checkDirectory(layoutDirectory_)) return false;
  cout << endl;

  askStandardDirectory();
  cout << endl;

  askMomenta();

  askTriggerMomenta();

  askThresholdProbabilities();

  ofstream configFile;
  configFile.open(configFileName.c_str(), ifstream::out);
  if (!configFile.good()) {
    cout << "Could not open " << configFileName << " for writing. I quit."  << endl;
    configFile.close();
    return false;
  } else {
    configFile << BINDIRECTORYDEFINITION << "=\"" << binDirectory_ << "\"" << endl;
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


bool MainConfigHandler::parseLine(const char* codeLine, string& parameter, string& value) {
  std::vector<string> tokens = split(codeLine, "=");
  if (tokens.empty()) {
    parameter = "";
    value = "";
    return false;
  } else if (tokens.size() < 2) { 
    cerr << "Cannot understand line: '" << codeLine << "' in the configuration file " << CONFIGURATIONFILENAME << endl;
    return false;
  } else {
    parameter = ctrim(tokens.at(0), " \"\n\t");
    value = ctrim(tokens.at(1), " \"\n\t");
    return true;
  }
}

vector<double> MainConfigHandler::parseDoubleList(string inString) {
  return split<double>(inString, ",");
}

bool MainConfigHandler::readConfigurationFile(string& configFileName) {
  char myLine[1024];
  string parameter, value;
  //bool styleFound=false;
  bool binFound=false;
  bool layoutFound=false;
  bool xmlFound=false;
  bool momentaFound=false;
  bool triggerMomentaFound=false;
  bool thresholdProbabilitiesFound=false;
  ifstream configFileIs;
  ofstream configFileOs;

  configFileIs.open(configFileName.c_str(), ifstream::in);

  // Parsing all the lines of the configuration file
  while (configFileIs.good()) {
    configFileIs.getline(myLine, 1024);
    if (parseLine(myLine, parameter, value)) {
      // Case insensitive: the configuration file should be a shell script too
      //std::transform(parameter.begin(), parameter.end(), parameter.begin(), ::tolower);
      // If he's defining the style directory, then we map it to styleDirectory_
      //if (parameter==STYLEDIRECTORYDEFINITION) {
      //styleDirectory_ = value;
      //styleFound=true;
      //} 
      if (parameter==BINDIRECTORYDEFINITION) {
        binDirectory_ = value;
        binFound = true;
      } else if (parameter==LAYOUTDIRECTORYDEFINITION) {
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
  configFileIs.close();

  //ask here for bin directory for backward compatibility of the install script
  if(!binFound&&layoutFound&&xmlFound&&momentaFound&&triggerMomentaFound&&thresholdProbabilitiesFound){
    askBinDirectory();
    if (checkDirectory(binDirectory_)) {
      configFileOs.open(configFileName.c_str(), ifstream::app);
      if (!configFileOs.good()) {
        cout << "Could not open " << configFileName << " for writing. I quit."  << endl;
        configFileOs.close();
        return false;
      } else {
        configFileOs << BINDIRECTORYDEFINITION << "=\"" << binDirectory_ << "\"" << endl;

        configFileOs.close();
        binFound = true;
      }
    }
    cout << endl;
  }

  //return (styleFound&&layoutFound&&xmlFound&&momentaFound);
  return (binFound&&layoutFound&&xmlFound&&momentaFound&&triggerMomentaFound&&thresholdProbabilitiesFound);
}

bool MainConfigHandler::getConfiguration(bool checkDirExists /* = true */) {
  return readConfiguration(checkDirExists);
}

bool MainConfigHandler::getConfiguration(string& layoutDirectory) {
  bool result = readConfiguration(true);
  if (result) {
    //styleDirectory = styleDirectory_;
    layoutDirectory = layoutDirectory_;
  }
  return result;
}

string MainConfigHandler::getConfigFileName() {
  char* specialConfigFile = getenv(CONFIGURATIONFILENAMEDEFINITION);
  if (specialConfigFile) return string(specialConfigFile);
  string homeDirectory = string(getenv(HOMEDIRECTORY));
  return homeDirectory+"/"+CONFIGURATIONFILENAME;
}

bool MainConfigHandler::readConfiguration( bool checkDirExists ) {
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
    configFile.close();
    // I will read the configuration out of that
    goodConfig = readConfigurationFile(configFileName);
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
      }
    } else { // not good config read
      cout << "Configuration file '" << configFileName << "' not properly formatted. You probably need to edit or delete it" << endl;
      return false;
    }
  }

  goodConfigurationRead_ = goodConfig;
  return goodConfig;
}

vector<double>& MainConfigHandler::getMomenta() {
  getConfiguration();
  return momenta_;
}

vector<double>& MainConfigHandler::getTriggerMomenta() {
  getConfiguration();
  return triggerMomenta_;
}

vector<double>& MainConfigHandler::getThresholdProbabilities() {
  getConfiguration();
  return thresholdProbabilities_;
}

string MainConfigHandler::getBinDirectory() {
  getConfiguration();
  return getBinDirectory_();
}

string MainConfigHandler::getLayoutDirectory() {
  getConfiguration();
  return getLayoutDirectory_();
}

string MainConfigHandler::getStandardDirectory() {
  getConfiguration();
  return getStandardDirectory_();
}

string MainConfigHandler::getStyleDirectory() {
  getConfiguration();
  return getStyleDirectory_();
}

string MainConfigHandler::getXmlDirectory() {
  getConfiguration();
  return getXmlDirectory_();
}

string MainConfigHandler::getMattabDirectory() {
  getConfiguration();
  return getMattabDirectory_();
}

string MainConfigHandler::getIrradiationDirectory() {
  getConfiguration();
  return getIrradiationDirectory_();
}

string MainConfigHandler::getDefaultMaterialsDirectory() {
  getConfiguration();
  return getDefaultMaterialsDirectory_();
}

string MainConfigHandler::getStandardIncludeDirectory() {
  getConfiguration();
  return getStandardIncludeDirectory_();
}

string MainConfigHandler::getGeometriesDirectory() {
  getConfiguration();
  return getGeometriesDirectory_();
}

string MainConfigHandler::getBinDirectory_() { return binDirectory_; }
string MainConfigHandler::getLayoutDirectory_() { return layoutDirectory_; }
string MainConfigHandler::getStandardDirectory_() { return standardDirectory_; }
string MainConfigHandler::getStyleDirectory_() { return layoutDirectory_+"/"+default_styledir; }
string MainConfigHandler::getXmlDirectory_() { return standardDirectory_+"/"+default_xmlpath; }
string MainConfigHandler::getMattabDirectory_() { return standardDirectory_+"/"+default_mattabdir; }
string MainConfigHandler::getIrradiationDirectory_() { return standardDirectory_+"/"+default_irradiationdir; }
string MainConfigHandler::getDefaultMaterialsDirectory_() { return standardDirectory_+"/"+default_materialsdir; }
string MainConfigHandler::getStandardIncludeDirectory_() { return standardDirectory_+"/"+default_configdir+"/"+default_stdincludedir; }
string MainConfigHandler::getGeometriesDirectory_() { return standardDirectory_+"/"+default_geometriesdir; }


std::set<string> MainConfigHandler::preprocessConfiguration(istream& is, ostream& os, const string& istreamid) {
  using namespace std;
  string line;
  int numLine = 1;
  std::set<string> includeSet;
  includeSet.insert(istreamid);
  while(getline(is, line).good()) {
    if (line.find("//") != string::npos) line = line.erase(line.find("//"));
    string trimmed = trim(line);
    int includeStart;
    if ((includeStart = trimmed.find("@include")) != string::npos) { //@include @include-std @include-weak @include-std-weak
      trimmed = trimmed.substr(includeStart);
      int quoteStart, quoteEnd;
      string filename;
      if ((quoteStart = trimmed.find_first_of("\"")) != string::npos && (quoteEnd = trimmed.find_last_of("\"")) != string::npos) {
        filename = ctrim(trimmed.substr(quoteStart, quoteEnd - quoteStart + 1), "\"");
      } else {
        auto tokens = split(trimmed, " ");
        filename = tokens.size() > 1 ? tokens[1] : "";
      }
      bool includeStdOld = trimmed.find("@includestd") != string::npos;  // both @includestd (deprecated) and @include-std (preferred) are supported 
      bool includeStdNew = trimmed.find("@include-std") != string::npos;
     // bool includeWeak = trimmed.find("@include-weak") != string::npos || trimmed.find("@include-std-weak") != string::npos || trimmed.find("@includestd-weak") != string::npos; // include weak command not supported for the moment
      string prefix = (includeStdOld || includeStdNew) ? getStandardIncludeDirectory()+"/" : std::string("");
      filename = prefix + filename;
      ifstream ifs(filename);
      //std::cout << "INFO: " << filename << std::endl;
      if (ifs) {
        stringstream ss;
        auto&& moreIncludes = preprocessConfiguration(ifs, ss, filename);
        includeSet.insert(moreIncludes.begin(), moreIncludes.end());
        string indent = line.substr(0, line.find_first_not_of(" \t"));
        while (getline(ss, line).good()) {
          os << indent << line << endl;   
        }
      } else {
        cerr << "WARNING: " << istreamid << ":" << numLine << ": Ignoring malformed @include or @includestd directive" << endl;
      }
    } else {
      os << line << endl;
    }
    numLine++;
  }
  return includeSet;
}

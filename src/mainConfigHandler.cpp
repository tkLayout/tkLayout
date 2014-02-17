#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

#include <stdio.h>
#include <stdlib.h>

#include <boost/filesystem/operations.hpp>

#include <sys/types.h>

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
  binDirectory_ = "";
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

  cout << "*** What is the bin directory where you want to" << endl
      << "    place your executables?" << endl
      << "    Example: " << getenv(HOMEDIRECTORY) << "/bin : ";
    cin >> binDirectory_;
    if (!checkDirectory(binDirectory_)) return false;
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
  triggerMomenta_ = parseDoubleList(tempString);

  tempString = "";
  tempString2 = "";
  cout << "*** Specify the list of trigger efficiency to be used for the" << endl
    << "    pt threshold find test (in percent: write 100 for full efficiency)" << endl
    << "    Example: 1, 50, 90, 95 : ";
  cin >> tempString;
  getline(cin,tempString2);
  tempString+=tempString2;
  thresholdProbabilities_ = parseDoubleList(tempString);
  for (vector<double>::iterator it = thresholdProbabilities_.begin();
       it!=thresholdProbabilities_.end(); ++it) (*it)/=100;

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
    parameter = ctrim(tokens.at(0), " \"\n\t");
    value = ctrim(tokens.at(1), " \"\n\t");
    return true;
  }
}

vector<double> mainConfigHandler::parseDoubleList(string inString) {
  return split<double>(inString, ",");
}

bool mainConfigHandler::readConfigurationFile(ifstream& configFile) {
  char myLine[1024];
  string parameter, value;
  //bool styleFound=false;
  bool binFound=false;
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

  //return (styleFound&&layoutFound&&xmlFound&&momentaFound);
  return (binFound&&layoutFound&&xmlFound&&momentaFound&&triggerMomentaFound&&thresholdProbabilitiesFound);
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

string mainConfigHandler::getBinDirectory() {
  getConfiguration();
  return getBinDirectory_();
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

string mainConfigHandler::getDefaultMaterialsDirectory() {
  getConfiguration();
  return getDefaultMaterialsDirectory_();
}

string mainConfigHandler::getStandardIncludeDirectory() {
  getConfiguration();
  return getStandardIncludeDirectory_();
}

string mainConfigHandler::getGeometriesDirectory() {
  getConfiguration();
  return getGeometriesDirectory_();
}

string mainConfigHandler::getBinDirectory_() { return binDirectory_; }
string mainConfigHandler::getLayoutDirectory_() { return layoutDirectory_; }
string mainConfigHandler::getStandardDirectory_() { return standardDirectory_; }
string mainConfigHandler::getStyleDirectory_() { return layoutDirectory_+"/"+insur::default_styledir; }
string mainConfigHandler::getXmlDirectory_() { return standardDirectory_+"/"+insur::default_xmlpath; } 
string mainConfigHandler::getMattabDirectory_() { return standardDirectory_+"/"+insur::default_mattabdir; }
string mainConfigHandler::getIrradiationDirectory_() { return standardDirectory_+"/"+insur::default_irradiationdir; }
string mainConfigHandler::getDefaultMaterialsDirectory_() { return standardDirectory_+"/"+insur::default_materialsdir; }
string mainConfigHandler::getStandardIncludeDirectory_() { return standardDirectory_+"/"+insur::default_configdir+"/"+insur::default_stdincludedir; }
string mainConfigHandler::getGeometriesDirectory_() { return standardDirectory_+"/"+insur::default_geometriesdir; }


std::set<string> mainConfigHandler::preprocessConfiguration(istream& is, ostream& os, const string& istreamid) {
  using namespace std;
  string line;
  int numLine = 1;
  std::set<string> includeSet;
  includeSet.insert(istreamid);
  while(getline(is, line).good()) {
    if (line.find("//") != string::npos) line = line.erase(line.find("//"));
    string trimmed = trim(line);
    int tstart;
    if ((tstart = trimmed.find("@include")) != string::npos) {
      trimmed = trimmed.substr(tstart);
      int qstart, qend;
      string filename;
      if ((qstart = trimmed.find_first_of("\"")) != string::npos && (qend = trimmed.find_last_of("\"")) != string::npos) {
        filename = ctrim(trimmed.substr(qstart, qend - qstart + 1), "\"");
      } else {
        auto tokens = split(trimmed, " ");
        filename = tokens.size() > 1 ? tokens[1] : "";
      }
      string prefix = (trimmed.find("@includestd") != string::npos ? getStandardIncludeDirectory()+"/" : std::string(""));
      filename = prefix + filename;
      ifstream ifs(filename);
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

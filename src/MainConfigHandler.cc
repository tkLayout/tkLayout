/**
 * @file MainConfigHandler.cc
 * @brief implementation of class that take care of creating and reading the configuration file .tkgeometryrc
 */
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

#include <MainConfigHandler.hh>
#include "Units.hh"

using namespace std;
using namespace boost;

template <class T> bool from_string(T& t, const std::string& s, 
                                    std::ios_base& (*f)(std::ios_base&)) {
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

mainConfigHandler::mainConfigHandler() {
  goodConfigurationRead_ = false;
  binDirectory_ = "";
  layoutDirectory_ = "";
  standardDirectory_ = "";
}

mainConfigHandler& mainConfigHandler::instance() {
  static mainConfigHandler instance_;
  return instance_;
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

void mainConfigHandler::askBinDirectory() {
  cout << "*** What is the bin directory where you want to" << endl
      << "    place your executables?" << endl
      << "    Example: " << getenv(HOMEDIRECTORY) << "/bin : ";
  cin >> binDirectory_;
}

void mainConfigHandler::askLayoutDirectory() {
  cout << "*** What is the web server directory where you want to" << endl
    << "    place your output?" << endl
    << "    Example: " << getenv(HOMEDIRECTORY) << "/www/layouts : ";
  cin >> layoutDirectory_;
}

void mainConfigHandler::askStandardDirectory() {
  cout << "*** What is the standard output directory?" << endl
      << "    xml files and other various output will be put here" << endl
      << "    Example: " << getenv(HOMEDIRECTORY) << "/tkgeometry : ";
  cin >> standardDirectory_;
}

void mainConfigHandler::askMomenta() {
  string tempString = "";
  string tempString2 = "";

  cout << "*** Specify the list of transverse momenta to be used for the" << endl
      << "    tracking performance test (in GeV/c)" << endl
      << "    Example: 1, 10, 100 : ";
  cin >> tempString;

  getline(cin,tempString2);
  tempString+=tempString2;
  momenta_ = parseDoubleList(tempString);
  for (double& iMom : momenta_) iMom *= Units::GeV;
}

void mainConfigHandler::askTriggerMomenta() {
  string tempString = "";
  string tempString2 = "";

  cout << "*** Specify the list of transverse momenta to be used for the" << endl
      << "    trigger efficiency performance test (in GeV/c)" << endl
      << "    Example: 1, 2, 5, 10 : ";
  cin >> tempString;
  getline(cin,tempString2);
  tempString+=tempString2;
  triggerMomenta_ = parseDoubleList(tempString);
  for (double& iMom : triggerMomenta_) iMom *= Units::GeV;
}

void mainConfigHandler::askThresholdProbabilities() {
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


bool mainConfigHandler::createConfigurationFileFromQuestions(string& configFileName) {

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
      configFile << std::fixed << std::setprecision(2) << (*it)/Units::GeV;
    }
    configFile << "\"" << std::endl;

    configFile << TRIGGERMOMENTADEFINITION << "=\"";
    for (std::vector<double>::iterator it = triggerMomenta_.begin(); it!=triggerMomenta_.end(); ++it) {
      if (it!=triggerMomenta_.begin()) configFile << ", ";
      configFile << std::fixed << std::setprecision(2) << (*it)/Units::GeV;
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
    cerr << "ERROR: Cannot understand line: '" << codeLine << "' in the configuration file " << CONFIGURATIONFILENAME << endl;
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

bool mainConfigHandler::readConfigurationFile(string& configFileName) {
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
      // std::transform(parameter.begin(), parameter.end(), parameter.begin(), ::tolower);
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
        for (double& iMom : momenta_) iMom *= Units::GeV;
        momentaFound = true;
      } else if (parameter==TRIGGERMOMENTADEFINITION) {
        triggerMomenta_ = parseDoubleList(value);
        for (double& iMom : triggerMomenta_) iMom *= Units::GeV;
        triggerMomentaFound = true;
      } else if (parameter==THRESHOLDPROBABILITIESDEFINITION) {
        thresholdProbabilities_ = parseDoubleList(value);
        thresholdProbabilitiesFound = true;
      } else {
        cerr << "ERROR: Unknown parameter " << parameter << " in the configuration file " << CONFIGURATIONFILENAME << endl;
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

  return (binFound&&layoutFound&&xmlFound&&momentaFound&&triggerMomentaFound&&thresholdProbabilitiesFound);
}

bool mainConfigHandler::getConfiguration(bool checkDirExists /* = true */) {
  return readConfiguration(checkDirExists);
}

bool mainConfigHandler::getConfiguration(string& layoutDirectory) {
  bool result = readConfiguration(true);
  if (result) {
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

string mainConfigHandler::getDetIdSchemesDirectory() {
  getConfiguration();
  return getDetIdSchemesDirectory_();
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
string mainConfigHandler::getDetIdSchemesDirectory_() { return standardDirectory_+"/"+insur::default_detidschemesdir; }
string mainConfigHandler::getStandardIncludeDirectory_() { return standardDirectory_+"/"+insur::default_configdir+"/"+insur::default_stdincludedir; }
string mainConfigHandler::getGeometriesDirectory_() { return standardDirectory_+"/"+insur::default_geometriesdir; }

string ConfigInputOutput::getIncludedFile(string fileName) {
  vector<string> result;
  string testPath;
  for (const auto& aPath : includePathList) {
    testPath = aPath+"/"+fileName;
    testPath = boost::filesystem::system_complete(testPath).string();
    if ( boost::filesystem::exists(testPath) ) {
      result.push_back(testPath);
    }
  }
  if (result.size()==1) return result.at(0);
  else if (result.size()==0) {
    cerr << "ERROR: cannot find file " << fileName << " in any of the include path directories: " << endl;
    for (const auto& aPath : includePathList)  cerr << "   * " << aPath << endl;
    return "";
  } else {
    cerr << "ERROR: file " << fileName << " found multiple times in the include path directories: " << endl;
    for (const auto& aFile : result)  cerr << "   * " << aFile << endl;
    return "";
  }
}

std::set<string> mainConfigHandler::preprocessConfiguration(ConfigInputOutput cfgInOut) {
  using namespace std;

  istream& is = cfgInOut.is;
  ostream& os = cfgInOut.os;
  const string& absoluteFileName = cfgInOut.absoluteFileName;
  const string& relativeFileName = cfgInOut.relativeFileName;
  bool& standardInclude = cfgInOut.standardInclude;
  string absoluteFileNameDirectory = boost::filesystem::canonical(absoluteFileName).parent_path().string();

  // For the first run and if it's not a standardInclude I have to att
  // ths file's path to the list of visited include paths
  if (cfgInOut.includePathList.size()==0) {
    string firstFileName = boost::filesystem::canonical(cfgInOut.absoluteFileName).string();
    addNodeUrl(firstFileName, "file://"+cfgInOut.absoluteFileName);
    setNodeLocal(firstFileName, true);
  }
  if ((cfgInOut.includePathList.size()==0)||(!standardInclude)) cfgInOut.includePathList.insert(absoluteFileNameDirectory);

  // Assing the file id to this file name
  int thisFileId = getFileId(absoluteFileName);
  setNodeLocal(absoluteFileName, ! cfgInOut.standardInclude);
  prepareNodeOutput(absoluteFileName, relativeFileName, cfgInOut.webOutput);

  // Avoid double-counting: files included from this one should be
  // counted only once
  clearGraphLinks(thisFileId);
  std::string full_path = boost::filesystem::system_complete(absoluteFileName).string();

  string line;
  int numLine = 1;
  std::set<string> includeSet;
  includeSet.insert(absoluteFileName);

  while(getline(is, line).good()) {
    if (line.find("//") != string::npos) line = line.erase(line.find("//"));
    string trimmed = trim(line);
    int includeStart;
    // Merging a spec file (adding the latest includepath to the filename)
    if ((includeStart = trimmed.find("tiltedLayerSpecFile")) != string::npos) { 
      string rightPart = trimmed.substr(includeStart + strlen("tiltedLayerSpecFile"));
      string leftPart = trimmed.substr(0, includeStart + strlen("tiltedLayerSpecFile"));
      int lastSpace;
      for (lastSpace=0; (rightPart[lastSpace]==' ')&&(lastSpace<rightPart.size()); lastSpace++);
      rightPart = rightPart.substr(lastSpace);
      line = leftPart + " " + absoluteFileNameDirectory + "/" + rightPart;
    }

    if ((includeStart = trimmed.find("@include")) != string::npos) { //@include @include-std
      trimmed = trimmed.substr(includeStart);
      int quoteStart, quoteEnd;
      string nextIncludeFileName;
      if ((quoteStart = trimmed.find_first_of("\"")) != string::npos && (quoteEnd = trimmed.find_last_of("\"")) != string::npos) {
        nextIncludeFileName = ctrim(trimmed.substr(quoteStart, quoteEnd - quoteStart + 1), "\"");
      } else {
        auto tokens = split(trimmed, " ");
        nextIncludeFileName = tokens.size() > 1 ? tokens[1] : "";
      }

      bool includeStdOld = trimmed.find("@includestd") != string::npos;  // both @includestd (deprecated) and @include-std (preferred) are supported 
      bool includeStdNew = trimmed.find("@include-std") != string::npos;

      string fullIncludedFileName;
      int includedFileId;
      if (includeStdOld || includeStdNew) {
	fullIncludedFileName = getStandardIncludeDirectory()+ "/" + nextIncludeFileName;
      } else {
	fullIncludedFileName = cfgInOut.getIncludedFile(nextIncludeFileName);
      }
      includedFileId = getFileId(fullIncludedFileName);
      
      ifstream ifs(fullIncludedFileName);

      if (ifs) {
        stringstream ss;
	ConfigInputOutput nextIncludeInputOutput(ifs, ss);
	nextIncludeInputOutput.includePathList=cfgInOut.includePathList;
	nextIncludeInputOutput.standardInclude=(includeStdOld || includeStdNew);
	nextIncludeInputOutput.absoluteFileName=fullIncludedFileName;
	nextIncludeInputOutput.relativeFileName=nextIncludeFileName;
	nextIncludeInputOutput.webOutput=cfgInOut.webOutput;
        auto&& moreIncludes = preprocessConfiguration(nextIncludeInputOutput);

	// Graph node links
	addGraphLink(thisFileId, includedFileId);
        includeSet.insert(moreIncludes.begin(), moreIncludes.end());
        string indent = line.substr(0, line.find_first_not_of(" \t"));

        while (getline(ss, line).good()) {
          os << indent << line << endl;   
        }
      } else {
        cerr << "ERROR: ignoring " << ( (includeStdOld||includeStdNew) ? "@include-std" : "@include" ) << " directive in " << absoluteFileName << ":" << numLine << " : could not open included file : " << nextIncludeFileName << endl;
      }
    } else {
      os << line << endl;
    }
    numLine++;
  }
  return includeSet;
}

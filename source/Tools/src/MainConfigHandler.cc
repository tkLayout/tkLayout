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

#include <sys/types.h>

#include "GraphVizCreator.h"
#include "MessageLogger.h"
#include "string_functions.h"
#include "Units.h"

using namespace std;
using namespace boost;

template <class T> bool from_string(T& t, const std::string& s, 
                                    std::ios_base& (*f)(std::ios_base&)) {
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

//
// Singleton constructor method -> must be private (singleton pattern)
//
MainConfigHandler::MainConfigHandler() {
  m_goodConfigurationRead = false;
  //m_styleDirectory = "";
  m_binDirectory      = "";
  m_layoutDirectory   = "";
  m_standardDirectory = "";
  m_projectName       = "CMS";
  m_resultsAuthor     = "Unknown author";
}

//
// MainConfigHandler access method -> get instance of singleton class MainConfigHandler
//
MainConfigHandler& MainConfigHandler::getInstance() {

  static MainConfigHandler s_instance;
  return s_instance;
}

//
// MainConfigHandler method to access GraphVizCreator member method
//
GraphVizCreator& MainConfigHandler::getGraphVizCreator() {

  static GraphVizCreator s_graphCreator;
  return s_graphCreator;
}

// Helper method to preprocess any input configuration file (is) (its full address specified by istreamid) ->
// recursively pass through all included files specified by @include & @includestd/@include-std pragma and
// include its content to one big configuration file (defined as os)
std::set<string> MainConfigHandler::preprocessConfiguration(istream& is, ostream& os, const string& istreamid, bool standardIncl) {

  using namespace std;

  // Create instance of configuration input/output helper class
  ConfigInputOutput cfgInOut(is, os, istreamid);

  bool& standardInclude = standardIncl;

  const string& absFileName = cfgInOut.absFileName;
  string absFileNameDirectory = boost::filesystem::canonical(absFileName).parent_path().string();

  // For the first run and if it's not a standardInclude I have to att
  // ths file's path to the list of visited include paths
  if (cfgInOut.includePathList.size()==0) {
    string firstFileName = boost::filesystem::canonical(cfgInOut.absFileName).string();
    getGraphVizCreator().addNodeUrl(firstFileName, "file://"+cfgInOut.absFileName);
    getGraphVizCreator().setNodeLocal(firstFileName, true);
  }
  if ((cfgInOut.includePathList.size()==0)||(!standardInclude)) cfgInOut.includePathList.insert(absFileNameDirectory);

  // Assign the file id to this file name
  int thisFileId = getGraphVizCreator().getFileId(cfgInOut.absFileName);
  getGraphVizCreator().setNodeLocal(cfgInOut.absFileName, ! cfgInOut.standardInclude);
  getGraphVizCreator().prepareNodeOutput(cfgInOut.absFileName, cfgInOut.relFileName, cfgInOut.webOutput);

  // Avoid double-counting: files included from this one should be
  // counted only once
  getGraphVizCreator().clearGraphLinks(thisFileId);
  std::string full_path = boost::filesystem::system_complete(cfgInOut.absFileName).string();

  string line;
  int numLine = 1;
  std::set<string> includeSet;
  includeSet.insert(cfgInOut.absFileName);

  while(getline(is, line).good()) {

    if (line.find("//") != string::npos) line = line.erase(line.find("//"));
    string trimmed = trim(line);
    int includeStart;

    // Merging a spec file (adding the latest includepath to the filename)
    if ((includeStart  = trimmed.find("tiltedLayerSpecFile")) != string::npos) {
      string rightPart = trimmed.substr(includeStart + strlen("tiltedLayerSpecFile"));
      string leftPart  = trimmed.substr(0, includeStart + strlen("tiltedLayerSpecFile"));
      int lastSpace;
      for (lastSpace=0; (rightPart[lastSpace]==' ')&&(lastSpace<rightPart.size()); lastSpace++);
      rightPart = rightPart.substr(lastSpace);
      line = leftPart + " " + absFileNameDirectory + "/" + rightPart;
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

      bool includeStdOld = trimmed.find("@includestd")  != string::npos;  // both @includestd (deprecated) and @include-std (preferred) are supported
      bool includeStdNew = trimmed.find("@include-std") != string::npos;

      string fullIncludedFileName;
      int includedFileId;
      if (includeStdOld || includeStdNew) {
        fullIncludedFileName = getStandardIncludeDirectory()+ "/" + nextIncludeFileName;
      } else {
        fullIncludedFileName = cfgInOut.getIncludedFile(nextIncludeFileName);
      }
      includedFileId = getGraphVizCreator().getFileId(fullIncludedFileName);

      ifstream ifs(fullIncludedFileName);
      if (ifs) {

        stringstream ss;
        auto&& moreIncludes = preprocessConfiguration(ifs, ss, fullIncludedFileName, includeStdOld || includeStdNew);

        // Graph node links
        getGraphVizCreator().addGraphLink(thisFileId, includedFileId);
        includeSet.insert(moreIncludes.begin(), moreIncludes.end());
        string indent = line.substr(0, line.find_first_not_of(" \t"));

        while (getline(ss, line).good()) {
          os << indent << line << endl;
        }
      } else {
        std::ostringstream message;
        message << "ERROR: ignoring "
                << ( (includeStdOld||includeStdNew) ? "@include-std" : "@include" )
                << " directive in "
                << absFileName << ":"
                << numLine
                << " : could not open included file : "
                << nextIncludeFileName << endl;
        logERROR(message);
      }
    } else {
      os << line << endl;
    }
    numLine++;
  }
  return includeSet;
}

//
// Check that given directory exists on the filesystem
//
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

//
// Ask user on standard input for bin directory, where to place executables
//
void MainConfigHandler::askBinDirectory() {
  cout << "*** What is the bin directory where you want to" << endl
       << "    place your executables?" << endl
       << "    Example: " << getenv(c_HOMEDIRECTORY) << "/bin : ";
  cin  >> m_binDirectory;
}

//
// Ask user on standard input for layout directory, where to place the www output
//
void MainConfigHandler::askLayoutDirectory() {
  cout << "*** What is the web server directory where you want to" << endl
       << "    place your output?" << endl
       << "    Example: " << getenv(c_HOMEDIRECTORY) << "/www/layouts : ";
  cin  >> m_layoutDirectory;
}

//
// Ask user on standard input for standard directory. xml files and other various output will be put here.
//
void MainConfigHandler::askStandardDirectory() {
  cout << "*** What is the standard output directory?" << endl
       << "    xml files and other various output will be put here" << endl
       << "    Example: " << getenv(c_HOMEDIRECTORY) << "/tkgeometry : ";
  cin  >> m_standardDirectory;
}

//
// Ask user on standard input for particle momenta to be studied
//
void MainConfigHandler::askMomenta() {
  string tempString = "";
  string tempString2 = "";

  cout << "*** Specify the list of transverse momenta to be used for the" << endl
       << "    tracking performance test (in GeV/c)" << endl
       << "    Example: 1, 10, 100 : ";
  cin  >> tempString;

  getline(cin,tempString2);
  tempString+=tempString2;
  m_momenta = parseDoubleList(tempString);
  for (double& iMom : m_momenta) iMom *= Units::GeV;
}

//
// Ask user on standard input for particle momenta to be studied by trigger
//
void MainConfigHandler::askTriggerMomenta() {
  string tempString = "";
  string tempString2 = "";

  cout << "*** Specify the list of transverse momenta to be used for the" << endl
       << "    trigger efficiency performance test (in GeV/c)" << endl
       << "    Example: 1, 2, 5, 10: ";
  cin  >> tempString;
  getline(cin,tempString2);
  tempString+=tempString2;
  m_triggerMomenta = parseDoubleList(tempString);
  for (double& iMom : m_triggerMomenta) iMom *= Units::GeV;
}

//
// Specify the list of trigger efficiency to be used for pt threshold find scan
//
void MainConfigHandler::askThresholdProbabilities() {
  string tempString = "";
  string tempString2 = "";

  cout << "*** Specify the list of trigger efficiency to be used for the" << endl
      << "    pt threshold find test (in percent: write 100 for full efficiency)" << endl
      << "    Example: 1, 50, 90, 95 : ";
  cin >> tempString;
  getline(cin,tempString2);
  tempString+=tempString2;
  m_thresholdProbabilities = parseDoubleList(tempString);
  for (double& threshold : m_thresholdProbabilities) threshold /=100;
}

//
// Create base configuration file for tkLayout .tkgeometryrc
//
bool MainConfigHandler::createConfigurationFileFromQuestions(string& configFileName) {

  // Clear screen
  //cout << "\033[2J"; // Clears the screen
  //cout << "\033[1;1H"; // Places cursor on line 1

  // I have no configuration, so I must create it
  cout << "Could not find the configuration file "  << configFileName << " maybe this is the first time you run with the new system." << endl;
  cout << "Answer to the following questions to have your configuration file automatically created." << endl;
  cout << "You will be later able to edit it manually, or you can just delete it and answer these questions again." << endl;
  cout << endl;

  askBinDirectory();
  if (!checkDirectory(m_binDirectory)) return false;
  cout << endl;

  askLayoutDirectory();
  if (!checkDirectory(m_layoutDirectory)) return false;
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
    configFile << c_BINDIRECTORYDEFINITION      << "=\"" << m_binDirectory      << "\"" << endl;
    configFile << c_LAYOUTDIRECTORYDEFINITION   << "=\"" << m_layoutDirectory   << "\"" << endl;
    configFile << c_STANDARDDIRECTORYDEFINITION << "=\"" << m_standardDirectory << "\"" << endl;

    configFile << c_MOMENTADEFINITION << "=\"";
    for (auto it = m_momenta.begin(); it!=m_momenta.end(); ++it) {
      if (it!=m_momenta.begin()) configFile << ", ";
      configFile << std::fixed << std::setprecision(2) << (*it)/Units::GeV;
    }
    configFile << "\"" << std::endl;

    configFile << c_TRIGGERMOMENTADEFINITION << "=\"";
    for (auto it = m_triggerMomenta.begin(); it!=m_triggerMomenta.end(); ++it) {
      if (it!=m_triggerMomenta.begin()) configFile << ", ";
      configFile << std::fixed << std::setprecision(2) << (*it)/Units::GeV;
    }
    configFile << "\"" << std::endl;

    configFile << c_THRESHOLDPROBABILITIESDEFINITION << "=\"";
    for (auto it = m_thresholdProbabilities.begin(); it!=m_thresholdProbabilities.end(); ++it) {
      if (it!=m_thresholdProbabilities.begin()) configFile << ", ";
      configFile << std::fixed << std::setprecision(2) << (*it)*100;
    }
    configFile << "\"" << std::endl;

    configFile.close();
  }

  return true;
}

//
// Helper method to read line
//
bool MainConfigHandler::parseLine(const char* codeLine, string& parameter, string& value) {

  std::vector<string> tokens = split(codeLine, "=");

  if (tokens.empty()) {
    parameter = "";
    value = "";
    return false;
  } else if (tokens.size() < 2) { 
    cerr << "Cannot understand line: '" << codeLine << "' in the configuration file " << c_CONFIGURATIONFILENAME << endl;
    return false;
  } else {
    parameter = ctrim(tokens.at(0), " \"\n\t");
    value = ctrim(tokens.at(1), " \"\n\t");
    return true;
  }
}

//
// Helper method to convert string list of doubles to std::vector<double>
//
vector<double> MainConfigHandler::parseDoubleList(string inString) {
  return split<double>(inString, ",");
}

//
// Read configuration from c_CONFIGURATIONFILENAME or ask user to define configuration on the fly if doesn't exist
//
bool MainConfigHandler::readConfiguration( bool checkDirExists ) {
  if (m_goodConfigurationRead) return true;

  ifstream configFile;
  string homeDirectory  = string(getenv(c_HOMEDIRECTORY));
  string configFileName = getConfigFileName();
  bool goodConfig=false;

  configFile.open(configFileName.c_str(), ifstream::in);
  if (!configFile.good()) {
    configFile.close();
    goodConfig = createConfigurationFileFromQuestions(configFileName);
  }
  else {
    configFile.close();
    // I will read the configuration out of that
    goodConfig = readConfigurationFile(configFileName);
    if (goodConfig) {
      if (checkDirExists) {
        // Check the basic configuration directories
        if (!checkDirectory(m_layoutDirectory)) {
          cout << "You probably need to edit or delete the configuration file " << c_CONFIGURATIONFILENAME << endl;
          return false;
        }
        if (!checkDirectory(m_standardDirectory)) {
          cout << "You probably need to edit or delete the configuration file " << c_CONFIGURATIONFILENAME << endl;
          return false;
        }
        // Check the mandatory subdirectories in the main
        if (!checkDirectory(getXmlDirectory_P())) return false;
        if (!checkDirectory(getMattabDirectory_P())) return false;
      }
    }
    else { // not good config read
      cout << "Configuration file '" << configFileName << "' not properly formatted. You probably need to edit or delete it" << endl;
      return false;
    }
  }

  m_goodConfigurationRead = goodConfig;
  return goodConfig;
}

// Helper method called by readConfiguration()
//
bool MainConfigHandler::readConfigurationFile(string& configFileName) {

  char     myLine[1024];
  string   parameter, value;
  //bool     styleFound=false;
  bool     binFound=false;
  bool     layoutFound=false;
  bool     xmlFound=false;
  bool     momentaFound=false;
  bool     triggerMomentaFound=false;
  bool     thresholdProbabilitiesFound=false;
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
      if (parameter==c_BINDIRECTORYDEFINITION) {
        m_binDirectory = value;
        binFound = true;
      }
      else if (parameter==c_LAYOUTDIRECTORYDEFINITION) {
        m_layoutDirectory = value;
        layoutFound = true;
      }
      else if (parameter==c_STANDARDDIRECTORYDEFINITION) {
        m_standardDirectory = value;
        xmlFound = true;
      }
      else if (parameter==c_MOMENTADEFINITION) {
        m_momenta = parseDoubleList(value);
        for (double& iMom : m_momenta) iMom *= Units::GeV;
        momentaFound = true;
      }
      else if (parameter==c_TRIGGERMOMENTADEFINITION) {
        m_triggerMomenta = parseDoubleList(value);
        for (double& iMom : m_triggerMomenta) iMom *= Units::GeV;
        triggerMomentaFound = true;
      }
      else if (parameter==c_THRESHOLDPROBABILITIESDEFINITION) {
        m_thresholdProbabilities = parseDoubleList(value);
        thresholdProbabilitiesFound = true;
      }
      else if (parameter==c_PROJECTNAME) {
        m_projectName = value;
      }
      else if (parameter==c_RESULTSAUTHOR) {
        m_resultsAuthor = value;
      }
      else {
        cerr << "Unknown parameter " << parameter << " in the configuration file " << c_CONFIGURATIONFILENAME << endl;
      }
    }
  }
  configFileIs.close();

  //ask here for bin directory for backward compatibility of the install script
  if(!binFound&&layoutFound&&xmlFound&&momentaFound&&triggerMomentaFound&&thresholdProbabilitiesFound){

    askBinDirectory();
    if (checkDirectory(m_binDirectory)) {
      configFileOs.open(configFileName.c_str(), ifstream::app);
      if (!configFileOs.good()) {
        cout << "Could not open " << configFileName << " for writing. I quit."  << endl;
        configFileOs.close();
        return false;
      }
      else {
        configFileOs << c_BINDIRECTORYDEFINITION << "=\"" << m_binDirectory << "\"" << endl;

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
    layoutDirectory = m_layoutDirectory;
  }
  return result;
}

string MainConfigHandler::getConfigFileName() {
  char* specialConfigFile = getenv(c_CONFIGURATIONFILENAMEDEFINITION);
  if (specialConfigFile) return string(specialConfigFile);
  string homeDirectory = string(getenv(c_HOMEDIRECTORY));
  return homeDirectory+"/"+c_CONFIGURATIONFILENAME;
}



vector<double>& MainConfigHandler::getMomenta() {
  getConfiguration();
  return m_momenta;
}

vector<double>& MainConfigHandler::getTriggerMomenta() {
  getConfiguration();
  return m_triggerMomenta;
}

vector<double>& MainConfigHandler::getThresholdProbabilities() {
  getConfiguration();
  return m_thresholdProbabilities;
}

string MainConfigHandler::getBinDirectory() {
  getConfiguration();
  return getBinDirectory_P();
}

string MainConfigHandler::getLayoutDirectory() {
  getConfiguration();
  return getLayoutDirectory_P();
}

string MainConfigHandler::getStandardDirectory() {
  getConfiguration();
  return getStandardDirectory_P();
}

string MainConfigHandler::getStyleDirectory() {
  getConfiguration();
  return getStyleDirectory_P();
}

string MainConfigHandler::getXmlDirectory() {
  getConfiguration();
  return getXmlDirectory_P();
}

string MainConfigHandler::getMattabDirectory() {
  getConfiguration();
  return getMattabDirectory_P();
}

string MainConfigHandler::getIrradiationDirectory() {
  getConfiguration();
  return getIrradiationDirectory_P();
}

string MainConfigHandler::getDefaultMaterialsDirectory() {
  getConfiguration();
  return getDefaultMaterialsDirectory_P();
}

string MainConfigHandler::getStandardIncludeDirectory() {
  getConfiguration();
  return getStandardIncludeDirectory_P();
}

string MainConfigHandler::getGeometriesDirectory() {
  getConfiguration();
  return getGeometriesDirectory_P();
}

std::string MainConfigHandler::getProjectName() {
  getConfiguration();
  return m_projectName;
}
std::string MainConfigHandler::getResultsAuthor() {
  getConfiguration();
  return m_resultsAuthor;
}

string MainConfigHandler::getBinDirectory_P()              { return m_binDirectory; }
string MainConfigHandler::getLayoutDirectory_P()           { return m_layoutDirectory; }
string MainConfigHandler::getStandardDirectory_P()         { return m_standardDirectory; }
string MainConfigHandler::getStyleDirectory_P()            { return m_layoutDirectory  +"/"+default_styledir; }
string MainConfigHandler::getXmlDirectory_P()              { return m_standardDirectory+"/"+default_xmlpath; }
string MainConfigHandler::getMattabDirectory_P()           { return m_standardDirectory+"/"+default_mattabdir; }
string MainConfigHandler::getIrradiationDirectory_P()      { return m_standardDirectory+"/"+default_irradiationdir; }
string MainConfigHandler::getDefaultMaterialsDirectory_P() { return m_standardDirectory+"/"+default_materialsdir; }
string MainConfigHandler::getStandardIncludeDirectory_P()  { return m_standardDirectory+"/"+default_configdir+"/"+default_stdincludedir; }
string MainConfigHandler::getGeometriesDirectory_P()       { return m_standardDirectory+"/"+default_geometriesdir; }

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
    std::ostringstream message;
    message << "ConfigInputOutput::getIncludedFile: Cannot find file " << fileName << " in any of the include path directories: " << endl;
    for (const auto& aPath : includePathList)  message << "   * " << aPath << endl;
    logERROR(message);
    return "";
  }
  else {
    std::ostringstream message;
    message << "ConfigInputOutput::getIncludedFile: File " << fileName << " found multiple times in the include path directories: " << endl;
    for (const auto& aFile : result)  message << "   * " << aFile << endl;
    logERROR(message);
    return "";
  }
}


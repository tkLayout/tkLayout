/**
 * @file mainConfigHandler.h
 * @brief Takes care of creating, reading and administering the configuration file .tkgeometryrc
 */

#ifndef _MAINCONFIGHANDLER_H__
#define _MAINCONFIGHANDLER_H__

#include <string>
#include <vector>
#include <set>
#include <global_constants.hh>
#include <global_funcs.hh>

#include "GraphVizCreator.hh"

using namespace std;

#define HOMEDIRECTORY "HOME"
#define CONFIGURATIONFILENAME ".tkgeometryrc"
#define CONFIGURATIONFILENAMEDEFINITION "TKGEOMETRYRC"
#define BINDIRECTORYDEFINITION "TKG_BINDIRECTORY"
#define LAYOUTDIRECTORYDEFINITION "TKG_LAYOUTDIRECTORY" 
#define STANDARDDIRECTORYDEFINITION "TKG_STANDARDDIRECTORY" 
#define MOMENTADEFINITION "TKG_MOMENTA" 
#define TRIGGERMOMENTADEFINITION "TKG_TRIGGERMOMENTA" 
#define THRESHOLDPROBABILITIESDEFINITION "TKG_THRESHOLD_PROB"

class ConfigInputOutput {
public:
  ConfigInputOutput(istream& newIs, ostream& newOs) : is(newIs) , os(newOs) {}
  istream& is;
  ostream& os;
  string absoluteFileName = "";
  string relativeFileName = "";
  bool standardInclude = false;
  set<string> includePathList;
  string getIncludedFile(string fileName);
  bool webOutput;
};

// This object wil read the configuration only once
// If the configuration file is not present, wuations will
// be asked directly through std::cin and the corresponding
// configuration file will be saved in the home directory

class mainConfigHandler : public GraphVizCreator {
private:
  mainConfigHandler();
public:
  ~mainConfigHandler() {};
  static mainConfigHandler& instance();
  bool getConfiguration(bool checkDirExists = true);
  //bool getConfiguration(string& layoutDirectory, string& xmlDirectory);
  bool getConfiguration(string& layoutDirectory);
  string getBinDirectory();
  string getLayoutDirectory();
  string getStandardDirectory();
  string getStyleDirectory();
  string getXmlDirectory();
  string getMattabDirectory();
  string getIrradiationDirectory();
  string getDefaultMaterialsDirectory();
  string getDetIdSchemesDirectory();
  string getStandardIncludeDirectory();
  string getGeometriesDirectory();
  string getConfigFileName();
  std::set<string> preprocessConfiguration(ConfigInputOutput);
  vector<double>& getMomenta();
  vector<double>& getTriggerMomenta();
  vector<double>& getThresholdProbabilities();
private:
  bool goodConfigurationRead_;
  //string styleDirectory_;
  string binDirectory_;
  string layoutDirectory_;
  //string xmlDirectory_;
  string standardDirectory_;
  vector<double> momenta_;
  vector<double> triggerMomenta_;
  vector<double> thresholdProbabilities_;
  bool checkDirectory(string dirName) ;
  void askBinDirectory();
  void askLayoutDirectory();
  void askStandardDirectory();
  void askMomenta();
  void askTriggerMomenta();
  void askThresholdProbabilities();
  bool createConfigurationFileFromQuestions(string& configFileName);
  bool parseLine(const char* codeLine, string& parameter, string& value);
  bool readConfigurationFile(string& configFileName);
  bool readConfiguration(bool checkDirExists);
  vector<double> parseDoubleList(string);
  string getBinDirectory_();
  string getLayoutDirectory_();
  string getStandardDirectory_();
  string getStyleDirectory_();
  string getXmlDirectory_();
  string getMattabDirectory_();
  string getIrradiationDirectory_();
  string getDefaultMaterialsDirectory_();
  string getDetIdSchemesDirectory_();
  string getStandardIncludeDirectory_();
  string getGeometriesDirectory_();

};

#endif

#ifndef _MAINCONFIGHANDLER_H__
#define _MAINCONFIGHANDLER_H__

#include <string>
#include <vector>
#include <global_constants.h>

using namespace std;

#define HOMEDIRECTORY "HOME"
#define CONFIGURATIONFILENAME ".tkgeometryrc"
#define CONFIGURATIONFILENAMEDEFINITION "TKGEOMETRYRC"
#define LAYOUTDIRECTORYDEFINITION "TKG_LAYOUTDIRECTORY" 
#define STANDARDDIRECTORYDEFINITION "TKG_STANDARDDIRECTORY" 
#define MOMENTADEFINITION "TKG_MOMENTA" 
#define TRIGGERMOMENTADEFINITION "TKG_TRIGGERMOMENTA" 
#define THRESHOLDPROBABILITIESDEFINITION "TKG_THRESHOLD_PROB"

// This object wil read the configuration only once
// If the configuration file is not present, wuations will
// be asked directly through std::cin and the corresponding
// configuration file will be saved in the home directory


class mainConfigHandler {
public:
  mainConfigHandler();
  ~mainConfigHandler() {};
  bool getConfiguration(bool checkDirExists = true);
  //bool getConfiguration(string& layoutDirectory, string& xmlDirectory);
  bool getConfiguration(string& layoutDirectory);
  string getLayoutDirectory();
  string getStandardDirectory();
  string getStyleDirectory();
  string getXmlDirectory();
  string getMattabDirectory();
  string getIrradiationDirectory();
  string getRootfileDirectory();
  string getGraphDirectory();
  string getSummaryDirectory();
  string getDefaultMaterialsDirectory();
  string getConfigFileName();
  vector<double>& getMomenta();
  vector<double>& getTriggerMomenta();
  vector<double>& getThresholdProbabilities();
private:
  bool goodConfigurationRead_;
  //string styleDirectory_;
  string layoutDirectory_;
  //string xmlDirectory_;
  string standardDirectory_;
  vector<double> momenta_;
  vector<double> triggerMomenta_;
  vector<double> thresholdProbabilities_;
  bool checkDirectory(string dirName) ;
  bool createConfigurationFileFromQuestions(string& configFileName);
  bool parseLine(const char* codeLine, string& parameter, string& value);
  bool readConfigurationFile(ifstream& configFile);
  bool readConfiguration(bool checkDirExists);
  vector<double> parseDoubleList(string);
  string getLayoutDirectory_();
  string getStandardDirectory_();
  string getStyleDirectory_();
  string getXmlDirectory_();
  string getMattabDirectory_();
  string getIrradiationDirectory_();
  string getRootfileDirectory_();
  string getGraphDirectory_();
  string getSummaryDirectory_();
  string getDefaultMaterialsDirectory_();
};

#endif

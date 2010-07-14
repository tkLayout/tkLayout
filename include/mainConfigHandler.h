#ifndef _MAINCONFIGHANDLER_H__
#define _MAINCONFIGHANDLER_H__

#include <string>
#include <vector>
#include <global_constants.h>

using namespace std;

#define HOMEDIRECTORY "HOME"
#define CONFIGURATIONFILENAME ".tkgeometryrc"
#define LAYOUTDIRECTORYDEFINITION "TKG_LAYOUTDIRECTORY" 
#define STANDARDDIRECTORYDEFINITION "TKG_STANDARDDIRECTORY" 
#define MOMENTADEFINITION "TKG_MOMENTA" 

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
  string getRootfileDirectory();
  string getGraphDirectory();
  string getSummaryDirectory();
  string getConfigFileName();
  vector<double> getMomenta();
 private:
  bool goodConfigurationRead_;
  //string styleDirectory_;
  string layoutDirectory_;
  //string xmlDirectory_;
  string standardDirectory_;
  vector<double> momenta_;
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
  string getRootfileDirectory_();
  string getGraphDirectory_();
  string getSummaryDirectory_();
};

#endif

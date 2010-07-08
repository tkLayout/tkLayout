#ifndef _MAINCONFIGHANDLER_H__
#define _MAINCONFIGHANDLER_H__

#include <string>
#include <vector>

using namespace std;

#define HOMEDIRECTORY "HOME"
#define CONFIGURATIONFILENAME ".tkgeometryrc"
#define STYLEDIRECTORYDEFINITION "styleDirectory"
#define STYLEDIRECTORYDEFINITIONLOWERCASE "styledirectory"
#define LAYOUTDIRECTORYDEFINITION "layoutDirectory" 
#define LAYOUTDIRECTORYDEFINITIONLOWERCASE "layoutdirectory" 
#define XMLDIRECTORYDEFINITION "xmlDirectory" 
#define XMLDIRECTORYDEFINITIONLOWERCASE "xmldirectory" 
#define MOMENTADEFINITION "momenta" 
#define MOMENTADEFINITIONLOWERCASE "momenta" 

// This object wil read the configuration only once
// If the configuration file is not present, wuations will
// be asked directly through std::cin and the corresponding
// configuration file will be saved in the home directory

class mainConfigHandler {
 public:
  mainConfigHandler();
  ~mainConfigHandler() {};
  bool getConfiguration();
  bool getConfiguration(string& styleDirectory, string& layoutDirectory, string& xmlDirectory);
  bool getConfiguration(string& styleDirectory, string& layoutDirectory);
  string getStyleDirectory() { getConfiguration(); return styleDirectory_;} ;
  string getLayoutDirectory() { getConfiguration(); return layoutDirectory_;} ;
  string getXmlDirectory() { getConfiguration(); return xmlDirectory_;} ;
  vector<double> getMomenta() { getConfiguration(); return momenta_; } ;
 private:
  bool goodConfigurationRead_;
  string styleDirectory_;
  string layoutDirectory_;
  string xmlDirectory_;
  vector<double> momenta_;
  bool checkDirectory(string dirName) ;
  bool createConfigurationFileFromQuestions(string& configFileName);
  bool parseLine(const char* codeLine, string& parameter, string& value);
  bool readConfigurationFile(ifstream& configFile);
  bool readConfiguration();
  vector<double> parseDoubleList(string);
};

#endif

#ifndef _MAINCONFIGHANDLER_H__
#define _MAINCONFIGHANDLER_H__

#include <string>

using namespace std;

#define HOMEDIRECTORY "HOME"
#define CONFIGURATIONFILENAME ".tkgeometryrc"
#define STYLEDIRECTORYDEFINITION "styleDirectory"
#define STYLEDIRECTORYDEFINITIONLOWERCASE "styledirectory"
#define LAYOUTDIRECTORYDEFINITION "layoutDirectory" 
#define LAYOUTDIRECTORYDEFINITIONLOWERCASE "layoutdirectory" 

class mainConfigHandler {
 public:
  mainConfigHandler() {} ;
  ~mainConfigHandler() {};
  bool getConfiguration(string& styleDirectory, string& layoutDirectory);
 private:
  bool checkDirectory(string dirName) ;
  bool createConfigurationFileFromQuestions(string& configFileName, string& styleDirectory, string& layoutDirectory);
  bool parseLine(const char* codeLine, string& parameter, string& value);
  bool readConfigurationFile(ifstream& configFile, string& styleDirectory, string& layoutDirectory);
};

#endif

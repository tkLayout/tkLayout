#ifndef _MAINCONFIGHANDLER_H__
#define _MAINCONFIGHANDLER_H__

#include <string>
#include <vector>
#include <set>

#define HOMEDIRECTORY "HOME"
#define CONFIGURATIONFILENAME ".tkgeometryrc"
#define CONFIGURATIONFILENAMEDEFINITION "TKGEOMETRYRC"
#define BINDIRECTORYDEFINITION "TKG_BINDIRECTORY"
#define LAYOUTDIRECTORYDEFINITION "TKG_LAYOUTDIRECTORY" 
#define STANDARDDIRECTORYDEFINITION "TKG_STANDARDDIRECTORY" 
#define MOMENTADEFINITION "TKG_MOMENTA" 
#define TRIGGERMOMENTADEFINITION "TKG_TRIGGERMOMENTA" 
#define THRESHOLDPROBABILITIESDEFINITION "TKG_THRESHOLD_PROB"

/*
 * @class MainConfigHandler
 * @brief Singleton class reading-in the tkLayout configuration (only once) from .tkgeometryrc file and
 * providing the information through out the software.
 * @details Singleton class reading-in the tkLayout configuration (only once) from .tkgeometryrc file and
 * providing the information through out the software. If the configuration file is missing, user will be
 * asked directly through std::cin on the fly. The corresponding configuration will be saved in a file in
 * the home directory.
 */
class MainConfigHandler {

 public:

  //! Destructor
  ~MainConfigHandler() {};

  //! MainConfigHandler access method -> get instance of singleton class MainConfigHandler
  static MainConfigHandler& getInstance();

  // Other methods
  bool getConfiguration(bool checkDirExists = true);
  //bool getConfiguration(std::string& layoutDirectory, std::string& xmlDirectory);
  bool getConfiguration(std::string& layoutDirectory);
  std::string getBinDirectory();
  std::string getLayoutDirectory();
  std::string getStandardDirectory();
  std::string getStyleDirectory();
  std::string getXmlDirectory();
  std::string getMattabDirectory();
  std::string getIrradiationDirectory();
  std::string getDefaultMaterialsDirectory();
  std::string getStandardIncludeDirectory();
  std::string getGeometriesDirectory();
  std::string getConfigFileName();
  std::set<std::string> preprocessConfiguration(std::istream& is, std::ostream& os, const std::string& istreamid);
  std::vector<double>& getMomenta();
  std::vector<double>& getTriggerMomenta();
  std::vector<double>& getThresholdProbabilities();

 private:

 //! Singleton constructor method
 MainConfigHandler();

 bool goodConfigurationRead_;
  //std::string styleDirectory_;
  std::string binDirectory_;
  std::string layoutDirectory_;
  //std::string xmlDirectory_;
  std::string standardDirectory_;
  std::vector<double> momenta_;
  std::vector<double> triggerMomenta_;
  std::vector<double> thresholdProbabilities_;
  bool checkDirectory(std::string dirName) ;
  void askBinDirectory();
  void askLayoutDirectory();
  void askStandardDirectory();
  void askMomenta();
  void askTriggerMomenta();
  void askThresholdProbabilities();
  bool createConfigurationFileFromQuestions(std::string& configFileName);
  bool parseLine(const char* codeLine, std::string& parameter, std::string& value);
  bool readConfigurationFile(std::string& configFileName);
  bool readConfiguration(bool checkDirExists);
  std::vector<double> parseDoubleList(std::string);
  std::string getBinDirectory_();
  std::string getLayoutDirectory_();
  std::string getStandardDirectory_();
  std::string getStyleDirectory_();
  std::string getXmlDirectory_();
  std::string getMattabDirectory_();
  std::string getIrradiationDirectory_();
  std::string getDefaultMaterialsDirectory_();
  std::string getStandardIncludeDirectory_();
  std::string getGeometriesDirectory_();
};

#endif

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
#include "MessageLogger.hh"
#include "GraphVizCreator.hh"

class ConfigInputOutput {
public:
  ConfigInputOutput(std::istream& newIs, std::ostream& newOs) : is(newIs) , os(newOs) {}
  std::istream& is;
  std::ostream& os;
  std::string absoluteFileName = "";
  std::string relativeFileName = "";
  bool standardInclude = false;
  std::set<std::string> includePathList;
  std::string getIncludedFile(std::string fileName);
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

  static const std::string HOMEDIRECTORY;
  static constexpr char CONFIGURATIONFILENAME[] = ".tkgeometryrc";
  static constexpr char CONFIGURATIONFILENAMEDEFINITION[] = "TKGEOMETRYRC";
  static constexpr char BINDIRECTORYDEFINITION[] = "TKG_BINDIRECTORY";
  static constexpr char LAYOUTDIRECTORYDEFINITION[] = "TKG_LAYOUTDIRECTORY";
  static constexpr char STANDARDDIRECTORYDEFINITION[] = "TKG_STANDARDDIRECTORY";
  static constexpr char MOMENTADEFINITION[] = "TKG_MOMENTA";
  static constexpr char TRIGGERMOMENTADEFINITION[] = "TKG_TRIGGERMOMENTA";
  static constexpr char THRESHOLDPROBABILITIESDEFINITION[] = "TKG_THRESHOLD_PROB";

  static mainConfigHandler& instance();
  bool getConfiguration(bool checkDirExists = true);
  bool getConfiguration(std::string& layoutDirectory);
  std::string getBinDirectory();
  std::string getLayoutDirectory();
  std::string getStandardDirectory();
  std::string getStyleDirectory();
  std::string getXmlDirectory();
  std::string getMattabDirectory();
  std::string getIrradiationDirectory();
  std::string getDefaultMaterialsDirectory();
  std::string getDetIdSchemesDirectory();
  std::vector<int> getDetIdScheme(std::string schemeName);
  std::string getStandardIncludeDirectory();
  std::string getGeometriesDirectory();
  std::string getConfigFileName();
  std::set<std::string> preprocessConfiguration(ConfigInputOutput);
  std::vector<double>& getMomenta();
  std::vector<double>& getTriggerMomenta();
  std::vector<double>& getThresholdProbabilities();

private:
  bool goodConfigurationRead_;
  std::map<std::string, std::vector<int> > detIdSchemes_;
  std::string binDirectory_;
  std::string layoutDirectory_;
  std::string standardDirectory_;
  std::vector<double> momenta_; // [MeV/c]
  std::vector<double> triggerMomenta_; // [MeV/c]
  std::vector<double> thresholdProbabilities_; // [0,1]
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
  std::string getDetIdSchemesDirectory_();
  void readDetIdSchemes();
  std::string getStandardIncludeDirectory_();
  std::string getGeometriesDirectory_();
};

#endif

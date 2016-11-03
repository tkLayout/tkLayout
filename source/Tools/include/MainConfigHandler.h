#ifndef _MAINCONFIGHANDLER_H__
#define _MAINCONFIGHANDLER_H__

#include <iostream>
#include <string>
#include <vector>
#include <set>

class GraphVizCreator;

/*
 * @class MainConfigHandler
 * @brief Singleton class reading-in the tkLayout configuration (only once) from .tkgeometryrc file and
 * providing the information through out the software.
 * @details Singleton class reading-in the tkLayout configuration (only once) from .tkgeometryrc file and
 * providing the information through out the software. If the configuration file is missing, user will be
 * directly asked through std::cin on the fly for required information. The corresponding configurations
 * will be saved in a file (defined by c_CONFIGURATIONFILENAME variable) in the home directory. The Handler
 * also provides a mechanism how preprocess any configuration file recursive using of @include/@includestd
 * pragmas -> used in geometry boost tree, when building full tracker geometry.
 */
class MainConfigHandler {

 public:

  //! Destructor
  ~MainConfigHandler() {};

  //! MainConfigHandler access method -> get instance of singleton class MainConfigHandler
  static MainConfigHandler& getInstance();

  //! MainConfigHandler method to access GraphVizCreator member (Used to create a graph description file illustrating recursive @include complexity)
  static GraphVizCreator& getGraphVizCreator();

  //! Helper method to preprocess any input configuration file (is) (its full address specified by istreamid) ->
  //! recursively pass through all included files specified by @include & @includestd/@include-std pragma and
  //! include its content to one big configuration file (defined as os)
  std::set<std::string> preprocessConfiguration(std::istream& is, std::ostream& os, const std::string& istreamid, bool standardIncl=false);

  // Getter methods - checking that configuration file read-in first by this singleton class
  bool        getConfiguration(bool checkDirExists = true);
  bool        getConfiguration(std::string& layoutDirectory);
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
  std::string getProjectName();
  std::string getResultsAuthor();

  //bool getConfiguration(std::string& layoutDirectory, std::string& xmlDirectory);

  std::vector<double>& getMomenta();
  std::vector<double>& getTriggerMomenta();
  std::vector<double>& getThresholdProbabilities();

 private:

  //! Singleton constructor method -> must be private (singleton pattern)
  MainConfigHandler();

  //! Check that given directory exists on the filesystem
  bool checkDirectory(std::string dirName) ;

  //! Ask user on standard input for bin directory, where to place executables
  void askBinDirectory();

  //! Ask user on standard input for layout directory, where to place the www output
  void askLayoutDirectory();

  //! Ask user on standard input for standard directory. xml files and other various output will be put here.
  void askStandardDirectory();

  //! Ask user on standard input for particle momenta to be studied
  void askMomenta();

  //! Ask user on standard input for particle momenta to be studied by trigger
  void askTriggerMomenta();

  //! Specify the list of trigger efficiency to be used for pt threshold find scan
  void askThresholdProbabilities();

  //! Create base configuration file for tkLayout .tkgeometryrc, if doesn't exist
  bool createConfigurationFileFromQuestions(std::string& configFileName);

  //! Helper method to read line
  bool parseLine(const char* codeLine, std::string& parameter, std::string& value);

  //! Helper method to convert string list of doubles to std::vector<double>
  std::vector<double> parseDoubleList(std::string);

  //! Read configuration from c_CONFIGURATIONFILENAME or ask user to define configuration on the fly if doesn't exist
  bool readConfiguration(bool checkDirExists);

  //! Helper method called by readConfiguration()
  bool readConfigurationFile(std::string& configFileName);

  // Private getter methods returning required configuration quantity
  std::string getBinDirectory_P();
  std::string getLayoutDirectory_P();
  std::string getStandardDirectory_P();
  std::string getStyleDirectory_P();
  std::string getXmlDirectory_P();
  std::string getMattabDirectory_P();
  std::string getIrradiationDirectory_P();
  std::string getDefaultMaterialsDirectory_P();
  std::string getStandardIncludeDirectory_P();
  std::string getGeometriesDirectory_P();

  bool m_goodConfigurationRead;          //!< Returns true if configuration read-in correctly

  std::string m_binDirectory;            //!< bin directory, where to place executables
  std::string m_layoutDirectory;         //!< layout directory, where to place the www output
  std::string m_standardDirectory;       //!< xml files and other various output will be put here.
  std::string m_projectName;             //!< software project used to produce tkLayout results
  std::string m_resultsAuthor;           //!< name of author producing tkLayout results
  //std::string m_styleDirectory;
  //std::string m_xmlDirectory;

  std::vector<double> m_momenta;                //!< Studied particle momenta
  std::vector<double> m_triggerMomenta;         //!< Studied particle momenta for trigger (i.e. strip detector, outer detector)
  std::vector<double> m_thresholdProbabilities; //!< Studied trigger efficiencies

  const char* c_HOMEDIRECTORY                   = "HOME";
  const char* c_CONFIGURATIONFILENAME           = ".tkgeometryrc";
  const char* c_CONFIGURATIONFILENAMEDEFINITION = "TKGEOMETRYRC";
  const char* c_BINDIRECTORYDEFINITION          = "TKG_BINDIRECTORY";
  const char* c_LAYOUTDIRECTORYDEFINITION       = "TKG_LAYOUTDIRECTORY";
  const char* c_STANDARDDIRECTORYDEFINITION     = "TKG_STANDARDDIRECTORY";
  const char* c_MOMENTADEFINITION               = "TKG_MOMENTA";
  const char* c_PROJECTNAME                     = "TKG_PROJECT";
  const char* c_RESULTSAUTHOR                   = "TKG_AUTHOR";
  const char* c_TRIGGERMOMENTADEFINITION        = "TKG_TRIGGERMOMENTA";
  const char* c_THRESHOLDPROBABILITIESDEFINITION= "TKG_THRESHOLD_PROB";

}; // Class

/*
 * Configuration input/output helper class
 */
class ConfigInputOutput {

public:

  ConfigInputOutput(std::istream& newIs, std::ostream& newOs, const std::string& name) :
    is(newIs) ,
    os(newOs),
    absFileName(name),
    relFileName(name) {}

  std::string getIncludedFile(std::string fileName);
  std::string getAbsFileNameDir() const;

  std::istream& is;
  std::ostream& os;
  std::string absFileName;
  std::string relFileName;
  bool standardInclude = false;
  std::set<std::string> includePathList;
  bool webOutput = false;
}; // Class

#endif

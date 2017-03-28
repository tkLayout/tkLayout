// 
// File:   Squid.h
// Author: ndemaio
//
// Created on May 29, 2009, 12:09 PM
//

/**
 * @file Squid.h
 * @brief This is the header file for the integrating wrapper class for the material budget operations of <i>tkgeometry</i>
 */

#ifndef _SQUID_H
#define	_SQUID_H

#include <string>
#include <MatParser.hh>
#include <InactiveSurfaces.hh>
#include <MaterialBudget.hh>
#include <Usher.hh>
#include <MatCalc.hh>
#include <Analyzer.hh>
#include <Vizard.hh>
#include <tk2CMSSW.hh>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/program_options/variables_map.hpp>
#include <RootWeb.hh>
#include <MainConfigHandler.hh>
#include <MessageLogger.hh>

#include <Tracker.hh>
#include <Support.hh>
#include "Materialway.hh"
#include "WeightDistributionGrid.hh"


using material::Materialway;
using material::WeightDistributionGrid;

namespace po = boost::program_options;
/**
 * A shorter alias for the filesystem library namespace
 */
namespace bfs = boost::filesystem;
namespace insur {
  /*
   * Error messages and warnings that may be reported.
   */
  static const std::string err_no_geomfile = "There is no recorded name for the geometry configuration file. Initialise the tracker first.";
  static const std::string err_no_matfile = "The provided material configuration file does not exist.";
  static const std::string err_no_matfile_pixel = "The material configuration file for the pixels does not exist.";
  static const std::string err_init_failed = "Initialisation of the material calculator failed.";
  static const std::string err_no_tracker = "The tracker object does not exist. The tracker must be created before calling this function.";
  static const std::string err_no_inacsurf = "The collection of inactive surfaces does not exist. It must be created before calling this function";
  static const std::string err_no_matbudget = "The material budget does not exist. It must be created before calling this function.";
  static const std::string err_no_triggerSummary = "Could not report on the trigger performance.";
  static const std::string warn_rootonly = "The collection of inactive surfaces does not exist. Only the .root file will be written.";
  static const std::string warn_custom_matfile = "A customized material file was used for the tracker";
  static const std::string warn_custom_matfile_pixel = "A customized material file was used for the pixel";

  static const std::string default_trackername = "defaultTrackerName";

  /**
   * @class Squid
   * @brief The Squid class integrates the components of the <i>tkmaterial</i> application and provides an
   * interface to its high-level functionality.
   *
   * It manages instances to all the necessary components internally. Its functions bundle the steps that are
   * required to carry out the main tasks of the application. The idea is to allow access to the underlying
   * via a series of well-defined pathways contained in a single class. This class acts as a sort of manager 
   * to the application and can be embedded into a main program easily. It maintains internal instances of
   * all the objects that are necessary to run the modelling, analysis and visualisation steps.
   */
  class Squid {
  public:
    Squid();
    virtual ~Squid();
    bool buildTracker();
    //bool dressTracker();
    //bool buildTrackerSystem();
    //bool irradiateTracker();
    bool buildInactiveSurfaces(bool verbose = false);
    bool buildMaterials(bool verbose = false);
    bool createMaterialBudget(bool verbose = false);
    //bool buildFullSystem(bool usher_verbose = false, bool mat_verbose = false);
    bool analyzeNeighbours(std::string graphout = "");
    bool translateFullSystemToXML(std::string xmlout = "");
    bool webOutput = false;

    // Functions using rootweb
    bool analyzeTriggerEfficiency(int tracks, bool detailed);
    bool pureAnalyzeGeometry(int tracks);
    bool pureAnalyzeMaterialBudget(int tracks, bool triggerRes, bool triggerPatternReco, bool debugResolution);
    bool reportGeometrySite(bool debugResolution);
    bool reportBandwidthSite();
    bool reportTriggerProcessorsSite();
    bool reportPowerSite();
    bool reportMaterialBudgetSite(bool debugServices);
    bool reportResolutionSite();
    bool reportPatternRecoSite();
    bool reportTriggerPerformanceSite(bool extended);
    bool reportNeighbourGraphSite();
    bool additionalInfoSite();
    bool makeSite(bool addLogPage = true);
    void setBasename(std::string newBaseName);
    void setGeometryFile(std::string geomFile);
    void setHtmlDir(std::string htmlDir);

    void simulateTracks(const po::variables_map& varmap, int seed);
    void setCommandLine(int argc, char* argv[]);
    //void pixelExtraction(std::string xmlout);
    void createAdditionalXmlSite(std::string xmlout);
  private:
    //std::string g;
    Tracker* tr;
    InactiveSurfaces* is;
    MaterialBudget* mb;
    Tracker* px;
    std::list<Support*> supports_;
    InactiveSurfaces* pi;
    MaterialBudget* pm;
    MatParser mp;
    Usher u;
    Materialway materialwayTracker;
    Materialway materialwayPixel;
    MatCalc tkMaterialCalc;
    MatCalc pxMaterialCalc;
    Analyzer a;
    Analyzer pixelAnalyzer;
    Vizard v;
    mainConfigHandler& mainConfiguration;
    tk2CMSSW t2c;
    bool fileExists(std::string filename);
    std::string extractFileName(const std::string& full);
    Squid(const Squid& s);
    Squid& operator=(const Squid& s);
    void resetVizard();
    std::string baseName_;
    std::string htmlDir_;
    std::string getGeometryFile();
    std::string getSettingsFile();
    std::string getMaterialFile();
    std::string getPixelMaterialFile();
    std::string myGeometryFile_;
    std::string mySettingsFile_;
    std::string myMaterialFile_;
    std::string myPixelMaterialFile_;
    std::set<std::string> includeSet_; // list of configuration files
    bool defaultMaterialFile;
    bool defaultPixelMaterialFile;

    RootWSite site;

    WeightDistributionGrid weightDistributionTracker;
    WeightDistributionGrid weightDistributionPixel;

    bool prepareWebsite();
    bool sitePrepared;
  };
}
#endif	/* _SQUID_H */


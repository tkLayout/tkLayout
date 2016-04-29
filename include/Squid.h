// 
// File:   Squid.h
// Author: ndemaio
//
// Created on May 29, 2009, 12:09 PM
// Rearranged in 2016 for FCC-hh project by Z. Drasal (CERN)
//

/**
 * @file Squid.h
 * @brief This is the header file for the integrating wrapper class for the material budget operations of <i>tkgeometry</i>
 */

#ifndef _SQUID_H
#define	_SQUID_H

#include <string>
#include <MatParser.h>
#include <InactiveSurfaces.h>
#include <MaterialBudget.h>
#include <Usher.h>
#include <MatCalc.h>
#include <Analyzer.h>
#include <Vizard.h>
#include <tk2CMSSW.h>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/program_options/variables_map.hpp>
#include <rootweb.hh>
#include <MainConfigHandler.h>
#include <MessageLogger.h>

#include <Tracker.h>
#include <Support.h>
#include "Materialway.h"
#include "WeightDistributionGrid.h"

class AnalyzerGeometry;

using material::Materialway;
using material::WeightDistributionGrid;

namespace po = boost::program_options;
/**
 * A shorter alias for the filesystem library namespace
 */
namespace bfs = boost::filesystem;
namespace insur {
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

    //Error messages and warnings that may be reported.
    static const std::string err_no_geomfile;
    static const std::string err_no_matfile;
    static const std::string err_no_matfile_pixel;
    static const std::string err_init_failed;
    static const std::string err_no_tracker;
    static const std::string err_no_inacsurf;
    static const std::string err_no_matbudget;
    static const std::string err_no_triggerSummary;
    static const std::string err_no_flukafile;
    static const std::string err_no_occupancy_failed;
    static const std::string warn_rootonly;
    static const std::string warn_custom_matfile;
    static const std::string warn_custom_matfile_pixel;
    static const std::string default_trackername;

    Squid();
    virtual ~Squid();

    // Build a bare-bones geometry of active tracker, i.e. active modules
    bool buildActiveTracker();
    // Build all pasive components (inactive surfaces) of active tracker
    bool buildPasiveTracker(bool verbose = false);

    bool createMaterialBudget(bool verbose = false);
    bool analyzeNeighbours(std::string graphout = "");
    bool translateFullSystemToXML(std::string xmlout = "");

    // Functions using rootweb
    //bool analyzeTriggerEfficiency(int tracks, bool detailed);

    // Check that tracker exists & analyze its layout geometry (no output through rootweb)
    bool analyzeGeometry(int tracks);


    bool pureAnalyzeMaterialBudget(int tracks);
    bool pureAnalyzeResolution(int tracks);
    bool reportGeometrySite();
    bool reportBandwidthSite();
    bool reportTriggerProcessorsSite();
    bool reportOccupancySite();
    bool reportPowerSite();
    bool reportMaterialBudgetSite(bool debugServices);
    bool reportResolutionSite();
    bool reportTriggerPerformanceSite(bool extended);
    bool reportNeighbourGraphSite();
    bool reportInfoSite();
    bool makeWebSite(bool addLogPage = true);
    void setBasename(std::string newBaseName);
    void setGeometryFile(std::string geomFile);
    void setHtmlDir(std::string htmlDir);

    // Other
    void simulateTracks(const po::variables_map& varmap, int seed);
    void setCommandLine(int argc, char* argv[]);

  private:

    // Trackers
    Tracker* m_ctrlPxd;
    Tracker* m_ctrlStd;
    Tracker* m_fwdPxd;
    Tracker* m_fwdStd;

    std::vector<Tracker*> m_trackers;

    // Inactive surfaces
    InactiveSurfaces* m_pasiveCtrlPxd;
    InactiveSurfaces* m_pasiveFwdPxd;
    InactiveSurfaces* m_pasiveCtrlStd;
    InactiveSurfaces* m_pasiveFwdStd;

    WeightDistributionGrid m_weightDistCtrlPxd;
    WeightDistributionGrid m_weightDistFwdPxd;
    WeightDistributionGrid m_weightDistCtrlStd;
    WeightDistributionGrid m_weightDistFwdStd;

    // Support structures
    std::list<Support*> m_supports;

    // Material budget containers
    MaterialBudget* m_mbBP;
    MaterialBudget* m_mbCtrlPxd;
    MaterialBudget* m_mbFwdPxd;
    MaterialBudget* m_mbCtrlStd;
    MaterialBudget* m_mbFwdStd;

    // Tracker analyzers
    Analyzer m_pixelAnalyzer;
    Analyzer m_stripAnalyzer;
    Analyzer m_fwdPxdAnalyzer;
    Analyzer m_fwdAnalyzer;
    Analyzer m_trackerAnalyzer;

    // Geometry layout analyzer
    AnalyzerGeometry* m_geomAnalyzer;

    // Visualization manager
    Vizard m_vizard;
    void   resetVizard();

    // Website container
    RootWSite m_webSite;        // General web container
    bool      prepareWebSite();
    bool      m_sitePrepared;

    Usher u;

    // Transform tkLayout geometry to CMSSW
    tk2CMSSW t2c;

    // File operations
    bool        fileExists(std::string filename);
    std::string extractFileName(const std::string& full);

    // tkLayout configuration & path variables
    std::string m_baseName;
    std::string m_htmlDir;
    std::string m_layoutName;
    std::string m_myGeometryFile;
    std::string m_mySettingsFile;

    std::string getGeometryFile();
    std::string getSettingsFile();
    std::string getMaterialFile();
    std::string getPixelMaterialFile();

    std::set<std::string> m_includeSet; // list of configuration files

  };
}
#endif	/* _SQUID_H */


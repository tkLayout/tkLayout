// System includes
#include <boost/program_options.hpp> 
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

// tkLayout includes
#include <global_constants.h>
#include <AnalysisManager.h>
#include <GeometryManager.h>
#include <SimParms.h>
#include <StopWatch.h>
#include <GitRevision.h>

namespace po = boost::program_options;

/*
 * \mainpage The tkLayout software
 */
int main(int argc, char* argv[]) {

  std::string usage("Usage: ");
  usage += argv[0];
  usage += " <geometry file> [options]";

  int geomTracks, matTracks;
  int verbosity;

  std::string geomFile, optFile;
  std::string expOption;

  // Program options - analysis
  po::options_description shown("Analysis options");
  shown.add_options()
    ("help,h"           , "Display help info.")
    ("opt-file"         , po::value<std::string>(&optFile)->implicit_value(""), "Specify an option file to parse the program options from (in addition to command line).")
    ("geometry-tracks,n", po::value<int>(&geomTracks)->default_value(default_n_tracks), "Number of tracks for geometry calculations.")
    ("material-tracks,N", po::value<int>(&matTracks)->default_value(default_n_tracks), "Number of tracks for material & resolution calculations.")
    ("occupancy,o"      , "Report occupancy studies based on Fluka data.")
    ("geometry,g"       , "Report geometry layout.")
    ("material,m"       , "Report material budget.")
    ("patternreco,p"    , "Report pattern recognition capabilities in given occupancy (using Fluka data).")
    ("resolution,r"     , "Report resolution studies.")
    ("all,a"            , "Report all studies.")
    ("extraction,e"     , po::value<std::string>(&expOption)->implicit_value("CMS"), "Extract tkLayout geometry to an XML file to be used in CMS/FCC SW frameworks. Supported values: CMS or FCC.")
    ("verbosity"        , po::value<int>(&verbosity)->default_value(1), "Verbosity level (Overridden by the option 'quiet').")
    ("quiet"            , "No output produced (Equivalent to verbosity 0, overrides the 'verbosity' option).")
    ("performance"      , "Outputs the CPU time needed for each computing step (Overrides the option 'quiet').")
    ;

  // Program options - other
  po::options_description otheropt("Other options");
  otheropt.add_options()("version,v", "Prints software version (GIT revision) and quits.");

  po::options_description hidden;
  hidden.add_options()("geom-file", po::value<std::string>(&geomFile));

  // Geometry/config file is defined without -- or - selector, that's defined as positional option
  po::positional_options_description posopt;
  posopt.add("geom-file", 1);

  po::options_description mainopt;
  mainopt.add(shown).add(hidden).add(otheropt);
  
  // Read user selected options
  po::variables_map progOptions;
  try {

    // Find all options including positional options
    po::store(po::command_line_parser(argc, argv).options(mainopt).positional(posopt).run(), progOptions);
    po::notify(progOptions);
    if (!optFile.empty()) {
      std::ifstream inFile(optFile.c_str());
      if (!inFile) throw po::error(("Options file \"" + optFile + "\" not found").c_str());

      // Parse also the optional input file
      po::store(po::parse_config_file(inFile, mainopt, true), progOptions);
      inFile.close();
    }
    po::notify(progOptions);

    // Check number of tracks
    if (geomTracks < 1) throw po::invalid_option_value("geometry-tracks");
    if (matTracks < 1)  throw po::invalid_option_value("material-tracks");
    if (!progOptions.count("geom-file") && !progOptions.count("help") && !progOptions.count("version")) throw po::error("Missing configuration/geometry file");

    // Check supported experimental option for full geometry extraction
    if (progOptions.count("extraction") && expOption!="CMS" && expOption!="FCC") throw po::invalid_option_value("extraction");

  } catch(po::error& e) {

    // Display the following options after program starts without configuration
    logERROR(e.what());
    std::cout << std::endl << usage << std::endl << shown << std::endl;
    return EXIT_FAILURE;
  }

  // Print help
  if (progOptions.count("help")) {
    std::cout << usage << std::endl << shown << otheropt << std::endl;
    return EXIT_SUCCESS;
  }

  // Print version
  if (progOptions.count("version")) {
    std::cout << "tklayout revision " << GitRevision::revisionNumber << std::endl;
    return EXIT_SUCCESS;
  }

  // Print geometry configuration file to be used
  if (progOptions.count("geom-file")) std::cout << "Reading configuration/geometry from: " << geomFile << std::endl;

  //
  // Set verbosity
  bool         verboseMaterial  = false;
  unsigned int verboseWatch     = verbosity;
  bool         performanceWatch = false ;

  if (progOptions.count("quiet"))       verboseWatch = 0;
  if (progOptions.count("performance")) verboseWatch = 1;

  StopWatch::getInstance().setVerbosity(verboseWatch, performanceWatch);

  //
  // Build full detector: tracker (active + passive parts) + beam pipe
  GeometryManager gManager(geomFile);
  bool activeDetOK  = gManager.buildActiveDetector();
  bool passiveDetOK = gManager.buildPassiveDetector();

  //
  // Set simulation generic parameters
  auto& simParms = SimParms::getInstance();
  simParms.setCommandLine(argc, argv);
  simParms.setLayoutName(gManager.getLayoutName());
  simParms.setBaseGeomFileName(gManager.getBaseGeomFileName());
  simParms.setRunDirPath(gManager.getRunDirPath());
  simParms.setWebDir(gManager.getWebDir());
  simParms.setListOfConfFiles(gManager.getListOfConfFiles());

  //
  // Analyze tracker - create analysis manager
  AnalysisManager aManager(gManager.getDetector());

  // Call individual analyzer units
  bool isAnalysisOK      = false;
  bool isVisualizationOK = false;
  if (activeDetOK) {

    // Geometry layout study
    if (progOptions.count("all") || progOptions.count("geometry") || progOptions.count("material") || progOptions.count("resolution") || progOptions.count("patternreco") || progOptions.count("occupancy")) {

      startTaskClock("Analyzing tracker geometry");
      aManager.initUnit(geomTracks, "AnalyzerGeometry");
      isAnalysisOK      = aManager.analyzeUnit("AnalyzerGeometry");
      isVisualizationOK = aManager.visualizeUnit("AnalyzerGeometry");
      stopTaskClock();
      if (!isAnalysisOK)      logERROR("Error in AnalyzerGeometry -> analysis failed!");
      if (!isVisualizationOK) logERROR("Error in AnalyzerGeometry -> visualization failed!");
    }

    // Material budget study
    if ((progOptions.count("all") || progOptions.count("material") || progOptions.count("resolution")) && passiveDetOK) {

      startTaskClock("Analyzing material budget");
      aManager.initUnit(matTracks, "AnalyzerMatBudget");
      isAnalysisOK      = aManager.analyzeUnit("AnalyzerMatBudget");
      isVisualizationOK = aManager.visualizeUnit("AnalyzerMatBudget");
      stopTaskClock();
      if (!isAnalysisOK)      logERROR("Error in AnalyzerMatBudget -> analysis failed!");
      if (!isVisualizationOK) logERROR("Error in AnalyzerMatBudget -> visualization failed!");
    }

    // Resolution study
    if ((progOptions.count("all") || progOptions.count("resolution")) && passiveDetOK) {

      startTaskClock("Analyzing resolution");
      aManager.initUnit(matTracks, "AnalyzerResolution");
      isAnalysisOK      = aManager.analyzeUnit("AnalyzerResolution");
      isVisualizationOK = aManager.visualizeUnit("AnalyzerResolution");
      stopTaskClock();
      if (!isAnalysisOK)      logERROR("Error in AnalyzerResolution -> analysis failed!");
      if (!isVisualizationOK) logERROR("Error in AnalyzerResolution -> visualization failed!");
    }

    // Pattern-recognition study
    if ((progOptions.count("all") || progOptions.count("patternreco")) && passiveDetOK) {

      startTaskClock("Analyzing pattern recognition");
      aManager.initUnit(matTracks, "AnalyzerPatternReco");
      isAnalysisOK      = aManager.analyzeUnit("AnalyzerPatternReco");
      isVisualizationOK = aManager.visualizeUnit("AnalyzerPatternReco");
      stopTaskClock();
      if (!isAnalysisOK)      logERROR("Error in AnalyzerPatternReco -> analysis failed!");
      if (!isVisualizationOK) logERROR("Error in AnalyzerPatternReco -> visualization failed!");
    }

    // Occupancy study
    if (progOptions.count("all") || progOptions.count("occupancy")) {

      startTaskClock("Analyzing occupancy");
      aManager.initUnit(matTracks, "AnalyzerOccupancy");
      isAnalysisOK      = aManager.analyzeUnit("AnalyzerOccupancy");
      isVisualizationOK = aManager.visualizeUnit("AnalyzerOccupancy");
      stopTaskClock();
      if (!isAnalysisOK)      logERROR("Error in AnalyzerOccupancy -> analysis failed!");
      if (!isVisualizationOK) logERROR("Error in AnalyzerOccupancy -> visualization failed!");
    }

    // Geometry extractor to CMSSW
    if (progOptions.count("extraction") && expOption=="CMS") {

      startTaskClock("Extracting geometry for use in CMSSW");
      aManager.initUnit(matTracks, "ExtractorCMSSW");
      isAnalysisOK      = aManager.analyzeUnit("ExtractorCMSSW");
      stopTaskClock();
      if (!isAnalysisOK) logERROR("Error in ExtractorCMSSW -> extraction failed!");
    }

    // Geometry extractor to FCCSW
    if (progOptions.count("extraction") && expOption=="FCC") {

      startTaskClock("Extracting geometry for use in FCCSW");
      aManager.initUnit(matTracks, "ExtractorFCCSW");
      isAnalysisOK      = aManager.analyzeUnit("ExtractorFCCSW");
      stopTaskClock();
      if (!isAnalysisOK) logERROR("Error in ExtractorFCCSW -> extraction failed!");
    }
  }

  //
  // Create html output with all results
  bool addInfoPage = true;
  bool addLogPage  = true;
  aManager.makeWebSite(addInfoPage, addLogPage);

  // Finish
  std::cout << std::endl;
  return EXIT_SUCCESS;
}





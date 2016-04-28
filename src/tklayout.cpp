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
#include <StopWatch.h>
#include <SvnRevision.h>

namespace po = boost::program_options;

int main(int argc, char* argv[]) {

  std::string usage("Usage: ");
  usage += argv[0];
  usage += " <geometry file> [options]";

  int geomTracks, matTracks;
  int verbosity;

  std::string geomFile, optFile;

  // Program options - analysis
  po::options_description shown("Analysis options");
  shown.add_options()
    ("help,h"           , "Display help info.")
    ("opt-file"         , po::value<std::string>(&optFile)->implicit_value(""), "Specify an option file to parse the program options from (in addition to command line).")
    ("geometry-tracks,n", po::value<int>(&geomTracks)->default_value(insur::default_n_tracks), "Number of tracks for geometry calculations.")
    ("material-tracks,N", po::value<int>(&matTracks)->default_value(insur::default_n_tracks), "Number of tracks for material & resolution calculations.")
    ("occupancy,o"      , "Report occupancy studies based on Fluka data.")
    ("geometry,g"       , "Report geometry layout.")
    ("material,m"       , "Report material buget.")
    ("resolution,r"     , "Report resolution studies.")
    ("all,a"            , "Report all studies.")
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

  } catch(po::error& e) {

    // Display the following options after program starts without configuration
    std::cerr << "\nERROR: " << e.what() << std::endl << std::endl;
    std::cout << usage << std::endl << shown << std::endl;
    return EXIT_FAILURE;
  }

  // Print help
  if (progOptions.count("help")) {
    std::cout << usage << std::endl << shown << otheropt << std::endl;
    return EXIT_SUCCESS;
  }

  // Print version
  if (progOptions.count("version")) {
    std::cout << "tklayout revision " << SvnRevision::revisionNumber << std::endl;
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

  StopWatch::instance()->setVerbosity(verboseWatch, performanceWatch);

  //
  // Build full tracker
  GeometryManager gManager(geomFile);
  bool activeTrackerOK = gManager.buildActiveTracker();

  //
  // Analyze tracker - create analysis manager
  AnalysisManager aManager(gManager.getLayoutName(), gManager.getWebDir(), gManager.getActiveTrackers(), gManager.getPasiveTrackers(), gManager.getTrackerSupports());
  aManager.setCommandLine(argc, argv);

  // Call individual analyzer modules
  bool isAnalysisOK      = false;
  bool isVisualizationOK = false;
  if (activeTrackerOK) {

    // Geometry layout study
    if (progOptions.count("all") || progOptions.count("geometry") || progOptions.count("material") || progOptions.count("resolution")) {

      startTaskClock("Analyzing tracker geometry");
      aManager.initModule(geomTracks, "AnalyzerGeometry");
      isAnalysisOK      = aManager.analyzeModule("AnalyzerGeometry");
      isVisualizationOK = aManager.visualizeModule("AnalyzerGeometry");
      stopTaskClock();
      if (!isAnalysisOK)      std::cerr << "\nERROR in AnalyzerGeometry -> analysis failed!"      << std::endl;
      if (!isVisualizationOK) std::cerr << "\nERROR in AnalyzerGeometry -> visualization failed!" << std::endl;
    }

    // Material budget study
    if (progOptions.count("all") || progOptions.count("material")) {

      startTaskClock("Analyzing material budget");
      aManager.initModule(matTracks, "AnalyzerMatBudget");
      isAnalysisOK      = aManager.analyzeModule("AnalyzerMatBudget");
      isVisualizationOK = aManager.visualizeModule("AnalyzerMatBudget");
      stopTaskClock();
      if (!isAnalysisOK)      std::cerr << "\nERROR in AnalyzerMatBudget -> analysis failed!"      << std::endl;
      if (!isVisualizationOK) std::cerr << "\nERROR in AnalyzerMatBudget -> visualization failed!" << std::endl;
    }

    // Resolution study
    if (progOptions.count("all") || progOptions.count("resolution")) {

      startTaskClock("Analyzing resolution");
      aManager.initModule(matTracks, "AnalyzerResolution");
      isAnalysisOK      = aManager.analyzeModule("AnalyzerResolution");
      isVisualizationOK = aManager.visualizeModule("AnalyzerResolution");
      stopTaskClock();
      if (!isAnalysisOK)      std::cerr << "\nERROR in AnalyzerResolution -> analysis failed!"      << std::endl;
      if (!isVisualizationOK) std::cerr << "\nERROR in AnalyzerResolution -> visualization failed!" << std::endl;
    }
  }

  //
  // Create html output with all results
  bool addInfoPage = false;
  bool addLogPage  = true;
  aManager.makeWebSite(addInfoPage, addLogPage);

  // Finish
  std::cout << std::endl;
  return EXIT_SUCCESS;
}





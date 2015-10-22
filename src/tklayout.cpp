#include <boost/program_options.hpp> 
#include <stdlib.h>
#include <iostream>
#include <string>
#include <Squid.h>
#include <global_constants.h>
#include "SvnRevision.h"

namespace po = boost::program_options;

int main(int argc, char* argv[]) {

  std::string usage("Usage: ");
  usage += argv[0];
  usage += " <geometry file> [options]";

  int geomTracks, matTracks;
  int verbosity;
  int randseed; 

  std::string geomFile, optFile, xmlDir, htmlDir;

  // Program options - analysis
  po::options_description shown("Analysis options");
  shown.add_options()
    ("help,h", "Display this help message.")
    ("opt-file", po::value<std::string>(&optFile)->implicit_value(""), "Specify an option file to parse program options from, in addition to the command line")
    ("geometry-tracks,n", po::value<int>(&geomTracks)->default_value(insur::default_n_tracks), "N. of tracks for geometry calculations.")
    ("material-tracks,N", po::value<int>(&matTracks)->default_value(insur::default_n_tracks), "N. of tracks for material & resolution calculations.")
    ("power,p", "Report irradiated power analysis.")
    ("bandwidth,b", "Report base bandwidth analysis.")
    ("bandwidth-cpu,B", "Report multi-cpu bandwidth analysis.\n\t(implies 'b')")
    ("material,m", "Report materials and weights analyses.")
    ("resolution,r", "Report resolution analysis.")
    ("trigger,t", "Report base trigger analysis.")
    ("trigger-ext,T", "Report extended trigger analysis.\n\t(implies 't')")
    ("debug-services,d", "Service additional debug info")
    ("all,a", "Report all analyses, except extended\ntrigger and debug page. (implies all other relevant\nreport options)")
    ("graph,g", "Build and report feeder/neighbour relations graph.")
    ("xml", po::value<std::string>(&xmlDir)->implicit_value(""), "Produce XML output files for materials.\nOptional arg specifies the subdirectory\nof the output directory (chosen via inst\nscript) where to create XML files.\nIf not supplied, the config file name (minus extension)\nwill be used as subdir.")
    ("html-dir", po::value<std::string>(&htmlDir), "Override the default html output dir\n(equal to the tracker name in the main\ncfg file) with the one specified.")
    ("verbosity", po::value<int>(&verbosity)->default_value(1), "Levels of details in the program's output (overridden by the option 'quiet').")
    ("quiet", "No output is produced, except the required messages (equivalent to verbosity 0, overrides the option 'verbosity')")
    ("performance", "Outputs the CPU time needed for each computing step (overrides the option 'quiet').")
    ("randseed", po::value<int>(&randseed)->default_value(0xcafebabe), "Set the random seed\nIf explicitly set to 0, seed is random")
    ;

  // Program options - simulation
  po::options_description trackopt("Track simulation options");
  trackopt.add_options()
    ("tracksim", "Switch to track sim mode, normal analysis disabled")
    ("num-events", po::value<std::string>(), "N. of events to simulate")
    ("num-tracks-ev", po::value<std::string>(), "N. of tracks per event")
    ("event-offset", po::value<std::string>(), "Start the event numbering from an offset value.")
    ("eta", po::value<std::string>(), "Particle eta")
    ("phi0", po::value<std::string>(), "Particle phi0")
    ("z0", po::value<std::string>(), "Particle z0")
    ("pt", po::value<std::string>(), "Particle transverse momentum, alternative to --invPt")
    ("invPt", po::value<std::string>(), "Specify pt in terms of its inverse, alternative to --pt")
    ("charge", po::value<std::string>(), "Particle charge")
    ("instance-id", po::value<std::string>(), "Id of the program instance, to tag the output file with")
    ("tracks-dir", po::value<std::string>(), "Override the default tracksim output dir.\nIf not supplied, the files will be saved in\nthe working dir")
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
  mainopt.add(shown).add(hidden).add(trackopt).add(otheropt);
  
  // Read user input
  po::variables_map vm;
  try {

    // Find all options including positional options
    po::store(po::command_line_parser(argc, argv).options(mainopt).positional(posopt).run(), vm);
    po::notify(vm);
    if (!optFile.empty()) {
      std::ifstream ifs(optFile.c_str());
      if (!ifs) throw po::error(("Options file \"" + optFile + "\" not found").c_str());

      // Parse also the optional input file
      po::store(po::parse_config_file(ifs, mainopt, true), vm);
      ifs.close();
    }
    po::notify(vm);

    // Check number of tracks
    if (geomTracks < 1) throw po::invalid_option_value("geometry-tracks");
    if (matTracks < 1)  throw po::invalid_option_value("material-tracks");
    if (!vm.count("geom-file") && !vm.count("help") && !vm.count("version")) throw po::error("Missing configuration/geometry file");

  } catch(po::error& e) {

    // Display the following options after program starts without configuration
    std::cerr << "\nERROR: " << e.what() << std::endl << std::endl;
    std::cout << usage << std::endl << shown << trackopt << std::endl;
    return EXIT_FAILURE;
  }

  if (vm.count("help")) {
    std::cout << usage << std::endl << shown << trackopt << otheropt << std::endl;
    return 0;
  }

  if (vm.count("version")) {
    std::cout << "tklayout revision " << SvnRevision::revisionNumber << std::endl;
    return 0;
  }

  if (vm.count("geom-file")) std::cout << "Reading configuration/geometry from: " << geomFile << std::endl;

  //
  // Start program using squid multi-interface
  insur::Squid squid;

  // Set configuration & other options
  squid.setCommandLine(argc, argv);

  // Set html directory
  if (htmlDir != "") squid.setHtmlDir(htmlDir);

  // Set geometry/configuration file
  squid.setGeometryFile(geomFile);

  // Set verbosity
  bool         verboseMaterial  = false;
  unsigned int verboseWatch     = verbosity;
  bool         performanceWatch = false ;
  if (vm.count("quiet")) verboseWatch = 0;
  if (vm.count("performance")) {
    if (verboseWatch==0) verboseWatch = 1;
  }
  StopWatch::instance()->setVerbosity(verboseWatch, performanceWatch);

  //
  // Build tracker: pixel & strip in barrel & forward region
  if (!squid.buildTracker()) return EXIT_FAILURE;

  // Do general analysis
  if (!vm.count("tracksim")) {

    // Analyze geometry
    if (!squid.pureAnalyzeGeometry(geomTracks)) return EXIT_FAILURE;

    //if ((vm.count("all") || vm.count("bandwidth") || vm.count("bandwidth-cpu")) && !squid.reportBandwidthSite())         return EXIT_FAILURE;
    //if ((vm.count("all") || vm.count("bandwidth-cpu"))                          && !squid.reportTriggerProcessorsSite()) return EXIT_FAILURE;
    //if ((vm.count("all") || vm.count("power"))                                  && !squid.reportPowerSite())             return EXIT_FAILURE;

    // If needed build material model & perform resolution simulation
    if ( vm.count("all") || vm.count("material") || vm.count("resolution") || vm.count("graph") || vm.count("xml") ) {
      if (squid.buildMaterials(verboseMaterial) && squid.createMaterialBudget(verboseMaterial)) {

        if ( vm.count("all") || vm.count("material") || vm.count("resolution") ) {

          // Perform material budget analysis & report it
          if ((vm.count("all") || vm.count("resolution") || vm.count("material")) && !squid.pureAnalyzeMaterialBudget(matTracks))                 return EXIT_FAILURE;
          if ((vm.count("all") || vm.count("resolution") || vm.count("material")) && !squid.reportMaterialBudgetSite(vm.count("debug-services"))) return EXIT_FAILURE;

          // Perform resolution analysis & report it
          if ((vm.count("all") || vm.count("resolution")) && !squid.pureAnalyzeResolution(matTracks)) return EXIT_FAILURE;
          if ((vm.count("all") || vm.count("resolution")) && !squid.reportResolutionSite())           return EXIT_FAILURE;
        }

        //TODO: Writes the feeder/neighbour relations in a collection of inactive surfaces to a file
        //if (vm.count("graph") && !squid.reportNeighbourGraphSite())     return EXIT_FAILURE;

        // Produce XML output files for materials
        if (vm.count("xml") && !squid.translateFullSystemToXML(xmlDir)) return (EXIT_FAILURE);
      }
    }

//    if ((vm.count("all") || vm.count("trigger") || vm.count("trigger-ext")) &&
//        ( !squid.analyzeTriggerEfficiency(matTracks, vm.count("trigger-ext")) || !squid.reportTriggerPerformanceSite(vm.count("trigger-ext"))) ) return EXIT_FAILURE;


    // Report tracker geometry
    if (!squid.reportGeometrySite()) return EXIT_FAILURE;
    if (!squid.reportInfoSite())     return EXIT_FAILURE;
    if (!squid.makeWebSite())        return EXIT_FAILURE;

  }
  // Do specific simulated tracks analysis
  else {
    //if (tracksim.size() < 1 || tracksim.size > 2) {
    //  std::cerr << "Wrong number of parameters. Syntax: --tracksim <num events> <num tracks/ev>" << std::endl;
    //  std::cerr << "                                    --tracksim parameterfile" << std::endl;
    //  std::cerr << "                                    --tracksim \"key1 = value1; key2 = value2 ...\"" << std::endl;
    //  return EXIT_FAILURE;
   // }
    if (!squid.pureAnalyzeGeometry(geomTracks)) return EXIT_FAILURE;
  
//    if (tracksim.size() == 2) {
//      vmtracks.insert(std::make_pair("num-events", po::variable_value(boost::any(tracksim[0]), false)));
//      vmtracks.insert(std::make_pair("num-tracks", po::variable_value(boost::any(tracksim[1]), false)));
//    }
    squid.simulateTracks(vm, randseed);

    //if (tracksim.size() == 2) { squid.simulateTracks(str2any<long int>(tracksim[0]), str2any<long int>(tracksim[1]), randseed, "", ""); }
    //else if (tracksim.size() == 1 && tracksim[0].at(0)=="\"") { squid.simulateTracks(0, 0, randseed, "", trim(tracksim[0], " \"")); }
    //else { squid.simulateTracks(0, 0, randseed, tracksim[0], ""); }
  }

  // Finish
  return EXIT_SUCCESS;
}





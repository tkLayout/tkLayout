
#include <boost/program_options.hpp> 
#include <stdlib.h>
#include <iostream>
#include <string>
#include <Squid.h>

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  std::string usage("Usage: ");
  usage += argv[0];
  usage += " <config basename> [options]";
  usage += "\n\n<config basename> is the user-specified part of the config file names.\nFull names are automatically inferred from it by appending appropriate suffixes\n";
  po::options_description shown("Allowed options");
  int geomtracks, mattracks;
  std::vector<int> tracksim;
  int verbosity;

  std::string basename, xmldir, htmldir;

  shown.add_options()
    ("help", "Display this help message.")
    ("geometry-tracks,n", po::value<int>(&geomtracks)->default_value(100), "N. of tracks for geometry calculations.")
    ("material-tracks,N", po::value<int>(&mattracks)->default_value(100), "N. of tracks for material calculations.")
    ("power,p", "Report irradiated power analysis.")
    ("bandwidth,b", "Report base bandwidth analysis.")
    ("bandwidth-cpu,B", "Report multi-cpu bandwidth analysis.\n\t(implies 'b')")
    ("material,m", "Report materials and weights analyses.")
    ("resolution,r", "Report resolution analysis.")
    ("trigger,t", "Report base trigger analysis.")
    ("trigger-ext,T", "Report extended trigger analysis.\n\t(implies 't')")
    ("all,a", "Report all analyses, except extended\ntrigger. (implies all other relevant\nreport options)")
    ("graph,g", "Build and report neighbour graph.")
    ("xml", po::value<std::string>(&xmldir)->implicit_value(""), "Produce XML output files for materials.\nOptional arg specifies the subdirectory\nof the output directory (chosen via inst\nscript) where to create XML files.\nIf not supplied, basename will be used\nas subdir.")
    ("html-dir", po::value<std::string>(&htmldir), "Override the default html output dir\n(equal to the tracker name in the main\ncfg file) with the one specified.")
    ("verbosity", po::value<int>(&verbosity)->default_value(1), "Levels of details in the program's output (overridden by the option 'quiet').")
    ("quiet", "No output is produced, except the required messages (equivalent to verbosity 0, overrides the option 'verbosity')")
    ("performance", "Outputs the CPU time needed for each computing step (overrides the option 'quiet').")
    ("tracksim", po::value<std::vector<int> >(&tracksim)->multitoken(), "Simulates tracks")
    ;

  po::options_description hidden;
  hidden.add_options()("base-name", po::value<std::string>(&basename));

  po::positional_options_description posopt;
  posopt.add("base-name", 1); 

  po::options_description allopt;
  allopt.add(shown).add(hidden);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(allopt).positional(posopt).run(), vm);

    po::notify(vm);

    if (geomtracks < 1) throw po::invalid_option_value("geometry-tracks");
    if (mattracks < 1) throw po::invalid_option_value("material-tracks");
    if (!vm.count("base-name")) throw po::error("Missing required config basename argument"); 

  } catch(po::error e) {
    std::cerr << e.what() << std::endl << std::endl;
    std::cout << usage << std::endl << shown << std::endl;
    return EXIT_FAILURE;
  }

  if (vm.count("help")) {
    std::cout << usage << std::endl << shown << std::endl;
    return 0;
  }

  insur::Squid squid;
  bool verboseMaterial = false;
  unsigned int verboseWatch = verbosity;
  if (vm.count("quiet")) verboseWatch=0;
  bool performanceWatch = vm.count("performance");
  if (performanceWatch) {
    if (verboseWatch==0) verboseWatch = 1;
  }
  StopWatch::instance()->setVerbosity(verboseWatch, performanceWatch);

  squid.setBasename(basename);
  if (htmldir != "") squid.setHtmlDir(htmldir);

  // The tracker (and possibly pixel) must be build in any case
  if (!squid.buildTracker()) return EXIT_FAILURE;

  // The tracker should pick the types here but in case it does not,
  // we can still write something
  if (squid.dressTracker()) {
    if (!squid.pureAnalyzeGeometry(geomtracks)) return EXIT_FAILURE;

    if ((vm.count("all") || vm.count("bandwidth") || vm.count("bandwidth-cpu")) && !squid.reportBandwidthSite()) return EXIT_FAILURE;
    if ((vm.count("all") || vm.count("bandwidth-cpu")) && (!squid.reportTriggerProcessorsSite()) ) return EXIT_FAILURE;
    if ((vm.count("all") || vm.count("power")) && (!squid.irradiateTracker() || !squid.reportPowerSite()) ) return EXIT_FAILURE;

    // If we need to have the material model, then we build it
    if ( vm.count("all") || vm.count("material") || vm.count("resolution") || vm.count("graph") || vm.count("xml") ) {
      if (squid.buildInactiveSurfaces(verboseMaterial) && squid.createMaterialBudget(verboseMaterial)) {
        if ( vm.count("all") || vm.count("material") || vm.count("resolution") ) {
          // TODO: the following call should know whether to compute resolution or material (or both)
          if (!squid.pureAnalyzeMaterialBudget(mattracks, true)) return EXIT_FAILURE;
          if ((vm.count("all") || vm.count("material"))  && !squid.reportMaterialBudgetSite()) return EXIT_FAILURE;
          if ((vm.count("all") || vm.count("resolution"))  && !squid.reportResolutionSite()) return EXIT_FAILURE;	  
        }
        if (vm.count("graph") && !squid.reportNeighbourGraphSite()) return EXIT_FAILURE;
        if (vm.count("xml") && !squid.translateFullSystemToXML(xmldir.empty() ? basename : xmldir, false)) return (EXIT_FAILURE); //TODO: take care of flag in a more intelligent way...
      }
    }

    if ((vm.count("all") || vm.count("trigger") || vm.count("trigger-ext")) &&
        ( !squid.analyzeTriggerEfficiency(mattracks, vm.count("trigger-ext")) || !squid.reportTriggerPerformanceSite(vm.count("trigger-ext"))) ) return EXIT_FAILURE;
  } else if (!squid.pureAnalyzeGeometry(geomtracks)) return EXIT_FAILURE;


  if (!squid.reportGeometrySite()) return EXIT_FAILURE;
  if (!squid.additionalInfoSite()) return EXIT_FAILURE;
  if (!squid.makeSite()) return EXIT_FAILURE;

  if (vm.count("tracksim") && tracksim.size()==2) {
    std::cout << "Track sim mode enabled" << std::endl;
    squid.simulateTracks(tracksim[0], tracksim[1]);
  }

  return EXIT_SUCCESS;
}





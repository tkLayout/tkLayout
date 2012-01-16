
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

    std::string basename, xmlname;

    shown.add_options()
        ("help", "Display this help message.")
        ("geometry-tracks,n", po::value<int>(&geomtracks)->default_value(50), "N. of tracks for geometry calculations.")
        ("material-tracks,N", po::value<int>(&mattracks)->default_value(2000), "N. of tracks for material calculations.")
        ("power,p", "Report irradiated power analysis.")
        ("bandwidth,b", "Report bandwidth analysis.")
        ("material,m", "Report materials and weights analyses.")
        ("resolution,r", "Report resolution analysis.")
        ("trigger,t", "Report base trigger analysis.")
        ("trigger-ext,T", "Report extended trigger analysis.\n\t(implies 't')")
        ("all,a", "Report all analyses, except extended\ntrigger. (implies all other relevant\nreport options)")
        ("graph,g", "Build and report neighbour graph.")
        ("xml", po::value<std::string>(&xmlname)->implicit_value(""), "Produce XML output files for materials.\nOptional arg specifies the subdirectory\nof the output directory (chosen via inst\nscript) where to create XML files.\nIf not supplied, basename will be used\nas subdir.")
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
        if (!vm.count("base-name")) throw po::too_few_positional_options_error("Missing required config basename argument"); 

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

    std::string geomfile = basename + ".cfg";
    std::string settingsfile = basename + "_Types.cfg";
    std::string matfile = basename + "_Materials.cfg";
    std::string pixmatfile = basename + "_Materials.cfg.pix";

    if (!squid.buildFullSystem(geomfile, settingsfile, matfile, pixmatfile, true, true)) return EXIT_FAILURE; // CUIDADO: true, true might need to be replaced by proper boolean switches 

    squid.pureAnalyzeGeometry(mattracks);

    std::cout << "Calling analyzer with " << geomtracks << std::endl;

    if (!squid.pureAnalyzeMaterialBudget(geomtracks)) return EXIT_FAILURE;

    if (!squid.reportGeometrySite()) return EXIT_FAILURE;
    if ((vm.count("all") || vm.count("bandwidth")) && !squid.reportBandwidthSite()) return EXIT_FAILURE;
    if ((vm.count("all") || vm.count("power"))     && !squid.reportPowerSite()) return EXIT_FAILURE;
    if ((vm.count("all") || vm.count("material"))  && !squid.reportMaterialBudgetSite()) return EXIT_FAILURE;
    if ((vm.count("all") || vm.count("trigger") || vm.count("trigger-ext")) && !squid.reportTriggerPerformanceSite(vm.count("trigger-ext"))) return EXIT_FAILURE;
    if (vm.count("graph") && !squid.reportNeighbourGraphSite()) return EXIT_FAILURE;
    if (!squid.additionalInfoSite(geomfile, settingsfile, matfile, pixmatfile)) return EXIT_FAILURE;
    if (!squid.makeSite()) return EXIT_FAILURE;
    
    if (vm.count("xml") && !squid.translateFullSystemToXML(xmlname.empty() ? basename : xmlname, false)) return (EXIT_FAILURE); //TODO: take care of flag in a more intelligent way...

    return EXIT_SUCCESS;
}





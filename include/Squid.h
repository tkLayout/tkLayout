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
#include <tracker.hh>
#include <configparser.hh>
#include <MatParser.h>
#include <InactiveSurfaces.h>
#include <MaterialBudget.h>
#include <Usher.h>
#include <MatCalc.h>
#include <Analyzer.h>
#include <Vizard.h>
#include <tk2CMSSW.h>
#include <boost/filesystem/operations.hpp>
#include <rootweb.hh>
#include <mainConfigHandler.h>

/**
 * A shorter alias for the filesystem library namespace
 */
namespace bfs = boost::filesystem;
namespace insur {
    /*
     * Error messages and warnings that may be reported.
     */
    static const std::string err_no_geomfile = "Error: there is no recorded name for the geometry configuration file. Initialise the tracker first.";
    static const std::string err_no_matfile = "Error: the provided material configuration file does not exist.";
    static const std::string err_init_failed = "Error: initialisation of the material calculator failed.";
    static const std::string err_no_tracker = "Error: the tracker object does not exist. The tracker must be created before calling this function.";
    static const std::string err_no_inacsurf = "Error: the collection of inactive surfaces does not exist. It must be created before calling this function";
    static const std::string err_no_matbudget = "Error: the material budget does not exist. It must be created before calling this function.";
    static const std::string warning_rootonly = "Warning: the collection of inactive surfaces does not exist. Only the .root file will be written.";
    static const std::string default_trackername = "defaultTrackerName";

    
    /**
     * @class Squid
     * @brief The Squid class integrates the components of the <i>tkmaterial</i> application and provides an
     * interface to its high-level functionality.
     *
     * It manages instances to all the necessary components internally. Its functions bundle the steps that are
     * required to carry out the main tasks of the application. The idea is to allow access to the underlying
     * via a series of well-defined pathways contained in a single class. This class acts as a sort of manager 
     * to the application and can be embedded into a main program easily.
     */
    class Squid {
    public:
        Squid();
        virtual ~Squid();
        bool buildTracker(std::string geomfile);
        bool dressTracker(std::string settingsfile);
        bool buildTrackerSystem(std::string geomfile, std::string settingsfile);
        bool buildInactiveSurfaces(bool verbose = false);
        bool buildInactiveSurfaces(std::string geomfile, bool verbose = false);
        bool buildInactiveSurfaces(std::string geomfile, std::string settingsfile, bool verbose = false);
        bool createMaterialBudget(std::string matfile, bool verbose = false);
        bool buildFullSystem(std::string geomfile, std::string settingsfile, std::string matfile, bool usher_verbose = false, bool mat_verbose = false);
        bool analyzeFullSystem(std::string htmlout = "", std::string rootout = "", std::string graphout = "", int tracks = 50, bool simplified = true);
        bool analyzeGeoMat(std::string htmlout = "", std::string rootout = "", int tracks = 50, bool simplified = true);
        bool analyzeGraphMat(std::string htmlout = "", std::string graphout = "", int tracks = 50);
        bool analyzeGeometry(std::string rootout = "", std::string graphout = "", bool simplified = true);
        bool analyzeGeometry(std::string rootout = "", bool simplified = true);
        bool analyzeNeighbours(std::string graphout = "");
        bool analyzeMaterialBudget(std::string htmlout = "", int tracks = 50);
        bool translateFullSystemToXML(std::string xmlout = "");
        bool trackerSummary(std::string configFileName, std::string dressFileName);
#ifdef USING_ROOTWEB
	// Functions using rootweb
	bool analyzeGeometrySite(int tracks = 1000);
	bool analyzeMaterialBudgetSite(int tracks = 50);
	bool additionalInfoSite(std::string& geomfile, std::string& settingsfile, std::string& matfile);
	bool makeSite();
#endif
    private:
        std::string g;
        Tracker* tr;
        InactiveSurfaces* is;
        MaterialBudget* mb;
        configParser cp;
        MatParser mp;
        Usher u;
        MatCalc c;
        Analyzer a;
        Vizard v;
        tk2CMSSW t2c;
        mainConfigHandler mainConfiguration;
        bool fileExists(std::string filename);
        std::string extractFileName(const std::string& full);
        Squid(const Squid& s);
        Squid& operator=(const Squid& s);
#ifdef USING_ROOTWEB
	RootWSite site;
	bool prepareWebsite();
	bool sitePrepared;
#endif
    };
}
#endif	/* _SQUID_H */


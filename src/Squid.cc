/**
 * void littleTestBench();
 * @file Squid.cc
 * @brief This implements the main interface between the tkgeometry library classes and the frontend
 */

#include <SvnRevision.h>
#include "Squid.h"
#include "AnalyzerOccupancy.h"
#include "StopWatch.h"
#include <boost/algorithm/string/split.hpp>

#include <AnalyzerGeometry.h>

namespace insur {

  // public
  const std::string Squid::err_no_geomfile           = "There is no recorded name for the geometry configuration file. Initialise the tracker first.";
  const std::string Squid::err_no_matfile            = "The provided material configuration file does not exist.";
  const std::string Squid::err_no_matfile_pixel      = "The material configuration file for the pixels does not exist.";
  const std::string Squid::err_init_failed           = "Initialisation failed.";
  const std::string Squid::err_no_tracker            = "The tracker object does not exist. The tracker must be created before calling this function.";
  const std::string Squid::err_no_inacsurf           = "The collection of inactive surfaces does not exist. It must be created before calling this function";
  const std::string Squid::err_no_matbudget          = "The material budget does not exist. It must be created before calling this function.";
  const std::string Squid::err_no_triggerSummary     = "Could not report on the trigger performance.";
  const std::string Squid::err_no_flukafile          = "The fluka charged hadrons or photons flux file doesn't exist.";
  const std::string Squid::err_no_occupancy_failed   = "Occupancy calculation or visualisation failed.";
  const std::string Squid::warn_rootonly             = "The collection of inactive surfaces does not exist. Only the .root file will be written.";
  const std::string Squid::warn_custom_matfile       = "A customized material file was used for the tracker";
  const std::string Squid::warn_custom_matfile_pixel = "A customized material file was used for the pixel";

  const std::string Squid::default_trackername       = "defaultTrackerName";

  /**
   * The constructor sets the internal pointers to <i>NULL</i>.
   */
  Squid::Squid() :
      t2c(mainConfigHandler::instance()),
      m_weightDistCtrlStd(0.1),
      m_weightDistFwdStd(0.1),
      m_weightDistCtrlPxd(0.1),
      m_weightDistFwdPxd(0.1) {

    m_ctrlPxd        = nullptr;
    m_ctrlStd        = nullptr;
    m_fwdPxd         = nullptr;
    m_fwdStd         = nullptr;

    m_pasiveCtrlPxd  = nullptr;
    m_pasiveFwdPxd   = nullptr;
    m_pasiveCtrlStd  = nullptr;
    m_pasiveFwdStd   = nullptr;

    m_mbBP           = nullptr;
    m_mbCtrlPxd      = nullptr;
    m_mbCtrlStd      = nullptr;
    m_mbFwdPxd       = nullptr;
    m_mbFwdStd       = nullptr;

    m_geomAnalyzer   = nullptr;

    //m_simParms       = nullptr;
    m_sitePrepared   = false;
    m_myGeometryFile = "";
    m_mySettingsFile = "";

    m_htmlDir        = "";
  }

  /**
   * The destructor deletes the heap-allocated internal objects if they exist.
   */
  Squid::~Squid() {

    if (m_ctrlPxd) delete m_ctrlPxd;
    if (m_ctrlStd) delete m_ctrlStd;
    if (m_fwdPxd)  delete m_fwdPxd;
    if (m_fwdStd)  delete m_fwdStd;

    if (m_pasiveCtrlPxd) delete m_pasiveCtrlPxd;
    if (m_pasiveCtrlStd) delete m_pasiveCtrlStd;

    if (m_mbBP)      delete m_mbBP;
    if (m_mbCtrlPxd) delete m_mbCtrlPxd;
    if (m_mbCtrlStd) delete m_mbCtrlStd;
    if (m_mbFwdPxd)  delete m_mbFwdPxd;
    if (m_mbFwdStd)  delete m_mbFwdStd;

    //if (m_simParms)  delete m_simParms;

    StopWatch::destroy();
  }

  /**
   * Build a bare-bones geometry of active modules. The resulting tracker object replaces the previously
   * registered one, if such an object existed. It remains in the squid until it is overwritten by a second call
   * to this function or by another one that creates a new tracker object. If the geometry file contains a
   * description of a pixel detector, that is created as well. It is treated exactly like the tracker with respect
   * to replacement and persistence.
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::buildActiveTracker() {

    if (m_ctrlPxd) {
      delete m_ctrlPxd;
      m_ctrlPxd = nullptr;
    }
    if (m_ctrlStd) {
      delete m_ctrlStd;
      m_ctrlStd = nullptr;
    }
    if (m_fwdPxd) {
      delete m_fwdPxd;
      m_fwdPxd = nullptr;
    }
    if (m_fwdStd) {
      delete m_fwdStd;
      m_fwdStd = nullptr;
    }

    // Geometry/config file not defined
    std::ifstream ifs(getGeometryFile());
    if (ifs.fail()) {
      std::cerr << "ERROR: cannot open geometry file " << getGeometryFile() << std::endl;
      return false;
    }

    // Build tracker
    startTaskClock("Building full tracker");

    // Parse config file & save everything in mainConfig
    std::stringstream ss;
    m_includeSet = mainConfigHandler::instance().preprocessConfiguration(ifs, ss, getGeometryFile());
    t2c.addConfigFile(tk2CMSSW::ConfigFile{getGeometryFile(), ss.str()});

    using namespace boost::property_tree;
    ptree pTree;
    info_parser::read_info(ss, pTree);

    /*
    class CoordExportVisitor : public ConstGeometryVisitor {
      std::ofstream barof, endof;
    public:
      CoordExportVisitor(std::string trid) : barof(trid + "_barrel_coords.txt"), endof(trid + "_endcap_coords.txt") {}
      void visit(const BarrelModule& m) override { barof << m.center().Z() << ", " << m.center().Rho() << ", " << m.center().Phi() << std::endl; }
      void visit(const EndcapModule& m) override { endof << m.center().Z() << ", " << m.center().Rho() << ", " << m.center().Phi() << std::endl; }
    };

    class ModuleDataVisitor : public ConstGeometryVisitor {
      std::ofstream of;
      const char sep = '\t';
    public:
      ModuleDataVisitor(std::string ofname) : of(ofname + "_mods.txt") {
         of << "cntName" << sep << "refZ" << sep << "refRho" << sep << "refPhi" << sep
            << "centerZ" << sep << "centerRho" << sep << "centerPhi" << sep
            << "dsDist" << sep << "thickn" << sep 
            << "minW" << sep << "maxW" << sep << "len" << sep
            << "type" << sep 
            << "resolRPhi" << sep << "resolY" << sep
            << "nStripA" << sep << "nSegmI" << sep << "nSegmO" << std::endl;
      }
      void visit(const DetectorModule& m) override {
        if (m.minZ() < 0.) return; // || m.posRef().phi != 1) return;
        of << m.cntName() << sep << (int)m.posRef().z << sep << (int)m.posRef().rho << sep << (int)m.posRef().phi << sep
           << m.center().Z() << sep << m.center().Rho() << sep << m.center().Phi() << sep
           << m.dsDistance() << sep << m.thickness() << sep 
           << m.minWidth() << sep << m.maxWidth() << sep << m.length() << sep
           << m.moduleType() << sep
           << m.resolutionLocalX() << sep << m.resolutionLocalY() << sep 
           << m.numStripsAcross() << sep << m.innerSensor().numSegments() << sep << m.outerSensor().numSegments() << std::endl;
      }
    };
    */

    try { 

      // Look for tag "Tracker" and build tracker
      auto childRange = getChildRange(pTree, "Tracker");
      std::for_each(childRange.first, childRange.second, [&](const ptree::value_type& kv) {

        Tracker* trk = new Tracker(kv.second);
        trk->build();

        // Distinguish between individual trackers
        if      (trk->myid()=="Inner"    || trk->myid()=="CtrlInner" ) {
          trk->setIsPixelType(true);
          m_trackers.push_back(trk);
          m_ctrlPxd    = trk;
        }
        else if (trk->myid()=="FwdInner" || trk->myid()=="FwdInner" || trk->myid()=="ForwardPixel") {
          trk->setIsPixelType(true);
          trk->setIsForwardType(true);
          m_trackers.push_back(trk);
          m_fwdPxd = trk;
        }
        else if (trk->myid()=="Outer"    || trk->myid()=="CtrlOuter" ) {
          trk->setIsStripType(true);
          m_trackers.push_back(trk);
          m_ctrlStd    = trk;
        }
        else if (trk->myid()=="FwdOuter" || trk->myid()=="FwdOuter" || trk->myid()=="ForwardStrip" ) {
          trk->setIsStripType(true);
          trk->setIsForwardType(true);
          m_trackers.push_back(trk);
          m_fwdStd = trk;
        }
        else {

          if (m_ctrlStd==nullptr) {

            std::string message = std::string("Squid::buildTracker: The following tracker "+trk->myid()+" set as default strip tracker! Its specific name not recognized!");
            logWARNING(message);
            m_ctrlStd = trk; // General sollution
          }
          else {

            std::string message = std::string("Squid::buildTracker: One tracker already assigned as a strip tracker! This one: "+trk->myid()+" will be thus omitted! Its specific name not recognized!");
            logWARNING(message);
            delete trk;
          }
        }
      });

      // Declare unmatched properties
      std::set<string> unmatchedProperties = PropertyObject::reportUnmatchedProperties();
      if (!unmatchedProperties.empty()) {

        std::ostringstream ss;
        ss << "The following unknown properties were ignored:" << std::endl;

        for (const string& s : unmatchedProperties) {
          ss << "  " << s << std::endl;
        }
        logERROR(ss);
      }

      // Look for tag "Support" not associated with a concrete Tracker and build supports
      childRange = getChildRange(pTree, "Support");
      std::for_each(childRange.first, childRange.second, [&](const ptree::value_type& kv) {

        Support* support = new Support();
        support->myid(kv.second.get_value(0));
        support->store(kv.second);
        support->build();
        m_supports.push_back(support);
      });

      // Create instance of simulation parameters
      SimParms* simParms = SimParms::getInstance();

      // Store data in a general tree & check
      simParms->store(getChild(pTree, "SimParms"));
      simParms->crosscheck();
    }
    catch (PathfulException& e) { 
      std::cerr << e.path() << " : " << e.what() << std::endl; 
      stopTaskClock();
      return false;
    }

    // Check that some tracker created
    if (m_trackers.size()==0) {

      std::cerr << "ERROR: Detector can't be built!!! No trackers defined ..." << std::endl;
      stopTaskClock();
      return false;
    }

    stopTaskClock();
    return true;
  }

  /*
   *  Build all pasive components of tracker & create database of overall material distribution in the tracker ...
   */
  bool Squid::buildPasiveTracker(bool verbose) {

    if (m_ctrlStd || m_fwdStd || m_ctrlPxd || m_fwdPxd) {

      startTaskClock("Building inactive materials of defined trackers + beam pipe");

      // Building pxd & strip inactive materials
      if (m_ctrlPxd) {
        if (!m_pasiveCtrlPxd) m_pasiveCtrlPxd = new InactiveSurfaces();
        Materialway materialwayPixel;
        materialwayPixel.build(*m_ctrlPxd, *m_pasiveCtrlPxd, m_weightDistCtrlPxd);
      }
      if (m_fwdPxd) {
        if (!m_pasiveFwdPxd) m_pasiveFwdPxd = new InactiveSurfaces();
         Materialway materialwayFwdPixel;
         materialwayFwdPixel.build(*m_fwdPxd, *m_pasiveFwdPxd, m_weightDistFwdPxd);
      }
      if (m_ctrlStd) {
        if (!m_pasiveCtrlStd) m_pasiveCtrlStd = new InactiveSurfaces();
        Materialway materialwayStrip;
        materialwayStrip.build(*m_ctrlStd, *m_pasiveCtrlStd, m_weightDistCtrlStd);
      }
      if (m_fwdStd) {
        if (!m_pasiveFwdStd) m_pasiveFwdStd = new InactiveSurfaces();
        Materialway materialwayFwdStrip;
        materialwayFwdStrip.build(*m_fwdStd, *m_pasiveFwdStd, m_weightDistFwdStd);
      }
    } else {

      std::string message=std::string("Squid::buildMaterials() :")+err_no_tracker;
      logERROR(message);
      stopTaskClock();
      return false;
    }

    stopTaskClock();
    return true;
  }

  /**
   * Calculate a material budget for the previously created tracker object and collection of inactive
   * surfaces. The resulting material budget replaces the previously registered one, if such an object
   * existed. It remains in the squid until it is overwritten by a second call to this function or by
   * another one that creates a new material budget. This function succeeds if either the tracker or
   * both the tracker and the inactive surfaces exist. If a pixel detector exists, it gets its own material
   * budget at the same time. This object is treated exactly the same as the one for the tracker.
   * @param verbose A flag that turns the final status summary of the material budget on or off
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::createMaterialBudget(bool verbose) {

    if (m_ctrlStd || m_fwdStd || m_ctrlPxd || m_fwdPxd) {

      startTaskClock("Calculating material budget of the tracker");
      if (m_ctrlStd && m_pasiveCtrlStd) {

        if (m_mbCtrlStd) delete m_mbCtrlStd;
        m_mbCtrlStd  = new MaterialBudget(*m_ctrlStd, *m_pasiveCtrlStd);
      }
      if (m_fwdStd && m_pasiveFwdStd) {

        if (m_mbFwdStd) delete m_mbFwdStd;
        m_mbFwdStd  = new MaterialBudget(*m_fwdStd, *m_pasiveFwdStd);
      }
      if (m_ctrlPxd && m_pasiveCtrlPxd) {

        if (m_mbCtrlPxd) delete m_mbCtrlPxd;
        m_mbCtrlPxd = new MaterialBudget(*m_ctrlPxd, *m_pasiveCtrlPxd);
      }
      if (m_fwdPxd && m_pasiveFwdPxd) {

        if (m_mbFwdPxd) delete m_mbFwdPxd;
        m_mbFwdPxd = new MaterialBudget(*m_fwdPxd, *m_pasiveFwdPxd);
      }
      stopTaskClock();
      return true;
    } else {

      std::string message = std::string("Squid::createMaterialBudget() :")+err_init_failed;
      logERROR(message);
      stopTaskClock();
      return false;
    }
  }

  /**
   * Translate an existing full tracker and material budget to a series of XML files that can be interpreted by CMSSW.
   * @param xmlout The name - without path - of the designated output subdirectory
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::translateFullSystemToXML(std::string xmlout) {

    MatCalc   matCalc;
    MatParser matParser;

    if (m_mbCtrlStd) {

      if (matCalc.initDone()) matCalc.reset();
      if (matParser.initMatCalc(matCalc, mainConfigHandler::instance().getMattabDirectory())) {
        m_mbCtrlStd->materialsAll(matCalc);
      }
      t2c.translate(matCalc.getMaterialTable(), *m_mbCtrlStd, xmlout.empty() ? m_baseName : xmlout, false); // false is setting a mysterious flag called wt which changes the way the XML is output. apparently setting it to true is of no use anymore.
      return true;
    }
    else {
      std::cout << "Squid::translateFullSystemToXML(): " << err_no_matbudget << std::endl;
      return false;
    }
  }

  // private
  /**
   * Check if a given configuration file actually exists.
   * @param filename The name of the configuration file
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::fileExists(std::string filename) {
    bfs::path p(filename);
    if (bfs::exists(p)) return true;
    return false;
  }

  /**
   * Extract the filename from a path. If there is nothing to extract, a copy of the input is returned.
   * @full The source string
   * @return A string containing everything after the last slash character, or the entire input if there were no slashes
   */
  std::string Squid::extractFileName(const std::string& full) {
    std::string::size_type idx = full.find_last_of("/");
    if (idx != std::string::npos)
      return full.substr(idx+1);
    else return full;
  }

  // private
  void Squid::resetVizard() {

    m_vizard.~Vizard();
    new ((void*) &m_vizard) Vizard();
  }

  /**
   * Prepare the website object (if not done yet) from the configuration file
   * it needs the tracker object to be already there
   * @return a boolean with the operation success
   */
  bool Squid::prepareWebSite() {

    if (m_sitePrepared) return true;

    if (m_htmlDir!="") m_webSite.setTargetDirectory(m_htmlDir);
    else return false;

    // Set layout title
    if (m_layoutName!="") m_webSite.setTitle(m_layoutName);
    else return false;

    m_webSite.setComment("Layouts");
    m_webSite.setCommentLink("../");
    m_webSite.addAuthor("Giovanni Bianchi");
    m_webSite.addAuthor("Nicoletta De Maio");
    m_webSite.addAuthor("Stefano Martina");
    m_webSite.addAuthor("Stefano Mersi");
    m_webSite.setRevision(SvnRevision::revisionNumber);
    return true;
  }

  /**
   * Actually creates the website where it was supposed to be
   * @return a boolean with the operation success
   */
  bool Squid::makeWebSite(bool addLogPage /* = true */) {

    ostringstream message;
    message << "Creating website in: " << m_htmlDir;

    startTaskClock(message.str());
    if (!prepareWebSite()) {

      logERROR("Problem in preparing website");
      stopTaskClock();
      return false;
    }

    // Add log webPage
    if (addLogPage) {
      m_vizard.makeLogPage(m_webSite);
    }

    bool result = m_webSite.makeSite(false);
    stopTaskClock();
    return result;
  }

  /**
   * Analyze tracker layout geometry with no output through rootweb.
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::analyzeGeometry(int nTracks) {

//    startTaskClock("Analyzing geometry");
//
//    if (m_trackers.size()>0) {
//
//      m_geomAnalyzer = new AnalyzerGeometry(m_trackers);
//
//      if (m_geomAnalyzer->init(nTracks)) {
//
//        bool isOK = m_geomAnalyzer->analyze();
//        stopTaskClock();
//        return isOK;
//      }
//      else {
//
//        std::string message = std::string("Squid::analyzeGeometry(): ")+err_init_failed;
//        logERROR(message);
//        return false;
//      }
//    }
//    else {
//
//      std::string message = std::string("Squid::analyzeGeometry(): ")+err_no_tracker;
//      logERROR(message);
//      return false;
//    }
  }

//  bool Squid::analyzeTriggerEfficiency(int tracks, bool detailed) {
//    // Call this before analyzetrigger if you want to have the map of suggested spacings
//    if (detailed) {
//      startTaskClock("Creating distance tuning plots");
//      m_stripAnalyzer.createTriggerDistanceTuningPlots(*m_ctrlStd, mainConfigHandler::instance().getTriggerMomenta());
//      stopTaskClock();
//    }
//    startTaskClock("Creating trigger efficiency plots");
//    m_stripAnalyzer.analyzeTriggerEfficiency(*m_ctrlStd,
//                               mainConfigHandler::instance().getTriggerMomenta(),
//                               mainConfigHandler::instance().getThresholdProbabilities(),
//                               tracks);
//    stopTaskClock();
//    return true;
//  }

  /**
   * Analyze the previously created material budget with no output.
   * @param tracks The number of tracks that should be fanned out across the analysed region
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::pureAnalyzeMaterialBudget(int nTracks) {

    if (m_mbCtrlStd || m_mbFwdStd || m_mbCtrlPxd || m_mbFwdPxd) {
      // TODO: insert the creation of sample tracks here, to compute intersections only once

      startTaskClock("Analyzing material budget" );
      std::vector<MaterialBudget*> trkMb;
      if (m_mbCtrlPxd) {
        trkMb.clear();
        trkMb.push_back(m_mbCtrlPxd);
        m_pixelAnalyzer.analyzeMaterialBudget(trkMb, nTracks);
      }
      if (m_mbCtrlStd) {
        trkMb.clear();
        trkMb.push_back(m_mbCtrlStd);
        m_stripAnalyzer.analyzeMaterialBudget(trkMb, nTracks);
      }
      if (m_mbFwdPxd) {
        trkMb.clear();
        if (m_mbFwdPxd) trkMb.push_back(m_mbFwdPxd);
        //if (m_mbFwdStd) trkMb.push_back(m_mbFwdStd);
        m_fwdPxdAnalyzer.analyzeMaterialBudget(trkMb, nTracks);
      }
      if (m_mbFwdStd) { //(m_mbFwdPxd || m_mbFwdStd) {
        trkMb.clear();
        //if (m_mbFwdPxd) trkMb.push_back(m_mbFwdPxd);
        if (m_mbFwdStd) trkMb.push_back(m_mbFwdStd);
        m_fwdAnalyzer.analyzeMaterialBudget(trkMb, nTracks);
      }
      if (m_mbCtrlPxd || m_mbFwdPxd || m_mbCtrlStd || m_mbFwdStd) {
        trkMb.clear();
        if (m_mbCtrlPxd) trkMb.push_back(m_mbCtrlPxd);
        if (m_mbFwdPxd)  trkMb.push_back(m_mbFwdPxd);
        if (m_mbCtrlStd) trkMb.push_back(m_mbCtrlStd);
        if (m_mbFwdStd)  trkMb.push_back(m_mbFwdStd);
        m_trackerAnalyzer.analyzeMaterialBudget(trkMb, nTracks);
      }
      trkMb.clear();
      stopTaskClock();

      //startTaskClock("Computing the weight summary");
      //m_stripAnalyzer.computeWeightSummary(*m_mbCtrlStd);
      //stopTaskClock();

      return true;
    } else {

      std::string message = std::string("Squid::pureAnalyzeMaterialBudget() :")+err_no_matbudget;
      logERROR(message);
      return false;
    }
  }

  /**
   * Analyze the tracker resolution with no output.
   * @param tracks The number of tracks that should be fanned out across the analysed region
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::pureAnalyzeResolution(int nTracks) {

    if (m_mbCtrlStd || m_mbFwdStd || m_mbCtrlPxd || m_mbFwdPxd) {
    // TODO: insert the creation of sample tracks here, to compute intersections only once

      startTaskClock("Analyzing tracking resolutions");
      std::vector<MaterialBudget*> trkMb;
      trkMb.clear();
      if (m_mbCtrlPxd) trkMb.push_back(m_mbCtrlPxd);
      if (m_mbFwdPxd)  trkMb.push_back(m_mbFwdPxd);
      if (m_mbCtrlStd) trkMb.push_back(m_mbCtrlStd);
      if (m_mbFwdStd)  trkMb.push_back(m_mbFwdStd);
      m_trackerAnalyzer.analyzeTaggedTracking(trkMb, mainConfigHandler::instance().getMomenta(), nTracks);

      trkMb.clear();
      stopTaskClock();
      return true;
    } else {

      std::string message = std::string("Squid::pureAnalyzeResolution() :")+err_no_matbudget;
      logERROR(message);
      return false;
    }
  }

  /**
   * Produces the output of the analysis of the geomerty
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::reportGeometrySite() {

    // Report geometry layout
    return m_geomAnalyzer->visualize(m_webSite);
  }

  bool Squid::reportBandwidthSite() {
    if (m_ctrlStd) {
      startTaskClock("Computing bandwidth and rates");
      m_stripAnalyzer.computeBandwidth(*m_ctrlStd);
      m_stripAnalyzer.computeTriggerFrequency(*m_ctrlStd);
      stopTaskClock();
      startTaskClock("Creating bandwidth and rates report");
      m_vizard.bandwidthSummary(m_stripAnalyzer, *m_ctrlStd, m_webSite);
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }
  }

  bool Squid::reportTriggerProcessorsSite() {
    if (m_ctrlStd) {
      startTaskClock("Computing multiple trigger tower connections");
      m_stripAnalyzer.computeTriggerProcessorsBandwidth(*m_ctrlStd);
      m_vizard.triggerProcessorsSummary(m_stripAnalyzer, *m_ctrlStd, m_webSite);
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }
  }
  /*
   * Calculate expected detector occupancy using Fluka simulation outputs:
   *  - charged particles map
   *  - photons map
   */
  bool Squid::reportOccupancySite() {
    if (m_trackers.size()!=0) {
      if (SimParms::getInstance()->chargedMapFile()!="" && SimParms::getInstance()->photonsMapFile()!="") {

        std::ostringstream message;
        startTaskClock("Computing occupancy");

        // Create occupancy analyzer & analyzer data
        std::string chargedMapFile          = mainConfigHandler::instance().getIrradiationDirectory() + "/" + SimParms::getInstance()->chargedMapFile();
        std::string chargedMapLowThFile     = SimParms::getInstance()->chargedMapLowThFile();
        std::string chargedMapBOffMatOnFile = SimParms::getInstance()->chargedMapBOffMatOnFile();
        std::string chargedMapBOnMatOffFile = SimParms::getInstance()->chargedMapBOnMatOffFile();
        std::string chargedMapBOffMatOffFile= SimParms::getInstance()->chargedMapBOffMatOffFile();
        std::string chargedMapBOffTrkOffFile= SimParms::getInstance()->chargedMapBOffTrkOffFile();
        std::string photonsMapFile          = mainConfigHandler::instance().getIrradiationDirectory() + "/" + SimParms::getInstance()->photonsMapFile();
        std::string photonsMapLowThFile     = SimParms::getInstance()->photonsMapLowThFile();
        std::string photonsMapBOffMatOnFile = SimParms::getInstance()->photonsMapBOffMatOnFile();
        std::string photonsMapBOnMatOffFile = SimParms::getInstance()->photonsMapBOnMatOffFile();
        std::string photonsMapBOffMatOffFile= SimParms::getInstance()->photonsMapBOffMatOffFile();
        std::string photonsMapBOffTrkOffFile= SimParms::getInstance()->photonsMapBOffTrkOffFile();

        std::string bFieldMapFile           = SimParms::getInstance()->bFieldMapFile();

//        AnalyzerOccupancy analyzerOccupancy(chargedMapFile, photonsMapFile, m_trackers);
//
//        // File existence is tested within the Analyzer class
//        if (!analyzerOccupancy.readMagFieldMap(mainConfigHandler::instance().getIrradiationDirectory(), bFieldMapFile)) {
//          logINFO("Couldn't read the mag. field map!");
//          logINFO(std::string(mainConfigHandler::instance().getIrradiationDirectory()+"/"+bFieldMapFile));
//        }
//        if (!analyzerOccupancy.readNoMagFieldIrradMap(mainConfigHandler::instance().getIrradiationDirectory(),chargedMapBOffMatOnFile, photonsMapBOffMatOnFile)) {
//          logINFO("Couldn't read the irradiation map with no mag. field!");
//          logINFO(std::string(mainConfigHandler::instance().getIrradiationDirectory()+"/"+chargedMapBOffMatOnFile+ " or "+mainConfigHandler::instance().getIrradiationDirectory()+"/"+photonsMapBOffMatOnFile));
//        }
//        if (!analyzerOccupancy.readNoMaterialIrradMap(mainConfigHandler::instance().getIrradiationDirectory(),chargedMapBOnMatOffFile, photonsMapBOnMatOffFile)) {
//          logINFO("Couldn't read the irradiation map with no material!");
//          logINFO(std::string(mainConfigHandler::instance().getIrradiationDirectory()+"/"+chargedMapBOnMatOffFile+ " or "+mainConfigHandler::instance().getIrradiationDirectory()+"/"+photonsMapBOnMatOffFile));
//        }
//        if (!analyzerOccupancy.readNoMagFieldNoMaterialIrradMap(mainConfigHandler::instance().getIrradiationDirectory(), chargedMapBOffMatOffFile, photonsMapBOffMatOffFile)) {
//          logINFO("Couldn't read the irradiation map with no mag. field & no material!");
//          logINFO(std::string(mainConfigHandler::instance().getIrradiationDirectory()+"/"+chargedMapBOffMatOffFile+ " or "+mainConfigHandler::instance().getIrradiationDirectory()+"/"+photonsMapBOffMatOffFile));
//        }
//        if (!analyzerOccupancy.readNoMagFieldNoTrackerIrradMap(mainConfigHandler::instance().getIrradiationDirectory(), chargedMapBOffTrkOffFile, photonsMapBOffTrkOffFile)) {
//          logINFO("Couldn't read the irradiation map with no mag. field & no tracker material!");
//          logINFO(std::string(mainConfigHandler::instance().getIrradiationDirectory()+"/"+chargedMapBOffTrkOffFile+ " or "+mainConfigHandler::instance().getIrradiationDirectory()+"/"+photonsMapBOffTrkOffFile));
//        }
//        if (!analyzerOccupancy.readLowThresholdIrradMap(mainConfigHandler::instance().getIrradiationDirectory(), chargedMapLowThFile, photonsMapLowThFile)) {
//          logINFO("Couldn't read the irradiation map with low thresholds set!");
//          logINFO(std::string(mainConfigHandler::instance().getIrradiationDirectory()+"/"+chargedMapLowThFile+ " or "+mainConfigHandler::instance().getIrradiationDirectory()+"/"+photonsMapLowThFile));
//        }
//
//        bool outCalc = analyzerOccupancy.analyze();
//        bool outVis  = analyzerOccupancy.visualize(m_webSite);

        stopTaskClock();
//        if (outCalc && outVis) return true;
//        else {
//          logERROR(err_no_flukafile);
//          return false;
//        }
        return true;
      }
      else {
        logERROR(err_no_flukafile);
        return false;
      }
    } else {
      logERROR(err_no_tracker);
      return false;
    }
  }

  bool Squid::reportPowerSite() {
    if (m_ctrlStd) {
      startTaskClock("Computing dissipated power");
      m_stripAnalyzer.analyzePower(*m_ctrlStd);
      stopTaskClock();
      startTaskClock("Creating power report");
      m_vizard.irradiatedPowerSummary(m_stripAnalyzer, *m_ctrlStd, m_webSite);
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }
  }

  /**
   * Produces the output of the analysis of the material budget analysis
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::reportMaterialBudgetSite(bool debugServices) {

    if (m_mbCtrlStd || m_mbCtrlPxd) {

      startTaskClock("Creating material budget report");
      std::vector<MaterialBudget*> materialBudgets;

      if (m_mbCtrlPxd) {

        materialBudgets.clear();
        materialBudgets.push_back(m_mbCtrlPxd);
        m_vizard.materialSummary(m_pixelAnalyzer, materialBudgets, false, m_webSite, "INNER");
      }
      if (m_mbCtrlStd) {

        materialBudgets.clear();
        materialBudgets.push_back(m_mbCtrlStd);
        m_vizard.materialSummary(m_stripAnalyzer, materialBudgets, false, m_webSite, "OUTER");
      }
      if (m_mbFwdStd) { //(m_mbFwdPxd || m_mbFwdStd) {

        materialBudgets.clear();
        //if (m_mbFwdPxd) materialBudgets.push_back(m_mbFwdPxd);
        if (m_mbFwdStd) materialBudgets.push_back(m_mbFwdStd);
        m_vizard.materialSummary(m_fwdAnalyzer, materialBudgets, false, m_webSite, "FWD");
      }
      if (m_mbFwdPxd) {

        materialBudgets.clear();
        if (m_mbFwdPxd) materialBudgets.push_back(m_mbFwdPxd);
        m_vizard.materialSummary(m_fwdPxdAnalyzer, materialBudgets, false, m_webSite, "FWDIN");
      }
      if (m_mbCtrlPxd || m_mbFwdPxd || m_mbCtrlStd || m_mbFwdStd) {

        materialBudgets.clear();
        if (m_mbCtrlPxd) materialBudgets.push_back(m_mbCtrlPxd);
        if (m_mbFwdPxd)  materialBudgets.push_back(m_mbFwdPxd);
        if (m_mbCtrlStd) materialBudgets.push_back(m_mbCtrlStd);
        if (m_mbFwdStd)  materialBudgets.push_back(m_mbFwdStd);
        m_vizard.materialSummary(m_trackerAnalyzer, materialBudgets, debugServices, m_webSite, "TRK");
      }
      //if (m_mbCtrlPxd) m_vizard.weigthSummart(m_pixelAnalyzer, weightDistributionPixel  , m_webSite, "INNER");
      //if (m_mbCtrlStd) m_vizard.weigthSummart(m_stripAnalyzer, weightDistributionTracker, m_webSite, "OUTER");

      materialBudgets.clear();
      stopTaskClock();
      return true;
    }
    else {
      logERROR(err_no_matbudget);
      return false;
    }
  }

  /**
   * Produces the output of the resolution measurement
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::reportResolutionSite() {

    if (m_ctrlStd || m_ctrlPxd) {

      startTaskClock("Creating resolution report");
      m_vizard.taggedErrorSummary(m_trackerAnalyzer, m_webSite);
      stopTaskClock();
      return true;
    }
    else {

      std::string message = std::string("Squid::reportResolutionSite() :")+err_no_matbudget;
      logERROR(message);
      return false;
    }
  }


  /**
   * Produces the output of the analysis of the material budget analysis
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::reportTriggerPerformanceSite(bool extended) {
    startTaskClock("Creating trigger summary report");
    if (m_vizard.triggerSummary(m_stripAnalyzer, *m_ctrlStd, m_webSite, extended)) {
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_triggerSummary);
      return false;
    }
  }

  bool Squid::reportNeighbourGraphSite() {

    bool outputOK = false;

    // Report PXD & STD separately
    if (m_ctrlPxd) {
      if (m_vizard.neighbourGraphSummary(*m_pasiveCtrlPxd, m_webSite, "INNER")) outputOK = true;
    }
    if (m_ctrlStd) {
      if (m_vizard.neighbourGraphSummary(*m_pasiveCtrlStd, m_webSite, "OUTER")) outputOK = true;
    }

    // Return
    if (!outputOK) {
      logERROR(err_no_inacsurf);
      return false;
    }
  }

  /*
   * Report additional information
   */
  bool Squid::reportInfoSite() {

    if (m_trackers.size()!=0) {

      startTaskClock("Creating additional info report");

      //m_vizard.additionalInfoSite(m_includeSet, getSettingsFile(), m_pixelAnalyzer, m_stripAnalyzer, m_fwdAnalyzer, m_trackers, m_webSite);
      stopTaskClock();
      return true;
    }
    else {

      std::string message = std::string("Squid::reportInfoSite() :")+err_no_tracker;
      logERROR(err_no_tracker);
      return false;
    }
  }

  void Squid::setBasename(std::string newBasename) {
    m_baseName = newBasename;
  }    

  void Squid::setGeometryFile(std::string geomFile) {

    m_myGeometryFile = geomFile;
    size_t pos = geomFile.find_last_of('.');
    if (pos != string::npos) { geomFile.erase(pos); }
    m_baseName = geomFile;

    // Set layout output htmlDir if not set explicite
    if (m_htmlDir=="") {

      m_htmlDir = m_baseName;
      pos = m_htmlDir.find_last_of('/');
      if (pos != string::npos) { m_htmlDir.erase(pos); }
      m_htmlDir+="/"+insur::default_htmldir;
    }

    // Set layout name according to geometry config file
    std::vector<std::string> info;
    boost::algorithm::split(info, geomFile, boost::algorithm::is_any_of("/"));
    if (info.size()!=0) m_layoutName = info[info.size()-1];
    else                m_layoutName = "";
  }

  void Squid::setHtmlDir(std::string htmlDir) {
    m_htmlDir = htmlDir;
  }


  std::string Squid::getGeometryFile() { 
    if (m_myGeometryFile == "") {
      m_myGeometryFile = m_baseName + suffix_geometry_file;
    }
    return m_myGeometryFile;
  }

  std::string Squid::getSettingsFile() { 
    if (m_mySettingsFile == "") {
      m_mySettingsFile = m_baseName + suffix_types_file;
    }
    return m_mySettingsFile;
  }

  void Squid::simulateTracks(const po::variables_map& varmap, int seed) { // CUIDADO not ported to coderev -- yet?
    startTaskClock("Shooting particles");
/*    TrackShooter ts;
    //std::ofstream ofs((outputfile + "." + any2str(getpid())).c_str());
    //ofs.rdbuf()->pubsetbuf(&(std::vector<char>(32768)[0]), 32768);
    //std::nounitbuf(ofs);
    ts.setOutput(std::cout);
    ts.setTrackerBoundaries(tr->getMaxR(), tr->getMaxBarrelZ(-1), tr->getMaxBarrelZ(+1));
    for (std::vector<Layer*>::const_iterator lit = tr->getLayers().begin(); lit != tr->getLayers().end(); ++lit) {
      ModuleVector* mods = (*lit)->getModuleVector();
      for (ModuleVector::const_iterator mit = mods->begin(); mit != mods->end(); ++mit) {
        ts.addModule(*mit);
      }
    }

    ts.shootTracks(varmap, seed);*/
    stopTaskClock();
  }

  //! TODO: Rework the global constant g!!! That's horrible approach
  void Squid::setCommandLine(int argc, char* argv[]) {
//    if (argc <= 1) return;
//    std::string cmdLine(argv[1]);
//    g=0; for (int i = 2; i < argc; i++) { if (argv[i] == "-"+std::string(1,103)) g=1; cmdLine += std::string(" ") + argv[i]; }
//    m_vizard.setCommandLine(cmdLine);
  }

} // Namespace






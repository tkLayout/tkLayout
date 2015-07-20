/**
 * void littleTestBench();
 * @file Squid.cc
 * @brief This implements the main interface between the tkgeometry library classes and the frontend
 */

#include <SvnRevision.h>
#include "Squid.h"
#include "StopWatch.h"

namespace insur {
  // public
  /**
   * The constructor sets the internal pointers to <i>NULL</i>.
   */
  Squid::Squid() :
      mainConfig_(mainConfigHandler::instance()),
      t2c(mainConfig_),
      weightDistributionTracker(0.1),
      weightDistributionPixel(0.1) {

    //bp_      = nullptr;
    pxd_      = nullptr;
    std_      = nullptr;
    fwdpxd_   = nullptr;
    fwdstd_   = nullptr;

    pxdPasive_= nullptr;
    stdPasive_= nullptr;

    bpMb_     = nullptr;
    pxdMb_    = nullptr;
    stdMb_    = nullptr;
    fwdpxdMb_ = nullptr;
    fwdstdMb_ = nullptr;

    simParms_                = nullptr;
    sitePrepared_            = false;
    myGeometryFile_          = "";
    mySettingsFile_          = "";
    myMaterialFile_          = "";
    myPixelMaterialFile_     = "";
    defaultMaterialFile      = false;
    defaultPixelMaterialFile = false;
  }

  /**
   * The destructor deletes the heap-allocated internal objects if they exist.
   */
  Squid::~Squid() {

    //if (bp_) delete bp_;
    if (pxd_)    delete pxd_;
    if (std_)    delete std_;
    if (fwdpxd_) delete fwdpxd_;
    if (fwdstd_) delete fwdstd_;

    if (pxdPasive_) delete pxdPasive_;
    if (stdPasive_) delete stdPasive_;

    if (bpMb_)  delete bpMb_;
    if (pxdMb_) delete pxdMb_;
    if (stdMb_) delete stdMb_;
    if (fwdpxdMb_) delete fwdpxdMb_;
    if (fwdstdMb_) delete fwdstdMb_;

    if (simParms_) delete simParms_;
    //if (pixelAnalyzer_) delete pixelAnalyzer_;
    StopWatch::destroy();
  }

  /**
   * Build a bare-bones geometry of active modules. The resulting tracker object replaces the previously
   * registered one, if such an object existed. It remains in the squid until it is overwritten by a second call
   * to this function or by another one that creates a new tracker object. If the geometry file contains a
   * description of a pixel detector, that is created as well. It is treated exactly like the tracker with respect
   * to replacement and persistence.
   * @param geomfile The name and - if necessary - path of the geometry configuration file
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::buildTracker() {

    if (pxd_) {
      delete pxd_;
      pxd_ = nullptr;
    }
    if (std_) {
      delete std_;
      std_ = nullptr;
    }
    if (fwdpxd_) {
      delete fwdpxd_;
      fwdpxd_ = nullptr;
    }
    if (fwdstd_) {
      delete fwdstd_;
      fwdstd_ = nullptr;
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
    includeSet_ = mainConfig_.preprocessConfiguration(ifs, ss, getGeometryFile());
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
        if      (trk->myid()=="Inner"          || trk->myid()=="BRL_Inner" || trk->myid()=="BarrelPixel"  ) {
          pxd_    = trk;
          pxd_->setIsPixelType(true);
        }
        else if (trk->myid()=="ForwardPixels"  || trk->myid()=="FWDPXD" || trk->myid()=="ForwarddPixel") {
          fwdpxd_ = trk;
          fwdpxd_->setIsPixelType(true);
          fwdpxd_->setIsForwardType(true);
        }
        else if (trk->myid()=="Outer"          || trk->myid()=="BRL_Outer" || trk->myid()=="BarrelStrip"  ) {
          std_    = trk;
          std_->setIsStripType(true);
        }
        else if (trk->myid()=="ForwardOuter"   || trk->myid()=="FWDSTD" || trk->myid()=="ForwardStrip" ) {
          fwdstd_ = trk;
          fwdstd_->setIsStripType(true);
          fwdstd_->setIsForwardType(true);
        }
        else std_ = trk; // General sollution
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

      // Read simulation parameters
      simParms_ = new SimParms();

      // Add default irradiation files to simParm
      for (auto singleIrradiationFile : insur::default_irradiationfiles) {
        simParms_->addIrradiationMapFile(mainConfig_.getIrradiationDirectory() + "/" + singleIrradiationFile);
      }

      // Read sim parameters file name from a general tree & build simParms
      simParms_->store(getChild(pTree, "SimParms"));
      simParms_->crosscheck();

      // Set sim parameters to analysis
      trackerAnalyzer_.simParms(simParms_);
      stripAnalyzer_.simParms(simParms_);
      pixelAnalyzer_.simParms(simParms_);

      // Look for tag support and build supports
      childRange = getChildRange(pTree, "Support");
      std::for_each(childRange.first, childRange.second, [&](const ptree::value_type& kv) {

        Support* support = new Support();
        support->myid(kv.second.get_value(0));
        support->store(kv.second);
        support->build();
        supports_.push_back(support);
      });
    }
    catch (PathfulException& e) { 
      std::cerr << e.path() << " : " << e.what() << std::endl; 
      stopTaskClock();
      return false;
    }

    stopTaskClock();
    return true;
  }

/**
 * Dress the previously created geometry with module options. The modified tracker object remains
 * in the squid as the current tracker until it is overwritten by a call to another function that creates
 * a new one. If the tracker object has not been created yet, the function fails. If a pixel detector was
 * also created in a previous step, it is dressed here as well. Just like the tracker object, the modified
 * pixel detector remains in the squid until it is replaced by a call to another function creating a new
 * one.
 * @param settingsfile The name and - if necessary - path of the module settings configuration file
 * @return True if there was an existing tracker to dress, false otherwise
 */
/*  bool Squid::dressTracker() {
    if (tr) {
      startTaskClock("Assigning module types to tracker and pixel"); 
      cp.dressTracker(tr, getSettingsFile());
      if (px) cp.dressPixels(px, getSettingsFile());
      stopTaskClock();
      return true;
    }
    else {
      logERROR(err_no_tracker);
      stopTaskClock();
      return false;
    }
  }
*/

  /**
   * Build up the inactive surfaces around the previously created tracker geometry. The resulting collection
   * of inactive surfaces replaces the previously registered one, if such an object existed. It remains in the
   * squid until it is overwritten by a second call to this function or by another one that creates a new
   * collection of inactive surfaces. If a pixel detector was also created in a previous step, it gets a new set
   * of inactive surfaces as well.
   * @param verbose A flag that turns the final status summary of the inactive surface placement algorithm on or off
   * @return True if there were no errors during processing, false otherwise
   */
/*  bool Squid::buildInactiveSurfaces(bool verbose) {
    startTaskClock("Building inactive surfaces");
    if (getGeometryFile()!="") {
      if (std_) {
        if (stdPasive_) delete stdPasive_;
        is = new InactiveSurfaces();
        u.arrange(*std_, *stdPasive_, supports_, verbose);
    }
    if (pxd_) {
        if (pxdPasive_) delete pxdPasive_;
          pxdPasive_ = new InactiveSurfaces();
          u.arrangePixels(*pxd_, *pxdPasive_, verbose);
        }
        stopTaskClock();
        return true;
      }
      else {
        logERROR(err_no_tracker);
        stopTaskClock();
        return false;
      }
    }
    else {
      logERROR(err_no_geomfile);
      stopTaskClock();
      return false;
    }
  } */

/*
 *  Build database of overall material distribution in the tracker ...
 */
bool Squid::buildMaterials(bool verbose) {

  if (std_ || pxd_) {

    startTaskClock("Building database of material distribution in tracker");
    if (std_) {
      if (!stdPasive_) stdPasive_ = new InactiveSurfaces();
      Materialway materialwayStrip;
      materialwayStrip.build(*std_, *stdPasive_, weightDistributionTracker);
    }
    if (pxd_) {
      if (!pxdPasive_) pxdPasive_ = new InactiveSurfaces();
      Materialway materialwayPixel;
      materialwayPixel.build(*pxd_, *pxdPasive_, weightDistributionPixel);
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

  if (std_ || pxd_) {

    startTaskClock("Calculating material budget of the tracker");
    if (std_) {

      if (!stdPasive_) stdPasive_ = new InactiveSurfaces();
      if (stdMb_) delete stdMb_;
      stdMb_  = new MaterialBudget(*std_, *stdPasive_);

      // Reste & calculate new
      if (stripMaterialCalc_.initDone()) stripMaterialCalc_.reset();
      if (matParser_.initMatCalc(stripMaterialCalc_, mainConfig_.getMattabDirectory())) {
        //stdMb_->materialsAll(stripMaterialCalc_);
        if (verbose) stdMb_->print();
      }
    }
    if (pxd_) {

      if (!pxdPasive_) pxdPasive_ = new InactiveSurfaces();
      if (pxdMb_) delete pxdMb_;
      pxdMb_ = new MaterialBudget(*pxd_, *pxdPasive_);

      // Reset & calculate new
      if (pixelMaterialCalc_.initDone()) pixelMaterialCalc_.reset();
	    if (matParser_.initMatCalc(pixelMaterialCalc_, mainConfig_.getMattabDirectory())) {

	      //pxdMb_->materialsAll(pixelMaterialCalc_);
	      if (verbose) pxdMb_->print();
	    }
    }
    stopTaskClock();
    return true;
  } else {
    if (stdMb_) {
      delete stdMb_;
      stdMb_ = nullptr;
    }
    if (pxdMb_) {
      delete pxdMb_;
      pxdMb_ = nullptr;
    }
    std::string message = std::string("Squid::createMaterialBudget() :")+err_init_failed;
    logERROR(message);
    stopTaskClock();
    return false;
  }
}

  /**
   * Build a full system consisting of tracker object, collection of inactive surfaces and material budget from the
   * given configuration files. All three objects replace the previously registered ones, if they existed. They remain
   * in the squid  until they are overwritten by a second call to this function or by another one that creates new
   * instances of them.
   * @param geomfile The name and - if necessary - path of the geometry configuration file
   * @param settingsfile The name and - if necessary - path of the module settings configuration file
   * @param usher_verbose A flag that turns the final status summary of the inactive surface placement algorithm on or off
   * @param mat_verbose A flag that turns the final status summary of the material budget on or off
   * @return True if there were no errors during processing, false otherwise
   */
  // bool Squid::buildFullSystem(bool usher_verbose, bool mat_verbose) {
  //   if (buildInactiveSurfaces(usher_verbose)) return createMaterialBudget(mat_verbose) && irradiateTracker();
  //   return false;
  // }

//  /**
//   * Build the feeder/neighbour graph of the previously created collection of inactive surfaces and
//   * save the results in a plain text file.
//   * @param graphout The name - without path - of the designated output file
//   * @return True if there were no errors during processing, false otherwise
//   */
//  bool Squid::analyzeNeighbours(std::string graphout) {
//    if (is) {
//      startTaskClock("Creating inactive materials hierarchy");
//      v.writeNeighbourGraph(*is, graphout);
//      stopTaskClock();
//      return true;
//    }
//    else {
//      std::cout << "Squid::analyzeNeighbours(): " << err_no_inacsurf << std::endl;
//      return false;
//    }
//  }

  /**
   * Translate an existing full tracker and material budget to a series of XML files that can be interpreted by CMSSW.
   * @param xmlout The name - without path - of the designated output subdirectory
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::translateFullSystemToXML(std::string xmlout) {
    if (stdMb_) {
      t2c.translate(stripMaterialCalc_.getMaterialTable(), *stdMb_, xmlout.empty() ? baseName_ : xmlout, false); // false is setting a mysterious flag called wt which changes the way the XML is output. apparently setting it to true is of no use anymore.
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

    vizard_.~Vizard();
    new ((void*) &vizard_) Vizard();
  }

/**
 * Prepare the website object (if not done yet) from the configuration file
 * it needs the tracker object to be already there
 * @return a boolean with the operation success
 */
bool Squid::prepareWebSite() {

  if (sitePrepared_) return true;

  string trackerName;
  if (htmlDir_ != "") trackerName = htmlDir_;
  else {

    if (std_) trackerName = baseName_;
    else trackerName = default_trackername;
  }

  string layoutDirectory;
  //styleDirectory=mainConfiguration.getStyleDirectory();
  layoutDirectory=mainConfig_.getLayoutDirectory();
  layoutDirectory+="/"+trackerName;
  if (layoutDirectory!="") webSite_.setTargetDirectory(layoutDirectory);
  else return false;

  webSite_.setTitle(trackerName);
  webSite_.setComment("Layouts");
  webSite_.setCommentLink("../");
  webSite_.addAuthor("Giovanni Bianchi");
  webSite_.addAuthor("Nicoletta De Maio");
  webSite_.addAuthor("Stefano Martina");
  webSite_.addAuthor("Stefano Mersi");
  webSite_.setRevision(SvnRevision::revisionNumber);
  return true;
}

/**
 * Actually creates the website where it was supposed to be
 * @return a boolean with the operation success
 */
bool Squid::makeWebSite(bool addLogPage /* = true */) {

  startTaskClock("Creating website");
  if (!prepareWebSite()) {

    logERROR("Problem in preparing website");
    stopTaskClock();
    return false;
  }

  // Add log webPage
  if (addLogPage) {
    vizard_.makeLogPage(webSite_);
  }

  bool result = webSite_.makeSite(false);
  stopTaskClock();
  return result;
}

/**
 * Analyze the previously created geometry and without no output  through rootweb.
 * @return True if there were no errors during processing, false otherwise
 */
bool Squid::pureAnalyzeGeometry(int tracks) {

  // Strip detector
  if (pxd_ || std_) {

    startTaskClock("Analyzing geometry");
    if (pxd_) pixelAnalyzer_.analyzeGeometry(*pxd_, tracks);
    if (std_) stripAnalyzer_.analyzeGeometry(*std_, tracks);
    stopTaskClock();
    return true;
  }
  else {
    std::string message = std::string("Squid::pureAnalyzeGeometry(): ")+err_no_tracker;
    logERROR(message);
    return false;
  }
}

  bool Squid::analyzeTriggerEfficiency(int tracks, bool detailed) {
    // Call this before analyzetrigger if you want to have the map of suggested spacings
    if (detailed) {
      startTaskClock("Creating distance tuning plots");
      stripAnalyzer_.createTriggerDistanceTuningPlots(*std_, mainConfig_.getTriggerMomenta());
      stopTaskClock();
    }
    startTaskClock("Creating trigger efficiency plots");
    stripAnalyzer_.analyzeTriggerEfficiency(*std_,
                               mainConfig_.getTriggerMomenta(),
                               mainConfig_.getThresholdProbabilities(),
                               tracks);
    stopTaskClock();
    return true;
  }

/**
 * Analyze the previously created material budget with no output.
 * @param tracks The number of tracks that should be fanned out across the analysed region
 * @return True if there were no errors during processing, false otherwise
 */
bool Squid::pureAnalyzeMaterialBudget(int nTracks, bool trkResolution) {

  if (stdMb_ || pxdMb_) {
    // TODO: insert the creation of sample tracks here, to compute intersections only once

    startTaskClock("Analyzing material budget" );
    if (stdMb_) stripAnalyzer_.analyzeMaterialBudget(stdMb_, nTracks);
    if (pxdMb_) pixelAnalyzer_.analyzeMaterialBudget(pxdMb_, nTracks);
    stopTaskClock();

    //startTaskClock("Computing the weight summary");
    //stripAnalyzer_.computeWeightSummary(*stdMb_);
    //stopTaskClock();

    if (trkResolution) {

      startTaskClock("Analyzing tracking resolutions");
      trackerAnalyzer_.analyzeTaggedTracking(stdMb_, pxdMb_,
                                             mainConfig_.getMomenta(),
                                             mainConfig_.getTriggerMomenta(),
                                             mainConfig_.getThresholdProbabilities(),
                                             nTracks);
      stopTaskClock();
    }
    return true;
  } else {

    std::string message = std::string("Squid::pureAnalyzeMaterialBudget() :")+err_no_matbudget;
    logERROR(message);
    return false;
  }
}

/**
 * Produces the output of the analysis of the geomerty
 * @return True if there were no errors during processing, false otherwise
 */
bool Squid::reportGeometrySite() {

  if (pxd_ || std_) {

    startTaskClock("Creating geometry report");
    if (std_) vizard_.geometrySummary(stripAnalyzer_, *std_, *simParms_, stdPasive_, webSite_,"STD");
    if (pxd_) vizard_.geometrySummary(pixelAnalyzer_, *pxd_, *simParms_, pxdPasive_, webSite_,"PXD");
    stopTaskClock();
    return true;
  } else {

    std::string message = std::string("Squid::reportGeometrySite(): ")+err_no_tracker;
    logERROR(message);
    return false;
  }
}

  bool Squid::reportBandwidthSite() {
    if (std_) {
      startTaskClock("Computing bandwidth and rates");
      stripAnalyzer_.computeBandwidth(*std_);
      stripAnalyzer_.computeTriggerFrequency(*std_);
      stopTaskClock();
      startTaskClock("Creating bandwidth and rates report");
      vizard_.bandwidthSummary(stripAnalyzer_, *std_, *simParms_, webSite_);
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }
  }

  bool Squid::reportTriggerProcessorsSite() {
    if (std_) {
      startTaskClock("Computing multiple trigger tower connections");
      stripAnalyzer_.computeTriggerProcessorsBandwidth(*std_);
      vizard_.triggerProcessorsSummary(stripAnalyzer_, *std_, webSite_);
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }

  }

  bool Squid::reportPowerSite() {
    if (std_) {
      startTaskClock("Computing dissipated power");
      stripAnalyzer_.analyzePower(*std_);
      stopTaskClock();
      startTaskClock("Creating power report");
      vizard_.irradiatedPowerSummary(stripAnalyzer_, *std_, webSite_);
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

    if (stdMb_ || pxdMb_) {

      startTaskClock("Creating material budget report");
      if (pxdMb_) vizard_.materialSummary(pixelAnalyzer_, *pxdMb_, debugServices, webSite_, "PXD");
      if (stdMb_) vizard_.materialSummary(stripAnalyzer_, *stdMb_, debugServices, webSite_, "STD");
      //if (pxdMb_) vizard_.weigthSummart(pixelAnalyzer_, weightDistributionPixel  , webSite_, "PXD");
      //if (stdMb_) vizard_.weigthSummart(stripAnalyzer_, weightDistributionTracker, webSite_, "STD");
      
      
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

  if (std_ || pxd_) {

    startTaskClock("Creating resolution report");
    vizard_.taggedErrorSummary(trackerAnalyzer_, webSite_);
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
    if (vizard_.triggerSummary(stripAnalyzer_, *std_, webSite_, extended)) {
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_triggerSummary);
      return false;
    }
  }

//  bool Squid::reportNeighbourGraphSite() {
//    if (v.neighbourGraphSummary(*is, webSite_)) return true;
//    else {
//      logERROR(err_no_inacsurf);
//      return false;
//    }
//  }

/*
 * Report additional information
 */
bool Squid::reportInfoSite() {

  if (std_ || pxd_) {

    startTaskClock("Creating additional info report");
    vizard_.additionalInfoSite(includeSet_, getSettingsFile(),
                           stripAnalyzer_, pixelAnalyzer_, *std_, *simParms_, webSite_);
    stopTaskClock();
    return true;
  }
  else {
    logERROR(err_no_tracker);
    return false;
  }
}

  void Squid::setBasename(std::string newBasename) {
    baseName_ = newBasename;
  }    

  void Squid::setGeometryFile(std::string geomFile) {
    myGeometryFile_ = geomFile;
    size_t pos = geomFile.find_last_of('.');
    if (pos != string::npos) { geomFile.erase(pos); }
    baseName_ = geomFile;
  }

  void Squid::setHtmlDir(std::string htmlDir) {
    htmlDir_ = htmlDir;
  }


  std::string Squid::getGeometryFile() { 
    if (myGeometryFile_ == "") {
      myGeometryFile_ = baseName_ + suffix_geometry_file;
    }
    return myGeometryFile_;
  }

  std::string Squid::getSettingsFile() { 
    if (mySettingsFile_ == "") {
      mySettingsFile_ = baseName_ + suffix_types_file;
    }
    return mySettingsFile_;
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

  void Squid::setCommandLine(int argc, char* argv[]) {
    if (argc <= 1) return;
    std::string cmdLine(argv[1]);
    g=0; for (int i = 2; i < argc; i++) { if (argv[i] == "-"+std::string(1,103)) g=1; cmdLine += std::string(" ") + argv[i]; }
    vizard_.setCommandLine(cmdLine);
  }
}






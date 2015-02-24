/**
 * void littleTestBench();
 * @file Squid.cc
 * @brief This implements the main interface between the tkgeometry library classes and the frontend
 */

#include "SvnRevision.h"
#include "Squid.h"
#include "StopWatch.h"

namespace insur {
  // public
  /**
   * The constructor sets the internal pointers to <i>NULL</i>.
   */
  Squid::Squid() :
      mainConfiguration(mainConfigHandler::instance()),
      t2c(mainConfiguration),
      weightDistributionTracker(0.1),
      weightDistributionPixel(0.1) {
    tr = NULL;
    is = NULL;
    mb = NULL;
    px = NULL;
    pi = NULL;
    pm = NULL;
    //pixelAnalyzer = NULL;
    sitePrepared = false;
    myGeometryFile_ = "";
    mySettingsFile_ = "";
    myMaterialFile_ = "";
    myPixelMaterialFile_ = "";
    defaultMaterialFile = false;
    defaultPixelMaterialFile = false;
  }

  /**
   * The destructor deletes the heap-allocated internal objects if they exist.
   */
  Squid::~Squid() {
    if (mb) delete mb;
    if (is) delete is;
    if (tr) delete tr;
    if (pm) delete pm;
    if (pi) delete pi;
    if (px) delete px;
    //if (pixelAnalyzer) delete pixelAnalyzer;    
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
    if (tr) delete tr;
    if (px) delete px;
    tr = NULL;
    px = NULL;

    std::ifstream ifs(getGeometryFile());
    if (ifs.fail()) {
      std::cerr << "ERROR: cannot open geometry file " << getGeometryFile() << std::endl;
      return false;
    }
    startTaskClock("Building tracker and pixel");
    std::stringstream ss;
    includeSet_ = mainConfiguration.preprocessConfiguration(ifs, ss, getGeometryFile());
    t2c.addConfigFile(tk2CMSSW::ConfigFile{getGeometryFile(), ss.str()});
    using namespace boost::property_tree;
    ptree pt;
    info_parser::read_info(ss, pt);

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
      auto childRange = getChildRange(pt, "Tracker");
      std::for_each(childRange.first, childRange.second, [&](const ptree::value_type& kv) {
        Tracker* t = new Tracker();
        t->setup();
        t->myid(kv.second.data());
        t->store(kv.second);
        t->build();
        //CoordExportVisitor v(t->myid());
        //ModuleDataVisitor v1(t->myid());
        //t->accept(v);
        //t->accept(v1);
        if (t->myid() == "Pixels") px = t;
        else tr = t;
      });

      std::set<string> unmatchedProperties = PropertyObject::reportUnmatchedProperties();
      if (!unmatchedProperties.empty()) {
        std::ostringstream ss;
        ss << "The following unknown properties were ignored:" << std::endl;
        for (const string& s : unmatchedProperties) {
          ss << "  " << s << std::endl;
        }
        logERROR(ss);
      }

      simParms_ = new SimParms();

      //iter between the default irradiation files vector and add each to simParm
      for (auto singleIrradiationFile : insur::default_irradiationfiles) {
        simParms_->addIrradiationMapFile(mainConfiguration.getIrradiationDirectory() + "/" + singleIrradiationFile);
      }
      //simParms_->irradiationMapFile(mainConfiguration.getIrradiationDirectory() + "/" + insur::default_irradiationfile);
      simParms_->store(getChild(pt, "SimParms"));
      simParms_->build();
      a.simParms(simParms_);
      pixelAnalyzer.simParms(simParms_);

      childRange = getChildRange(pt, "Support");
      std::for_each(childRange.first, childRange.second, [&](const ptree::value_type& kv) {
        Support* s = new Support();
        s->myid(kv.second.get_value(0));
        s->store(kv.second);
        s->build();
        supports_.push_back(s);
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

 /*
  bool Squid::buildNewTracker() {
    boost::ptree pt;
    info_parser::read_info(getGeometryFile(), pt);
  }
*/
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
  bool Squid::buildInactiveSurfaces(bool verbose) {
    startTaskClock("Building inactive surfaces");
    if (getGeometryFile()!="") {
      if (tr) {
        if (is) delete is;
        is = new InactiveSurfaces();
        u.arrange(*tr, *is, supports_, verbose);
        if (px) {
          if (pi) delete pi;
          pi = new InactiveSurfaces();
          u.arrangePixels(*px, *pi, verbose);
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
  }

  bool Squid::buildMaterials(bool verbose) {
    startTaskClock("Building materials");

    if (tr) {
        std::string trackm = getMaterialFile();
        if (trackm=="") { stopTaskClock(); return false; }
        if (!is) is = new InactiveSurfaces();
        //if (mb) delete mb;
        //mb  = new MaterialBudget(*tr, *is);
        //if (tkMaterialCalc.initDone()) tkMaterialCalc.reset(); // TODO: obsolete these
        //if (pxMaterialCalc.initDone()) pxMaterialCalc.reset(); // TODO: obsolete these

        //if (mp.initMatCalc(trackm, tkMaterialCalc, mainConfiguration.getMattabDirectory())) {
        materialwayTracker.build(*tr, *is, weightDistributionTracker);

          // mb->materialsAll(tkMaterialCalc);
          // if (verbose) mb->print();

          if (px) {
            std::string pixm = getPixelMaterialFile();
            if (pixm!="") {
              //if (mp.initMatCalc(pixm, pxMaterialCalc, mainConfiguration.getMattabDirectory())) {
                if (!pi) pi = new InactiveSurfaces();
                //if (pm) delete pm;
                //pm = new MaterialBudget(*px, *pi);
                materialwayPixel.build(*px, *pi, weightDistributionPixel);

                //pm->materialsAll(pxMaterialCalc);
                //if (verbose) pm->print();
                //}
            }
          }

          /*
        } else {
          //if (mb) delete mb;
          //mb = NULL;
          //if (pm) delete pm;
          //pm = NULL;
          logERROR(err_init_failed);
          return false;
        }
        */
      } else {
        logERROR(err_no_tracker);
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
    if (tr) {
      std::string trackm = getMaterialFile();
      if (trackm=="") return false;
      if (!is) is = new InactiveSurfaces();
      if (mb) delete mb;
      mb  = new MaterialBudget(*tr, *is);
      if (tkMaterialCalc.initDone()) tkMaterialCalc.reset();
      if (pxMaterialCalc.initDone()) pxMaterialCalc.reset();
      if (mp.initMatCalc(trackm, tkMaterialCalc, mainConfiguration.getMattabDirectory())) {
        mb->materialsAll(tkMaterialCalc);
        if (verbose) mb->print();

        if (px) {
          std::string pixm = getPixelMaterialFile();
          if (pixm!="") {
            if (mp.initMatCalc(pixm, pxMaterialCalc, mainConfiguration.getMattabDirectory())) {
              if (!pi) pi = new InactiveSurfaces();
              if (pm) delete pm;
              pm = new MaterialBudget(*px, *pi);
              pm->materialsAll(pxMaterialCalc);
              if (verbose) pm->print();
            }
          }
        }
        return true;
      } else {
        if (mb) delete mb;
        mb = NULL;
        if (pm) delete pm;
        pm = NULL;
        logERROR(err_init_failed);
        return false;
      }
    } else {
      logERROR(err_no_tracker);
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

  /**
   * Build the feeder/neighbour graph of the previously created collection of inactive surfaces and
   * save the results in a plain text file.
   * @param graphout The name - without path - of the designated output file
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::analyzeNeighbours(std::string graphout) {
    if (is) {
      startTaskClock("Creating inactive materials hierarchy");
      v.writeNeighbourGraph(*is, graphout);
      stopTaskClock();
      return true;
    }
    else {
      std::cout << "Squid::analyzeNeighbours(): " << err_no_inacsurf << std::endl;
      return false;
    }
  }

  /**
   * Translate an existing full tracker and material budget to a series of XML files that can be interpreted by CMSSW.
   * @param xmlout The name - without path - of the designated output subdirectory
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::translateFullSystemToXML(std::string xmlout) {
    if (mb) {
      t2c.translate(tkMaterialCalc.getMaterialTable(), *mb, xmlout.empty() ? baseName_ : xmlout, false); // false is setting a mysterious flag called wt which changes the way the XML is output. apparently setting it to true is of no use anymore.
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
    v.~Vizard();
    new ((void*) &v) Vizard();
  }

  /**
   * Prepare the website object (if not done yet) from the configuration file
   * it needs the tracker object to be already there
   * @return a boolean with the operation success
   */
  bool Squid::prepareWebsite() {
    if (sitePrepared) return true;
    string trackerName;
    if (htmlDir_ != "") trackerName = htmlDir_;
    else {
      if (tr) trackerName = baseName_;
      else trackerName = default_trackername;
    }
    string layoutDirectory;
    //styleDirectory=mainConfiguration.getStyleDirectory();
    layoutDirectory=mainConfiguration.getLayoutDirectory();
    layoutDirectory+="/"+trackerName;
    if (layoutDirectory!="") site.setTargetDirectory(layoutDirectory);
    else return false;
    site.setTitle(trackerName);
    site.setComment("layouts");
    site.setCommentLink("../");
    site.addAuthor("Giovanni Bianchi");
    site.addAuthor("Nicoletta De Maio");
    site.addAuthor("Stefano Martina");
    site.addAuthor("Stefano Mersi");
    site.setRevision(SvnRevision::revisionNumber);
    return true;
  }


  /**
   * Actually creates the website where it was supposed to be
   * @return a boolean with the operation success
   */
  bool Squid::makeSite(bool addLogPage /* = true */) {
    startTaskClock("Creating website");
    if (!prepareWebsite()) {
      logERROR("Problem in preparing website");
      stopTaskClock();
      return false;
    }

    if (addLogPage) {
      v.makeLogPage(site);
    }

    bool result = site.makeSite(false);
    stopTaskClock();
    return result;
  }

  /**
   * Analyze the previously created geometry and without no output  through rootweb.
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::pureAnalyzeGeometry(int tracks) {
    if (tr) {
      startTaskClock("Analyzing geometry");
      a.analyzeGeometry(*tr, tracks);
      if (px) pixelAnalyzer.analyzeGeometry(*px, tracks);
      stopTaskClock();
      return true; // TODO: this return value is not really meaningful
    } else {
      std::cout << "Squid::pureAnalyzeGeometry(): " << err_no_tracker << std::endl;
      return false;
    }
  }

  bool Squid::analyzeTriggerEfficiency(int tracks, bool detailed) {
    // Call this before analyzetrigger if you want to have the map of suggested spacings
    if (detailed) {
      startTaskClock("Creating distance tuning plots");
      a.createTriggerDistanceTuningPlots(*tr, mainConfiguration.getTriggerMomenta());
      stopTaskClock();
    }
    startTaskClock("Creating trigger efficiency plots");
    a.analyzeTriggerEfficiency(*tr,
                               mainConfiguration.getTriggerMomenta(),
                               mainConfiguration.getThresholdProbabilities(),
                               tracks);
    stopTaskClock();
    return true;
  }

  /**
   * Analyze the previously created material budget with no output.
   * @param tracks The number of tracks that should be fanned out across the analysed region
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::pureAnalyzeMaterialBudget(int tracks, bool triggerResolution) {
    if (mb) {
//      startTaskClock(!trackingResolution ? "Analyzing material budget" : "Analyzing material budget and estimating resolution");
      // TODO: insert the creation of sample tracks here, to compute intersections only once
      startTaskClock("Analyzing material budget" );
      a.analyzeMaterialBudget(*mb, mainConfiguration.getMomenta(), tracks, pm);
      stopTaskClock();
      if (pm) {
        startTaskClock("Analyzing pixel material budget");
        pixelAnalyzer.analyzeMaterialBudget(*pm, mainConfiguration.getMomenta(), tracks, NULL);
        stopTaskClock();
      }
      startTaskClock("Computing the weight summary");
      a.computeWeightSummary(*mb);
      stopTaskClock();
      if (triggerResolution) {
        startTaskClock("Estimating tracking resolutions");
        a.analyzeTaggedTracking(*mb,
                                mainConfiguration.getMomenta(),
                                mainConfiguration.getTriggerMomenta(),
                                mainConfiguration.getThresholdProbabilities(),
                                tracks, pm);
        stopTaskClock();
      }
      return true;
    } else {
      logERROR(err_no_matbudget);
      return false;
    }
  }

  /**
   * Produces the output of the analysis of the geomerty analysis
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::reportGeometrySite() {
    if (tr) {
      startTaskClock("Creating geometry report");
      v.geometrySummary(a, *tr, *simParms_, is, site);
      if (px) v.geometrySummary(pixelAnalyzer, *px, *simParms_, pi, site, "pixel");
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }
  }

  bool Squid::reportBandwidthSite() {
    if (tr) {
      startTaskClock("Computing bandwidth and rates");
      a.computeBandwidth(*tr);
      a.computeTriggerFrequency(*tr);
      stopTaskClock();
      startTaskClock("Creating bandwidth and rates report");
      v.bandwidthSummary(a, *tr, *simParms_, site);
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }
  }

  bool Squid::reportTriggerProcessorsSite() {
    if (tr) {
      startTaskClock("Computing multiple trigger tower connections");
      a.computeTriggerProcessorsBandwidth(*tr);
      v.triggerProcessorsSummary(a, *tr, site);
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }

  }

  bool Squid::reportPowerSite() {
    if (tr) {
      startTaskClock("Computing dissipated power");
      a.analyzePower(*tr);
      stopTaskClock();
      startTaskClock("Creating power report");
      v.irradiatedPowerSummary(a, *tr, site);
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
  bool Squid::reportMaterialBudgetSite() {
    if (mb) {
      startTaskClock("Creating material budget report");
      v.histogramSummary(a, site, "outer");
      if (pm) v.histogramSummary(pixelAnalyzer, site, "pixel");
      v.weigthSummart(a, weightDistributionTracker, site, "outer");
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
    if (mb) {
      startTaskClock("Creating resolution report");
      v.errorSummary(a, site, "", false);
#ifdef NO_TAGGED_TRACKING
      v.errorSummary(a, site, "trigger", true);
#else
      v.taggedErrorSummary(a, site);
#endif
      stopTaskClock();
      return true;
    }
    else {
      logERROR(err_no_matbudget);
      return false;
    }
  }


  /**
   * Produces the output of the analysis of the material budget analysis
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::reportTriggerPerformanceSite(bool extended) {
    startTaskClock("Creating trigger summary report");
    if (v.triggerSummary(a, *tr, site, extended)) {
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_triggerSummary);
      return false;
    }
  }

  bool Squid::reportNeighbourGraphSite() {
    if (v.neighbourGraphSummary(*is, site)) return true;
    else {
      logERROR(err_no_inacsurf);
      return false;
    }
  }

  bool Squid::additionalInfoSite() {
    if (!tr) {
      logERROR(err_no_tracker);
      return false;
    } else {
      getMaterialFile();
      getPixelMaterialFile();
      startTaskClock("Saving additional information");
      v.additionalInfoSite(includeSet_, getSettingsFile(),
                           getMaterialFile(), getPixelMaterialFile(),
                           defaultMaterialFile, defaultPixelMaterialFile,
                           a, pixelAnalyzer, *tr, *simParms_, site);
      stopTaskClock();
      return true;
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

  std::string Squid::getMaterialFile() {
    if (myMaterialFile_=="") {
      myMaterialFile_ = baseName_ + suffix_tracker_material_file; 
      if (fileExists(myMaterialFile_)) {
        logWARNING(warn_custom_matfile);
        defaultMaterialFile = false;
      } else {
        myMaterialFile_ = mainConfiguration.getDefaultMaterialsDirectory() + "/" + default_tracker_materials_file;
        defaultMaterialFile = true;
        if (!fileExists(myMaterialFile_)) {
          logERROR(err_no_matfile);
          myMaterialFile_ = ""; // TODO: put an "undefined" here to mark your passage
        }
      }
    }
    return myMaterialFile_;
  }

  std::string Squid::getPixelMaterialFile() {
    if (myPixelMaterialFile_=="") {
      myPixelMaterialFile_ = baseName_ + suffix_pixel_material_file;         
      if (fileExists(myPixelMaterialFile_)) {
        logWARNING(warn_custom_matfile_pixel);
        defaultPixelMaterialFile = false;
      } else {
        myPixelMaterialFile_ = mainConfiguration.getDefaultMaterialsDirectory() + "/" + default_pixel_materials_file;
        defaultPixelMaterialFile = true;
        if (!fileExists(myPixelMaterialFile_)) {
          logERROR(err_no_matfile_pixel);
          myPixelMaterialFile_ = ""; // TODO: put an "undefined" here to mark your passage
        }
      }
    }
    return myPixelMaterialFile_;
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
    v.setCommandLine(cmdLine);
  }
}






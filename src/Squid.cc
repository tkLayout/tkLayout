/**
 * @file Squid.cc
 * @brief This implements the main interface between the tkgeometry library classes and the frontend
 */

#include <Squid.h>
namespace insur {
  // public
  /**
   * The constructor sets the internal pointers to <i>NULL</i>.
   */
  Squid::Squid() : t2c(mainConfiguration) {
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
    startTaskClock("Building tracker and pixel");
    tr = cp.parseFile(getGeometryFile(), getSettingsFile());
    if (tr) {
      if (px) delete px;
      px = cp.parsePixelsFromFile(getGeometryFile());
      if (px) px->setZError(tr->getZError());
      stopTaskClock();
      return true;
    }
    stopTaskClock();
    return false;
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
  bool Squid::dressTracker() {
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

  /**
   * Irradiate a previously created tracker.
   * @return True if there was an existing tracker to irradiate, false otherwise
   */
  bool Squid::irradiateTracker() {
    if (tr) {
      startTaskClock("Evaluating modules irradiation");
      cp.irradiateTracker(tr, mainConfiguration.getIrradiationDirectory() + "/" + insur::default_irradiationfile);
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }
  }

    
  /**
   * Build a geometry of active modules using both geometry constraints and module settings. If there
   * was an existing tracker object, it is destroyed and replaced by a new one as described in the geometry
   * and settings configuration files. If there is a pixel detector, it gets the same treatment and may not
   * exist anymore afterwards, depending on whether the new geometry file describes one or not.
   * @param geomfile The name and - if necessary - path of the geometry configuration file
   * @param settingsfile The name and - if necessary - path of the module settings configuration file
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::buildTrackerSystem() {
    if ( buildTracker() ) return dressTracker();
    else return false;
  }
  
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
	u.arrange(*tr, *is, getGeometryFile(), verbose);
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
  bool Squid::translateFullSystemToXML(std::string xmlout, bool wt) {
    if (mb) {
      t2c.translate(tkMaterialCalc.getMaterialTable(), *mb, xmlout, wt);
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
        if (tr) trackerName = tr->getName();
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
    site.addAuthor("Stefano Mersi");
#ifdef REVISIONNUMBER
    site.setRevision(REVISIONNUMBER);
#endif
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
      startTaskClock("Analyzing material budget and estimating resolution");
      a.analyzeMaterialBudget(*mb, mainConfiguration.getMomenta(), tracks, pm, true);
      stopTaskClock();
      if (pm) {
        startTaskClock("Analyzing pixel material budget");
        pixelAnalyzer.analyzeMaterialBudget(*pm, mainConfiguration.getMomenta(), tracks, NULL, false);
        stopTaskClock();
      }
      a.computeWeightSummary(*mb);
      if (triggerResolution) {
	startTaskClock("Estimating tracking resolution of track-trigger");
	a.analyzeTrigger(*mb,
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
      v.geometrySummary(a, *tr, site);
      if (px) v.geometrySummary(pixelAnalyzer, *px, site, "pixel");
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
      v.bandwidthSummary(a, *tr, site);
      stopTaskClock();
      return true;
    } else {
      logERROR(err_no_tracker);
      return false;
    }
  }

bool Squid::reportTriggerProcessorsSite() {
    if (tr) {
      startTaskClock("Computing dissipated power");
        a.computeTriggerProcessorsBandwidth(*tr);
        v.triggerProcessorsSummary(a, *tr, site);
        return true;
    } else {
    std::cerr << "Squid::reportTriggerProcessorsSite(): " << err_no_tracker << std::endl;
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
      v.weigthSummart(a, site, "outer");
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
      v.errorSummary(a, site, "trigger", true);
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
    if (v.triggerSummary(a, site, extended)) {
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
      v.additionalInfoSite(getGeometryFile(), getSettingsFile(),
			   getMaterialFile(), getPixelMaterialFile(),
			   defaultMaterialFile, defaultPixelMaterialFile,
			   a, *tr, site);
      stopTaskClock();
      return true;
    }
  }

  void Squid::setBasename(std::string newBasename) {
    baseName_ = newBasename;
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

}

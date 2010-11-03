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
    Squid::Squid() {
        tr = NULL;
        is = NULL;
        mb = NULL;
        px = NULL;
        pi = NULL;
        pm = NULL;
	pixelAnalyzer = NULL;
#ifdef USING_ROOTWEB
        sitePrepared = false;
#endif
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
	if (pixelAnalyzer) delete pixelAnalyzer;
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
    bool Squid::buildTracker(std::string geomfile) {
        if (tr) delete tr;
        tr = cp.parseFile(geomfile);
        if (tr) {
            if (px) delete px;
            px = cp.parsePixelsFromFile(geomfile);
            if (px) px->setZError(tr->getZError());
            g = geomfile;
            return true;
        }
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
    bool Squid::dressTracker(std::string settingsfile) {
        if (tr) {
            cp.dressTracker(tr, settingsfile);
            if (px) cp.dressPixels(px, settingsfile);
            return true;
        }
        else {
            std::cout << "Squid::dressTracker(): " << err_no_tracker << std::endl;
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
    bool Squid::buildTrackerSystem(std::string geomfile, std::string settingsfile) {
        if (buildTracker(geomfile)) return dressTracker(settingsfile);
        else {
            g = "";
            return false;
        }
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
        if (!g.empty()) {
            if (tr) {
                if (is) delete is;
                is = new InactiveSurfaces();
                u.arrange(*tr, *is, g, verbose);
                if (px) {
                    if (pi) delete pi;
                    pi = new InactiveSurfaces();
                    u.arrangePixels(*px, *pi, verbose);
                }
                return true;
            }
            else {
                std::cout << "Squid::buildInactiveSurfaces(): " << err_no_tracker << std::endl;
                return false;
            }
        }
        else {
            std::cout << "Squid::buildInactiveSurfaces(): " << err_no_geomfile << std::endl;
            return false;
        }
    }
    
    /**
     * Build up a bare-bones geometry of active modules, then arrange the inactive surfaces around it. The
     * resulting tracker object and collection of inactive surfaces replace the previously registered ones, if
     * such objects existed. They remain in the squid until they are overwritten by a second call to this
     * function or by another one that creates a new tracker object and collection of inactive surfaces. If a
     * pixel detector is specified in the geometry and settings files, it recieves the same treatment as the
     * tracker.
     * @param geomfile The name and - if necessary - path of the geometry configuration file
     * @param verbose A flag that turns the final status summary of the inactive surface placement algorithm on or off
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::buildInactiveSurfaces(std::string geomfile, bool verbose) {
        if (is) delete is;
        if (buildTracker(geomfile)) {
            is = new InactiveSurfaces();
            u.arrange(*tr, *is, geomfile, verbose);
            if (px) {
                if (pi) delete pi;
                pi = new InactiveSurfaces();
                u.arrangePixels(*px, *pi, verbose);
            }
            return true;
        }
        is = NULL;
        return false;
    }
    
    /**
     * Build a geometry of active modules using both geometry constraints and module settings, then
     * arrange the inactive surfaces around it. The resulting tracker object and collection of inactive
     * surfaces replace the previously registered ones, if such objects existed. They remain in the squid
     * until they are overwritten by a second call to this function or by another one that creates a new
     * tracker object and collection of inactive surfaces. If a pixel detector was also created in a previous
     * step, its inactive surfaces are treated exactly the same as those of the tracker.
     * @param geomfile The name and - if necessary - path of the geometry configuration file
     * @param settingsfile The name and - if necessary - path of the module settings configuration file
     * @param verbose A flag that turns the final status summary of the inactive surface placement algorithm on or off
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::buildInactiveSurfaces(std::string geomfile, std::string settingsfile, bool verbose) {
        if (is) delete is;
        if (buildTrackerSystem(geomfile, settingsfile)) {
            is = new InactiveSurfaces();
            u.arrange(*tr, *is, geomfile, verbose);
            if (px) {
                if (pi) delete pi;
                pi = new InactiveSurfaces();
                u.arrangePixels(*px, *pi, verbose);
            }
            return true;
        }
        is = NULL;
        return false;
    }

   /**
     * Calculate a material budget for the previously created tracker object and collection of inactive
     * surfaces. The resulting material budget replaces the previously registered one, if such an object
     * existed. It remains in the squid until it is overwritten by a second call to this function or by
     * another one that creates a new material budget. This function succeeds if either the tracker or
     * both the tracker and the inactive surfaces exist. If a pixel detector exists, it gets its own material
     * budget at the same time. This object is treated exactly the same as the one for the tracker.
     * @param matfile The name and - if necessary - path of the materials configuration file
     * @param pixmatfile The name and - if necessary - path of the materials configuration file fro the pixels
     * @param verbose A flag that turns the final status summary of the material budget on or off
     * @return True if there were no errors during processing, false otherwise
     */
  bool Squid::createMaterialBudget(std::string matfile, bool verbose) {   
     std::string dummyString;
     return createMaterialBudget(matfile, dummyString, verbose);
  }
 
    /**
     * Calculate a material budget for the previously created tracker object and collection of inactive
     * surfaces. The resulting material budget replaces the previously registered one, if such an object
     * existed. It remains in the squid until it is overwritten by a second call to this function or by
     * another one that creates a new material budget. This function succeeds if either the tracker or
     * both the tracker and the inactive surfaces exist. If a pixel detector exists, it gets its own material
     * budget at the same time. This object is treated exactly the same as the one for the tracker.
     * @param matfile The name and - if necessary - path of the materials configuration file
     * @param pixmatfile The name and - if necessary - path of the materials configuration file fro the pixels
     * @param verbose A flag that turns the final status summary of the material budget on or off
     * @return True if there were no errors during processing, false otherwise
     */
  bool Squid::createMaterialBudget(std::string matfile, std::string& pixmatfile, bool verbose) {
        if (fileExists(matfile)) {
            if (tr) {
                if (!is) is = new InactiveSurfaces();
                if (mb) delete mb;
                mb  = new MaterialBudget(*tr, *is);
                if (tkMaterialCalc.initDone()) tkMaterialCalc.reset();
                if (pxMaterialCalc.initDone()) pxMaterialCalc.reset();
                if (mp.initMatCalc(matfile, tkMaterialCalc, mainConfiguration.getMattabDirectory())) {
                    mb->materialsAll(tkMaterialCalc);
                    if (verbose) mb->print();
                    if (px) {
		      if (fileExists(pixmatfile)) {
                        if (mp.initMatCalc(pixmatfile, pxMaterialCalc, mainConfiguration.getMattabDirectory())) {
			  if (!pi) pi = new InactiveSurfaces();
			  if (pm) delete pm;
			  pm = new MaterialBudget(*px, *pi);
			  pm->materialsAll(pxMaterialCalc);
			  if (verbose) pm->print();
                        }
		      } else {
			std::cout << "Squid::createMaterialBudget(): " << err_no_pixmatfile << std::endl;
			return false;
		      }
                    } else pixmatfile="";
                    return true;
                }
                else {
                    if (mb) delete mb;
                    mb = NULL;
                    if (pm) delete pm;
                    pm = NULL;
                    std::cout << "Squid::createMaterialBudget(): " << err_init_failed << std::endl;
                    return false;
                }
            }
            else {
                std::cout << "Squid::createMaterialBudget(): " << err_no_tracker << std::endl;
                return false;
            }
        }
        else {
            std::cout << "Squid::createMaterialBudget(): " << err_no_matfile << std::endl;
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
     * @param matfile The name and - if necessary - path of the materials configuration file
     * @param usher_verbose A flag that turns the final status summary of the inactive surface placement algorithm on or off
     * @param mat_verbose A flag that turns the final status summary of the material budget on or off
     * @return True if there were no errors during processing, false otherwise
     */
  bool Squid::buildFullSystem(std::string geomfile, std::string settingsfile, std::string matfile, bool usher_verbose, bool mat_verbose) {
    std::string dummyString;
    return buildFullSystem(geomfile, settingsfile, matfile, dummyString, usher_verbose, mat_verbose);
  }
 
    /**
     * Build a full system consisting of tracker object, collection of inactive surfaces and material budget from the
     * given configuration files. All three objects replace the previously registered ones, if they existed. They remain
     * in the squid  until they are overwritten by a second call to this function or by another one that creates new
     * instances of them.
     * @param geomfile The name and - if necessary - path of the geometry configuration file
     * @param settingsfile The name and - if necessary - path of the module settings configuration file
     * @param matfile The name and - if necessary - path of the materials configuration file
     * @param pixmatfile The name and - if necessary - path of the materials configuration file for the pixels
     * @param usher_verbose A flag that turns the final status summary of the inactive surface placement algorithm on or off
     * @param mat_verbose A flag that turns the final status summary of the material budget on or off
     * @return True if there were no errors during processing, false otherwise
     */
  bool Squid::buildFullSystem(std::string geomfile, std::string settingsfile, std::string matfile, std::string& pixmatfile, bool usher_verbose, bool mat_verbose) {
    if (buildInactiveSurfaces(geomfile, settingsfile, usher_verbose)) return createMaterialBudget(matfile, pixmatfile, mat_verbose);
        return false;
    }
    
    /**
     * Analyse the previously created full system, writing the full series of output files: HTML for the histograms that
     * were filled during analysis of the material budget, ROOT for the geometry visualisation and plain text for the
     * feeder/neighbour relations. The tracker object, the collection of inactive surfaces and the material budget must
     * all exist already for this function to succeed.
     *
     * WARNING: do <i>not</i> turn the <i>simplified</i> flag off unless you have a very small number of modules,
     * or another very good reason! In most cases, loading a geometry with individual modules into one of the geometry
     * viewers will simply crash ROOT because the number of volumes is too large. So once again - use with caution!
     *
     * @param htmlout The name - without path - of the designated HTML output file
     * @param rootout The name - without path - of the designated ROOT output file
     * @param graphout The name - without path - of the designated plain text output file
     * @param tracks The number of tracks that should be fanned out across the analysed region
     * @param simplified A flag that turns bounding boxes instead of individual modules in the ROOT output file on or off
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeFullSystem(std::string htmlout, std::string rootout, std::string graphout, int tracks, bool simplified) {
        if (analyzeGeometry(rootout, graphout, simplified)) return analyzeMaterialBudget(htmlout, tracks);
        return false;
    }
    
    /**
     * Analyse the previously created full system, writing the result to an HTML file and the geometry
     * visualisation to a ROOT file.
     * @param htmlout The name - without path - of the designated HTML output file
     * @param rootout The name - without path - of the designated ROOT output file
     * @param tracks The number of tracks that should be fanned out across the analysed region
     * @param simplified A flag that turns bounding boxes instead of individual modules in the ROOT output file on or off
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeGeoMat(std::string htmlout, std::string rootout, int tracks, bool simplified) {
        if (analyzeGeometry(rootout, simplified)) return analyzeMaterialBudget(htmlout, tracks);
        return false;
    }
    
    /**
     * Analyse the previously created full system, writing the result to an HTML file and the feeder/neighbour
     * graph to a plain text file.
     * @param htmlout The name - without path - of the designated HTML output file
     * @param graphout The name - without path - of the designated plain text output file
     * @param tracks The number of tracks that should be fanned out across the analysed region
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeGraphMat(std::string htmlout, std::string graphout, int tracks) {
        if (analyzeNeighbours(graphout)) return analyzeMaterialBudget(htmlout, tracks);
        return false;
    }
    
    /**
     * Build a ROOT representation of a partial or complete tracker geometry that can be visualised
     * in a ROOT viewer later. In addition, build the feeder/neighbour graph of the collection of inactive
     *  surfaces if it exists, and save both results to file. This function succeeds if either the tracker or
     * both the tracker and the inactive surfaces exist, but the graph file will only be created for a full
     * geometry.
     * @param rootout The name - without path - of the designated ROOT output file
     * @param graphout The name - without path - of the designated plain text output file
     * @param simplified A flag that turns bounding boxes instead of individual modules in the ROOT output file on or off
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeGeometry(std::string rootout, std::string graphout, bool simplified) {
        if (tr) {
            if (is) {
                if (v.needsReset()) resetVizard();
                v.display(*tr, *is, rootout, simplified);
                v.writeNeighbourGraph(*is, graphout);
            }
            else {
                std::cout << "Squid::analyzeGeometry(): " << warning_rootonly << std::endl;
                is = new InactiveSurfaces();
                if (v.needsReset()) resetVizard();
                v.display(*tr, *is, rootout, simplified);
                delete is;
                is = NULL;
            }
            return true;
        }
        else {
            std::cout << "Squid::analyzeGeometry(): " << err_no_tracker << std::endl;
            return false;
        }
    }
    
    /**
     * Build a ROOT representation of a partial or complete tracker geometry that can be visualised
     * in a ROOT viewer later, and save the result to a ROOT file. This function succeeds if either the
     * tracker or both the tracker and the inactive surfaces exist.
     * @param rootout The name - without path - of the designated output file
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeGeometry(std::string rootout, bool simplified) {
        if (tr) {
            if (is) v.display(*tr, *is, rootout, simplified);
            else {
                is = new InactiveSurfaces();
                if (v.needsReset()) resetVizard();
                v.display(*tr, *is, rootout, simplified);
                delete is;
                is = NULL;
            }
            return true;
        }
        else {
            std::cout << "Squid::analyzeGeometry(): " << err_no_tracker << std::endl;
            return false;
        }
    }
    
    /**
     * Build the feeder/neighbour graph of the previously created collection of inactive surfaces and
     * save the results in a plain text file.
     * @param graphout The name - without path - of the designated output file
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeNeighbours(std::string graphout) {
        if (is) {
            v.writeNeighbourGraph(*is, graphout);
            return true;
        }
        else {
            std::cout << "Squid::analyzeNeighbours(): " << err_no_inacsurf << std::endl;
            return false;
        }
    }
    
    /**
     * Analyze the previously created material budget and save the results in an HTML file. If a material budget
     * for an associated pixel detector exists, it is included in the calculations.
     * @param htmlout The name - without path - of the designated output file
     * @param tracks The number of tracks that should be fanned out across the analysed region
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeMaterialBudget(std::string htmlout, int tracks) {
        if (mb) {
            a.analyzeMaterialBudget(*mb, mainConfiguration.getMomenta(), tracks, pm);
            v.histogramSummary(a, htmlout);
            return true;
        }
        else {
            std::cout << "Squid::analyzeMaterialBudget(): " << err_no_matbudget << std::endl;
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
    
    /**
     * Create the geometry summary page for an existing tracker. Unless specified otherwise in an <i>Output</i>
     * block in the geometry file, the page will be written to a subfolder based on the tracker name.
     * @configFileName The name - preferably without the path - of the geometry configuration file
     * @dressFileName The name - preferably without the path - of the module settings configuration file
     * @return True if there is an existing tracker to summarise, false otherwise.
     */
  /*    bool Squid::trackerSummary(std::string configFileName, std::string dressFileName) {
        if (tr) {
            std::string myDirectory;
            std::string destConfigFile;
            std::string destDressFile;
            // Optical transmission
            tr->createGeometry(true);
            tr->computeBandwidth();
            
            // Analysis
            tr->analyze(2000);
            
            // Summary and save
            tr->writeSummary(true, extractFileName(configFileName), extractFileName(dressFileName));
            //tr->printBarrelModuleZ();
            tr->save();
            
            myDirectory = tr->getActiveDirectory();
            destConfigFile = myDirectory + "/" + extractFileName(configFileName);
            destDressFile = myDirectory + "/" + extractFileName(dressFileName);
            
            bfs::remove(destConfigFile);
            bfs::remove(destDressFile);
            bfs::copy_file(configFileName, destConfigFile);
            bfs::copy_file(dressFileName, destDressFile);
            return true;
        }
        else {
            std::cout << "Squid::trackerSummary(): " << err_no_tracker << std::endl;
            return false;
        }
	} */
    
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
    
#ifdef USING_ROOTWEB
    /**
     * Prepare the website object (if not done yet) from the configuration file
     * it needs the tracker object to be already there
     * @return a boolean with the operation success
     */
    bool Squid::prepareWebsite() {
        if (sitePrepared) return true;
        string trackerName;
        if (tr) trackerName = tr->getName();
        else trackerName = default_trackername;
        string layoutDirectory;
        //styleDirectory=mainConfiguration.getStyleDirectory();
        layoutDirectory=mainConfiguration.getLayoutDirectory();
        layoutDirectory+="/"+trackerName;
        if (layoutDirectory!="") site.setTargetDirectory(layoutDirectory);
        else return false;
        site.setTitle(trackerName);
        site.setComment("layouts");
        site.setCommentLink("../");
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
        if (!prepareWebsite()) return false;
        if (addLogPage) {
            v.makeLogPage(site);
        }
        
        return site.makeSite();
    }
    
    /**
     * Analyze the previously created material budget and save the results through rootweb. If a material budget
     * for an associated pixel detector exists, it will be included in the calculations as well.
     * @param tracks The number of tracks that should be fanned out across the analysed region
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeMaterialBudgetSite(int tracks) {
        if (mb) {
            a.analyzeMaterialBudget(*mb, mainConfiguration.getMomenta(), tracks, pm);
            v.histogramSummary(a, site, "outer");
            if (pm) {
	      // TODO: make this much neater!
	      Analyzer* pixelAnalyzer = new Analyzer;
	      pixelAnalyzer->analyzeMaterialBudget(*pm, mainConfiguration.getMomenta(), tracks);
	      v.histogramSummary(*pixelAnalyzer, site, "pixel");
            }
	    a.computeWeightSummary(*mb);
	    v.weigthSummart(a, site, "outer");
            v.errorSummary(a, site);
            return true;
        }
        else {
            std::cout << "Squid::analyzeMaterialBudgetSite(): " << err_no_matbudget << std::endl;
            return false;
        }
    }
    
    /**
     * Analyze the previously created geometry and produces the html ouput through rootweb. If a material budget
     * for an associated pixel detector exists, it will be included in the calculations as well.
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeGeometrySite(int tracks /* = 1000 */ ) {
        if (tr) {
            a.analyzeGeometry(*tr, tracks);
            a.createGeometryLite(*tr);
            v.geometrySummary(a, *tr, site);
            a.computeBandwidth(*tr);
            v.bandwidthSummary(a, *tr, site);
            return true; // TODO: is not really meaningful
        } else {
            std::cout << "Squid::analyzeGeometrySite(): " << err_no_tracker << std::endl;
            return false;
        }
    }

  /**
   * Analyze the previously created geometry and without no output  through rootweb.
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::pureAnalyzeGeometry(int tracks) {
    if (tr) {
      a.analyzeGeometry(*tr, tracks);
      a.createGeometryLite(*tr);
      a.computeBandwidth(*tr);
      return true; // TODO: is not really meaningful
    } else {
      std::cout << "Squid::pureAnalyzeGeometry(): " << err_no_tracker << std::endl;
      return false;
    }
  }

  /**
   * Analyze the previously created material budget with no output.
   * @param tracks The number of tracks that should be fanned out across the analysed region
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::pureAnalyzeMaterialBudget(int tracks) {
    if (mb) {
      a.analyzeMaterialBudget(*mb, mainConfiguration.getMomenta(), tracks, pm);
      if (pm) {
	// TODO: make this much neater!
	if (pixelAnalyzer) delete pixelAnalyzer;
	pixelAnalyzer = new Analyzer;
	pixelAnalyzer->analyzeMaterialBudget(*pm, mainConfiguration.getMomenta(), tracks);
      }
      a.computeWeightSummary(*mb);
      return true;
    } else {
      std::cout << "Squid::pureAnalyzeMaterialBudget(): " << err_no_matbudget << std::endl;
      return false;
    }
  }

  /**
   * Produces the output of the analysis of the geomerty analysis
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::reportGeometrySite() {
    if (tr) {
      v.geometrySummary(a, *tr, site);
      v.bandwidthSummary(a, *tr, site);
      return true;
    } else {
      std::cout << "Squid::reportGeometrySite(): " << err_no_tracker << std::endl;
      return false;
    }
  }

  /**
   * Produces the output of the analysis of the material budget analysis
   * @return True if there were no errors during processing, false otherwise
   */
  bool Squid::reportMaterialBudgetSite() {
    if (mb) {
      v.histogramSummary(a, site, "outer");
      if ((pm)&&(pixelAnalyzer)) {
	v.histogramSummary(*pixelAnalyzer, site, "pixel");
      }
      v.weigthSummart(a, site, "outer");
      v.errorSummary(a, site);
      return true;
    }
    else {
      std::cout << "Squid::reportMaterialBudgetSite(): " << err_no_matbudget << std::endl;
      return false;
    }
  }
    
  bool Squid::additionalInfoSite(std::string& geomfile, std::string& settingsfile, std::string& matfile, std::string& pixmatfile) {
        if (!tr) {
            std::cout << "Squid::additionalInfoSite(): " << err_no_tracker << std::endl;
            return false;
        } else {
            v.additionalInfoSite(geomfile, settingsfile, matfile, pixmatfile, a, *tr, site);
            return true;
        }
    }
    
    
#endif
    
}

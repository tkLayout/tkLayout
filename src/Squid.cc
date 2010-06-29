/**
 * @file Squid.cc
 * @brief
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
    }
    
    /**
     * The destructor deletes the heap-allocated internal objects if they exist.
     */
    Squid::~Squid() {
        if (mb) delete mb;
        if (is) delete is;
        if (tr) delete tr;
    }
    
    /**
     * Build a bare-bones geometry of active modules. The resulting tracker object replaces the previously
     * registered one, if such an object existed. It remains in the squid until it is overwritten by a second call
     * to this function or by another one that creates a new tracker object.
     * @param geomfile The name and - if necessary - path of the geometry configuration file
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::buildTracker(std::string geomfile) {
        if (tr) delete tr;
        tr = cp.parseFile(geomfile);
        if (tr) {
            g = geomfile;
            return true;
        }
        return false;
    }
    
    /**
     * Dress the previously created geometry with module options. The resulting tracker object replaces the
     * previously registered one, if such an object existed. It remains in the squid until it is overwritten by a
     * second call to this function or by another one that creates a new tracker object.
     * @param settingsfile The name and - if necessary - path of the module settings configuration file
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::dressTracker(std::string settingsfile) {
        if (tr) {
            cp.dressTracker(tr, settingsfile);
            return true;
        }
        else {
            std::cout << "Squid::dressTracker(): " << err_no_tracker << std::endl;
            return false;
        }
    }
    
    /**
     * Build a geometry of active modules using both geometry constraints and module settings.
     * @param geomfile The name and - if necessary - path of the geometry configuration file
     * @param settingsfile The name and - if necessary - path of the module settings configuration file
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::buildTrackerSystem(std::string geomfile, std::string settingsfile) {
        if (tr) delete tr;
        if ((tr = cp.parseFile(geomfile))) {
            cp.dressTracker(tr, settingsfile);
            g = geomfile;
            return true;
        }
        else {
            g = "";
            return false;
        }
    }
    
    /**
     * Build up the inactive surfaces around the previously created tracker geometry. The resulting collection
     * of inactive surfaces replaces the previously registered one, if such an object existed. It remains in the
     * squid until it is overwritten by a second call to this function or by another one that creates a new
     * collection of inactive surfaces.
     * @param verbose A flag that turns the final status summary of the inactive surface placement algorithm on or off
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::buildInactiveSurfaces(bool verbose) {
        if (!g.empty()) {
            if (tr) {
                if (is) delete is;
                is = new InactiveSurfaces();
                u.arrange(*tr, *is, g, verbose);
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
     * function or by another one that creates a new tracker object and collection of inactive surfaces.
     * @param geomfile The name and - if necessary - path of the geometry configuration file
     * @param verbose A flag that turns the final status summary of the inactive surface placement algorithm on or off
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::buildInactiveSurfaces(std::string geomfile, bool verbose) {
        if (is) delete is;
        if (buildTracker(geomfile)) {
            is = new InactiveSurfaces();
            u.arrange(*tr, *is, geomfile, verbose);
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
     * tracker object and collection of inactive surfaces.
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
     * both the tracker and the inactive surfaces exist.
     * @param matfile The name and - if necessary - path of the materials configuration file
     * @param verbose A flag that turns the final status summary of the material budget on or off
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::createMaterialBudget(std::string matfile, bool verbose) {
        if (fileExists(matfile)) {
            if (tr) {
                if (!is) is = new InactiveSurfaces();
                if (mb) delete mb;
                mb  = new MaterialBudget(*tr, *is);
                if (c.initDone()) c.reset();
                if (mp.initMatCalc(matfile, c)) {
                    mb->materialsAll(c);
                    if (verbose) mb->print();
                    return true;
                }
                else {
                    delete mb;
                    mb = NULL;
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
        if (buildInactiveSurfaces(geomfile, settingsfile, usher_verbose)) return createMaterialBudget(matfile, mat_verbose);
        return false;
    }
    
    /**
     * Analyse the previously created full system, writing to the full series of output files: HTML for the histograms
     * that were filled during analysis of the material budget, ROOT for the geometry visualisation and plain text for
     * the feeder/neighbour relations. The tracker object, the collection of inactive surfaces and the material budget
     * must all exist already for this function to succeed.
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
     * in a ROOT viewer later as well as the feeder/neighbour graph of the collection of inactive surfaces
     *  if it exists, and save both results to file. This function succeeds if either the tracker or both the
     * tracker and the inactive surfaces exist, but the graph file will only be created for a full geometry.
     * @param rootout The name - without path - of the designated ROOT output file
     * @param graphout The name - without path - of the designated plain text output file
     * @param simplified A flag that turns bounding boxes instead of individual modules in the ROOT output file on or off
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeGeometry(std::string rootout, std::string graphout, bool simplified) {
        if (tr) {
            if (is) {
                v.display(*tr, *is, rootout, simplified);
                v.writeNeighbourGraph(*is, graphout);
            }
            else {
                std::cout << "Squid::analyzeGeometry(): " << warning_rootonly << std::endl;
                is = new InactiveSurfaces();
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
     * in a ROOT viewer later and save the result to a ROOT file. This function succeeds if either the
     * tracker or both the tracker and the inactive surfaces exist.
     * @param rootout The name - without path - of the designated output file
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeGeometry(std::string rootout, bool simplified) {
        if (tr) {
            if (is) v.display(*tr, *is, rootout, simplified);
            else {
                is = new InactiveSurfaces();
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
     * Analyze the previously created material budget and save the results in an HTML file.
     * @param htmlout The name - without path - of the designated output file
     * @param tracks The number of tracks that should be fanned out across the analysed region
     * @return True if there were no errors during processing, false otherwise
     */
    bool Squid::analyzeMaterialBudget(std::string htmlout, int tracks) {
        if (mb) {
            a.analyzeMaterialBudget(*mb, tracks);
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
    bool Squid::translateFullSystemToXML(std::string xmlout) {
        if (mb) {
            t2c.translate(c.getMaterialTable(), *mb, xmlout);
            return true;
        }
        else {
            std::cout << "Squid::translateFullSystemToXML(): " << err_no_matbudget << std::endl;
            return false;
        }
    }
    
    bool Squid::trackerSummary(std::string configFileName, std::string dressFileName) {
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
    
    std::string Squid::extractFileName(const std::string& full) {
        std::string::size_type idx = full.find_last_of("/");
        if (idx != std::string::npos)
            return full.substr(idx+1);
        else return full;
    }
}

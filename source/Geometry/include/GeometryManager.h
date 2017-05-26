/*
 * GeometryManager.h
 *
 *  Created on: 27. 4. 2016
 *      Author: Drasal (CERN)
 */

#ifndef INCLUDE_GEOMETRYMANAGER_H_
#define INCLUDE_GEOMETRYMANAGER_H_

#include <boost/property_tree/ptree_fwd.hpp>
#include <memory>
#include <set>
#include <string>
#include <vector>

// Forward declarations
class Detector;

/*
 * @class GeometryManager
 * @brief The core geometry class, building the overall tracker with its inactive components & beam-pipe, so called Detector.
 * @details The core geometry class, building the overall tracker with its inactive components & beam-pipe, so called Detector.
 * The tracker consists of individual sub-trackers and corresponding inactive parts. In addition various support structures,
 * services, etc. are built as well. As an interface for analysis, it provides getDetector() methods to obtain geometrical info:
 * built active, trackers, their passive parts (passive trackers), supports etc. This class takes over the role of previous
 * Squid class and its geometry related methods ...
 */
class GeometryManager {

 public:

  //! Constructor - initializes manager -> reads-in all geometry configuration files (using the mainConfigHandler)
  //! @param[in]  baseGeomFile Base geometry file with fill address
  GeometryManager(std::string baseGeomFile);

  //! Destructor - clean out the memory
  ~GeometryManager();

  //! Build active part of the detector: individual trackers: pixel, strip, central, forward, ... (a bare-bones geometry of active modules).
  //! @return True if there were no errors during processing, false otherwise
  bool buildActiveDetector();

  //! Build all passive components related to individual active sub-trackers and/or beam-pipe.
  //! @return True if there were no errors during processing, false otherwise
  bool buildPassiveDetector();

  //! Get detector - container for active/passive parts of trackers, beam-pipe etc.
  const Detector& getDetector() const {return *m_detector;}

  //! Get geometry layout name
  //! @return layout name
  std::string getLayoutName() const {return m_layoutName;}

  //! Get default html directory, where all results are saved
  //! @return directory address
  std::string getWebDir() const {return m_htmlDir;}

  //! Get name of main geometry config file
  //! @return file name
  std::string getBaseGeomFileName() const {return m_geomFile;}

  //! Get path of tkLayout run directory
  //! @return run directory path
  std::string getRunDirPath() const {return m_baseDir;}

  //! Get list of all configuration files obtained from the base geometry file using @include command
  //! @return set of all include files (strings)
  std::set<std::string> getListOfConfFiles() const {return m_includeSet;}

 private:

  //! Set base geometry file, directory, layout name, ...
  void setGeometryFile(std::string geomFile);

  bool        m_initOK;      //!< Initialization OK

  std::string m_geomFile;    //!< Name of main geometry config file
  std::string m_baseDir;     //!< Base tkLayout run directory
  std::string m_htmlDir;     //!< Default html directory, where all results are saved
  std::string m_layoutName;  //!< Geometry layout name

  std::unique_ptr<boost::property_tree::ptree> m_geomTree;   //!< Geometry property tree (built from the tree structure of geometry configuration files)
  std::set<std::string>                        m_includeSet; //!< List of all configuration files obtained from the base geometry file using @include command

  std::unique_ptr<Detector>                    m_detector;   //!< The overall detector geometry: trackers, beam-pipe, etc.

  // Constants
  const std::string c_defaultHtmlDir = "results";

}; // Class

#endif /* INCLUDE_GEOMETRYMANAGER_H_ */

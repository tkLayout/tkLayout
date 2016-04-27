/*
 * GeometryManager.h
 *
 *  Created on: 27. 4. 2016
 *      Author: Drasal (CERN)
 */

#ifndef INCLUDE_GEOMETRYMANAGER_H_
#define INCLUDE_GEOMETRYMANAGER_H_

#include <vector>
#include <set>
#include <string>
#include <boost/property_tree/ptree_fwd.hpp>

class Tracker;
namespace insur {
  class InactiveSurfaces;
}
class Support;

/*
 * The core geometry class, building the overall tracker with its inactive components. The tracker consists of individual sub-trackers and corresponding
 * inactive parts. In addition various support structures, services, etc. are built as well. For the analysis, it provides several get methods to provide
 * a vector of built active trackers, their pasive parts (pasive trackers) etc. This class takes over the role of previous Squid methods related to geometry.
 */
class GeometryManager {

 public:

  //! Constructor - initializes manager -> reads-in all geometry configuration files (using the mainConfigHandler)
  //! @param[in]  baseGeomFile Base geometry file with fill address
  GeometryManager(std::string baseGeomFile);

  //! Destructor - clean out the memory
  ~GeometryManager();

  //! Build active part of the tracker (a bare-bones geometry of active modules). The resulting tracker object consists of individual sub-trackers
  //! (pixel, strip, central, forward, ...). This procedure replaces the previously registered tracker (applying correct memory managment), if such an object
  //! existed.
  //! @return True if there were no errors during processing, false otherwise
  bool buildActiveTracker();

  //! Build tracker support structures, which are independent on individual sub-trackers (the supports directly related to sub-trackers are built as pasive
  //! components of active sub-tracker). This procedure replaces the previously registered support (applying correct memory managment), if such an object
  //! existed.
  //! @return True if there were no errors during processing, false otherwise
  bool buildTrackerSupport();

  //! Build all pasive components related to individual active sub-trackers. This procedure replaces the previously registered support (applying correct
  //! memory managment), if such an object existed.
  //! @return True if there were no errors during processing, false otherwise
  bool buildPasiveTracker();

  //! Get active sub-trackers
  //! @return std::vector<Tracker*>
  std::vector<const Tracker*> getActiveTrackers() const;

  //! Get tracker supports, which are independent on individual sub-trackers
  //! @return std::vector<Support*>
  std::vector<const Support*> getTrackerSupports() const;

  //! Get pasive components related to active sub-trackers
  //! @return std::vector<InactiveSurfaces*>
  std::vector<const insur::InactiveSurfaces*> getPasiveTrackers() const;

 private:

  //! Set base geometry file, directory, layout name, ...
  void setGeometryFile(std::string geomFile);

  bool        m_initOK;      //!< Initialization OK

  std::string m_geomFile;    //!< Base geometry file name
  std::string m_geomDir;     //!< Base geometry directory
  std::string m_htmlDir;     //!< Default html directory, where all results are saved
  std::string m_layoutName;  //!< Geometry layout name

  boost::property_tree::ptree*   m_geomTree;       //!< Geometry property tree (built from the tree structure of geometry configuration files)
  std::set<std::string>          m_includeSet;     //!< List of all configuration files obtained from the base geometry file using @include command

  std::vector<Tracker*>                 m_activeTrackers; //!< Vector of active sub-trackers
  std::vector<Support*>                 m_supports;       //!< Vector of independent support structures not directly related to active sub-trackers
  std::vector<insur::InactiveSurfaces*> m_pasiveTrackers; //!< Vector of pasive sub-trackers

  // Constants
  const std::string c_defaultHtmlDir = "results";

}; // Class

#endif /* INCLUDE_GEOMETRYMANAGER_H_ */

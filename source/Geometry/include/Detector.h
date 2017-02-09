/*
 * Detector.h
 *
 *  Created on: 25. 7. 2016
 *      Author: Z.Drasal (CERN)
 */
#ifndef INCLUDE_DETECTOR_H_
#define INCLUDE_DETECTOR_H_

#include <memory>
#include <string>
#include <vector>
#include <boost/property_tree/ptree_fwd.hpp>

#include "Visitable.h"

// Forward declaration
class BeamPipe;
class Tracker;

/*
 * @class Detector
 * @details Detector class holds all information about individual trackers systems & beam-pipe.
 * It represents an object at the top of geometrical hierarchy to be passed by reference from
 * GeometryManager to AnalysisManager, etc.
 */
class Detector : public Visitable {

 public:

  //! Default constructor
  Detector(std::string name);

  //! Default destructor
  ~Detector();

  //! Build active parts of the tracker (a bare-bones geometry of active modules). The resulting tracker object consists of individual sub-trackers
  //! (pixel, strip, central, forward, ...). This procedure replaces the previously registered tracker (applying correct memory managment), if such an object
  //! existed.
  //! @return True if there were no errors during processing, false otherwise
  bool buildActiveTracker(boost::property_tree::ptree& geomTree);

  //! Build & route all passive components (as related to individual active sub-trackers) & calculate material assigned to them.
  //! @return True if there were no errors during processing, false otherwise
  bool buildPassiveTracker();

  //! Build beam pipe. This procedure replaces the previously registered beam pipe
  //! @return True if there were no errors during processing, false otherwise
  bool buildBeamPipe(boost::property_tree::ptree& geomTree);

  //! GeometryVisitor pattern -> detector visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> detector visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  //! Get references to all sub-trackers
  //! @return vector of pointers to active trackers
  std::vector<const Tracker*> getTrackers() const;

  //! Get reference to beam pipe
  //! @return InactiveTube
  const BeamPipe* getBeamPipe() const;

 private:

  std::string m_name; //!< Detector configuration name

  std::vector<std::unique_ptr<Tracker>> m_trackers; //!< Vector of sub-trackers
  std::unique_ptr<BeamPipe>             m_beamPipe; //!< Passive surface (tube) simulating beam pipe

}; // Class


#endif /* INCLUDE_DETECTOR_H_ */

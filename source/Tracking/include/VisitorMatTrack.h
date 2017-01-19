/*
 * VisitorMatTrack.h
 *
 *  Created on: 17. 1. 2017
 *      Author: drasal
 */

#ifndef INCLUDE_VISITORMATTRACK_H_
#define INCLUDE_VISITORMATTRACK_H_

#include <Visitor.h>

class Track;
namespace insur {
  class InactiveElement;
}

//! Visitor pattern class - checks that modules, beam-pipe etc. hit by a track -> assigns all active or passive hits to the track
class VisitorMatTrack : public ConstGeometryVisitor {

public:

  //! Constructor
  //! As a parameter provide a track, under which the material is studied
  VisitorMatTrack(Track& matTrack);

  //! Destructor
  ~VisitorMatTrack();

  //! Visit BeamPipe
  void visit(const BeamPipe& bp) override;

  //! Visit Barrel
  void visit(const Barrel& b) override;

  //! Visit BarrelModule (no limits on Rods, Layers or Barrels)
  void visit(const BarrelModule& m) override;

  //! Visit EndcapModule (no limits on Rings, Discs or Endcaps)
  void visit(const EndcapModule& m) override;

  //! Visit Support structure
  void visit(const SupportStructure& s) override;

private:

  //! Analyze if module crossed by given track & how much material is on the way
  void analyzeModuleMB(const DetectorModule& m);

  //! Helper method - analyse inactive element & estimate how much material is in the way
  void analyzeInactiveElement(const insur::InactiveElement& e);

  Track& m_matTrack;   //!< Shooting direction + origin encapsulated in a track class -> update track with hits once found

}; // Class


#endif /* INCLUDE_VISITORMATTRACK_H_ */

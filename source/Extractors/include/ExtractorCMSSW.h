/*
 * ExtractorCMSSW.h
 *
 *  Created on: 1. 9. 2016
 */

#ifndef INCLUDE_EXTRACTORCMSSW_H_
#define INCLUDE_EXTRACTORCMSSW_H_

// System libraries
#include <string>
#include <vector>
#include <map>
#include <memory>

#include <Visitor.h>
#include <AnalyzerUnit.h>

/*
 * @class ExtractorCMSSW
 * @brief A specialized analysis module extracting full tkLayout geometry description for use in CMSSW.
 * @details TODO:
 */
class ExtractorCMSSW : public AnalyzerUnit {

 public:

  //! Constructor
  ExtractorCMSSW(const Detector& detector);

  //! Destructor
  virtual ~ExtractorCMSSW();

  //! Initialize - mostly histograms & other containers, nGeomTracks is assumed to be zero
  //! @return True if OK
  bool init(int nGeomTracks=0);

  //! Gather information about the geometry layout (if init OK) & output the data to XML file
  //! @return True if OK
  bool analyze();

  //! Visualization method remains non-used in case of geometry extraction
  //! @return True if OK
  bool visualize(RootWSite& webSite) {m_isVisOK=true; return m_isVisOK;}

 private:

}; // Class

#endif /* INCLUDE_EXTRACTORCMSSW_H_ */

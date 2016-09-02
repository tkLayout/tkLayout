/*
 * ExtractorFCCSW.h
 *
 *  Created on: 1. 9. 2016
 *      Author: Z.Drasal (CERN)
 */

#ifndef INCLUDE_EXTRACTORFCCSW_H_
#define INCLUDE_EXTRACTORFCCSW_H_

// System libraries
#include <string>
#include <vector>
#include <map>
#include <memory>

#include <Visitor.h>
#include <AnalyzerUnit.h>

// Forward declaration
namespace tinyxml2
{
  class XMLDocument;
  class XMLNode;
}

/*
 * @class ExtractorFCCSW
 * @brief A specialized analysis module extracting full tkLayout geometry description to an XML file, for use in FCCSW.
 * @details The FCCSW requires DD4Hep objects as basic geometry building blocks. To provide an interface between tkLayout geometry description
 * and FCCSW geometry description (DD4Hep) a specialized analysis module, ExtractorFCCSW, analyses the full geometry using a Visitor pattern
 * and outputs experiment related description to an XML file (through TinyXML library). XML file is supposed to be used by DD4Hep creator to
 * build the same tracker geometry as optimized by tkLayout for full/fast simulation studies within FCCSW.
 */
class ExtractorFCCSW : public AnalyzerUnit {

 public:

  //! Constructor
  ExtractorFCCSW(const Detector& detector);

  //! Destructor
  virtual ~ExtractorFCCSW();

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

  std::unique_ptr<tinyxml2::XMLDocument> m_xmlDoc;      //!< An XML base container (pointing to XML file created on disc)
  tinyxml2::XMLNode*                     m_xmlNodeRoot; //!< XML root node

}; // Class

#endif /* INCLUDE_EXTRACTORFCCSW_H_ */

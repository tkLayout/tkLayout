/*
 * ExtractorFCCSW.cc
 *
 *  Created on: 1. 9. 2016
 *      Author: Z.Drasal (CERN)
 */
#include "ExtractorFCCSW.h"

// Include files
#include "tinyxml2.h"
#include "MainConfigHandler.h"
#include "MessageLogger.h"
#include "SimParms.h"

// Used namespaces
using namespace tinyxml2;

//
// Constructor
//
ExtractorFCCSW::ExtractorFCCSW(const Detector& detector) :
 AnalyzerUnit("ExtractorFCCSW", detector)
{}

//
// Destructor
//
ExtractorFCCSW::~ExtractorFCCSW()
{}

//
// Initialize - mostly histograms & other containers
// @return True if OK
bool ExtractorFCCSW::init(int nGeomTracks)
{
  // Create an empty XML document
  m_xmlDoc = std::unique_ptr<XMLDocument>(new XMLDocument());

  // Create XML declaration
  XMLNode* xmlDeclare = m_xmlDoc->NewDeclaration("xml version=\"1.0\" encoding=\"UTF-8\"");
  m_xmlDoc->InsertFirstChild(xmlDeclare);

  // Create DD4Hep standard Root node: lcdd
  XMLNode* m_xmlNodeRoot = m_xmlDoc->NewElement("lcdd");
  m_xmlNodeRoot->ToElement()->SetAttribute("xmlns:compact","http://www.lcsim.org/schemas/compact/1.0");
  m_xmlNodeRoot->ToElement()->SetAttribute("xmlns:xs","http://www.w3.org/2001/XMLSchema");
  m_xmlNodeRoot->ToElement()->SetAttribute("xs:noNamespaceSchemaLocation","http://www.lcsim.org/schemas/compact/1.0/compact.xsd");

  m_xmlDoc->InsertAfterChild(xmlDeclare, m_xmlNodeRoot);

  // Add info to Root node
  auto xmlInfo = m_xmlDoc->NewElement("info");
  xmlInfo->SetAttribute("name", "TkLayoutTracker");
  xmlInfo->SetAttribute("title", "TkLayoutTracker");
  xmlInfo->SetAttribute("author", MainConfigHandler::getInstance().getResultsAuthor().c_str());
  xmlInfo->SetAttribute("url", "http://fcc-tklayout.web.cern.ch/fcc-tklayout");
  xmlInfo->SetAttribute("status", "optimization");
  xmlInfo->SetAttribute("version", SimParms::getInstance().getLayoutName().c_str());

  // Add comment to info element
  auto xmlInfoComment = m_xmlDoc->NewElement("comment");
  xmlInfoComment->SetText(std::string("The tracker geometry as designed and optimized by tkLayout software - design version: "+SimParms::getInstance().getLayoutName()).c_str());
  xmlInfo->InsertEndChild(xmlInfoComment);

  m_xmlNodeRoot->InsertEndChild(xmlInfo);

  m_isInitOK = true;
  return m_isInitOK;
}

//
// Gather information about the geometry layout (if init OK) & output the data to XML file
// @return True if OK
bool ExtractorFCCSW::analyze()
{

  // Save XML file
  std::string xmlDocName = SimParms::getInstance().getWebDir()+"/"+SimParms::getInstance().getLayoutName()+".xml";
  auto eResult = m_xmlDoc->SaveFile(xmlDocName.c_str());

  if (eResult==XML_SUCCESS) {

    logINFO(std::string("Geometry extracted to the following XML file: "+xmlDocName));
    m_isAnalysisOK = true;
  }
  else {

    logERROR("ExtractorFCCSW::analyze(): A problem occured when saving an XML file!");
    m_isAnalysisOK = false;
  }

  return m_isAnalysisOK;
}

/*
 * ExtractorFCCSW.cc
 *
 *  Created on: 1. 9. 2016
 *      Author: Z.Drasal (CERN)
 */
#include "ExtractorFCCSW.h"

// System include files
#include <boost/range/adaptor/reversed.hpp>
#include <map>

// Project include files
#include "Barrel.h"
#include "Endcap.h"
#include "MainConfigHandler.h"
#include "MaterialTab.h"
#include "MessageLogger.h"
#include "SimParms.h"
#include "tinyxml2.h"
#include "Tracker.h"
#include "Units.h"

// Used namespaces
using namespace tinyxml2;

//
// Constructor
//
ExtractorFCCSW::ExtractorFCCSW(const Detector& detector) :
 AnalyzerUnit("ExtractorFCCSW", detector),
 m_xmlNodeRoot(nullptr)
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
  m_xmlDefinitionsDoc = std::unique_ptr<XMLDocument>(new XMLDocument());

  // Create XML declaration
  XMLNode* xmlDeclare = m_xmlDoc->NewDeclaration("xml version=\"1.0\" encoding=\"UTF-8\"");
  m_xmlDoc->InsertFirstChild(xmlDeclare);

  // Create DD4Hep standard Root node: lcdd
  m_xmlNodeRoot = m_xmlDoc->NewElement("lccdd");
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
  // workaround: hardcoded ids for the names
  std::map<std::string, std::string> detNames2IDs;
  detNames2IDs.insert(std::make_pair("InnerBRL_0", "BarTrackerInner_id"));
  detNames2IDs.insert(std::make_pair("InnerBRL_1", "BarTrackerInner_id + 16"));
  detNames2IDs.insert(std::make_pair("InnerECAPpos", "EndCapTrackerInner_id"));
  detNames2IDs.insert(std::make_pair("InnerECAPneg", "EndCapTrackerInner_id + 16"));
  detNames2IDs.insert(std::make_pair("OuterBRL", "BarTrackerOuter_id"));
  detNames2IDs.insert(std::make_pair("OuterECAPpos", "EndCapTrackerOuter_id"));
  detNames2IDs.insert(std::make_pair("OuterECAPneg", "EndCapTrackerOuter_id + 16"));
  detNames2IDs.insert(std::make_pair("FwdECAPpos", "27"));
  detNames2IDs.insert(std::make_pair("FwdECAPneg", "28"));
  detNames2IDs.insert(std::make_pair("FwdHelpECAPpos", "29"));
  detNames2IDs.insert(std::make_pair("FwdHelpECAPneg", "30"));

  // Save XML file
  std::string xmlDocBaseName = SimParms::getInstance().getWebDir();
  std::string xmlDocName = SimParms::getInstance().getLayoutName()  + ".xml";
  std::string xmlDefinitionsDocName = SimParms::getInstance().getLayoutName() + "_Definitions.xml";
  // Access all geometry parameters through visitor pattern
  //GeomExtractVisitor geomExtractVisitor(m_xmlDoc);

  // Access material database
  auto materialTab = MaterialTab::getInstance();

  // Readouts
  auto xmlDetReadouts = m_xmlDoc->NewElement("readouts");
  m_xmlNodeRoot->InsertEndChild(xmlDetReadouts);

  // Brl read-out
  auto xmlDetBrlReadout = m_xmlDoc->NewElement("readout");
  xmlDetBrlReadout->SetAttribute("name", c_defaultBrlReadout);
  xmlDetReadouts->InsertEndChild(xmlDetBrlReadout);

  auto xmlDetBrlSeg = m_xmlDoc->NewElement("segmentation");
  xmlDetBrlSeg->SetAttribute("type","CartesianGridXZ");
  xmlDetBrlSeg->SetAttribute("grid_size_x", printWithUnit(c_readoutGridX, c_precision, "mm").c_str());
  xmlDetBrlSeg->SetAttribute("grid_size_z", printWithUnit(c_readoutGridZ, c_precision, "mm").c_str());
  xmlDetBrlReadout->InsertEndChild(xmlDetBrlSeg);

  auto xmlDetBrlId = m_xmlDoc->NewElement("id");
  xmlDetBrlId->SetText("system:5,layer:5,module:16,component:4,x:-15,z:-15");
  xmlDetBrlReadout->InsertEndChild(xmlDetBrlId);

  // Ecap read-out
  auto xmlDetEcapReadout = m_xmlDoc->NewElement("readout");
  xmlDetEcapReadout->SetAttribute("name", c_defaultEcapReadout);
  xmlDetReadouts->InsertEndChild(xmlDetEcapReadout);

  auto xmlDetEcapSeg = m_xmlDoc->NewElement("segmentation");
  xmlDetEcapSeg->SetAttribute("type","CartesianGridXZ");
  xmlDetEcapSeg->SetAttribute("grid_size_x", printWithUnit(c_readoutGridX, c_precision, "mm").c_str());
  xmlDetEcapSeg->SetAttribute("grid_size_z", printWithUnit(c_readoutGridZ, c_precision, "mm").c_str());
  xmlDetEcapReadout->InsertEndChild(xmlDetEcapSeg);

  auto xmlDetEcapId = m_xmlDoc->NewElement("id");
  xmlDetEcapId->SetText("system:5,disc:5,module:16,component:4,x:-15,z:-15");
  xmlDetEcapReadout->InsertEndChild(xmlDetEcapId);


  // Initialize detector id
  short detId = c_defaultTrackerId;

  auto xmlIncludeDefinitions = m_xmlDoc->NewElement("include");
  xmlIncludeDefinitions->SetAttribute("ref", xmlDefinitionsDocName.c_str());
  m_xmlNodeRoot->InsertEndChild(xmlIncludeDefinitions);
  
  auto xmlSubdetectorAssemblies = m_xmlDoc->NewElement("detectors");
  m_xmlNodeRoot->InsertEndChild(xmlSubdetectorAssemblies);

  auto xmlSubdetectorBeampipe = m_xmlDoc->NewElement("detector");
  xmlSubdetectorBeampipe->SetAttribute("id", 0);
  xmlSubdetectorBeampipe->SetAttribute("name", "beampipe");
  xmlSubdetectorBeampipe->SetAttribute("type", "BeamTube");
  auto xmlBeamPipeDimensions = m_xmlDoc->NewElement("dimensions");
  xmlBeamPipeDimensions->SetAttribute("rmin", "CentralBeamTube_rmin");
  xmlBeamPipeDimensions->SetAttribute("rmax", "CentralBeamTube_rmax");
  xmlBeamPipeDimensions->SetAttribute("z", "CentralBeamTube_dz");
  xmlBeamPipeDimensions->SetAttribute("material", "Beryllium");
  xmlBeamPipeDimensions->SetAttribute("vis", "violet");
  xmlSubdetectorBeampipe->InsertEndChild(xmlBeamPipeDimensions);
  xmlSubdetectorAssemblies->InsertEndChild(xmlSubdetectorBeampipe);

  auto xmlSubdetectorInnermostBarrel = m_xmlDoc->NewElement("detector");
  xmlSubdetectorInnermostBarrel->SetAttribute("id", 100); // this id is not used in readout, just in dd4hep plugin
  xmlSubdetectorInnermostBarrel->SetAttribute("name", "FCChhInner0");
  xmlSubdetectorInnermostBarrel->SetAttribute("type", "DD4hep_SubdetectorAssembly");
  xmlSubdetectorInnermostBarrel->SetAttribute("vis", "BlueVisTrans");
  auto xmlInnermostComposite = m_xmlDoc->NewElement("composite");
  xmlInnermostComposite->SetAttribute("name", "InnerBRL_0");
  xmlSubdetectorInnermostBarrel->InsertEndChild(xmlInnermostComposite);
  xmlSubdetectorAssemblies->InsertEndChild(xmlSubdetectorInnermostBarrel);
  
  auto xmlSubdetectorInnerBarrel = m_xmlDoc->NewElement("detector");
  xmlSubdetectorInnerBarrel->SetAttribute("id", 101); // this id is not used in readout, just in dd4hep plugin
  xmlSubdetectorInnerBarrel->SetAttribute("name", "FCChhInner");
  xmlSubdetectorInnerBarrel->SetAttribute("type", "DD4hep_SubdetectorAssembly");
  xmlSubdetectorInnerBarrel->SetAttribute("vis", "BlueVisTrans");
  auto xmlInnerComposite = m_xmlDoc->NewElement("composite");
  xmlInnerComposite->SetAttribute("name", "InnerBRL_1");
  xmlSubdetectorInnerBarrel->InsertEndChild(xmlInnerComposite);
  auto xmlInnerCompositeEcapPos = m_xmlDoc->NewElement("composite");
  xmlInnerCompositeEcapPos->SetAttribute("name", "InnerECAPpos");
  xmlSubdetectorInnerBarrel->InsertEndChild(xmlInnerCompositeEcapPos);
  auto xmlInnerCompositeEcapNeg = m_xmlDoc->NewElement("composite");
  xmlInnerCompositeEcapNeg->SetAttribute("name", "InnerECAPneg");
  xmlSubdetectorInnerBarrel->InsertEndChild(xmlInnerCompositeEcapNeg);
  xmlSubdetectorAssemblies->InsertEndChild(xmlSubdetectorInnerBarrel);

  auto xmlSubdetectorOuterBarrel = m_xmlDoc->NewElement("detector");
  xmlSubdetectorOuterBarrel->SetAttribute("id", 102); // this id is not used in readout, just in dd4hep plugin
  xmlSubdetectorOuterBarrel->SetAttribute("name", "FCChhOuter");
  xmlSubdetectorOuterBarrel->SetAttribute("type", "DD4hep_SubdetectorAssembly");
  xmlSubdetectorOuterBarrel->SetAttribute("vis", "BlueVisTrans");
  auto xmlOuterComposite = m_xmlDoc->NewElement("composite");
  xmlOuterComposite->SetAttribute("name", "OuterBRL");
  xmlSubdetectorOuterBarrel->InsertEndChild(xmlOuterComposite);
  auto xmlOuterCompositeEcapPos = m_xmlDoc->NewElement("composite");
  xmlOuterCompositeEcapPos->SetAttribute("name", "OuterECAPpos");
  xmlSubdetectorOuterBarrel->InsertEndChild(xmlOuterCompositeEcapPos);
  auto xmlOuterCompositeEcapNeg = m_xmlDoc->NewElement("composite");
  xmlOuterCompositeEcapNeg->SetAttribute("name", "OuterECAPneg");
  xmlSubdetectorOuterBarrel->InsertEndChild(xmlOuterCompositeEcapNeg);
  xmlSubdetectorAssemblies->InsertEndChild(xmlSubdetectorOuterBarrel);


  // Go through the whole geometry hierarchy & build corresponding XML node
  auto xmlDetectors = m_xmlDefinitionsDoc->NewElement("detectors");
  m_xmlDefinitionsDoc->InsertFirstChild(xmlDetectors);

  auto xmlDetComment = m_xmlDoc->NewElement("comment");
  std::string text   = "The tracker geometry as described here follows the concept of CMS tracker design. Each sub-tracker (inner, outer, etc. tracker) is either of barrel type or end-cape type. ";
  text              += "In case of barrel-type, individual modules are arranged in layers, layers consist of individual rods (ladders) and each rod contains several silicon modules. ";
  text              += "Within a layer, even/odd rods are offset by +/-bigDelta parameter with respect to the layer average radius (rods are arranged in so-called staggered structure). ";
  text              += "Further on, within each rod the modules are positioned in Z in a way, that they hermetically cover the primary vertex (taking also into account the Gaussian shape of p-p beam-spot). ";
  text              += "Similar concept is hence used in positioning in Z, even/odd modules are offset by +-smallDelta parameter with respect to the rod average radius. ";
  text              += "In case of end-cap type, individual modules are arranged in discs, discs consist of individual rings and each ring contains several silicon modules. ";
  text              += "The modules hermetically cover the primary vertex up-to required eta. ";
  text              += "In the same manner as for the rods, the rings are offset by +-/bigDelta parameter and as for the modules within a rod, the modules are offset by +/-smallDelta parameter. ";
  text              += "Each sub-tracker is assumed to have a phi-symmetry, therefore the data are saved in the XML file in a way to avoid their duplication. ";
  text              += "Each layer is described only by first two rods: an even and odd rod. ";
  text              += "Similarly, each ring is described by first two modules: an even and odd module. ";
  text              += "The symmetry considerations are also applied in a description of individual end-cap sub-trackers, where +/-Z symmetry is naturally assumed. ";
  text              += "Therefore, only information about discs placed at positive Z are written out. ";
  text              += "The end-cap tracker is designed in a way that the geometry of individual discs is replicated accross the whole end-cap sub-detector. ";
  text              += "The XML file describes only one disc for each end-cap detector. ";
  text              += "Additionally, the modules material budget as described by moduleProperties tag (as a sequence of material components) is assumed to be for simplicity uniformaly distributed across the module surface. ";
  text              += "Individual material properties are specified as defined within tkLayout, i.e. by their density, radiation and interaction lengths. ";
  text              += "Finally, positioning of individual modules have been optimized in the tkLayout software. Hence, the relevant geometry parameters are directly the X, Y, Z coordinates of modules ";
  text              += "centre position plus modules rotation: in phi (phiTilt) and theta (thetaTilt=0*deg for barrel modules, thetaTilt=90*deg for end-cap modules. ";
  text              += "Other parameters, based on which tkLayout algorithms have found the optimal positions of tracker modules are given just for completeness. ";

  xmlDetComment->SetText(text.c_str());
  xmlDetectors->InsertEndChild(xmlDetComment);

  auto xmlDetComment2 = m_xmlDoc->NewElement("comment");
  text                = "The tracker numbering scheme applied in the XML is designed as follows: all sub-trackers, barrels/discs, layers/rings etc. are numbered by increasing ID from inside, i.e. IP, out. ";
  text               += "Individual subtrackers' ID start with 10, subtracker components' ID with 1. In addition, no global numbering scheme accross the overall tracker is applied. Instead, ";
  text               += "the hierarchical scheme is used, e.g. individual rings numbering starts with ID=1 for each disc, modules numbering with ID=1 for each ring etc. ";
  text               += "As for the sub-components arranged in R-Phi, ID increases along the increasing phi angle (rotation along Z-axis in right-handed coordinate scheme), starting with ID=1 ";
  text               += "for a component positioned at phi0. As for the sub-components arranged in Z (e.g. modules in a rod), ID increases from -Z to +Z.";

  xmlDetComment2->SetText(text.c_str());
  xmlDetectors->InsertEndChild(xmlDetComment2);

  for (auto iTrk : m_trackers) {
    for (const auto& iBrl : iTrk->barrels()) {

      // Create detector element assigned to given barrel
      auto xmlBrlDet = m_xmlDefinitionsDoc->NewElement("detector");
      std::string brlName = std::string(iTrk->myid()+iBrl.myid());
      xmlBrlDet->SetAttribute("name", brlName.c_str());
      xmlBrlDet->SetAttribute("id", detNames2IDs[brlName].c_str());
      xmlBrlDet->SetAttribute("type", c_defaultBrlGeoCreator);
      xmlBrlDet->SetAttribute("readout", c_defaultBrlReadout);
      xmlDetectors->InsertEndChild(xmlBrlDet);

      // Add detailed info about given barrel
      auto xmlBrlDim = m_xmlDefinitionsDoc->NewElement("dimensions");

      xmlBrlDim->SetAttribute("rmin", printWithUnit(iBrl.minRAllMat(), c_precision, "mm").c_str());
      xmlBrlDim->SetAttribute("rmax", printWithUnit(iBrl.maxRAllMat(), c_precision, "mm").c_str());
      xmlBrlDim->SetAttribute("zmin", printWithUnit(iBrl.minZ(), c_precision, "mm").c_str());
      xmlBrlDim->SetAttribute("zmax", printWithUnit(iBrl.maxZ(), c_precision, "mm").c_str());
      xmlBrlDet->InsertEndChild(xmlBrlDim);

      // Add sensitive type
      auto xmlBrlSens = m_xmlDefinitionsDoc->NewElement("sensitive");
      xmlBrlSens->SetAttribute("type", c_defaultSensDet);
      xmlBrlDet->InsertEndChild(xmlBrlSens);

      // Add detailed info about layers
      auto xmlBrlLayers = m_xmlDefinitionsDoc->NewElement("layers");
      xmlBrlLayers->SetAttribute("repeat" , iBrl.numLayers());
      for (const auto& iLayer : iBrl.layers()) {

        if (iLayer.myid()==1) {
          xmlBrlDet->InsertEndChild(xmlBrlLayers);
        }

        // Add detailed info about individual layer
        auto xmlBrlIthLayer = m_xmlDefinitionsDoc->NewElement("layer");
        xmlBrlIthLayer->SetAttribute("id",            iLayer.myid());
        xmlBrlIthLayer->SetAttribute("radius",        printWithUnit(iLayer.avgBuildRadius(), c_precision, "mm").c_str());
        xmlBrlIthLayer->SetAttribute("rmin",          printWithUnit(iLayer.minRAllMat()    , c_precision, "mm").c_str());
        xmlBrlIthLayer->SetAttribute("rmax",          printWithUnit(iLayer.maxRAllMat()    , c_precision, "mm").c_str());
        xmlBrlIthLayer->SetAttribute("bigDelta",      printWithUnit(iLayer.bigDelta(),       c_precision, "mm").c_str());
        xmlBrlIthLayer->SetAttribute("phi0",          printWithUnit(iLayer.layerRotation(),2*c_precision, "rad").c_str());
        xmlBrlLayers->InsertEndChild(xmlBrlIthLayer);

        // Rods general info
        auto xmlBrlRods = m_xmlDefinitionsDoc->NewElement("rods");

        // Add detailed info about straight rods
        // TODO: Add tilted rods
        for (const auto& iRod : iLayer.flatRods()) {

          // Odd rods
          if (iRod.myid()==1) {

            // Fill general info
            xmlBrlRods->SetAttribute("repeat",        iLayer.numRods());
            xmlBrlRods->SetAttribute("smallDelta",    printWithUnit(iLayer.smallDelta(), c_precision, "mm").c_str());
            xmlBrlRods->SetAttribute("rPhiOverlap",   printWithUnit(iLayer.phiOverlap(), c_precision, "mm").c_str());
            xmlBrlRods->SetAttribute("nRPhiSegments", iLayer.phiSegments());
            try {
              xmlBrlRods->SetAttribute("zOverlap", printWithUnit((dynamic_cast<const RodPairStraight&>(iRod)).zOverlap(), c_precision, "mm").c_str());
            }
            catch(std::bad_cast& e) {}
            xmlBrlRods->SetAttribute("nModules", iRod.numModules());

            xmlBrlIthLayer->InsertEndChild(xmlBrlRods);

            // Odd rod
            auto xmlBrlRodOdd = m_xmlDefinitionsDoc->NewElement("rodOdd");
            xmlBrlRodOdd->SetAttribute("id", iRod.myid());
            xmlBrlRods->InsertEndChild(xmlBrlRodOdd);

            // Add all modules info (negative Z modules first, then positive)
            auto xmlBrlModules = m_xmlDefinitionsDoc->NewElement("modules");
            xmlBrlRodOdd->InsertEndChild(xmlBrlModules);

            auto xmlBrlModProperties = m_xmlDefinitionsDoc->NewElement("moduleProperties");
            xmlBrlRodOdd->InsertEndChild(xmlBrlModProperties);

            auto xmlBrlSensorProperties = m_xmlDefinitionsDoc->NewElement("sensorProperties");
            xmlBrlRodOdd->InsertEndChild(xmlBrlSensorProperties);

            int idMod = 1;

            for (const auto& iMod : boost::adaptors::reverse(iRod.modules().second)) {

              // Module properties
              if (iMod.myid()==1) {
                xmlBrlModProperties->SetAttribute("modLength"      , printWithUnit(iMod.physicalLength(), c_precision, "mm").c_str());
                xmlBrlModProperties->SetAttribute("modWidth"       , printWithUnit(iMod.meanWidth(),      c_precision, "mm").c_str());

                xmlBrlSensorProperties->SetAttribute("sensorLength"   , printWithUnit(iMod.length(),   c_precision, "mm").c_str());
                xmlBrlSensorProperties->SetAttribute("sensorWidth"    , printWithUnit(iMod.meanWidth(),c_precision, "mm").c_str());
                xmlBrlSensorProperties->SetAttribute("sensorThickness", printWithUnit(iMod.thickness(),c_precision, "mm").c_str());
                xmlBrlSensorProperties->SetAttribute("resRPhi"        , printWithUnit(iMod.resolutionLocalX(), 1, "um").c_str());
                xmlBrlSensorProperties->SetAttribute("resZ"           , printWithUnit(iMod.resolutionLocalY(), 1, "um").c_str());

                auto xmlBrlModComponents = m_xmlDefinitionsDoc->NewElement("components");
                xmlBrlModProperties->InsertEndChild(xmlBrlModComponents);

                double modThick = 0;

                for (auto& elem : iMod.materialObject().getLocalElements()) {

                  auto xmlBrlModComponent = m_xmlDefinitionsDoc->NewElement("component");
                  xmlBrlModComponent->SetAttribute("name",      elem->componentName().c_str());
                  xmlBrlModComponent->SetAttribute("thickness", printWithUnit(elem->quantity(), c_precision, elem->unit()).c_str());
                  xmlBrlModComponent->SetAttribute("material",  elem->elementName().c_str());
                  // TODO: Density, rad & length have fixed units in tkLayout -> don't apply unit coef.
                  xmlBrlModComponent->SetAttribute("density",   printWithUnit(materialTab.density(elem->elementName()),          c_precision, "g/cm3").c_str());
                  xmlBrlModComponent->SetAttribute("radLength", printWithUnit(materialTab.radiationLength(elem->elementName()),  c_precision, "g/cm2").c_str());
                  xmlBrlModComponent->SetAttribute("intLength", printWithUnit(materialTab.interactionLength(elem->elementName()),c_precision, "g/cm2").c_str());
                  if (elem->quantity()==iMod.thickness()) xmlBrlModComponent->SetAttribute("sensitive", "true");
                  else                                    xmlBrlModComponent->SetAttribute("sensitive", "false");
                  xmlBrlModComponents->InsertEndChild(xmlBrlModComponent);

                  modThick += elem->quantity();
                }
                xmlBrlModProperties->SetAttribute("modThickness", printWithUnit(modThick, c_precision, "mm").c_str());
              }

              // negative Z position & rotation
              auto xmlBrlMod = m_xmlDefinitionsDoc->NewElement("module");
              xmlBrlMod->SetAttribute("id",       idMod++);
              xmlBrlMod->SetAttribute("X",        printWithUnit(iMod.center().X(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("Y",        printWithUnit(iMod.center().Y(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("Z",        printWithUnit(iMod.center().Z(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("phiTilt",  printWithUnit(iMod.skewAngle(), 2*c_precision, "rad").c_str());
              xmlBrlMod->SetAttribute("thetaTilt",printWithUnit(iMod.tiltAngle(), 2*c_precision, "rad").c_str());
              xmlBrlModules->InsertEndChild(xmlBrlMod);
            }
            for (const auto& iMod : iRod.modules().first) {

              // positive Z position & rotation
              auto xmlBrlMod = m_xmlDefinitionsDoc->NewElement("module");
              xmlBrlMod->SetAttribute("id",       idMod++);
              xmlBrlMod->SetAttribute("X",        printWithUnit(iMod.center().X(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("Y",        printWithUnit(iMod.center().Y(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("Z",        printWithUnit(iMod.center().Z(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("phiTilt",  printWithUnit(iMod.skewAngle(), 2*c_precision, "rad").c_str());
              xmlBrlMod->SetAttribute("thetaTilt",printWithUnit(iMod.tiltAngle(), 2*c_precision, "rad").c_str());
              xmlBrlModules->InsertEndChild(xmlBrlMod);
            }
          }
          // Even rods
          else if (iRod.myid()==2) {

            auto xmlBrlRodEven = m_xmlDefinitionsDoc->NewElement("rodEven");
            xmlBrlRodEven->SetAttribute("id", iRod.myid());
            xmlBrlRods->InsertEndChild(xmlBrlRodEven);

            // Add all modules info (negative Z modules first, then positive)
            auto xmlBrlModules = m_xmlDefinitionsDoc->NewElement("modules");
            xmlBrlRodEven->InsertEndChild(xmlBrlModules);

            auto xmlBrlModProperties = m_xmlDefinitionsDoc->NewElement("moduleProperties");
            xmlBrlRodEven->InsertEndChild(xmlBrlModProperties);

            auto xmlBrlSensorProperties = m_xmlDefinitionsDoc->NewElement("sensorProperties");
            xmlBrlRodEven->InsertEndChild(xmlBrlSensorProperties);

            int idMod = 1;

            for (const auto& iMod : boost::adaptors::reverse(iRod.modules().second)) {

              // Module properties
              if (iMod.myid()==1) {
                xmlBrlModProperties->SetAttribute("modLength"      , printWithUnit(iMod.physicalLength(), c_precision, "mm").c_str());
                xmlBrlModProperties->SetAttribute("modWidth"       , printWithUnit(iMod.meanWidth(),      c_precision, "mm").c_str());

                xmlBrlSensorProperties->SetAttribute("sensorLength"   , printWithUnit(iMod.length(),   c_precision, "mm").c_str());
                xmlBrlSensorProperties->SetAttribute("sensorWidth"    , printWithUnit(iMod.meanWidth(),c_precision, "mm").c_str());
                xmlBrlSensorProperties->SetAttribute("sensorThickness", printWithUnit(iMod.thickness(),c_precision, "mm").c_str());
                xmlBrlSensorProperties->SetAttribute("resRPhi"        , printWithUnit(iMod.resolutionLocalX(), 1, "um").c_str());
                xmlBrlSensorProperties->SetAttribute("resZ"           , printWithUnit(iMod.resolutionLocalY(), 1, "um").c_str());

                auto xmlBrlModComponents = m_xmlDefinitionsDoc->NewElement("components");
                xmlBrlModProperties->InsertEndChild(xmlBrlModComponents);

                double modThick = 0;

                for (auto& elem : iMod.materialObject().getLocalElements()) {

                  auto xmlBrlModComponent = m_xmlDefinitionsDoc->NewElement("component");
                  xmlBrlModComponent->SetAttribute("name",      elem->componentName().c_str());
                  xmlBrlModComponent->SetAttribute("thickness", printWithUnit(elem->quantity(), c_precision, elem->unit()).c_str());
                  xmlBrlModComponent->SetAttribute("material",  elem->elementName().c_str());
                  // TODO: Density, rad & length have fixed units in tkLayout -> don't apply unit coef.
                  xmlBrlModComponent->SetAttribute("density",   printWithUnit(materialTab.density(elem->elementName()),          c_precision, "g/cm3").c_str());
                  xmlBrlModComponent->SetAttribute("radLength", printWithUnit(materialTab.radiationLength(elem->elementName()),  c_precision, "g/cm2").c_str());
                  xmlBrlModComponent->SetAttribute("intLength", printWithUnit(materialTab.interactionLength(elem->elementName()),c_precision, "g/cm2").c_str());
                  if (elem->quantity()==iMod.thickness()) xmlBrlModComponent->SetAttribute("sensitive", "true");
                  else                                    xmlBrlModComponent->SetAttribute("sensitive", "false");
                  xmlBrlModComponents->InsertEndChild(xmlBrlModComponent);

                  modThick += elem->quantity();
                }
                xmlBrlModProperties->SetAttribute("modThickness", printWithUnit(modThick, c_precision, "mm").c_str());
              }

              // negative Z position & rotation
              auto xmlBrlMod = m_xmlDefinitionsDoc->NewElement("module");
              xmlBrlMod->SetAttribute("id",       idMod++);
              xmlBrlMod->SetAttribute("X",        printWithUnit(iMod.center().X(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("Y",        printWithUnit(iMod.center().Y(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("Z",        printWithUnit(iMod.center().Z(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("phiTilt",  printWithUnit(iMod.skewAngle(), 2*c_precision, "rad").c_str());
              xmlBrlMod->SetAttribute("thetaTilt",printWithUnit(iMod.tiltAngle(), 2*c_precision, "rad").c_str());
              xmlBrlModules->InsertEndChild(xmlBrlMod);
            }
            for (const auto& iMod : iRod.modules().first) {

              // positive Z position & rotation
              auto xmlBrlMod = m_xmlDefinitionsDoc->NewElement("module");
              xmlBrlMod->SetAttribute("id",       idMod++);
              xmlBrlMod->SetAttribute("X",        printWithUnit(iMod.center().X(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("Y",        printWithUnit(iMod.center().Y(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("Z",        printWithUnit(iMod.center().Z(), c_precision, "mm").c_str());
              xmlBrlMod->SetAttribute("phiTilt",  printWithUnit(iMod.skewAngle(), 2*c_precision, "rad").c_str());
              xmlBrlMod->SetAttribute("thetaTilt",printWithUnit(iMod.tiltAngle(), 2*c_precision, "rad").c_str());
              xmlBrlModules->InsertEndChild(xmlBrlMod);
            }

          }
          else {
            break;
          } // Only first odd & even rod read-in
        } // Rods
      } // Layers
    } // Barrels
    std::array<std::string, 2> posnegStringArray {{"pos", "neg"}};
    for (auto posnegString: posnegStringArray) {
    for (const auto& iEcap : iTrk->endcaps()) {

      // Create detector element assigned to given endcap
      auto xmlEcapDet = m_xmlDefinitionsDoc->NewElement("detector");
      std::string ecapName = std::string(iTrk->myid()+iEcap.myid()+posnegString);
      xmlEcapDet->SetAttribute("name", ecapName.c_str());
      xmlEcapDet->SetAttribute("id", detNames2IDs[ecapName].c_str());
      xmlEcapDet->SetAttribute("type", c_defaultEcapGeoCreator);
      xmlEcapDet->SetAttribute("readout", c_defaultEcapReadout);
      xmlDetectors->InsertEndChild(xmlEcapDet);

      // Add detailed info about given endcap
      auto xmlEcapDim = m_xmlDefinitionsDoc->NewElement("dimensions");

      xmlEcapDim->SetAttribute("rmin", printWithUnit(iEcap.minR(), c_precision, "mm").c_str());
      xmlEcapDim->SetAttribute("rmax", printWithUnit(iEcap.maxR(), c_precision, "mm").c_str());
      xmlEcapDim->SetAttribute("zmin", printWithUnit(iEcap.minZAllMat(), c_precision, "mm").c_str());
      xmlEcapDim->SetAttribute("zmax", printWithUnit(iEcap.maxZAllMat(), c_precision, "mm").c_str());
      if ("pos" == posnegString) {
        xmlEcapDim->SetAttribute("reflect", "false");
      } else {
        xmlEcapDim->SetAttribute("reflect", "true");
      }
      xmlEcapDet->InsertEndChild(xmlEcapDim);

      // Add sensitive type
      auto xmlEcapSens = m_xmlDefinitionsDoc->NewElement("sensitive");
      xmlEcapSens->SetAttribute("type", c_defaultSensDet);
      xmlEcapDet->InsertEndChild(xmlEcapSens);

      // Add detailed info about discs
      auto xmlEcapDiscs = m_xmlDefinitionsDoc->NewElement("discs");
      xmlEcapDiscs->SetAttribute("repeat" , iEcap.numDisks());
      xmlEcapDet->InsertEndChild(xmlEcapDiscs);

      for (const auto& iDisc : iEcap.disks()) {

        //// Print out only pos. Z
        if (iDisc.averageZ()>0) {

          // Add detailed info about individual disc
          auto xmlEcapIthDisc = m_xmlDefinitionsDoc->NewElement("discZPls");
          xmlEcapIthDisc->SetAttribute("id",       iDisc.myid());
          xmlEcapIthDisc->SetAttribute("z",        printWithUnit(iDisc.averageZ(),   c_precision, "mm").c_str());
          xmlEcapIthDisc->SetAttribute("zmin",     printWithUnit(iDisc.minZAllMat(), c_precision, "mm").c_str());
          xmlEcapIthDisc->SetAttribute("zmax",     printWithUnit(iDisc.maxZAllMat(), c_precision, "mm").c_str());
          xmlEcapIthDisc->SetAttribute("bigDelta", printWithUnit(iDisc.bigDelta(),   c_precision, "mm").c_str());
          xmlEcapDiscs->InsertEndChild(xmlEcapIthDisc);

          // Rods general info
          auto xmlEcapRings = m_xmlDefinitionsDoc->NewElement("rings");

          // Add detailed info about rods
          for (const auto& iRing : iDisc.rings()) {

            if (iRing.myid()==1) {

              // Fill general info
              xmlEcapRings->SetAttribute("repeat",        iDisc.numRings());
              xmlEcapRings->SetAttribute("smallDelta",    printWithUnit(iRing.smallDelta(), c_precision, "mm").c_str());
              xmlEcapRings->SetAttribute("rPhiOverlap",   printWithUnit(iRing.phiOverlap(), c_precision, "mm").c_str());
              xmlEcapRings->SetAttribute("nRPhiSegments", iRing.phiSegments());
              xmlEcapRings->SetAttribute("rOverlap",      printWithUnit(iDisc.rOverlap(), c_precision, "mm" ).c_str());

              xmlEcapIthDisc->InsertEndChild(xmlEcapRings);
            }

            // Discs are symmetrically shifted - fill only the first disc
            if (iDisc.myid()==1) {
              auto xmlEcapRing = m_xmlDefinitionsDoc->NewElement("ring");
              xmlEcapRing->SetAttribute("id",       iRing.myid());
              xmlEcapRing->SetAttribute("phi0",     printWithUnit(iRing.zRotation(),2*c_precision, "rad").c_str());
              xmlEcapRing->SetAttribute("nModules", iRing.numModules());
              xmlEcapRings->InsertEndChild(xmlEcapRing);

              // Add all modules info (negative Z modules first, then positive)
              auto xmlEcapModules = m_xmlDefinitionsDoc->NewElement("modules");
              xmlEcapRing->InsertEndChild(xmlEcapModules);

              auto xmlEcapModProperties = m_xmlDefinitionsDoc->NewElement("moduleProperties");
              xmlEcapRing->InsertEndChild(xmlEcapModProperties);

              auto xmlEcapSensorProperties = m_xmlDefinitionsDoc->NewElement("sensorProperties");
              xmlEcapRing->InsertEndChild(xmlEcapSensorProperties);

              for (const auto& iMod : iRing.modules()) {

                // Odd module -> use symmetry to build the ring
                if (iMod.myid()==1) {

                  auto xmlOddMod = m_xmlDefinitionsDoc->NewElement("moduleOdd");
                  xmlOddMod->SetAttribute("id",       iMod.myid());
                  xmlOddMod->SetAttribute("X",        printWithUnit(iMod.center().X(), c_precision, "mm").c_str());
                  xmlOddMod->SetAttribute("Y",        printWithUnit(iMod.center().Y(), c_precision, "mm").c_str());
                  xmlOddMod->SetAttribute("Z",        printWithUnit(iMod.center().Z(), c_precision, "mm").c_str());
                  xmlOddMod->SetAttribute("phiTilt",  printWithUnit(iMod.skewAngle(), 2*c_precision, "rad").c_str());
                  xmlOddMod->SetAttribute("thetaTilt",printWithUnit(iMod.tiltAngle(), 2*c_precision, "rad").c_str());
                  xmlEcapModules->InsertEndChild(xmlOddMod);

                  double modLength = iMod.physicalLength();
                  if (modLength==0.0) modLength = iMod.length(); // Parameter not specified, to build module wafer the radius was used instead -> hence use direcly length parameter
                  xmlEcapModProperties->SetAttribute("modLength"      , printWithUnit(modLength,       c_precision, "mm").c_str());
                  xmlEcapModProperties->SetAttribute("modWidthMin"    , printWithUnit(iMod.minWidth(), c_precision, "mm").c_str());
                  xmlEcapModProperties->SetAttribute("modWidthMax"    , printWithUnit(iMod.maxWidth(), c_precision, "mm").c_str());

                  xmlEcapSensorProperties->SetAttribute("sensorLength"   , printWithUnit(iMod.length(),   c_precision, "mm").c_str());
                  xmlEcapSensorProperties->SetAttribute("sensorWidthMin" , printWithUnit(iMod.minWidth(), c_precision, "mm").c_str());
                  xmlEcapSensorProperties->SetAttribute("sensorWidthMax" , printWithUnit(iMod.maxWidth(), c_precision, "mm").c_str());
                  xmlEcapSensorProperties->SetAttribute("sensorThickness", printWithUnit(iMod.thickness(),c_precision, "mm").c_str());
                  xmlEcapSensorProperties->SetAttribute("resRPhi"        , printWithUnit(iMod.resolutionLocalX(), 1, "um").c_str());
                  xmlEcapSensorProperties->SetAttribute("resZ"           , printWithUnit(iMod.resolutionLocalY(), 1, "um").c_str());

                  auto xmlEcapModComponents = m_xmlDefinitionsDoc->NewElement("components");
                  xmlEcapModProperties->InsertEndChild(xmlEcapModComponents);

                  double modThick = 0;

                  for (auto& elem : iMod.materialObject().getLocalElements()) {

                    auto xmlEcapModComponent = m_xmlDefinitionsDoc->NewElement("component");
                    xmlEcapModComponent->SetAttribute("name",      elem->componentName().c_str());
                    xmlEcapModComponent->SetAttribute("thickness", printWithUnit(elem->quantity(), c_precision, elem->unit()).c_str());
                    xmlEcapModComponent->SetAttribute("material",  elem->elementName().c_str());
                    // TODO: Density, rad & length have fixed units in tkLayout -> don't apply unit coef.
                    xmlEcapModComponent->SetAttribute("density",   printWithUnit(materialTab.density(elem->elementName()),          c_precision, "g/cm3").c_str());
                    xmlEcapModComponent->SetAttribute("radLength", printWithUnit(materialTab.radiationLength(elem->elementName()),  c_precision, "g/cm2").c_str());
                    xmlEcapModComponent->SetAttribute("intLength", printWithUnit(materialTab.interactionLength(elem->elementName()),c_precision, "g/cm2").c_str());
                    if (elem->quantity()==iMod.thickness()) xmlEcapModComponent->SetAttribute("sensitive", "true");
                    else                                    xmlEcapModComponent->SetAttribute("sensitive", "false");
                    xmlEcapModComponents->InsertEndChild(xmlEcapModComponent);

                    modThick += elem->quantity();
                  }
                  xmlEcapModProperties->SetAttribute("modThickness", printWithUnit(modThick, c_precision, "mm").c_str());
                }
                // Even module -> use symmetry to build the ring
                else if (iMod.myid()==2) {

                  auto xmlEvenMod = m_xmlDefinitionsDoc->NewElement("moduleEven");
                  xmlEvenMod->SetAttribute("id",       iMod.myid());
                  xmlEvenMod->SetAttribute("X",        printWithUnit(iMod.center().X(), c_precision, "mm").c_str());
                  xmlEvenMod->SetAttribute("Y",        printWithUnit(iMod.center().Y(), c_precision, "mm").c_str());
                  xmlEvenMod->SetAttribute("Z",        printWithUnit(iMod.center().Z(), c_precision, "mm").c_str());
                  xmlEvenMod->SetAttribute("phiTilt",  printWithUnit(iMod.skewAngle(), 2*c_precision, "rad").c_str());
                  xmlEvenMod->SetAttribute("thetaTilt",printWithUnit(iMod.tiltAngle(), 2*c_precision, "rad").c_str());
                  xmlEcapModules->InsertEndChild(xmlEvenMod);

                }
                else {
                  break;
                }
              } // Modules
            } // If pos. disk
          } // Rings
        }
      } // Discs
    } // Endcaps
  }

  } // Trackers

  auto eResult = m_xmlDoc->SaveFile((xmlDocBaseName + "/" + xmlDocName).c_str());
  auto eResultDefinitions = m_xmlDefinitionsDoc->SaveFile((xmlDocBaseName + "/" + xmlDefinitionsDocName).c_str());

  if (eResult==XML_SUCCESS && eResultDefinitions==XML_SUCCESS) {

    logINFO(std::string("Geometry extracted to the following XML file: "+xmlDocName));
    m_isAnalysisOK = true;
  }
  else {

    logERROR("ExtractorFCCSW::analyze(): A problem occured when saving an XML file!");
    m_isAnalysisOK = false;
  }

  return m_isAnalysisOK;
}

//
// Helper method printing out number with unit & given precision
//
std::string ExtractorFCCSW::printWithUnit(double value, int precision, std::string unit)
{
  float unitCoef = 1;

  if      (unit=="mm")   unitCoef = Units::mm;
  else if (unit=="cm")   unitCoef = Units::cm;
  else if (unit=="m")    unitCoef = Units::m;
  else if (unit=="um")   unitCoef = Units::um;
  else if (unit=="rad")  unitCoef = 1;
  else if (unit=="g/cm3")unitCoef = Units::g/Units::cm3;
  else if (unit=="g/cm2")unitCoef = Units::g/Units::cm2;
  else if (unit=="g/cm") unitCoef = Units::g/Units::cm;
  else if (unit=="g/mm3")unitCoef = Units::g/Units::mm3;
  else if (unit=="g/mm2")unitCoef = Units::g/Units::mm2;
  else if (unit=="g/mm") unitCoef = Units::g/Units::mm;
  else {

    logERROR(std::string("Extractor::printWithUnit - Unsupported unit: ")+unit);
  }

  // Avoid minus zero values for the required precision
  if (fabs(value/unitCoef*pow(10,precision))<0.1) value = 0;

  std::stringstream myNum;
  myNum.clear();
  myNum << std::dec << std::fixed << std::setprecision(precision) << value/unitCoef;

  std::string numWithUnit = myNum.str()+"*"+unit;

  return numWithUnit;


}

////
//// Helper class constructer
////
//GeomExtractVisitor::GeomExtractVisitor(tinyxml2::XMLDocument& xmlDoc) :
// m_xmlDoc(xmlDoc),
// m_trkID(0),
//{}
//
//void GeomExtractVisitor::visit(const Tracker& t)
//{
//
//
//}
////  void visit(const Barrel& b) override;
////  void visit(const Endcap& e) override;
////  void visit(const Layer& l) override;
////  void visit(const Disk& d) override;
////  void visit(const Ring& r) override;
////  void visit(const DetectorModule& m) override;
////  void visit(const EndcapModule& m) override;
////  void postVisit();

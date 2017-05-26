#include "AnalyzerVisitors/MaterialBillAnalyzer.h"
#include "InactiveElement.h"
#include "MaterialBudget.h"
#include "ModuleCap.h"

using namespace insur;

/*
 * Constructor
 */
MaterialBillAnalyzer::MaterialBillAnalyzer() {

  // Set label part of the csv container
  m_matBillCsv.addCsvElement("Label", "Material bill");
  m_matBillCsv.addCsvEOL("Label");
  m_matBillCsv.addCsvEOL("Label");
}

/*
 *  Test if text block exists
 */
bool MaterialBillAnalyzer::existCsvText(std::string id) {

  if (m_matBillCsv.existCsvText(id)) return true;
  else                               return false;
}

/*
 * Calculate material bill for given tracker and fill internal Csv table with results
 */
void MaterialBillAnalyzer::inspectTracker(insur::MaterialBudget& mb) {

  // Tracker id
  std::string id = mb.getTracker().myid();

  // Go through all barrel and end-cap modules
  inspectModules(mb.getBarrelModuleCaps());
  inspectModules(mb.getEndcapModuleCaps());

  // Set label for barrel/end-cap modules
  m_matBillCsv.addCsvElement(id, "Material in barrel/end-cap layers");
  m_matBillCsv.addCsvElement(id, "Layer name");
  m_matBillCsv.addCsvElement(id, "Material");
  m_matBillCsv.addCsvElement(id, "Weight [g]");
  m_matBillCsv.addCsvEOL(id);
  m_matBillCsv.addCsvElement(id, std::string("Tracker: "+id));
  m_matBillCsv.addCsvEOL(id);

  // Set csv formatted layer content
  for (const auto &it : m_layerMaterialMap) {
    for (const auto &layMats : it.second) {

      m_matBillCsv.addCsvElement(id, "");
      m_matBillCsv.addCsvElement(id, it.first);
      m_matBillCsv.addCsvElement(id, layMats.first);
      m_matBillCsv.addCsvElement(id, layMats.second);
      m_matBillCsv.addCsvEOL(id);
    }
  } // For layers

  // Set services & supports content
  InactiveSurfaces& is = mb.getInactiveSurfaces();
  std::vector<insur::InactiveElement>& barrelServices = is.getBarrelServices();
  std::vector<insur::InactiveElement>& endcapServices = is.getEndcapServices();
  std::vector<insur::InactiveElement>& supports       = is.getSupports();

  m_matBillCsv.addCsvElement(id, "Material in services & supports");
  m_matBillCsv.addCsvElement(id, "R_in [mm]");
  m_matBillCsv.addCsvElement(id, "R_out [mm]");
  m_matBillCsv.addCsvElement(id, "Z_in [mm]");
  m_matBillCsv.addCsvElement(id, "Z_out [mm]");
  m_matBillCsv.addCsvElement(id, "Material");
  m_matBillCsv.addCsvElement(id, "Weight [g]");
  m_matBillCsv.addCsvEOL(id);
  m_matBillCsv.addCsvElement(id, std::string("Tracker: "+id));
  m_matBillCsv.addCsvEOL(id);

  inspectInactiveElements(id, barrelServices);
  inspectInactiveElements(id, endcapServices);
  inspectInactiveElements(id, supports);
}

/*
 * Internal method getting material information about detector modules
 */
void MaterialBillAnalyzer::inspectModules(std::vector<std::vector<insur::ModuleCap> >& tracker) {

  // Loop over layers
  for (auto layerIt : tracker ) {
    // Loop over modules
    for (auto moduleIt : layerIt ) {

      auto myModuleCap = (moduleIt);
      auto myModule    = &(myModuleCap.getModule());

      // Fill layer material map
      MaterialVisitor visitor;
      myModule->accept(visitor);

      MaterialMap& layerMaterial = m_layerMaterialMap[visitor.getId()];
      const std::map<std::string, double>& localMasses   = myModuleCap.getLocalMasses();
      const std::map<std::string, double>& exitingMasses = myModuleCap.getExitingMasses();
      for (const auto &it : localMasses)   layerMaterial[it.first]+=it.second;
      for (const auto &it : exitingMasses) layerMaterial[it.first]+=it.second;
    }
  }
}

/*
 * Internal method getting material information about services
 */
void MaterialBillAnalyzer::inspectInactiveElements(std::string id, const std::vector<insur::InactiveElement>& inactiveElements) {

  // Loop over services
  for (const auto& it : inactiveElements) {

    const std::map<std::string, double>& localMasses   = it.getLocalMasses();
    const std::map<std::string, double>& exitingMasses = it.getExitingMasses();

    std::map<std::string, double> totalMasses;
    for (auto massIt : localMasses)   totalMasses[massIt.first]+=massIt.second;
    for (auto massIt : exitingMasses) totalMasses[massIt.first]+=massIt.second;
    for (auto massIt : totalMasses) {

      // Set csv formatted services content
      m_matBillCsv.addCsvElement(id, "");
      m_matBillCsv.addCsvElement(id, any2str(it.getInnerRadius()) );
      m_matBillCsv.addCsvElement(id, any2str(it.getInnerRadius()+it.getRWidth()) );
      m_matBillCsv.addCsvElement(id, any2str(it.getZOffset()) );
      m_matBillCsv.addCsvElement(id, any2str(it.getZOffset()+it.getZLength()) );
      m_matBillCsv.addCsvElement(id, massIt.first );
      m_matBillCsv.addCsvElement(id, massIt.second );
      m_matBillCsv.addCsvEOL(id);
    }
  } // For services
}






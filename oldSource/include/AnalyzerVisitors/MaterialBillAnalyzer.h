#ifndef MATERIAL_BILL_ANALYZER_H
#define MATERIAL_BILL_ANALYZER_H

#include <string>
#include <map>
#include <vector>

#include "global_funcs.h"
#include "CsvTextBuilder.h"
#include "DetectorModule.h"
#include "Visitor.h"

namespace insur {

class InactiveElement;
class MaterialBudget;

/*
 * Calculate all contributions to the material budget for a given tracker, i.e.
 * build material bill.
 */
class MaterialBillAnalyzer {

  class ServiceElement;
  typedef std::map<std::string, double>      MaterialMap;
  typedef std::map<std::string, MaterialMap> LayerMaterialMap;
  typedef std::vector<ServiceElement>        ServicesMaterialVector;

  class ServiceElement {
  public: 
    double zmin, zmax, rmin, rmax;
    MaterialMap materialMap;
  };

  class MaterialVisitor : public ConstGeometryVisitor {

   public:
    virtual ~MaterialVisitor() {;}

    std::string getId()               { return m_id; }
    void visit(const BarrelModule& m) { m_id = m.cntName() + "_Layer" + any2str(m.layer()); }
    void visit(const EndcapModule& m) { m_id = m.cntName() + "_Disk"  + any2str(m.disk()); }

   private:
    std::string m_id;
  };

 public:
  // Constructor
  MaterialBillAnalyzer();
  // Calculate material bill for given tracker and fill internal Csv table with results
  void inspectTracker(insur::MaterialBudget&);

  // Test if text block exists
  bool existCsvText(std::string id);
  // Read csv text block content
  std::string getCsvText(std::string id) { return m_matBillCsv.getCsvText(id); }
  // Get the final Csv formatted date through iterators (ids) - "label" = label table, "tracker id" = table for each tracker
  std::map<std::string, std::string>::iterator getCsvTextBegin() { return m_matBillCsv.getCsvTextBegin();}
  std::map<std::string, std::string>::iterator getCsvTextEnd()   { return m_matBillCsv.getCsvTextEnd();}

 private:
  // Csv formatted output container
  CsvTextBuilder m_matBillCsv;
  // Layer material container
  LayerMaterialMap m_layerMaterialMap;
  // Services material container
  ServicesMaterialVector m_servicesMaterialVector;

  // Internal method getting material information about detector modules
  void inspectModules(std::vector<std::vector<insur::ModuleCap> >& tracker);
  // Internal method getting material information about services
  void inspectInactiveElements(std::string id, const std::vector<insur::InactiveElement>& inactiveElements);
}; // Class

}; // Namespace

#endif // MATERIAL_BILL_ANALYZER_H

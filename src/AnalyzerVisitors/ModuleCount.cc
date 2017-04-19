#include "AnalyzerVisitors/ModuleCount.hh"

#include <Tracker.hh>
#include <Barrel.hh>
#include <Endcap.hh>
#include <DetectorModule.hh>
#include <RootWeb.hh>


void ModuleCountVisitor::sortTypesAndDetectors() {
  int iType=1;
  int iSub=1;
  for (auto aType : moduleTypes) aType.second=iType++;
  for (auto aType : subDetectors) aType.second=iSub++;
}

RootWTable* ModuleCountVisitor::makeTable() {
  RootWTable* moduleCountTable = new RootWTable();
  int iType;
  int iSub;
  std::map<std::string, int> typeTotal;
  std::map<std::string, int> subTotal;
  int grandTotal=0;
  // Top left corner
  moduleCountTable->setContent(0, 0, "Module types");
  iType = 0;
  iSub = 0;
  // Row titles & column titles
  for (auto aType : moduleTypes) moduleCountTable->setContent(++iType, 0, aType.first);
  for (auto aSub : subDetectors) moduleCountTable->setContent(0, ++iSub, aSub.first);
  // The table central part
  iType=0;
  for (auto aType : moduleTypes) {
    iType++;
    iSub=0;	
    for (auto aSub : subDetectors) {
      iSub++;
      int cellContent = count_TypeSub[make_pair(aType.first, aSub.first)];
      if (cellContent!=0) moduleCountTable->setContent(iType, iSub, cellContent);
      typeTotal[aType.first]+=cellContent;
      subTotal[aSub.first]+=cellContent;
      grandTotal+=cellContent;
    }
  }
  // Total column and row titles
  moduleCountTable->setContent(++iType, 0, "Total");
  moduleCountTable->setContent(0, ++iSub, "Total");
  int totalType = iType;
  int totalSub = iSub;
  iType = 0;
  iSub = 0;
  // Total columns, rows and corner
  for (auto aType : moduleTypes) moduleCountTable->setContent(++iType, totalSub, typeTotal[aType.first]);
  for (auto aSub : subDetectors) moduleCountTable->setContent(totalType, ++iSub, subTotal[aSub.first]);
  moduleCountTable->setContent(totalType, totalSub, grandTotal);
      
  return moduleCountTable;
}

void ModuleCountVisitor::preVisit(const Tracker& tracker) {
  firstVisit = true;
  tracker.accept(*this);
  sortTypesAndDetectors();
  firstVisit = false;
}
    
void ModuleCountVisitor::visit(const Barrel& barrel) {
  currentSubdetector = barrel.myid();
  if (firstVisit) {
    subDetectors[currentSubdetector]=1;
  }
}

void ModuleCountVisitor::visit(const Endcap& endcap) {
  currentSubdetector = endcap.myid();
  if (firstVisit) {
    subDetectors[currentSubdetector]=1;
  }
}

void ModuleCountVisitor::visit(const DetectorModule& m) {
  std::string modType = m.summaryType();
  if (firstVisit) {
    moduleTypes[modType]=1;
  } else {
    auto index = make_pair(modType, currentSubdetector);
    count_TypeSub[index]++;
  }
}




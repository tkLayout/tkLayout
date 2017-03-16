#include <ReportModuleCount.hh>
#include <Tracker.h>
#include <Module.h>
#include <rootweb.hh>

#include <map>
#include <string>

std::string VisitorModuleCount::moduleSummaryType(const DetectorModule& m) const  {
  std::string result;
  result+=m.moduleType();
  if (m.dsDistance()!=0) result+=" "+any2str(m.dsDistance(), 1)+" mm";
  return result;
};

void VisitorModuleCount::sortTypesAndDetectors() {
  int iType=1;
  int iSub=1;
  for (auto aType : moduleTypes) aType.second=iType++;
  for (auto aType : subDetectors) aType.second=iSub++;
}

RootWTable* VisitorModuleCount::makeTable() {
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

void VisitorModuleCount::preVisit(const Tracker& tracker) {
  firstVisit = true;
  tracker.accept(*this);
  sortTypesAndDetectors();
  firstVisit = false;
}
    
void VisitorModuleCount::visit(const Barrel& barrel) {
  currentSubdetector = barrel.myid();
  if (firstVisit) {
    subDetectors[currentSubdetector]=1;
  }
}

void VisitorModuleCount::visit(const Endcap& endcap) {
  currentSubdetector = endcap.myid();
  if (firstVisit) {
    subDetectors[currentSubdetector]=1;
  }
}

void VisitorModuleCount::visit(const DetectorModule& m) {
  std::string modType = moduleSummaryType(m);
  if (firstVisit) {
    moduleTypes[modType]=1;
  } else {
    auto index = make_pair(modType, currentSubdetector);
    count_TypeSub[index]++;
  }
}


void ReportModuleCount::analyze(const Tracker& tracker) {
  moduleCounter_.preVisit(tracker);
  tracker.accept(moduleCounter_);  
}

void ReportModuleCount::visualizeTo(RootWContent& myContent) {
  myContent.addItem(moduleCounter_.makeTable());
}



#include <Module.hh>
#include <ReportModuleCount.hh>
#include <RootWeb.hh>
#include <Tracker.hh>

#include <map>
#include <string>

void ReportModuleCount::analyze(const Tracker &tracker) {
  moduleCounter_.preVisit(tracker);
  tracker.accept(moduleCounter_);
}

void ReportModuleCount::visualizeTo(RootWContent &myContent) {
  myContent.addItem(moduleCounter_.makeTable());
}

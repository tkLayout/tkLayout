#ifndef MODULESTOELINKSCONNECTOR_HH
#define MODULESTOELINKSCONNECTOR_HH

#include <global_constants.hh>
#include "global_funcs.hh"
#include "ITCabling/inner_cabling_functions.hh"
#include "ITCabling/ELink.hh"

#include "Module.hh"


/* This class is used to CONNECT MODULES TO ELINKS.
*/
class ModulesToELinksConnector : public GeometryVisitor {
public:
  std::map<std::string, ELink*> getELinks() { return eLinks_; }

  void visit(Barrel& b);
  void visit(Layer& l);     
  void visit(BarrelModule& m);

  void visit(Endcap& e);
  void visit(Ring& r);
  void visit(EndcapModule& m);

private:
  // BUILDING
  const int computeNumELinksPerModule(const std::string subdetectorName, const int layerOrRingNumber) const;

  void buildELink(DetectorModule& m, std::map<std::string, ELink*>& eLinks, const int eLinkIndex);
  const std::string computeELinkId(DetectorModule& m, const int eLinkIndex) const;
  ELink* createAndStoreELink(std::map<std::string, ELink*>& eLinks, DetectorModule& m, const std::string eLinkId);
  void connectModuleToELink(DetectorModule& m, ELink* eLink) const;

  std::map<std::string, ELink*> eLinks_;

  std::string barrelName_;
  int layerNumber_;

  std::string endcapName_;
  int ringNumber_;
};


#endif  // MODULESTOELINKSCONNECTOR_HH

#ifndef BAG_H
#define BAG_H

#include <map>
#include <string>
#include <set>

#include <TProfile.h>
#include <TGraph.h>
#include <TH2.h>

using std::string;

  class graphBag {
  public:
    static const double Triggerable;
    static const int RhoGraph;
    static const int PhiGraph;
    static const int DGraph;
    static const int CtgthetaGraph;
    static const int Z0Graph;
    static const int PGraph;
    static const int TriggeredGraph;
    static const int IdealGraph;
    static const int RealGraph;
    static const int TriggerGraph;
    static const int StandardGraph;
    std::map<double, TGraph>& getGraphs(const int& attribute);
    std::map<double, TGraph>& getTaggedGraphs(int attribute, const string& tag);
    const std::set<string>& getTagSet() const { return tagSet_; }
    int clearTriggerGraphs();
    int clearStandardGraphs();
    static int buildAttribute(bool ideal, bool isTrigger);
    //static std::pair<double, double> splitMomenta(double momentum);
    //static double joinMomenta(double momentum1, double momentum2);
  private:
    std::map<int, std::map<double, TGraph> > graphMap_;
    std::map<std::pair<int, string>, std::map<double, TGraph> > taggedGraphMap_;
    std::set<string> tagSet_;
    int clearGraphs(const int& attributeMask);
  };

  /**
   * @class mapBag
   * @brief A bag of graphs sorted by variable, scope and track's pt
   */
  class mapBag {
  public:
    static const int efficiencyMap;
    static const int thresholdMap;
    static const int thicknessMap;
    static const int windowMap;
    static const int suggestedSpacingMap;
    static const int suggestedSpacingMapAW;
    static const int nominalCutMap;
    static const int irradiatedPowerConsumptionMap;
    static const int totalPowerConsumptionMap;
    static const int moduleConnectionEtaMap;
    static const int moduleConnectionPhiMap;
    static const int moduleConnectionEndcapPhiMap;
    static const double dummyMomentum;
    std::map<double, TH2D>& getMaps(const int& attribute);
    int clearMaps(const int& attributeMask);
  private:
    std::map<int, std::map<double, TH2D> > mapMap_;
  };

  /**
   * @class profileBag
   * @brief A bag of TProfiles sorted by a double variable and scope
   */
  class profileBag {
  public:
    static const double Triggerable;
    static const int TriggeredProfile;
    static const int TriggerProfile;
    static const int TriggeredFractionProfile;
    static const int TriggerPurityProfile;
    static const std::string TriggerProfileName;
    static const std::string TriggerProfileNameWindow;
    static const std::string TurnOnCurveName;
    std::map<double, TProfile>& getProfiles(const int& attribute);
    int clearTriggerProfiles();
    int clearTriggerNamedProfiles();
    std::map<double, TProfile>& getNamedProfiles(const std::string& name);
    std::vector<std::string> getProfileNames(const std::string& name);
  private:
    int clearProfiles(const int& attributeMask);
    int clearNamedProfiles(const std::string& name);
    std::map<int, std::map<double, TProfile> > profileMap_;
    std::map<std::string, std::map<double, TProfile> > namedProfileMap_;
  };     

#endif

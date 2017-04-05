#ifndef _GEOMETRICINFO_HH
#define	_GEOMETRICINFO_HH

#include <global_constants.hh>
#include <VizardTools.hh>
#include <Tracker.hh>
#include <SimParms.hh>
#include <TagMaker.hh>
#include <RootWeb.hh>


    //**************************************//
    //*             Visitor                *//
    //*  Layers and disks : global info    *//
    //*                                    *//
    //**************************************//

    // Build the module type maps
    // with a pointer to a sample module
class LayerDiskSummaryVisitor : public ConstGeometryVisitor {
public:
  // info
  RootWTable* layerTable = new RootWTable();
  RootWTable* diskTable = new RootWTable();
  std::vector<RootWTable*> diskNames;
  std::vector<RootWTable*> ringTables;
  std::map<std::string, std::set<std::string> > tagMapPositions;
  std::map<std::string, int> tagMapCount;
  std::map<std::string, long> tagMapCountChan;
  std::map<std::string, double> tagMapMaxStripOccupancy;
  std::map<std::string, double> tagMapAveStripOccupancy;
  std::map<std::string, double> tagMapMaxHitOccupancy;
  std::map<std::string, double> tagMapAveHitOccupancy;
  std::map<std::string, double> tagMapAveRphiResolution;
  std::map<std::string, double> tagMapAveRphiResolutionRmse;
  std::map<std::string, double> tagMapSumXResolution;
  std::map<std::string, double> tagMapSumSquaresXResolution;
  std::map<std::string, double> tagMapCountXResolution;
  std::map<std::string, double> tagMapIsParametrizedXResolution;
  std::map<std::string, double> tagMapAveYResolution;
  std::map<std::string, double> tagMapAveYResolutionRmse;
  std::map<std::string, double> tagMapSumYResolution;
  std::map<std::string, double> tagMapSumSquaresYResolution;
  std::map<std::string, double> tagMapCountYResolution;
  std::map<std::string, double> tagMapIsParametrizedYResolution;
  std::map<std::string, double> tagMapAveRphiResolutionTrigger;
  std::map<std::string, double> tagMapAveYResolutionTrigger;
  std::map<std::string, double> tagMapSensorPowerAvg;
  std::map<std::string, double> tagMapSensorPowerMax;
  std::map<std::string, const DetectorModule*> tagMap;

  // counters
  int nBarrelLayers=0;
  int nEndcaps=0;
  int nDisks=0;
  int nRings=0;
  int totalBarrelModules = 0;
  int totalEndcapModules = 0;

  double totArea = 0;
  int totCountMod = 0;
  int totCountSens = 0;
  long totChannel = 0;
  double totalSensorPower = 0;

  double nMB;

  void preVisit();
  void visit(const Barrel& b) override;
  void visit(const Layer& l) override;
  void visit(const Endcap& e) override;
  void visit(const Disk& d) override;
  void visit(const Ring& r) override;
  void visit(const Module& m) override;
  void visit(const EndcapModule& m) override;
  void postVisit();
};


    //***************************************//
    //*                Visitor              *//
    //* Automatic-placement tilted layers : *//
    //*            Additional info          *//
    //*                                     *//
    //***************************************//

class TiltedLayersVisitor : public ConstGeometryVisitor {
public:
  // tilted info
  std::vector<RootWTable*> tiltedLayerNames;
  std::vector<RootWTable*> flatPartNames;
  std::vector<RootWTable*> tiltedPartNames;
  std::vector<RootWTable*> flatPartTables;
  std::vector<RootWTable*> tiltedPartTables;

  // counter
  int numTiltedLayers = 0;

  void visit(const Layer& l) override;     
};


    //************************************//
    //*               Visitor             //
    //*            AllModulesCsv          //
    //*                                   //
    //************************************//
class TrackerVisitor : public ConstGeometryVisitor {
  std::stringstream output_;
  string sectionName_;
  int layerId_;

public:
  void preVisit();
  void visit(const Barrel& b);
  void visit(const Endcap& e);
  void visit(const Layer& l);
  void visit(const Disk& d);
  void visit(const Module& m);
  std::string output() const { return output_.str(); }
};


    //************************************//
    //*               Visitor             //
    //*            BarrelModulesCsv       //
    //*                                   //
    //************************************//
class BarrelVisitor : public ConstGeometryVisitor {
  std::stringstream output_;
  string barName_;
  int layId_;
  int numRods_;

public:
  void preVisit();
  void visit(const Barrel& b);
  void visit(const Layer& l);
  void visit(const BarrelModule& m);
  std::string output() const;
};


    //************************************//
    //*               Visitor             //
    //*            EndcapModulesCsv       //
    //*                                   //
    //************************************//
class EndcapVisitor : public ConstGeometryVisitor {
  std::stringstream output_;
  string endcapName_;
  int diskId_;

public:
  void preVisit();
  void visit(const Endcap& e);
  void visit(const Disk& d);
  void visit(const EndcapModule& m);

  std::string output() const;
};


    //************************************//
    //*               Visitor             //
    //*            Sensors DetIds         //
    //*                                   //
    //************************************//
class TrackerSensorVisitor : public SensorGeometryVisitor {
  std::stringstream output_;
  string sectionName_;
  int layerId_;
  int moduleRing_;

public:
  void visit(Barrel& b);
  void visit(Endcap& e);
  void visit(Layer& l);
  void visit(Disk& d);
  void visit(Module& m);
  void visit(Sensor& s);

  std::string output() const;
};

#endif // _GEOMETRICINFO_HH

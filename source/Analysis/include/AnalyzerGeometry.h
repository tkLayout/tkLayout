/*
 * AnalyzerGeometry.h
 *
 *  Created on: 20. 4. 2016
 *      Author: Drasal (CERN)
 */

#ifndef INCLUDE_ANALYZERGEOMETRY_H_
#define INCLUDE_ANALYZERGEOMETRY_H_

// System libraries
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <set>

#include <Visitor.h>
#include <AnalyzerUnit.h>

// Forward declaration
class Barrel;
class BeamPipe;
class ConstGeometryVisitor;
class Detector;
class Disk;
class Endcap;
class EndcapModule;
class Layer;
class VisitorLayerName;
class DetectorModule;
class TCanvas;
class Tracker;
class TH2D;
class TH2I;
class TProfile;
class Ring;
class RootWTable;

/*
 * @class AnalyzerGeometry
 * Analyze geometry layout, vizualize data and print them out in a html formatted output.
 * Unique name defined as "AnalyzerGeometry".
 */
class AnalyzerGeometry : public AnalyzerUnit {

 public:

  //! Constructor
  AnalyzerGeometry(const Detector& detector);

  //! Destructor
  virtual ~AnalyzerGeometry();

  //! Initialize - mostly histograms & other containers
  //! @return True if OK
  bool init(int nGeomTracks);

  //! Inspect geometry layout (if init OK) -> collect data to histograms & tables
  //! @return True if OK
  bool analyze();

  //! Visualize geometry layout (if init & analysis OK) -> add html page with collected tables & created histograms
  //! @return True if OK
  bool visualize(RootWSite& webSite);

  //! Get number of used tracks
  int getNGeomTracks() const { return m_nTracks;}

 private:

  //! Draw beam-pipe in RZ to given canvas
  void drawBeamPipeRZ(TCanvas& canvas, double maxZ);
  //! Draw beam-pipe in XY to given canvas
  void drawBeamPipeXY(TCanvas& canvas);

  int    m_nTracks;     //!< Number of geometry tracks to be used in the analysis
  double m_etaSpan;     //!< Eta interval to be analyzed
  double m_etaMin;      //!< Minimum eta value;
  double m_etaMax;      //!< Maximum eta value

  const float c_etaSafetyMargin = 0.01;
  const int   c_nBinsProfile    = 100;

  std::unique_ptr<VisitorLayerName> m_layerNamesVisitor; //! Visitor pattern to be used to find layer names

  // Histogram output
  std::map<std::string,TH2D>      m_hitMapPhiEta;            //!< Number of hits - map in phi & eta for each tracker with given name
  std::map<std::string,TH2I>      m_trackMapPhiEta;          //!< Number of tracks - map in phi & eta for each tracker with given name

  std::map<std::string,TProfile>  m_moduleHitEtaProfile;     //!< Number of hits in modules versus eta for each tracker with given name
  std::map<std::string,TProfile>  m_sensorHitEtaProfile;     //!< Number of hits in sensors versus eta for each tracker with given name
  std::map<std::string,TProfile>  m_stubHitEtaProfile;       //!< Number of hits in stubs versus eta for each tracker with given name

  std::map<std::string,std::set<std::string>>           m_moduleTypes;             //!< List of module types in a given tracker
  std::map<std::string,std::map<std::string,int>>       m_moduleTypeColor;         //!< Find color to given module type for a given tracker
  std::map<std::string,std::map<std::string,TProfile>>  m_moduleTypeHitEtaProfile; //!< Number of hits in various types of modules versus eta for each tracker with given name
  std::map<std::string,std::map<std::string,TProfile>>  m_moduleTypeStubEtaProfile;//!< Number of stubs in various types of modules versus eta for each tracker with given name

  std::map<std::string,std::map<std::string, TProfile>> m_layerEtaCoverProfile;    //!< For each tracker: Eta coverage profile of individual layers defined by unique name
  std::map<std::string,std::map<std::string, TProfile>> m_layerStubEtaCoverProfile;//!< For each tracker: Eta coverage profile of individual stub layers defined by unique name

}; // Class

/*
 * Helper class: Layer name visitor (visitor pattern) -> get names of individual layers
 */
class VisitorLayerName : public ConstGeometryVisitor {

  std::string m_idBRLorEC; //!< Barrel/Endcap id number
  std::string m_idTRK;     //!< Tracker name

  std::map<std::string, std::set<std::string>> m_data; //!< Layer/Disk names for given tracker

 public:

  VisitorLayerName(std::vector<const Tracker*>& trackers);
  virtual ~VisitorLayerName() {};

  //! Fill container with layer names for defined tracker if tracker exists
  bool getLayerNames(std::string trkName, std::set<std::string>& layerNames);

  void visit(const Barrel& b);
  void visit(const Endcap& e);
  void visit(const Layer& l);
  void visit(const Disk& d);
}; // Helper Class

/*
 *  Helper class: Layer/disk summary visitor (visitor pattern) - gather information for geometry tables
 */
class VisitorLayerDiscSummary : public ConstGeometryVisitor {

 public:

  virtual ~VisitorLayerDiscSummary();

  void preVisit();
  void visit(const Layer& l) override;
  void visit(const Disk& d) override;
  void visit(const Ring& r) override;
  void visit(const DetectorModule& m) override;
  void visit(const EndcapModule& m) override;
  void postVisit();

  std::unique_ptr<RootWTable> m_layerTable; //!< Web table containing info about layers
  std::unique_ptr<RootWTable> m_diskTable;  //!< Web table containing info about disks
  std::unique_ptr<RootWTable> m_ringTable;  //!< Web table containing info about rings
  std::unique_ptr<RootWTable> m_moduleTable;//!< Web table containing info about modules

  // Counters
  int m_nBarrelLayers      = 0; //!< Number of barrel layers
  int m_nDisks             = 0; //!< Number of disks
  int m_nRings             = 0; //!< Number of rings
  int m_totalBarrelModules = 0; //!< Total number of barrel modules
  int m_totalEndcapModules = 0; //!< Total number of end-cap modules

  double m_totalArea       = 0; //!< Total tracker area
  int    m_totalModules    = 0; //!< Total number of modules
  int    m_totalSensors    = 0; //!< Total number of sensors
  long   m_totalChannels   = 0; //!< Total number of channels
  double m_totalSensorPower= 0; //!< Total power needed

  std::map<std::string, std::set<std::string> > m_moduleTagToPositionsMap;

  std::map<std::string, const DetectorModule*> m_modulePtrMap;            //!< Module (by tag) to module pointer map -> to get module properties
  std::map<std::string, int>                   m_moduleCount;             //!< Number of modules of given module type (by tag)
  std::map<std::string, long>                  m_moduleChannels;          //!< Number of channels of given module type
  std::map<std::string, long>                  m_moduleAvgChannelsRPhi;   //!< Number of channels in R-Phi of given module type
  std::map<std::string, long>                  m_moduleAvgChannelsZ;      //!< Number of channels in Z of given module type
  std::map<std::string, long>                  m_moduleAvgROCs;           //!< Number of read-out chips per module type
  std::map<std::string, double>                m_moduleMaxStripOccupancy; //!< Maximum strip occupancy of given module type
  std::map<std::string, double>                m_moduleAvgStripOccupancy; //!< Average hit occupancy of given module type
  std::map<std::string, double>                m_moduleMaxHitOccupancy;   //!< Maximum hit occupancy of given module type
  std::map<std::string, double>                m_moduleAvgHitOccupancy;   //!< Average hit occupancy of given module type
  std::map<std::string, double>                m_moduleMinRphiResolution; //!< Minimum R-phi resolution of given module type
  std::map<std::string, double>                m_moduleAvgRphiResolution; //!< Average R-phi resolution of given module type
  std::map<std::string, double>                m_moduleMaxRphiResolution; //!< Maximum R-phi resolution of given module type
  std::map<std::string, double>                m_moduleMinZResolution;    //!< Minimum Z resolution of given module type
  std::map<std::string, double>                m_moduleAvgZResolution;    //!< Average Z resolution of given module type
  std::map<std::string, double>                m_moduleMaxZResolution;    //!< Maximum Z resolution of given module type
  //std::map<std::string, double>                m_moduleMinRphiPitch;      //!< Minimum R-phi pitch of given module type
  //std::map<std::string, double>                m_moduleMaxRphiPitch;      //!< Maximum R-phi pitch of given module type
  //std::map<std::string, double>                m_moduleMinZPitch;         //!< Minimum Z pitch of given module type
  //std::map<std::string, double>                m_moduleMaxZPitch;         //!< Maximum Z pitch of given module type
  std::map<std::string, double>                m_moduleAvgPower;          //!< Average power for given module type
  std::map<std::string, double>                m_moduleMaxPower;          //!< Maximum power for given module type

  std::map<int, const EndcapModule*> m_ringModuleMap; //!< Ring (by id) module map -> to get ring properties
  std::vector<int>                   m_ringNModules;  //!< Number of modules in a given ring

 private:
  const int   c_coordPrecision       = 1;
  const int   c_areaPrecision        = 1;
  const int   c_occupancyPrecision   = 1;
  const int   c_resolutionPrecision  = 1;
  const int   c_channelPrecision     = 2;
}; // Helper Class


/*
 *  Helper class: Visit tilted layers for info on website.
 */
class TiltedLayersVisitor : public ConstGeometryVisitor {
 
 public:

  void visit(const Tracker& t) override;
  void visit(const Layer& l) override;

  // tilted info containers
  std::vector<std::unique_ptr<RootWTable>> m_tiltedLayerNames;
  std::vector<std::unique_ptr<RootWTable>> m_flatPartTables;
  std::vector<std::unique_ptr<RootWTable>> m_tiltedPartTables;

  // counters
  int m_nTiltedLayers;
  int m_nLayers;

 private:
  const int   c_tiltedCoordPrecision = 2;
  const int   c_zOverlapPrecision    = 3;
  const int   c_anglePrecision       = 1;
}; // Helper Class

#endif /* INCLUDE_ANALYZERGEOMETRY_H_ */

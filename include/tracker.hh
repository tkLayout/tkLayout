#ifndef _TRACKER_HH_
#define _TRACKER_HH_

#include <vector>
#include "Math/Vector3D.h"
#include <string>
#include "TGeoManager.h"
#include "module.hh"
#include "moduleType.hh"
#include "layer.hh"
#include <Palette.h>

#include "TCanvas.h"
#include "TRandom3.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"

#include <messageLogger.h>

#define RANDOM_SEED 0xcaffe

#define DEFAULTBARRELNAME "defaultBarrel"
#define DEFAULTENDCAPNAME "defaultEndcap"
#define COLOR_PLOT_BACKGROUND kWhite
#define COLOR_BACKGROUND kGray
#define COLOR_GRID kGreen-10
#define COLOR_HARD_GRID kGray

#define COLOR_INVALID_MODULE kGray+1

// TODO: add slanted gap between barrel and end-cap

using namespace ROOT::Math;

typedef std::vector<Layer*> LayerVector;
typedef std::map<std::string, LayerVector> SectionMap;
typedef std::map<int, double> SpecialDelta;

class Tracker : public MessageLogger {
public:
  enum {TypeBarrel, TypeEndcap};
protected:
  LayerVector layerSet_;
  LayerVector barrelLayerSet_;
  LayerVector endcapLayerSet_;
  ModuleVector endcapSample_;
  SectionMap sectionMap_;
  double nMB_;

  std::map<int, double> mapTypeToCost_;
  std::map<std::string, ModuleType> mapType_;
  std::map<std::pair<int, int>, double> mapIrradiation_;

  SpecialDelta specialSmallDelta_;
  SpecialDelta specialBigDelta_;

  double maxL_;
  double maxR_;

  std::string arguments_;

  // Default variables (distances in mm)
  static const double defaultZError_ = 70.;    // Vertex displacement sigma in z (7 cm)
  //static const double defaultRError_ = 0.02;   // Vertex displacement sigma in r (20 um)
  static const double defaultRError_ = 5.00;   // Vertex displacement sigma in r (5 mm)
  static const double defaultSmallDelta_ = 2.; // Space between overlapping modules
  static const double defaultBigDelta_ = 12.;  // Space between different faces of the same structure
  static const double defaultOverlap_ = 1.;    // Safety overlap between modules
  static const double defaultNMB_ = 200;
  static const bool defaultUseIPConstraint_ = true; // in the trigger fit use the ip constraint by default
  
private:
  void setDefaultParameters();
  void shapeVolume();
  void shapeLayerVolumes();

  void shapeModuleVolumes(bool lite = false, int section = Layer::NoSection);
  void shapeModuleVolumesEndcapSample(bool lite = false);

  void placeModule(Module* aModule);
  void placeModuleLite(Module* aModule);

  // Color picking according to type
  //std::map<std::string, Color_t> colorPickMap_; // obsolete
  //Color_t lastPickedColor_;                     // obsolete

  int iModule_;
  double efficiency_;
  double pixelEfficiency_;

  std::string comment_;

  TRandom3 myDice_; // obsolete

  TGeoVolume* myVolume_;
  TGeoMedium* myMed_;
  TGeoManager* myGeom_;

  TCanvas* geomLite_;   // should be made obsolete
  TCanvas* geomLiteXY_; // should be made obsolete
  TCanvas* geomLiteYZ_; // should be made obsolete
  TCanvas* geomLiteEC_; // should be made obsolete
  TCanvas* etaProfileCanvas_; // obsolete
  TCanvas* bandWidthCanvas_;

  TH1F* bandWidthDist_;
  TH1F* bandWidthDistSp_;
  TH1F* chanHitDist_;

  TH2D* mapPhiEta_; // obsolete

  std::vector<TObject* > savingV_; // (partly?) obsolete: check

  double rError_;
  double zError_;
  double useIPConstraint_;
  double smallDelta_;
  double bigDelta_;
  double overlap_;
  double etaCut_;
  int phiSegments_;

  double numInvFemtobarns_;
  double operatingTemp_;
  double referenceTemp_;
  int chargeDepletionVoltage_;
  double alphaParam_;

  std::string summaryDirectory_;
  std::string storeDirectory_;
  std::string trackerName_;
  std::string activeDirectory_;

  std::map<int, int> ringDirectives_;
  std::map<int,double> layerDirectives_;
  std::map<int,LayerOption> layerOptions_;

  // Color picking
  //Color_t colorPicker(std::string); // obsolete

  // Internal stats
  std::vector<int> LpB_, DpE_;
  std::vector<double> rMinpB_, rMaxpB_, dZpB_, rMinpE_, rMaxpE_, dZpE_;

  // Geometry validation functions
  ModuleVector trackHit(const XYZVector& origin, const XYZVector& direction, ModuleVector* properModules); // obsolete
  void resetTypeCounter(std::map <std::string, int> &modTypes); // obsolete
  int createResetCounters(std::map <std::string, int> &modTypes); // obsolete
  std::pair <XYZVector, double > shootDirection(double minEta, double maxEta); // obsolete
  std::pair <XYZVector, double > shootFixedDirection(double phi, double eta);  // obsolete
  std::pair <XYZVector, double > shootDirectionFixedPhi(double minEta, double maxEta); // obsolete
  
  // Formatted output
  void printHtmlTableRow(std::ofstream *output, std::vector<std::string> myRow);
  void printHtmlTableRow(std::ofstream *output, std::vector<double> myRow, int coordPrecision = 0, bool skimZero = false);
public:
  void compressBarrelLayers(LayerVector aLayerSet, bool oneSided, double destZ = 0);
private:
  void createDirectories();

  enum {ViewSectionXY=3, ViewSectionYZ=1, ViewSectionXZ=2}; // obsolete
  void drawTicks(TView* myView, double maxL, double maxR, int noAxis=1, double spacing = 100., Option_t* option = "same"); // shold become obsolete
  void drawGrid(double maxL, double maxR, int noAxis=1, double spacing = 100., Option_t* option = "same"); // shold become obsolete
  void drawSummary(double maxZ, double maxRho, std::string fileName); // obsolete
  void drawLayout(double maxZ, double maxRho, std::string fileName);  // obsolete

public:
  ~Tracker();
  Tracker();
  Tracker(std::string trackerName);

  std::vector<int>& getLayersPerBarrel() { return LpB_; }
  std::vector<int>& getDiscsPerEndcap() { return DpE_; }
  std::vector<double>& getMinRPerBarrel() { return rMinpB_; }
  std::vector<double>& getMaxRPerBarrel() { return rMaxpB_; }
  std::vector<double>& getdZPerBarrel() { return dZpB_; }
  std::vector<double>& getMinRPerEndcap() { return rMinpE_; }
  std::vector<double>& getMaxRPerEndcap() { return rMaxpE_; }
  std::vector<double>& getdZPerEndcap() { return dZpE_; }

  // Standard barrel builder
  LayerVector buildBarrel(int nLayer, double minRadius, double maxRadius,
		   int nModules, BarrelModule* sampleModule, std::string barrelName, int section = Layer::NoSection,
		   bool compressed = false, bool shortBarrel = false );
  // Barrel builder for backwards compatibility with the command-line version
/*  void buildBarrel(int nLayer, double minRadius, double maxRadius,
		   int nModules, BarrelModule* sampleModule, int section = Layer::NoSection,
		   bool compressed = false ); */

  // Adjustment for short barrels
  void alignShortBarrels();

  // Sort barrel and endcap layer collections by rho and z, respectively
  void sortLayers();


  // Standard endcap builder
  void buildEndcaps(int nDisks, double minZ, double maxZ, double minRadius, double maxRadius,
		    std::map<int, EndcapModule*> sampleModule, std::string endcapName, int diskParity,
		    bool oddSegments, bool alignEdges,
		    int sectioned = Layer::NoSection ); 
  // Gets the minimum radius from eta and minimum Z
  void buildEndcapsAtEta(int nDisks, double minZ, double maxZ, double maxEta, double maxRadius,
			 std::map<int, EndcapModule*> sampleModule, std::string endcapName, int diskParity,
			 bool oddSegments, bool alignEdges,
			 int sectioned = Layer::NoSection ); 

  // Endcap ring remover
  void removeDiskRings(std::string sectionName, int iDisk, int iRing, bool directionOuter);


  // Access to parameters
  void setRError(const double& newError) { rError_ = newError; }
  void setZError(const double& newError) { zError_ = newError; }
  void setUseIPConstraint(bool newUse) { useIPConstraint_ = newUse; }
  void setBigDelta(const double& newDelta) { bigDelta_ = newDelta; }
  void setSmallDelta(const double& newDelta) { smallDelta_ = newDelta; }
  void setSpecialSmallDelta(const int& specialIndex, const double& newSpecialSmallDelta) { specialSmallDelta_[specialIndex]=newSpecialSmallDelta; };
  void setSpecialBigDelta(const int& specialIndex, const double& newSpecialBigDelta) { specialBigDelta_[specialIndex]=newSpecialBigDelta; };
  double getSpecialSmallDelta(const int& specialIndex) { return specialSmallDelta_[specialIndex]; };
  double getSpecialBigDelta(const int& specialIndex) { return specialBigDelta_[specialIndex]; };
  void resetSpecialSmallDelta() { while (!specialSmallDelta_.empty()) specialSmallDelta_.erase(specialSmallDelta_.begin()); };
  void resetSpecialBigDelta() { while (!specialBigDelta_.empty()) specialBigDelta_.erase(specialBigDelta_.begin()); };
  void resetSpecialDeltas() { resetSpecialSmallDelta(); resetSpecialBigDelta(); };
  void setOverlap(const double& newOverlap) { overlap_ = newOverlap; };
  void setEtaCut(const double& newEta) { etaCut_ = newEta; };

  void setPhiSegments(const int& newPhiSegments) { phiSegments_ = newPhiSegments; };

  void setStoreDirectory(const std::string newDir) { storeDirectory_ = newDir; };
  void setSummaryDirectory(const std::string newDir) { summaryDirectory_ = newDir; };
  void setActiveDirectory(const std::string newDir) { activeDirectory_ = newDir; };
  void setTrackerName(const std::string newName) { trackerName_ = newName; }; // deprecated: TODO remove it
  void setName(const std::string newName) { trackerName_ = newName; }; // deprecated: TODO remove it
  void setLayerDirectives(const std::map<int, double> newDirectives ) { layerDirectives_=newDirectives; };
  void setLayerOptions(const std::map<int, LayerOption> newOptions ) { layerOptions_=newOptions; };
  void setRingDirectives(const std::map<int, int> newDirectives ) { ringDirectives_=newDirectives; };
  void setArguments(const std::string &newArgs) {arguments_=newArgs;};
  void setComment(const std::string &newComment) {comment_=newComment;};
  void setNMB(const double nMB) {nMB_=nMB;};
  void setSparsifiedHeaderBits(string typeIndex, int bits)  { mapType_[typeIndex].setSparsifiedHeaderBits(bits);  }
  void setSparsifiedPayloadBits(string typeIndex, int bits) { mapType_[typeIndex].setSparsifiedPayloadBits(bits); }
  void setTriggerDataHeaderBits(string typeIndex, int bits)  { mapType_[typeIndex].setTriggerDataHeaderBits(bits);  }
  void setTriggerDataPayloadBits(string typeIndex, int bits) { mapType_[typeIndex].setTriggerDataPayloadBits(bits); }
  void setSensorThickness(string typeIndex, double thickness) { mapType_[typeIndex].setSensorThickness(thickness); }
  void setInefficiencyType(string typeIndex, ptError::InefficiencyType type) { mapType_[typeIndex].setInefficiencyType(type); }
  void setNumInvFemtobarns(double numInvFemtobarns) { numInvFemtobarns_ = numInvFemtobarns; }
  void setOperatingTemp(double operatingTemp) { operatingTemp_ = operatingTemp; }
  void setReferenceTemp(double referenceTemp) { referenceTemp_ = referenceTemp; }
  void setChargeDepletionVoltage(int chargeDepletionVoltage) { chargeDepletionVoltage_ = chargeDepletionVoltage; }
  void setAlphaParam(double alphaParam) { alphaParam_ = alphaParam; }


  // Summary parameters
  double getCost(const int& type) { return(mapTypeToCost_[type]); }; // should be made obsolete (should go into Analyzer)
  //double getPower(const int& type) { return(mapTypeToPower_[type]); }; // should be made obsolete (should go into Analyzer)
  void setCost(const int& type, const double& newCost) { mapTypeToCost_[type]=newCost; }; // should be made obsolete (should go into Analyzer)
  //void setPower(const int& type, const double& newPower) { mapTypeToPower_[type]=newPower; }; // should be made obsolete (should go into Analyzer)
  void setPower(string typeIndex, int powerIndex, double newPower) { mapType_[typeIndex].setPowerPerStrip(newPower, powerIndex); }
  void setPowerFixed(string typeIndex, int powerIndex, double newPower) { mapType_[typeIndex].setPowerPerModule(newPower, powerIndex); }
  std::map<std::string, ModuleType>& getTypes() { return mapType_ ; }
  ModuleType& getModuleType(std::string typeName) { return mapType_[typeName] ; }


  // Trigger module types error increase
  void setTriggerErrorX(string typeIndex, double errorIncrease) { mapType_[typeIndex].setTriggerErrorX(errorIncrease); }
  void setTriggerErrorY(string typeIndex, double errorIncrease) { mapType_[typeIndex].setTriggerErrorY(errorIncrease); }

  // Overlaps / error
  double getRError() { return rError_; }
  double getZError() { return zError_; }
  bool getUseIPConstraint() { return useIPConstraint_; }
  double getSmallDelta() { return smallDelta_; }
  double getBigDelta() { return bigDelta_; }
  double getSmallDelta(const int& index);
  double getBigDelta(const int& index);
  double getOverlap() { return overlap_; }
  double getEtaCut() { return etaCut_; }

  // Module efficiency
  double getEfficiency() {return efficiency_ ; }
  void setEfficiency(double newEfficiency) { efficiency_ = newEfficiency ; }

  // Pixel module efficiency
  double getPixelEfficiency() {return pixelEfficiency_ ; }
  void setPixelEfficiency(double newEfficiency) { pixelEfficiency_ = newEfficiency ; }

  // Access to layer vectors
  LayerVector* getBarrelLayers() { return &barrelLayerSet_; }
  LayerVector* getEndcapLayers() { return &endcapLayerSet_; }
  LayerVector& getLayers() { return layerSet_; }
  ModuleVector& getEndcapSample() { return endcapSample_ ; }

  // Other
  std::string getStoreDirectory() { return storeDirectory_; };
  std::string getSummaryDirectory() { return summaryDirectory_; };
  std::string getActiveDirectory() { return activeDirectory_; };
  std::string getTrackerName() { return trackerName_; }; // deprecated (TODO: remove it)
  std::string getName() { return trackerName_; }; // deprecated (TODO: remove it)  
  std::string getArguments() {return arguments_;};
  std::string getComment() {return comment_;};
  double getNMB() {return nMB_;};
  double getMaxL() {return maxL_;};
  double getMaxR() {return maxR_;};
  TCanvas* getGeomLite() {return geomLite_;};
  TCanvas* getGeomLiteXY() {return geomLiteXY_;};
  TCanvas* getGeomLiteYZ() {return geomLiteYZ_;};
  TCanvas* getGeomLiteEC() {return geomLiteEC_;};

  int  getSparsifiedHeaderBits(string typeIndex)  { return mapType_[typeIndex].getSparsifiedHeaderBits();  }
  int  getSparsifiedPayloadBits(string typeIndex) { return mapType_[typeIndex].getSparsifiedPayloadBits(); }
  int  getTriggerDataHeaderBits(string typeIndex)  { return mapType_[typeIndex].getTriggerDataHeaderBits();  }
  int  getTriggerDataPayloadBits(string typeIndex) { return mapType_[typeIndex].getTriggerDataPayloadBits(); }
  double  getSensorThickness(string typeIndex) { return mapType_[typeIndex].getSensorThickness(); }
  ptError::InefficiencyType getInefficiencyType(string typeIndex) { return mapType_[typeIndex].getInefficiencyType(); }
  std::map<std::pair<int,int>, double>& getIrradiationMap() { return mapIrradiation_; } 
  double getNumInvFemtobarns() const { return numInvFemtobarns_; }
  double getOperatingTemp() const { return operatingTemp_; }
  double getReferenceTemp() const { return referenceTemp_; }
  int getChargeDepletionVoltage() const { return chargeDepletionVoltage_; }  
  double getAlphaParam() const { return alphaParam_; }


  void addLayer(Layer* aLayer, std::string sectionName, int type = TypeBarrel) {
    layerSet_.push_back(aLayer);
    if (type==TypeBarrel) {
      barrelLayerSet_.push_back(aLayer);
    }
    if (type==TypeEndcap) {
      endcapLayerSet_.push_back(aLayer);
    }

    sectionMap_[sectionName].push_back(aLayer);
  };
  
  // 3D geometry output preparation
  void createGeometry(bool lite = false); // should be made obsolete (currently used by analyzer)  or be restricted in use: no canvases

  // Summary output
  //void writeSummary(bool configFiles, std::string configFile, std::string dressFile, std::string fileType = "html",
  //                std::string barrelModuleCoordinatesFile ="",
  //                std::string endcapModuleCoordinatesFile ="");
  //void writeSummary(std::string fileType = "html"); // obsolete
  //void createPackageLayout(std::string dirName);
  void printBarrelModuleZ(ostream& outfile); // to be made obsolete
  void printEndcapModuleRPhiZ(ostream& outfile); // to be made obsolete

  // Save everything
  void save();

  // Geometry validation
  std::pair<double, double> getEtaMinMax();
  //void analyze(int nTracks = 1000, int section = Layer::NoSection); // obsolete
  int cutOverEta(double etaCut);
  double getMaxBarrelZ(int direction);
  //  void compressBarrelLayers();

  // Module adjustments
  void changeRingModules(std::string diskName, int ringN, std::string newtype, Color_t newColor);

  //void setModuleTypes();
  void setModuleTypes(std::string sectionName,
		      std::map<int, int> nStripsAcross,
		      std::map<int, int> nFaces,
		      std::map<int, int> nSegments,
		      std::map<int, std::string> myType, 
                      std::map<int, double> dsDistance,
                      std::map<int, int> triggerWindow,
                      std::map<int, double> dsRotation,
		      std::map<int, int> divideBack,
		      std::map<int, double> xResolution,
		      std::map<int, double> yResolution,
		      std::map<std::pair<int, int>, int> nStripsAcrossSecond,
		      std::map<std::pair<int, int>, int> nFacesSecond,
		      std::map<std::pair<int, int>, int> nSegmentsSecond,
		      std::map<std::pair<int, int>, std::string> myTypeSecond,
                      std::map<std::pair<int, int>, double> dsDistanceSecond,
                      std::map<std::pair<int, int>, int> triggerWindowSecond,
                      std::map<std::pair<int, int>, double> dsRotationSecond,
		      std::map<std::pair<int, int>, int> divideBackSecond,
		      std::map<std::pair<int, int>, bool> specialSecond);

};


#endif

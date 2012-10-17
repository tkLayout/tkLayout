#ifndef _MODULE_HH_
#define _MODULE_HH_

// Standard stuff
#include <vector>
#include <time.h>
#include <string>

// ROOT objects
#include "Math/Vector3D.h"
#include "TGeoManager.h"
#include "TPolyLine3D.h"

// Our objects
#include <messageLogger.h>
#include <moduleType.hh>
#include <ptError.h>

/*****************************/
/*                           */
/* Measurements are  in mm   */
/* It's simple to remember:  */
/* it's an  IS measurement!  */
/*                           */
/*****************************/

// TODO: add the power consumption estimate
// constant (digital: per module) + 
// scaling (analogue? : per strip / strip size..?) - Mark Raymond

#define maxNFaces 2

using namespace ROOT::Math;

typedef std::pair<double, int> edge;

template<class T> struct InitableProperty {
  T value;
  bool inited;
  InitableProperty() : inited(false) {}
  InitableProperty(const T& v) : value(v), inited(true) {}
};


struct PosRef {
  int8_t cnt;
  int8_t z;
  int8_t rho;
  int8_t phi;
};



class Module {
  friend ostream& operator<<(ostream& output, const Module& m);

protected:
  int shape_;
  double aspectRatio_;

  double boundaryMinPhi_;
  double boundaryMaxPhi_;
  double boundaryMinEta_;
  double boundaryMaxEta_;
  std::map<int, double> optimalSpacing_;
  bool isAcrossPi_;
  bool computedBoundaries_;

  double waferDiameter_;
  
  // Shape-specific paramters
  double height_;
  double thickness_;
  double moduleThickness_;
  double area_;
  double stripArea_;
  double stereodist_;
  int triggerWindow_;
  double stereorot_;
  double dphideta_;
  double phiWidth_, etaWidth_;
  int nSegmentsFace_[maxNFaces];
  int nStripAcross_;
  int nFaces_;
  int readoutType_;
  int readoutMode_;
  double resolutionRphi_;
  double resolutionY_;

  double irradiatedPowerConsumption_;

  int processorConnectionsEta_, processorConnectionsPhi_;

  int inSection_;

  int ring_;
  
  int phiIndex_;

  int containerId_;
  int zSide_; // whether in the positive or negative Z side of the detector

  std::string containerName_;
  
  //  point corner_[4];
  XYZVector corner_[4];

  double phi_; // A che serve?
  double r_;
  int nHits_;

  // Default variables
  static const double  defaultWaferDiameter_ = 131.; // Wafer diameter 131 mm
  static const double  defaultAspectRatio_ = 1.;
  static const double  defaultThickness_ = 0.3; // Wafer thickness: 300 um
  static const double  defaultModuleThickness_ = 1.0; //total module thickness
  static const Color_t defaultColor_ = kBlack;
  static const int     defaultNHits_ = 0;
  static const int     defaultHeight_ = 0;
  static const int     defaultArea_ = 0;
  static const int     defaultStereoDist_ = 2; // default distance between sensors = 2 mm
  static const int     defaultStereoRot_ = 0;
  static const int     defaultInSection_ = 0;
  static const int     defaultChannelsPerFace_ = 1;
  static const int     defaultSegments_ = 1;
  static const int     defaultStripAcross_ = 1;
  static const int     defaultFaces_ = 1;
  static const double  defaultResolutionRphi_ = 0.0;
  static const double  defaultResolutionY_ = 0.0;
  
  std::string id_;   // Ids of the module
  std::string tag_;  // Tags the module
  std::string type_; // Specifies the module type
  Color_t color_;
  
  // Returns the gamma parameter of the line if hits (complete line)
  double triangleCross(const XYZVector& P1, // Triangle points
		       const XYZVector& P2,
		       const XYZVector& P3,
		       const XYZVector& PL, // Base line point
		       const XYZVector& PU);

  TGeoVolume* thisVolume_;

  edge getEdgeRho(int direction);

  //void computeStripArea(int nFace);
  void computeMaxDphiDeta();

  const ModuleType* moduleType_;

  ptError myPtError;

/* mutable property bank. values are cached for optimized performance and simply returned instead of being recalculated once geometry is locked */
  mutable InitableProperty<double> minTheta_, maxTheta_, meanTheta_;
  mutable InitableProperty<double> minRho_, maxRho_;
  mutable InitableProperty<double> minZ_, maxZ_;
  mutable InitableProperty<double> minPhi_, maxPhi_;
  mutable InitableProperty<XYZVector> meanPoint_;
  mutable InitableProperty<int> octant_;

 private:
  void setDefaultParameters();
  int findMaxSegmentsFace_(); // starting from 0 !!
  int findMinSegmentsFace_(); // starting from 0 !!
    
 public:
  virtual ~Module();
  Module();
  Module(double waferDiameter);

  virtual int getSubdetectorType() const { return Undefined; };
  int getShape() const { return shape_; };
  
  void translate(XYZVector Delta);
  void rotatePhi(double phi);
  void rotateX(double phi);
  void rotateY_PI();
  void reflectZ();
  void shiftRho(double Delta);

  XYZVector* marginBorder(double widthMargin, double lengthMargin, int corner );
  XYZVector* marginBorderSide(double margin, int side);

  void print();
  
  void setOptimalSpacing(const int& windowSize, const double& newSpacing) { optimalSpacing_[windowSize] = newSpacing ; }
  const double& getOptimalSpacing(const int& windowSize) { return optimalSpacing_[windowSize] ; }

  int projectSideZ(int side, double destZ, double displaceZ = 0);
  int projectSideRho(int side, double destRho, double displaceZ = 0);
  
  // Returns the distance if hits or -1
  double trackCross (const XYZVector& PL,  // Base line point
		     const XYZVector& PU); // Base direction

  void setId(const std::string newId) {id_=newId;}
  std::string getId() {return id_;}

  void setTag(const std::string newTag) {tag_=newTag;}
  std::string getTag() {return tag_;}

  void setType(const std::string newType) {type_=newType;}
  std::string getType() {return type_;}

  void setModuleType(const ModuleType* newType) { moduleType_ = newType; }
  const ModuleType* getModuleType() const { return moduleType_; }

  void setStereoDistance(double sdist) { stereodist_=sdist; }
  double getStereoDistance() const { return (nFaces_ -1) * stereodist_; }; // TODO: check if this creates any problem around!

  void setTriggerWindow(const int& newWindow) { triggerWindow_ = newWindow ; }
  const int& getTriggerWindow() { return triggerWindow_ ; }

  void setStereoRotation(double srot) { stereorot_=srot; }
  double getStereoRotation() const { return stereorot_; };

  void setColor(const Color_t newColor) {color_ = newColor;}
  Color_t getColor() const {return color_;}

  void resetNHits() { nHits_ = defaultNHits_ ; };
  int getNHits() const { return nHits_ ;};

  void shapeVolume(TGeoVolume* container,
		   TGeoMedium* medium,
		   TGeoManager* geom);

  TPolyLine3D* getContour();

  double getHeight() const {return height_;};
  double getArea() const { return area_;};
  double getDiameter() const {return waferDiameter_; };
  double getThickness() const { return thickness_; }; // TODO: important: Check the use of "thickness" everywhere. it might have been confused with 'stereodistance'  !!!DEPRECATED!!! USE getSensorThickness in moduleType
  double getModuleThickness() const { return moduleThickness_; };
  const XYZVector& getCorner(int index) const { return corner_[index]; };

  edge getEdgeRhoSide(int direction);
  int setEdgeRho(double destRho, int direction);
  int setEdgeRhoSide(double destRho, int direction);

  virtual double getMinTheta() const;
  virtual double getMaxTheta() const;
  virtual double getMeanTheta() const;
  virtual double getMaxRho() const;
  virtual double getMinRho() const;
  virtual double getMaxZ() const;
  virtual double getMinZ() const;
  virtual const XYZVector& getMeanPoint() const;

  virtual double getMaxPhi() const { return 0.; }; 
  virtual double getMinPhi() const { return 0.; };

  virtual std::string getSensorTag() {  return std::string("");  };
  virtual std::string getSensorGeoTag() {  return std::string("");  };
  virtual std::string getPositionTag() {  return std::string("");  };

  int getRing() const {return ring_;};
  void setRing(const int& newRing) {ring_ = newRing;};

  int getSection() const { return inSection_ ;};
  void setSection(const int newSection) {inSection_ = newSection; };

  int getPhiIndex() const { return phiIndex_; }
  void setPhiIndex(int phiIndex) { phiIndex_ = phiIndex; }

  int getContainerId() const { return containerId_; }
  void setContainerId(int containerId) { containerId_ = containerId; }

  int getZSide() const { return zSide_; }
  void setZSide(int zSide) { zSide_ = zSide; }

  int getNChannels() const;
  int getNChannelsFace(int nFace) const;
  int getNMaxChannelsFace();
  int getNMinChannelsFace();

  int getNStripAcross()     { return nStripAcross_ ;};
  int getNStripsAcross()    { return nStripAcross_ ;};

  int getNSegments(int nFace);
  int getNMaxSegments();
  double getNMeanSegments();
  int getNMinSegments();

  void setNSegments(const int& face, const int& newNSegments);
  int getNFaces()           { return nFaces_ ;};
  int getReadoutType()      { return readoutType_ ;};
  int getReadoutMode()      { return readoutMode_ ;};

  void setNStripsAcross(const int& newN) { nStripAcross_=newN; }
  void setNSegments(const int& newN);
  void setNFaces(const int& newN);
  void setReadoutType(const int& newN) { readoutType_=newN; }; // TODO: check validity
  void setReadoutMode(const int& newN) { readoutMode_=newN; }; // TODO: check validity

  virtual double getLowPitch();
  virtual double getHighPitch();

  void setIrradiatedPowerConsumption(double irradiatedPowerConsumption) { irradiatedPowerConsumption_ = irradiatedPowerConsumption; }
  double getIrradiatedPowerConsumption() const { return irradiatedPowerConsumption_; }

  // R-Phi resolution 
  void setResolutionRphi(const double& newRes ) { resolutionRphi_ = newRes; };
  void setResolutionRphi(){};
  virtual double getResolutionRphi() { std::cerr << "BAD ERROR: this shouldn't happen" << std::endl; return 0; };
  virtual double getResolutionRphiTrigger()  { std::cerr << "BAD ERROR: this shouldn't happen" << std::endl; return 0; };

  // Z resolution 
  void setResolutionY(const double& py) { resolutionY_ = py; };
  void setResolutionY(){};
  double getResolutionY();
  double getResolutionYTrigger();

  //virtual double getOccupancyPerEvent();
  virtual double getStripOccupancyPerEvent();
  virtual double getHitOccupancyPerEvent();

  // Boundaries to ease the computation of tracks
  void computeBoundaries(double zError);
  bool couldHit(double eta, double phi);

  enum { Strip, Pixel, Pt,   // sensor types
	 Barrel, Endcap,     // module subdetector type
	 Rectangular, Wedge, // sensor shapes
	 Undefined };
  enum { Binary, Cluster }; // readout modes
  static const double BoundaryEtaSafetyMargin = 5. ; // track origin shift in units of zError to compute boundaries

  void setContainerName(const std::string& newName) {containerName_ = newName;};
  std::string getContainerName() {return containerName_ ;};

  virtual int getLayer() const { return 0; };
  virtual int getDisk() const { return 0;};

  virtual PosRef getPositionalReference() const { return (PosRef){ 0, 0, 0, 0 }; } 

  double getPtThreshold(const double& myEfficiency);
  double getTriggerProbability(const double& trackPt, const double& stereoDistance = 0, const int& triggerWindow=0);
  double getPtCut();

  //double getTriggerFrequencyTruePerEvent();
  double getTriggerFrequencyTruePerEventAbove(const double& myCut);
  double getTriggerFrequencyTruePerEventBelow(const double& myCut);
  double getParticleFrequencyPerEventAbove(const double& myCut);
private:
  double getTriggerFrequencyTruePerEventBetween(double myLowCut, double myHighCut);
  double getParticleFrequencyPerEventBetween(double myLowCut, double myHighCut);
public:
  double getTriggerFrequencyFakePerEvent();

  int getProcessorConnectionsEta() const { return processorConnectionsEta_; }
  void setProcessorConnectionsEta(int processorConnectionsEta) { processorConnectionsEta_ = processorConnectionsEta; }
  int getProcessorConnectionsPhi() const { return processorConnectionsPhi_; }
  void setProcessorConnectionsPhi(int processorConnectionsPhi) { processorConnectionsPhi_ = processorConnectionsPhi; }
  int getProcessorConnections() const { return processorConnectionsEta_*processorConnectionsPhi_; }

  void setProperty(std::string name, double value) { properties_[name] = value; }
  double getProperty(std::string name) const { return hasProperty(name) ? properties_.at(name) : 0; }
  bool hasProperty(std::string name) const { return properties_.find(name) != properties_.end(); }

  ptError* getPtError() { return &myPtError; }

  void lockGeometry() { geometryLocked_ = true; setPterrorParameters(); }
  void unlockGeometry() { geometryLocked_ = false; }
  bool geometryLocked() const { return geometryLocked_; }

private:
  void setPterrorParameters();

  std::map<std::string, double> properties_;

  bool geometryLocked_;
};



// Barrel Module

class BarrelModule : public Module {

 private:
  void setSensorRectGeometry(double heightOverWidth);
  edge getEdgePhi(int direction, double margin = 0);
public:
  edge getEdgeZ(int direction, double margin = 0);
private:
  double width_;
  int layer_;
  int zIndex_; // index of a module in the rod
 public:
  ~BarrelModule();
  BarrelModule(double waferDiameter, double heightOverWidth);
  BarrelModule(double heightOverWidth=1);

  double getLowPitch();
  double getHighPitch();

  double getResolutionRphi();
  double getResolutionRphiTrigger();

  virtual int getSubdetectorType() const { return Barrel; };
  edge getEdgeZSide(int direction, double margin = 0);
  int setEdgeZ(double newZ, int direction);

  edge getEdgePhiSide(int direction, double margin = 0);
  int setEdgePhi(double newPhi, int direction);

  double getMaxPhi() const; 
  double getMinPhi() const;
  double getWidth() const {return width_;};

  std::string getSensorTag();
  std::string getSensorGeoTag();
  std::string getPositionTag();

  PosRef getPositionalReference() const { return (PosRef){ containerId_, (getZSide() > 0 ? ring_ : -ring_), layer_, phiIndex_+1 }; }

  int getLayer() const {return layer_;};
  void setLayer(const int& newLayer) {layer_ = newLayer;};

 // double getOccupancyPerEvent();
  virtual double getStripOccupancyPerEvent();
  virtual double getHitOccupancyPerEvent();
};


class EndcapModule : public Module {

private:
  // Remember to update the copyconstructor when adding member variables
  //double phiWidth_;
  double widthLo_;
  double widthHi_;
  double dist_;
  bool cut_;
  double lost_; // lost millimeters in height because of cut
  void setSensorWedgeGeometry(double alpha, double d, double maxRho = -1);
  void setSensorRectGeometry(double heightOverWidth, double d);
  int disk_;  

public:
  ~EndcapModule();
  double getLowPitch();
  double getHighPitch();

  double getResolutionRphi();
  double getResolutionRphiTrigger();

  // Rectangular-shaped detectors
  // EndcapModule(double waferDiameter, double heightOverWidth); // TODO: treat the wafer's diameter properly
  EndcapModule(double heightOverWidth);
  EndcapModule(double waferDiameter, double heightOverWidth);

  // Wedge- or square-shaped (no real geometry setting here)
  //  EndcapModule(const Module& aModule);

  // Generating wedge-shaped modules
  EndcapModule(double alpha, double d, double maxRho); //use -1 as default for maxRho
  EndcapModule(const Module& sampleModule, double alpha, double d, double maxRho); //use -1 as default for maxRho
  EndcapModule(const EndcapModule& sampleModule, double alpha, double d, double maxRho); //use -1 as default for maxRho

  // Generating square-shaped modules
  EndcapModule(const Module& aModule, double d);

  EndcapModule(int shape /*=wedge*/);

  virtual int getSubdetectorType() const { return Endcap; };
  std::string getSensorTag();
  std::string getSensorGeoTag();
  std::string getPositionTag();

  bool wasCut() const { return cut_ ; };
  double getLost() const { if (cut_) return lost_; return 0;};
  double getDist() const { return dist_;};

  double getMaxPhi() const; 
  double getMinPhi() const;
  double getWidthLo() const { return widthLo_ ; };
  double getWidthHi() const { return widthHi_ ; };

  int getDisk() const { return disk_;};
  void setDisk(const int& newDisk) { disk_ = newDisk;};

  void setLayer(int newLayer) { disk_ = newLayer; }  // Layer == Disk for EndCap modules!
  int getLayer() const { return disk_; }

  PosRef getPositionalReference() const { return (PosRef){ containerId_, (getZSide() > 0 ? disk_ : -disk_), ring_, phiIndex_+1 }; }

  //double getOccupancyPerEvent();
  virtual double getStripOccupancyPerEvent();
  virtual double getHitOccupancyPerEvent();
};

#endif

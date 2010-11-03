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

using namespace ROOT::Math;

typedef std::pair<double, int> edge;

class Module {

 protected:
  int shape_;
  double aspectRatio_;

  double boundaryMinPhi_;
  double boundaryMaxPhi_;
  double boundaryMinEta_;
  double boundaryMaxEta_;
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
  double stereorot_;
  double dphideta_;
  int nChannelsPerFace_;
  int nSegments_;
  int nStripAcross_;
  int nFaces_;
  int readoutType_;
  int readoutMode_;
  double resolutionRphi_;
  double resolutionY_;

  int inSection_;

  int ring_;

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
  static const int     defaultStereoDist_ = 0;
  static const int     defaultStereoRot_ = 0;
  static const int     defaultInSection_ = 0;
  static const int     defaultChannelsPerFace_ = 1;
  static const int     defaultSegments_ = 1;
  static const int     defaultStripAcross_ = 1;
  static const int     defaultFaces_ = 1;
  static const double defaultResolutionRphi_ = 0.0;
  static const double defaultResolutionY_ = 0.0;
  
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

  void computeStripArea();
  void computeDphiDeta();

 private:
  void setDefaultParameters();
    
 public:
  virtual ~Module();
  Module();
  Module(double waferDiameter);

  virtual int getSubdetectorType() { return Undefined; };
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

  void setStereoDistance(double sdist) { stereodist_=sdist; }
  double getStereoDistance() { return stereodist_; };

  void setStereoRotation(double srot) { stereorot_=srot; }
  double getStereoRotation() { return stereorot_; };

  void setColor(const Color_t newColor) {color_ = newColor;}
  Color_t getColor() {return color_;}

  void resetNHits() { nHits_ = defaultNHits_ ; };
  int getNHits() { return nHits_ ;};

  void shapeVolume(TGeoVolume* container,
		   TGeoMedium* medium,
		   TGeoManager* geom);

  TPolyLine3D* getContour();

  double getHeight() {return height_;};
  double getArea() { return area_;};
  double getDiameter() {return waferDiameter_; };
  double getThickness() { return thickness_; };
  double getModuleThickness() { return moduleThickness_; };
  XYZVector getCorner(int index) { return corner_[index]; };

  edge getEdgeRhoSide(int direction);
  int setEdgeRho(double destRho, int direction);
  int setEdgeRhoSide(double destRho, int direction);

  const double getMinTheta() const;
  const double getMaxTheta() const;
  const double getMeanTheta() const;
  const double getMaxRho() const;
  const double getMinRho() const;
  const double getMaxZ() const;
  const double getMinZ() const;
  const XYZVector getMeanPoint() const;

  virtual std::string getSensorTag() {  return std::string("");  };

  int getRing() {return ring_;};
  void setRing(const int& newRing) {ring_ = newRing;};

  int getSection() { return inSection_ ;};
  void setSection(const int newSection) {inSection_ = newSection; };

  int getNChannels()        { return nChannelsPerFace_ * nFaces_ ;};
  int getNChannelsPerFace() { return nChannelsPerFace_ ;};
  int getNPerFace()         { return nChannelsPerFace_ ;};
  int getNStripAcross()     { return nStripAcross_ ;};
  int getNStripsAcross()    { return nStripAcross_ ;};
  int getNSegments()        { return nSegments_ ;};
  int getNFaces()           { return nFaces_ ;};
  int getReadoutType()      { return readoutType_ ;};
  int getReadoutMode()      { return readoutMode_ ;};

  void setNStripAcross(const int& newN) { nStripAcross_=newN; nChannelsPerFace_=nStripAcross_*nSegments_;  };
  void setNStripsAcross(const int& newN) { setNStripAcross(newN); }
  void setNSegments(const int& newN) { nSegments_=newN; nChannelsPerFace_=nStripAcross_*nSegments_;  };
  void setNFaces(const int& newN) { nFaces_=newN; };
  void setReadoutType(const int& newN) { readoutType_=newN; }; // TODO: check validity
  void setReadoutMode(const int& newN) { readoutMode_=newN; }; // TODO: check validity

  virtual double getLowPitch();
  virtual double getHighPitch();

  // R-Phi resolution 
  void setResolutionRphi(const double& newRes ) { resolutionRphi_ = newRes; };
  void setResolutionRphi(){};
  virtual double getResolutionRphi() { std::cerr << "BAD ERROR: this shouldn't happen" << std::endl; return 0; };

  // Z resolution 
  void setResolutionY(const double& py) { resolutionY_ = py; };
  void setResolutionY(){};
  double getResolutionY();

  virtual double getOccupancyPerEvent();

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


};



// Barrel Module

class BarrelModule : public Module {

 private:
  void setSensorRectGeometry(double heightOverWidth);
  edge getEdgeZ(int direction, double margin = 0);
  edge getEdgePhi(int direction, double margin = 0);
  double width_;
  int layer_;
  
 public:
  ~BarrelModule();
  BarrelModule(double waferDiameter, double heightOverWidth);
  BarrelModule(double heightOverWidth=1);

  double getLowPitch();
  double getHighPitch();

  double getResolutionRphi();

  virtual int getSubdetectorType() { return Barrel; };
  edge getEdgeZSide(int direction, double margin = 0);
  int setEdgeZ(double newZ, int direction);

  edge getEdgePhiSide(int direction, double margin = 0);
  int setEdgePhi(double newPhi, int direction);

  double getMaxPhi();
  double getMinPhi();
  double getWidth() {return width_;};

  std::string getSensorTag();

  int getLayer() const {return layer_;};
  void setLayer(const int& newLayer) {layer_ = newLayer;};

  double getOccupancyPerEvent();
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

  virtual int getSubdetectorType() { return Endcap; };

  // The constructor needs 
  // alpha: the angle covered by the module
  // d: the dstance of the module's base from the z axis
  EndcapModule(int shape=Wedge);

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

  std::string getSensorTag();

  bool wasCut() const { return cut_ ; };
  double getLost() const { if (cut_) return lost_; return 0;};
  double getDist() const { return dist_;};

  double getWidthLo() const { return widthLo_ ; };
  double getWidthHi() const { return widthHi_ ; };

  int getDisk() const { return disk_;};
  void setDisk(const int& newDisk) {disk_ = newDisk;};

  double getOccupancyPerEvent();
};

#endif

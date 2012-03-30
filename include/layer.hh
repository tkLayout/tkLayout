#ifndef _LAYER_HH_
#define _LAYER_HH_

// Standard stuff
#include <vector>
#include <string>

// ROOT stuff
#include <Math/Vector3D.h>
#include <TGeoManager.h>

// Our stuff
#include <module.hh>
#include <messageLogger.h>

using namespace ROOT::Math;

typedef std::vector<Module* > ModuleVector;
typedef std::pair<int,double> LayerOption;

class Layer : public MessageLogger {

protected:
  ModuleVector moduleSet_;
  std::string layerName_;
  std::string containerName_;
  int layerIndex_;
  
private:
  virtual void setDefaultParameters();
  
public:
  virtual ~Layer();
  Layer();
  void translate(XYZVector Delta);
  void rotatePhi(double phi) {/*TODO*/};
  void shapeVolume(TGeoVolume* container, TGeoMedium* medium, TGeoManager* geom);
  void shapeModuleVolumes(TGeoVolume* container, TGeoMedium* medium, TGeoManager* geom);
  ModuleVector* getModuleVector() { return &moduleSet_; }
  virtual Module* getSampleModule() { return NULL; }

  std::string getName() {return layerName_; };
  int getIndex() {return layerIndex_; };
  std::string getContainerName() {return containerName_; };
  void setName(const std::string& newName, const int& newIndex) { layerName_ = newName; layerIndex_ = newIndex; setObjectName(layerName_); };
  void setContainerName(const std::string& newName ) { containerName_ = newName; };
  
  double getMaxZ();
  double getMinZ();
  double getMaxRho();
  double getMinRho();
  virtual double getTilt() { return 0.0; }
  virtual double getStartAngle() { return 0.0; }
  virtual double getMaxModuleThickness() { return 0.0; }

  virtual int cutOverEta(double etaCut) { return 0; }

  virtual void decreaseModCount(int ring) {}

  enum {NoSection = 0x0,
	XYSection = 0x1,
	YZSection = 0x2,
	Forward   = 0x4};

  // Directives
  enum {SHRINK  = -1,
	FIXED   = -2,
	ENLARGE = -3,
	AUTO    = -4};

  // Options
  enum {Stacked = 1};

  // Other constants
  enum {InvalidRadius = -1 }; // It needs to be negative

protected:
  int round(const double& x, const bool& odd);

};


class BarrelLayer : public Layer {
private:
  int nOfRods_, nModsOnString_;
  BarrelModule* sampleModule_;
  BarrelModule* getSampleModule() { return sampleModule_; }
  double averageRadius_;
  void setDefaultParameters(){ averageRadius_=0;};

  std::pair<double, int> computeRadius(const double& x,    // radius of the inner module
				       const double& g,    // Gap between inner and outer
				       const double& o,    // Needed overlap in mm
				       const double& b,    // Module's half width
				       const int& optimal, // wether to shrink or enlarge, fix or auto
				       const int& base );

  int buildString(ModuleVector& thisModuleSet,
		  double stringAverageRadius,
		  double smallDelta, // Half the distance between inner and outer modules
		  double zOverlap,
		  double safetyOrigin,
		  double maxZ,
		  BarrelModule* sampleModule);
  int buildString(ModuleVector& thisModuleSet,
		  double stringAverageRadius,
		  double smallDelta, // Half the distance between inner and outer modules
		  double zOverlap,
		  double safetyOrigin,
		  int nModules,
		  BarrelModule* sampleModule,
		  double minZ = 0);
  
  int buildMezzanineString(ModuleVector& thisModuleSet,
		  double stringAverageRadius,
		  double smallDelta, // Half the distance between inner and outer modules
		  double zOverlap,
		  double safetyOrigin,
		  int nModules,
		  BarrelModule* sampleModule,
		  double farthestZ = 0);

double computeListZ(  // New-fangled kickass function
    std::vector<double>& listZ,
    double startZ,
    double radiiRatio,  // r_max / r_min
    double modLengthZ,
    double originDeltaZ,
    double baseOverlapZ,
    int numModules,
    int parity,
    int direction,
    bool looseStartZ = false); 

  void buildStringPair(  // New-fangled kickass function
          ModuleVector& thisModuleSet,
          double averageRadius,
          double smallDelta,
          double baseOverlap,
          double zDelta,
          double startZ,
          int numModules,
          int smallParity,
          BarrelModule* sampleModule);
  void buildMezzanineStringPair(  // New-fangled kickass function
          ModuleVector& thisModuleSet,
          double averageRadius,
          double smallDelta,
          double baseOverlap,
          double zDelta,
          double startZ,
          int numModules,
          int smallParity,
          BarrelModule* sampleModule);
/*  void buildStringPair(ModuleVector& thisModuleSet, // CUIDADO: TBR
		       double stringAverageRadius,
		       double smallDelta, // Half the distance between inner and outer modules
		       double zOverlap,
		       double safetyOrigin,
		       double maxZ,
		       BarrelModule* sampleModule); */
  void buildStringPair(ModuleVector& thisModuleSet, // CUIDADO: maybe TBR
		       double stringAverageRadius,
		       double smallDelta, // Half the distance between inner and outer modules
		       double zOverlap,
		       double safetyOrigin,
		       int nModules,
		       BarrelModule* sampleModule);
  double layerRadius(const double& nMod,  // number of strings in a layer
		     const double& g,     // Gap between inner and outer
		     const double& o,     // Needed overlap in mm
		     const double& b);    // Module's half width
  // This function was tuned "by hand" (see code for comments)
  std::pair<double, int> layerPhi(double, double, double, double, double, int, int, bool);

  
public:
  ~BarrelLayer();
  BarrelLayer();
  BarrelLayer(BarrelLayer& sampleLayer);
  BarrelLayer(double waferDiameter, double heightOverWidth);
  BarrelLayer(double heightOverWidth);
  BarrelLayer(const BarrelModule& mySample);
  BarrelLayer(BarrelModule* mySample);
  
  // An optimization function to run prior of settling the geometry
  void neededModulesPlot(double smallDelta, // Half distance between modules in the same string
			 double bigDelta,   // String half gap
			 double o,          // overlap
			 double b,          // Module half width
			 int base);


  void buildLayer (double averageRadius,
		   double smallDelta, 
		   double bigDelta, 
		   double overlap, 
		   double safetyOrigin, 
		   double maxZ, 
		   int pushDirection, 
		   int base,
		   bool stringSameParity,
		   BarrelModule* sampleModule,
		   int sectioned=NoSection) { std::cout << "DEPRECATED" << std::endl; } ; // TODO: remove this

  void buildLayer (double averageRadius,
		   double smallDelta, 
		   double bigDelta, 
		   double overlap, 
		   double safetyOrigin, 
		   int nModules, 
		   int pushDirection, 
		   int base,
		   bool stringSameParity,
		   BarrelModule* sampleModule,
		   int sectioned = NoSection,
		   double minZ = 0.);

  int cutOverEta(double etaCut);

  int getRods() { return nOfRods_; }
  int getModulesOnRod() { return nModsOnString_; }

  double getMaxZ(int direction);
  void compressToZ(double newMaxZ);
  void compressExceeding(double newMaxZ, double newMinZ);

  double computeAverageRadius();
  double getAverageRadius() { if (averageRadius_<=0) return computeAverageRadius(); return averageRadius_;};
  void rotateY_PI();
  void reflectZ();
  void shiftRho(double Delta);

  double getMaxModuleThickness() { if (sampleModule_ != NULL) return sampleModule_->getModuleThickness(); return 0.0; }
};


class EndcapLayer : public Layer {
private:
  int nOfRings_;
  std::vector<int> nModsOnRing_;
  EndcapModule* sampleModule_;
  EndcapModule* getSampleModule() { return sampleModule_; }
  double averageZ_;
  void setDefaultParameters(){averageZ_=0;};

  double solvex(double y);
  double gamma1(double x, double y, double r);
  double gamma2(double x, double y, double r);
  double Area(double x, double y, double r);
  double compute_l(double x, double y, double d);
  double compute_d(double x, double y, double l);

  
public:
  ~EndcapLayer();
  EndcapLayer();
  EndcapLayer(EndcapLayer& sampleLayer);
  EndcapLayer(const Module& mySample, double alpha, double d);
  EndcapLayer(double alpha, double d);
  EndcapLayer(const EndcapModule& mySample);
  EndcapLayer(const EndcapModule& mySample, double alpha, double d);

  void translateZ(const double& zShift);

  void buildSingleDisk(double minRadius,
		       double maxRadius,
		       double smallDelta, 
		       double bigDelta,
		       double diskZ, 
		       double overlap, 
		       double zError,
		       int phiSegments, 
		       bool oddSegments, bool alignEdges,
		       std::map<int, EndcapModule*> sampleModule, 
		       std::map<int, int> ringDirectives, 
		       int diskParity,
		       int sectioned = NoSection);

  double buildRing(double minRadius,
		   double smallDelta, 
		   double bigDelta, 
		   double diskZ, 
		   double overlap, 
		   int phiSegments,
		   bool oddSegments, bool alignEdges,
		   int nearDirection, 
		   EndcapModule* sampleModule,
		   double maxRadius = -1,
		   int addModules = 0,
		   int sectioned = NoSection);

  int cutOverEta(double etaCut);

  int getRings() { return nOfRings_; }
  std::vector<int>& getModulesOnRing() { return nModsOnRing_; }
  virtual void decreaseModCount(int ring) { nModsOnRing_.at(ring)--; }

  void rotateY_PI();
  double getAverageZ() {return averageZ_;};
  
  double getMaxModuleThickness();

};


#endif







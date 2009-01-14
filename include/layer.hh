#ifndef _LAYER_HH_
#define _LAYER_HH_

#include <vector>
#include "Math/Vector3D.h"
#include <string>
#include "TGeoManager.h"
#include "module.hh"

using namespace ROOT::Math;

typedef std::vector<Module* > ModuleVector;

class Layer {
protected:
  ModuleVector moduleSet_;
  std::string layerName_;
  
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

  std::string getName() {return layerName_; };
  void setName(const std::string newName ) { layerName_ = newName; };
  
  double getMaxZ();
  double getMinZ();
  double getMaxRho();
  double getMinRho();

  int cutOverEta(double etaCut);

  enum {NoSection = 0x0,
	XYSection = 0x1,
	YZSection = 0x2,
	Forward   = 0x4};

  enum {SHRINK  = -1,
	FIXED   = -2,
	ENLARGE = -3,
	AUTO    = -4};
};


class BarrelLayer : public Layer {
private:
  BarrelModule* sampleModule_;
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
  
  void buildStringPair(ModuleVector& thisModuleSet,
		       double stringAverageRadius,
		       double smallDelta, // Half the distance between inner and outer modules
		       double zOverlap,
		       double safetyOrigin,
		       double maxZ,
		       BarrelModule* sampleModule);
  void buildStringPair(ModuleVector& thisModuleSet,
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

  double getMaxZ(int direction);
  void compressToZ(double newMaxZ);
  void compressExceeding(double newMaxZ, double newMinZ);

  double getAverageRadius() {return averageRadius_;};
  void rotateY_PI();
  void reflectZ();

};


class EndcapLayer : public Layer {
private:
  EndcapModule* sampleModule_;
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
		       int base, 
		       EndcapModule* sampleModule, 
		       std::map<int, int> ringDirectives, 
		       int diskParity = -1,
		       int sectioned = NoSection);

  double buildRing(double minRadius,
		   double smallDelta, 
		   double bigDelta, 
		   double diskZ, 
		   double overlap, 
		   int base,
		   int nearDirection, 
		   EndcapModule* sampleModule,
		   double maxRadius = -1,
		   int addModules = 0,
		   int sectioned = NoSection);


  void rotateY_PI();
  double getAverageZ() {return averageZ_;};
  
};


#endif







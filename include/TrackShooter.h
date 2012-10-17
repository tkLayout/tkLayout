#ifndef TRACK_SHOOTER_H
#define TRACK_SHOOTER_H

#include <ostream>
#include <map>
#include <functional>
#include <list>
#include <vector>

#include <TRandom3.h>
#include <TCanvas.h>
#include <TImage.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TView.h>
#include <TSystem.h> // CUIDADO: WTF is this?
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <Math/GenVector/RotationZ.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/Transform3D.h>

#include <global_funcs.h>
#include <global_constants.h>
#include <module.hh>
#include <PlotDrawer.h>
#include <Palette.h>


template<typename T> int getPointOctant(const T& x, const T& y, const T& z) {
  return (int(x < T(0)) << 2) || (int(x < T(0)) << 1) || int(x < T(0));
}


std::set<int> getModuleOctants(const Module* mod);


template<int NumSides, class Coords, class Random>
class AbstractPolygon {  // any number of dimensions, any number of sides, convex or concave. it only has 2 properties
protected:
  double area_;
  std::vector<Coords> v_; // vertices of the polygon
  virtual void computeProperties() = 0;
public:
  AbstractPolygon() { v_.reserve(NumSides); }
  virtual AbstractPolygon<NumSides, Coords, Random>& operator<<(const Coords& vertex);
  virtual AbstractPolygon<NumSides, Coords, Random>& operator<<(const std::vector<Coords>& vertices);
  const std::vector<Coords>& getVertices() const;
  const Coords& getVertex(int index) const;
  int getNumSides() const;
  bool isComplete() const;
  double getDoubleArea() const; // we save a division by 2 by returning the double area. for our purposes it's the same
  virtual bool isPointInside(const Coords& p) const = 0;
  virtual Coords generateRandomPoint(Random* die) const = 0;
};

template<class Polygon>
struct PolygonLess {
  bool operator()(const Polygon& t1, const Polygon& t2) const {
    return t1.getDoubleArea() < t2.getDoubleArea();
  }
};

template<int NumSides>
class Polygon3d : public AbstractPolygon<NumSides, XYZVector, TRandom> { // no checks are made on the convexity, but the algorithms in the class only work for convex polygons, so beware!
public:
  typedef std::multiset<Polygon3d<3>, PolygonLess<Polygon3d<3> > > TriangleSet;
protected:
  TriangleSet trianglesByArea_;
  void computeProperties();
public:
  const TriangleSet& getTriangulation() const;
  bool isPointInside(const XYZVector& p) const;
  XYZVector generateRandomPoint(TRandom* die) const;
};

template<>
class Polygon3d<3> : public AbstractPolygon<3, XYZVector, TRandom> {  // a triangle can be defined with more than 3 vertices. no error checking is made, simply the additional vertices besides the 3rd one are ignored
  void computeProperties();
public:
  bool isPointInside(const XYZVector& p) const;
  XYZVector generateRandomPoint(TRandom* die) const;
};

typedef Polygon3d<3> Triangle3d;




class ParticleGenerator {
public:
  static const int MAX_ENTRIES = 12000;
  static const int Z0_SMEAR_MM = 70;

  TRandom& die;

  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent; //!current Tree number in a TChain

  Long64_t numEntries;

  // Declaration of leaf types
  Float_t         qScale;
  Float_t         pThat;
  Int_t           NumPart;
  Int_t           pdgId[MAX_ENTRIES];   //[NumPart]
  Int_t           status[MAX_ENTRIES];   //[NumPart]
  Int_t           charge[MAX_ENTRIES];   //[NumPart]
  Float_t         Eta[MAX_ENTRIES];   //[NumPart]
  Float_t         Phi[MAX_ENTRIES];   //[NumPart]
  Float_t         Pt[MAX_ENTRIES];   //[NumPart]
  Float_t         E[MAX_ENTRIES];   //[NumPart]

  // List of branches
  TBranch        *b_qScale;   //!
  TBranch        *b_pThat;   //!
  TBranch        *b_NumPart;   //!
  TBranch        *b_pdgId;   //!
  TBranch        *b_status;   //!
  TBranch        *b_charge;   //!
  TBranch        *b_Eta;   //!
  TBranch        *b_Phi;   //!
  TBranch        *b_Pt;   //!
  TBranch        *b_E;   //!

  struct Particle { 
    double pt; 
    double eta; 
    double phi; 
    double z0;
  };


  ParticleGenerator(TRandom& aDie);
  virtual ~ParticleGenerator();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual Bool_t   Notify();

  Particle getParticle();
  Long64_t getNumEntries() const;
};





class TrackShooter {
  static const float HIGH_PT_THRESHOLD = 2.;

  std::ostream* output_;
  const char *FS, *LS;

  std::vector<Module*> allMods_;
  typedef std::list<BarrelModule*> BarrelModules;   // accessed sequentially. last step in hit localization
  typedef std::list<EndcapModule*> EndcapModules;
  typedef std::vector<BarrelModules> BarrelOctants; // these vectors are accessed randomly when the octant of the hit is known (function of the radius/z discovered in the previous step)
  typedef std::vector<EndcapModules> EndcapOctants;
  typedef std::map<double, BarrelOctants> BarrelRadii; // first step in hit localization. we fix radius/z and we calculate the rest of the coordinates of a hit
  typedef std::map<double, EndcapOctants> EndcapZs;
  BarrelRadii barrelModsByRadius_;
  EndcapZs endcapModsByZ_;

public:
  void setOutput(std::ostream& output, const char* fieldSeparator = "\t", const char* lineSeparator = "\n", bool synced = false);
  void addModule(Module* module);
  void shootTracks(long int numEvents, long int numTracksPerEvent);

  void manualPolygonTestBench();
  void moduleTestBench();
};


#endif

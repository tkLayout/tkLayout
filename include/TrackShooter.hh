#ifndef TRACK_SHOOTER_H
#define TRACK_SHOOTER_H

#include <ostream>
#include <map>
#include <functional>
#include <algorithm>
#include <list>
#include <vector>
#include <memory>


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
#include <Math/GenVector/DisplacementVector2D.h>
#include <Math/GenVector/RotationZ.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/Transform3D.h>
#include <boost/program_options/variables_map.hpp>

#include <global_funcs.hh>
#include <global_constants.hh>
#include <module.hh>
#include <PlotDrawer.hh>
#include <Palette.hh>

namespace po = boost::program_options;

template<typename T> int getPointOctant(const T& x, const T& y, const T& z) {
  return (int(x < T(0)) << 2) || (int(x < T(0)) << 1) || int(x < T(0));
}

typedef ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<double>,ROOT::Math::DefaultCoordinateSystemTag> XYVector; // CUIDADO The version of ROOT tkLayout is linked with misses this typedef

std::set<int> getModuleOctants(const Module* mod);


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

struct ModuleData {
  double x, y, z;
  double rho, phi;
  double widthlo, widthhi, height;
  double stereo; 
  double pitchlo, pitchhi;
  double striplen;
  double yres;
  char inefftype;
  char refcnt, refz, refrho, refphi; // positional reference
  char type;
};

struct Tracks { // holder struct for TTree export
  std::vector<unsigned> eventn;
  std::vector<unsigned> trackn;
  std::vector<double> eta;
  std::vector<double> phi0;
  std::vector<double> z0;
  std::vector<double> pt;
  std::vector<char> nhits;
  const std::string name;
  Tracks(const std::string& name_) : name(name_) {}

  void clear() {
    eventn.clear();
    trackn.clear();
    eta.clear();
    phi0.clear();
    z0.clear();
    pt.clear();
    nhits.clear();
  }
  void push_back(int eventn_, int trackn_, double eta_, double phi0_, double z0_, double pt_, int8_t nhits_) {
    eventn.push_back(eventn_);
    trackn.push_back(trackn_);
    eta.push_back(eta_);
    phi0.push_back(phi0_);
    z0.push_back(z0_);
    pt.push_back(pt_);
    nhits.push_back(nhits_);
  }
  void setupBranches(TTree& tree) {
    tree.Branch((name + ".eventn").c_str(), &eventn);
    tree.Branch((name + ".trackn").c_str(), &trackn);
    tree.Branch((name + ".eta").c_str(), &eta);
    tree.Branch((name + ".phi0").c_str(), &phi0);
    tree.Branch((name + ".z0").c_str(), &z0);
    tree.Branch((name + ".pt").c_str(), &pt);
    tree.Branch((name + ".nhits").c_str(), &nhits);
  }
};

struct Hits { // holder struct for TTree export
  std::vector<double> glox, gloy, gloz;
  std::vector<double> locx, locy;
  std::vector<float> pterr, hitprob;
  std::vector<float> deltas;
  std::vector<char> cnt, z, rho, phi;
//  std::vector<double> distx, disty;
//  std::vector<double> eta;
  const std::string name;
  Hits(const std::string& name_) : name(name_) {}

  void clear() {
    glox.clear(); gloy.clear(); gloz.clear();
    locx.clear(); locy.clear();
    pterr.clear(); hitprob.clear();
    deltas.clear();
    cnt.clear(); z.clear(); rho.clear(); phi.clear();
//    distx.clear(); disty.clear();
//    eta.clear();
  }
  void push_back(double glox_, double gloy_, double gloz_, double locx_, double locy_, float pterr_, float hitprob_, float deltas_, int8_t cnt_, int8_t z_, int8_t rho_, int8_t phi_/*, double distx_, double disty_, double eta_*/) {
    glox.push_back(glox_);
    gloy.push_back(gloy_);
    gloz.push_back(gloz_);
    locx.push_back(locx_);
    locy.push_back(locy_);
    pterr.push_back(pterr_);
    hitprob.push_back(hitprob_);
    deltas.push_back(deltas_);
    cnt.push_back(cnt_);
    z.push_back(z_);
    rho.push_back(rho_);
    phi.push_back(phi_);
//    distx.push_back(distx_);
//    disty.push_back(disty_);
//    eta.push_back(eta_);
  }
  int size() const { return glox.size(); }

  void setupBranches(TTree& tree) {
    tree.Branch((name + ".glox").c_str(), &glox);
    tree.Branch((name + ".gloy").c_str(), &gloy);
    tree.Branch((name + ".gloz").c_str(), &gloz);
    tree.Branch((name + ".locx").c_str(), &locx);
    tree.Branch((name + ".locy").c_str(), &locy);
    tree.Branch((name + ".pterr").c_str(), &pterr);
    tree.Branch((name + ".hitprob").c_str(), &hitprob);
    tree.Branch((name + ".deltas").c_str(), &deltas);
    tree.Branch((name + ".cnt").c_str(), &cnt);
    tree.Branch((name + ".z").c_str(), &z);
    tree.Branch((name + ".rho").c_str(), &rho);
    tree.Branch((name + ".phi").c_str(), &phi);
//    tree.Branch((name + ".distx").c_str(), &distx);
//    tree.Branch((name + ".disty").c_str(), &disty);
//    tree.Branch((name + ".eta").c_str(), &eta);
  }
};



struct Helix {
  double h, k;
  double R, Rs; // Rs = Radius with sign depending on particle direction
  double pt;
  double eta, theta, L;
  double phi0;
  double z0;
  double d0;
  int dir;

  Helix(double pt_, double eta_, double phi0_, double z0_, double d0_)
    : pt(pt_), eta(eta_), phi0(phi0_), z0(z0_), d0(d0_) {
      theta = 2*atan(exp(-eta));
      R = fabs(pt)/(0.3*insur::magnetic_field) * 1e3;
      Rs = R*dir + d0;
      dir = signum(pt);
      h = -Rs*sin(phi0);
      k =  Rs*cos(phi0);
      L = tan(M_PI/2 - theta);
  }
};



template<class T>
class Value {
public:
  virtual T get() = 0;
  virtual std::string toString() const = 0;
};

template<class T>
class ConstValue : public Value<T> {
  T value_;
public:
  ConstValue(T value) : value_(value) {}
  T get() { return value_; }
  std::string toString() const { return any2str(value_); } 
};

template<class T>
class UniformValue : public Value<T> {
  TRandom& die_;
  T min_, max_;
public:
  UniformValue(TRandom& die, T min, T max) : die_(die), min_(min), max_(max) {}
  T get() { return die_.Uniform(min_, max_); }
  std::string toString() const { return any2str(min_) + ":" + any2str(max_); }
};

template<class T>
class BinaryValue : public Value<T> {
  TRandom& die_;
  T value0_, value1_;
public:
  BinaryValue(TRandom& die, T value0, T value1) : die_(die), value0_(value0), value1_(value1) {}
  T get() { return die_.Integer(2) ? value1_ : value0_; }
  std::string toString() const { return any2str(value0_) + "," + any2str(value1_); }
};

/*
template<>
class BinaryValue<bool> : public Value<bool> {
  TRandom& die_;
  BinaryValue(TRandom& die) : die_(die) {}
  bool get() { return die_.Integer(2); }
  std::string toString() const { return "true|false"; }
};
*/

template<class T>
std::auto_ptr<Value<T> > valueFromString(TRandom& die, const std::string& str) {
  std::vector<std::string> values = split(str, ":,");
  if (str.find(":") != std::string::npos && values.size() == 2) return std::auto_ptr<Value<T> >(new UniformValue<T>(die, str2any<T>(values[0]), str2any<T>(values[1])));
  else if (str.find(",") != std::string::npos && values.size() == 2) return std::auto_ptr<Value<T> >(new BinaryValue<T>(die, str2any<T>(values[0]), str2any<T>(values[1])));
  else return std::auto_ptr<Value<T> >(new ConstValue<T>(str2any<T>(values[0])));
}



class TrackShooter {
  static const float HIGH_PT_THRESHOLD = 2.;

  double trackerMaxRho_;
  double barrelMinZ_, barrelMaxZ_;

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

  TRandom3 die_;
  std::auto_ptr<Value<double> > eta_, phi0_, z0_, pt_, invPt_;  // auto_ptr because I don't want to be bothered with deletion
  std::auto_ptr<Value<int> > charge_;
  bool useInvPt_;

  long int numEvents_, numTracksEv_, eventOffset_;
  std::string instanceId_;
  std::string tracksDir_;

  int detectCollisionSlanted(const Helix& helix, const Polygon3d<4>& poly, std::vector<XYZVector>& collisions);
  int detectCollisionBarrel(const Helix& helix, const Polygon3d<4>& poly, std::vector<XYZVector>& collisions);
  int detectCollisionEndcap(const Helix& helix, const Polygon3d<4>& poly, std::vector<XYZVector>& collisions);

  XYVector convertToLocalCoords(const XYZVector& globalHit, const BarrelModule* mod) const;
  XYVector convertToLocalCoords(const XYZVector& globalHit, const EndcapModule* mod) const;

  void shootTracks();

  void setDefaultParameters();

  void printParameters();
public:
  TrackShooter() : trackerMaxRho_(std::numeric_limits<double>::max()), barrelMinZ_(0.), barrelMaxZ_(0.) { setDefaultParameters(); }
  void setOutput(std::ostream& output, const char* fieldSeparator = "\t", const char* lineSeparator = "\n", bool synced = false);
  void setTrackerBoundaries(double trackerMaxRho, double barrelMinZ, double barrelMaxZ);  // if not set the tracker is considered infinite in the rho direction and no particle can escape without ever curving back
  void addModule(Module* module);
  void shootTracks(long int numEvents, long int numTracksPerEvent, int seed);
  void shootTracks(const po::variables_map& varmap, int seed);
  void exportGeometryData();
  //void manualPolygonTestBench();
  //void moduleTestBench();
};

#endif

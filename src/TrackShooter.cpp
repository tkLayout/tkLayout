#include <TrackShooter.h>


std::set<int> getModuleOctants(const Module* mod) {
  std::set<int> octants;
  octants.insert(getPointOctant(mod->getCorner(0).X(), mod->getCorner(0).Y(), mod->getCorner(0).Z())); 
  octants.insert(getPointOctant(mod->getCorner(1).X(), mod->getCorner(1).Y(), mod->getCorner(1).Z())); 
  octants.insert(getPointOctant(mod->getCorner(2).X(), mod->getCorner(2).Y(), mod->getCorner(2).Z())); 
  octants.insert(getPointOctant(mod->getCorner(3).X(), mod->getCorner(3).Y(), mod->getCorner(3).Z())); 
  return octants;
}




template<int NumSides, class Coords, class Random, class FloatType> 
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::operator<<(const Coords& vertex) {
  if (!isComplete()) v_.push_back(vertex);
  if (isComplete()) computeProperties();
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> 
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::operator<<(const std::vector<Coords>& vertices) {
  for(typename std::vector<Coords>::const_iterator it = vertices.begin(); (v_.size() < NumSides) && (it != vertices.end()); ++it) {
    *this << *it;
  }
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType> 
inline bool AbstractPolygon<NumSides, Coords, Random, FloatType>::isComplete() const { 
  return v_.size() == NumSides; 
}

template<int NumSides, class Coords, class Random, class FloatType> 
inline const std::vector<Coords>& AbstractPolygon<NumSides, Coords, Random, FloatType>::getVertices() const { 
  return v_; 
}

template<int NumSides, class Coords, class Random, class FloatType> 
inline const Coords& AbstractPolygon<NumSides, Coords, Random, FloatType>::getVertex(int index) const { 
  return v_[index]; 
}

template<int NumSides, class Coords, class Random, class FloatType> 
double AbstractPolygon<NumSides, Coords, Random, FloatType>::getDoubleArea() const {
  return area_;
}


template<int NumSides, class Coords, class Random, class FloatType> 
const Coords& AbstractPolygon<NumSides, Coords, Random, FloatType>::getCenter() const {
  if (dirty_) {
    for (int i=0; i<NumSides; i++)
      center_ += v_[i];
    center_ /= NumSides;
    dirty_ = false;
  }
  return center_;
}

template<int NumSides, class Coords, class Random, class FloatType>
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::translate(const Coords& vector) {
  std::transform(v_.begin(), v_.end(), std::bind2nd(std::plus<Coords>(), vector));
  dirty_ = true;
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType>
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::rotateX(FloatType angle) {
  for (typename std::vector<Coords>::iterator it = v_.begin(); it != v_.end(); ++it) {
    *it = Coords(it->X(),
                 it->Y()*cos(angle) - it->Z()*sin(angle),
                 it->Y()*sin(angle) + it->Z()*cos(angle));
  }
  dirty_ = true;
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType>
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::rotateY(FloatType angle) {
  for (typename std::vector<Coords>::iterator it = v_.begin(); it != v_.end(); ++it) {
    *it = Coords(it->Z()*sin(angle) + it->X()*cos(angle),
                 it->Y(),
                 it->Z()*cos(angle) - it->X()*sin(angle));
  }
  dirty_ = true;
  return *this;
}

template<int NumSides, class Coords, class Random, class FloatType>
AbstractPolygon<NumSides, Coords, Random, FloatType>& AbstractPolygon<NumSides, Coords, Random, FloatType>::rotateZ(FloatType angle) {
  for (typename std::vector<Coords>::iterator it = v_.begin(); it != v_.end(); ++it) {
    *it = Coords(it->X()*cos(angle) - it->Y()*sin(angle),
                 it->X()*sin(angle) + it->Y()*cos(angle),
                 it->Z());
  }
  dirty_ = true;
  return *this;
}



template<int NumSides>
bool Polygon3d<NumSides>::isPointInside(const XYZVector& p) const {
  double sum = 0;
  for (typename std::vector<XYZVector>::const_iterator it = this->v_.begin(); it != this->v_.end(); ++it) {
    Triangle3d t;
    t << p << *it << (it+1 < this->v_.end() ? *(it+1) : *this->v_.begin());
    sum += t.getDoubleArea();
    if (sum - this->getDoubleArea() > 1e-4) return false;  // early quit if sum area is already bigger
  }
  return fabs(this->getDoubleArea() - sum) < 1e-4;
}



template<int NumSides>
void Polygon3d<NumSides>::computeProperties() {
  this->area_ = 0;
  for (typename std::vector<XYZVector>::const_iterator it = this->v_.begin()+1; it != this->v_.end()-1; ++it) {
    Triangle3d t;
    t << *this->v_.begin() << *it << *(it+1);
    double triangleArea = t.getDoubleArea();
    trianglesByArea_.insert(t);
    this->area_ += triangleArea;
  }
}

template<int NumSides>
inline const std::multiset<Triangle3d, PolygonLess<Triangle3d> >& Polygon3d<NumSides>::getTriangulation() const {
  return trianglesByArea_;
}

template<int NumSides> 
XYZVector Polygon3d<NumSides>::generateRandomPoint(TRandom* die) const {
  double binPicker = die->Uniform(this->getDoubleArea());
  double prevBin = 0;
  TriangleSet::const_iterator it;
  for (it = getTriangulation().begin(); it != getTriangulation().end(); it++) {
    if (it->getDoubleArea() + prevBin >= binPicker) break;  // found triangle
    prevBin += it->getDoubleArea();
  }
  return it->generateRandomPoint(die);
}

void Triangle3d::computeProperties() {
   this->area_ = sqrt((this->v_[1] - this->v_[0]).Cross(this->v_[2] - this->v_[0]).Mag2());
}

XYZVector Triangle3d::generateRandomPoint(TRandom* die) const {
  double a = die->Rndm(), b = die->Rndm();
  if (a + b > 1) { a = 1 - a; b = 1 - b; }
  return this->v_[0] + a*(this->v_[1]-this->v_[0]) + b*(this->v_[2]-this->v_[0]); // seriously C++??? you've been around for a while now, time to fix this horrid syntax??
}

bool Triangle3d::isPointInside(const XYZVector& p) const { // note: this is exactly the same as the non specialized case
  double sum = 0;
  for (std::vector<XYZVector>::const_iterator it = this->v_.begin(); it != this->v_.end(); ++it) {
    Triangle3d t;
    t << p << *it << (it+1 < this->v_.end() ? *(it+1) : *this->v_.begin());
    sum += t.getDoubleArea();
    if (sum - this->getDoubleArea() > 1e-4) return false;  // early quit if sum area is already bigger
  }
  return fabs(this->getDoubleArea() - sum) < 1e-4;
}




ParticleGenerator::ParticleGenerator(TRandom& aDie) : die(aDie) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("BunchX_PhaseIISLHC_rootuple.root");
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MinBias12k_ppAt14TeV_rootuple.root");
  if (!f) {
    //f = new TFile("BunchX_PhaseIISLHC_rootuple.root");
    f = new TFile("MinBias12k_ppAt14TeV_rootuple.root");
  }
  TTree* tree = (TTree*)gDirectory->Get("particles");

  //myPtz=NULL;
  Init(tree);
}


ParticleGenerator::~ParticleGenerator() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}


Int_t ParticleGenerator::GetEntry(Long64_t entry) {
   // Read contents of entry.
  if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t ParticleGenerator::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void ParticleGenerator::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
  fChain->SetBranchAddress("pThat", &pThat, &b_pThat);
  fChain->SetBranchAddress("NumPart", &NumPart, &b_NumPart);
  fChain->SetBranchAddress("pdgId", pdgId, &b_pdgId);
  fChain->SetBranchAddress("status", status, &b_status);
  fChain->SetBranchAddress("charge", charge, &b_charge);
  fChain->SetBranchAddress("Eta", Eta, &b_Eta);
  fChain->SetBranchAddress("Phi", Phi, &b_Phi);
  fChain->SetBranchAddress("Pt", Pt, &b_Pt);
  fChain->SetBranchAddress("E", E, &b_E);
  Notify();

  numEntries = fChain->GetEntriesFast();
}

Bool_t ParticleGenerator::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}


ParticleGenerator::Particle ParticleGenerator::getParticle() {
    if (fChain == 0 || getNumEntries() == 0) return (Particle){ 0, 0, 0, 0 };
    Long64_t entryIndex;
    Int_t partIndex; 
    do {
      entryIndex = die.Integer(getNumEntries()); // event entry index
      fChain->GetEntry(entryIndex); // retrieves the entry, appropriately setting the class variables
      partIndex = die.Integer(NumPart); // particle index inside the event entry
    } while(charge[partIndex] == 0.);

    return (Particle){ Pt[partIndex]/charge[partIndex], Eta[partIndex], Phi[partIndex], die.Uniform(-Z0_SMEAR_MM, Z0_SMEAR_MM) };
}

Long64_t ParticleGenerator::getNumEntries() const { return numEntries; }




void TrackShooter::setOutput(ostream& output, const char* fieldSeparator, const char* lineSeparator, bool) {
  output_ = &output;
  FS = fieldSeparator;
  LS = lineSeparator;
}

void TrackShooter::setTrackerBoundaries(double trackerMaxRho, double barrelMinZ, double barrelMaxZ) {
  trackerMaxRho_ = trackerMaxRho;
  barrelMinZ_ = barrelMinZ;
  barrelMaxZ_ = barrelMaxZ;
}  

void TrackShooter::addModule(Module* module) { // also locks modules geometry so that only cached values for geometry properties are returned from now on
  BarrelModule* bmod; EndcapModule* emod;
  module->lockGeometry();
  allMods_.push_back(module);
  
  const std::set<int>& octants = getModuleOctants(module);

  if ((bmod=dynamic_cast<BarrelModule*>(module))) {
    BarrelOctants& barrelOctants = barrelModsByRadius_[module->getMeanPoint().Rho()-module->getStereoDistance()/2]; // the radius is the radius of the center point of the lower module
    barrelOctants.resize(8); // prepare the octants vector. it's going to be called for each module, but after the first time the call will do nothing
    for (std::set<int>::const_iterator oit = octants.begin(); oit != octants.end(); ++oit) {
      barrelOctants[*oit].push_back(bmod);
    }
  } else if ((emod=dynamic_cast<EndcapModule*>(module))) { 
    EndcapOctants& endcapOctants = endcapModsByZ_[module->getMeanPoint().Z()-(module->getZSide()*module->getStereoDistance()/2)]; // getZSide() changes the sign of the addition because if mod is in Z+ then the first sensor to be hit is the one at lower Z, viceversa for Z-
    endcapOctants.resize(8);
    for (std::set<int>::const_iterator oit = octants.begin(); oit != octants.end(); ++oit) {
      endcapOctants[*oit].push_back(emod);
    }
  }

}

XYVector TrackShooter::convertToLocalCoords(const XYZVector& globalHit, const BarrelModule* mod) const {
  XYZVector locv = RotationZ(-mod->getMeanPoint().Phi())*(globalHit - mod->getMeanPoint()); // we translate the hit back to the origin, we rotate it like the module at phi=0 was hit (which now looks like a vertical segment in the XY plane), then we drop the irrelevant coordinate (X)
  return XYVector(locv.Y(), locv.Z());
}

XYVector TrackShooter::convertToLocalCoords(const XYZVector& globalHit, const EndcapModule* mod) const {
  XYZVector locv = RotationZ(-mod->getMeanPoint().Phi())*(globalHit - mod->getMeanPoint()); // we translate the hit back to the origin, we rotate it like the module at phi=0 was hit (which now looks like a polygon in the XY plane with the local axes XY reverse-matched with the global YX), then we drop the irrelevant coordinate (Z)
  return XYVector(locv.Y(), locv.X());
}

/*
int detectCollisionSlanted(XYVector& helixAxis, double R, double z0, double theta, const Polygon3d<4>& poly, std::vector<XYZVector>& collisions) {
  std::vector<XYZVector> bcol, ecol;
  detectCollisionBarrel(helixAxisOrigin, poly, bcol);
  detectCollisionEndcap(helixAxisOrigin, poly, ecol);

  double phi0 = origin.Phi() - M_PI/2;
  double h = helixAxis.X();
  double k = helixAxis.Y();

  double L = tan(M_PI/2 - theta);
  if (bcol.size() > 0 || ecol.size() > 0) {
    for (int i = 0; i < bcol.size(); i++) {
      double zz = (ecol.size() > 0) ? ecol[0].Z()*sin(mod.getGamma()) + bcol[i].Z()*cos(mod.getGamma()) : bcol[i].Z();
      XYZVector& P = mod.getCorner(0);
      XYZVector& N = mod.getNormal();
      XYZVector Q;
      do {
        double tt = (zz - z0)/(L*R);
        XYZVector H(h + R*sin(phi0 + tt), k - R*cos(phi0 + tt), zz);
        Q = H - N.Dot(H - P) * N;
        zz = Q.Z(); 
      } while (N.Dot(Q - P) > 1e-3)
      collisions.push_back(Q);
    }
  }  
  return collisions.size();
}
*/

int TrackShooter::detectCollisionEndcap(double h, double k, double pt, double z0, double phi0, double theta, const Polygon3d<4>& poly, std::vector<XYZVector>& collisions) {
  // get z from module
  // put z in the helix equations, solve for t and obtain x and y
  // check whether x and y are inside the polygon of the module (only 1 collision is possible)
  
  double R = pt/(0.3*insur::magnetic_field) * 1e3;

  double z = poly.getVertex(0).Z();
  
  double L = tan(M_PI/2 - theta);

  double t = (z - z0)/(L*R);
  double x = h + R*sin(phi0 + t);
  double y = k - R*cos(phi0 + t);
  
  XYZVector hit(x, y, z);
  if (poly.isPointInside(hit)) {
    collisions.push_back(hit);
    return 1;
  } else {
    return 0; 
  }
}

int TrackShooter::detectCollisionBarrel(double h, double k, double pt, double z0, double phi0, double theta, const Polygon3d<4>& poly, std::vector<XYZVector>& collisions) {
  // this method solves the equations of the intersection between a circle (particle track in the XY plane) and a line segment (module in XY)
  // both equation are expressed in their parametric forms, to aid in boundary checking.
  // In case of the mocule the parameter u must be between 0 and 1 for the hit to be inside the module
  // In case of the track the parameter t helps in understanding which of the two intersections is the one to consider given the rotation direction of the track (determined by its charge)
  // (a circle can also intersect a line in 0 or just 1 point, but those cases are dealt with separately)
  //
  // the resulting system of equations is:
  // x = h + R*cos(t)
  // y = k + R*sin(t)
  // x = v + m*u
  // y = w + n*u
  //
  // where: h, k are the center of the track circle, R is its radius
  //        v, w is the starting point of the module segment, m, n are its lengths along the x and y axes (m = endpointX - v, n = endpointY - w)
  //
  // the first two equations are put together with the other two over x and y
  // 
  // The method returns the x,y,t of the intersections (2, 1 or possibly a null intersection if the track never hits the module) 

  double R = pt/(0.3*insur::magnetic_field) * 1e3;
  if (fabs(2*R) < poly.getCenter().Rho()) return 0; // not even a tangent collision is possible (this is only valid inasmuch the center of the module is the closest point to the beam axis - always true for barrel modules)

  double v = poly.getVertex(0).X(),     w = poly.getVertex(0).Y();
  double m = poly.getVertex(1).X() - v, n = poly.getVertex(1).Y() - w;

  double minZ = MIN(poly.getVertex(0).Z(), poly.getVertex(2).Z());
  double maxZ = MAX(poly.getVertex(0).Z(), poly.getVertex(2).Z());

  if (theta < M_PI/2 && maxZ < z0) return 0; // this filters out modules which will nevel be hit because they lie in the opposite direction as the particle's trajectory
  else if (theta > M_PI/2 && minZ > z0) return 0; // ditto here

  double a = m*m + n*n;
  double b = 2*(m*v + n*w - h*m - k*n);
  double c = h*h + k*k + v*v + w*w - 2*h*v - 2*k*w - R*R;
  double delta = b*b - 4*a*c;
  if (delta < 0.) return 0; // track can never collide with module

  double u1 = (-b + sqrt(delta))/(2*a);
  double u2 = (-b - sqrt(delta))/(2*a);
 
  double L = tan(M_PI/2 - theta);

  int initialCollisionSize = collisions.size();
  //if (phi0 < 0) phi0 += 2*M_PI;
  if (u1 >= 0. && u1 <= 1.) {
    double x = v + m*u1;
    double y = w + n*u1;
    double t1_1 = fmod(asin((v + m*u1 - h)/R) - phi0, 2*M_PI); // each u root, plugged in the first parametric circle eq, results in 2 candidate t roots (because of sine, which in [0,2pi] always has two angles resulting in the same sine value)
    t1_1 = t1_1 - (pt < 0 && t1_1 > 0 ? 2*M_PI : 0.); // wrapping the root in 0,2pi, if pt is < 0 we want negative numbers
    double t1_2 = fmod(M_PI - asin((v + m*u1 - h)/R) - phi0, 2*M_PI);
    t1_2 = t1_2 - (pt < 0 && t1_2 > 0 ? 2*M_PI : 0.);
    double y1 = k - R*cos(phi0 + t1_1);
    double y2 = k - R*cos(phi0 + t1_2);
    y2 = y2;
    double t1 = fabs(y1 - y) < 1e-3 ? t1_1 : t1_2; // to figure out which one we want, we plug the first candidate in the second parametric circle eq and we check whether the resulting y has the same sign of the y of the intersection point. If signs differ, the good t root is the other one. This is because the 2 two roots result in cosines with same modulo but opposite sign.
    double z = z0 + L*R*t1;
    if (fabs(pt) >= HIGH_PT_THRESHOLD) {
      if (z >= minZ && z <= maxZ && fabs(t1) < M_PI/2) collisions.push_back(XYZVector(x, y, z)); // the condition on t1 weeds out tracks curving back in the detector (in a coarse way)
    } else {
      double imin = (minZ - z0 - L*R*t1)/(L*R*t1*2*M_PI);
      double imax = (maxZ - z0 - L*R*t1)/(L*R*t1*2*M_PI);
      if (signum(pt) == signum(imax)) {
        for (int i = ceil(imin); i <= floor(imax); i++) {
          collisions.push_back(XYZVector(x, y, z0 + L*R*t1*(1 + 2*M_PI*i)));
        }
      }
    }
  }
  if (u2 >= 0. && u2 <= 1.) {
    double x = v + m*u2;
    double y = w + n*u2;
    double t2_1 = fmod(asin((v + m*u2 - h)/R) - phi0, 2*M_PI);
    t2_1 = t2_1 - (pt < 0 && t2_1 > 0 ? 2*M_PI : 0.); // each u root, plugged in the first parametric circle eq, results in 2 candidate t roots (because of sine, which in [0,2pi] always has two angles resulting in the same sine value)
    double t2_2 = fmod(M_PI - asin((v + m*u2 - h)/R) - phi0, 2*M_PI);
    t2_2 = t2_2 - (pt < 0 && t2_2 > 0 ? 2*M_PI : 0.);
    double y1 = k - R*cos(phi0 + t2_1);
    double y2 = k - R*cos(phi0 + t2_2);
    y2 = y2;
    double t2 = fabs(y1 - y) < 1e-3 ? t2_1 : t2_2; // to figure out which one we want, we plug the first candidate in the second parametric circle eq and we check whether the resulting y has the same sign of the y of the intersection point. If signs differ, the good t root is the other one. This is because the 2 two roots result in cosines with same modulo but opposite sign.
    double z = z0 + L*R*t2;
    if (fabs(pt) >= HIGH_PT_THRESHOLD) {
      if (z >= minZ && z <= maxZ && fabs(t2) < M_PI/2) collisions.push_back(XYZVector(x, y, z));
    } else { 
      double imin = (minZ - z0 - L*R*t2)/(L*R*t2*2*M_PI);
      double imax = (maxZ - z0 - L*R*t2)/(L*R*t2*2*M_PI);
      if (signum(pt)*signum(imax)) {
        for (int i = ceil(imin); i <= floor(imax); i++) {
          collisions.push_back(XYZVector(x, y, z0 + L*R*t2*(1 + 2*M_PI*i)));
        }
      }
    }
  }

  return collisions.size() - initialCollisionSize;
}



//#define FAKE_HITS

void TrackShooter::shootTracks(long int numEvents, long int numTracksEv, int seed) {

  numEvents_ = numEvents;
  numTracksEv_ = numTracksEv;

  die_.SetSeed(seed);

  shootTracks();
}


void TrackShooter::setDefaultParameters() {

  numEvents_ = 1000;
  numTracksEv_ = 1; 
  eventOffset_ = 0;

  eta_ = std::auto_ptr<Value<double> >(new UniformValue<double>(die_, -2, 2));
  phi0_ = std::auto_ptr<Value<double> >(new UniformValue<double>(die_, M_PI, -M_PI));
  pt_ = std::auto_ptr<Value<double> >(new UniformValue<double>(die_, 2, 50));
  z0_ = std::auto_ptr<Value<double> >(new UniformValue<double>(die_, 0.01, 0.5));
  charge_ = std::auto_ptr<Value<int> >(new BinaryValue<int>(die_, -1, 1));

  instanceId_ = any2str(getpid()) + "_" + any2str(time(NULL)); 
  tracksDir_ = ".";

  useInvPt_ = false;
}


void TrackShooter::shootTracks(const po::variables_map& varmap, int seed) {

//  std::string line;

//  while (line = std::getline(cfg, line, ";")) {
//    if (line.find("//") != string::npos) line.erase(line.find("//")); // remove comments
//    line = trim(line);
//    if (line.empty()) continue;
//    std::vector<string> tokens = split(line, "=,");
  for (po::variables_map::const_iterator it = varmap.begin(); it != varmap.end(); ++it) {
    std::string key(it->first);
    if (key == "eta") eta_ = valueFromString<double>(die_, it->second.as<std::string>());
    else if (key == "phi0") phi0_ = valueFromString<double>(die_, it->second.as<std::string>());
    else if (key == "z0") z0_ = valueFromString<double>(die_, it->second.as<std::string>());
    else if (key == "pt") {
      useInvPt_ = false; // only either pt or invPt can be specified
      pt_ = valueFromString<double>(die_, it->second.as<std::string>());
    } else if (key == "invPt") {
      useInvPt_ = true;
      invPt_ =  valueFromString<double>(die_, it->second.as<std::string>());
    } else if (key == "charge") charge_ = valueFromString<int>(die_, it->second.as<std::string>());
    else if (key == "num-events") numEvents_ = str2any<long int>(it->second.as<std::string>()); 
    else if (key == "num-tracks-ev") numTracksEv_ = str2any<long int>(it->second.as<std::string>());
    else if (key == "event-offset") eventOffset_ = str2any<long int>(it->second.as<std::string>());
    else if (key == "instance-id") {
      static const std::string pidtag = "%PID%";
      static const std::string timetag = "%TIME%";
      size_t pos;
      instanceId_ = it->second.as<std::string>();
      if ((pos = instanceId_.find(timetag)) != std::string::npos) instanceId_.replace(pos, timetag.size(), any2str(time(NULL)));
      if ((pos = instanceId_.find(pidtag)) != std::string::npos) instanceId_.replace(pos, pidtag.size(), any2str(getpid()));
    } else if (key == "tracks-dir") tracksDir_ = it->second.as<std::string>();
  
  }

  die_.SetSeed(seed);

  shootTracks();
}


void TrackShooter::printParameters() {
  std::cout << "\nSimulation parameters summary" << std::endl;
  std::cout << "num-events = " << numEvents_ << std::endl;
  std::cout << "num-tracks-ev = " << numTracksEv_ << std::endl;
  std::cout << "event-offset = " << eventOffset_ << std::endl;
  std::cout << "eta = " << eta_->toString() << std::endl;
  std::cout << "phi0 = " << phi0_->toString() << std::endl;
  std::cout << "z0 = " << z0_->toString() << std::endl;
  std::cout << "pt = " << (!useInvPt_ ? pt_->toString() : "n/a") << std::endl;
  std::cout << "inv-pt = " << (useInvPt_ ? invPt_->toString() : "n/a") << std::endl;
  std::cout << "charge = " << charge_->toString() << std::endl;
  std::cout << "instance-id = " << instanceId_ << std::endl;
  std::cout << "tracks-dir = " << tracksDir_ << std::endl;
  std::cout << "rand-seed = " << die_.GetSeed() << std::endl;
}


void TrackShooter::exportGeometryData() {

  ModuleData mdata;

  TTree* tree = new TTree("geomdata", "Geometry data");
  tree->Branch("mdata", &mdata, "x/D:y:z:rho:phi:widthlo:widthhi:height:stereo:pitchlo:pitchhi:striplen:yres:inefftype/B:refcnt:refz:refrho:refphi:type"); 

  for (std::vector<Module*>::const_iterator it = allMods_.begin(); it != allMods_.end(); ++it) {
    Module* mod = (*it);
    PosRef posref = mod->getPositionalReference();
    XYZVector center = mod->getMeanPoint();
    mdata = (ModuleData){ center.X(), center.Y(), center.Z(),
                          center.Rho(), center.Phi(),
                          mod->getWidthLo(), mod->getWidthHi(), mod->getHeight(),
                          mod->getStereoDistance(),
                          mod->getLowPitch(), mod->getHighPitch(),
                          mod->getStripLength(), 
                          mod->getResolutionYTrigger(),
                          mod->getModuleType()->getInefficiencyType(),
                          posref.cnt, posref.z, posref.rho, posref.phi,
                          mod->getSubdetectorType() };

    tree->Fill();
  }
}

void TrackShooter::shootTracks() {

  Tracks tracks("tracks");
  Hits hits("hits");
  Hits plhits("plhits");

  printParameters();

  std::string outfileName = tracksDir_ + "/tracks_" + instanceId_ + ".root";

  TFile* outfile = new TFile(outfileName.c_str(), "recreate");
  if (outfile->IsZombie()) {
    std::cerr << "Failed opening file \"" << outfileName << "\" for writing. Simulation aborted." << std::endl;
    return;
  }

  exportGeometryData();

  TTree* tree = new TTree("trackhits", "TTree containing hits of simulated tracks");

  gROOT->ProcessLine("#include <vector>");

  tracks.setupBranches(*tree);
  hits.setupBranches(*tree);
  plhits.setupBranches(*tree);


  int numpl = 0, numcyl = 0;
  // build ordered maps
  for (long int i=eventOffset_, totTracks = eventOffset_*numTracksEv_; i<numEvents_+eventOffset_; i++) {
    // new event
    //
    for (long int j=0; j<numTracksEv_; j++, totTracks++) {
      double eta = eta_->get();
      double phi0 = phi0_->get();
      double z0 = z0_->get();
      double pt = charge_->get() * (!useInvPt_ ? pt_->get() : 1./invPt_->get());

#ifdef FAKE_HITS
      if ((totTracks % 500) == 0) cout << "Track " << totTracks << " of " << (numEvents_+eventOffset_)*numTracksEv_ << std::endl;
#else
      if ((totTracks % 5000) == 0) cout << "Track " << totTracks << " of " << (numEvents_+eventOffset_)*numTracksEv_ << std::endl;
#endif
/*      
      // or pick from file
      //
      //
      ParticleGenerator::Particle particle = partGen.getParticle();
      
      double pt = particle.pt;
      double eta = particle.eta;
      double phi0 = particle.phi;
      double z0 = particle.z0;
  */    

      int dir = signum(pt);
      double theta = 2*atan(exp(-eta));
      double R = fabs(pt)/(0.3*insur::magnetic_field) * 1e3;
      double B = tan(theta)/R;

      double helixCenterX = -dir*R*sin(phi0);
      double helixCenterY =  dir*R*cos(phi0);

      std::vector<XYZVector> collisions;
      collisions.clear();

      for (BarrelRadii::const_iterator rit = barrelModsByRadius_.begin(); rit != barrelModsByRadius_.end(); ++rit) {
        double r = rit->first;
        if (r >= 2*R) continue;
        double z = 1/B*acos(1-r*r/(2*R*R)) + z0;
        double x = R*sin(B*(z-z0));
        double y = dir*R*(1-cos(B*(z-z0)));
        double xrot = x*cos(phi0) - y*sin(phi0);
        double yrot = x*sin(phi0) + y*cos(phi0);
        double phirot = atan2(y,x) + phi0;
        if (phirot > M_PI) phirot -= 2*M_PI;

        const BarrelOctants& octants = rit->second;
        const BarrelModules& bmods = octants[getPointOctant(xrot, yrot, z)]; // jump to the octant of the point (CUIDADO octant is determined used non-planar mods)

        for (BarrelModules::const_iterator mit = bmods.begin(); mit != bmods.end(); ++mit) {
          BarrelModule* mod = (*mit);

          Polygon3d<4> poly;
          poly << mod->getCorner(0) << mod->getCorner(1) << mod->getCorner(2) << mod->getCorner(3);
          double xpl;
          double ypl;
          double zpl;
          bool planarcoll = false;
          if (detectCollisionBarrel(helixCenterX, helixCenterY, pt, z0, phi0, theta, poly, collisions)) { // planar collisions 
            xpl = collisions[0].X();
            ypl = collisions[0].Y();
            zpl = collisions[0].Z(); 
            ptError* modPtError = mod->getPtError();
            float pterr = modPtError->computeError(pt);
            float hitprob = mod->getTriggerProbability(pt);
            float deltaStrips = modPtError->pToStrips(pt);
            XYVector locv = convertToLocalCoords(collisions[0], mod); 
            PosRef posref = mod->getPositionalReference();
            plhits.push_back(xpl, ypl, zpl, locv.X(), locv.Y(), pterr, hitprob, deltaStrips, posref.cnt, posref.z, posref.rho, posref.phi, -1., -1.);
            collisions.clear();
            planarcoll = true;
            numpl++;
            xpl = locv.X();
            ypl = locv.Y();
          }

          double minPhi = mod->getMinPhi();
          double maxPhi = mod->getMaxPhi();
          if (mod->getMinZ() <= z && z <= mod->getMaxZ() &&
              minPhi <= phirot && phirot <= maxPhi) {
            ptError* modPtError = mod->getPtError();
            float pterr = modPtError->computeError(pt);
            float hitprob = mod->getTriggerProbability(pt);
            float deltaStrips = modPtError->pToStrips(pt);
            mod->setProperty("tracksimHits", mod->getProperty("tracksimHits")+1);
            XYVector locv = convertToLocalCoords(XYZVector(xrot, yrot, z), mod); 
            PosRef posref = mod->getPositionalReference();
           // double dist = planarcoll ? sqrt((xrot - xpl)*(xrot - xpl) + (yrot - ypl)*(yrot - ypl) + (z - zpl)*(z - zpl)) : -1.;
            double distx = planarcoll ? fabs(locv.X() - xpl) : -1.;
            double disty = planarcoll ? fabs(locv.Y() - ypl) : -1.;
            hits.push_back(xrot, yrot, z, locv.X(), locv.Y(), pterr, hitprob, deltaStrips, posref.cnt, posref.z, posref.rho, posref.phi, distx, disty);

            numcyl++;
            //if (fabs(pt) >= HIGH_PT_THRESHOLD) break; // high pT particles never curve back inside the detector so after a layer/disk has been hit it makes no sense to look for more hits in modules in the same layer/disk
            // CUIDADO this optimization has been disabled for debug
          }
        }
      }

      if (2*R > trackerMaxRho_) { // track will at some point escape the tracker
        double z = 1/B*acos(1-(trackerMaxRho_*trackerMaxRho_)/(2*R*R)) + z0;
        if (barrelMinZ_ <= z && z <= barrelMaxZ_) { // check whether the track will escape from the barrel volume
          tracks.push_back(i, j, eta, phi0, z0, pt, hits.size());
          tree->Fill();
          hits.clear();
          tracks.clear();
          continue; // particle has escaped the detector from the barrel volume, we don't want the endcaps to see escaped particles curving back into the tracker, so we skip on
        }
      }

      for (EndcapZs::const_iterator zit = endcapModsByZ_.begin(); zit != endcapModsByZ_.end(); ++zit) {
        double z = zit->first;
        if (signum(eta)*signum(z) <= 0) continue; // we need to skip the negative tracker section if the particle is headed towards positive Z's and viceversa, as the trajectory is an infinite sinusoid
        double x = R*sin(B*(z-z0));
        double y = dir*R*(1-cos(B*(z-z0)));
        double xrot = x*cos(phi0) - y*sin(phi0);
        double yrot = x*sin(phi0) + y*cos(phi0);

        const EndcapOctants& octants = zit->second;
        const EndcapModules& emods = octants[getPointOctant(xrot, yrot, z)];


        for (EndcapModules::const_iterator mit = emods.begin(); mit != emods.end(); ++mit) {
          EndcapModule* mod = (*mit);

          Polygon3d<4> poly;
          poly << (mod->getCorner(0) - XYZVector(0, 0, mod->getZSide()*mod->getStereoDistance()/2))
            << (mod->getCorner(1) - XYZVector(0, 0, mod->getZSide()*mod->getStereoDistance()/2))
            << (mod->getCorner(2) - XYZVector(0, 0, mod->getZSide()*mod->getStereoDistance()/2))
            << (mod->getCorner(3) - XYZVector(0, 0, mod->getZSide()*mod->getStereoDistance()/2)); 

          if (detectCollisionEndcap(helixCenterX, helixCenterY, pt, z0, phi0, theta, poly, collisions)) {
            ptError* modPtError = mod->getPtError();
            float pterr = modPtError->computeError(pt);
            float hitprob = mod->getTriggerProbability(pt);
            float deltaStrips = modPtError->pToStrips(pt);
            PosRef posref = mod->getPositionalReference();
            double xpl = collisions[0].X();
            double ypl = collisions[0].Y();
            double zpl = collisions[0].Z();
            XYVector locv = convertToLocalCoords(collisions[0], mod);
            plhits.push_back(xpl, ypl, zpl, locv.X(), locv.Y(), pterr, hitprob, deltaStrips, posref.cnt, posref.z, posref.rho, posref.phi, -1., -1.);
            collisions.clear();
          }

          XYZVector hitv(xrot, yrot, z);
          // check if hit lies in module

          if (poly.isPointInside(hitv)) {
            // module was hit
            ptError* modPtError = mod->getPtError();
            float pterr = modPtError->computeError(pt);
            float hitprob = mod->getTriggerProbability(pt);
            float deltaStrips = modPtError->pToStrips(pt);
            mod->setProperty("tracksimHits", mod->getProperty("tracksimHits")+1);
            PosRef posref = mod->getPositionalReference();
            XYVector locv = convertToLocalCoords(XYZVector(xrot, yrot, z), mod);

            hits.push_back(xrot, yrot, z, locv.X(), locv.Y(), pterr, hitprob, deltaStrips, posref.cnt, posref.z, posref.rho, posref.phi, -1., -1.);

            //if (fabs(pt) >= HIGH_PT_THRESHOLD) break; // high pT particles never curve back inside the detector so after a layer/disk has been hit it makes no sense to look for more hits in modules in the same layer/disk
          }      
        }

      }
      tracks.push_back(i, j, eta, phi0, z0, pt, hits.size());

      tree->Fill();

      hits.clear();
      tracks.clear();
    }     
#ifdef FAKE_HITS
    // NOW BROKEN! ENABLING THEM IS USELESS!!
    for (std::vector<Module*>::iterator mit = allMods_.begin(); mit != allMods_.end(); ++mit) {
      Module* mod = (*mit);
      // loop over all the modules again
      // fake hits
      Polygon3d<4> poly;
      poly << mod->getCorner(0)-mod->getMeanPoint() << mod->getCorner(1)-mod->getMeanPoint() << mod->getCorner(2)-mod->getMeanPoint() << mod->getCorner(3)-mod->getMeanPoint();

      int numFake = die.Poisson(mod->getTriggerFrequencyFakePerEvent());
      for (int k = 0; k < numFake; k++) {
        XYZVector fakehit = poly.generateRandomPoint(&die);
        float deltaStrips = 1 + die.Integer(mod->getTriggerWindow()/2);
        float fakept = mod->getPtError()->stripsToP(deltaStrips);
        float pterr = mod->getPtError()->computeError(fakept);
        PosRef posref = mod->getPositionalReference();
      }
    }
#endif
  }
  std::cout << "numpl=" << numpl << " numcyl=" << numcyl << std::endl;

  outfile->Write();
  outfile->Close();
  delete outfile;

  std::cout << "Output written to file " << outfileName << std::endl;
}



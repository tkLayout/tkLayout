#include <TrackShooter.h>


std::set<int> getModuleOctants(const Module* mod) {
  std::set<int> octants;
  octants.insert(getPointOctant(mod->getCorner(0).X(), mod->getCorner(0).Y(), mod->getCorner(0).Z())); 
  octants.insert(getPointOctant(mod->getCorner(1).X(), mod->getCorner(1).Y(), mod->getCorner(1).Z())); 
  octants.insert(getPointOctant(mod->getCorner(2).X(), mod->getCorner(2).Y(), mod->getCorner(2).Z())); 
  octants.insert(getPointOctant(mod->getCorner(3).X(), mod->getCorner(3).Y(), mod->getCorner(3).Z())); 
  return octants;
}




template<int NumSides, class Coords, class Random>
AbstractPolygon<NumSides, Coords, Random>& AbstractPolygon<NumSides, Coords, Random>::operator<<(const Coords& vertex) {
  if (!isComplete()) v_.push_back(vertex);
  if (isComplete()) computeProperties();
  return *this;
}

template<int NumSides, class Coords, class Random>
AbstractPolygon<NumSides, Coords, Random>& AbstractPolygon<NumSides, Coords, Random>::operator<<(const std::vector<Coords>& vertices) {
  for(typename std::vector<Coords>::const_iterator it = vertices.begin(); (v_.size() < NumSides) && (it != vertices.end()); ++it) {
    *this << *it;
  }
  return *this;
}

template<int NumSides, class Coords, class Random>
inline bool AbstractPolygon<NumSides, Coords, Random>::isComplete() const { 
  return v_.size() == NumSides; 
}

template<int NumSides, class Coords, class Random>
inline const std::vector<Coords>& AbstractPolygon<NumSides, Coords, Random>::getVertices() const { 
  return v_; 
}

template<int NumSides, class Coords, class Random>
inline const Coords& AbstractPolygon<NumSides, Coords, Random>::getVertex(int index) const { 
  return v_[index]; 
}

template<int NumSides, class Coords, class Random>
double AbstractPolygon<NumSides, Coords, Random>::getDoubleArea() const {
  return area_;
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
    BarrelOctants& barrelOctants = barrelModsByRadius_[module->getMeanPoint().Rho()];
    barrelOctants.resize(8); // prepare the octants vector. it's going to be called for each module, but after the first time the call will do nothing
    for (std::set<int>::const_iterator oit = octants.begin(); oit != octants.end(); ++oit) {
      barrelOctants[*oit].push_back(bmod);
    }
  } else if ((emod=dynamic_cast<EndcapModule*>(module))) { 
    EndcapOctants& endcapOctants = endcapModsByZ_[module->getMeanPoint().Z()];
    endcapOctants.resize(8);
    for (std::set<int>::const_iterator oit = octants.begin(); oit != octants.end(); ++oit) {
      endcapOctants[*oit].push_back(emod);
    }
  }

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


void TrackShooter::shootTracks() {

  Tracks tracks;
  Hits hits;

  printParameters();

  std::string outfileName = tracksDir_ + "/tracks_" + instanceId_ + ".root";

  TFile* outfile = new TFile(outfileName.c_str(), "recreate");
  if (outfile->IsZombie()) {
    std::cerr << "Failed opening file \"" << outfileName << "\" for writing. Simulation aborted." << std::endl;
    return;
  }
  TTree* tree = new TTree("trackhits", "TTree containing hits of simulated tracks");

  gROOT->ProcessLine("#include <vector>");

  tree->Branch("tracks.eventn", &tracks.eventn);
  tree->Branch("tracks.trackn", &tracks.trackn);
  tree->Branch("tracks.eta", &tracks.eta);
  tree->Branch("tracks.phi0", &tracks.phi0);
  tree->Branch("tracks.z0", &tracks.z0);
  tree->Branch("tracks.pt", &tracks.pt);

  tree->Branch("hits.glox", &hits.glox);
  tree->Branch("hits.gloy", &hits.gloy);
  tree->Branch("hits.gloz", &hits.gloz);
  tree->Branch("hits.locx", &hits.locx);
  tree->Branch("hits.locy", &hits.locy);
  tree->Branch("hits.pterr", &hits.pterr);
  tree->Branch("hits.hitprob", &hits.hitprob);
  tree->Branch("hits.deltas", &hits.deltas);
  tree->Branch("hits.cnt", &hits.cnt);
  tree->Branch("hits.z", &hits.z);
  tree->Branch("hits.rho", &hits.rho);
  tree->Branch("hits.phi", &hits.phi);


  // build ordered maps
  for (long int i=eventOffset_, totTracks = eventOffset_*numTracksEv_; i<numEvents_+eventOffset_; i++) {
    // new event
    //
    for (long int j=0; j<numTracksEv_; j++, totTracks++) {
      // randomly generate particle with eta = [-3,+3], phi = [0,2pi], Pt = [0.6,15]
      //double eta = die.Uniform(-3, 3);
      double eta = eta_->get(); // die.Uniform(-2, 2);
      double phi0 = phi0_->get(); //die.Uniform(0.38, 1.634);
      double z0 = z0_->get(); //die.Uniform(-20, 20);
      double pt = charge_->get() * (!useInvPt_ ? pt_->get() : 1./invPt_->get()); //(die.Integer(2) ? 1 : -1) * 1./die.Uniform(0.01, 0.5);

      // or pick from file
      //
      //
#ifdef FAKE_HITS
      if ((totTracks % 500) == 0) cout << "Track " << totTracks << " of " << (numEvents_+eventOffset_)*numTracksEv_ << std::endl;
#else
      if ((totTracks % 5000) == 0) cout << "Track " << totTracks << " of " << (numEvents_+eventOffset_)*numTracksEv_ << std::endl;
#endif
/*      
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

      for (BarrelRadii::const_iterator rit = barrelModsByRadius_.begin(); rit != barrelModsByRadius_.end(); ++rit) {
        double r = rit->first;
        if (r >= 2*R) continue;
        double z = 1/B*acos(1-r*r/(2*R*R)) + z0;
        double x = R*sin(B*(z-z0));
        double y = dir*R*(1-cos(B*(z-z0)));
        double xrot = x*cos(phi0) - y*sin(phi0);
        double yrot = x*sin(phi0) + y*cos(phi0);
        double phirot = atan2(y,x) + phi0;

        const BarrelOctants& octants = rit->second;
        const BarrelModules& bmods = octants[getPointOctant(xrot, yrot, z)]; // jump to the octant of the point 

        for (BarrelModules::const_iterator mit = bmods.begin(); mit != bmods.end(); ++mit) {
          BarrelModule* mod = (*mit);
          double minPhi = mod->getMinPhi();
          double maxPhi = mod->getMaxPhi();
          if (mod->getMinZ() <= z && z <= mod->getMaxZ() &&
              minPhi <= phirot && phirot <= maxPhi) {
            ptError* modPtError = mod->getPtError();
            float pterr = modPtError->computeError(pt);
            float hitprob = mod->getTriggerProbability(pt);
            int deltaStrips = modPtError->pToStrips(pt);
            mod->setProperty("tracksimHits", mod->getProperty("tracksimHits")+1);
            PosRef posref = mod->getPositionalReference();
            XYZVector locv = RotationZ(-mod->getMeanPoint().Phi())*(XYZVector(xrot, yrot, z) - mod->getMeanPoint()); // we translate the hit back to the origin, we rotate it like the module at phi=0 was hit, then we drop the irrelevant coordinate (X)
            
            hits.push_back(xrot, yrot, z, locv.Y(), locv.Z(), pterr, hitprob, deltaStrips, posref.cnt, posref.z, posref.rho, posref.phi);

            if (fabs(pt) >= HIGH_PT_THRESHOLD) break; // high pT particles never curve back inside the detector so after a layer/disk has been hit it makes no sense to look for more hits in modules in the same layer/disk
          }
        }
      }

      if (2*R > trackerMaxRho_) { // track will at some point escape the tracker
        double z = 1/B*acos(1-(trackerMaxRho_*trackerMaxRho_)/(2*R*R)) + z0;
        if (barrelMinZ_ <= z && z <= barrelMaxZ_) { // check whether the track will escape from the barrel volume
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
          Module* mod = (*mit);
          XYZVector hitv(xrot, yrot, z);
          // check if hit lies in module
          Polygon3d<4> poly;
          poly << mod->getCorner(0) << mod->getCorner(1) << mod->getCorner(2) << mod->getCorner(3); 
          if (poly.isPointInside(hitv)) {
            // module was hit
            ptError* modPtError = mod->getPtError();
            float pterr = modPtError->computeError(pt);
            float hitprob = mod->getTriggerProbability(pt);
            int deltaStrips = modPtError->pToStrips(pt);
            mod->setProperty("tracksimHits", mod->getProperty("tracksimHits")+1);
            PosRef posref = mod->getPositionalReference();
            XYZVector locv = RotationZ(-mod->getMeanPoint().Phi())*(hitv - mod->getMeanPoint());  // we translate the hit back to the origin, we rotate it like the module at phi=0 was hit, then we drop the irrelevant coordinate (Z)

            hits.push_back(xrot, yrot, z, locv.Y(), locv.X(), pterr, hitprob, deltaStrips, posref.cnt, posref.z, posref.rho, posref.phi);

            if (fabs(pt) >= HIGH_PT_THRESHOLD) break; // high pT particles never curve back inside the detector so after a layer/disk has been hit it makes no sense to look for more hits in modules in the same layer/disk
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
        int deltaStrips = 1 + die.Integer(mod->getTriggerWindow()/2);
        float fakept = mod->getPtError()->stripsToP(deltaStrips);
        float pterr = mod->getPtError()->computeError(fakept);
        PosRef posref = mod->getPositionalReference();
      }
    }
#endif
  }

  outfile->Write();
  outfile->Close();
  delete outfile;

  std::cout << "Output written to file " << outfileName << std::endl;
}



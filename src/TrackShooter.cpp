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
const std::multiset<Triangle3d, PolygonLess<Triangle3d> >& Polygon3d<NumSides>::getTriangulation() const {
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

//  PosRef p = module->getPositionalReference();
  //cout << "Mod posref= " << (int)p.cnt << "," << (int)p.z << "," << (int)p.rho << "," << (int)p.phi << " coords=" << module->getMeanPoint().Z() << "," << module->getMeanPoint().Rho() << "," << module->getMeanPoint().Phi() << endl; 
}

//#define VDEBUG

//#define FAKE_HITS

void TrackShooter::shootTracks(long int numEvents, long int numTracksEv) {
#ifdef VDEBUG
  TCanvas* Canvas3d = new TCanvas("Canvas3d", "Badass 3d canvas", 1200, 1200);
  TView *view = TView::CreateView(1);
  view->SetRange(0,0,0,1300,600,600);

  std::vector<float> modXYZ;

  TCanvas* RZCanvas = new TCanvas("RZCanvasTS", "RZView Canvas", 1200, 600 );
  RZCanvas->cd();
  PlotDrawer<YZ, Property, Sum> yzDrawer(0, 0, Property("tracksimHits"));

  TCanvas* XYCanvas = new TCanvas("XYCanvasTS", "XYView Canvas", 600, 600 );
  XYCanvas->cd();
  PlotDrawer<XY, Property, Sum> xyBarrelDrawer(0, 0, Property("tracksimHits"));

  TCanvas* XYCanvasEC = new TCanvas("XYCanvasECTS", "XYView Canvas (Endcap)", 600, 600 );
  XYCanvasEC->cd();
  PlotDrawer<XY, Property, Sum> xyEndcapDrawer(0, 0, Property("tracksimHits")); 
#endif
  //*output_ << LS;
  //*output_ << "P" << FS << "ev#" << FS << "tr#" << FS << "eta" << FS << "phi0" << FS << "z0" << FS << "pt" << LS;
  //*output_ << "T/F" << FS << "ev#" << FS << "tr#" << FS << "eta" << FS << "phi0" << FS << "z0" << FS << "pt" << FS << "x" << FS << "y" << FS << "z" << FS << "pterr" << FS << "hitprob" << LS;
  //*output_ << LS;

  TRandom3 die(1);

  struct Hit { 
    long int nevent;
    long int ntrack;
    double eta;
    double phi0;
    double z0;
    double pt; 
    double locx;
    double locy;
    float pterr;
    float hitprob;
    int deltas;
    int8_t cntref;
    int8_t zref;
    int8_t rhoref;
    int8_t phiref;
  } hit;

  TFile* outfile = new TFile(("tracks." + any2str(getpid()) + ".root").c_str(), "recreate");
  TTree* tree = new TTree("trackhits", "TTree containing hits of simulated tracks");
  tree->Branch("hit", &hit, "nevent/L:ntrack/L:eta/D:phi0/D:z0/D:pt/D:locx/D:locy/D:pterr/F:hitprob/F:deltas/I:cntref/B:zref/B:rhoref/B:phiref/B");
  //tree->SetAutoSave(100000);


  // build ordered maps
  for (long int i=0, totTracks = 0; i<numEvents; i++) {
    // new event
    for (long int j=0; j<numTracksEv; j++, totTracks++) {
      // randomly generate particle with eta = [-3,+3], phi = [0,2pi], Pt = [0.6,15]
      //double eta = die.Uniform(-3, 3);
      double eta = die.Uniform(-2, 2);
      double phi0 = die.Uniform(0.38, 1.634);
      double pt = 1./die.Uniform(0.01, 0.5);
      double z0 = die.Uniform(-20, 20);
      //double pt  = die.Uniform(0.6, 15); // everything linear just for debug
      //double eta = 1;
      //double phi0 = M_PI/3+0.2;
      //double pt = 0.6;

      // or pick from file
      //
      //
#ifdef FAKE_HITS
      if ((totTracks % 500) == 0) cout << "Track " << totTracks << " of " << numEvents*numTracksEv << std::endl;
#else
      if ((totTracks % 5000) == 0) cout << "Track " << totTracks << " of " << numEvents*numTracksEv << std::endl;
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

      //*output_ << "P" << FS << i << FS << j << FS << eta << FS << phi0 << FS << z0 << FS << pt << LS;

#ifdef VDEBUG
      double polylineZ[100], polylineR[100];
      double polylineX[100], polylineY[100];
      double polylineXec[100], polylineYec[100]; 
      int nPtsZR = 0, nPtsXY = 0, nPtsXYec = 0;

      double polyline3dX[100], polyline3dY[100], polyline3dZ[100];
      int nPts3d = 0;
#endif


      for (BarrelRadii::const_iterator rit = barrelModsByRadius_.begin(); rit != barrelModsByRadius_.end(); ++rit) {
        double r = rit->first;
        if (r >= 2*R) continue;
        double z = 1/B*acos(1-r*r/(2*R*R)) + z0;
        double x = dir*R*sin(B*(z-z0));
        double y = dir*R*(1-cos(B*(z-z0)));
        double xrot = x*cos(phi0) - y*sin(phi0);
        double yrot = x*sin(phi0) + y*cos(phi0);
        double phirot = atan2(y,x) + phi0;

#ifdef VDEBUG
        *output_ << i << FS << j << FS << " when r=" << r << ", z=" << z << ", x=" << xrot << ", y=" << yrot << LS;
        if (z > 0 && z < 1300 && y > 0) {
          polylineZ[nPtsZR] = z; polylineR[nPtsZR++] = yrot;
          polylineX[nPtsXY] = xrot; polylineY[nPtsXY++] = yrot;
          polyline3dX[nPts3d] = z; polyline3dY[nPts3d] = yrot; polyline3dZ[nPts3d++] = xrot;
        }
#endif

        const BarrelOctants& octants = rit->second;
        const BarrelModules& bmods = octants[getPointOctant(xrot, yrot, z)]; // jump to the octant of the point 

        for (BarrelModules::const_iterator mit = bmods.begin(); mit != bmods.end(); ++mit) {
          BarrelModule* mod = (*mit);
#ifdef VDEBUG 
          if (mod->getMeanPoint().Z() > 0 && mod->getLayer() == 2) {
            modXYZ.push_back(mod->getCorner(0).Z());
            modXYZ.push_back(mod->getCorner(0).Y());
            modXYZ.push_back(mod->getCorner(0).X());
            modXYZ.push_back(mod->getCorner(1).Z());
            modXYZ.push_back(mod->getCorner(1).Y());
            modXYZ.push_back(mod->getCorner(1).X());
            modXYZ.push_back(mod->getCorner(2).Z());
            modXYZ.push_back(mod->getCorner(2).Y());
            modXYZ.push_back(mod->getCorner(2).X());
            modXYZ.push_back(mod->getCorner(3).Z());
            modXYZ.push_back(mod->getCorner(3).Y());
            modXYZ.push_back(mod->getCorner(3).X());
            modXYZ.push_back(mod->getCorner(0).Z());
            modXYZ.push_back(mod->getCorner(0).Y());
            modXYZ.push_back(mod->getCorner(0).X());
          }
#endif
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
            hit = (Hit){ i, j, eta, phi0, z0, pt, locv.Y(), locv.Z(), pterr, hitprob, deltaStrips, posref.cnt, posref.z, posref.rho, posref.phi };
            tree->Fill();
            //*output_ << "T" << FS << i << FS << j << FS << eta << FS << phi0 << FS << z0 << FS << pt << FS << xrot << FS << yrot << FS << z << FS << pterr << FS << hitprob << LS;
            if (fabs(pt) >= HIGH_PT_THRESHOLD) break; // high pT particles never curve back inside the detector so after a layer/disk has been hit it makes no sense to look for more hits in modules in the same layer/disk
          }
        }
      }


      for (EndcapZs::const_iterator zit = endcapModsByZ_.begin(); zit != endcapModsByZ_.end(); ++zit) {
        double z = zit->first;
        double x = dir*R*sin(B*(z-z0));
        double y = dir*R*(1-cos(B*(z-z0)));
        double xrot = x*cos(phi0) - y*sin(phi0);
        double yrot = x*sin(phi0) + y*cos(phi0);
#ifdef VDEBUG
        *output_ << i << FS << j << FS << " when z=" << z << ", x=" << xrot << ", y=" << yrot << LS;

        if (z > 0 && y > 0) {
          polylineZ[nPtsZR] = z; polylineR[nPtsZR++] = yrot;
          polylineXec[nPtsXYec] = xrot; polylineYec[nPtsXYec++] = yrot;
        }
#endif

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
            hit = (Hit){ i, j, eta, phi0, z0, pt, locv.Y(), locv.X(), pterr, hitprob, deltaStrips, posref.cnt, posref.z, posref.rho, posref.phi }; // Y and X are inverted, because the wide side of the module lies along the Y axis
            tree->Fill();
            //*output_ << "T" << FS << i << FS << j << FS << eta << FS << phi0 << FS << z0 << FS << pt << FS << xrot << FS << yrot << FS << z << FS << pterr << FS << hitprob << LS;
            if (fabs(pt) >= HIGH_PT_THRESHOLD) break; // high pT particles never curve back inside the detector so after a layer/disk has been hit it makes no sense to look for more hits in modules in the same layer/disk
          }      
        }
      }
#ifdef VDEBUG

      yzDrawer.addModulesType(allMods_.begin(), allMods_.end(), Module::Barrel | Module::Endcap);
      yzDrawer.drawFrame<HistogramFrameStyle>(*RZCanvas);
      yzDrawer.drawModules<ContourStyle>(*RZCanvas);

      xyBarrelDrawer.addModulesType(allMods_.begin(), allMods_.end(), Module::Barrel);
      xyBarrelDrawer.drawFrame<HistogramFrameStyle>(*XYCanvas);
      xyBarrelDrawer.drawModules<ContourStyle>(*XYCanvas);

      xyEndcapDrawer.addModulesType(allMods_.begin(), allMods_.end(), Module::Endcap);
      xyEndcapDrawer.drawFrame<HistogramFrameStyle>(*XYCanvasEC);
      xyEndcapDrawer.drawModules<ContourStyle>(*XYCanvasEC);

      TPolyLine* polylineZR = new TPolyLine(nPtsZR, polylineZ, polylineR);
      TPolyLine* polylineXY = new TPolyLine(nPtsXY, polylineX, polylineY);
      TPolyLine* polylineXYec = new TPolyLine(nPtsXYec, polylineXec, polylineYec);
      RZCanvas->cd();
      polylineZR->SetLineColor(Palette::color(j)+2);
      polylineZR->SetLineWidth(2);
      polylineZR->Draw("P");
      XYCanvas->cd();
      polylineXY->SetLineColor(Palette::color(j)+2);
      polylineXY->SetLineWidth(2);
      polylineXY->Draw("P");
      XYCanvasEC->cd();
      polylineXYec->SetLineColor(Palette::color(j)+2);
      polylineXYec->SetLineWidth(2);
      polylineXYec->Draw("P");

      TPolyLine3D* mod3d = new TPolyLine3D(modXYZ.size()/3, &modXYZ[0]);
      Canvas3d->cd();
      mod3d->Draw();

      TPolyLine3D* polyline3d = new TPolyLine3D(nPts3d, polyline3dX, polyline3dY, polyline3dZ);
      polyline3d->SetLineColor(Palette::color(j)+2);
      polyline3d->SetLineWidth(2);
      polyline3d->Draw();

      Canvas3d->SaveAs("Canvas3d.root");

      gSystem->ProcessEvents();

      TImage* img = TImage::Create();
      img->FromPad(RZCanvas);
      img->WriteImage("RZCanvas.png");

      img->FromPad(XYCanvas);
      img->WriteImage("XYCanvas.png");

      img->FromPad(XYCanvasEC);
      img->WriteImage("XYCanvasEC.png");

      //img->FromPad(Canvas3d);
      //img->WriteImage("Canvas3d.png");


      delete mod3d;
      delete polylineZR;
      delete polylineXY;
      delete polylineXYec;
      delete polyline3d;
      delete img;
#endif
    }     
#ifdef FAKE_HITS
    for (std::vector<Module*>::iterator mit = allMods_.begin(); mit != allMods_.end(); ++mit) {
      Module* mod = (*mit);
      // loop over all the modules again
      // fake hits
      Polygon3d<4> poly;
      //poly << mod->getCorner(0) << mod->getCorner(1) << mod->getCorner(2) << mod->getCorner(3);  
      poly << mod->getCorner(0)-mod->getMeanPoint() << mod->getCorner(1)-mod->getMeanPoint() << mod->getCorner(2)-mod->getMeanPoint() << mod->getCorner(3)-mod->getMeanPoint();

      int numFake = die.Poisson(mod->getTriggerFrequencyFakePerEvent());
      for (int k = 0; k < numFake; k++) {
        XYZVector fakehit = poly.generateRandomPoint(&die);
        int deltaStrips = 1 + die.Integer(mod->getTriggerWindow()/2);
        float fakept = mod->getPtError()->stripsToP(deltaStrips);
        float pterr = mod->getPtError()->computeError(fakept);
        PosRef posref = mod->getPositionalReference();
        hit = (Hit){ i, -1, 0., 0., 0., fakept, fakehit.X(), fakehit.Y(), fakehit.Z(), pterr, 1.0, deltaStrips, posref.cnt, posref.z, posref.rho, posref.phi };
        tree->Fill();
        //*output_ << "F" << FS << i << FS << k << FS << 0. << FS << 0. << FS << 0. << FS << fakept << FS << fakehit.X() << FS << fakehit.Y() << FS << fakehit.Z() << FS << pterr << FS << 1.0 << LS;
      }
    }
#endif
  }

  outfile->Write();
  outfile->Close();
  //delete tree;
  delete outfile;
}



void TrackShooter::manualPolygonTestBench() {
  // generate poly
  TRandom3 die;
  while (true) {
    cout << std::endl << "---- New polygon (CTRL-C to exit) ----" << std::endl;
    double x, y, z;
    Polygon3d<4> poly;
    for (int i = 0; i < 4; i++) {
      std::cout << "Insert vertex " << i << " coords separated by spaces: ";
      std::cin >> x >> y >> z;
      XYZVector v(x, y, z);
      poly << v;
    }
    cout << "Polygon has vertices: " << std::endl;
    for (int i = 0; i < 4; i++) {
      XYZVector v = poly.getVertex(i);
      std::cout << "    " << v.X() << " " << v.Y() << " " << v.Z() << std::endl;
    }
    //cout << "Random points test: : " << std::endl;
    //for (int i = 0; i < 10; i++) {
    //  XYZVector randPoint = poly.generateRandomPoint(&die);
    //  bool isInside = poly.isPointInside(randPoint); // self-correctness check
    //  std::cout << "    Random point at: " << randPoint.X() << " " << randPoint.Y() << " " << randPoint.Z() << "... Is point inside? " << (isInside ? "yes" : "no") << std::endl; 
    //}
    while (true) {
      cout << "Manual point input (42 42 42 to exit): ";
      std::cin >> x >> y >> z;
      if (x == 42. && y == 42. && z == 42.) break;
      XYZVector myPoint(x, y, z);
      bool isInside = poly.isPointInside(myPoint);
      std::cout << "    User point at: " << myPoint.X() << " " << myPoint.Y() << " " << myPoint.Z() << "... Is point inside? " << (isInside ? "yes" : "no") << std::endl; 
    }
  }
}

void TrackShooter::moduleTestBench() {
  TRandom3 die;
  int a;
  while (true) {
    int index = die.Integer(allMods_.size());
    cout << std::endl << "---- Random module pick: #" << index << " ----" << std::endl;
    Module* m = allMods_[index];
    cout << "Type: " << (dynamic_cast<BarrelModule*>(m) ? "barrel" : "endcap") << std::endl;
    cout << "Ring: " << m->getRing() << " Layer/Disk: " << m->getLayer() << " Phi Idx: " << m->getPhiIndex() << std::endl;
    cout << "Center Point (x, y, z):     " << m->getMeanPoint().X() << " " << m->getMeanPoint().Y() << " " << m->getMeanPoint().Z() << std::endl;
    cout << "             (z, rho, phi): " << m->getMeanPoint().Z() << " " << m->getMeanPoint().Rho() << " " << m->getMeanPoint().Phi() << std::endl;
    Polygon3d<4> poly;
    poly << m->getCorner(0) << m->getCorner(1) << m->getCorner(2) << m->getCorner(3);  
    cout << "Corners (x, y, z):" << std::endl;
    for (int i = 0; i < 4; i++) {
      cout << "    " << m->getCorner(i).X() << " " << m->getCorner(i).Y() << " " << m->getCorner(i).Z() << std::endl;
    }
    cout << "Random points self-correctness test: " << std::endl;
    for (int i = 0; i < 10; i++) {
      XYZVector randPoint = poly.generateRandomPoint(&die);
      bool isInside = poly.isPointInside(randPoint); // self-correctness check
      std::cout << "    Random point at: " << randPoint.X() << " " << randPoint.Y() << " " << randPoint.Z() << "... Is point inside? " << (isInside ? "yes" : "no") << std::endl; 
    }
    while (true) {
      double x, y, z;
      cout << "Manual point input (42 42 42 to exit): ";
      std::cin >> x >> y >> z;
      if (x == 42. && y == 42. && z == 42.) break;
      XYZVector myPoint(x, y, z);
      bool isInside = poly.isPointInside(myPoint);
      std::cout << "    User point at: " << myPoint.X() << " " << myPoint.Y() << " " << myPoint.Z() << "... Is point inside? " << (isInside ? "yes" : "no") << std::endl; 
    }

    cout << "More?" << std::endl;
    cin >> a;
  }
}




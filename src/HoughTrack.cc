#include <HoughTrack.hh>


HoughTrack::~HoughTrack() {
}

double HoughTrack::calcPhi0(double x, double y, double pt) {
  /*
   * rotate hit to phi=0 (for ease of calculation)
  R^2 = (x - h)^2 + (y - k)^2
    R^2 = (0 - h)^2 + (0 - k)^2
    R = trajectory circle radius, depends on pT
    x,y = hit coordinates
    0,0 = particle origin
    find h1,k1. h2,k2
    2 solutions: for pt > 0 and < 0 
    phi1 = atan2(k1,h1)  phi2=atan2(k2,h2)
    rotate the angles back by phi and add a pi/2 correction factor (why?)
    */

  double phi = atan2(y,x);
  double xrot = x*cos(-phi) - y*sin(-phi);
  double yrot = x*sin(-phi) + y*cos(-phi);

  x = xrot;
  y = yrot;

  double R = fabs(pt)/(0.3*insur::magnetic_field) * 1e3; // particle curve radius

  double sqdelta = sqrt(4*R*R*pow(x,4) + 4*R*R*x*x*y*y - pow(x,6) - 2*pow(x,4)*y*y - x*x*pow(y,4));

  double k1 = (x*x*y + pow(y,3) + sqdelta)/(2*(x*x + y*y));
  double k2 = (x*x*y + pow(y,3) - sqdelta)/(2*(x*x + y*y));

  double h1 = sqrt(R*R - k1*k1);
  double h2 = sqrt(R*R - k2*k2);


  double phi1 = atan2(k1,h1) + phi - M_PI/2;
  double phi2 = atan2(k2,h2) + phi + M_PI/2;

  return pt > 0 ? phi1 : phi2;
}


double myatan2(double y, double x) {
  return 2*atan((sqrt(x*x + y*y)-x)/y);
}


double HoughTrack::calcTheta(double x, double y, double z, double z0, double pt) {
  
  double r = sqrt(x*x + y*y);

  double R = fabs(pt)/(0.3*insur::magnetic_field) * 1e3;

  double theta = myatan2(R*acos(1-r*r/(2*R*R)),(z-z0));
//  double theta = atan(R*acos(1-r*r/(2*R*R))/(z-z0));

  //double theta = atan2(r, z-z0);

  return theta;
}



double HoughTrack::rectangularSmear(double mean, double sigma, int nsteps, int step) {
  return mean - sigma + step*2*sigma/nsteps;
}


void HoughTrack::processHit(int evid, int hitid, double x, double y, double z, double pt, double ptError, double yResolution) { 
  const double sigmaZ0 = 70;
  const double sigmaZ = yResolution*sqrt(12)/2;
  double sigmaInvPt = 3*ptError*1/fabs(pt);
  //double invPt = die_.Gaus(1/pt, ptError); 
  double invPt = die_.Uniform(1/pt - sigmaInvPt, 1/pt + sigmaInvPt);
  int nSamplesPt = 2*sigmaInvPt/histo_.getWbins(H_K); 
  z = die_.Uniform(z-sigmaZ, z+sigmaZ);
  for (int k = 0; k < nSamplesPt; k++) {
    double invPtSample = rectangularSmear(invPt, sigmaInvPt, nSamplesPt, k);
    double phi0 = calcPhi0(x, y, 1/invPtSample);
    int nSamplesZ0 = 2*sigmaZ0/histo_.getWbins(H_Z0);
    for (int l = 0; l < nSamplesZ0; l++) {
      double z0Sample = rectangularSmear(0, sigmaZ0, nSamplesZ0, l);
      int nSamplesZ = 2*sigmaZ/histo_.getWbins(H_Z0);
      for (int m = 0; m < nSamplesZ; m++) {
        double zSample = rectangularSmear(z, sigmaZ, nSamplesZ, m);
        double theta = calcTheta(x, y, zSample, z0Sample, 1/invPtSample);
        //histo_.fill(seq<4>(invPtSample)(phi0)(z0Sample)(theta));
        histo_[invPtSample][phi0][z0Sample][theta] += SmartBin(1, evid, 1 << hitid);
      }
    }
  }
}

void HoughTrack::loadGeometryData(TFile* infile) {
  TTree* tree;
  ModuleData mdata;

  infile->GetObject("geomdata", tree);

  tree->SetBranchAddress("mdata", &mdata); 

  long int nentries = tree->GetEntriesFast();
  for (int i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    mods_[mdata.refcnt][mdata.refz][mdata.refrho][mdata.refphi] = mdata;
  }

}


void printTHisto(TH1* th, const char* drawopts = "", bool logx = false, bool logy = false) {
  TCanvas* canvas = new TCanvas("histo_canvas", "Histo Canvas", 1200, 1200 );
  canvas->cd();
  if (logx) canvas->SetLogx();
  if (logy) canvas->SetLogy();
  th->Draw(drawopts);
  TImage* img = TImage::Create();
  img->FromPad(canvas);
  img->WriteImage((th->GetName()+std::string(".png")).c_str());
  delete canvas;
}


//#define GENERATE_HIT_MAP

void HoughTrack::processTree(std::string filename, long int startev, long int howmany) {

  TracksP tracks;
  HitsP hits;

  ptError pterror;

  TFile* infile = new TFile(filename.c_str(), "read");
  if (infile->IsZombie()) {
    std::cerr << "Failed opening file \"" << filename << "\" for reading. Processing aborted." << std::endl;
    return;
  }
  TTree* tree;

  gROOT->ProcessLine("#include <vector>");

  loadGeometryData(infile);

  infile->GetObject("trackhits", tree);

  //TBranch *eventn, *trackn, *eta, *phi0, *z0, *pt, *glox, *gloy, *gloz, *locx, *locy, *pterr, *hitprob, *deltas, *cnt, *z, *rho, *phi;

  tree->SetBranchAddress("tracks.eventn", &tracks.eventn);//, &eventn);
  tree->SetBranchAddress("tracks.trackn", &tracks.trackn);
  tree->SetBranchAddress("tracks.eta", &tracks.eta);
  tree->SetBranchAddress("tracks.phi0", &tracks.phi0);
  tree->SetBranchAddress("tracks.z0", &tracks.z0);
  tree->SetBranchAddress("tracks.pt", &tracks.pt);
  tree->SetBranchAddress("tracks.nhits", &tracks.nhits);

  tree->SetBranchAddress("hits.glox", &hits.glox);
  tree->SetBranchAddress("hits.gloy", &hits.gloy);
  tree->SetBranchAddress("hits.gloz", &hits.gloz);
  tree->SetBranchAddress("hits.locx", &hits.locx);
  tree->SetBranchAddress("hits.locy", &hits.locy);
  tree->SetBranchAddress("hits.pterr", &hits.pterr);
  tree->SetBranchAddress("hits.hitprob", &hits.hitprob);
  tree->SetBranchAddress("hits.deltas", &hits.deltas);
  tree->SetBranchAddress("hits.cnt", &hits.cnt);
  tree->SetBranchAddress("hits.z", &hits.z);
  tree->SetBranchAddress("hits.rho", &hits.rho);
  tree->SetBranchAddress("hits.phi", &hits.phi);

  long int nevents = tree->GetEntriesFast();
#ifdef GENERATE_HIT_MAP
  //TH2I* trackHitsHisto = new TH2I("track_hits", "track hits;pt;eta", 100, -50, 50, 100, -2, 2);
  Histo<2, int> invPtEtaTrackCount(Seq<2,int>   (100)(100),
                                   Seq<2,double>(-.6)(-2.2),
                                   Seq<2,double> (.6)(2.2));
  Histo<2, double> invPtEtaAverageHits(Seq<2,int>   (100) (100),
                                       Seq<2,double>(-.6)(-2.2),
                                       Seq<2,double> (.6) (2.2));
#endif
  int minHits = 100, maxHits = 0;
  float minAvgHits = 100, maxAvgHits = 0;
  for (long int i = startev; i < howmany+startev && i < nevents; i++) {
    tree->GetEntry(i);
    /*if ((i-startev)%2 == 0)*/ std::cout << "Event " << i+1 << " of " << MIN(nevents,howmany+startev) << std::endl;
    for (unsigned int j = 0; j < tracks.trackn->size(); j++) {
      minHits = tracks.nhits->at(j) < minHits ? hits.cnt->size() /*tracks.nhits->at(j)*/ : minHits;
      maxHits = tracks.nhits->at(j) > maxHits ? hits.cnt->size() /*tracks.nhits->at(j)*/ : maxHits;
#ifndef GENERATE_HIT_MAP
      for (size_t k = 0; k < hits.cnt->size(); k++) {
        ModuleData& mdata = mods_[hits.cnt->at(k)][hits.z->at(k)][hits.rho->at(k)][hits.phi->at(k)];
        pterror.setPitch((mdata.pitchhi+mdata.pitchlo)/2);
        pterror.setStripLength(mdata.striplen);
        pterror.setZ(mdata.z);
        pterror.setR(mdata.rho);
        pterror.setDistance(mdata.stereo);
        pterror.setHeight(mdata.height);
        pterror.setInefficiencyType((ptError::InefficiencyType)mdata.inefftype);
        pterror.setModuleType(mdata.type);
        processHit(i, k, hits.glox->at(k), hits.gloy->at(k), hits.gloz->at(k), tracks.pt->at(j), hits.pterr->at(k), mdata.yres);
        //if (tracks.nhits->at(j) > 10) 
        //  std::cout << "  Mod z, rho, phi: " << mdata.z << "," << mdata.rho << "," << mdata.phi << " Hit invPt, pterr: " << 1/pt << "," << hits.pterr->at(k) << std::endl;
      }
#else
      double invpt = 1/tracks.pt->at(j), eta = tracks.eta->at(j);
      invPtEtaTrackCount.fill(seq<2>(invpt)(eta));
      double avg = invPtEtaAverageHits.get(seq<2>(invpt)(eta));
      invPtEtaAverageHits.fill(seq<2>(invpt)(eta), (hits.cnt->size() - avg)/invPtEtaTrackCount.get(seq<2>(invpt)(eta)));
#endif
    }
  }
  cout << "Transform done. Histo size: " << histo_.size() << " entries. " << histo_.size()*(sizeof(BinKey<4, uint16_t>) + sizeof(SmartBin))/1048576 << " MB (no overhead)" << std::endl;

#ifdef GENERATE_HIT_MAP
  std::ofstream hout("pt_eta_average_hits_3million.hst");
  invPtEtaAverageHits.serialize(hout);
  hout.close();
#else
  std::ifstream hin("pt_eta_average_hits_3million.hst", ios::in);
  Histo<2, double> invPtEtaAverageHits(hin); 
  hin.close();
#endif

  cout << "track with the fewest hits has: " << minHits << " hits" << std::endl;
  cout << "track with the most hits has: " << maxHits << " hits" << std::endl;

  int minCell = 100, maxCell = 0;
  int maxStacked = 0;
  TH1I* cellLoadHisto = new TH1I("cell_load", "cell load over theoretical number of hits;C/N", 50, 0, 2);
  for (HistoType::const_iterator it = histo_.begin(); it != histo_.end(); ++it) {
    if (it->first.overflow() || it->first.underflow()) continue;
    double avgload = invPtEtaAverageHits.get( seq<2>(it->first.at(0))(-log(tan(it->first.at(3)/2))) );
    if (avgload == 0)
      cout << "unmapped value at invpt,eta,theta: " << it->first.at(0) << "," << -log(tan(it->first.at(3)/2)) << "," << it->first.at(3) << std::endl;
    minAvgHits = MIN(minAvgHits, avgload);
    maxAvgHits = MAX(maxAvgHits, avgload);
    minCell = MIN(minCell, (int)it->second);
    maxCell = MAX(maxCell, (int)it->second);
    maxStacked = MAX(maxStacked, it->second.stacked);
    cellLoadHisto->Fill( ((double)it->second) / floor(avgload) );
  }
  
  cout << "minimum cell value for track streak: " << minCell << " hits" << std::endl;
  cout << "maximum cell value for track streak: " << maxCell << " hits" << std::endl;
  cout << "minimum average hits from the hit map: " << minAvgHits << " hits" << std::endl;
  cout << "maximum average hits from the hit map: " << maxAvgHits << " hits" << std::endl;
  cout << "maximum stacked events: " << maxStacked << " events" << std::endl;

  TH2I* thisto1 = new TH2I();
  TH2I* thisto2 = new TH2I();
  TH2I* thisto3 = new TH2I();
  TH2I* thisto4 = new TH2I();
  TH2D* thistom = new TH2D();

  thisto1->SetNameTitle("invpt_phi0", "histo1;invpt;phi0");
  toTH2(histo_, *thisto1, H_K, H_PHI0); 

  thisto2->SetNameTitle("z0_theta", "histo2;z0;theta");
  toTH2(histo_, *thisto2, H_Z0, H_THETA);

//  thisto3->SetNameTitle("invpt_theta", "histo3;invpt;theta");
//  toTH2(histo_, *thisto3, 0, 3);

//  thisto4->SetNameTitle("invpt_z0", "histo4;invpt;z0");
//  toTH2(histo_, *thisto4, 0, 2);

  gStyle->SetOptStat(0);

#ifdef GENERATE_HIT_MAP
  thistom->SetNameTitle("invpt_eta_average_hits", "avg hits;invpt;eta");
  toTH2(invPtEtaAverageHits, *thistom);
  printTHisto(thistom, "colz");
#endif
  tree->ResetBranchAddresses();


  printTHisto(thisto1, "colz");
  printTHisto(thisto2, "colz");
  printTHisto(thisto3, "colz");
  printTHisto(thisto4, "colz");


  printTHisto(cellLoadHisto, "", false, true);


  delete thisto1;
  delete thisto2;
  delete thisto3;
  delete thisto4;
  delete thistom;
}



int main(int argc, char* argv[]) {

  HoughTrack ht;
  ht.processTree(argv[1], str2any<long int>(argv[2]), str2any<long int>(argv[3]));
  
  return 0;
}




#include <TFile.h>
#include <TCanvas.h>
#include <string>
#include <TList.h>
#include <TKey.h>
#include <TObject.h>
#include <TProfile.h>
#include <TAxis.h>
#include <string>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>


TFile* inputFile = NULL;
TCanvas* aCanvas = NULL;
TCanvas* anotherCanvas = NULL;

std::map<double, TProfile*> ptProfiles;

std::ostream* osp = &std::cout;

#define VNAME(x) #x
#ifdef DEBUG
#define VDUMP(x) std::cerr << #x << " = " << x << std::endl
#else
#define VDUMP(x)
#endif

void browse(std::string fileNameString) {
  const char* fileName = fileNameString.c_str();
  inputFile = TFile::Open(fileName);
  if (!inputFile) {
    std::cerr << "File '" << fileNameString << "' does not exist" << std::endl;
    return;
  }
  // new TBrowser;
  TKey* lll = (TKey*)inputFile->GetListOfKeys()->At(0);
  if (lll) {
    TObject* myObject = lll->ReadObj();
    if (myObject) {
      std::string aClass=myObject->ClassName();
      if (aClass=="TCanvas") {
        aCanvas = (TCanvas*) myObject;
        TList* aList=aCanvas->GetListOfPrimitives();
        for (int i=0; i<aList->GetSize(); ++i) {
          std::string aClassName = aList->At(i)->ClassName();
          if (aClassName=="TProfile") {
            TProfile* myProfile = (TProfile*)aList->At(i);
            // std::cerr << "TProfile:" << myProfile->GetName() << std::endl;
            double aMomentum;
            if (sscanf(myProfile->GetName(), "pt_vs_eta%lftracker_profile", &aMomentum)==1) {
              // std::cerr << "aMomentum=" << aMomentum << std::endl;
              ptProfiles[aMomentum]=(TProfile*) myProfile->Clone();
            }
          }
        }
      }
    }
  }
  //inputFile->Close();

  std::string plotOpt = "";
  aCanvas = new TCanvas();
  aCanvas->cd();
  aCanvas->SetLogy();
  //for (std::map<int, TProfile*>::iterator it = ptProfiles.begin(); it!=ptProfiles.end(); ++it) {
  for (auto it : ptProfiles) {
    // std::cerr << "aMomentum=" << it.first << std::endl;
    it.second->Draw(plotOpt.c_str());
    plotOpt="same";
  }
  aCanvas->SaveAs(Form("%s/ptResolution.png", getenv("TKG_DELPHES_TARGET_DIRECTORY")));
  aCanvas->SaveAs(Form("%s/ptResolution.root", getenv("TKG_DELPHES_TARGET_DIRECTORY")));
}

void printEtaRange(double lowEta, double highEta) {
  (*osp) << Form("   (abs(eta) >= %.04f && abs(eta) < %.04f) * ",
                    lowEta, highEta);
}

void printPtRange(double lowPt, double highPt) {
  if (highPt>lowPt) {
    (*osp) << Form("(pt >= %.04f && pt < %.04f) * ",
                      lowPt, highPt);
  } else {
    (*osp) << Form("(pt >= %.04f) * ",
                      lowPt);
  }
}

void printNewLine(bool last) {
  if (!last) (*osp) << " + \\" << std::endl;
  else (*osp) << std::endl << "}" << std::endl;;
}

void printResolutionScale(double resolution1,
                          double pt1,
                          double resolution2,
                          double pt2) {
  if (pt1==pt2) {
    (*osp) << Form("(%.8f)", resolution1);
  } else { 
    // Linear scaling!
    (*osp) << Form("(%f + "                 // resolution1
                      "(pt-%f)"                // pt1
                      "* %f)",                 // (res2-res1) / (pt2-pt1)
                      resolution1, pt1,
                      (resolution2-resolution1)/(pt2-pt1) );
  }
}

void printResolutionStandardWorsen(double myResolution, double myPt) {
  (*osp) << Form("(%f*pt/%f)", myResolution, myPt);
}


bool delphize(std::string rootFileName, double etaSlice=.2) {
  browse(rootFileName.c_str());
  auto firstProfile = ptProfiles.begin();
  if (firstProfile==ptProfiles.end()) {
    std::cerr << "Error: the collection of profiles is empty" << std::endl;
    return false;
  }
  TProfile* exampleProfile = firstProfile->second;
  double maxEta = exampleProfile->GetXaxis()->GetXmax();
  VDUMP(maxEta);
  double minEta = exampleProfile->GetXaxis()->GetXmin();
  VDUMP(minEta);
  double originalEtaStep = exampleProfile->GetXaxis()->GetBinWidth(1);
  VDUMP(originalEtaStep);
  int rebinScale = etaSlice / originalEtaStep;
  VDUMP(rebinScale);
  for (auto it : ptProfiles) it.second->Rebin(rebinScale);
  int newBins = exampleProfile->GetXaxis()->GetNbins();
  VDUMP(newBins);

  // OK, for each eta slice (rebinned bins) we will compute
  // an interpolation formula. Keep your finger crossed...
  
  // Looping on eta range first:
  double lowEta, highEta;
  double myResolution, nextResolution;
  double myPt, nextPt;
  (*osp) << "# set ResolutionFormula for layout " << getenv("SELECTED_LAYOUT_NAME") << std::endl;
  (*osp) << "set ResolutionFormula { ";
  for (int iBin=1; iBin<=newBins; ++iBin) {
    // Looping on pT
    lowEta = exampleProfile->GetXaxis()->GetBinLowEdge(iBin);
    highEta = exampleProfile->GetXaxis()->GetBinUpEdge(iBin);
    for (auto it = ptProfiles.begin(); it!=ptProfiles.end(); ++it) {
      auto nextItem = (++it)--;
      TProfile* myProfile = it->second;
      myPt = it->first;
      myResolution = myProfile->GetBinContent(iBin)/100.;
      printEtaRange(lowEta, highEta);
      if (it==ptProfiles.begin()) {
        // From here down the resolution will be always the same
        printPtRange(0, myPt);
        printResolutionScale(myResolution, myPt, myResolution, myPt);
        printNewLine(false);
        printEtaRange(lowEta, highEta);
      }
      if (nextItem==ptProfiles.end()) {
        // From here up the resolution will just scale with pT
        printPtRange(myPt, -1);
        printResolutionStandardWorsen(myResolution, myPt);
        printNewLine(iBin==newBins); // true if it is the last
      } else {
        TProfile* nextProfile = nextItem->second;
        nextPt = nextItem->first;
        nextResolution = nextProfile->GetBinContent(iBin)/100.;
        printPtRange(myPt, nextPt);
        printResolutionScale(myResolution, myPt, nextResolution, nextPt);
        printNewLine(false);
      }
    }
  }

  return true;
}

int main(int argc, char* argv[]) {
  if (argc<2) {
    std::cerr << Form("Syntax: %s filename.root [delphesConfig.txt]", argv[0]) << std::endl;
    return -1;
  }
  if (argc==3) {
    std::ofstream* fout = new std::ofstream();
    fout->open(argv[2]);
    if (fout->is_open()) {
      osp = fout;
    } else {
      std::cerr << "ERROR: could not open " << argv[2] << " writing config to the screen" << std::endl;
    }
  }
  std::string rootFileName=argv[1];
  double etaSlice=.5;
  if (delphize(rootFileName, etaSlice))  {
    return 0;
  } else return -1;
}

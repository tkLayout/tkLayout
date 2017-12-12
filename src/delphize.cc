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
#include <time.h>
#include <iomanip>
#include <boost/program_options.hpp>

std::ostream* osp = &std::cout;

#define VNAME(x) #x
#ifdef DEBUG
#define VDUMP(x) std::cerr << #x << " = " << x << std::endl
#else
#define VDUMP(x)
#endif

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

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

std::map<double, TProfile*> ptProfiles;

void addProfile(TObject* anObject) {
  char buffer[1024];
  char buffer2[1024];

  std::string aClassName = anObject->ClassName();
  if (aClassName=="TProfile") {
    TProfile* myProfile = (TProfile*)anObject;;
    // std::cerr << "TProfile: " << myProfile->GetName() << std::endl;
    double aMomentum;
    if (sscanf(myProfile->GetName(), "pt_vs_eta%lftracker_profile", &aMomentum)==1) {
      std::cerr << "Momentum [GeV]: " << aMomentum << std::endl;
      ptProfiles[aMomentum]=(TProfile*) myProfile->Clone();
    }
    else if (sscanf(myProfile->GetName(), "Total_pt_vs_eta%lftracker_profile", &aMomentum)==1) {
      std::cerr << "Momentum [GeV]: " << aMomentum << std::endl;
      ptProfiles[aMomentum]=(TProfile*) myProfile->Clone();
    }
    else if (sscanf(myProfile->GetName(), "Resolution_(%[^)]).ptres_tracker_MS_Pt.pt_vs_eta%lftracker_profile", buffer, &aMomentum)==2) {
      std::cerr << "Momentum [GeV]: " << aMomentum << " buffer: " << buffer << " (" << myProfile->GetName() << ")" << std::endl;
      ptProfiles[aMomentum]=(TProfile*) myProfile->Clone();
    }
  }
}

int main(int argc, char* argv[]) {

  //
  // Define variables
  std::string usage("Usage: ");
  usage += argv[0];
  usage += " [options]";

  std::string rootFileName;
  std::string layoutName;
  TFile* inputFile       = nullptr;
  std::ofstream* outFile = nullptr;

  double      etaSlice = 0.0;
  std::string author   = "";

  //
  // Program options
  boost::program_options::options_description help("Program options");
  help.add_options()
    ("help,h"        , "Display help")
    ("input-file,i"  , boost::program_options::value<std::string>()                                        , "Specify name of input root file (from tkLayout), containing canvas with pT profiles")
    ("layout-name,n" , boost::program_options::value<std::string>()                                        , "Specify given layout name")
    ("author,a"      , boost::program_options::value<std::string>(&author)->default_value("Unknown author"), "Specify author of DELPHES file (optional: default = Unknown author)")
    ("eta-slicing,s" , boost::program_options::value<double>(&etaSlice)->default_value(0.2)                , "Specify eta slicing <0.0?; 1.0>, e.g. 0.1, 0.2, ... (optional: default = 0.2)")
    ;

  // Read user input
  boost::program_options::variables_map varMap;
  try {

    // Parse user defined options
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(help).run(), varMap);
    boost::program_options::notify(varMap);

    // Check
    if      (etaSlice <= 0.0 || etaSlice > 1.0)                     throw boost::program_options::invalid_option_value("eta-slicing");
    else if (!varMap.count("input-file")  && !varMap.count("help")) throw boost::program_options::error("Forgot to define input root file???");
    else if (!varMap.count("layout-name") && !varMap.count("help")) throw boost::program_options::error("Forgot to define layout name???");
    else {

      // Write help
      if (varMap.count("help")) {
        std::cout << usage << std::endl << help << std::endl;
        return -1;
      }

      // Check that the provided ROOT input file corresponds to an existing file.
      rootFileName = varMap["input-file"].as<std::string>();
      inputFile = new TFile(rootFileName.c_str(),"READ");
      if (!inputFile) {
        std::cerr << "File '" << rootFileName << "' does not exist" << std::endl;
        return 0;
      }

      // Delphes output file
      outFile = new std::ofstream();
      layoutName = varMap["layout-name"].as<std::string>();
      std::string delphesName  = layoutName+"_Delphes.conf";
      outFile->open(delphesName);
      if (outFile->is_open()) {
        osp = outFile;
      } else {
        std::cerr << "ERROR: could not open " << outFile << " -> writing config to the screen" << std::endl;
        return -1;
      }
    }
  }catch(boost::program_options::error& e) {

    // Display error type
    std::cerr << "\nERROR: " << e.what() << std::endl << std::endl;
    std::cout << usage                   << std::endl << help << std::endl;
    return -1;
  }


  // Read pt profiles > browse;
  // This works if the user passes a TCanvas root file
  TKey* key = (TKey*)inputFile->GetListOfKeys()->At(0);
  if (key) {
    // Object
    TObject* myObject = key->ReadObj();
    if (myObject) {
      // Canvas
      std::string aClass=myObject->ClassName();
      if (aClass=="TCanvas") {

        TCanvas* aCanvas = (TCanvas*) myObject;
        TList* aList=aCanvas->GetListOfPrimitives();
        for (int i=0; i<aList->GetSize(); ++i) {
	  // Add the graph (if it's an appropriate TProfile)
	  addProfile(aList->At(i));
        }
      }
    }
  }
  // Now try to get the correct file from the summary.root file (instead)
  TIter nextItem(inputFile->GetListOfKeys());
  while ((key = (TKey*)nextItem())) {
    addProfile(key->ReadObj());
  }
  //inputFile->Close();

  // Check that non-zero profile map
  auto itPtProfiles = ptProfiles.begin();
  
  if (itPtProfiles==ptProfiles.end()) {
    std::cerr << "Error: the collection of profiles is empty" << std::endl;
    return false;
  }
  
  // First profile -> get eta, rescale bins, ...
  TProfile* exampleProfile = itPtProfiles->second;

  double maxEta = exampleProfile->GetXaxis()->GetXmax();
  VDUMP(maxEta);
  double minEta = exampleProfile->GetXaxis()->GetXmin();
  VDUMP(minEta);
  double originalEtaStep = exampleProfile->GetXaxis()->GetBinWidth(1);
  VDUMP(originalEtaStep);
  int rebinScale = etaSlice / originalEtaStep;
  VDUMP(rebinScale);
  for (auto it = ptProfiles.begin(); it!=ptProfiles.end(); ++it) it->second->Rebin(rebinScale);
  int newBins = exampleProfile->GetXaxis()->GetNbins();
  VDUMP(newBins);

  // OK, for each eta slice (rebinned bins) we will compute
  // an interpolation formula. Keep your finger crossed...
  
  // Looping on eta range first:
  double lowEta, highEta;
  double myResolution, nextResolution;
  double myPt, nextPt;
  (*osp) << "#" << std::endl;
  (*osp) << "# Automatically generated tracker resolution formula for layout: " << layoutName << std::endl;
  (*osp) << "#" << std::endl;
  (*osp) << "#  By " << author << " on: " << currentDateTime() << std::endl;
  (*osp) << "#" << std::endl;
  (*osp) << "set ResolutionFormula { ";

  for (int iBin=1; iBin<=newBins; ++iBin) {

    // Get eta
    lowEta  = exampleProfile->GetXaxis()->GetBinLowEdge(iBin);
    highEta = exampleProfile->GetXaxis()->GetBinUpEdge(iBin);

    // Loop over pT
    for (auto it = ptProfiles.begin(); it!=ptProfiles.end(); ++it) {

      auto nextItem = (++it)--;

      double aMomentum = it->first;
      if (iBin==1) std::cout << "Using momentum [GeV]: " << std::setw(6) << aMomentum << std::endl;

      TProfile* myProfile = it->second;
      myPt         = it->first;
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
      }
      else {
        TProfile* nextProfile = nextItem->second;
        nextPt = nextItem->first;
        nextResolution = nextProfile->GetBinContent(iBin)/100.;
        printPtRange(myPt, nextPt);
        printResolutionScale(myResolution, myPt, nextResolution, nextPt);
        printNewLine(false);
      }
    }
  }

  // Close file
  inputFile->Close();
  delete inputFile;
  outFile->close();
  delete outFile;

  return true;
}

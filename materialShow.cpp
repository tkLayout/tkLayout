// Copyright Vladimir Prus 2002-2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/* Shows how to use both command line and config file. */

#include <TCanvas.h>
#include <TFile.h>
#include <TFrame.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TAxis.h>
#include <TH1.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <TLegend.h>

#include <Palette.h>


#include <mainConfigHandler.h>
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

mainConfigHandler mainConfiguration;
TCanvas* plotCanvas = NULL;
TLegend* myLegend;
TH1D* lastPixel;
double range;

TFile* openFile(std::string dirName, std::string fileName) {
  std::string completeFileName = dirName + "/" + fileName;
  return new TFile(completeFileName.c_str());
}

TCanvas* getCanvasFile(string dirName, string canvasName) {
  TCanvas* result = NULL;
  string fileName = canvasName + ".root";
  TFile* myFile = openFile(dirName, fileName);
  if (myFile) {
    result = (TCanvas*)myFile->FindObjectAny(canvasName.c_str());
  }
  return result;
}

void plotSumMaterials(TCanvas* myTkMat, TCanvas* myPxMat, string myPlot, int iPad, int kColor, string layoutName) {
  TH1D* myGraph = (TH1D*)myTkMat->GetPad(iPad)->FindObject(myPlot.c_str());
  if (!myGraph) {
    cerr << myPlot << " not found, will not continue" << endl;
    return;
  }
  TH1D* tempGraph = new TH1D(*myGraph);
  TH1D* myGraph2 = NULL;

  if (myPxMat) {
    myGraph2 = (TH1D*)myPxMat->GetPad(iPad)->FindObject(myPlot.c_str());
    if (myGraph2) {
      tempGraph->Add(myGraph2);
      myGraph2->SetFillColor(kBlack);
      myGraph2->SetFillStyle(3001);
      lastPixel = myGraph2;
    }
  }

  plotCanvas->cd();
  tempGraph->SetFillColor(kColor);
  tempGraph->Draw("same"); // total
  myLegend->AddEntry(tempGraph, layoutName.c_str(), "F");

  if (myGraph2) 
    myGraph2->Draw("same");  // pixel

}

void plotMaterial(string dirName, int kColor, bool usePixel, string layoutName) {
  TCanvas* myTkMat=NULL;
  TCanvas* myPxMat=NULL;
  myTkMat = getCanvasFile(dirName, "matFull000");
  
  if (usePixel)
    myPxMat = getCanvasFile(dirName, "matFull001");

  //cerr << "myTkMat = " << myTkMat << endl; // debug
  //cerr << "myPxMat = " << myPxMat << endl; // debug

  plotSumMaterials(myTkMat, myPxMat, "rfullvolume", 1, kColor, layoutName);
}

void prepareCanvas() {
  plotCanvas = new TCanvas();
  plotCanvas->SetFillColor(kWhite);
  plotCanvas->cd();
  plotCanvas->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  TH1D* ranger = new TH1D("ranger","", 100, 0, 3);
  ranger->SetMaximum(range);
  TAxis* myAxis;
  myAxis = ranger->GetXaxis();
  myAxis->SetTitle("#eta");
  myAxis = ranger->GetYaxis();
  myAxis->SetTitle("Material amount [x/X_{0}]");
  myLegend = new TLegend(0.15, 0.75, .35, .99);
  ranger->Draw();
}

void printCanvas(string fileName) {
  string pngFileName=fileName+".png";
  string rootFileName=fileName+".root";
  myLegend->Draw();
  plotCanvas->Print(pngFileName.c_str());
  plotCanvas->Print(rootFileName.c_str());
}


// A helper function to simplify the main part.
template<class T> ostream& operator<<(ostream& os, const vector<T>& v) {
    copy(v.begin(), v.end(), ostream_iterator<T>(cout, " "));
    return os;
}

void printUsage(po::options_description& visible, string progName) {
  cout << progName << " creates a plot (png and root) with the materials of the selected layouts" << endl;
  cout << "The layouts are taken from the default output directory" << endl;
  cout << "(currently " << mainConfiguration.getLayoutDirectory() << ")" << endl;
  cout << endl;
  cout << "Syntax: " << endl;
  cout << "  " << progName << " [options] layout1" << endl;
  cout << "  " << progName << " [options] layout1 layout2" << endl;
  cout << "" << endl;
  cout << "  Shows the material of one layout or compares two layouts" << endl << endl;
  cout << visible << "\n";
}

int main(int ac, char* av[]) {
  try {
    range = 2;
    string outputName = "output";
    vector<string> fileV;
    lastPixel = NULL;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic");
    generic.add_options()
      ("help", "produce help message")
      ;

    // Declare optional options
    po::options_description optional("Material plot");
    optional.add_options()
      ("trackerOnly,t", "takes only into account the outer tracker")
      ("maxRange,r", po::value<double>(&range), "sets the vertical range for the material plot")
      ("output,o", po::value<string>(&outputName), "output file name (.png will be added)")
      ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", po::value< vector<string> >(&fileV), "input file")
      ;

    // Create the command line option set
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(optional).add(hidden);

    // Create the list to be shown
    po::options_description visible("Options");
    visible.add(generic).add(optional);

    // Input file is the positional option
    po::positional_options_description p;
    p.add("input-file", -1);

    // Actually read the options from command line
    po::variables_map vm;
    store(po::command_line_parser(ac, av).
	  options(cmdline_options).positional(p).run(), vm);
    notify(vm);

    // List the options
    if (vm.count("help")) {
      printUsage(visible, av[0]);
      return 0;
    }
    
    bool usePixel = (vm.count("trackerOnly")==0);
    int myColor;

    // Look into the input-file list
    if (vm.count("input-file")) {
      string dirName;
      prepareCanvas();
      Palette::skipColors(100);
      Palette::prepare(fileV.size(),210,0.75,0.8);
      for (unsigned int iMat=0 ; iMat<fileV.size(); ++iMat) {
	dirName = (mainConfiguration.getLayoutDirectory()
		   + "/" + fileV.at(iMat));
	cerr << dirName << endl; // debug
	myColor = Palette::color(iMat);
	plotMaterial(dirName, myColor, usePixel, fileV.at(iMat));
      }
      if (lastPixel) {
	myLegend->AddEntry(lastPixel, "pixel", "F");
      }
      //Palette::prepare(10, 210, 0.9, 0.4);
      printCanvas(outputName);
    } else {
      cout << "Error: at least a layout name should be given" << endl << endl;

      printUsage(visible, av[0]);
      return -1;
    }

  } catch(exception& e) {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}
  

#include "TStyle.h"
#include "stylePlot.h"

#include "configparser.hh"
#include "tracker.hh"
#include <string>
#include <boost/filesystem/operations.hpp>
#include <iostream>
#include <fstream>
using namespace boost::filesystem;

using namespace std;


// This function only takes a geometry configuration file and creates
// the tracker object with the rules described there. No module "dressing"
// is done here.
Tracker* createGeometry(string configFileName) {
  configParser myParser;

  return myParser.parseFile(configFileName);

}

// This function only takes a geometry configuration files and creates
// the tracker object with the rules described there. No module "dressing"
// is done here.
Tracker* createDressedGeometry(string configFileName, string dressFileName) {
  configParser myParser;
  Tracker* myTracker;

  if ( (myTracker=myParser.parseFile(configFileName)) ) {
    myParser.dressTracker(myTracker, dressFileName);
  } 

  return myTracker;
}

// This function creates a geometry package from a given configuration file
// The output directory will be res/TrackerName/
// TrackerName is taken form the configuration file itself
Tracker* createGeometryPackage(string configFileName) {
  Tracker* myTracker = createGeometry(configFileName);
  string dirName;
  string trackerName;
  string destConfigFileName;

  if (myTracker) {
    myTracker->createGeometry(true);
    trackerName=myTracker->getName();
    dirName = "./res/"+myTracker->getName();
    // TODO: treat the error properly
    create_directory(dirName);
    // This creates the layout and description files
    myTracker->createPackageLayout(dirName);
    std::string descriptionFile = dirName + "/desc.txt";
    if (!exists(descriptionFile)) {
      std::ofstream myfile;
      myfile.open(descriptionFile.c_str());
      myfile << myTracker->getComment() << std::endl;
      myfile.close();
    }
    destConfigFileName =  dirName+"/geometry.cfg";
    std::cout << "Copying " << configFileName <<" into " << destConfigFileName <<endl;
    remove(destConfigFileName);
    copy_file(configFileName, destConfigFileName);

    std::string settingsFile = dirName + "/settings.cfg";
    if (!exists(settingsFile)) {
      std::ofstream myfile;
      myfile.open(settingsFile.c_str());
      myfile << std::endl;
      myfile.close();
    }
    
  }

  return myTracker;
}

std::string extractFileName(const std::string& full) {
  std::string::size_type idx = full.find_last_of("/");
  if (idx != std::string::npos)
    return full.substr(idx+1);
  else return full;
}

Tracker* analyzeGeometryPackage(string configFileName, string dressFileName) {
  Tracker* myTracker = createDressedGeometry(configFileName, dressFileName);

  std::string myDirectory;
  std::string destConfigFile;
  std::string destDressFile;
  std::string endcapModuleCoordinatesFile = "endcapCoordinates.csv";

  if (myTracker) {
    myTracker->createGeometry(true);

    // Optical transmission
    myTracker->computeBandwidth();

    // Analysis
    //myTracker->analyze(2000, Layer::YZSection);
    myTracker->analyze(2000);

    // Summary and save
    //myTracker->writeSummary(true, extractFileName(configFileName), extractFileName(dressFileName), endcapModuleCoordinatesFile);
    myTracker->writeSummary(true, extractFileName(configFileName), extractFileName(dressFileName), "html", "", endcapModuleCoordinatesFile);
    myTracker->save();

    myDirectory = myTracker->getActiveDirectory();
    destConfigFile = myDirectory + "/" + extractFileName(configFileName);
    destDressFile = myDirectory + "/" + extractFileName(dressFileName);

    // File(s) with module positions
    endcapModuleCoordinatesFile = myDirectory + "/" + extractFileName(endcapModuleCoordinatesFile);
    ofstream endcapFile;
    endcapFile.open(endcapModuleCoordinatesFile.c_str());
    myTracker->printEndcapModuleRPhiZ(endcapFile);
    endcapFile.close();
    //myTracker->printBarrelModuleZ(std::cout);

    remove(destConfigFile);
    remove(destDressFile);
    copy_file(configFileName, destConfigFile);
    copy_file(dressFileName, destDressFile);
  }

  return myTracker;
}

void displaySyntax(char* programName) {
  cout << "Syntax: " << endl;
  cout << programName << " geometry.cfg" << endl;
  cout << "\t Creates a geometry package according to geometry definition config file 'geometry.cfg'" << endl;
  cout << programName << " geometry.cfg moduleTypes.cfg" << endl;
  cout << "\t Creates a summary of a complete geometry defined by config files 'geometry.cfg' and 'moduleTypes.cfg'" << endl;
}

int main(int argc, char* argv[]) {
  //configParser myParser;
  //Tracker* myTracker;

  if (argc==3) {
    string geometryConfig(argv[1]);
    string dressingConfig(argv[2]);
    analyzeGeometryPackage(geometryConfig, dressingConfig);
  } else if (argc==2) {
    string geometryConfig(argv[1]);
    createGeometryPackage(geometryConfig);
  } else {
    displaySyntax(argv[0]);
  }

  return 0;
}

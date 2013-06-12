#ifndef _CONFIGPARSER_HH_
#define _CONFIGPARSER_HH_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <exception>
#include "tracker.hh"
#include "global_funcs.h"
#include <StopWatch.h>

#define NMAXCHARS 100 // Max number of characters for parameterName[index]

using namespace std;

class configParser {
public:
  // Constructor and destructor
  configParser();
  ~configParser();

  // Parse the geometry config file and creates a Tracker object
  Tracker* parseFile(string configFileName, string typesFileName = "");
  Tracker* parsePixelsFromFile(string configFileName);

  // Parse the module type config file and "dresses" a Tracker object
  bool dressTracker(Tracker* aTracker, string configFileName);
  bool dressPixels(Tracker* aTracker, string configFileName);
  // Parse an irradiation map so that modules can be assigned radiation levels
  bool irradiateTracker(Tracker* aTracker, string irrFileName);
  bool flukaGridSteps(Tracker* aTracker, string irrFileName);

  // Extract the user-defined support structures from the geometry config file
  std::list<std::pair<int, double> >* parseSupportsFromFile(std::string fileName);

private:

  // Temporary streams and objects
  ifstream rawConfigFile_;
  stringstream configFile_;
  Tracker* myTracker_;

  std::map<std::string, TiltedBarrelSpecs> tiltedBarrelSpecs_;

  // Generic parsing functions
  string getTill(istream &inStream, char delimiter, bool singleWord, bool allowNothing /* = false */);
  bool parseParameter(string& name, string& value, istream &inStream);
  
  // Parsing functions for the tracker building
  bool parseTracker(string myName, istream &inStream);
  bool parseBarrel(string myName, istream &inStream);
  bool parseEndcap(string myName, istream &inStream);
  bool parsePixels(string myName, istream &inStream);
  bool parseObjectType(string myType);

  // Parsing functions for the tracker dressing
  bool parseDressType(string myType);
  bool breakParameterName(string& parameterName, int& ringIndex, int& diskIndex);
  bool breakParameterName(string& parameterName, string& stringIndex);
  bool parseBarrelType(string myName, istream &inStream);
  bool parseEndcapType(string myName, istream &inStream);
  bool parsePixelType(istream& inStream);
  bool parseAnyType(string myName, istream &inStream);
  bool parseOutput(istream &inStream);

  // Parsing function for inactive surfactes
  bool parseSupportParameters(std::istream& inStream, std::list<std::pair<int, double> >& list);

  bool peekTypes(string typesFileName); // goes through the types file to get the dsDistances (sensor spacing) to use for the geometry before tracker is dressed (or even created!)

  bool parseTilted(const std::string& fileName, const std::string& barrelName);

  std::map<std::string, std::map<int, double> > geometryDsDistance_; // CUIDADO: not pretty but it will do the job
  std::map<std::string, std::map<std::pair<int, int>, double> > geometryDsDistanceSecond_;
};

// Definition of the parsing exception
class parsingException: public exception
{
  virtual const char* what() const throw()
  {
    return "Error parsing the config file.";
  }
  
};

#endif

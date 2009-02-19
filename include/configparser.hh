#ifndef _CONFIGPARSER_HH_
#define _CONFIGPARSER_HH_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <exception>
#include "tracker.hh"

#define NMAXCHARS 100 // Max number of characters for parameterName[index]

using namespace std;

class configParser {
public:
  // Constructor and destructor
  configParser();
  ~configParser();

  // Parse the geometry config file and creates a Tracker object
  Tracker* parseFile(string fileName);

  // Parse the module type config file and "dresses" a Tracker object
  bool dressTracker(Tracker* aTracker, string fileName);

  // Extract the user-defined support structures from the geometry config file
  list<double>* parseSupportsFromFile(string fileName);

private:

  // Temporary streams and objects
  ifstream rawConfigFile_;
  stringstream configFile_;
  Tracker* myTracker_;

  // Generic parsing functions
  string getTill(istream &inStream, char delimiter, bool singleWord, bool allowNothing /* = false */);
  bool parseParameter(string& name, string& value, istream &inStream);
  
  // Parsing functions for the tracker building
  bool parseTracker(string myName, istream &inStream);
  bool parseBarrel(string myName, istream &inStream);
  bool parseEndcap(string myName, istream &inStream);
  bool parseObjectType(string myType);

  // Parsing functions for the tracker dressing
  bool parseDressType(string myType);
  bool breakParameterName(string& parameterName, int& ringIndex, int& diskIndex);
  bool parseBarrelType(string myName, istream &inStream);
  bool parseEndcapType(string myName, istream &inStream);
  bool parseAnyType(string myName, istream &inStream);
  bool parseOutput(istream &inStream);

  // Parsing function for inactive surfactes
  bool parseSupportParameters(istream& inStream, list<double>& list);

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

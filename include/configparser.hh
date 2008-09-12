#ifndef _CONFIGPARSER_HH_
#define _CONFIGPARSER_HH_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <exception>
#include "tracker.hh"


using namespace std;

class configParser {
public:
  configParser();
  ~configParser();
  bool parseFile(string fileName);
  Tracker* getTracker() { return myTracker_;};

private:
  ifstream rawConfigFile_;
  stringstream configFile_;
  bool parseType(string myType);

  string getTill(istream &inStream, char delimiter, bool singleWord, bool allowNothing /* = false */);
  
  bool parseTracker(string myName, istream &inStream);
  bool parseBarrel(string myName, istream &inStream);
  bool parseEndcap(string myName, istream &inStream);
  bool parseParameter(string& name, string& value, istream &inStream);

  // Tracker objects
  Tracker* myTracker_;

};


class parsingException: public exception
{
  virtual const char* what() const throw()
  {
    return "Error parsing the config file.";
  }
  
};

#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
//#include "tracker.hh"

using namespace std;

class configParser {
public:
  configParser();
  ~configParser();
  bool parseFile(string fileName);
private:
  ifstream rawConfigFile_;
  stringstream configFile_;
  bool parseType(string myType);

  string getTill(istream &inStream, char delimiter, bool singleWord, bool allowNothing /* = false */);
  
  bool parseBarrel(string myName, istream &inStream);
  bool parseEndcap(string myName, istream &inStream);
  bool parseParameter(string& name, string& value, istream &inStream);
};


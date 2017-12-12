/**
 * @file setup.cc
 * @brief called by the install.sh script, use mainConfigHandler for reading config file
 */
#include <global_constants.hh>
#include <MainConfigHandler.hh>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
  bool result = false;
  mainConfigHandler& m = mainConfigHandler::instance();
  // No arguments: just read the configuration file
  // or create it if necessary
  if (argc==1) {
    result = m.getConfiguration(false);
  } else if (argc==2) {
    string myArg = argv[1];
    if (myArg=="--checkDir") result = m.getConfiguration(true);
    if (myArg=="--dirNames") {
      m.getConfiguration(false);
      cout << "# Main config file" << endl;
      cout << "TKG_CONFIGFILE=\"" << m.getConfigFileName() << "\"" << endl;
      cout << "# Directories needing copy during install" << endl;
      cout << "TKG_SOURCE_MATTAB=\"" << insur::default_mattabdir << "\"" << endl;
      cout << "TKG_SOURCE_XML=\"" << insur::default_xmlpath << "\"" << endl;
      cout << "TKG_SOURCE_STYLE=\"" << insur::default_styledir << "\"" << endl;
      cout << "TKG_SOURCE_STDINCLUDE=\"" << insur::default_stdincludedir << "\"" << endl;
      cout << "TKG_SOURCE_GEOMETRIES=\"" << insur::default_geometriesdir << "\"" << endl;
      cout << "TKG_DESTINATION_MATTAB=\"" << m.getMattabDirectory() << "\"" << endl;
      cout << "TKG_DESTINATION_STDINCLUDE=\"" << m.getStandardIncludeDirectory() << "\"" << endl;
      cout << "TKG_DESTINATION_GEOMETRIES=\"" << m.getGeometriesDirectory() << "\"" << endl;
      cout << "TKG_DESTINATION_XML=\"" << m.getXmlDirectory() << "\"" << endl;
      cout << "TKG_DESTINATION_STYLE=\"" << m.getStyleDirectory() << "\"" << endl;
      result=true;
    }
  }

  if (result) return 0;
  else return -1;
}

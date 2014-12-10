#include <global_constants.h>
#include <mainConfigHandler.h>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
  bool result = false;
  mainConfigHandler m;
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
      cout << "TKG_DESTINATION_MATTAB=\"" << m.getMattabDirectory() << "\"" << endl;
      cout << "TKG_DESTINATION_XML=\"" << m.getXmlDirectory() << "\"" << endl;
      cout << "TKG_DESTINATION_STYLE=\"" << m.getStyleDirectory() << "\"" << endl;
      cout << "# Directoriesto be created during install" << endl;
      cout << "TKG_DESTINATION_ROOT=\"" << m.getRootfileDirectory() << "\"" << endl;
      cout << "TKG_DESTINATION_GRAPH=\"" << m.getGraphDirectory() << "\"" << endl;
      cout << "TKG_DESTINATION_SUMMARY=\"" << m.getSummaryDirectory() << "\"" << endl;
      result=true;
    }
  }

  if (result) return 0;
  else return -1;
}

#include "SimParms.h"


void SimParms::build() {
  try {
    check();

    std::ifstream filein(irradiationMapFile().c_str());
    if (!filein.is_open()) { 
      logERROR("Failed opening irradiation map file!");
    }
    std::string line;
    while(std::getline(filein, line)) {
      if (line.find_first_of("#//;")==0 || line=="") continue;
      std::stringstream ss(line);
      double z, r = -1.0, fluence = -1.0, error = -1.0; // set to -1.0 to check for parsing errors
      ss >> z;
      ss >> r;
      ss >> fluence;
      ss >> error;
      if (r < 0.0 || fluence < 0.0) {
        std::ostringstream tempSS;
        tempSS << "Error while parsing irradiation map line: " << z << " ," << r << " ," << fluence << " ," << error;
        logERROR(tempSS);
      }
      irradiationMap_[std::make_pair(int(z/2.5),int(r/2.5))] = fluence;
    }
    cleanup();
    builtok(true);
  } catch (PathfulException& pe) { pe.pushPath("SimParms"); throw; }
}

#include "SimParms.h"


void SimParms::crosscheck() {
  try {
    check();

    //iter between irradiation map file names and feed the irradiationMapsManager
    for(std::vector<std::string>::const_iterator iterMapFile = irradiationMapFiles.begin(); iterMapFile != irradiationMapFiles.end(); ++ iterMapFile) {
      irradiationMapsManager_.addIrradiationMap((*iterMapFile).c_str());
    }

    builtok(true);
    cleanup();
  } catch (PathfulException& pe) { pe.pushPath("SimParms"); throw; }
}

void SimParms::addIrradiationMapFile(std::string path) {
  irradiationMapFiles.appendString(path);
}


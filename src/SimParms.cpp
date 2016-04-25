#include "SimParms.h"

#include <mainConfigHandler.h>


//
// Static instance of this class
//
SimParms * SimParms::s_instance = nullptr;

//
// Method constructing static instance of this class
//
SimParms * SimParms::getInstance()
{
  if (s_instance == nullptr) s_instance = new SimParms;
  return s_instance;
}

//
// Check that sim parameters read-in correctly
//
void SimParms::crosscheck() {
  try {
    check();
    builtok(true);
    cleanup();
  } catch (PathfulException& pe) { pe.pushPath("SimParms"); throw; }
}

//
// Read-in irradiation maps
//
void SimParms::readIrradiationMaps()
{
  for (auto file : insur::default_irradiationfiles) {

    std::string path = mainConfigHandler::instance().getIrradiationDirectory() + "/" + file;
    irradiationMapFiles.appendString(path);
    m_irradiationMapsManager.addIrradiationMap(path.c_str());
  }
}

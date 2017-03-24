/*
 * MaterialTab.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: smartina
 */

// System include files
#include <fstream>
#include <sstream>
#include <stdexcept>

// tkLayout include files
#include "MaterialTab.h"
#include "global_constants.h"
#include "MainConfigHandler.h"
#include "MessageLogger.h"
#include "Units.h"


const std::string MaterialTab::msg_no_mat_file        = "Material tab file does not exist.";
const std::string MaterialTab::msg_no_mat_file_entry1 = "Material '";
const std::string MaterialTab::msg_no_mat_file_entry2 = "' not found in Material tab file.";

//
// Singleton private constructor -> reads-in default file, where all materials are defined, converts quantities into tkLayout natural units
//
MaterialTab::MaterialTab()
{
  std::string mattabFile(MainConfigHandler::getInstance().getMattabDirectory() + "/" + default_mattabfile);
  std::ifstream mattabStream(MainConfigHandler::getInstance().getMattabDirectory() + "/" + default_mattabfile);
  std::string line;
  std::string material;
  std::istringstream lineStream;

  double density, radLength, intLength;

  if (mattabStream.good()) {

    while (!mattabStream.eof()) {
      std::getline(mattabStream, line);
      lineStream.str(line);
      lineStream >> material;

      //check if is a comment
      if (material[0] != '#') {
        lineStream >> density >> radLength >> intLength;

        // Convert to natural tkLayout units
        density   *= Units::g/Units::cm3; // Quantities expected in g/cm3 on the input
        radLength *= Units::g/Units::cm2; // Quantities expected in g/cm2 on the input
        intLength *= Units::g/Units::cm2; // Quantities expected in g/cm2 on the input

        insert(make_pair(material, make_tuple(density, radLength, intLength)));
      }

      lineStream.clear();
    }
  } else {
    logERROR(msg_no_mat_file);
  }
}

//
// Material access method -> get instance of singleton class SimParms
//
const MaterialTab& MaterialTab::getInstance() {

  static MaterialTab s_instance;
  return s_instance;
}

//
// Get material info - density
//
double MaterialTab::density(std::string material) const {
  double val = 0;
  try {
    val = std::get<0>(at(material));
  } catch (const std::out_of_range& ex) {
    logERROR(msg_no_mat_file_entry1 + material + msg_no_mat_file_entry2);
    val = -1;
  }
  return val;
}

//
// Get material info - rad. length
//
double MaterialTab::radiationLength(std::string material) const {
  double val = 0;
  try {
    val = std::get<1>(at(material));
  } catch (const std::out_of_range& ex) {
    logERROR(msg_no_mat_file_entry1 + material + msg_no_mat_file_entry2);
    val = -1;
  }
  return val;
}

//
// Get material info - rad. length
//
double MaterialTab::interactionLength(std::string material) const {
  double val = 0;
  try {
    val = std::get<2>(at(material));
  } catch (const std::out_of_range& ex) {
    logERROR(msg_no_mat_file_entry1 + material + msg_no_mat_file_entry2);
    val = -1;
  }
  return val;
}

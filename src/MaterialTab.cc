/*
 * MaterialTab.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: smartina
 */


#include <fstream>
#include <sstream>
#include "MaterialTab.hh"
#include "global_constants.hh"
#include "MainConfigHandler.hh"
#include <MessageLogger.hh>
#include <stdexcept>

namespace material {

  const std::string MaterialTab::msg_no_mat_file = "Material tab file does not exist.";
  const std::string MaterialTab::msg_no_mat_file_entry1 = "Material '";
  const std::string MaterialTab::msg_no_mat_file_entry2 = "' not found in Material tab file.";

  MaterialTab::MaterialTab() {
    std::string mattabFile(mainConfigHandler::instance().getMattabDirectory() + "/" + insur::default_mattabfile);
    std::ifstream mattabStream(mainConfigHandler::instance().getMattabDirectory() + "/" + insur::default_mattabfile);
    std::string line;
    std::string material;
    std::istringstream lineStream;
    double density, radiationLength, interactionLength;

    if (mattabStream.good()) {
      while (!mattabStream.eof()) {
        std::getline(mattabStream, line);
        lineStream.str(line);
        lineStream >> material;

        //check if is a comment
        if (material[0] != '#') {
          lineStream >> density >> radiationLength >> interactionLength;
          density /= 1000; // convert g/cm3 in g/mm3
          insert(make_pair(material, make_tuple(density, radiationLength, interactionLength)));
        }

        lineStream.clear();
      }
    } else {
      logERROR(msg_no_mat_file);
    }
  }

  const MaterialTab& MaterialTab::instance() {
    static MaterialTab instance_;
    return instance_;
  }

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
} /* namespace material */

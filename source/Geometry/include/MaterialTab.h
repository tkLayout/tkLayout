/**
 * @file MaterialTab.h
 *
 * @date Jul 29, 2014
 * @author Stefano Martina
 */

#ifndef INCLUDE_MATERIALTAB_H_
#define INCLUDE_MATERIALTAB_H_

#include <map>
#include <tuple>

//! Define typedef: material name -> density, radLength, interactLength
typedef std::map<std::string, std::tuple<double, double, double> > MaterialTabType;

/*
 * @brief Singleton class providing interface for all materials as defined by user in default_mattabfile file.
 * @details User defines all materials in a default_mattabfile. For each material density, rad. length and int.
 * length are defined in predefined units. These quantities are translated into tkLayout natural units and
 * accessible through the overall program by calling instance of this class and using a unique name of the
 * material.
 */
class MaterialTab : public MaterialTabType {

 public:

  //! Material access method -> get instance of singleton class SimParms
  static const MaterialTab& getInstance();

  //! Get material info - density
  double density(std::string material) const;

  //! Get material info - rad. length
  double radiationLength(std::string material) const;

  //! Get material info - int. length
  double interactionLength(std::string material) const;

 private:

  //! Singleton private constructor -> reads-in default file, where all materials are defined, converts quantities into tkLayout natural units
  MaterialTab();

  static const std::string msg_no_mat_file;
  static const std::string msg_no_mat_file_entry1;
  static const std::string msg_no_mat_file_entry2;

}; // Class

#endif /* INCLUDE_MATERIALTAB_H_ */

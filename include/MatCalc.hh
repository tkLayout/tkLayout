//
// File:   MatCalc.h
// Author: ndemaio
//
// Created on October 24, 2008, 11:24 AM
//

/**
 * @file MatCalc.h
 * @brief This is the base class header file for the algorithms that assign materials to the elements of the tracker geometry
 */

#ifndef _MATCALC_H
#define	_MATCALC_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <list>
#include <string>
#include <stdexcept>
#include <ModuleCap.hh>
#include <MaterialTable.hh>
#include <InactiveElement.hh>
namespace insur {
  /**
   * Messages that may appear during the calculations
   */
  static const std::string err_no_such_type = "Error: the requested type is not on the list.";
  static const std::string err_no_material = "Error: no material with the specified tag was found.";
  static const std::string err_no_service = "Error: no service material with the specified properties was found.";
  static const std::string err_no_support = "Error: no support material with the specified properties was found.";
  static const std::string err_unknown_type = "Error: unknown module type.";
  static const std::string err_up_general = "Error updating parameter entry.";
  static const std::string err_conversion = "Error: material unit not recognised during conversion.";
  static const std::string err_matadd_weird = "Something weird happened when trying to add an entry to one of the vectors for material parameters...";
  static const std::string msg_negative_area = "Warning: module surface is negative.";
  static const std::string msg_abort = "Aborting function.";
  static const std::string msg_ignore_tag = " Ignoring module material labelled ";
  /**
   * @class MatCalc
   * @brief The MatCalc class provides the core material assignment algorithm for a given tracker geometry.
   *
   * Once its internal data structures have been initialised from the material config file by the <i>MatParser</i> class,
   * it uses that information, combined with the geometry and position of an individual tracker element, to set
   * the local materials vector that element before getting it to calculate its overall mass, its radiation length
   * and its interaction length. The tracker geometry is ready for further study afterwards. Some of the access functions
   * for various internal list elements may return an exception if the requested element does not exist on the list.
   */
  class MatCalc {
  public:
    /**
     * @enum Matunit The allowed measurement units for the quantities of material that the config file defines
     */
    enum Matunit { gr, mm3, mm, grpm };
    MatCalc() { init_done = false; }
    virtual ~MatCalc() {}
    bool initDone();
    void initDone(bool yes);
    void reset();
    MaterialTable& getMaterialTable();
    bool init_done;
    MaterialTable mt;
  };
}
#endif	/* _MATCALC_H */


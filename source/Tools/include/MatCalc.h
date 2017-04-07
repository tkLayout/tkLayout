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
#include <ModuleCap.h>
#include <MaterialTable.h>
#include <InactiveElement.h>
#include <DetectorModule.h>
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
   * the local and exiting materials vector that element before getting it to calculate its overall mass, its radiation length
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
    int getStripsAcross(std::string type); // throws exception
    int getSegmentsAlong(std::string type); // throws exception
    std::pair<int, int> getDefaultDimensions(std::string type); // throws exception
    void addTypeInfo(std::string type, int strips, int segments);
    void updateTypeInfoStrips(std::string type, int strips);
    void updateTypeInfoSegments(std::string type, int segments);
    void addModuleParameters(std::string tag, std::string type, std::string comp,
                             double A, Matunit uA, double B, Matunit uB, double C, Matunit uC, double D, Matunit uD, bool local);
    void addServiceParameters(std::string tag, double Q, Matunit uQ);
    void addServiceParameters(std::string tagIn, double In, Matunit uIn, std::string tagOut, double Out, Matunit uOut, bool local);
    void addSupportParameters(std::string tag, double M, Matunit uM, MaterialProperties::Category cM);
    void clearModVectors();
    void clearModVector(std::string type);
    void copyContents(std::string source, std::string dest);
    void appendContents(std::string source, std::string dest);
    bool typeRegistered(std::string type);
    unsigned int registeredTypes();
    MaterialTable& getMaterialTable();
    virtual bool calculateBarrelMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps);
    virtual bool calculateEndcapMaterials(std::vector<std::vector<ModuleCap> >& endcapcaps);
    virtual bool calculateBarrelServiceMaterials(std::vector<std::vector<ModuleCap> >& barrelcaps,
                                                 std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices);
    virtual bool calculateEndcapServiceMaterials(
      std::vector<std::vector<ModuleCap> >& endcapcaps,
      std::vector<InactiveElement>& barrelservices, std::vector<InactiveElement>& endcapservices);
    virtual bool calculateSupportMaterials(std::vector<InactiveElement>& supports);
    void printInternals();
  protected:
    /**
     * @struct TypeInfo
     * @brief Instances of this struct contain the default dimensions for modules of a certain type.
     * @param type One of <i>rphi</i>, <i>stereo</i> or <i>pt</i> indicating the module type
     * @param strips_across The number of chips across the module
     * @param segments_along The number of segments along the module
     */
    struct TypeInfo {
      std::string type;
      int strips_across;
      int segments_along;
    };
    /**
     * @struct SingleMod
     * @brief This struct contains the values for a single material that contributes to the mix within the modules.
     * @param tag A string that identifies the material unambiguously
     * @param comp A string that identifies the component originating the material
     * @param A The component that depends on the position of the module on the rod as well as its layout of chips and segments
     * @param B The component that depends on the layout of chips and segments on an individual module
     * @param C The component that depends on the position of a module on a rod
     * @param D The constant component local to every module in the tracker
     * @param uA The unit of component A
     * @param uB The unit of component B
     * @param uC The unit of component C
     * @param uD The unit of component D
     * @param is_local A flag that states whether the described material is local to the module or exits the layer
     */
    struct SingleMod {
      std::string tag;
      std::string comp;
      double A, B, C, D;
      Matunit uA, uB, uC, uD;
      bool is_local;
    };
    /**
     * @struct SingleSerLocal
     * @brief This struct contains the values for the local components of a service material at the layer-service boundary.
     * @param tag A string that identifies the material unambiguously
     * @param Q The amount of material that needs to be scaled to the actual size of the service volume
     * @param uQ The unit of component Q
     */
    struct SingleSerLocal {
      std::string tag;
      double Q;
      Matunit uQ;
    };
    /**
     * @struct SingleSerExit
     * @brief This struct contains one mapping from layer material to service material as it is used to describe the layer-service boundary.
     * @param tagIn A string that identifies the layer material unambiguously
     * @param tagOut A string that identifies the resulting service material unambiguously
     * @param In The amount of layer material that needs to be mapped to the service
     * @param Out The amount of service material that results from the given amount of layer material
     * @param uIn The unit of component In
     * @param uOut The unit of component Out
     * @param is_local A flag that states whether the resulting material is local to the service or continues outward to the next volume
     */
    struct SingleSerExit {
      std::string tagIn, tagOut;
      double In, Out;
      Matunit uIn, uOut;
      bool is_local;
    };
    /**
     * @struct SingleSup
     * @brief This struct contains one material that is used for a certain category of support structure
     * @param tag A string that identifies the supporting material unambiguously
     * @param M The amount of material that needs to be scaled to the actual size of the support volume
     * @param uM The unit of component M
     * @param cM The category of support structure component M is used for
     */
    struct SingleSup {
      std::string tag;
      double M;
      Matunit uM;
      MaterialProperties::Category cM;
    };
    /**
     * @struct MatInfo
     * @brief An instance of this struct bundles the material information for all types of volumes as it comes out of the config file.
     * @param modinfo A nested arrangement of vectors and pairs listing the module types and their associated parameters and materials
     * @param serlocalinfo A vector describing the local materials in services at the layer-service boundary
     * @param serexitinfo A vector describing the material mappings at the layer-service boundary
     * @param supinfo A vector describing the materials that are found in the various categories of support structures
     */
    struct MatInfo {
      std::vector<std::pair<TypeInfo, std::vector<SingleMod> > > modinfo;
      std::vector<SingleSerLocal> serlocalinfo;
      std::vector<SingleSerExit> serexitinfo;
      std::vector<SingleSup> supinfo;
    };
    bool init_done;
    MaterialTable mt;
    MatInfo internals;
    std::vector<SingleMod>& getModVector(std::string type); // throws exception
    TypeInfo& getTypeInfoByType(std::string type); // throws exception
    SingleMod& getSingleMod(std::string tag, std::string type, std::string comp, Matunit uA, Matunit uB, Matunit uC, Matunit uD, bool local); // throws exception
    SingleSerLocal& getSingleSer(std::string tag, Matunit u); // throws exception
    SingleSerExit& getSingleSer(std::string tag1, std::string tag2, Matunit u1, Matunit u2, bool local); // throws exception
    SingleSup& getSingleSup(std::string tag, Matunit uM, MaterialProperties::Category cM); // throws exception
  private:
    bool entryExists(std::string type);
    bool entryExists(std::string tag, std::string type, std::string comp, Matunit uA, Matunit uB, Matunit uC, Matunit uD, bool local);
    bool entryExists(std:: string tag, Matunit uQ);
    bool entryExists(std::string tag1, std::string tag2, Matunit uIn, Matunit uOut, bool local);
    bool entryExists(std::string tag, Matunit uM, MaterialProperties::Category cM);
    int findBarrelRods(std::vector<std::vector<ModuleCap> >& caps, int layer);
    int findEndcapRods(std::vector<std::vector<ModuleCap> >& caps, int layer);
    double convert(double value, Matunit unit, double densityorlength, double surface = 0); // throws exception
    void adjacentDifferentCategory(std::vector<ModuleCap>& source, InactiveElement& dest, int r, double l, double s);
    void adjacentSameCategory(InactiveElement& source, InactiveElement& dest);
  };
#endif	/* _MATCALC_H */


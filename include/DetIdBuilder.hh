#ifndef _DETIDBUILDER_HH
#define	_DETIDBUILDER_HH

#include "global_constants.hh"
#include "global_funcs.hh"
#include <Tracker.hh>

/** See https://github.com/tkLayout/tkLayout/wiki/DetIds-in-tkLayout-for-the-entire-Tracker .
 * All geometry hierarchy levels (for example OT Barrel, or Ring) are assigned a size (number of bits) and an Id.
 * Each module or Sensor DetId is a 32-bit integer, which gathers all this information.

 * The process to calculate the DetIds is as follows :
 * 1) The detIdScheme, specified in the cfg files, provides the sizes (in number of bits) associated to each geometry hierarchy level (geometryHierarchySizes).
 * 2) This allows to calculate the geometryHierarchyIds : the Ids associated to each geometry hierarchy level.
 * This is done here, with BarrelDetIdBuilder and EndcapDetIdBuilder.
 * 3) Lastly, one can easily calculate the DetId for each Module and Sensor. 
 * This is done using DetIdentifiable::buildDetId(geometryHierarchySizes, geometryHierarchyIds) method.
 */


    //************************************//
    //*               Visitor             //
    //*          BarrelDetIdBuilder       //
    //*                                   //
    //************************************//
class BarrelDetIdBuilder : public SensorGeometryVisitor {

public:
  BarrelDetIdBuilder(bool isPixelTracker, std::vector<int> geometryHierarchySizes);

  void visit(Barrel& b);
  void visit(Layer& l);
  void visit(RodPair& r);
  void visit(BarrelModule& m);
  void visit(Sensor& s);

private:
  bool isPixelTracker_;
  std::vector<int> geometryHierarchySizes_;  // Size (in number of bits) associated to each level in the geometry hierarchy.
  std::map< std::pair<std::string, int>, int > sortedLayersIds_;  // Full subdetector layers sorted by radius.
  // For example, TBPS Layer 1 is asigned number 1, TBPS Layer 2 is asigned number 2, TB2S Layer 1 is assigned number 4.

  std::map<int, uint32_t> geometryHierarchyIds_;  // WHAT IS CALCULATED HERE !! Id associated to each level in the geometry hierarchy.

  bool isTiltedLayer_;
  int numRods_;
  int numFlatRings_;
  int numRings_;
  bool isCentered_;
  uint32_t phiRef_;
};



    //************************************//
    //*               Visitor             //
    //*          EndcapDetIdBuilder       //
    //*                                   //
    //************************************//
class EndcapDetIdBuilder : public SensorGeometryVisitor {

public:
  EndcapDetIdBuilder(bool isPixelTracker, std::vector<int> geometryHierarchySizes);

  void visit(Endcap& e);
  void visit(Disk& d);
  void visit(Ring& r);
  void visit(EndcapModule& m);
  void visit(Sensor& s);

private:
  bool isPixelTracker_;
  std::vector<int> geometryHierarchySizes_;  // Size (in number of bits) associated to each level in the geometry hierarchy.
  std::map< std::tuple<std::string, int, bool>, int > sortedDisksIds_;  // Full subdetector disks sorted by Z.
  // For example, TEDD_1 Disk 1 is asigned number 1, TEDD_1 Disk 2 is asigned number 2, TEDD_2 Disk 1 is assigned number 3.

  std::map<int, uint32_t> geometryHierarchyIds_;  // WHAT IS CALCULATED HERE !! Id associated to each level in the geometry hierarchy.

  int numEmptyRings_;
  int numModules_;
};


#endif // _DETIDBUILDER_HH

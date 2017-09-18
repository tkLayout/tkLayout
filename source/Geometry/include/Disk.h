#ifndef INCLUDE_DISK_H_
#define INCLUDE_DISK_H_

#include <memory>
#include <vector>
#include <string>

#include "Property.h"
#include "Ring.h"
#include "Layer.h"
#include "Visitable.h"
#include "MaterialObject.h"

namespace material {
  class ConversionStation;
}

using material::MaterialObject;
using material::ConversionStation;

// Typedefs
typedef PtrVector<Disk> Disks;

/*
 * @class Disk
 * @details Disk class holds information about individual endcap disks. It's building procedure is executed automatically via
 * Endcap class. Similarly, all its components (rings -> modules) are recursively build through build() method. Depending on
 * the configuration, building algorithm uses either the bottom-up or up-down approach. i.e. from the closest ring to the
 * beam-pipe to the furthest one or vice-versa. The positioning algorithm always tests two extreme cases to fully cover
 * eta region and simultaneously position individual rings within the disk. Hence, the leftmost and rightmost edges of
 * the given endcap are investigated.
 */
class Disk : public PropertyObject, public Buildable, public Identifiable<int>, public Visitable {

 public:

  //! Constructor - specify unique id, const pointer to parent module & parse geometry config file using boost property tree & read-in module parameters
  Disk(int id, double zOffset, double zECCentre,  double zECHalfLength, const PropertyNode<int>& nodeProperty, const PropertyTree& treeProperty);

  //! Build recursively individual subdetector systems: rings -> modules & conversion stations
  void build(const vector<double>& buildDsDistances);
   
  //! Position newly individual rings of given disc -> clone (+Z) side disks and rotate them around FCC_Y with angle Pi.
  // One does not want to build different disks for both sides, hence no 'mirror' should be considered.
  void buildMirror(int id) { myid(id); rotateToNegativeZSide(); }

  //! Position newly individual rings of given disc -> offset them in Z after cloning
  void buildClone(int id, double offset) { myid(id); translateZ(offset); }

  //! Setup: link lambda functions to various layer related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Cross-check parameters provided from geometry configuration file
  void check() override;

  //! Limit disk geometry by eta cut
  void cutAtEta(double eta);

  //! Return disk rings
  const Rings& rings()           const { return m_rings; }
  const RingIndexMap& ringsMap() const { return m_ringIndexMap; }

  //! Return disk as a material object
  const MaterialObject& materialObject()       const { return m_materialObject; }

  //! Return first order conversion station -> TODO: check if needed to be updatable, use PtrVector instead
  ConversionStation* flangeConversionStation() const {return m_flangeConversionStation; }

  //! Return vector of second order conversion stations -> TODO: check if needed to be updatable, use reference instead
  ConversionStations secondConversionStations() const {return m_secondConversionStations; }

  //! GeometryVisitor pattern -> layer visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> layer visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  //! Return average disc Z position
  double averageZ()  const { return m_averageZ; }

  //! Return disc thickness
  double thickness() const { return bigDelta()*2 + maxRingThickness(); }

  Property<int   , NoDefault> numRings; //!< Required number of rings in the disk -> TODO: Compression as for barrel rods
  Property<double, Default>   zError;   //!< When positioning modules take into account beam spot spread in Z (read-in through SimParms interface)

  ReadonlyProperty<double, Computable> minZ;       //!< Disk minimum Z position
  ReadonlyProperty<double, Computable> maxZ;       //!< Disk maximum Z position
  ReadonlyProperty<double, Computable> minZAllMat; //!< Disk minimum Z position taking into account all material structures
  ReadonlyProperty<double, Computable> maxZAllMat; //!< Disk maximum Z position taking into account all material structures
  ReadonlyProperty<double, Computable> minR;       //!< Disk minimum radius
  ReadonlyProperty<double, Computable> maxR;       //!< Disk maximum radius
  ReadonlyProperty<int   , Computable> totalModules;     //!< Total number of modules
  ReadonlyProperty<double, Computable> maxRingThickness; //!< Maximum ring thickness

  Property<double, NoDefault> bigDelta;    //!< Ring versus another ring are positioned by +-bigDelta from the central Z pos.
  Property<double, Default>   rOverlap;    //!< Required ring overlap in radius

 private:

  //! Use top to bottom approach when building rings -> internally called by build method
  void buildTopDown(const vector<double>& buildDsDistances);

  //! Use bottom to top approach when building rings -> internally called by build method
   void buildBottomUp(const vector<double>& buildDsDistances);

  //! Helper method translating Disc z position by given offset
  void translateZ(double z);

  //! Helper method mirroring the whole Disc from zPos to -zPos or vice versa
  void rotateToNegativeZSide();

  Rings          m_rings;             //!< Disk rings
  RingIndexMap   m_ringIndexMap;
  MaterialObject m_materialObject;

  ConversionStation*              m_flangeConversionStation;  //!< First order layer conversion unit
  std::vector<ConversionStation*> m_secondConversionStations; //!< Vector of second order layer conversion units

  Property<double, NoDefault> m_innerRadius; //!< Disc innermost radius
  Property<double, NoDefault> m_outerRadius; //!< Disc outermost radius
  Property<int   , Default>   m_bigParity;   //!< Use +bigDelta or -bigDelta as starting value in the positioning algorithm

  double m_zEndcapHalfLength; //!< Z halflength of endcap in which the disc is to bebuilt
  double m_zEndcapCentre;     //!< Z central position of endcap in which the disc is to bebuilt
  double m_zOffset;           //!< Relative position of the disc with respect to endcap central position
  double m_averageZ = 0;      //!< Average Z position of the given disc

  PropertyNode<int>               m_ringNode;     //!< Property tree node for ring (to grab properties for specific ring modules)
  PropertyNodeUnique<std::string> m_stationsNode; //!< Property tree nodes for conversion stations (included geometry config file)

  inline double getDsDistance(const vector<double>& buildDsDistances, int rindex) const;
}; // Class

#endif /* INCLUDE_DISK_H_ */

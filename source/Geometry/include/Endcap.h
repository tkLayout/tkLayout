#ifndef INCLUDE_ENDCAP_H_
#define INCLUDE_ENDCAP_H_

#include <memory>
#include <vector>
#include <string>

#include "Disk.h"
#include "Property.h"
#include "Visitable.h"

// Forward declaration
class SupportStructure;

// Typedefs
typedef PtrVector<Endcap>           Endcaps;
typedef PtrVector<SupportStructure> EndcapSupportStructures;

/*
 * @class Endcap
 * @details Endcap class holds information about tracker endcap system. It's building procedure is executed automatically via
 * Tracker class. Similarly, all its components (disks -> rings -> modules) are recursively called through private build()
 * method.
 */
class Endcap : public PropertyObject, public Buildable, public Identifiable<std::string>, public Visitable {

 public:

  //! Constructor - parse geometry config file using boost property tree & read-in Disk, Support nodes
  Endcap(double barrelOuterZ, const std::string& name, const PropertyTree& nodeProperty, const PropertyTree& treeProperty);

  //! Build recursively individual subdetector systems: Disks -> rings -> modules
  void build();

  //! Setup: link lambda functions to various endcap related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Limit endcap geometry by eta cut
  void cutAtEta(double eta);

  //! Return endcap disks
  const Disks& disks() const { return m_disks; }

  //! Return endcap supports
  const EndcapSupportStructures& supports() const { return m_supportStructures; }

  //! GeometryVisitor pattern -> endcap visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> endcap visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  Property<        int   , NoDefault>  numDisks;     //!< Total number of discs to be built between innerZ and outerZ positions
  Property<        double, NoDefault>  innerZ;       //!< Start building disks at innerZ
  Property<        double, NoDefault>  outerZ;       //!< End building disks at outerZ
  ReadonlyProperty<double, Computable> minZ;         //!< Minimum Z position of an endcap
  ReadonlyProperty<double, Computable> maxZ;         //!< Maximum Z position of an endcap
  ReadonlyProperty<double, Computable> minZAllMat;   //!< Minimum Z position of an endcap taking into account all material structures
  ReadonlyProperty<double, Computable> maxZAllMat;   //!< Maximum Z position of an endcap taking into account all material structures
  ReadonlyProperty<double, Computable> minR;         //!< Minimum radius of an endcap
  ReadonlyProperty<double, Computable> maxR;         //!< Maximum radius of an endcap
  ReadonlyProperty<bool  , Default>    skipServices; // TODO: Comment

 private:

  //! Drill down into the property tree to find the maximum dsDistances -> used by build algorithm to optimize disk positions
  vector<double> findMaxDsDistances();

  Disks                   m_disks;              //!< Disks of given endcap
  EndcapSupportStructures m_supportStructures;  //!< Endcap supports

  Property<double, NoDefault>     m_barrelOuterZ; //!< Start building disks in Z given by maximum Z distance of all barrels + m_barrelGap
  Property<double, NoDefault>     m_barrelGap;    //!< Gap between the end-point in Z of barrel detector and starting point in Z of end-cap detector

  PropertyNode<int>               m_diskNode;     //!< Property tree nodes for disks (included geometry config file)
  PropertyNodeUnique<std::string> m_supportNode;  //!< Property tree nodes for endcap supports (included geometry config file)

}; // Class

#endif /* INCLUDE_ENDCAP_H_ */

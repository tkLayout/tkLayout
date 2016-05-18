#ifndef ENDCAP_H
#define ENDCAP_H

#include <vector>
#include <string>
#include <memory>

#include <boost/ptr_container/ptr_vector.hpp>

#include "Property.h"
#include "Visitable.h"

#include "Disk.h"

namespace material {
  class SupportStructure;
}

// Typedefs
typedef PtrVector<Disk>                       Disks;
typedef PtrVector<material::SupportStructure> EndcapSupportStructures;

/*
 * @class Endcap
 * @details Endcap class holds information about tracker endcap system. It's building procedure is executed automatically via
 * Tracker class. Similarly, all its components (disks -> rings -> modules) are recursively called through private build()
 * method.
 */
class Endcap : public PropertyObject, public Buildable, public Identifiable<std::string>, public Visitable {

 public:

  //! Constructor - parse geometry config file using boost property tree & read-in Disk, Support nodes
  Endcap(double brlMaxZ, const std::string& name, const PropertyTree& nodeProperty, const PropertyTree& treeProperty);

  //! Limit endcap geometry by eta cut
  void cutAtEta(double eta);

  //! Return endcap disks
  const Disks& disks() const { return m_disks; }

  //! Return endcap supports which can be updated
  EndcapSupportStructures& supportStructures() { return m_supportStructures; }

  //! GeometryVisitor pattern -> endcap visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> endcap visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  Property<        int   , NoDefault>  numDisks;     //!< Total number of discs to be built between minZ (innerZ) and maxZ (outerZ) positions
  Property<        double, NoDefault>  barrelMaxZ;   //!< Start building disks in Z given by maximum Z distance of all barrels + m_barrelGap
  Property<        double, NoDefault>  innerZ;       //!< Start building disks at innerZ (i.e. minimum Z)
  Property<        double, NoDefault>  outerZ;       //!< End building disks at outerZ (i.e. maximum Z)
  ReadonlyProperty<double, Computable> minZ;         //!< Minimum Z position of an endcap
  ReadonlyProperty<double, Computable> maxZ;         //!< Maximum Z position of an endcap
  ReadonlyProperty<double, Computable> minR;         //!< Minimum radius of an endcap
  ReadonlyProperty<double, Computable> maxR;         //!< Maximum radius of an endcap
  ReadonlyProperty<bool  , Default>    skipServices; // TODO: Comment

 private:

  //! Build recursively individual subdetector systems: Disks -> rings -> modules -> private method called by constructor
  void build();

  //! Calculate various endcap related properties -> private method called by constructor
  void setup();

  //! Drill down into the property tree to find the maximum dsDistances -> used by build algorithm to optimize disk positions
  vector<double> findMaxDsDistances();

  Disks                   m_disks;              //!< Disks of given endcap
  EndcapSupportStructures m_supportStructures;  //!< Endcap supports

  Property<double, NoDefault>     m_barrelGap;  //!< Gap between the end-point in Z of barrel detector and starting point in Z of end-cap detector

  PropertyNode<int>               m_diskNode;   //!< Property tree nodes for disks (included geometry config file)
  PropertyNodeUnique<std::string> m_supportNode;//!< Property tree nodes for endcap supports (included geometry config file)

};



#endif

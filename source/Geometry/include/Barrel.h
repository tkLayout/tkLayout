#ifndef INCLUDE_BARREL_H_
#define INCLUDE_BARREL_H_

#include <vector>
#include <string>
#include <memory>

#include <boost/ptr_container/ptr_vector.hpp>

#include "Property.h"
#include "Visitable.h"

#include "Layer.h"

// Forward declaration
namespace material {
  class SupportStructure;
}

// Typedefs
typedef PtrVector<Layer>                      Layers;
typedef PtrVector<material::SupportStructure> BarrelSupportStructures;

/*
 * @class Barrel
 * @details Barrel class holds information about tracker barrel system. It's building procedure is expected to be executed
 * automatically via Tracker class build() method. Similarly, all barrel components (layers -> rods -> modules) are
 * recursively built through barrel build() method. Call setup() method after the geometry is built to assign (lambda)
 * functions to various barrel related properties.
 */
class Barrel : public PropertyObject, public Buildable, public Identifiable<string>, Clonable<Barrel>, public Visitable {

 public:
  
  //! Constructor - parse geometry config file using boost property tree & read-in Layer, Support nodes
  Barrel(const std::string& name, const PropertyTree& nodeProperty, const PropertyTree& treeProperty);

  //! Build recursively individual subdetector systems: Layers -> rods -> modules
  void build();

  //! Setup: link lambda functions to various barrel related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Limit barrel geometry by eta cut
  void cutAtEta(double eta);

  //! Return barrel layers
  const Layers& layers() const { return m_layers; }

  //! Return barrel supports which can be updated
  BarrelSupportStructures& supportStructures() { return m_supportStructures; }

  //! GeometryVisitor pattern -> barrel visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> barrel visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  Property<        int   , NoDefault>  numLayers;   //!< Number of layers in a barrel
  ReadonlyProperty<double, Computable> minZ;        //!< Minimum Z position of a barrel
  ReadonlyProperty<double, Computable> maxZ;        //!< Maximum Z position of a barrel
  ReadonlyProperty<double, Computable> minR;        //!< Minimum radius of a barrel
  ReadonlyProperty<double, Computable> maxR;        //!< Maximum radius of a barrel
  ReadonlyProperty<bool  , Default>    skipServices;// TODO: Comment

 private:

  Layers                  m_layers;                 //!< Layers of given barrel
  BarrelSupportStructures m_supportStructures;      //!< Barrel supports

  Property<double, NoDefault> m_innerRadius;        //!< Starting barrel inner radius (algorithm may optimize its value)
  Property<double, NoDefault> m_outerRadius;        //!< Starting barrel outer radius (algorithm may optimize its value)
  Property<double, Default>   m_barrelRotation;     //!< Start rod (ladder) positioning in R-Phi at Phi=barrelRotation [rad]
  Property<double, Default>   m_supportMarginOuter; // TODO: Comment
  Property<double, Default>   m_supportMarginInner; // TODO: Comment
  Property<bool  , Default>   m_innerRadiusFixed;   //!< Is inner radius fixed or floating -> internal algorithm finds optimal radius
  Property<bool  , Default>   m_outerRadiusFixed;   //!< Is outer radius fixed or floating -> internal algorithm finds optimal radius

  PropertyNode<int>               m_layerNode;      //!< Property tree nodes for layers (included geometry config file)
  PropertyNodeUnique<std::string> m_supportNode;    //!< Property tree nodes for barrel supports (included geometry config file)

}; // Class

#endif /* INCLUDE_BARREL_H_ */

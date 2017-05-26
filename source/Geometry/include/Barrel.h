#ifndef INCLUDE_BARREL_H_
#define INCLUDE_BARREL_H_

#include <memory>
#include <vector>
#include <string>

#include "Layer.h"
#include "Property.h"
#include "Visitable.h"

// Forward declaration
class InactiveElement;
class SupportStructure;

// Typedefs
typedef PtrVector<Barrel>           Barrels;
typedef PtrVector<SupportStructure> BarrelSupportStructures;
typedef PtrVector<InactiveElement>  BarrelServices;

/*
 * @class Barrel
 * @details Barrel class holds information about tracker barrel system. It's building procedure is expected to be executed
 * automatically via Tracker class build() method. Similarly, all barrel components (layers -> rods -> modules) are
 * recursively built through barrel build() method.
 */
class Barrel : public PropertyObject, public Buildable, public Identifiable<string>, Clonable<Barrel>, public Visitable {

 public:
  
  //! Constructor - parse geometry config file using boost property tree & read-in Layer, Support nodes
  Barrel(const std::string& name, const PropertyTree& nodeProperty, const PropertyTree& treeProperty);

  //! Build recursively individual subdetector systems: Layers -> rods -> modules
  void build();

  //! Setup: link lambda functions to various barrel related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Cross-check parameters provided from geometry configuration file
  void check() override;

  //! Limit barrel geometry by eta cut
  void cutAtEta(double eta);

  //! Return barrel layers
  const Layers& layers() const { return m_layers; }

  //! Return barrel supports
  const BarrelSupportStructures& supports() const { return m_supportStructures; }

  //! Add barrel service line
  void addServiceLine(InactiveElement* service);

  //! Return barrel services
  const BarrelServices& services() const { return m_services; }

  //! Update barrel services
  BarrelServices& updateServices() { return m_services; }

  //! GeometryVisitor pattern -> barrel visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> barrel visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  Property<        int   , NoDefault>  numLayers;   //!< Number of layers in a barrel
  ReadonlyProperty<double, Computable> minZ;        //!< Minimum Z position of a barrel
  ReadonlyProperty<double, Computable> maxZ;        //!< Maximum Z position of a barrel
  ReadonlyProperty<double, Computable> minR;        //!< Minimum radius of a barrel
  ReadonlyProperty<double, Computable> maxR;        //!< Maximum radius of a barrel
  ReadonlyProperty<double, Computable> minRAllMat;  //!< Minimum radius of a barrel taking into account all material structures
  ReadonlyProperty<double, Computable> maxRAllMat;  //!< Maximum radius of a barrel taking into account all material structures
  ReadonlyProperty<bool  , Default>    skipServices;// TODO: Comment

 private:

  Layers                  m_layers;                 //!< Layers of given barrel
  BarrelSupportStructures m_supportStructures;      //!< Barrel supports
  BarrelServices          m_services;               //!< Barrel services

  Property<double, NoDefault> m_innerRadius;        //!< Starting barrel inner radius (algorithm may optimize its value)
  Property<double, NoDefault> m_outerRadius;        //!< Starting barrel outer radius (algorithm may optimize its value)
  Property<double, Default>   m_barrelRotation;     //!< Start rod (ladder) positioning in R-Phi at Phi=barrelRotation [rad]
  Property<bool  , Default>   m_innerRadiusFixed;   //!< Is inner radius fixed or floating -> internal algorithm finds optimal radius
  Property<bool  , Default>   m_outerRadiusFixed;   //!< Is outer radius fixed or floating -> internal algorithm finds optimal radius
  Property<bool  , Default>   m_sameRods;           //!< Build the same rods in the whole barrel -> advantage from the engineering point of view


  PropertyNode<int>               m_layerNode;      //!< Property tree nodes for layers (included geometry config file)
  PropertyNodeUnique<std::string> m_supportNode;    //!< Property tree nodes for barrel supports (included geometry config file)

}; // Class

#endif /* INCLUDE_BARREL_H_ */

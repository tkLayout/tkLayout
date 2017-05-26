/*
 * BeamPipe.h
 *
 *  Created on: 11. 5. 2016
 *      Author: Z.Drasal (CERN)
 */

#ifndef BEAMPIPE_H_
#define BEAMPIPE_H_

#include <memory>

#include "Property.h"
#include "Visitable.h"

// Forward declaration
class ConstGeometryVisitor;
class GeometryVisitor;
class InactiveTube;

/*
 * @class BeamPipe
 * @brief Geometry & material object holding information about the beam pipe
 */
class BeamPipe : public PropertyObject, public Identifiable<string>, public Buildable, Visitable {

 public:

  //! Constructor -> object set through build method
  BeamPipe(const PropertyTree& treeProperty);

  //! Destructor
  ~BeamPipe();

  //! Build method - setting all parameters
  void build();

  //! Setup: link lambda functions to various beampipe related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup() {};

  //! GeometryVisitor pattern -> beam pipe visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> beam pipe visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  //! Get beam pipe material information
  //! @return InactiveTube, i.e. returns beam pipe as inactive element
  const InactiveTube* getMaterial() const;

  // Beam pipe radius, thickness, thickness in rad. length, in int. length
  ReadonlyProperty<double, Default>   radius;           //!< Beam pipe radius
  ReadonlyProperty<double, Default>   thickness;        //!< Beam pipe thickness
  ReadonlyProperty<double, Default>   maxZ;             //!< Beam pipe extends from minimum Z to maximum Z
  ReadonlyProperty<double, Default>   radLength;        //!< Beam pipe rad. length @ 90 deg
  ReadonlyProperty<double, Default>   intLength;        //!< Beam pipe int. lenght @ 90 deg

 private:

  std::unique_ptr<InactiveTube> m_tube; //!< Beam pipe as inactive tube;

}; // Class

#endif /* BEAMPIPE_H_ */

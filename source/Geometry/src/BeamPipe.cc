/*
 * BeamPipe.cc
 *
 *  Created on: 11. 5. 2016
 *      Author: drasal
 */

#include "BeamPipe.h"

#include "InactiveTube.h"
#include "Visitor.h"

//
// Constructor -> object set through build method
//
BeamPipe::BeamPipe(const PropertyTree& treeProperty) :
 radius(   "radius"   , parsedAndChecked(), 0.0),
 thickness("thickness", parsedAndChecked(), 0.0),
 maxZ(     "maxZ"     , parsedAndChecked(), 0.0),
 radLength("radLength", parsedAndChecked(), 0.0),
 intLength("intLength", parsedAndChecked(), 0.0),
 m_tube(nullptr)
{
  // Set the geometry config parameters
  this->myid(treeProperty.data());
  this->store(treeProperty);
}

BeamPipe::~BeamPipe()
{
  m_tube.reset();
}

//
// Buidl method - setting all parameters
//
void BeamPipe::build()
{
  try {
    check();

    // Build beam pipe as inactive material
    double zLength = 2*maxZ();
    double zOffset = 0.0;

    m_tube.reset(new InactiveTube(-zLength/2.,zLength,radius(),thickness()));
    m_tube->setRadiationLength(radLength());
    m_tube->setInteractionLength(intLength());
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }

  cleanup();
  builtok(true);
}

//
// GeometryVisitor pattern -> beam pipe visitable
//
void BeamPipe::accept(GeometryVisitor& v)
{
  v.visit(*this);
}

//
// GeometryVisitor pattern -> beam pipe visitable (const. option)
//
void BeamPipe::accept(ConstGeometryVisitor& v) const
{
  v.visit(*this);
}

//
// Get beam pipe as inactive tube
//
const InactiveTube* BeamPipe::getMaterial() const
{
  const InactiveTube* tube = m_tube.get();
  return tube;
}

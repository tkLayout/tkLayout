#include "Barrel.h"

#include "InactiveElement.h"
#include "Layer.h"
#include "MessageLogger.h"
#include "SupportStructure.h"

//
// Constructor - parse geometry config file using boost property tree & read-in Layer, Support nodes
//
Barrel::Barrel(const std::string& name, const PropertyTree& nodeProperty, const PropertyTree& treeProperty) :
 numLayers(           "numLayers"         , parsedAndChecked()),
 minZ(         string("minZ")             ),
 maxZ(         string("maxZ")             ),
 minR(         string("minR")             ),
 maxR(         string("maxR")             ),
 skipServices(        "skipServices"      , parsedOnly(), false), // broken, do not use
 m_innerRadius(       "innerRadius"       , parsedAndChecked()),
 m_outerRadius(       "outerRadius"       , parsedAndChecked()),
 m_innerRadiusFixed(  "innerRadiusFixed"  , parsedAndChecked(), true),
 m_outerRadiusFixed(  "outerRadiusFixed"  , parsedAndChecked(), true),
 m_sameRods(          "sameRods"          , parsedAndChecked(), false),
 m_barrelRotation(    "barrelRotation"    , parsedOnly(), 0.),
 m_layerNode(         "Layer"             , parsedOnly()),
 m_supportNode(       "Support"           , parsedOnly())
{
  // Set the geometry config parameters
  this->myid(name);
  this->store(treeProperty);
  this->store(nodeProperty);
}

//
// Limit barrel geometry by eta cut
//
void Barrel::cutAtEta(double eta)
{
  for (auto& l : m_layers) l.cutAtEta(eta);
  m_layers.erase_if([](const Layer& l) { return l.numRods() == 0; });
  numLayers(m_layers.size());
}

//
// Build recursively individual subdetector systems: Layers -> rods -> modules
//
void Barrel::build()
{
  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    // Build layers
    for (auto i=1; i <= numLayers(); i++) {

      Layer* layer = GeometryFactory::make<Layer>(i, numLayers(), m_sameRods(), m_innerRadiusFixed(), m_outerRadiusFixed(), m_barrelRotation(), m_layerNode, propertyTree());
      layer->build(numLayers(), m_innerRadius(), m_outerRadius());

      m_layers.push_back(layer);
    }

  }
  catch (PathfulException& pe) {

    pe.pushPath(fullid(*this));
    throw;
  }

  // Supports defined within a Barrel
  for (auto& mapel : m_supportNode) {
    SupportStructure* s = new SupportStructure();
    s->store(propertyTree());
    s->store(mapel.second);
    s->buildInBarrel(*this);
    m_supportStructures.push_back(s);
  }

  cleanup();
  builtok(true);
}

//
// Setup: link lambda functions to various barrel related properties (use setup functions for ReadOnly Computable properties)
//
void Barrel::setup()
{
  maxR.setup([&]()       { double max = 0;                                  for (const auto& l : m_layers) { max = MAX(max, l.maxR()); }       return max; });
  minR.setup([&]()       { double min = std::numeric_limits<double>::max(); for (const auto& l : m_layers) { min = MIN(min, l.minR()); }       return min; });
  maxRAllMat.setup([&]() { double max = 0;                                  for (const auto& l : m_layers) { max = MAX(max, l.maxRAllMat()); } return max; });
  minRAllMat.setup([&]() { double min = std::numeric_limits<double>::max(); for (const auto& l : m_layers) { min = MIN(min, l.minRAllMat()); } return min; });
  maxZ.setup([&]()       { double max =-std::numeric_limits<double>::max(); for (const auto& l : m_layers) { max = MAX(max, l.maxZ()); }       return max; });
  minZ.setup([&]()       { double min = std::numeric_limits<double>::max(); for (const auto& l : m_layers) { min = MIN(min, l.minZ()); }       return min; });
}

//
// Cross-check parameters provdied from geometry configuration file
//
void Barrel::check() {
  PropertyObject::check();

  if (m_sameRods() && (!m_innerRadiusFixed() || !m_outerRadiusFixed())) throw PathfulException("Rod building algorithm assumes that inner&outer barrel radii are fixed if same rods required to be built within whole barrel!");
}

//
// Add barrel service line
//
void Barrel::addServiceLine(InactiveElement* service)
{
  m_services.push_back(service);
}

//
// GeometryVisitor pattern -> barrel visitable
//
void Barrel::accept(GeometryVisitor& v)
{
  v.visit(*this);
  for (auto& l : m_layers) { l.accept(v); }
  for (auto& s : m_supportStructures) { s.accept(v); }
}

//
// GeometryVisitor pattern -> tracker visitable (const. option)
//
void Barrel::accept(ConstGeometryVisitor& v) const
{
  v.visit(*this);
  for (const auto& l : m_layers) { l.accept(v); }
  for (const auto& s : m_supportStructures) { s.accept(v); }
}

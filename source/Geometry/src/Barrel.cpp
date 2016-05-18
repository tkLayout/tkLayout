#include "Barrel.h"

#include "Layer.h"
#include "MessageLogger.h"
#include "SupportStructure.h"

using material::SupportStructure;

//
// Constructor - parse geometry config file using boost property tree & read-in Layer, Support nodes
//
Barrel::Barrel(const std::string& name, const PropertyTree& nodeProperty, const PropertyTree& treeProperty) :
 numLayers(           "numLayers"         , parsedAndChecked()),
 skipServices(        "skipServices"      , parsedOnly(), false), // broken, do not use
 m_innerRadius(       "innerRadius"       , parsedAndChecked()),
 m_outerRadius(       "outerRadius"       , parsedAndChecked()),
 m_innerRadiusFixed(  "innerRadiusFixed"  , parsedAndChecked(), true),
 m_outerRadiusFixed(  "outerRadiusFixed"  , parsedAndChecked(), true),
 m_sameRods(          "sameRods"          , parsedAndChecked(), false),
 m_barrelRotation(    "barrelRotation"    , parsedOnly(), 0.),
 m_supportMarginOuter("supportMarginOuter", parsedOnly(), 2.),
 m_supportMarginInner("supportMarginInner", parsedOnly(), 2.),
 m_layerNode(         "Layer"             , parsedOnly()),
 m_supportNode(       "Support"           , parsedOnly())
{
  // Set the geometry config parameters
  this->myid(name);
  this->store(treeProperty);
  this->store(nodeProperty);

  // Build & setup tracker
  this->build();
  this->setup();
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
// Build recursively individual subdetector systems: Layers -> rods -> modules -> private method called by constructor
//
void Barrel::build()
{
  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    for (int i = 1; i <= numLayers(); i++) {
      Layer* layer = GeometryFactory::make<Layer>();
      layer->myid(i);

      if      (i == 1)           { if (m_innerRadiusFixed()) layer->radiusMode(Layer::FIXED); layer->placeRadiusHint(m_innerRadius()); }
      else if (i == numLayers()) { if (m_outerRadiusFixed()) layer->radiusMode(Layer::FIXED); layer->placeRadiusHint(m_outerRadius()); }
      else                       { layer->placeRadiusHint(m_innerRadius() + (m_outerRadius()-m_innerRadius())/(numLayers()-1)*(i-1)); }

      if (m_sameRods()) {
        layer->minBuildRadius(m_innerRadius());
        layer->maxBuildRadius(m_outerRadius());
        layer->sameParityRods(true);
      }

      layer->store(propertyTree());
      if (m_layerNode.count(i) > 0) layer->store(m_layerNode.at(i));
      layer->build();
      layer->rotateZ(m_barrelRotation());
      layer->rotateZ(layer->layerRotation());
      m_layers.push_back(layer);
    }

  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

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
// Calculate various barrel related properties -> private method called by constructor
//
void Barrel::setup()
{
  maxR.setup([this]() { double max = 0;                                  for (const auto& l : m_layers) { max = MAX(max, l.maxR()); } return max; });
  minR.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& l : m_layers) { min = MIN(min, l.minR()); } return min; });
  maxZ.setup([this]() { double max =-std::numeric_limits<double>::max(); for (const auto& l : m_layers) { max = MAX(max, l.maxZ()); } return max; });
  minZ.setup([this]() { double min = std::numeric_limits<double>::max(); for (const auto& l : m_layers) { min = MIN(min, l.minZ()); } return min; });
}

//
// GeometryVisitor pattern -> barrel visitable
//
void Barrel::accept(GeometryVisitor& v)
{
  v.visit(*this);
  for (auto& l : m_layers) { l.accept(v); }
}

//
// GeometryVisitor pattern -> tracker visitable (const. option)
//
void Barrel::accept(ConstGeometryVisitor& v) const
{
  v.visit(*this);
  for (const auto& l : m_layers) { l.accept(v); }
}

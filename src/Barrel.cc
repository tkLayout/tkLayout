#include "Barrel.hh"
#include "MessageLogger.hh"
#include "SupportStructure.hh"

using material::SupportStructure;

void Barrel::cutAtEta(double eta) { 
  for (auto& l : layers_) l.cutAtEta(eta); 
  layers_.erase_if([](const Layer& l) { return l.numRods() == 0; });
  numLayers(layers_.size()); 
}

void Barrel::build() {
  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    for (int i = 1; i <= numLayers(); i++) {
      Layer* layer = GeometryFactory::make<Layer>(myid());
      layer->myid(i);

      // Compute layer's radii
      if      (i == 1)           { if (innerRadiusFixed()) layer->radiusMode(Layer::FIXED); layer->placeRadiusHint(innerRadius()); } 
      else if (i == numLayers()) { if (outerRadiusFixed()) layer->radiusMode(Layer::FIXED); layer->placeRadiusHint(outerRadius()); } 
      else                       { layer->placeRadiusHint(innerRadius() + (outerRadius()-innerRadius())/(numLayers()-1)*(i-1)); }
      if (sameRods()) { 
        layer->minBuildRadius(innerRadius()); 
        layer->maxBuildRadius(outerRadius()); 
        layer->sameParityRods(true);
      }

      // Store all layer properties
      layer->store(propertyTree());
      if (layerNode.count(i) > 0) layer->store(layerNode.at(i)); // TO DO: WARNING!! layer->placeRadiusHint is reassigned here!!

      // Build the layer
      layer->build();

      // Additional barrel rotation
      layer->rotateZ(barrelRotation());
      if (rotateBarrelByHalfPi()) { layer->rotateZ(M_PI / 2.); }

      // Additional layer rotation
      layer->rotateZ(layer->layerRotation());

      // Store the layer
      layers_.push_back(layer);
    }

  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  // Supports defined within a Barrel
  for (auto& mapel : supportNode) {
    SupportStructure* s = new SupportStructure();
    s->store(propertyTree());
    s->store(mapel.second);
    s->buildInBarrel(*this);
    supportStructures_.push_back(s);
  }

  cleanup();
  builtok(true);
}


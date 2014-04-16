#include "Barrel.h"

void Barrel::cutAtEta(double eta) { 
  for (auto& l : layers_) l.cutAtEta(eta); 
  layers_.erase_if([](const Layer& l) { return l.numRods() == 0; });
  numLayers(layers_.size()); 
}

void Barrel::build() {
  try {
    std::cout << ">>> Building " << fullid(*this) << " <<<" << std::endl;
    check();

    for (int i = 1; i <= numLayers(); i++) {
      Layer* layer = GeometryFactory::make<Layer>();
      layer->myid(i);

      if      (i == 1)           { layer->radiusMode(Layer::FIXED); layer->placeRadiusHint(innerRadius()); } 
      else if (i == numLayers()) { layer->radiusMode(Layer::FIXED); layer->placeRadiusHint(outerRadius()); } 
      else                       { layer->placeRadiusHint(innerRadius() + (outerRadius()-innerRadius())/(numLayers()-1)*(i-1)); }

      if (sameRods()) { 
        layer->minBuildRadius(innerRadius()); 
        layer->maxBuildRadius(outerRadius()); 
        layer->sameParityRods(true);
      }

      layer->store(propertyTree());
      if (layerNode.count(i) > 0) layer->store(layerNode.at(i));
      layer->build();
      layer->rotateZ(barrelRotation());
      layers_.push_back(layer);
    }

  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  cleanup();
  builtok(true);
}

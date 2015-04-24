/**
 * @file SupportStructure.cpp
 *
 * @date 25/Nov/2014
 * @author Stefano Martina
 */

#include <set>
#include "SupportStructure.h"
#include "MaterialTab.h"
#include "messageLogger.h"
#include "MaterialProperties.h"
#include "InactiveElement.h"
#include "InactiveTube.h"
#include "InactiveRing.h"
#include "InactiveSurfaces.h"
#include "Barrel.h"

using insur::InactiveTube;
using insur::InactiveRing;

namespace material {
  //=============== begin class SupportStructure
  const std::map<std::string, SupportStructure::Type> SupportStructure::typeStringMap = {
    {"custom", CUSTOM},
    {"auto", AUTO},
    {"top", TOP},
    {"bottom", BOTTOM}
  };
  const std::map<std::string, SupportStructure::Direction> SupportStructure::directionStringMap = {
    {"horizontal", HORIZONTAL},
    {"vertical", VERTICAL}
  };


  SupportStructure::SupportStructure() :
    componentsNode("Component", parsedOnly()),
    type("type", parsedAndChecked()),
    autoPosition("autoPosition", parsedOnly()),
    customZMin("customZMin", parsedOnly()),
    customRMin("customRMin", parsedOnly()),
    customLength("customLength", parsedOnly()),
    customDir("customDir", parsedOnly())
  {}
  
  void SupportStructure::buildInTracker() {
    InactiveElement* inactiveElement;

    buildBase();

    try {
      supportType_ = typeStringMap.at(type());
    }  catch (const std::out_of_range& ex) {
      logERROR("Unrecognized value " + type() + ".");
      return;
    }

    if (supportType_ == CUSTOM) {
      if (customZMin.state() && customRMin.state() && customLength.state() && customDir.state()) {
        try {
          direction_ = directionStringMap.at(customDir());
        }  catch (const std::out_of_range& ex) {
          logERROR("Unrecognized value " + customDir() + ".");
          return;
        }

        buildInactiveElementPair(direction_, customZMin(), customRMin(), customLength());
      } else {
        logERROR("Property customZMin, customRMin, customLength, or customDir not set.");
        return;
      }
    } else {
      logERROR("Support type \"" + type() + "\" defined in wrong place.");
      return;
    }

    cleanup();
  }

  void SupportStructure::buildInBarrel(Barrel& barrel) {
    InactiveElement* inactiveElement;

    buildBase();

    try {
      supportType_ = typeStringMap.at(type());
    }  catch (const std::out_of_range& ex) {
      logERROR("Unrecognized value " + type() + ".");
      return;
    }

    switch(supportType_) {
    case AUTO :
      if (autoPosition.state()) {
        direction_ = VERTICAL; //AUTO supports are only for barrels, and vertical

        std::set<double> layerRadiuses; //use set for ordering safe check

        // class LayerVisitor : public ConstGeometryVisitor {
        // private:
        //   std::set<double>& layerRadiuses_;
        //   const double& autoLayerMarginUpper_;
        //   const double& autoLayerMarginLower_;
          
        // public:
        //   LayerVisitor(std::set<double>& layerRadiuses, const double& autoLayerMarginUpper, const double& autoLayerMarginLower) : 
        //     layerRadiuses_(layerRadiuses),
        //     autoLayerMarginUpper_(autoLayerMarginUpper),
        //     autoLayerMarginLower_(autoLayerMarginLower)
        //   {};
          
        //   void visit(Layer& layer) {
        //     layerRadiuses_.insert(layer.minR() - autoLayerMarginUpper_);
        //     layerRadiuses_.insert(layer.maxR() + autoLayerMarginLower_);
        //   }
        // };

        // LayerVisitor visitor(layerRadiuses, autoLayerMarginUpper, autoLayerMarginLower);
        // barrel.accept(visitor);

        //get radiuses around layers
        for(const Layer& layer : barrel.layers()) {
          layerRadiuses.insert(layer.minR() - autoLayerMarginUpper);
          layerRadiuses.insert(layer.maxR() + autoLayerMarginLower);
        }

        if(layerRadiuses.size() < 4) {
          logERROR("Barrel with only one or zero layers. Auto support impossible to build.");
          return;
        }

        //delete minimum and maximum useless radiuses (support structure only in the inner spaces)
        layerRadiuses.erase(layerRadiuses.begin());
        layerRadiuses.erase(std::prev(layerRadiuses.end()));

        //build inactiveElements inside spaces
        for(std::set<double>::iterator minIter = layerRadiuses.begin(), maxIter = ++ layerRadiuses.begin(); minIter != layerRadiuses.end(); std::advance(minIter,2), std::advance(maxIter,2)) {
          buildInactiveElementPair(direction_, autoPosition(), *minIter, *maxIter - *minIter);
        }
      } else {
        logERROR("Property autoPosition not set.");
        return;
      }
      break;

    case TOP :
      direction_ = HORIZONTAL; //TOP and BOTTOM supports are only for barrels, and horizontal
      buildInactiveElementPair(direction_,
                               MAX(0, barrel.minZ()),
                               barrel.maxR() + autoLayerMarginLower,
                               barrel.maxZ() - MAX(0, barrel.minZ()));
      break;

    case BOTTOM :
      direction_ = HORIZONTAL;
      buildInactiveElementPair(direction_,
                               MAX(0, barrel.minZ()),
                               barrel.minR() - inactiveElementWidth - autoLayerMarginUpper,
                               barrel.maxZ() - MAX(0, barrel.minZ()));
      break;

    default :
      logERROR("Support type \"" + type() + "\" defined in wrong place.");
      return;
    }

    cleanup();
  }

  void SupportStructure::updateInactiveSurfaces(InactiveSurfaces& inactiveSurfaces) {
    for(InactiveElement* inactiveElement : inactiveElements) {
      inactiveSurfaces.addSupportPart(*inactiveElement);
    }
  }

  void SupportStructure::buildBase() {
    for (auto& currentComponentNode : componentsNode) {
      Component* newComponent = new Component();
      newComponent->store(propertyTree());
      newComponent->store(currentComponentNode.second);
      newComponent->check();
      newComponent->build();

      components_.push_back(newComponent);
    }
  }

  void SupportStructure::populateMaterialProperties(MaterialProperties& materialProperties) const {
    for (const Component* currComponent : components_) {
      currComponent->populateMaterialProperties(materialProperties);
    }
  }

  void SupportStructure::buildInactiveElementPair(Direction direction, double startZ, double startR, double length) {
    InactiveElement* zPositiveElement;
    InactiveElement* zNegativeElement;

    if(direction == HORIZONTAL) {
      zPositiveElement = new InactiveTube;
      zPositiveElement->setZLength(length);
      zPositiveElement->setZOffset(startZ);
      zPositiveElement->setInnerRadius(startR);
      zPositiveElement->setRWidth(inactiveElementWidth);
      zPositiveElement->setFinal(true);
      //zPositiveElement->setCategory(insur::MaterialProperties::o_sup);      
      zPositiveElement->setCategory(insur::MaterialProperties::b_sup);      

      zNegativeElement = new InactiveTube;
      zNegativeElement->setZLength(length);
      zNegativeElement->setZOffset(-1 * startZ - length);
      zNegativeElement->setInnerRadius(startR);
      zNegativeElement->setRWidth(inactiveElementWidth);
      zNegativeElement->setFinal(true);
      //zNegativeElement->setCategory(insur::MaterialProperties::o_sup);      
      zNegativeElement->setCategory(insur::MaterialProperties::b_sup);      
    } else {
      zPositiveElement = new InactiveRing;
      zPositiveElement->setZLength(inactiveElementWidth);
      zPositiveElement->setZOffset(startZ);
      zPositiveElement->setInnerRadius(startR);
      zPositiveElement->setRWidth(length);
      zPositiveElement->setFinal(true);
      //zPositiveElement->setCategory(insur::MaterialProperties::u_sup);
      zPositiveElement->setCategory(insur::MaterialProperties::b_sup);

      zNegativeElement = new InactiveRing;
      zNegativeElement->setZLength(inactiveElementWidth);
      zNegativeElement->setZOffset(-1 * startZ - inactiveElementWidth);
      zNegativeElement->setInnerRadius(startR);
      zNegativeElement->setRWidth(length);
      zNegativeElement->setFinal(true);
      //zNegativeElement->setCategory(insur::MaterialProperties::u_sup);
      zNegativeElement->setCategory(insur::MaterialProperties::b_sup);
    }

      populateMaterialProperties(*zPositiveElement);
      populateMaterialProperties(*zNegativeElement);
      inactiveElements.push_back(zPositiveElement);
      inactiveElements.push_back(zNegativeElement);
  }


  //=============== end class SupportStructure

  //=============== begin class SupportStructure::Component
  SupportStructure::Component::Component() :
    componentsNode("Component", parsedOnly()),
    elementsNode("Element", parsedOnly()) {}

  void SupportStructure::Component::build() {
    //sub components
    for (auto& currentComponentNode : componentsNode) {
      Component* newComponent = new Component();
      newComponent->store(propertyTree());
      newComponent->store(currentComponentNode.second);
      newComponent->check();
      newComponent->build();

      components_.push_back(newComponent);
    }
    //elements
    for (auto& currentElementNode : elementsNode) {
      Element* newElement = new Element();
      newElement->store(propertyTree());
      newElement->store(currentElementNode.second);
      newElement->check();
      newElement->cleanup();

      elements_.push_back(newElement);
    }
    cleanup();
  }
  void SupportStructure::Component::populateMaterialProperties(MaterialProperties& materialProperties) const {
    for (const Component* currComponent : components_) {
      currComponent->populateMaterialProperties(materialProperties);
    }
    for (const Element* currElement : elements_) {
      currElement->populateMaterialProperties(materialProperties);
    }
  }

  //=============== end class SupportStructure::Component
  
  //=============== begin class SupportStructure::Element
  SupportStructure::Element::Element() :
    componentName ("componentName", parsedOnly()),
    elementName ("elementName", parsedAndChecked()),
    quantity ("quantity", parsedAndChecked()),
    unit ("unit", parsedAndChecked()),
    debugInactivate ("debugInactivate", parsedOnly(), false),
    materialTab_ (MaterialTab::instance()) {}
    
  const std::string SupportStructure::Element::msg_no_valid_unit = "No valid unit: ";

  const std::map<std::string, SupportStructure::Element::Unit> SupportStructure::Element::unitStringMap = {
      {"g", GRAMS},
      {"mm", MILLIMETERS},
      {"g/m", GRAMS_METER}
  };

  double SupportStructure::Element::quantityInGrams(double length, double surface) const {
    double returnVal;
    try {
      switch (unitStringMap.at(unit())) {
      case Element::GRAMS:
        returnVal = quantity();
        break;

      case Element::GRAMS_METER:
        returnVal = length * quantity() / 1000.0;
        break;

      case Element::MILLIMETERS:
        std::string elementNameString = elementName();
        double elementDensity = materialTab_.density(elementNameString);
        returnVal = elementDensity * surface * quantity();
        break;
      }
    } catch (const std::out_of_range& ex) {
      logERROR(msg_no_valid_unit + unit());
    }

    return returnVal;
  }

  void SupportStructure::Element::populateMaterialProperties(MaterialProperties& materialProperties) const {
    double quantity;
    
    if(debugInactivate() == false) {
      quantity = quantityInGrams(materialProperties.getLength(), materialProperties.getSurface());
      materialProperties.addLocalMass(elementName(), componentName(), quantity);
    }
  }
  
  //=============== end class SupportStructure::Element
}


/**
 * @file SupportStructure.cpp
 *
 * @date 25/Nov/2014
 * @author Stefano Martina
 */

#include <set>
#include "SupportStructure.h"

#include "MaterialTab.h"
#include "MessageLogger.h"
#include "MaterialProperties.h"
#include "InactiveElement.h"
#include "InactiveTube.h"
#include "InactiveRing.h"
#include "Barrel.h"
#include "Endcap.h"
#include "Units.h"

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
    componentsNode("Component"   , parsedOnly()),
    type(          "type"        , parsedAndChecked()),
    autoPosition(  "autoPosition", parsedOnly()),
    customZMin(    "customZMin"  , parsedOnly()),
    customRMin(    "customRMin"  , parsedOnly()),
    customLength(  "customLength", parsedOnly()),
    customDir(     "customDir"   , parsedOnly())
  {}
  
  void SupportStructure::buildInTracker() {
    //InactiveElement* inactiveElement;

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
        logINFO("Building custom support structure, positioned at Z: "+any2str(customZMin())+ " R: "+any2str(customRMin())+" of length: "+any2str(customLength()));
      } else {
        logERROR("Property customZMin, customRMin, customLength, or customDir not set.");
        return;
      }
    } else {
      logERROR("Support type \"" + type() + "\" not supported at the Tracker level.");
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

        std::set<double> layerRadii; //use set for ordering safe check

         /* class LayerVisitor : public ConstGeometryVisitor {
         private:
           std::set<double>& layerRadiuses_;
           const double& autoLayerMarginUpper_;
           const double& autoLayerMarginLower_;
          
         public:
           LayerVisitor(std::set<double>& layerRadiuses, const double& autoLayerMarginUpper, const double& autoLayerMarginLower) :
             layerRadiuses_(layerRadiuses),
             autoLayerMarginUpper_(autoLayerMarginUpper),
             autoLayerMarginLower_(autoLayerMarginLower)
           {};
          
           void visit(Layer& layer) {
             layerRadiuses_.insert(layer.minR() - autoLayerMarginUpper_);
             layerRadiuses_.insert(layer.maxR() + autoLayerMarginLower_);
           }
         };

         LayerVisitor visitor(layerRadiuses, autoLayerMarginUpper, autoLayerMarginLower);
         barrel.accept(visitor); */

        //get radiuses around layers
        for(const Layer& layer : barrel.layers()) {
          layerRadii.insert(layer.minR() - autoLayerMarginUpper);
          layerRadii.insert(layer.maxR() + autoLayerMarginLower);
        }
        if(layerRadii.size() < 4) {
          logERROR("Barrel with only one or zero layers. Auto support impossible to build.");
          return;
        }

        //delete minimum and maximum useless radiuses (support structure only in the inner spaces)
        layerRadii.erase(layerRadii.begin());
        layerRadii.erase(std::prev(layerRadii.end()));

        //build inactiveElements inside spaces
        for(std::set<double>::iterator minIter = layerRadii.begin(), maxIter = ++ layerRadii.begin(); minIter != layerRadii.end(); std::advance(minIter,2), std::advance(maxIter,2)) {
          buildInactiveElementPair(direction_, autoPosition(), *minIter, *maxIter - *minIter);
        }
        logINFO("Building barrel support structure vertically oriented, positioned at Z: "+any2str(autoPosition()));
      } else {
        logERROR("Property autoPosition not set.");
        return;
      }
      break;
    case TOP :
      direction_ = HORIZONTAL;
      buildInactiveElementPair(direction_,
                               MAX(0, barrel.minZ()),
                               barrel.maxR() + autoLayerMarginLower,
                               barrel.maxZ() - MAX(0, barrel.minZ()));
      //std::cout << ">>Top barrel support>> " << "maxR: " << barrel.maxR() + autoLayerMarginLower << " minZ: " << MAX(0, barrel.minZ()) << " maxZ: " << barrel.maxZ() << std::endl;
      logINFO("Building barrel top support structure horizontally oriented");
      break;
    case BOTTOM :
      direction_ = HORIZONTAL;
      buildInactiveElementPair(direction_,
                               MAX(0, barrel.minZ()),
                               barrel.minR() - inactiveElementWidth - autoLayerMarginUpper,
                               barrel.maxZ() - MAX(0, barrel.minZ()));
      //std::cout << ">>Bottom barrel support>> " << "minR: " << barrel.minR() - inactiveElementWidth - autoLayerMarginUpper << " minZ: " << MAX(0, barrel.minZ()) << " maxZ: " << barrel.maxZ() << std::endl;
      logINFO("Building barrel bottom support structure horizontally oriented");
      break;
    default :
      logERROR("Support type \"" + type() + "\" not supported at the barrel level.");
      return;
    }
    cleanup();
  }

  void SupportStructure::buildInEndcap(Endcap& endcap) {

    InactiveElement* inactiveElement;
    buildBase();

    try {
      supportType_ = typeStringMap.at(type());
    }  catch (const std::out_of_range& ex) {
      logERROR("Unrecognized value " + type() + ".");
      return;
    }

    switch(supportType_) {
    case TOP :

      direction_ = HORIZONTAL;
      buildInactiveElementPair(direction_,
                               MAX(0, endcap.minZ()),
                               endcap.maxR() + autoLayerMarginLower,
                               endcap.maxZ() - MAX(0, endcap.minZ()));
      //std::cout << ">>Top end-cap support>> " << "maxR: " << endcap.maxR() + geom_support_margin_top << " minZ: " << MAX(0, endcap.minZ()) << " maxZ: " << endcap.maxZ() << std::endl;
      logINFO("Building end-cap top support structure horizontally oriented");
      break;
    case BOTTOM :

      direction_ = HORIZONTAL;
      buildInactiveElementPair(direction_,
                               MAX(0, endcap.minZ()),
                               endcap.minR() - inactiveElementWidth - autoLayerMarginUpper,
                               endcap.maxZ() - MAX(0, endcap.minZ()));
      //std::cout << ">>Bottom end-cap support>> " << "minR: " << endcap.minR() - geom_support_margin_bottom << " minZ: " << MAX(0, endcap.minZ()) << " maxZ: " << endcap.maxZ() << std::endl;
      logINFO("Building end-cap bottom support structure horizontally oriented");
      break;
    default :
      logERROR("Support type \"" + type() + "\" not supported at the end-cap level.");
      return;
    }
    cleanup();
  }

  //
  // GeometryVisitor pattern -> support structure visitable
  //
  void SupportStructure::accept(GeometryVisitor& v)
  {
    v.visit(*this);
  }

  //
  // GeometryVisitor pattern -> support structure visitable (const. option)
  //
  void SupportStructure::accept(ConstGeometryVisitor& v) const {
    v.visit(*this);
  }

//  void SupportStructure::updateInactiveSurfaces(InactiveSurfaces& inactiveSurfaces) {
//    for(auto& inactiveElement : inactiveElements) {
//      inactiveSurfaces.addSupportPart(inactiveElement);
//    }
//  }

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
      zPositiveElement = new InactiveTube(startZ,length,startR,inactiveElementWidth);
      //zPositiveElement->setCategory(MaterialProperties::o_sup);
      zPositiveElement->setCategory(MaterialProperties::b_sup);

      zNegativeElement = new InactiveTube(-1*startZ-length,length,startR,inactiveElementWidth);
      //zNegativeElement->setCategory(MaterialProperties::o_sup);
      zNegativeElement->setCategory(MaterialProperties::b_sup);
    } else {
      zPositiveElement = new InactiveRing(startZ,inactiveElementWidth,startR,length);
      //zPositiveElement->setCategory(MaterialProperties::u_sup);
      zPositiveElement->setCategory(MaterialProperties::b_sup);

      zNegativeElement = new InactiveRing(-1*startZ-inactiveElementWidth,inactiveElementWidth,startR,length);
      //zNegativeElement->setCategory(MaterialProperties::u_sup);
      zNegativeElement->setCategory(MaterialProperties::b_sup);
    }

      populateMaterialProperties(*zPositiveElement);
      populateMaterialProperties(*zNegativeElement);
      m_inactiveElements.push_back(zPositiveElement);
      m_inactiveElements.push_back(zNegativeElement);
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
    m_materialTab (MaterialTab::getInstance()) {}
    
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
        returnVal = length * quantity() * Units::g/Units::m;
        break;

      case Element::MILLIMETERS:
        std::string elementNameString = elementName();
        double elementDensity = m_materialTab.density(elementNameString);
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

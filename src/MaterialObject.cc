/**
 * @file MaterialObject.cc
 *
 * @date 19/jun/2014
 * @author Stefano Martina
 */

//#include "Materialway.hh"
#include "MaterialObject.hh"
#include "ConversionStation.hh"
#include "global_constants.hh"
#include "MaterialTab.hh"
//#include "InactiveElement.hh"
#include "MaterialProperties.hh"
#include "DetectorModule.hh"
#include "MessageLogger.hh"
#include <stdexcept>


namespace material {
  const std::map<MaterialObject::Type, const std::string> MaterialObject::typeString = {
      {MODULE, "module"},
      {ROD, "rod"},
      {LAYER, "layer"}
  };

  MaterialObject::MaterialObject(Type materialType) :
      materialType_ (materialType),
      type_ ("type", parsedOnly()),
      destination_ ("destination", parsedOnly()),
      debugInactivate_ ("debugInactivate", parsedOnly(), false),
      materialsNode_ ("Materials", parsedOnly()),
      // sensorNode_ ("Sensor", parsedOnly()),
      materials_ (nullptr) {}

  MaterialObject::MaterialObject(const MaterialObject& other) :
    MaterialObject(other.materialType_) {
    materials_ = other.materials_;
    serviceElements_ = other.serviceElements_; //do shallow copies
  }

  MaterialObject::~MaterialObject() {
    // if (materials_ != nullptr) {
    //   delete materials_;
    //   materials_ = nullptr;
    // }
  }

  const std::string MaterialObject::getTypeString() const {
    auto mapIter = typeString.find(materialType_);
    if (mapIter != typeString.end()) {
      return mapIter->second;
    } else {
      return "";
    }
  }

  double MaterialObject::totalGrams(double length, double surface) const {
    double result = 0.0;
    for (const Element* currElement : serviceElements_) {
      result += currElement->totalGrams(length, surface);
    }
    if (materials_ != nullptr) {
      result += materials_->totalGrams(length, surface);
    }
    return result;
  }

  void MaterialObject::build() {
    check();
    if (!debugInactivate_()) {
      // for (auto& currentSensor : sensorNode_) {
      //   ReferenceSensor temporarySensor;
      //   temporarySensor.store(currentSensor.second);
      //   temporarySensor.check();
      //   temporarySensor.cleanup();

      //   std::cout << "[" << currentSensor.first << "]=" << temporarySensor.numChannels() << "; ";
      //   sensorChannels[currentSensor.first] = temporarySensor.numChannels();
      // }
      // std::cout << "}" << std::endl;
      

      static std::map<MaterialObjectKey, Materials*> materialsMap_; //for saving memory
      for (auto& currentMaterialNode : materialsNode_) {
        store(currentMaterialNode.second);

        check();
        if (type_().compare(getTypeString()) == 0) {
          MaterialObjectKey myKey(currentMaterialNode.first, sensorChannels, destination_.state()? destination_() : std::string(""));
          if (materialsMap_.count(myKey) == 0) {
            Materials * newMaterials  = new Materials(materialType_);
            newMaterials->store(currentMaterialNode.second);

            //pass destination to newMaterials
            if(destination_.state()) {
              PropertyTree destinationPt;
              destinationPt.add(destination_.name(), destination_());
              newMaterials->store(destinationPt);
            }

            newMaterials->build(sensorChannels);
            materialsMap_[myKey] = newMaterials;
          }
          materials_ = materialsMap_[myKey];

          break;
        }
      }

    }
    cleanup();
  }

  void MaterialObject::deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices /*= false */, double gramsMultiplier /*= 1.*/) const {
    for(const Element * currElement : serviceElements_) {
      currElement->deployMaterialTo(outputObject, unitsToDeploy, onlyServices, gramsMultiplier);
    }
    
    if (materials_ != nullptr) {
      materials_->deployMaterialTo(outputObject, unitsToDeploy, onlyServices, gramsMultiplier);
    }    
  }

  void MaterialObject::addElement(const MaterialObject::Element* element) {
    if(element != nullptr) {
      serviceElements_.push_back(element);
    }
  }

  void MaterialObject::populateMaterialProperties(MaterialProperties& materialProperties) const {
    double quantity = 0;
    
    for (const Element* currElement : serviceElements_) {
      //currElement.populateMaterialProperties(materialProperties);
      //populate directly because need to skip the control if is a service
      //TODO: check why componentName is not present in no Element
      
      if (currElement->debugInactivate() == false) {
        quantity = currElement->totalGrams(materialProperties);

        if (currElement->componentName.state()) {
	  /*if (currElement->componentName() == "Sensor HV line") {
	    std::cout << "currElement->componentName()" << currElement->componentName() << "currElement->elementName() = " << currElement->elementName() << "quantity = " << quantity << std::endl;
	    }*/
          materialProperties.addLocalMass(currElement->elementName(), currElement->componentName(), quantity);
        } else {
          materialProperties.addLocalMass(currElement->elementName(), quantity);
        }
      }
    }

    if (materials_ != nullptr) {
      materials_->populateMaterialProperties(materialProperties);
    }
  }

  ElementsVector& MaterialObject::getLocalElements() const {
    ElementsVector* elementsList = new ElementsVector;
    if (materials_ != nullptr) {
      materials_->getLocalElements(*elementsList);
    }

    return *elementsList;
  }

  bool MaterialObject::isPopulated() const {
    return (materials_ != nullptr);
  }


  //void MaterialObject::chargeTrain(Materialway::Train& train) const {
  //  materials_->chargeTrain(train);
  //}

  MaterialObject::Materials::Materials(MaterialObject::Type newMaterialType) :
    componentsNode_ ("Component", parsedOnly()),
    materialType_(newMaterialType) {};

  MaterialObject::Materials::~Materials() {}

  double MaterialObject::Materials::totalGrams(double length, double surface) const {
    double result = 0.0;
    for (const Component* currentComponentNode : components_) {
      result += currentComponentNode->totalGrams(length, surface);
    }
    return result;
  }

  void MaterialObject::Materials::build(const std::map<int, int>& newSensorChannels) {
    check();        
    for (auto& currentComponentNode : componentsNode_) {
      Component* newComponent = new Component(materialType_);
      newComponent->store(propertyTree());
      newComponent->store(currentComponentNode.second);
      newComponent->check();
      newComponent->build(newSensorChannels);

      components_.push_back(newComponent);
    }
    cleanup();
  }

  void MaterialObject::Materials::deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices, double gramsMultiplier /*= 1.*/) const {
    for (const Component* currComponent : components_) {
      currComponent->deployMaterialTo(outputObject, unitsToDeploy, onlyServices, gramsMultiplier);
    }
  }

  void MaterialObject::Materials::populateMaterialProperties(MaterialProperties& materialProperties) const {
    for (const Component* currComponent : components_) {
      currComponent->populateMaterialProperties(materialProperties);
    }
  }

  void MaterialObject::Materials::getLocalElements(ElementsVector& elementsList) const {
    for (const Component* currComponent : components_) {
      currComponent->getLocalElements(elementsList);
    }
  }

  MaterialObject::Component::Component(MaterialObject::Type& newMaterialType) :
    //componentName ("componentName", parsedAndChecked()),
    componentsNode_ ("Component", parsedOnly()),
    elementsNode_ ("Element", parsedOnly()),
    materialType_(newMaterialType) {};

  MaterialObject::Component::~Component() { }

  double MaterialObject::Component::totalGrams(double length, double surface) const {
    double result = 0.0;
    for (Component* currentComponentNode : components_) {
      result += currentComponentNode->totalGrams(length, surface);
    }
    for  (const Element* currentElementNode : elements_) {
      result += currentElementNode->totalGrams(length, surface);
    }
    return result;
  }

  void MaterialObject::Component::build(const std::map<int, int>& newSensorChannels) {
    check();
    //std::cout << "COMPONENT " << componentName() << std::endl;

    //sub components
    for (auto& currentComponentNode : componentsNode_) {
      Component* newComponent = new Component(materialType_);
      newComponent->store(propertyTree());
      newComponent->store(currentComponentNode.second);
      newComponent->check();
      newComponent->build(newSensorChannels);

      components_.push_back(newComponent);
    }
    //elements
    for (auto& currentElementNode : elementsNode_) {
      Element* newElement = new Element(materialType_);
      newElement->store(propertyTree());
      newElement->store(currentElementNode.second);
      newElement->check();
      newElement->cleanup();
      newElement->build(newSensorChannels);
      //bool test1 = newElement->componentName.state();
      //bool test2 = newElement->nSegments.state();

      elements_.push_back(newElement);
    }
    cleanup();
  }

  void MaterialObject::Component::deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices, double gramsMultiplier /*= 1.*/) const {
    for(const Element* currElement : elements_) {
      currElement->deployMaterialTo(outputObject, unitsToDeploy, onlyServices, gramsMultiplier);
    }
    for (const Component* currComponent : components_) {
      currComponent->deployMaterialTo(outputObject, unitsToDeploy, onlyServices, gramsMultiplier);
    }
  }

  void MaterialObject::Component::populateMaterialProperties(MaterialProperties& materialProperties) const {
    for (const Component* currComponent : components_) {
      currComponent->populateMaterialProperties(materialProperties);
    }
    for (const Element* currElement : elements_) {
      currElement->populateMaterialProperties(materialProperties);
    }
  }

  void MaterialObject::Component::getLocalElements(ElementsVector& elementsList) const {
    for (const Component* currComponent : components_) {
      currComponent->getLocalElements(elementsList);
    }
    for (const Element* currElement : elements_) {
      currElement->getLocalElements(elementsList);
    }
  }


  /*
  const std::map<MaterialObject::Type, const std::string> MaterialObject::Element::unitString = {
      {GRAMS, "g"},
      {MILLIMETERS, "mm"},
      {GRAMS_METER, "gm"}
  };
  */

  MaterialObject::Element::Element(MaterialObject::Type& newMaterialType) :
    componentName ("componentName", parsedOnly()),
    //numStripsAcrossEstimate("numStripsAcrossEstimate", parsedOnly()),
    //numSegmentsEstimate("numSegmentsEstimate", parsedOnly()),
    //nStripsAcross("nStripsAcross", parsedOnly()),
    //nSegments("nSegments", parsedOnly()),
    referenceSensorNode ("ReferenceSensor", parsedOnly()),
    elementName ("elementName", parsedAndChecked()),
    service ("service", parsedOnly(), false),
    scaleOnSensor ("scaleOnSensor", parsedOnly(), false),
    quantity ("quantity", parsedAndChecked()),
    unit ("unit", parsedAndChecked()),
    debugInactivate ("debugInactivate", parsedOnly(), false),
    destination ("destination", parsedOnly()),
    targetVolume ("targetVolume", parsedOnly(), 0),
    materialTab_ (MaterialTab::instance()),
    materialType_(newMaterialType) {};

  MaterialObject::Element::Element(const Element& original, double multiplier) : Element(original.materialType_) {
    if(original.destination.state())
      destination(original.destination());
    if(original.componentName.state())
      componentName(original.componentName());
    elementName(original.elementName());
    service(original.service());
    quantity(original.quantity() * original.scalingMultiplier() * multiplier); //apply the scaling in the copied object
    scaleOnSensor(0);
    unit(original.unit());
    debugInactivate(original.debugInactivate());
  }
  
  MaterialObject::Element::~Element() { }

  const std::string MaterialObject::Element::msg_no_valid_unit = "No valid unit: ";

  const std::map<std::string, MaterialObject::Element::Unit> MaterialObject::Element::unitStringMap = {
      {"g", GRAMS},
      {"mm", MILLIMETERS},
      {"g/m", GRAMS_METER}
  };

  void MaterialObject::Element::deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices /*= false*/, double gramsMultiplier /*= 1.*/) const {
    const Element* elementToDeploy = this;
    bool valid = false;
    if ((! onlyServices) || (onlyServices && (service() == true))) {
      if((unit().compare("g") == 0) && (service() == true)) { 
        logERROR(err_service1 + elementName() + err_service2);
      } else {
        valid = true;
      }
      // if (materialType_ == STATION) {
      //   if(unit().compare("g") == 0) { 
      //     logERROR(err_service1 + elementName() + err_service2);
      //   } else {
      //     valid = true;
      //   }
      // } else if (materialType_ == ROD) {
      //   if((unit().compare("g") == 0) && (service() == true) { 
      //     logERROR(err_service1 + elementName() + err_service2);
      //   } else {
      //     valid = true;
      //   }
      // } else if (materialType_ == MODULE) {
      //   if (service() == true) {
      //     valid = true;
      //   }      
      // }
        
      if(valid) {
        for(const std::string& unitToDeploy : unitsToDeploy) {
          if (unit().compare(unitToDeploy) == 0) {
            if (((materialType_ == ROD) || (materialType_ == MODULE)) && (service()==true) && (unit().compare("mm") == 0)) {
              logUniqueWARNING("Definition of services in \"mm\" is deprecated");
            }
            if (unit().compare("g") == 0) {
              elementToDeploy = new Element(*this, gramsMultiplier);
            }
            outputObject.addElement(elementToDeploy);
            break;
          }
        }
      }
    }
  }

  double MaterialObject::Element::quantityInGrams(const DetectorModule& module) const {
    return quantityInUnit("g", module.length(), module.area());
  }

  double MaterialObject::Element::quantityInGrams(const MaterialProperties& materialProperties) const {
    return quantityInUnit("g", materialProperties.getLength(), materialProperties.getSurface());
  }

  double MaterialObject::Element::quantityInGrams(const double length, const double surface) const {
    return quantityInUnit("g", length, surface);
  }

  double MaterialObject::Element::quantityInUnit(const std::string desiredUnit, const MaterialProperties& materialProperties) const {
    return quantityInUnit(desiredUnit, materialProperties.getLength(), materialProperties.getSurface());
  }

  /**
   * return the quantity in the desired unit quantity
   * @param desiredUnit the desired unit, one between 'g', 'g/m', 'mm'
   * @param length the length in mm
   * @param surface the surface in mm^2
   */
     
  double MaterialObject::Element::quantityInUnit(const std::string desiredUnit, const double length, const double surface) const {
    double returnVal = 0;
    double density = materialTab_.density(elementName());
    bool invert;
    Unit desiredUnitVal, elementUnitVal, tempUnit;

    //Conversion matrix:
    //            g              g/m                 mm
    //      __________________________________________________
    //  g  |      1             l/1000            rho*S       |
    // g/m |    1000/l            1           (rho*S*1000)/l  |
    //  mm |   1/(rho*S)    l/(rho*S*1000)          1         |
    //
    // rows:    desired unit
    // columns: original unit
    // l:       length
    // rho:     density
    // S:       surface

    /*
    std::map<std::pair<Unit, Unit>, double> conversionMatrix = {
      {{GRAMS, GRAMS}, 1}, {{GRAMS, GRAMS_METER}, length/1000}, {{GRAMS, MILLIMETERS}, density*surface},
      {{GRAMS_METER, GRAMS}, 1000/length}, {{GRAMS_METER, GRAMS_METER}, 1}, {{GRAMS_METER, MILLIMETERS}, (density*surface*1000)/length},
      {{MILLIMETERS, GRAMS}, 1/(density*surface)}, {{MILLIMETERS, GRAMS_METER}, length/(density*surface*1000)}, {{MILLIMETERS, MILLIMETERS}, 1}
    };

    try {
      returnVal = quantity() * conversionMatrix.at({unitStringMap.at(desiredUnit), unitStringMap.at(unit())});
    } catch (const std::out_of_range& ex) {
      logERROR(msg_no_valid_unit + unit() + ", " + desiredUnit + ".");
    }
    */

    try {
      desiredUnitVal = unitStringMap.at(desiredUnit);
      elementUnitVal = unitStringMap.at(unit());
      
      if (desiredUnitVal == elementUnitVal) {
        return quantity();
      } else if (desiredUnitVal > elementUnitVal) {
        invert = true;
        tempUnit = desiredUnitVal;
        desiredUnitVal = elementUnitVal;
        elementUnitVal = tempUnit;
      } else {
        invert = false;
      }
      
      if      ((desiredUnitVal == GRAMS) && (elementUnitVal == GRAMS_METER))
        returnVal = quantity() * length / 1000.;
      else if ((desiredUnitVal == GRAMS) && (elementUnitVal == MILLIMETERS))
        returnVal = quantity() * density * surface;
      else if ((desiredUnitVal == GRAMS_METER) && (elementUnitVal == MILLIMETERS))
        returnVal = quantity() * (density * surface * 1000.) / length;

      if (invert)
        returnVal = 1 / returnVal;
    } catch (const std::out_of_range& ex) {
      logERROR(msg_no_valid_unit + unit() + ", " + desiredUnit + ".");
    }
    return returnVal;
  }

  double MaterialObject::Element::totalGrams(const DetectorModule& module) const {
    return totalGrams(module.length(), module.area());
  }

  double MaterialObject::Element::totalGrams(const MaterialProperties& materialProperties) const {
    return totalGrams(materialProperties.getLength(), materialProperties.getSurface());
  }
  
  double MaterialObject::Element::totalGrams(double length, double surface) const {
    return quantityInGrams(length, surface) * scalingMultiplier();
  }

  double MaterialObject::Element::scalingMultiplier() const {
    int sensorIndex = scaleOnSensor();

    if(sensorIndex != 0) {
      if(materialType_ == MODULE) {
        try {
          return sensorChannels_.at(sensorIndex) / referenceSensors_.at(sensorIndex)->numChannels();
        } catch (std::out_of_range& e) {
          std::stringstream error;
          error << "Sensor " << sensorIndex << " don't exists." << std::endl;
          logERROR(error.str());
        }
      } else {
        logUniqueERROR("Is not possible to define scaling for materials not in module.");
      }
    }
    return 1.;
  }

  void MaterialObject::Element::build(const std::map<int, int>& newSensorChannels) {
    check();
    // if(destination.state())
    //   std::cout << "DESTINATION " << destination() << " for " << elementName() << std::endl;
    for (const auto& aSensorChannel : newSensorChannels ) {
      sensorChannels_[aSensorChannel.first] = aSensorChannel.second;
    }

    for (const auto& currentSensorNode : referenceSensorNode ) {      
      ReferenceSensor* newReferenceSensor = new ReferenceSensor();
      newReferenceSensor->store(currentSensorNode.second);
      newReferenceSensor->check();
      newReferenceSensor->cleanup();
      referenceSensors_[currentSensorNode.first] = newReferenceSensor;
    }
    /*
    std::cout << "  ELEMENT " << elementName() << std::endl;
    std::cout << "    DATA "
        << " componentName " << (componentName.state() ? componentName() : "NOT_SET")
        << " nSegments " << (nSegments.state() ? std::to_string(nSegments()) : "NOT_SET")
        << " service " << service()
        << " scaleOnSensor " << scaleOnSensor()
        << " quantity " << quantity()
        << " unit " << unit()
        << " station " << (destination.state() ? destination() : "NOT_SET")
        << std::endl;
    */
  }

//  void MaterialObject::Element::chargeTrain(Materialway::Train& train) const {
//    if (service()) {
//      train.addWagon(elementName(), )
//  }

  void MaterialObject::Element::populateMaterialProperties(MaterialProperties& materialProperties) const {
    double quantity;

    if(debugInactivate() == false) {
      if(service() == false) {
        quantity = totalGrams(materialProperties);
        materialProperties.addLocalMass(elementName(), componentName(), quantity);
      }
    }
  }

  void MaterialObject::Element::getLocalElements(ElementsVector& elementsList) const {
    if(service() == false) {
      elementsList.push_back(this);
    }
  }


} /* namespace material */

void MaterialObject::ReferenceSensor::check() {
  PropertyObject::check();
  
  if (!numStripsAcross.state() && !pitchEstimate.state()) throw PathfulException("At least one between numStripsAcross and pitchEstimate must be specified");
  if (numStripsAcross.state() && pitchEstimate.state()) throw PathfulException("Only one between numStripsAcross and pitchEstimate can be specified");
  if (!numSegments.state() && !stripLengthEstimate.state()) throw PathfulException("At least one between numSegments and stripLengthEstimate must be specified");
  if (numSegments.state() && stripLengthEstimate.state()) throw PathfulException("Only one between numSegments and stripLengthEstimate can be specified");
}


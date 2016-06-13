/**
 * @file Station.cpp
 *
 * @date Jul 21, 2014
 * @author Stefano Martina
 */

#include <set>
#include "ConversionStation.h"
//#include "MaterialObject.h"
#include "InactiveElement.h"
#include "MessageLogger.h"

namespace material {

  const std::map<std::string, ConversionStation::Type> ConversionStation::typeString = {
    {"flange", FLANGE},
    {"second", SECOND}
  };

  void ConversionStation::build() {
    try {
      stationType_ = typeString.at(type_());
      buildConversions();
    } catch (const std::out_of_range& ex) {
      logERROR("Station type \"" + type_() + "\" not recognized.");
    }
    cleanup();
  }

  void ConversionStation::routeConvertedElements(MaterialObject& localOutput, MaterialObject& serviceOutput, InactiveElement& inactiveElement) {
    MaterialObject::Element* inputElement;
    //double totalGrams = 0.0;
    double multiplier = 0.0;
    bool converted = false;
    std::set<std::string> warningMaterials;
    
    for (const MaterialObject::Element* currElement : serviceElements_) { //inputElements) {
      converted = false;
      //if the material need to be converted (flange station, or endcap station with right destination)
      if ((stationType_ == FLANGE) || (stationType_ == SECOND && currElement->destination.state() && currElement->destination().compare(stationName_()) == 0)) {
        for (const Conversion* currConversion : conversions) {
          inputElement = currConversion->input->elements[0];
          if (inputElement->elementName().compare(currElement->elementName()) == 0) {
            converted = true;

            multiplier = currElement->quantityInUnit(inputElement->unit(), inactiveElement) / 
              inputElement->quantityInUnit(inputElement->unit(), inactiveElement);
          
            for (const MaterialObject::Element* outputElement : currConversion->outputs->elements) {
              MaterialObject::Element * newElement = new MaterialObject::Element(*outputElement, multiplier);
              if(currElement->debugInactivate()) {  //apply the inactivation also to converteds
                newElement->debugInactivate(true);
              }
              //TODO: check if is ok to do this
              if(currElement->destination.state() && !newElement->destination.state()) {  //apply the same destination of converted element (only if not defined in output conversion rule)
                newElement->destination(currElement->destination());
              }
              if (newElement->service()) {
                if (newElement->unit().compare("g") != 0) {
                  serviceOutput.addElement(newElement);
                } else {
                  logERROR(err_service1 + newElement->elementName() + err_service2);
                }        
              } else {
                localOutput.addElement(newElement);
              }
            }
          }
        }
      }
      if (!converted) {
        serviceOutput.addElement(currElement);
        warningMaterials.insert(currElement->elementName());
      }
    }

    for(auto& warningMaterial : warningMaterials) {
      logWARNING("Element \"" + warningMaterial + "\" ignored by station \"" + stationName_() + "\".");
    }    

    /*
    for (const Conversion* currConversion : conversions) {
      inputElement = currConversion->input->elements[0];
      totalGrams = 0.0;

      for (const MaterialObject::Element* currElement : inputElements) {
        if (inputElement->elementName().compare(currElement->elementName()) == 0) {
          totalGrams += currElement->quantityInGrams(inactiveElement);
        }
      }

      multiplier = totalGrams / inputElement->quantityInGrams(inactiveElement);

      for (const MaterialObject::Element* outputElement : currConversion->outputs->elements) {
    
        MaterialObject::Element * newElement = new MaterialObject::Element(*outputElement, multiplier);

        if (newElement->service()) {
          serviceOutput.addElement(newElement);
        } else {
          localOutput.addElement(newElement);
        }
      }
    }
    */
  }

  /*
  void ConversionStation::routeConvertedServicesTo(MaterialObject& outputObject) const {

  }

  void ConversionStation::routeConvertedLocalsTo(MaterialObject& outputObject) const {

  }
  */


  ConversionStation::Type ConversionStation::stationType() const {
    return stationType_;
  }

  void ConversionStation::buildConversions() {
    //std::cout << "STATION" << std::endl;

    for (auto& currentConversionNode : conversionsNode_) {
      Conversion* newConversion = new Conversion();
      newConversion->store(propertyTree());
      newConversion->store(currentConversionNode.second);
      newConversion->check();
      newConversion->build();

      conversions.push_back(newConversion);
    }
  }

  void ConversionStation::Conversion::build() {
    //std::cout << "  CONVERSION" << std::endl;

    if (inputNode_.size() > 0) {
      input = new Inoutput();
      input->store(propertyTree());
      input->store(inputNode_.begin()->second);
      input->check();
      input->build();

      if(input->elements[0]->unit().compare("g") == 0) {
        logWARNING("Converted element unit is 'g' in conversion rule.");
      }
    }

    if (outputNode_.size() > 0) {
      outputs = new Inoutput();
      outputs->store(propertyTree());
      outputs->store(outputNode_.begin()->second);
      outputs->check();
      outputs->build();
    }
    cleanup();
  }

  void ConversionStation::Inoutput::build() {
    //std::cout << "    INPUT/OUTPUT" << std::endl;

    for  (auto& currentElementNode : elementsNode_) {
      MaterialObject::Element* newElement = new MaterialObject::Element(elementMaterialType);
      newElement->store(propertyTree());
      newElement->store(currentElementNode.second);
      newElement->check();
      //newElement->build();
      newElement->cleanup();

      elements.push_back(newElement);
    }
    cleanup();
  }

  /*
  void ConversionStation::Element::build() {
    std::cout << "      ELEMENT -> "
        << " elementName " << elementName()
        << "; quantity " << quantity()
        << "; unit " << unit()
        << "; service " << (service.state() ? std::to_string(service()) : "NOT_SET" )
        << std::endl;
  }
  */

} /* namespace material */

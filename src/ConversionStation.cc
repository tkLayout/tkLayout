/**
 * @file ConversionStation.cc 
 *
 * @date Jul 21, 2014
 * @author Stefano Martina
 */

#include <set>
#include "ConversionStation.hh"
//#include "MaterialObject.hh"
#include "InactiveElement.hh"
#include "MessageLogger.hh"

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
    std::set<std::string> infoMaterials;
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
	      // newElement is the converted Element.
	      // newElement is assigned the properties from OUTPUT Element.
	      // Notably, in the rare case where outputElement has a componentName (see conversion station cfg file), 
	      // that name is assigned to newElement.
              MaterialObject::Element * newElement = new MaterialObject::Element(*outputElement, multiplier);

	      // Here, new Element is assigned additional info from the INPUT Elements.
	      // newElement is assigned the same componentName as in input
	      // (except if outputElement has a componentName).
	      if (currElement->componentName.state() && !outputElement->componentName.state()) {
		newElement->componentName(currElement->componentName());
	      }
	      // If we end up with a converted element belonging to no component, there is a problem !
	      if (!newElement->componentName.state()) {
		logWARNING("Element " + newElement->elementName() + "which is at output of station" + stationName_() + "has no assigned componentName.");
	      }
              if(currElement->debugInactivate()) {  //apply the inactivation also to converteds
                newElement->debugInactivate(true);
              }
              //TODO: check if is ok to do this
              if(currElement->destination.state() && !newElement->destination.state()) {  //apply the same destination of converted element (only if not defined in output conversion rule)
                newElement->destination(currElement->destination());
              }

	      // We now have our newElement :)
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
        // We may want a warning or an info here
        bool doWarning = false;
        bool doInfo = false;
        if (stationType_ != SECOND) {
          // Always throw a warning for flange stations
          doWarning=true;
        } else {
          // it's a SECOND-level station: let's not overwarn
          // If I was the target of the failed conversion, then it's a warning, otherwise it's an info
          if (currElement->destination.state() && currElement->destination().compare(stationName_()) == 0) doWarning=true;
          else doInfo = true;
        }
        if (doWarning) warningMaterials.insert(currElement->elementName());
        if (doInfo) infoMaterials.insert(currElement->elementName());
      }
    }

    for(auto& warningMaterial : warningMaterials) {
      logWARNING("Element \"" + warningMaterial + "\" ignored by station \"" + stationName_() + "\".");
    }    
    for(auto& infoMaterial : infoMaterials) {
      logINFO("Element \"" + infoMaterial + "\" ignored by station \"" + stationName_() + "\".");
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

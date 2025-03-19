#include <set>
#include "ConversionStation.hh"
#include "InactiveElement.hh"
#include "MessageLogger.hh"

namespace material {

  /////////////////////////////////////////////////////
  /* CONVERSION STATION: 
   * VOLUME AT A SPECIFIC LOCATION, WICH PRODUCES CONVERSIONS.
   * CAN ALSO HAVE CONVERSION-INDEPENDENT MATERIALS.
   *///////////////////////////////////////////////////

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
    double multiplier = 0.0;
    std::set<std::string> infoMaterials;
    std::set<std::string> warningMaterials;
    
    // CONVERSIONS
    // Loop on all incoming materials
    for (const MaterialObject::Element* currElement : serviceElements_) {
      bool converted = false;
      // Find whether the material need to be converted (flange station, or endcap station with right destination)
      if ( (stationType_ == FLANGE) 
	  || (stationType_ == SECOND && currElement->destination.state() && currElement->destination().compare(stationName_()) == 0)
	  ) {

	// Loop on all conversions of the station
        for (const std::unique_ptr<Conversion>& currConversion : conversions_) {

	  // Find the input material of each conversion
	  Inoutput* myInput = currConversion->input();
	  MaterialObject::Element* inputElement = nullptr;
	  if (myInput) {
	    const std::vector<std::unique_ptr<MaterialObject::Element> >& myInputElements = myInput->elements();
	    if (myInputElements.size() == 1) { inputElement = myInputElements.front().get(); }
	  }

	  if (inputElement && currElement) {
	    // If incoming material and conversion input element are the same: do conversion!
	    if (inputElement->elementName().compare(currElement->elementName()) == 0) {
	      converted = true;

	      // Calculate conversion factor
	      multiplier = currElement->quantityInUnit(inputElement->unit(), inactiveElement) / 
		inputElement->quantityInUnit(inputElement->unit(), inactiveElement);

	      // Compute the output and store result
	      Inoutput* myOutput = currConversion->outputs();
	      if (myOutput) {
		for (const std::unique_ptr<MaterialObject::Element>& outputElementUPtr : myOutput->elements()) {
		  MaterialObject::Element* outputElement = outputElementUPtr.get();
		  if (outputElement) {
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

		    // We now have our converted newElement :)
		    // Amount of its material is proportional to the input.
		    // We add it to the relevant volumes (these volumes depend whether the newElement is a service or local).
		    addOutputElements(newElement, localOutput, serviceOutput);
		  }
		}
	      }
	    }
	  }
        }
      }
      // WARNINGS / INFOS on materials that were assigned in cfg to be converted, but that were eventually not converted!!
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


    // MATERIALS INDEPENDENT FROM THE CONVERSIONS
    for (const std::unique_ptr<NonConvertedMaterials>& it : nonConvertedMaterials_) {
      for (const std::unique_ptr<MaterialObject::Element>& nonConvertedOutputElementUPtr : it->elements()) {

	// We now have our non-converted output Element :)
	// Amount of its material is INDEPENDENT from the input.
	// We add it to the relevant volumes (these volumes depend whether the output Element is a service or local).
	MaterialObject::Element* outputElement = nonConvertedOutputElementUPtr.get();
	addOutputElements(outputElement, localOutput, serviceOutput);
      }
    }

  }


  /*
   * Build conversion station.
   * This stores all the station's conversions.
   * It also stores all the materials which are independent from the conversion.
   */
  void ConversionStation::buildConversions() {

    // Stores all conversions.
    for (auto& currentConversionNode : conversionsNode_) {
      std::unique_ptr<Conversion> newConversion (new Conversion(subdetectorName_));
      newConversion->store(propertyTree());
      newConversion->store(currentConversionNode.second);
      newConversion->check();
      newConversion->build();

      conversions_.push_back(std::move(newConversion));
    }

    // Stores conversion-independent info.
    for (const auto& node : nonConvertedMaterialsNode_) {
      std::unique_ptr<NonConvertedMaterials> newNonConverted (new NonConvertedMaterials(subdetectorName_));
      newNonConverted->store(propertyTree());
      newNonConverted->store(node.second);
      newNonConverted->check();
      newNonConverted->build();

      nonConvertedMaterials_.push_back(std::move(newNonConverted));
    }
  
  }


  /*
   * Add the output of the conversion station (whether in converted or non-converted category, local or service) to the relevant volumes.
   * The 'relevant volumes" are:
   * local: MaterialObject& localOutput (assigned to the station volume).
   * service: MaterialObject& serviceOutput (routed from the station volume).
   */
  void ConversionStation::addOutputElements(MaterialObject::Element* outputElement, MaterialObject& localOutput, MaterialObject& serviceOutput) {

    // OUTPUT SERVICES 
    if (outputElement->service()) {
      if (outputElement->unit().compare("g/m") == 0) {
	serviceOutput.addElement(outputElement);
      } 
      else { 
	logERROR(any2str("Services (") 
		 + any2str(outputElement->elementName())
		 + any2str(") in ") 
		 + any2str(outputElement->unit())
		 + any2str(" is not supported. Please use g/m !!")
		 );
      }

      // OUTPUT LOCAL MATERIALS  
    } else {
      localOutput.addElement(outputElement);
    }

  }



  /////////////////////////////////////////////////////
  /* CONVERSIONS: 
   * DEFINE THE ASSOCIATED OUTPUTS TO A GIVEN INPUT
   *///////////////////////////////////////////////////

  ConversionStation::Conversion::Conversion(const std::string subdetectorName) 
    : inputNode_ ("Input", parsedAndChecked()),
      outputNode_ ("Output", parsedAndChecked()),
      subdetectorName_(subdetectorName)
  {};

 
  /*
   * Takes info from cfg files.
   */
  void ConversionStation::Conversion::build() {

    // Takes info from Input node.
    if (inputNode_.size() > 0) {
      input_.reset(new Inoutput(subdetectorName_));
      input_->store(propertyTree());
      input_->store(inputNode_.begin()->second);
      input_->check();
      input_->build();

      const std::vector<std::unique_ptr<MaterialObject::Element> >& myInputElements = input_->elements();
      if (myInputElements.size() != 1) { 
	logERROR("Conversion with " + any2str(myInputElements.size()) + " elements in input, only 1 element is supported."); 
      }
      else {
	if(myInputElements.front()->unit().compare("g") == 0) {
	  logWARNING("Converted element unit is 'g' in conversion rule.");
	}
      }
    }

    // Takes info from Output node.
    if (outputNode_.size() > 0) {
      outputs_.reset(new Inoutput(subdetectorName_));
      outputs_->store(propertyTree());
      outputs_->store(outputNode_.begin()->second);
      outputs_->check();
      outputs_->build();
    }
    cleanup();
  }



  ////////////////////////////////////////////////////////////////////
  /* INOUTPUTS: 
   * DEFINE A LIST OF ELEMENTS MAKING UP THE INPUT OR THE OUTPUT NODE
   *//////////////////////////////////////////////////////////////////

  ConversionStation::Inoutput::Inoutput(const std::string subdetectorName) 
    : elementsNode_ ("Element", parsedOnly()),
      elementMaterialType_(MaterialObject::Type::STATION),
      subdetectorName_(subdetectorName)
  {};


  /*
   * Takes info from cfg files.
   */
  void ConversionStation::Inoutput::build() {

    // Takes info from Element nodes.
    for (auto& currentElementNode : elementsNode_) {
      std::unique_ptr<MaterialObject::Element> newElement (new MaterialObject::Element(elementMaterialType_, subdetectorName_));
      newElement->store(propertyTree());
      newElement->store(currentElementNode.second);
      newElement->check();
      newElement->cleanup();

      elements_.push_back(std::move(newElement));
    }
    cleanup();
  }



  /////////////////////////////////////////////////////
  /* MATERIALS INDEPENDENT FROM THE CONVERSIONS.
   *///////////////////////////////////////////////////

  ConversionStation::NonConvertedMaterials::NonConvertedMaterials(const std::string subdetectorName) 
    : elementsNode_ ("Element", parsedOnly()),
      elementMaterialType_(MaterialObject::Type::STATION),
      subdetectorName_(subdetectorName)
  {};


  /*
   * Takes info from cfg files.
   */
  void ConversionStation::NonConvertedMaterials::build() {

    // Takes info from Element nodes.
    for (const auto& currentElementNode : elementsNode_) {
      std::unique_ptr<MaterialObject::Element> newElement (new MaterialObject::Element(elementMaterialType_, subdetectorName_));
      newElement->store(propertyTree());
      newElement->store(currentElementNode.second);
      newElement->check();     
      newElement->cleanup();

      elements_.push_back(std::move(newElement));
    }
    cleanup();
  }
 

} /* namespace material */

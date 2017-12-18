/**
 * @file MaterialObject.h
 *
 * @date 19/Jun/2014
 * @author Stefano Martina
 */

#ifndef MATERIALOBJECT_H_
#define MATERIALOBJECT_H_

#include "Property.hh"
//#include "Materialway.hh"

class DetectorModule;

namespace insur {
  class MaterialProperties;
}

using insur::MaterialProperties;

namespace material {
  static const std::string err_service1 = "Impossible to use 'g' as unit for service materials, element '";
  static const std::string err_service2 = "' ignored.";

  class MaterialTab;
  class ConversionStation;

  class MaterialObject : public PropertyObject {
  public:
    class Element; //forward declaration for getElementIfService(Element& inputElement)
    class Component;

    typedef std::vector<Component*> ComponentsVector;
    typedef std::vector<const Element*> ElementsVector;
    static const bool ONLY_SERVICES = true;
    static const bool SERVICES_AND_LOCALS = false;

    enum Type {MODULE, ROD, LAYER, SERVICE, STATION};

    MaterialObject(Type materialType);
    MaterialObject(const MaterialObject& other);
    virtual ~MaterialObject();

    double totalGrams(double length, double surface) const;

    virtual void build();
    
    void deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices = false, double gramsMultiplier = 1.) const;
    void addElement(const MaterialObject::Element* element);
    void populateMaterialProperties(MaterialProperties& materialProperties) const;

    ElementsVector& getLocalElements() const;

    bool isPopulated() const;

    //TODO: do methods for interrogate/get materials

  private:
    static const std::map<Type, const std::string> typeString;
    Type materialType_;
    ReadonlyProperty<std::string, NoDefault> type_;
    ReadonlyProperty<std::string, NoDefault> destination_;
    ReadonlyProperty<bool, Default> debugInactivate_;
    PropertyNodeUnique<std::string> materialsNode_; //TODO:see if is better PropertyNode instead
    //    PropertyNode<int> sensorNode_;

    const std::string getTypeString() const;

  public:
    // The sensor channel count
    std::map<int, int> sensorChannels;
    class ReferenceSensor : public PropertyObject {
    public:
      ReadonlyProperty<int, NoDefault> numStripsAcross;
      ReadonlyProperty<double, NoDefault> pitchEstimate;
      ReadonlyProperty<int, NoDefault> numSegments;
      ReadonlyProperty<double, NoDefault> stripLengthEstimate;
      Property<double, NoDefault> width;
      Property<double, NoDefault> length;
    ReferenceSensor() :
      numStripsAcross("numStripsAcross", parsedOnly()),
	pitchEstimate("pitchEstimate", parsedOnly()),
	numSegments("numSegments", parsedOnly()),
	stripLengthEstimate("stripLengthEstimate", parsedOnly()),
	length("length", parsedOnly()), // not checked because a custom checker is defined for RectangularModules
	width("width", parsedOnly())  // same here
	  {}
      int numStripsAcrossEstimate() const {
	if (numStripsAcross.state()) return numStripsAcross();
	else return floor(width() / pitchEstimate() + 0.5);
      }
      int numSegmentsEstimate() const {
	if (numSegments.state()) return numSegments();
	else return floor(length() / stripLengthEstimate() + 0.5);
      }
      int numChannels() const { return numStripsAcrossEstimate() * numSegmentsEstimate(); }

      void check() override;
    };

    // The index of module types, including channel count per sensor
    class MaterialObjectKey {
    public:
      // TODO embed sensorChannel map into a fancier object?...?
      MaterialObjectKey(const std::string& newName, std::map<int, int> newSensorChannels, const std::string& newDestination) : 
        name(newName), 
        sensorChannels(newSensorChannels),
        destination(newDestination)  {}

      bool operator<(const MaterialObjectKey& r ) const {
        if (this->sensorChannels < r.sensorChannels) return true;
        if (this->sensorChannels > r.sensorChannels) return false;
        if (this->name < r.name) return true;
        if (this->name > r.name) return false;
        if (this->destination < r.destination) return true;
        return false;
      }
    private:
      std::string name;
      std::map<int, int> sensorChannels;
      std::string destination;
    };

    class Element : public PropertyObject {
    public:
      enum Unit{GRAMS, MILLIMETERS, GRAMS_METER};
      //static const std::map<Unit, const std::string> unitString;
      static const std::map<std::string, Unit> unitStringMap;
      Property<std::string, NoDefault> componentName; //only the inner component's name
      Property<std::string, NoDefault> elementName;
      Property<bool, Default> service;
      Property<int, Default> scaleOnSensor;
      Property<double, NoDefault> quantity;
      Property<std::string, NoDefault> unit;
      Property<bool, Default> debugInactivate;
      Property<std::string, NoDefault> destination;
      Property<int, Default> targetVolume;
      PropertyNode<int> referenceSensorNode;

      Element(MaterialObject::Type& newMaterialType);
      Element(const Element& original, double multiplier = 1.0);
      //Element(const Element& originElement);
      std::map<int, ReferenceSensor*> referenceSensors_;

      virtual ~Element();
      void build(const std::map<int, int>& newSensorChannels);
      void deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices = false, double gramsMultiplier = 1.) const;
      double quantityInGrams(const DetectorModule& module) const;
      double quantityInGrams(const MaterialProperties& materialProperties) const;
      double quantityInGrams(const double length, const double surface) const;
      double quantityInUnit(const std::string desiredUnit, const MaterialProperties& materialProperties) const;
      double quantityInUnit(const std::string desiredUnit, const double length, const double surface) const;
      double totalGrams(const DetectorModule& module) const;
      double totalGrams(const MaterialProperties& materialProperties) const;
      double totalGrams(const double length, const double surface) const;
      double scalingMultiplier() const;
      void populateMaterialProperties(MaterialProperties& materialProperties) const;
      void getLocalElements(ElementsVector& elementsList) const;
      std::map<int, int> sensorChannels_;

    private:
      const MaterialTab& materialTab_;
      static const std::string msg_no_valid_unit;
      MaterialObject::Type& materialType_;
      //static const std::map<std::string, Materialway::Train::UnitType> unitTypeMap;
    };

    class Component : public PropertyObject {
    public:
      //Property<std::string, NoDefault> componentName;
      PropertyNodeUnique<std::string> componentsNode_;
      PropertyNodeUnique<std::string> elementsNode_;
      Component(MaterialObject::Type& newMaterialType);
      virtual ~Component();
      double totalGrams(double length, double surface) const;
      void build(const std::map<int, int>& newSensorChannels);
      void deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices = false, double gramsMultiplier = 1.) const;
      void populateMaterialProperties(MaterialProperties& materialPropertie) const;
      void getLocalElements(ElementsVector& elementsList) const;

      ComponentsVector components_;
      ElementsVector elements_;

      MaterialObject::Type& materialType_;
    };

    class Materials : public PropertyObject {
    public:
      PropertyNodeUnique<std::string> componentsNode_;
      //Property<double, Computable> radiationLength, interactionLenght;
      Materials(MaterialObject::Type newMaterialType);
      virtual ~Materials();
      double totalGrams(double length, double surface) const;
      void build(const std::map<int, int>& newSensorChannels);
      void deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices = false, double gramsMultiplier = 1.) const;
      void populateMaterialProperties(MaterialProperties& materialProperties) const;
      void getLocalElements(ElementsVector& elementsList) const;

      ComponentsVector components_;

      MaterialObject::Type materialType_;
    };

    //ATTENTION: Materials objects of the same structure are shared between MaterialObject objects
    //   of the modules/layer/etc.. (for containing memory use).
    //   This is not for service routing objects.
    Materials * materials_;

    ElementsVector serviceElements_; //used for MaterialObject not from config file (service routing)
    
  };

  typedef std::vector<const MaterialObject::Element*> ElementsVector;


} /* namespace material */

#endif /* MATERIALOBJECT_H_ */

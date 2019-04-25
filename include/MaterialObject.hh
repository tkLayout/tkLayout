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


#include "MaterialTab.hh"


class DetectorModule;

namespace insur {
  class MaterialProperties;
}

using insur::MaterialProperties;

namespace material {

  class MaterialTab;
  class ConversionStation;

  enum Position { COMMON, DEE, EXTERNAL };

  class MaterialObject : public PropertyObject {
  public:
    class Element; //forward declaration for getElementIfService(Element& inputElement)
    class Component;

    typedef std::vector<Component*> ComponentsVector;
    typedef std::vector<const Element*> ElementsVector;
    static const bool ONLY_SERVICES = true;
    static const bool SERVICES_AND_LOCALS = false;

    enum Type {MODULE, ROD, LAYER, SERVICE, STATION};

    MaterialObject(Type materialType, const std::string subdetectorName);
    MaterialObject(const MaterialObject& other);
    virtual ~MaterialObject();

    const std::string subdetectorName() const { return subdetectorName_; }

    double totalGrams(double length, double surface) const;

    virtual void build();
    
    void deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices = false, 
			  double gramsMultiplier = 1.,
			  Position requestedPosition = COMMON
			  ) const;
    void addElement(const MaterialObject::Element* element);
    void populateMaterialProperties(MaterialProperties& materialProperties) const;

    ElementsVector& getLocalElements() const;

    bool isPopulated() const;

    //TODO: do methods for interrogate/get materials

  private:
    static const std::map<Type, const std::string> typeString;
    std::string subdetectorName_;
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
        width("width", parsedOnly()),  // same here
        length("length", parsedOnly()) // not checked because a custom checker is defined for RectangularModules
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
      MaterialObjectKey(const std::string subdetectorName, const std::string& newName, std::map<int, int> newSensorChannels, const std::string& newDestination) : 
	subdetectorName_(subdetectorName),
        name(newName), 
        sensorChannels(newSensorChannels),
        destination(newDestination)  {}

      bool operator<(const MaterialObjectKey& r ) const {
	if (this->subdetectorName_ < r.subdetectorName_) return true;
	if (this->subdetectorName_ > r.subdetectorName_) return false;
        if (this->sensorChannels < r.sensorChannels) return true;
        if (this->sensorChannels > r.sensorChannels) return false;
        if (this->name < r.name) return true;
        if (this->name > r.name) return false;
        if (this->destination < r.destination) return true;
        return false;
      }
    private:
      std::string subdetectorName_;
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
      Property<Position, Default> position;
      PropertyNode<int> referenceSensorNode;

      Element(MaterialObject::Type& newMaterialType, const std::string subdetectorName);
      Element(const Element& original, double multiplier = 1.0);
      //Element(const Element& originElement);
      std::map<int, ReferenceSensor*> referenceSensors_;

      virtual ~Element();
      void build(const std::map<int, int>& newSensorChannels);
      void deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices = false, double gramsMultiplier = 1., Position requestedPosition = COMMON) const;
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
      const std::string subdetectorName() const { return subdetectorName_; }
      std::map<int, int> sensorChannels_;

    private:
      //const MaterialTab& materialTab_;
      std::string subdetectorName_;
      const MaterialsTable& materialsTable_;
      static const std::string msg_no_valid_unit;
      MaterialObject::Type& materialType_;
      //static const std::map<std::string, Materialway::Train::UnitType> unitTypeMap;
    };

    class Component : public PropertyObject {
    public:
      PropertyNodeUnique<std::string> componentsNode_;
      PropertyNodeUnique<std::string> elementsNode_;
      Component(MaterialObject::Type& newMaterialType, const std::string subdetectorName);
      virtual ~Component();
      double totalGrams(double length, double surface) const;
      void build(const std::map<int, int>& newSensorChannels);
      void deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices = false, double gramsMultiplier = 1., Position requestedPosition = COMMON) const;
      void populateMaterialProperties(MaterialProperties& materialPropertie) const;
      void getLocalElements(ElementsVector& elementsList) const;

      const std::string subdetectorName() const { return subdetectorName_; }

      ComponentsVector components_;
      ElementsVector elements_;

      MaterialObject::Type& materialType_;

    private:
      std::string subdetectorName_;
    };

    class Materials : public PropertyObject {
    public:
      PropertyNodeUnique<std::string> componentsNode_;
      Materials(MaterialObject::Type newMaterialType, const std::string subdetectorName);
      virtual ~Materials();
      double totalGrams(double length, double surface) const;
      void build(const std::map<int, int>& newSensorChannels);
      void deployMaterialTo(MaterialObject& outputObject, const std::vector<std::string>& unitsToDeploy, bool onlyServices = false, double gramsMultiplier = 1., Position requestedPosition = COMMON) const;
      void populateMaterialProperties(MaterialProperties& materialProperties) const;
      void getLocalElements(ElementsVector& elementsList) const;

      const std::string subdetectorName() const { return subdetectorName_; }

      ComponentsVector components_;

      MaterialObject::Type materialType_;

    private:
      std::string subdetectorName_;
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

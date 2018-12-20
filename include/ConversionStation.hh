/**
 * @file ConversionStation.h
 *
 * @date Jul 21, 2014
 * @author Stefano Martina
 */

#ifndef CONVERSIONSTATION_H_
#define CONVERSIONSTATION_H_

#include "Property.hh"
#include "MaterialObject.hh"

namespace insur {
  class InactiveElement;
}

using insur::InactiveElement;

namespace material {
  //class MaterialObject;
  
  class ConversionStation :public MaterialObject {
  public:
    enum Type {ERROR, FLANGE, SECOND};
    
    ConversionStation(const std::string subdetectorName) :
      MaterialObject(MaterialObject::Type::STATION, subdetectorName),
      stationName_ ("stationName", parsedAndChecked()),
      type_ ("type", parsedAndChecked()),
      minZ_ ("minZ", parsedOnly()),
      maxZ_ ("maxZ", parsedOnly()),
      subdetectorName_(subdetectorName),
      stationType_ (ERROR),
      conversionsNode_ ("Conversion", parsedOnly()),
      nonConvertedLocalMaterialsNode_ ("NonConvertedLocalMaterials", parsedOnly())
    {};
    virtual ~ConversionStation() {};

    void build();
    void routeConvertedElements(MaterialObject& localOutput, MaterialObject& serviceOutput, InactiveElement& inactiveElement);
    //void routeConvertedServicesTo(MaterialObject& outputObject) const;
    //void routeConvertedLocalsTo(MaterialObject& outputObject) const;
    Type stationType() const;

    ReadonlyProperty<std::string, NoDefault> stationName_;
    ReadonlyProperty<std::string, NoDefault> type_;
    ReadonlyProperty<double, NoDefault> minZ_;
    ReadonlyProperty<double, NoDefault> maxZ_;
    const double meanZ() const { return (minZ_() + maxZ_()) / 2.; }

  private:
    std::string subdetectorName_;
    static const std::map<std::string, Type> typeString;
    Type stationType_;
    bool valid_;
    PropertyNodeUnique<std::string> conversionsNode_;
    PropertyNodeUnique<std::string> nonConvertedLocalMaterialsNode_;

    void buildConversions();

    /*
    class Element : public PropertyObject {
    public:
      ReadonlyProperty<std::string, NoDefault> elementName;
      ReadonlyProperty<long, NoDefault> quantity;
      ReadonlyProperty<std::string, NoDefault> unit;
      ReadonlyProperty<bool, Default> service;

      Element() :
        elementName ("elementName", parsedAndChecked()),
        quantity ("quantity", parsedAndChecked()),
        unit ("unit", parsedAndChecked()),
        service ("service", parsedOnly(), false) {};
      virtual ~Element() {};

      //void build();
    };
    */


    class NonConvertedLocalMaterials : public PropertyObject {
    public:
      NonConvertedLocalMaterials(const std::string subdetectorName) :
        elementsNode_ ("Element", parsedOnly()),
	elementMaterialType_(MaterialObject::Type::STATION),
	subdetectorName_(subdetectorName)
      {};
      virtual ~NonConvertedLocalMaterials() {};
      void build();
      const std::vector<MaterialObject::Element*> elements() const { return elements_; }

      PropertyNodeUnique<std::string> elementsNode_;

    private:
      std::vector<MaterialObject::Element*> elements_;
      MaterialObject::Type elementMaterialType_;
      std::string subdetectorName_;
    };


    class Inoutput : public PropertyObject {
    public:
      PropertyNodeUnique<std::string> elementsNode_;

      Inoutput(const std::string subdetectorName) :
        elementsNode_ ("Element", parsedOnly()),
        elementMaterialType(MaterialObject::Type::STATION),
	subdetectorName_(subdetectorName)
      {};

        virtual ~Inoutput() {};

      void build();

      std::vector<MaterialObject::Element*> elements;
      MaterialObject::Type elementMaterialType;

    private:
      std::string subdetectorName_;
    };

    class Conversion : public PropertyObject {
    public:
      PropertyNode<std::string> inputNode_;
      PropertyNode<std::string> outputNode_;

      Conversion(const std::string subdetectorName) :
        inputNode_ ("Input", parsedAndChecked()),
        outputNode_ ("Output", parsedAndChecked()),
	subdetectorName_(subdetectorName)
      {};
      virtual ~Conversion() {};

      void build();

      Inoutput* input;
      Inoutput* outputs;

    private:
      std::string subdetectorName_;
    };

    std::vector<Conversion *> conversions;
    std::vector<NonConvertedLocalMaterials*> nonConvertedLocalMaterials_;
  };

} /* namespace material */

#endif /* CONVERSIONSTATION_H_ */

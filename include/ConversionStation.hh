#ifndef CONVERSIONSTATION_H_
#define CONVERSIONSTATION_H_

#include "Property.hh"
#include "MaterialObject.hh"

namespace insur {
  class InactiveElement;
}

using insur::InactiveElement;

namespace material {
  
  /////////////////////////////////////////////////////
  /* CONVERSION STATION: 
   * VOLUME AT A SPECIFIC LOCATION, WICH PRODUCES CONVERSIONS.
   * CAN ALSO HAVE CONVERSION-INDEPENDENT LOCAL MATERIALS.
   *///////////////////////////////////////////////////

  class ConversionStation : public MaterialObject {
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
    Type stationType() const { return stationType_; }

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

    class Conversion;
    class Inoutput;
    class NonConvertedLocalMaterials;
   
    std::vector<Conversion*> conversions;
    std::vector<NonConvertedLocalMaterials*> nonConvertedLocalMaterials_;
  };



  /////////////////////////////////////////////////////
  /* CONVERSIONS: 
   * DEFINE THE ASSOCIATED OUTPUTS TO A GIVEN INPUT
   *///////////////////////////////////////////////////
  class ConversionStation::Conversion : public PropertyObject {
  public:
    PropertyNode<std::string> inputNode_;
    PropertyNode<std::string> outputNode_;

    Conversion(const std::string subdetectorName);
    virtual ~Conversion() {};

    void build();

    Inoutput* input;
    Inoutput* outputs;

  private:
    std::string subdetectorName_;
  };



  ////////////////////////////////////////////////////////////////////
  /* INOUTPUTS: 
   * DEFINE A LIST OF ELEMENTS MAKING UP THE INPUT OR THE OUTPUT NODE
   *//////////////////////////////////////////////////////////////////
  class ConversionStation::Inoutput : public PropertyObject {
  public:
    PropertyNodeUnique<std::string> elementsNode_;

    Inoutput(const std::string subdetectorName);
    virtual ~Inoutput() {};

    void build();

    std::vector<MaterialObject::Element*> elements;
    MaterialObject::Type elementMaterialType;

  private:
    std::string subdetectorName_;
  };



  /////////////////////////////////////////////////////
  /* LOCAL MATERIALS INDEPENDENT FROM THE CONVERSIONS.
   * These are the materials that should not be scaled as a function of the input. They are independent from the conversions!!
   * Whether the input is nothing, 10 g of Cu, or 50 g of ALN, the non converted materials are the same!
   *///////////////////////////////////////////////////
  class ConversionStation::NonConvertedLocalMaterials : public PropertyObject {
  public:
    NonConvertedLocalMaterials(const std::string subdetectorName);
    virtual ~NonConvertedLocalMaterials() {};
    void build();
    const std::vector<std::unique_ptr<MaterialObject::Element> >& elements() const { return elements_; }

  private:
    PropertyNodeUnique<std::string> elementsNode_;
    MaterialObject::Type elementMaterialType_;
    std::string subdetectorName_;
    std::vector<std::unique_ptr<MaterialObject::Element> > elements_;
  };



} /* namespace material */

#endif /* CONVERSIONSTATION_H_ */

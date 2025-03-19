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
   * CAN ALSO HAVE CONVERSION-INDEPENDENT MATERIALS.
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
      nonConvertedMaterialsNode_ ("NonConvertedMaterials", parsedOnly())
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
    PropertyNodeUnique<std::string> nonConvertedMaterialsNode_;

    void buildConversions();
    void addOutputElements(MaterialObject::Element* outputElement, MaterialObject& localOutput, MaterialObject& serviceOutput);

    class Conversion;
    class Inoutput;
    class NonConvertedMaterials;
   
    std::vector<std::unique_ptr<Conversion> > conversions_;
    std::vector<std::unique_ptr<NonConvertedMaterials> > nonConvertedMaterials_;
  };



  /////////////////////////////////////////////////////
  /* CONVERSIONS: 
   * DEFINE THE ASSOCIATED OUTPUTS TO A GIVEN INPUT
   *///////////////////////////////////////////////////
  class ConversionStation::Conversion : public PropertyObject {
  public:
    Conversion(const std::string subdetectorName);
    virtual ~Conversion() {};

    void build();
    Inoutput* input() { return input_.get(); }
    Inoutput* outputs() { return outputs_.get(); }

  private:
    PropertyNode<std::string> inputNode_;
    PropertyNode<std::string> outputNode_;
    std::string subdetectorName_;
    std::unique_ptr<Inoutput> input_;
    std::unique_ptr<Inoutput> outputs_;
  };



  ////////////////////////////////////////////////////////////////////
  /* INOUTPUTS: 
   * DEFINE A LIST OF ELEMENTS MAKING UP THE INPUT OR THE OUTPUT NODE
   *//////////////////////////////////////////////////////////////////
  class ConversionStation::Inoutput : public PropertyObject {
  public:
    Inoutput(const std::string subdetectorName);
    virtual ~Inoutput() {};

    void build();
    const std::vector<std::unique_ptr<MaterialObject::Element> >& elements() const { return elements_; }

  private:
    PropertyNodeUnique<std::string> elementsNode_;
    MaterialObject::Type elementMaterialType_;
    std::string subdetectorName_;
    std::vector<std::unique_ptr<MaterialObject::Element> > elements_;
  };



  /////////////////////////////////////////////////////
  /* MATERIALS INDEPENDENT FROM THE CONVERSIONS.
   * These are the materials that should not be scaled as a function of the input. They are independent from the conversions!!
   * The non converted materials are of the same quantity, whether the input is nothing, 10 g of Cu, or 50 g of ALN!
   * These materials can be local (assigned to the station volume) or services (routed from the station volume).
   *///////////////////////////////////////////////////
  class ConversionStation::NonConvertedMaterials : public PropertyObject {
  public:
    NonConvertedMaterials(const std::string subdetectorName);
    virtual ~NonConvertedMaterials() {};

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

/**
 * @file SupportStructure.h
 *
 * @date 25/Nov/2014
 * @author Stefano Martina
 */

#ifndef SUPPORTSTRUCTURE_H_
#define SUPPORTSTRUCTURE_H_

#include <vector>
#include <string>
#include "Property.h"
#include "global_constants.h"

namespace insur {
  class MaterialProperties;
  class InactiveElement;
  class InactiveSurfaces;
}
using insur::MaterialProperties;
using insur::InactiveElement;
using insur::InactiveSurfaces;

class Barrel;

namespace material {
  class MaterialTab;

  class SupportStructure : public PropertyObject {
  private:
    class Component;
    class Element;

  public:
    enum Type {CUSTOM, AUTO, TOP, BOTTOM};
    enum Direction {HORIZONTAL, VERTICAL};

    static const std::map<std::string, Type> typeStringMap;
    static const std::map<std::string, Direction> directionStringMap;

    typedef std::vector<const Component*> ComponentsVector;
    typedef std::vector<const Element*> ElementsVector;

    Property<std::string, NoDefault> type;
    Property<double, NoDefault> autoPosition;
    Property<double, NoDefault> customZMin;
    Property<double, NoDefault> customRMin;
    Property<double, NoDefault> customLength;
    Property<std::string, NoDefault> customDir;
    
    SupportStructure();
    virtual ~SupportStructure() {};
    void buildInTracker();
    void buildInBarrel(Barrel& barrel);

    void updateInactiveSurfaces(InactiveSurfaces& inactiveSurfaces);
    
  private:
    const double inactiveElementWidth = insur::volume_width;
    const double autoLayerMarginUpper = 1.; //margins for the auto barrel support
    const double autoLayerMarginLower = 2. + insur::volume_width; //upper is upper for the support (is lower for the layer) and viceversa

    PropertyNodeUnique<std::string> componentsNode;

    ComponentsVector components_;
    std::vector<InactiveElement*> inactiveElements;
    Type supportType_;
    Direction direction_;

    void buildBase();
    void populateMaterialProperties(MaterialProperties& materialPropertie) const;
    void buildInactiveElementPair(Direction direction, double zStart, double rStart, double length);


    class Component : public PropertyObject {
    public:
      PropertyNodeUnique<std::string> componentsNode;
      PropertyNodeUnique<std::string> elementsNode;
      Component();
      virtual ~Component() {};

      void build();
      void populateMaterialProperties(MaterialProperties& materialPropertie) const;

      ComponentsVector components_;
      ElementsVector elements_;
    };
    
    class Element : public PropertyObject {
    public:
      enum Unit{GRAMS, MILLIMETERS, GRAMS_METER};
      static const std::map<std::string, Unit> unitStringMap;
      Property<std::string, NoDefault> componentName; //only the inner component's name
      Property<std::string, NoDefault> elementName;
      Property<double, NoDefault> quantity;
      Property<std::string, NoDefault> unit;
      Property<bool, Default> debugInactivate;

      Element();
      virtual ~Element() {};

      double quantityInGrams(double length, double surface) const;
      void populateMaterialProperties(MaterialProperties& materialPropertie) const;
    private:
      const MaterialTab& materialTab_;
      static const std::string msg_no_valid_unit;
    };
  };
}

#endif // SUPPORTSTRUCTURE_H_

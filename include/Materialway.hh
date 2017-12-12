/**
 * @file Materialway.h
 *
 * @date 31/mar/2014
 * @author Stefano Martina
 */

#ifndef MATERIALWAY_H_
#define MATERIALWAY_H_

#include <map>
#include <vector>
#include <utility>
#include <set>
#include <string>
#include "MaterialObject.hh"
//#include "global_constants.hh"

class DetectorModule;
class Tracker;
class Barrel;
class Endcap;
class Visitable;
class Disk;
class Layer;

namespace insur {
  class InactiveSurfaces;
  class InactiveTube;
  class InactiveRing;
  class InactiveElement;
  class MatCalc;
}

using insur::InactiveSurfaces;
using insur::InactiveTube;
using insur::InactiveRing;
using insur::InactiveElement;
using insur::MatCalc;


namespace material {
  static const std::string inactiveElementError = "Section without inactiveElement.";

  class MaterialObject;
  class ConversionStation;
  class WeightDistributionGrid;

  /**
   * @class Materialway
   * @brief Represents a track where the materials are
   * routed from the modules to the appropriate end
   *
   * The materialway is make up from single elements (Element),
   * every element point to the possible single next element.
   * Every module is linked to a materialway element where
   * it puts its materials to be routed, the material is then deposited
   * on the element and forwarded to the next element, passing through the
   * chain from an element to the subsequent one until the material reach
   * its designated end.
   * Every materialway element can perform some transformation on the routed
   * material before forward it to the next element.
   */
  class Materialway {
  private:
    enum Direction { HORIZONTAL, VERTICAL };

    class Section;
    class Station;      //because Train need to use Station and Section, and Station and Section need to use Train

    class RodSectionsStation {
    private:
      std::vector<Section*> sections_;
      Section* station_;
    public:
      RodSectionsStation();
      ~RodSectionsStation();
      void addSection(Section* section);
      void setStation(Section* station);
      std::vector<Section*>& getSections();
      Section* getStation();
    };

    //typedef std::map<const Layer*, std::pair<std::vector<Section*>, Section*> > LayerRodSectionsMap;
    //typedef std::map<const Disk*, std::pair<std::vector<Section*>, Section*> > DiskRodSectionsMap;

    typedef std::map<const Layer*, RodSectionsStation> LayerRodSectionsMap;
    typedef std::map<const Disk*, RodSectionsStation> DiskRodSectionsMap;

  private:

    /**
     * @class Section
     * @brief Represents a single element of the materialway
     */
    class Section {
    public:
      Section(int minZ, int minR, int maxZ, int maxR, Direction bearing, Section* nextSection, bool debug);
      Section(int minZ, int minR, int maxZ, int maxR, Direction bearing, Section* nextSection);
      Section(int minZ, int minR, int maxZ, int maxR, Direction bearing);
      Section(const Section& other);
      virtual ~Section();

      int isHit(int z, int r, int end, Direction aDirection) const;
      void minZ(int minZ);
      void minR(int minR);
      void maxZ(int maxZ);
      void maxR(int maxR);
      void bearing(Direction bearing);
      void nextSection(Section* nextSection);
      int minZ() const;
      int minR() const;
      int maxZ() const;
      int maxR() const;
      int lenght() const;
      Direction bearing() const;
      Section* nextSection() const;
      bool hasNextSection() const;
      MaterialObject& materialObject();
      void inactiveElement(InactiveElement* inactiveElement);
      InactiveElement* inactiveElement() const;
      //Section* appendNewSection

      virtual void getServicesAndPass(const MaterialObject& source);
      virtual void getServicesAndPass(const MaterialObject& source, const std::vector<std::string>& unitsToPass);

      bool debug_;
    private:
      int minZ_, minR_, maxZ_, maxR_;
      Section* nextSection_;
      Direction bearing_;
      MaterialObject materialObject_;
      InactiveElement* inactiveElement_; /**< The InactiveElement for hooking up to the existing infrastructure */
    protected:
      const std::vector<std::string> unitsToPass_ = {"g/m", "mm"};
    }; //class Section

    class Station : public Section {
    public:
      Station(int minZ, int minR, int maxZ, int maxR, Direction bearing, ConversionStation& conversionStation, Section* nextSection);
      Station(int minZ, int minR, int maxZ, int maxR, Direction bearing, ConversionStation& conversionStation);
      virtual ~Station();

      // enum Type { LAYER, TERMINUS };
      /*struct ConversionRule {

        };*/

      virtual void getServicesAndPass(const MaterialObject& source);
      virtual void getServicesAndPass(const MaterialObject& source, const std::vector<std::string>& unitsToPass);

      ConversionStation& conversionStation();
      MaterialObject& outgoingMaterialObject();

    private:
      size_t labelHash;
      ConversionStation& conversionStation_;
      MaterialObject outgoingMaterialObject_;
    };

    typedef std::vector<Section*> SectionVector;
    typedef std::vector<Station*> StationVector;
    typedef std::map<const DetectorModule*, Section*> ModuleSectionMap;

    /**
     * @class Boundary
     * @brief Represents a boundary where the services are routed around
     */
    class Boundary {
    public:
      Boundary();
      //Boundary(const Visitable* containedElement, int minZ, int minR, int maxZ, int maxR);
      Boundary(const Visitable* containedElement, const bool isBarrel, int minZ, int minR, int maxZ, int maxR);
      virtual ~Boundary();

      int isHit(int z, int r, Direction aDirection) const;

      void minZ(int minZ);
      void minR(int minR);
      void maxZ(int maxZ);
      void maxR(int maxR);
      void outgoingSection(Section* outgoingSection);
      int minZ() const;
      int minR() const;
      int maxZ() const;
      int maxR() const;
      Section* outgoingSection();
      const bool isBarrel() const { return isBarrel_; }

      bool operator==(Boundary* const& b) {
        return (minZ_ == b->minZ() && minR_ == b->minR() && maxZ_ == b->maxZ() && maxR_ == b->maxR());
      }

    private:
      int minZ_, minR_, maxZ_, maxR_;
      Section* outgoingSection_;
      const Visitable* containedElement_;
      const bool isBarrel_;
    }; //class Boundary

    typedef std::map<const Barrel*, Boundary*> BarrelBoundaryMap;
    typedef std::map<const Endcap*, Boundary*> EndcapBoundaryMap;

    struct BoundaryComparator {
      bool operator()(Boundary* const& one, Boundary* const& two) {
        //return ((one->maxZ() > two->maxZ()) || (one->maxR() > two->maxR()));
        return ((one->maxZ() + one->maxR()) > (two->maxZ() + two->maxR()));
      }
    };

    typedef std::set<Boundary*, BoundaryComparator> BoundariesSet;

    /**
     * @class OuterUsher
     * @brief Is the core of the functionality that builds sections across boundaries
     * starting from a point and ending to another section or to the upper right angle
     */
    class OuterUsher {
     public:
      OuterUsher(SectionVector& sectionsList, BoundariesSet& boundariesList);
      virtual ~OuterUsher();

      void go(Boundary* boundary, const Tracker& tracker);         /**< start the process of section building, returns pointer to the first */
    private:
      SectionVector& sectionsList_;
      BoundariesSet& boundariesList_;

      Direction buildDirection(const int& startZ, const int& startR, const bool& hasStepInEndcapsOuterRadius, const int& numBarrels);
      bool findBoundaryCollision(int& collision, int& border, int startZ, int startR, const Tracker& tracker, Direction direction);
      bool findSectionCollision(std::pair<int,Section*>& sectionCollision, int startZ, int startR, int end, Direction direction);
      bool buildSection(Section*& firstSection, Section*& lastSection, int& startZ, int& startR, int end, Direction direction);
      bool buildSectionPair(Section*& firstSection, Section*& lastSection, int& startZ, int& startR, int collision, int border, Direction direction);
      Section* splitSection(Section* section, int collision, Direction direction);
      Direction inverseDirection(Direction direction) const;
      void updateLastSectionPointer(Section* lastSection, Section* newSection);
    }; //class OuterUsher

    /**
     * @class InnerUsher
     * @brief Build the internal sections of a boundary
     */
    class InnerUsher {
    public:
      InnerUsher(
          SectionVector& sectionsList,
          StationVector& stationListFirst,
          StationVector& stationListSecond,
          BarrelBoundaryMap& barrelBoundaryAssociations,
          EndcapBoundaryMap& endcapBoundaryAssociations,
          ModuleSectionMap& moduleSectionAssociations,
          LayerRodSectionsMap& layerRodSections,
          DiskRodSectionsMap& diskRodSections);
      virtual ~InnerUsher();

      void go(Tracker& tracker);
    private:
      SectionVector& sectionsList_;
      StationVector& stationListFirst_;
      StationVector& stationListSecond_;
      BarrelBoundaryMap& barrelBoundaryAssociations_;
      EndcapBoundaryMap& endcapBoundaryAssociations_;
      ModuleSectionMap& moduleSectionAssociations_;
      LayerRodSectionsMap& layerRodSections_;
      DiskRodSectionsMap& diskRodSections_;
    }; //class InnerUsher

    /**
     * @class ModuleUsher
     * @brief Add materials to the modules
     */
    class ModuleUsher {
    public:
      ModuleUsher();
      void go(Tracker& tracker);
    }; //class ModuleUsher


  public:
    Materialway();
    virtual ~Materialway();

    bool build(Tracker& tracker, InactiveSurfaces& inactiveSurface, WeightDistributionGrid& weightDistribution);

    static const double gridFactor;                                     /**< the conversion factor for using integers in the algorithm (helps finding collisions),
                                                                            actually transforms millimiters in microns */
    static const int sectionWidth;     /**< the width of a section */
    static const int safetySpace;           /**< the safety space between sections */
    //static const double globalMaxZ_mm;                     /**< the Z coordinate of the end point of the sections */
    //static const double globalMaxR_mm;                   /**< the rho coordinate of the end point of the sections */
    //static const int globalMaxZ;
    //static const int globalMaxR;
    static const int boundaryPaddingBarrel;
    static const int boundaryPaddingEndcaps;
    static const int boundaryPrincipalPaddingBarrel;
    static const int boundaryPrincipalPaddingEndcaps;
    static const int globalMaxZPadding;
    static const int globalMaxRPadding;
    static const int layerSectionMargin;
    static const int diskSectionMargin;
    static const int layerSectionRightMargin;
    static const int diskSectionUpMargin;
    static const int sectionTolerance;
    static const int layerStationLenght;
    static const int layerStationWidth;
    static const double radialDistribError;

    static int discretize(double input);
    static double undiscretize(int input);

  private:
    BoundariesSet boundariesList_;       /**< Vector for storing all the boundaries */
    SectionVector sectionsList_;         /**< Vector for storing all the sections (also stations)*/
    StationVector stationListFirst_;         /**< Pointers to first step stations*/
    StationVector stationListSecond_;         /**< Pointers to second step stations*/

    OuterUsher outerUsher;
    InnerUsher innerUsher;

    bool buildBoundaries(const Tracker& tracker);             /**< build the boundaries around barrels and endcaps */
    void buildExternalSections(const Tracker& tracker);       /**< build the sections outside the boundaries */
    void buildInternalSections(Tracker& tracker);                             /**< build the sections inside the boundaries */
    void buildInactiveElements();
    void routeServices(const Tracker& tracker);
    //void routeModuleServices();
    //void routeRodMaterials();
    void firstStepConversions();
    void secondStepConversions();
    void createModuleCaps(Tracker& tracker);
    void duplicateSections();
    void populateAllMaterialProperties(Tracker& tracker, WeightDistributionGrid& weightDistribution);
    //void calculateMaterialValues(Tracker& tracker);
    void buildInactiveSurface(Tracker& tracker, InactiveSurfaces& inactiveSurface);
    void calculateMaterialValues(InactiveSurfaces& inactiveSurface, Tracker& tracker);
    //InactiveElement* buildOppositeInactiveElement(InactiveElement* inactiveElement);


    BarrelBoundaryMap barrelBoundaryAssociations_;
    EndcapBoundaryMap endcapBoundaryAssociations_;
    ModuleSectionMap moduleSectionAssociations_; /**< Map that associate each module with the section that it feeds */
    LayerRodSectionsMap layerRodSections_;      /**< maps for sections of the rods */
    DiskRodSectionsMap diskRodSections_;
    //std::map<Boundary&, Section*> boundarySectionAssociations;         /**< Map that associate each boundary with the outgoing section (for the construction) */
  };

} /* namespace material */

#endif /* MATERIALWAY_H_ */

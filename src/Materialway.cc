/**
 * @file Materialway.cc
 *
 * @date 31/mar/2014
 * @author Stefano Martina
 */

#include "Materialway.hh"
#include "DetectorModule.hh"
#include "Tracker.hh"
#include "Visitable.hh"
#include "InactiveSurfaces.hh"
#include "InactiveTube.hh"
#include "InactiveRing.hh"
#include "InactiveElement.hh"
#include "MatCalc.hh"
#include "MaterialObject.hh"
#include "ConversionStation.hh"
#include "Barrel.hh"
#include "Endcap.hh"
#include "Disk.hh"
#include "Layer.hh"
#include "WeightDistributionGrid.hh"
#include "StopWatch.hh"

#include <ctime>


namespace material {

  Materialway::RodSectionsStation::RodSectionsStation() {}

  Materialway::RodSectionsStation::~RodSectionsStation() {}

  void Materialway::RodSectionsStation::addSection(Materialway::Section* section) {
    sections_.push_back(section);
  }

  void Materialway::RodSectionsStation::setStation(Materialway::Section* station) {
    station_ = station;
  }

  std::vector<Materialway::Section*>& Materialway::RodSectionsStation::getSections() {
    return sections_;
  }

  Materialway::Section* Materialway::RodSectionsStation::getStation() {
    return station_;
  }


  //=================================================================================
  //START Materialway::Boundary
  Materialway::Boundary::Boundary(const Visitable* containedElement, const bool isBarrel, int minZ, int minR, int maxZ, int maxR, bool hasCShapeOutgoingServices, int cShapeMinZ) :
      containedElement_(containedElement),
      isBarrel_(isBarrel),
      minZ_(minZ),
      maxZ_(maxZ),
      minR_(minR),
      maxR_(maxR),
      hasCShapeOutgoingServices_(hasCShapeOutgoingServices),
      cShapeMinZ_(cShapeMinZ),
      outgoingSection_(nullptr) {}

  Materialway::Boundary::Boundary() :
    Boundary(nullptr, true, 0, 0, 0, 0, false, std::numeric_limits<int>::max()){}
  Materialway::Boundary::~Boundary() {}


  int Materialway::Boundary::isHit(int z, int r, Direction aDirection) const {
    if (aDirection==HORIZONTAL) {
      if ((minR_<r)&&(maxR_>r)) {
        if (minZ_>z) return minZ_;
        else if (maxZ_>z) return -1;
      }
    } else {
      if ((minZ_<z)&&(maxZ_>z)) {
        if (minR_>r) return minR_;
        else if (maxR_>r) return -1;
      }
    }
    return 0;
  }

  void Materialway::Boundary::minZ(int minZ) {
    minZ_ = minZ;
  }
  void Materialway::Boundary::minR(int minR) {
    minR_ = minR;
  }
  void Materialway::Boundary::maxZ(int maxZ) {
    maxZ_ = maxZ;
  }
  void Materialway::Boundary::maxR(int maxR) {
    maxR_ = maxR;
  }
  void Materialway::Boundary::outgoingSection(Section* outgoingSection) {
    outgoingSection_ = outgoingSection;
  }
  int Materialway::Boundary::minZ() const {
    return minZ_;
  }
  int Materialway::Boundary::minR() const {
    return minR_;
  }
  int Materialway::Boundary::maxZ() const {
    return maxZ_;
  }
  int Materialway::Boundary::maxR() const {
    return maxR_;
  }
  Materialway::Section* Materialway::Boundary::outgoingSection() {
    return outgoingSection_;
  }
  //END Materialway::Boundary
  //=================================================================================
  //START Materialway::Section
  Materialway::Section::Section(int minZ, int minR, int maxZ, int maxR, Direction bearing, Section* nextSection, bool debug) :
    debug_(debug),
    minZ_(minZ),
    minR_(minR),
    maxZ_(maxZ),
    maxR_(maxR),
    myLocation_(Location::EXTERNAL),
    bearing_(bearing),
    nextSection_(nextSection),
    inactiveElement_(nullptr), 
    materialObject_ (MaterialObject::SERVICE, "") 
  {}

  Materialway::Section::Section(int minZ, int minR, int maxZ, int maxR, Direction bearing, Section* nextSection) :
    Section(minZ,
            minR,
            maxZ,
            maxR,
            bearing,
            nextSection,
            false) 
  {}

  Materialway::Section::Section(int minZ, int minR, int maxZ, int maxR, Direction bearing, Section* nextSection, const Location& location) :
    Section(minZ,
            minR,
            maxZ,
            maxR,
            bearing,
            nextSection,
            false) 
  {
    myLocation_ = location;
  }

  Materialway::Section::Section(int minZ, int minR, int maxZ, int maxR, Direction bearing) :
    Section(minZ,
            minR,
            maxZ,
            maxR,
            bearing,
            nullptr) 
  {}

  Materialway::Section::Section(const Section& other) :
    debug_(other.debug_),
    minZ_(other.minZ_),
    minR_(other.minR_),
    maxZ_(other.maxZ_),
    maxR_(other.maxR_),
    bearing_(other.bearing_),
    nextSection_(nullptr),
    inactiveElement_(nullptr), //attention, nextSection and inactiveElement are not copied
    materialObject_(other.materialObject_) {

    double zLength = undiscretize(other.maxZ() - other.minZ());
    double zOffset = undiscretize(other.minZ());
    double innerRadius = undiscretize(other.minR());
    double rWidth = undiscretize(other.maxR() - other.minR());

    if (other.bearing() == HORIZONTAL) {
      InactiveTube* tube = new InactiveTube;
      tube->setZLength(zLength);
      tube->setZOffset(zOffset);
      tube->setInnerRadius(innerRadius);
      tube->setRWidth(rWidth);
      tube->setFinal(true);
      inactiveElement(tube);
    } else {
      InactiveRing* ring = new InactiveRing;
      ring->setZLength(zLength);
      ring->setZOffset(zOffset);
      ring->setInnerRadius(innerRadius);
      ring->setRWidth(rWidth);
      ring->setFinal(true);
      inactiveElement(ring);
    }
  }


  int Materialway::Section::isHit(int z, int r, int end, Direction aDirection) const {
    if (aDirection==HORIZONTAL) {
      if ((minR() - sectionWidth - safetySpace < r)&&(maxR() + sectionWidth + safetySpace > r)) {
        if (minZ_>z){
          if (minZ() <= end + safetySpace) {
            return minZ();
          }
        } else if (maxZ()>z) {
          return -1;
        }
      }
    } else {
      if ((minZ() - sectionWidth - safetySpace < z)&&(maxZ() + sectionWidth + safetySpace > z)) {
        if (minR()>r) {
          if (minR() <= end + safetySpace) {
            return minR();
          }
        } else if (maxR()>r) {
          return -1;
        }
      }
    }
    return 0;
  }

  void Materialway::Section::minZ(int minZ){
    minZ_ = minZ;
  }
  void Materialway::Section::minR(int minR){
    minR_ = minR;
  }
  void Materialway::Section::maxZ(int maxZ){
    maxZ_ = maxZ;
  }
  void Materialway::Section::maxR(int maxR){
    maxR_ = maxR;
  }
  void Materialway::Section::bearing(Direction bearing) {
    bearing_ = bearing;
  }
  void Materialway::Section::nextSection(Section* nextSection) {
    nextSection_ = nextSection;
  }

  int Materialway::Section::minZ() const {
    return minZ_;
  }
  int Materialway::Section::minR() const {
    return minR_;
  }
  int Materialway::Section::maxZ() const {
    return maxZ_;
  }
  int Materialway::Section::maxR() const {
    return maxR_;
  }
  int Materialway::Section::lenght() const {
    if (bearing() == HORIZONTAL) {
      return maxZ() - minZ();
    } else {
      return maxR() - minR();
    }
  }
  // NB: Very bad idea to have lengths in integers (even if all lengths are multiplied by 1000, ie precision is 1 um).
  // Do what we can with present codebase...
  const double Materialway::Section::getVolume() const {
    const double crossSection =  M_PI * (pow((double)maxR(), 2.) - pow((double)minR(), 2.));
    const double zLength = (double)maxZ() - (double)minZ();
    const double volume = crossSection * zLength;
    return volume;
  }
  Materialway::Direction Materialway::Section::bearing() const {
    return bearing_;
  }
  Materialway::Section* Materialway::Section::nextSection() const {
    return nextSection_;
  }
  bool Materialway::Section::hasNextSection() const {
    return (nextSection_ == nullptr)? false : true;
  }

  MaterialObject& Materialway::Section::materialObject() {
    return materialObject_;
  }

  void Materialway::Section::inactiveElement(InactiveElement* inactiveElement) {
    inactiveElement_ = inactiveElement;
  }
  InactiveElement* Materialway::Section::inactiveElement() const {
    return inactiveElement_;
  }
  void Materialway::Section::getServicesAndPass(const MaterialObject& source) {
    source.deployMaterialTo(materialObject(), unitsToPass_, MaterialObject::ONLY_SERVICES);
    if(hasNextSection()) {
      nextSection()->getServicesAndPass(source);
    }
  }

  void Materialway::Section::getServicesAndPass(const MaterialObject& source, const std::vector<std::string>& unitsToPass) {
    source.deployMaterialTo(materialObject(), unitsToPass, MaterialObject::ONLY_SERVICES);
    if(hasNextSection()) {
      nextSection()->getServicesAndPass(source);
    }
  }

  //END Materialway::Section
  //=================================================================================
  //START Materialway::Station

  Materialway::Station::Station(int minZ, int minR, int maxZ, int maxR, Direction bearing, ConversionStation& conversionStation, Section* nextSection) :
      Section(minZ, minR, maxZ, maxR, bearing, nextSection),
      conversionStation_ (conversionStation),
      outgoingMaterialObject_ (MaterialObject::SERVICE, "") {}

  Materialway::Station::Station(int minZ, int minR, int maxZ, int maxR, Direction bearing, ConversionStation& conversionStation) :
      Station(minZ, minR, maxZ, maxR, bearing, conversionStation, nullptr) {}
  Materialway::Station::~Station() {}

  void Materialway::Station::getServicesAndPass(const MaterialObject& source) {
    source.deployMaterialTo(conversionStation_, unitsToPass_, MaterialObject::ONLY_SERVICES);
    //don't pass
  }

  void Materialway::Station::getServicesAndPass(const MaterialObject& source, const std::vector<std::string>& unitsToPass) {
    source.deployMaterialTo(conversionStation_, unitsToPass, MaterialObject::ONLY_SERVICES);
    //don't pass
  }

  ConversionStation& Materialway::Station::conversionStation() {
    return conversionStation_;
  }

  MaterialObject& Materialway::Station::outgoingMaterialObject() {
    return outgoingMaterialObject_;
  }


  //END Materialway::Station
  //=================================================================================
  //START Materialway::OuterUsher
  Materialway::OuterUsher::OuterUsher(SectionVector& sectionsList, BoundariesSet& boundariesList) :
    sectionsList_(sectionsList),
    boundariesList_(boundariesList) {}
  Materialway::OuterUsher::~OuterUsher() {}

  void Materialway::OuterUsher::go(Boundary* boundary, const Tracker& tracker) {
    int startZ, startR, collision, border;
    Direction direction;
    bool wasRoutedAlongCShape = false;
    bool foundBoundaryCollision, noSectionCollision;
    Section* lastSection = nullptr;
    Section* firstSection = nullptr;

    startZ = boundary->maxZ();
    startR = boundary->maxR();

    bool going=true;
    while (going) {
      // compute direction along which external sections will be built
      direction = buildDirection(startZ, startR, tracker.hasStepInEndcapsOuterRadius(), tracker.barrels().size());

      // If requested, build the 'C-shape' outgoing services.
      // Example: this is used on the outgoing services over TBPX.
      // The services leave from barrel's (maxZ, MaxR). They go horizontally towards Z = cShapeMinZ, then back horizontally towards bigger |Z|.
      if (boundary->hasCShapeOutgoingServices() && !wasRoutedAlongCShape) {
	routeOutgoingServicesAlongCShape(firstSection, lastSection, startR, direction, startZ, boundary->cShapeMinZ());	
	wasRoutedAlongCShape = true;
      }
 
      // Look for a possible collision between the services to be routed and all boundaries
      foundBoundaryCollision = findBoundaryCollision(collision, border, startZ, startR, tracker, direction);

      // Lastly, Build the corresponding "elbow shaped" services (pair of 2 services)
      noSectionCollision = buildSectionPair(firstSection, lastSection, startZ, startR, collision, border, direction);

      going = foundBoundaryCollision && noSectionCollision;
    }

    boundary->outgoingSection(firstSection);
  }


 /**
   * This computes the direction along which the external sections will be built.
   * @param startZ : Z coordinate of the starting point
   * @param hasStepInEndcapsOuterRadius : Has the considered tracker a step in its outer radius, in the endcaps region ?
   * @return direction : The direction along which the external sections will be built
   */
  Materialway::Direction Materialway::OuterUsher::buildDirection(const int& startZ, const int& startR, const bool& hasStepInEndcapsOuterRadius, const int& numBarrels) {
    // By default, direction is vertical
    Direction direction = VERTICAL;
    if (boundariesList_.size() > 0) {

      // Special case where there is a step in the outer radius in the endcaps
      if (hasStepInEndcapsOuterRadius && (startZ <= (*boundariesList_.cbegin())->minZ())) {

	// set direction to horizontal for the barrel and all the first endcaps
	direction = HORIZONTAL;

	// In case where there are several barrels, direction horizontal should only be set for the top of the outermost barrel.
	if (numBarrels > 1) {
	  // Get the outermost Barrel boundary (The Barrel which has the biggest value = (boundary->maxR() + boundary->maxZ()))
	  auto biggestBarrelBoundary = std::find_if(boundariesList_.begin(), boundariesList_.end(), 
						    [&](const Boundary* b) { return b->isBarrel(); } );
	  if (startR < (*biggestBarrelBoundary)->maxR()) direction = VERTICAL;
	}    	
      }
    }
    else logERROR("Try to build direction for routing services along boundaries, but no boundary !");
    return direction;
  }


  /**
   * This is used to build the 'C-shape' outgoing services.
   * Example: this is used on the outgoing services out of TBPX.
   * The services leave from barrel's (maxZ, MaxR). They go horizontally towards Z = cShapeMinZ, then back horizontally towards bigger |Z|.
   */
  void Materialway::OuterUsher::routeOutgoingServicesAlongCShape(Section*& firstSection, Section*& lastSection, int& startR, const Direction direction, int startZ, int cShapeMinZ) {
    if (direction == VERTICAL) {
      logERROR("Route outgoing services along C-shape: supported on horizontal direction only. Otherwise, it's a U-shape!");
    }

    // HORIZONTAL DIRECTION ONLY
    else {
      // FIRST, GO TOWARDS SMALLER |Z|.
      // Z: startZ -> cShapeMinZ
      const bool isTowardsBiggerZ = false;
      int towardsSmallerZCablingMinR = startR + safetySpace;
      bool noCollision = buildSection(firstSection, lastSection, startZ, towardsSmallerZCablingMinR, cShapeMinZ, direction, isTowardsBiggerZ);
      if (!noCollision) logERROR("Found a collision, which is not supported.");

      // THEN, GO BACK TOWARDS BIGGER |Z|.
      // Z: cShapeMinZ -> startZ
      int towardsBiggerZCablingMinR = towardsSmallerZCablingMinR + sectionWidth + layerStationLenght + safetySpace;
      noCollision = buildSection(firstSection, lastSection, cShapeMinZ, towardsBiggerZCablingMinR, startZ, direction, !isTowardsBiggerZ);
      if (!noCollision) logERROR("Found a collision, which is not supported.");
      

      startR += sectionWidth + safetySpace;

      // firstSection has been updated and now points to the Section which is the 'towards smaller |Z|' branch of the C.
      // lastSection has been updated and now points to the Section which is the 'towards bigger |Z|' branch of the C.
      // startR is updated, so that the services will leave at a reasonable radius.
      // startZ is not updated.
    }
  }


  /**
   * Look for the nearest collision in every defined boundary, starting from the (startZ, startR) point.
   * If no collision is found return the external perimeter coordinates
   * @param collision is a reference for the return value of the coordinate of the collision point (Z if horizontal, rho if vertical)
   * @param border is a reference for the return value of the coordinate of the border of boundary (rho if horizontal, Z if vertical)
   * @param startZ is the Z coordinate of the starting point
   * @param startR is the rho coordinate of the starting point
   * @param direction
   * @return true if a collision is found, false otherwise
   */
  bool Materialway::OuterUsher::findBoundaryCollision(int& collision, int& border, int startZ, int startR, const Tracker& tracker, Direction direction) {
    int hitCoord;
    int globalMaxZ = discretize(tracker.maxZwithHybrids()) + globalMaxZPadding;
    int globalMaxR = discretize(tracker.maxRwithHybrids()) + globalMaxRPadding;
    std::map<int, BoundariesSet::const_iterator> hitBoundariesCoords;
    bool foundCollision = false;

    //test all the boundaries for an hit (keep it on a map ordered for key = collision coords, so the first is the nearest collision)
    //for (Boundary& currBoundary : boundariesList_) {
    for(BoundariesSet::const_iterator it = boundariesList_.cbegin(); it != boundariesList_.cend(); ++it) {
      hitCoord = (*it)->isHit(startZ, startR, direction);
      if (hitCoord>0) {
        hitBoundariesCoords[hitCoord] = it;
      }
    }

    if (hitBoundariesCoords.size()) {
      collision = hitBoundariesCoords.begin()->first;
      if(direction == HORIZONTAL) {
        border = (*hitBoundariesCoords.begin()->second)->maxR();
      } else {
        border = (*hitBoundariesCoords.begin()->second)->maxZ();
      }
      foundCollision = true;
    } else {
      if(direction == HORIZONTAL) {
        collision = globalMaxZ;
        border = globalMaxR + safetySpace;
      } else {
        collision = globalMaxR;
        border = globalMaxZ + safetySpace;
      }
    }

    return foundCollision;
  }

  /**
   * Look for the nearest collision in every defined section, starting from the (startZ, startR) point.
   * @param sectionCollision is a reference for the return value of the coordinate of the collision point (Z if horizontal, rho if vertical)
   * and the pointer to the collided section
   * @param startZ is the Z coordinate of the starting point
   * @param startR is the rho coordinate of the starting point
   * @param end is the maximum range of searching for a collision, around (startZ, startR)
   * @param direction
   * @return true if a collision is found, false otherwise
   */
  bool Materialway::OuterUsher::findSectionCollision(std::pair<int,Section*>& sectionCollision, int startZ, int startR, int end, Direction direction) {
    int hitCoord;
    std::map<int, Section*> hitSectionsCoords;

    //test all the sections for an hit (keep it on a map ordered for key = collision coords, so the first is the nearest collision)
    for (Section* currSection : sectionsList_) {
      hitCoord = currSection->isHit(startZ, startR, end, direction);
      if (hitCoord>0) {
        hitSectionsCoords[hitCoord] = currSection;
      }
    }

    if (hitSectionsCoords.size()) {
      sectionCollision = std::make_pair(hitSectionsCoords.begin()->first, hitSectionsCoords.begin()->second);
      return true;
    }
    return false;
  }

  /**
   * Build a section segment from the start point to the end in given direction. Also update the value of the start
   * with the new start point
   * @param firstSection is a pointer for returning the first built section
   * @param lastSection is a pointer to the previous section, is used for updating the nextSection pointer to this
   * @param startZ is a reference to the Z start coordinate
   * @param startR is a reference to the rho start coordinate
   * @param end is a reference to the end coordinate (Z if horizontal, rho if vertical)
   * @param direction is the direction
   * @param isTowardsBiggerZ
   * @param isTowardsBiggerR  NB: !isTowardsBiggerR not used at all so far, added for flexibility purposes!
   * @return true if no section collision found, false otherwise
   */
  bool Materialway::OuterUsher::buildSection(Section*& firstSection, Section*& lastSection, int& startZ, int& startR, int end, Direction direction, 
					     bool isTowardsBiggerZ, bool isTowardsBiggerR) {
    int minZ, minR, maxZ, maxR;
    int trueEnd = end;
    int cutCoordinate;
    Section* newSection;
    bool foundSectionCollision = false;
    bool returnValue = true;
    std::pair<int,Section*> sectionCollision;
    int N_subsections;

    //search for collisions
    // Dont handle this if goes towards smaller |Z|, or smaller R. 
    // These are corner cases which are rare and handled independently (routeOutgoingServicesAlongCShape, or routeOutgoingServicesAlongUShape).
    if (isTowardsBiggerZ && isTowardsBiggerR) {
      foundSectionCollision = findSectionCollision(sectionCollision, startZ, startR, end, direction);
    }

    //if section collision found update the end
    if(foundSectionCollision) {
      trueEnd = sectionCollision.first - safetySpace;
      returnValue = false;
    }

    //set coordinates
    // HORIZONTAL DIRECTION
    if (direction == HORIZONTAL) {
      // TOWARDS BIGGER |Z|
      if (isTowardsBiggerZ) {
	minZ = startZ;      
	maxZ = trueEnd;
	startZ = trueEnd + safetySpace;
      }
      // TOWARDS SMALLER |Z|
      else {
	minZ = trueEnd;      
	maxZ = startZ;
	// Do not silently update startZ. Here we will just create a section with the desired coordinates.
      }
      minR = startR;
      maxR = startR + sectionWidth;    
    } 

    // VERTICAL DIRECTION
    else {
      // TOWARDS BIGGER R
      if (isTowardsBiggerR) {
	minR = startR;
	maxR = trueEnd;
	startR = trueEnd + safetySpace;
      }
      // TOWARDS SMALLER R
      else {
	minR = trueEnd;
	maxR = startR;
	// Do not silently update startR. Here we will just create a section with the desired coordinates.
      }
      minZ = startZ;    
      maxZ = startZ + sectionWidth;   
    }

    //build section and update last section's nextSection pointer
    if ((maxZ > minZ) && (maxR > minR)) { 
      //calculate number of subsections. If direction == HORIZONTAL, no sectioning
      N_subsections = (direction == HORIZONTAL) ? 1 : ceil( log( double(maxR) / minR ) / log( (1 + radialDistribError) / (1 - radialDistribError) ) );

      if (N_subsections == 1){    
	newSection = new Section(minZ, minR, maxZ, maxR, direction);
	sectionsList_.push_back(newSection);
	updateLastSectionPointer(lastSection, newSection);
	lastSection = newSection;
	//if firstSection is not yet set, set it
	if(firstSection == nullptr) {
	  firstSection = newSection;
	}
      }

      else {
	double ratioR = pow( double(minR) / maxR, 1 / double(N_subsections) );
	double maxR_temp = maxR; 
	double minR_temp = maxR * ratioR;
	for (int i=0; i<N_subsections; i++){
	  newSection = new Section(minZ, minR_temp, maxZ, maxR_temp, direction);
	  sectionsList_.push_back(newSection);
	  updateLastSectionPointer(lastSection, newSection);
	  lastSection = newSection;
	  maxR_temp = minR_temp;
	  minR_temp *= ratioR;
	  //if firstSection is not yet set, set it
	  if(firstSection == nullptr && i==0) {
	    firstSection = newSection;
	  }
	}
      } 
    }
    else { logERROR( Form( "While building sections of services, I found minZ=%d maxZ=%d minR=%d maxR=%d", minZ, maxZ, minR, maxR ) ); }

    //if found a collision and the section is perpendicular split the collided section and update nextSection pointer
    if(foundSectionCollision) {
      if (sectionCollision.second->bearing() != direction) {
        if (direction == HORIZONTAL) {
          cutCoordinate = minR;
        } else {
          cutCoordinate = minZ;
        }
        newSection = splitSection(sectionCollision.second, cutCoordinate, direction);
        //sectionsList_.push_back(newSection);
        updateLastSectionPointer(lastSection, newSection);
      } else {
        //set directly the pointer
        updateLastSectionPointer(lastSection, sectionCollision.second);
      }
    }

    return returnValue;
  }

  /**
   * Build a pair of section around a border in the shape of an L, the first section from the start point to the boundary collision,
   * the second section from the collision to the boundary border
   * @param firstSection is a pointer for returning the first built section
   * @param lastSection is a pointer to the previous section, is used for updating the nextSection pointer to this
   * @param startZ is a reference to the Z start coordinate
   * @param startR is a reference to the rho start coordinate
   * @param collision is a coordinate for the first collision point with boundary (Z if horizontal, rho if vertical)
   * @param border is a coordinate for the second point, the border of the boundary (rho if horizontal, Z if vertical)
   * @param direction is the direction
   * @return true if no section collision found, false otherwise
   */
  bool Materialway::OuterUsher::buildSectionPair(Section*& firstSection, Section*& lastSection, int& startZ, int& startR, int collision, int border, Direction direction) {

    bool noCollision = true;

    //build first section and if no collision happened build also the second section
    noCollision = buildSection(firstSection, lastSection, startZ, startR, collision - sectionWidth - safetySpace, direction);
    if (noCollision) {
      noCollision = buildSection(firstSection, lastSection, startZ, startR, border - safetySpace, inverseDirection(direction));
    }

    return noCollision;
  }

  Materialway::Section* Materialway::OuterUsher::splitSection(Section* section, int collision, Direction direction) {
    Section* retValue = nullptr;
    Section* useless = nullptr;
    int secMinZ = section->minZ();
    int secMinR = section->minR();
    int secCollision = collision;

      if (direction == HORIZONTAL) {
        buildSection(useless, retValue, secMinZ, secCollision, section->maxR(), inverseDirection(direction));

        section->maxR(collision - safetySpace);
        updateLastSectionPointer(section, retValue);
      } else {
        buildSection(useless, retValue, secCollision, secMinR, section->maxZ(), inverseDirection(direction));

        section->maxZ(collision - safetySpace);
        updateLastSectionPointer(section, retValue);
      }
      return retValue;
  }

  Materialway::Direction Materialway::OuterUsher::inverseDirection(Direction direction) const {
    return (direction == HORIZONTAL)?VERTICAL:HORIZONTAL;
  }

  void Materialway::OuterUsher::updateLastSectionPointer(Section* lastSection, Section* newSection) {
    if(lastSection != nullptr) {
      lastSection->nextSection(newSection);
    }
  }
  //END Materialway::OuterUsher
  //=================================================================================
  //START Materialway::InnerUsher
  Materialway::InnerUsher::InnerUsher(
      SectionVector& sectionsList,
      StationVector& stationListFirst,
      StationVector& stationListSecond,
      BarrelBoundaryMap& barrelBoundaryAssociations,
      EndcapBoundaryMap& endcapBoundaryAssociations,
      ModuleSectionMap& moduleSectionAssociations,
      LayerRodSectionsMap& layerRodSections,
      DiskRodSectionsMap& diskRodSections) :
        sectionsList_(sectionsList),
        stationListFirst_(stationListFirst),
        stationListSecond_(stationListSecond),
        barrelBoundaryAssociations_(barrelBoundaryAssociations),
        endcapBoundaryAssociations_(endcapBoundaryAssociations),
        moduleSectionAssociations_(moduleSectionAssociations),
        layerRodSections_(layerRodSections),
        diskRodSections_(diskRodSections) {}
  Materialway::InnerUsher::~InnerUsher() {}


  void Materialway::InnerUsher::go(Tracker& tracker) {


    //Visitor for the barrel layers
    class MultipleVisitor : public GeometryVisitor {
    public:
      MultipleVisitor(
          SectionVector& sectionsList,
          StationVector& stationListFirst,
          StationVector& stationListSecond,
          BarrelBoundaryMap& barrelBoundaryAssociations,
          EndcapBoundaryMap& endcapBoundaryAssociations,
          ModuleSectionMap& moduleSectionAssociations,
          LayerRodSectionsMap& layerRodSections,
          DiskRodSectionsMap& diskRodSections) :
        sectionsList_(sectionsList),
        stationListFirst_(stationListFirst),
        stationListSecond_(stationListSecond),
        barrelBoundaryAssociations_(barrelBoundaryAssociations),
        endcapBoundaryAssociations_(endcapBoundaryAssociations),
        moduleSectionAssociations_(moduleSectionAssociations),
        layerRodSections_(layerRodSections),
        diskRodSections_(diskRodSections),
        //startLayerZMinus(nullptr),
        startBarrel(nullptr),
	startLayer(nullptr),
        startEndcap(nullptr),
        startDisk(nullptr),
        currLayer_(nullptr),
        currDisk_(nullptr),
	hasCShapeOutgoingServices_(false) {}
        //currEndcapPosition(POSITIVE) {}//, splitCounter(0) {}
      virtual ~MultipleVisitor() {}


      //For the barrels ----------------------------------------------------
      void visit(Barrel& barrel) {
        //build section right to layers
        Boundary* boundary = barrelBoundaryAssociations_[&barrel];

        int minZ = discretize(barrel.maxZwithHybrids()) + layerSectionRightMargin + safetySpace + layerStationLenght;
        int minR = discretize(barrel.minRwithHybrids());
        int maxZ = minZ + sectionWidth;
        int maxR = boundary->maxR() - safetySpace;

        startBarrel = new Section(minZ, minR, maxZ, maxR, VERTICAL, boundary->outgoingSection());

	hasCShapeOutgoingServices_ = boundary->hasCShapeOutgoingServices();
      }


      void visit(Layer& layer) {
        ConversionStation* flangeConversionStation = layer.flangeConversionStation();
        const std::vector<ConversionStation*>& secondConversionStations = layer.secondConversionStations();

        currLayer_ = &layer;
        Section* section = nullptr;
        Station* station = nullptr;

        int attachPoint = 0;
        int sectionMinZ = 0;
        int sectionMinR = 0;
        int sectionMaxZ = 0;
        int sectionMaxR = 0;
        int stationMinZ = 0;
        int stationMinR = 0;
        int stationMaxZ = 0;
        int stationMaxR = 0;


        //split the right section
        section = startBarrel;
        attachPoint = discretize(layer.maxRwithHybrids()) + layerSectionMargin;        //discretize(layer.minRwithHybrids());

	
	while(section->maxR() < attachPoint + sectionTolerance) {
	  if(!section->hasNextSection()) {
	    logERROR("Error in finding attach point during construction of services");
	    return;
	  }
	  section = section->nextSection();
	}

	if (section->minR() < attachPoint - sectionTolerance) {
	  section = splitSection(section, attachPoint);
	}

	//built main sections above the layer

	sectionMinZ = safetySpace;
	sectionMinR = attachPoint;
	//sectionMaxZ = discretize(layer.maxZwithHybrids()) + layerSectionRightMargin;
	sectionMaxZ = section->minZ() - safetySpace - layerStationLenght;
	sectionMaxR = sectionMinR + sectionWidth;

	stationMinZ = sectionMaxZ + safetySpace;
	stationMinR = sectionMinR - (layerStationWidth / 2);
	stationMaxZ = sectionMaxZ + safetySpace + layerStationLenght;
	stationMaxR = sectionMinR + (layerStationWidth / 2);

	if(flangeConversionStation != nullptr) {
	  station = new Station(stationMinZ, stationMinR, stationMaxZ, stationMaxR, VERTICAL, *flangeConversionStation, section); //TODO: check if is ok VERTICAL
	  sectionsList_.push_back(station);
	  stationListFirst_.push_back(station);
	  startLayer = new Section(sectionMinZ, sectionMinR, sectionMaxZ, sectionMaxR, HORIZONTAL, station);
	  layerRodSections_[currLayer_].setStation(station);
	} else {
	  startLayer = new Section(sectionMinZ, sectionMinR, sectionMaxZ, sectionMaxR, HORIZONTAL, section);
	  layerRodSections_[currLayer_].setStation(section);
	  logERROR("Flange conversion not defined, bypassed.");
	}
	sectionsList_.push_back(startLayer);
	layerRodSections_[currLayer_].addSection(startLayer);

        
        //==========second level conversion station
	buildSecondConversionStations(section, secondConversionStations);
      }


      void visit(BarrelModule& module) {
        Section* section;
        int attachPoint;

        //if the module is in z plus
        if(module.maxZ() > 0) {
          if(module.maxZ() <= currLayer_->maxZwithHybrids()) {
            attachPoint = discretize(module.maxZ());
          } else {
            attachPoint = discretize(currLayer_->maxZwithHybrids());
          }
          section = startLayer;
          while (section->maxZ() < attachPoint + sectionTolerance) {
            if(!section->hasNextSection()) {
              logERROR("No section above the module found.");
              return;
            }
            section = section->nextSection();
          }
          if (section->minZ() < attachPoint - sectionTolerance) {
            section = splitSection(section, attachPoint);
            layerRodSections_[currLayer_].addSection(section);
          }
          moduleSectionAssociations_[&module] = section;
        }
        /*
        else {
          attachPoint = discretize(module.minZ());
          section = startLayerZMinus;
          while (section->minZ() > attachPoint - sectionTolerance) {
            if(!section->hasNextSection()) {
              //TODO: messaggio di errore
              return;
            }
            section = section->nextSection();
          }
          if (section->maxZ() > attachPoint + sectionTolerance) {
            section = splitSection(section, attachPoint, false);
          }
          moduleSectionAssociations_[&module] = section;
        }
        */
      }


      //For the endcaps ----------------------------------------------------
      void visit(Endcap& endcap) {
        if (endcap.minZwithHybrids() >= 0) {
          //currEndcapPosition = POSITIVE;

          //build section above the disks
          Boundary* boundary = endcapBoundaryAssociations_[&endcap];
          int minZ = discretize(endcap.minZwithHybrids());
          int minR = discretize(endcap.maxRwithHybrids())  + diskSectionUpMargin + safetySpace;
          int maxZ = boundary->maxZ() - safetySpace;
          int maxR = minR + sectionWidth;
          startEndcap = new Section(minZ, minR, maxZ, maxR, HORIZONTAL, boundary->outgoingSection());

	  hasCShapeOutgoingServices_ = boundary->hasCShapeOutgoingServices();
        }
      }


      void visit(Disk& disk) {
        Section* section = nullptr;
        Station* station = nullptr;

        ConversionStation* flangeConversionStation = disk.flangeConversionStation();        
        const std::vector<ConversionStation*>& secondConversionStations = disk.secondConversionStations();

        currDisk_ = &disk;
        //if(currEndcapPosition == POSITIVE) {
        if (disk.minZwithHybrids() >= 0) {
          //split the right section
          section = startEndcap;
          int attachPoint = discretize(disk.maxZwithHybrids()) + diskSectionMargin;

          while(section->maxZ() < attachPoint + sectionTolerance) {
            if(!section->hasNextSection()) {
              logERROR("No section over the disk module found.");
              return;
            }
            section = section->nextSection();
          }

          if (section->minZ() < attachPoint - sectionTolerance) {
            section = splitSection(section, attachPoint);
          }

          //built two main sections above the layer

	  // EXTERNAL SECTION ON THE RIGHT OF THE DISK
          int sectionMinZ = attachPoint;
          int sectionMinR = discretize(disk.minRwithHybrids());
          int sectionMaxZ = sectionMinZ + sectionWidth;
          //int sectionMaxR = discretize(disk.maxRwithHybrids()) + diskSectionUpMargin;
          int sectionMaxR = section->minR() - safetySpace - layerStationLenght;

	  // DEE VOLUMES (1 disk = 2 dees)
	  // TO DO: move this to tuneable disk properties.
	  const int deeZThickness = sectionWidth;
	  const int deeRadialDistanceToEdge = discretize(0.); // = 2 mm

	  // small |Z| disk
	  const int smallAbsZDiskMinZ = discretize(disk.smallAbsZDeeCenterZ()) - deeZThickness / 2.;
	  const int smallAbsZDiskMaxZ = discretize(disk.smallAbsZDeeCenterZ()) + deeZThickness / 2.;
          const int smallAbsZDiskMinR = discretize(disk.minRwithHybrids()) - deeRadialDistanceToEdge;
	  const int smallAbsZDiskMaxR = discretize(disk.maxRwithHybrids()) + deeRadialDistanceToEdge;

	  // big |Z| disk
	  const int bigAbsZDiskMinZ = discretize(disk.bigAbsZDeeCenterZ()) - deeZThickness / 2.;
	  const int bigAbsZDiskMaxZ = discretize(disk.bigAbsZDeeCenterZ()) + deeZThickness / 2.;
          const int bigAbsZDiskMinR = smallAbsZDiskMinR;
	  const int bigAbsZDiskMaxR = smallAbsZDiskMaxR;


	  // FLANGE STATION OVER THE DISK
          int stationMinZ = sectionMinZ - (layerStationWidth / 2);
          int stationMinR = sectionMaxR + safetySpace;
          int stationMaxZ = sectionMinZ + (layerStationWidth / 2);
          int stationMaxR = sectionMaxR + safetySpace + layerStationLenght;

	  Section* smallAbsZDisk = nullptr;
	  Section* bigAbsZDisk = nullptr;
          if(flangeConversionStation != nullptr) {
            station = new Station(stationMinZ, stationMinR, stationMaxZ, stationMaxR, HORIZONTAL, *flangeConversionStation, section);
            sectionsList_.push_back(station);
            stationListFirst_.push_back(station);
            diskRodSections_[currDisk_].setStation(station);

	    startDisk = new Section(sectionMinZ, sectionMinR, sectionMaxZ, sectionMaxR, VERTICAL, station, Location::EXTERNAL);
	    smallAbsZDisk = new Section(smallAbsZDiskMinZ, smallAbsZDiskMinR, smallAbsZDiskMaxZ, smallAbsZDiskMaxR, VERTICAL, station, Location::DEE);
	    bigAbsZDisk = new Section(bigAbsZDiskMinZ, bigAbsZDiskMinR, bigAbsZDiskMaxZ, bigAbsZDiskMaxR, VERTICAL, station, Location::DEE);
          } 
	  else {         
            diskRodSections_[currDisk_].setStation(section);

	    startDisk = new Section(sectionMinZ, sectionMinR, sectionMaxZ, sectionMaxR, VERTICAL, section, Location::EXTERNAL);
	    smallAbsZDisk = new Section(smallAbsZDiskMinZ, smallAbsZDiskMinR, smallAbsZDiskMaxZ, smallAbsZDiskMaxR, VERTICAL, section, Location::DEE);
	    bigAbsZDisk = new Section(bigAbsZDiskMinZ, bigAbsZDiskMinR, bigAbsZDiskMaxZ, bigAbsZDiskMaxR, VERTICAL, section, Location::DEE);
          }

          sectionsList_.push_back(startDisk);
	  sectionsList_.push_back(smallAbsZDisk);
	  sectionsList_.push_back(bigAbsZDisk);

          diskRodSections_[currDisk_].addSection(startDisk);
	  diskRodSections_[currDisk_].addSection(smallAbsZDisk);
	  diskRodSections_[currDisk_].addSection(bigAbsZDisk);
        

	  //==========second level conversion station
	  buildSecondConversionStations(section, secondConversionStations);
        }
      }


      void visit(EndcapModule& module) {
        //if(currEndcapPosition == POSITIVE) {
        if (module.minZ() >= 0) {
          Section* section;
          int attachPoint;

          attachPoint = discretize(module.maxR());
          section = startDisk;
          while (section->maxR() < attachPoint + sectionTolerance) {
            if(!section->hasNextSection()) {
              logERROR("No section next to module found.");
              return;
            }
            section = section->nextSection();
          }
          if (section->minR() < attachPoint - sectionTolerance) {
            section = splitSection(section, attachPoint);
            diskRodSections_[currDisk_].addSection(section);
          }

          moduleSectionAssociations_[&module] = section;
        }
      }


    private:
      /*
       * This is used to build the conversion stations volumes + make the succession of sections and stations correctly point to each other.
       */
      void buildSecondConversionStations(Section* section, const std::vector<ConversionStation*>& secondConversionStations) {

	// Loop on all cfg conversion stations
        for(ConversionStation* stationFromCfg : secondConversionStations) {
          bool isStationAlreadyBuilt = false;
          
          // Check whether the station is already built
          for (auto& existentStation : stationListSecond_) {
            if (existentStation->conversionStation().stationName_().compare(stationFromCfg->stationName_()) == 0) {
              isStationAlreadyBuilt = true;
              break;
            }
          }


	  // Only do sth if conversion station does not exist yet.
          if (!isStationAlreadyBuilt) {   

	    // STANDARD CASE: THERE IS NO C-SHAPE LOOP OVER THE BARREL / ENDCAP.
	    // STATIONS ARE BUILT TOWARDS BIGGER |Z|.
	    if (!hasCShapeOutgoingServices_) {

	      // PART A: FIND THE CABLING SECTION LOCATED BELOW THE CONVERSION STATION.

	      // Check whether station is placed towards bigger |Z|.
	      if (discretize(stationFromCfg->maxZ_()) < section->minZ()) {
		logERROR(any2str("Want to place conversion station (type second) ") + any2str(stationFromCfg->stationName_()) 
			 + any2str(" towards bigger |Z|. ")
			 + any2str("But we are at minZ = ") + any2str(section->minZ()/1000.) + any2str(" mm")
			 + any2str(" and station maxZ = ") + any2str(stationFromCfg->maxZ_()) + any2str(" mm.")
			 );
		continue;
	      }

	      const int attachPoint = discretize(stationFromCfg->meanZ());
        
	      // If we are at a too low |Z|, we go on from section to section, towards bigger |Z|.
	      while(section->maxZ() < attachPoint + sectionTolerance) {
		if(!section->hasNextSection()) {
		  logERROR(any2str("Want to place conversion station (type second) ") + any2str(stationFromCfg->stationName_()) 
			   + any2str(" towards bigger |Z|. ")
			   + any2str("But we are at maxZ = ") + any2str(section->maxZ()/1000.) + any2str(" mm")
			   + any2str(" and we have no next section (routing is stopping here).")
			   );
		  return;
		}
		section = section->nextSection();
	      }
     
	      // Now we are correctly placed at section = the section below the conversion station. :)

	      // PART B: BUILD THE 2 CABLING SECTIONS BELOW THE STATION.
	      if (section->minZ() < attachPoint - sectionTolerance) {
		splitSection(section, attachPoint); 
		// Split the exisiting section under the conversion station into 2 sections, around Z = attachPoint.
		// The section at smaller |Z| points towards the section at bigger |Z|.
	      }
	    }


	    // C-SHAPE LOOP OVER THE BARREL / ENDCAP.
	    // For example, C-shape loop built over TBPX.
	    // STATIONS ARE BUILT TOWARDS SMALLER |Z|.
	    // THEY ARE BUILT ALONG THE RADIALLY LOWER HORIZONTAL BRANCH OF THE C.
	    else {
	      const bool isTowardsBiggerZ = false;

	      // PART A: FIND THE CABLING SECTION LOCATED BELOW THE CONVERSION STATION.

	      // Check whether station is placed towards smaller |Z|.
	      if (discretize(stationFromCfg->minZ_()) > section->maxZ()) {
		logERROR(any2str("Want to place conversion station (type second) ") + any2str(stationFromCfg->stationName_()) 
			 + any2str(" towards smaller |Z|. ")
			 + any2str("But we are at maxZ = ") + any2str(section->maxZ()/1000.) + any2str(" mm")
			 + any2str(" and station minZ = ") + any2str(stationFromCfg->minZ_()) + any2str(" mm.")
			 );
		continue;
	      }

	      const int attachPoint = discretize(stationFromCfg->meanZ());
        
	      // If we are at a too big |Z|, we go on from section to section, towards smaller |Z|.
	      while(section->minZ() > attachPoint - sectionTolerance) {
		if(!section->hasNextSection()) {
		  logERROR(any2str("Want to place conversion station (type second) ") + any2str(stationFromCfg->stationName_()) 
			   + any2str(" towards smaller |Z|. ")
			   + any2str("But we are at minZ = ") + any2str(section->minZ()/1000.) + any2str(" mm")
			   + any2str(" and we have no next section (routing is stopping here).")
			   );
		  return;
		}
		section = section->nextSection();
	      }
	   
	      // Now we are correctly placed at section = the section below the conversion station. :)

	      // PART B: BUILD THE 2 CABLING SECTIONS BELOW THE STATION.
	      if (section->maxZ() > attachPoint + sectionTolerance) {
		splitSection(section, attachPoint, isTowardsBiggerZ);
		// Split the exisiting section under the conversion station into 2 sections, around Z = attachPoint.
		// The section at bigger |Z| points towards the section at smaller |Z|.
	      }
	    }


	    // PART C: BUILD CONVERSION STATION
	    const int stationMinZ = discretize(stationFromCfg->minZ_());
	    const int stationMaxZ = discretize(stationFromCfg->maxZ_());
	    const int stationMinR = section->maxR() + safetySpace;
	    const int stationMaxR = stationMinR + layerStationLenght;
	   
	    Station* myStation = nullptr;
            if(section->hasNextSection()) {
              myStation = new Station(stationMinZ, stationMinR, stationMaxZ, stationMaxR, HORIZONTAL, *stationFromCfg, section->nextSection());
            } else {
              myStation = new Station(stationMinZ, stationMinR, stationMaxZ, stationMaxR, HORIZONTAL, *stationFromCfg);
            }

            section->nextSection(myStation);


	    // PART D: STORE CONVERSION STATION.
            sectionsList_.push_back(myStation);
            stationListSecond_.push_back(myStation);
          }
        }
      }


      /*
       * Split section into 2 sections, at the collision coordinate.
       * NB: !isTowardsBiggerR not used at all so far, just added for flexibility purposes!
       */
      Section* splitSection(Section* section, const int collision, bool isTowardsBiggerZ = true, bool isTowardsBiggerR = true, bool debug = false) {
        Section* newSection = nullptr;

	// HORIZONTAL DIRECTION
        if (section->bearing() == HORIZONTAL) { 
	  // A: Create new section
	  // TOWARDS BIGGER |Z|
	  if (isTowardsBiggerZ) { 
	    newSection = new Section(collision, section->minR(), section->maxZ(), section->maxR(), section->bearing(), section->nextSection(), debug); 
	  }
	  // TOWARDS SMALLER |Z|
	  else { 
	    newSection = new Section(section->minZ(), section->minR(), collision, section->maxR(), section->bearing(), section->nextSection(), debug);
	  }
          sectionsList_.push_back(newSection);

	  // B: Update previous section
	  if (isTowardsBiggerZ) { section->maxZ(collision - safetySpace); } 
	  else { section->minZ(collision + safetySpace); }
          section->nextSection(newSection);
        } 

	// VERTICAL DIRECTION
	else {
	  // A: Create new section
	  // TOWARDS BIGGER R
	  if (isTowardsBiggerR) {
	    newSection = new Section(section->minZ(), collision, section->maxZ(), section->maxR(), section->bearing(), section->nextSection(), debug);	    
	  }
	  // TOWARDS SMALLER R
	  else {
	    newSection = new Section(section->minZ(), section->minR(), section->maxZ(), collision, section->bearing(), section->nextSection(), debug);
	  }
	  sectionsList_.push_back(newSection);

	  // B: Update previous section
          if (isTowardsBiggerR) { section->maxR(collision - safetySpace); }
	  else { section->minR(collision + safetySpace); }
          section->nextSection(newSection);
        }
        return newSection;
      }
      

      enum ZPosition {POSITIVE, NEGATIVE};
      //int splitCounter;
      SectionVector& sectionsList_;
      StationVector& stationListFirst_;
      StationVector& stationListSecond_;
      BarrelBoundaryMap& barrelBoundaryAssociations_;
      EndcapBoundaryMap& endcapBoundaryAssociations_;
      ModuleSectionMap& moduleSectionAssociations_;
      LayerRodSectionsMap& layerRodSections_;  //map each layer with corresponding sections vector and station
      DiskRodSectionsMap& diskRodSections_;
      Section* startBarrel;
      Section* startLayer;
      //Section* startLayerZMinus;
      Section* startEndcap;
      Section* startDisk;     
      //ZPosition currEndcapPosition;
      Layer* currLayer_;
      Disk* currDisk_;
      bool hasCShapeOutgoingServices_;
    }; //END class MultipleVisitor

    MultipleVisitor visitor (sectionsList_, stationListFirst_, stationListSecond_, barrelBoundaryAssociations_, endcapBoundaryAssociations_, moduleSectionAssociations_, layerRodSections_, diskRodSections_);
    tracker.accept(visitor);
  }

  //END Materialway::InnerUsher
  //=================================================================================
  //START Materialway

  const double Materialway::gridFactor = 1000.0;                                     /**< the conversion factor for using integers in the algorithm (helps finding collisions),
                                                                              actually transforms millimiters in microns */
  const int Materialway::sectionWidth = discretize(insur::geom_inactive_volume_width);     /**< the width of a section */
  const int Materialway::safetySpace = discretize(insur::geom_epsilon);           /**< the safety space between sections */
  //const double Materialway::globalMaxZ_mm = insur::max_length;                     /**< the Z coordinate of the end point of the sections */
  //const double Materialway::globalMaxR_mm = insur::outer_radius;                   /**< the rho coordinate of the end point of the sections */
  //const int Materialway::globalMaxZ = discretize(globalMaxZ_mm);
  //const int Materialway::globalMaxR = discretize(globalMaxR_mm);

  const int Materialway::boundaryPaddingBarrel = discretize(12.0);             /**< the space between the barrel/endcap and the containing box (for routing services) */
  const int Materialway::boundaryPaddingEndcaps = discretize(10.0); 
  const int Materialway::boundaryPrincipalPaddingBarrel = discretize(21.0);       /**< the space between the barrel/endcap and the containing box only right for the barrel, up for endcap */  
  const int Materialway::globalMaxZPadding = discretize(100.0);          /**< the space between the tracker and the right limit (for routing services) */
  const int Materialway::globalMaxRPadding = discretize(30.0);          /**< the space between the tracker and the upper limit (for routing services) */

  const int Materialway::layerSectionMargin = discretize(2.0);          /**< the space between the layer and the service sections over it */
  const int Materialway::diskSectionMargin = discretize(2.0);          /**< the space between the disk and the service sections right of it */
  const int Materialway::layerSectionRightMargin = discretize(5.0);     /**< the space between the end of the layer (on right) and the end of the service sections over it */
  const int Materialway::diskSectionUpMargin = discretize(5.0);     /**< the space between the end of the disk (on top) and the end of the service sections right of it */
  const int Materialway::sectionTolerance = discretize(1.0);       /**< the tolerance for attaching the modules in the layers and disk to the service section next to it */
  const int Materialway::cShapeMinZTolerance = 1.1 * sectionTolerance; /** Distance in Z between the minZ of the cabling along C-shape, 
									   and the center of the conversion station (of type second). 
									   If there are several conversion stations (of type second) on the boundary, 
									   consider the station located at the lowest |Z|.
									   WARNING: By construction, cShapeMinZTolerance need to be > sectionTolerance. */
  const int Materialway::layerStationLenght = discretize(insur::geom_conversion_station_width);  /**< the lenght of the converting station on right of the layers */
  const int Materialway::layerStationWidth = discretize(20.0);         /**< the width of the converting station on right of the layers */
  const double Materialway::radialDistribError = 0.05;                 /**< 5% max error in the material radial distribution */

  Materialway::Materialway() :
    boundariesList_(),
    outerUsher(sectionsList_, boundariesList_),
    innerUsher(sectionsList_, stationListFirst_, stationListSecond_, barrelBoundaryAssociations_, endcapBoundaryAssociations_, moduleSectionAssociations_, layerRodSections_, diskRodSections_)
  {}
  Materialway::~Materialway() {}

  int Materialway::discretize(double input) {
    return int(input * gridFactor);
  }
  double Materialway::undiscretize(int input) {
    return double(input / gridFactor);
  }

  bool Materialway::build(Tracker& tracker, InactiveSurfaces& inactiveSurface, WeightDistributionGrid& weightDistribution) {
    bool retValue = false;

    //int startTime = time(0);
    startTaskClock("Building boundaries");
    if (buildBoundaries(tracker)) {
      stopTaskClock();
      startTaskClock("Building external sections"); buildExternalSections(tracker); stopTaskClock();
      startTaskClock("Building internal sections"); buildInternalSections(tracker); stopTaskClock();
    } else stopTaskClock();

    startTaskClock("Building inactive elements"); buildInactiveElements(); stopTaskClock();
    startTaskClock("Rounting services"); routeServices(tracker); stopTaskClock();
    startTaskClock("First step conversions"); firstStepConversions(); stopTaskClock();
    startTaskClock("Second step conversions"); secondStepConversions(); stopTaskClock();
    startTaskClock("Creating ModuleCaps"); createModuleCaps(tracker); stopTaskClock();
    startTaskClock("Duplicating sections"); duplicateSections(); stopTaskClock();
    startTaskClock("Populating MaterialProperties"); populateAllMaterialProperties(tracker, weightDistribution); stopTaskClock();
    startTaskClock("Building inactive surfaces"); buildInactiveSurface(tracker, inactiveSurface); stopTaskClock();
    startTaskClock("Computing material amounts"); calculateMaterialValues(inactiveSurface, tracker); stopTaskClock();
    return retValue;
  }

  bool Materialway::buildBoundaries(const Tracker& tracker) {
    bool retValue = false;

    class BarrelVisitor : public ConstGeometryVisitor {
    public:
      BarrelVisitor(BoundariesSet& boundariesList, BarrelBoundaryMap& barrelBoundaryAssociations, EndcapBoundaryMap& endcapBoundaryAssociations) :
        boundariesList_(boundariesList),
        barrelBoundaryAssociations_(barrelBoundaryAssociations),
        endcapBoundaryAssociations_(endcapBoundaryAssociations) {}
      virtual ~BarrelVisitor() {}

      void visit(const Barrel& barrel) {
        int boundMinZ = discretize(barrel.minZwithHybrids()) - boundaryPaddingBarrel;
        int boundMinR = discretize(barrel.minRwithHybrids()) - boundaryPaddingBarrel;
        int boundMaxZ = discretize(barrel.maxZwithHybrids()) + boundaryPrincipalPaddingBarrel;
        int boundMaxR = discretize(barrel.maxRwithHybrids()) + boundaryPaddingBarrel;

	const int barrelSecondStationsMinMeanZ = computeBarrelSecondStationsMinMeanZ(barrel);
	const bool hasCShapeOutgoingServices = (barrelSecondStationsMinMeanZ < boundMaxZ); // KEY POINT: Need to route services along a C-shape.
	const int cShapeMinZ = barrelSecondStationsMinMeanZ - cShapeMinZTolerance;
	Boundary* newBoundary = new Boundary(&barrel, true, boundMinZ, boundMinR, boundMaxZ, boundMaxR, hasCShapeOutgoingServices, cShapeMinZ);

        boundariesList_.insert(newBoundary);
        barrelBoundaryAssociations_.insert(std::make_pair(&barrel, newBoundary));
      }

      void visit(const Endcap& endcap) {
        const int boundMinZ = discretize(endcap.minZwithHybrids()) - boundaryPaddingEndcaps;
        const int boundMinR = discretize(endcap.minRwithHybrids()) - boundaryPaddingEndcaps;
        const int boundMaxZ = discretize(endcap.maxZwithHybrids()) + discretize(endcap.distanceFromEndcapsModulesMaxZtoRoutedMaterial());
        const int boundMaxR = discretize(endcap.maxRwithHybrids()) + discretize(endcap.distanceFromEndcapsModulesMaxRtoRoutedMaterial());
        Boundary* newBoundary = new Boundary(&endcap, false, boundMinZ, boundMinR, boundMaxZ, boundMaxR);

        boundariesList_.insert(newBoundary);
        endcapBoundaryAssociations_.insert(std::make_pair(&endcap, newBoundary));
      }


    private:
      /*
       * This is used to compute the minimum meanZ of the barrel second conversion stations.
       */
      const int computeBarrelSecondStationsMinMeanZ(const Barrel& barrel) {
	int barrelSecondStationsMinMeanZ = std::numeric_limits<int>::max();
	
	for (const auto& layer: barrel.layers()) {
	  const std::vector<ConversionStation*> secondConversionStations = layer.secondConversionStations();
	  for (const auto& station : secondConversionStations) {
	    barrelSecondStationsMinMeanZ = MIN(barrelSecondStationsMinMeanZ, 
					       discretize(station->meanZ())
					       );
	  }
	}
	return barrelSecondStationsMinMeanZ;
      }

      BoundariesSet& boundariesList_;
      BarrelBoundaryMap& barrelBoundaryAssociations_;
      EndcapBoundaryMap& endcapBoundaryAssociations_;
    };


    BarrelVisitor visitor (boundariesList_, barrelBoundaryAssociations_, endcapBoundaryAssociations_);
    tracker.accept(visitor);

    if (boundariesList_.size() > 0) {

      // For a given tracker, in case of a step in the Outer Radius in the Endcaps : 
      // This sets the maxR boundary of the outermost barrel to the value of the maxR of the closest-to-barrel Endcap.
      // As a result, the services will be routed properly, directly from the outermost barrel to the closest-to-barrel Endcap.
      if (tracker.hasStepInEndcapsOuterRadius() && barrelBoundaryAssociations_.size() > 0 && endcapBoundaryAssociations_.size() > 0) {
	
	// In boundariesList_, the boundaries are sorted by decreasing order of the value of (boundary->maxR() + boundary->maxZ()).
	// As a result, one can define a comparaison relation between barrels (or between endcaps) : the comparaison relation that exists between their corresponding boundaries.
	// This allows to get the "biggest barrel" and "smallest endcap" in the sense of this relation of comparaison. As a result, we have access to the outermost barrel, and the closest-to-barrel endcap, respectively.
	
	// This is TEMPORARY. A relation of comparaison for barrels and endcaps should be introduced much upper in the code.
	// + Ugly to have boundariesList_, barrelBoundaryAssociations_, endcapBoundaryAssociations_ all storing the same data !
	// Convenient because allows to have both a relation of comparaison (1 set) and a sortage by key (2 maps), but should be optimized !!


	// Get the closest-to-barrel Endcap boundary (The Endcap which has the smallest value = (boundary->maxR() + boundary->maxZ()))
	auto smallestEndcapBoundary = std::find_if(boundariesList_.rbegin(), boundariesList_.rend(), 
										  [&](const Boundary* b) { return !b->isBarrel(); } );
	int smallestEndcapMaxR = (*smallestEndcapBoundary)->maxR(); // takes its maxR

	// Get the outermost Barrel boundary (The Barrel which has the biggest value = (boundary->maxR() + boundary->maxZ()))
	auto biggestBarrelBoundary = std::find_if(boundariesList_.begin(), boundariesList_.end(), 
									 [&](const Boundary* b) { return b->isBarrel(); } );

	// set the maxR for the outermost Barrel boundary, which is in barrelBoundaryAssociations_ 
	for (auto& it : barrelBoundaryAssociations_) {
	  if (*biggestBarrelBoundary == it.second) it.second->maxR(smallestEndcapMaxR);
	}
	// set the maxR for the outermost Barrel boundary, which is also in boundariesList_ ...
	(*biggestBarrelBoundary)->maxR(smallestEndcapMaxR);
      }
      retValue = true;
    }

    return retValue;
  }


  void Materialway::buildExternalSections(const Tracker& tracker) {
    for(BoundariesSet::iterator it = boundariesList_.begin(); it != boundariesList_.end(); ++it) {
      outerUsher.go(const_cast<Boundary*>(*it), tracker);
    }
  }


  void Materialway::buildInternalSections(Tracker& tracker) {
   innerUsher.go(tracker);
  }

  void Materialway::buildInactiveElements() {
    double zLength, zOffset, innerRadius, rWidth;

    for(Section* section : sectionsList_) {
    //for(SectionVector::iterator sectionIter = sectionsList_.begin(); sectionIter != sectionsList_.end(); ++sectionIter) {
      if(section->debug_ == false) {
        zLength = undiscretize(section->maxZ() - section->minZ());
        zOffset = undiscretize(section->minZ());
        innerRadius = undiscretize(section->minR());
        rWidth = undiscretize(section->maxR() - section->minR());

        if (section->bearing() == HORIZONTAL) {
          InactiveTube* tube = new InactiveTube;
          tube->setZLength(zLength);
          tube->setZOffset(zOffset);
          tube->setInnerRadius(innerRadius);
          tube->setRWidth(rWidth);
          tube->setFinal(true);
          section->inactiveElement(tube);
          //tube->addLocalMass("Steel", 1000.0*zLength);
          //inactiveSurface.addBarrelServicePart(tube);
        } else {
          InactiveRing* ring = new InactiveRing;
          ring->setZLength(zLength);
          ring->setZOffset(zOffset);
          ring->setInnerRadius(innerRadius);
          ring->setRWidth(rWidth);
          ring->setFinal(true);
          section->inactiveElement(ring);
          //ring->addLocalMass("Steel", 1000.0*rWidth);
          //inactiveSurface.addBarrelServicePart(ring);
        }
      }
    }
  }

  void Materialway::routeServices(const Tracker& tracker) {
    class ServiceVisitor : public ConstGeometryVisitor {
    private:
      ModuleSectionMap& moduleSectionAssociations_;  /**< Map that associate each module with the section that it feeds */
      LayerRodSectionsMap& layerRodSections_;      /**< maps for sections of the rods */
      DiskRodSectionsMap& diskRodSections_;

      const Layer* currLayer_;
      const Disk* currDisk_;

      bool printGuard;
      int printCounter;
      bool firstRod;
      bool firstRing;
      int rodSectionsSize;
      
      
      const std::vector<std::string> unitsToPassRodGGM = {"g", "g/m"};
      const std::vector<std::string> unitsToPassRodGM = {"g/m"};
      const std::vector<std::string> unitsToPassRodMM = {"mm"};
      const std::vector<std::string> unitsToPassLayer = {"g", "g/m", "mm"};
      const std::vector<std::string> unitsToPassLayerServ = {"g/m", "mm"};
    public:
      ServiceVisitor(ModuleSectionMap& moduleSectionAssociations, LayerRodSectionsMap& layerRodSections, DiskRodSectionsMap& diskRodSections) :
        moduleSectionAssociations_(moduleSectionAssociations),
        layerRodSections_(layerRodSections),
        diskRodSections_(diskRodSections), printGuard(true), printCounter(0), firstRing(false), rodSectionsSize(0) {}

      void visit(const Layer& layer) {
        currLayer_ = &layer;
        firstRod = true;
        if(layer.maxZwithHybrids() > 0.) {
          int totalLength = 0;
          for (Section* currSection : layerRodSections_.at(currLayer_).getSections()) {
            totalLength += currSection->maxZ() - currSection->minZ();
          }

          for (Section* currSection : layerRodSections_.at(currLayer_).getSections()) {
            layer.materialObject().deployMaterialTo(currSection->materialObject(), unitsToPassLayer, MaterialObject::SERVICES_AND_LOCALS, double(currSection->maxZ()-currSection->minZ()) / totalLength);
          }
          layerRodSections_.at(currLayer_).getStation()->getServicesAndPass(layer.materialObject(), unitsToPassLayerServ);
        }
      }

      void visit(const RodPair& rod) {
        if (firstRod) {
          rodSectionsSize = 0;
          for (Section* currSection : layerRodSections_.at(currLayer_).getSections()) {
            rod.materialObject().deployMaterialTo(currSection->materialObject(), unitsToPassRodMM);
            rodSectionsSize += currSection->maxZ() - currSection->minZ();
          }
          firstRod = false;
          layerRodSections_.at(currLayer_).getStation()->getServicesAndPass(rod.materialObject(), unitsToPassRodMM);
        }
        for (Section* currSection : layerRodSections_.at(currLayer_).getSections()) {
          rod.materialObject().deployMaterialTo(currSection->materialObject(), unitsToPassRodGGM, MaterialObject::SERVICES_AND_LOCALS, double(currSection->maxZ()-currSection->minZ()) / rodSectionsSize);
        }
        layerRodSections_.at(currLayer_).getStation()->getServicesAndPass(rod.materialObject(), unitsToPassRodGM);
      }

      void visit(const BarrelModule& module) {
        if(module.maxZ() > 0) {
          moduleSectionAssociations_.at(&module)->getServicesAndPass(module.materialObject());

          return;

          if (printGuard) {
            std::cout << "MODULO:\t\t< "
                      << std::setw(9) << module.minZ() << ";\tv " 
                      << std::setw(9) << module.minR() << ";\t> " 
                      << std::setw(9) << module.maxZ() << ";\t^ " 
                      << std::setw(9) << module.maxR() << std::endl;
            std::cout << "popolo sez mod:\t< " 
                      << std::setw(9) << undiscretize(moduleSectionAssociations_.at(&module)->minZ()) << ";\tv " 
                      << std::setw(9) << undiscretize(moduleSectionAssociations_.at(&module)->minR()) << ";\t> " 
                      << std::setw(9) << undiscretize(moduleSectionAssociations_.at(&module)->maxZ()) << ";\t^ " 
                      << std::setw(9) << undiscretize(moduleSectionAssociations_.at(&module)->maxR()) << std::endl;
          }

          //route layer rod services

          /* NO
          for (
               currSectionIter = layerRodSections_.at(currLayer_).getSections().begin();
               currSectionIter != layerRodSections_.at(currLayer_).getSections().end();
               ++currSectionIter) {
            currLayer_->materialObject().copyServicesTo((*currSectionIter)->materialObject());
            currLayer_->materialObject().copyLocalsTo((*currSectionIter)->materialObject());
          }
          */
         
          /* SI
          RodSectionsStation& rodSectionsStation = layerRodSections_.at(currLayer_);
          std::vector<Section*>& layerSections = rodSectionsStation.getSections();
          for (currSectionIter = layerSections.begin(); currSectionIter != layerSections.end(); ++currSectionIter) {
            Section* currSection = *currSectionIter;
            currLayer_->materialObject().copyServicesTo(currSection->materialObject());
            currLayer_->materialObject().copyLocalsTo(currSection->materialObject());
          }
          */

            /* BOH
          for (Section* currSection : layerRodSections_.at(currLayer_).getSections()) {
            currLayer_->materialObject().copyServicesTo(currSection->materialObject());
            currLayer_->materialObject().copyLocalsTo(currSection->materialObject());
            if (printGuard) {
              std::cout << "popolo sez rod:\t< " 
                        << std::setw(9) << undiscretize(currSection->minZ()) << ";\tv " 
                        << std::setw(9) << undiscretize(currSection->minR()) << ";\t> " 
                        << std::setw(9) << undiscretize(currSection->maxZ()) << ";\t^ " 
                        << std::setw(9) << undiscretize(currSection->maxR()) << std::endl;
            }
          }
           */

          /* NO
          RodSectionsStation& rodSectionsStation = layerRodSections_.at(currLayer_);
          std::vector<Section*>& layerSections = rodSectionsStation.getSections();
          for (Section* currSection : layerSections) {
            currLayer_->materialObject().copyServicesTo(currSection->materialObject());
            currLayer_->materialObject().copyLocalsTo(currSection->materialObject());
          }
          */

          //layerRodSections_.at(currLayer_).getStation()->getServicesAndPass(currLayer_->materialObject());
          
          if (printGuard) {
            std::cout << "popolo staz:\t< " 
                      << std::setw(9) << undiscretize(layerRodSections_.at(currLayer_).getStation()->minZ()) << ";\tv " 
                      << std::setw(9) << undiscretize(layerRodSections_.at(currLayer_).getStation()->minR()) << ";\t> " 
                      << std::setw(9) << undiscretize(layerRodSections_.at(currLayer_).getStation()->maxZ()) << ";\t^ " 
                      << std::setw(9) << undiscretize(layerRodSections_.at(currLayer_).getStation()->maxR()) << std::endl;
            
            if (++ printCounter == 2) {
              printGuard = false;
            }
          }
        }
      }

      void visit(const Disk& disk) {
        currDisk_ = &disk;

        if(disk.maxZwithHybrids() > 0) {
          firstRing = true;

	  // COMPUTE TOTAL VOLUME OF A DOUBLE-DISK
	  double totalDeesVolume = 0; 
          double totalExternalVolume = 0; 
	  // This loops on all the sections making up the double disk.
          for (Section* currSection : diskRodSections_.at(currDisk_).getSections()) {
	    if (currSection->getLocation() == Location::DEE) {
	      totalDeesVolume += currSection->getVolume();
	    }
	    else if (currSection->getLocation() == Location::EXTERNAL) {
	      totalExternalVolume += currSection->getVolume();
	    }
	    else {
	      logERROR("Double-Disk: created a section of unsupported location type.");
	    }
	  }

	  // ASSIGN MATERIALS OF THE DOUBLE-DISK (DEFINED IN CFG) TO THE SECTIONS (DEFINED IN MARTINA'S CODE).
	  // This loops on all the sections making up the double disk.
	  // The full point here is that the sections have shapes already, but no material affected to them yet!!
	  for (Section* currSection : diskRodSections_.at(currDisk_).getSections()) {
	    // Ratio of volume of the section divided by the total volume of the double-disk.
	    // This is a fix to Martina's code!!
	    // We want to use volume ratio, because we want to keep density uniform within the double-disk.
	    const double totalVolume = (currSection->getLocation() == Location::DEE ? totalDeesVolume : totalExternalVolume);

	    const double massRatio = currSection->getVolume() / totalVolume;
	    disk.materialObject().deployMaterialTo(currSection->materialObject(), 
						   unitsToPassLayer, 
						   MaterialObject::SERVICES_AND_LOCALS, 
						   massRatio,
						   currSection->getLocation()
						   );
	    // disk.materialObject() is the materials info from cfg files (all double-disk).
	    // currSection->materialObject() is the small section to which we want to assign a fraction of the materials.
	  }

	  // ROUTE SERVICES FROM DOUBLE-DISK   
	  diskRodSections_.at(currDisk_).getStation()->getServicesAndPass(disk.materialObject(), 
									  unitsToPassLayerServ);
	}

      }

      void visit(const Ring& ring) {
        if(firstRing) {
          //for the mm
          // for (Section* currSection : diskRodSections_.at(currDisk_).getSections()) {
          //   ring.materialObject().deployMaterialTo(currSection->materialObject(), unitsToPassRodMM);
          // }
          //diskRodSections_.at(currDisk_).getStation()->getServicesAndPass(ring.materialObject(), unitsToPassRodMM);
          //for the g/m
       
          // for(int i=0; i<ring.modules().size(); ++i) {
          //   for (Section* currSection : diskRodSections_.at(currDisk_).getSections()) {
          //     ring.materialObject().deployMaterialTo(currSection->materialObject(), unitsToPassRodGM);
          //   }
          //   diskRodSections_.at(currDisk_).getStation()->getServicesAndPass(ring.materialObject(), unitsToPassRodGM);
          // }

          //TODO: ADD A WARNING IF MATERIAL DEFINED FOR RODS
      
          firstRing = false;
        }
      }
      /*
      void visit(const Ring& ring) {
        //route disk rod services
        if(ring.minZwithHybrids() > 0) {
          for (Section* currSection : diskRodSections_.at(currDisk_).getSections()) {
            ring.materialObject().copyServicesTo(currSection->materialObject());
            ring.materialObject().copyLocalsTo(currSection->materialObject());
          }
          diskRodSections_.at(currDisk_).getStation()->getServicesAndPass(ring.materialObject());
        }
      }
      */

      void visit(const EndcapModule& module) {
        if (module.minZ() >= 0) {
          //route module services
          moduleSectionAssociations_.at(&module)->getServicesAndPass(module.materialObject());

          /*
          //route disk rod services
          for (Section* currSection : diskRodSections_.at(currDisk_).getSections()) {
            currDisk_->materialObject().copyServicesTo(currSection->materialObject());
            currDisk_->materialObject().copyLocalsTo(currSection->materialObject());
          }
          diskRodSections_.at(currDisk_).getStation()->getServicesAndPass(currDisk_->materialObject());
          */
        }
      }
    };

    ServiceVisitor v(moduleSectionAssociations_, layerRodSections_, diskRodSections_);
    tracker.accept(v);
  }

  /*

  void Materialway::routeModuleServices() {
    for (std::pair<const DetectorModule* const, Section*>& pair : moduleSectionAssociations_) {
      //pair.first->materialObject().routeServicesTo(pair.second->materialObject());
      pair.second->getServicesAndPass(pair.first->materialObject());
    }
  }

  void Materialway::routeRodMaterials() {
    bool first = true;
    for (std::pair<const Layer* const, RodSectionsStation> currAssociation : layerRodSections_) {
      first = true;
      for (Section* currSection : currAssociation.second.getSections()) {
        //currAssociation.first->materialObject().copyServicesTo(currSection->materialObject());
        if (first) {
          currSection->getServicesAndPass(currAssociation.first->materialObject());
          first = false;
        }
        currAssociation.first->materialObject().copyLocalsTo(currSection->materialObject());
      }
      currAssociation.second.getStation()->getServicesAndPass(currAssociation.first->materialObject());
    }

    for (std::pair<const Disk* const, RodSectionsStation> currAssociation : diskRodSections_) {
      for (Section* currSection : currAssociation.second.getSections()) {
        currAssociation.first->materialObject().copyServicesTo(currSection->materialObject());
        currAssociation.first->materialObject().copyLocalsTo(currSection->materialObject());
      }
      currAssociation.second.getStation()->getServicesAndPass(currAssociation.first->materialObject());
    }
  }

  */

  void Materialway::firstStepConversions() {
    for (Station* station : stationListFirst_) {
      if ((station->nextSection() != nullptr) && (station->inactiveElement() != nullptr)) {
        //put local materials in the materialObject of the station
        //put routed services in the materialObject of the adiacent section of station
        station->conversionStation().routeConvertedElements(station->materialObject(), station->outgoingMaterialObject(), *station->inactiveElement());
        if (station->nextSection() != nullptr) {
          //route converted materials
          station->nextSection()->getServicesAndPass(station->outgoingMaterialObject());
        }
      }
    }
  }

  void Materialway::secondStepConversions() {
    for (Station* station : stationListSecond_) {
      if ((station->nextSection() != nullptr) && (station->inactiveElement() != nullptr)) {
        //put local materials in the materialObject of the station
        //put routed services in the materialObject of the adiacent section of station
        station->conversionStation().routeConvertedElements(station->materialObject(), station->outgoingMaterialObject(), *station->inactiveElement());
        if (station->nextSection() != nullptr) {
          //route converted materials
          station->nextSection()->getServicesAndPass(station->outgoingMaterialObject());
        }
      }
    }
  }

  void Materialway::createModuleCaps(Tracker& tracker) {

    class CapsVisitor : public GeometryVisitor {
    private:
      int m_layerID;
      int m_discID;
      std::string m_detName;
    public:
      //std::map<BarrelModule*, int> mappaB;
      //std::map<EndcapModule*, int> mappaE;
      //int totB = 0;
      //int totE = 0;
      CapsVisitor() {m_layerID=-1; m_discID=-1;}
      virtual ~CapsVisitor() {};
      void visit(Barrel& b)
      {
        m_detName = b.myid();
      }
      void visit(Endcap& e)
      {
        m_detName = e.myid();
      }
      void visit(Layer& l)
      {
        m_layerID = l.myid();
        m_discID  = -1;
      }
      void visit(Disk& d)
      {
        m_layerID = -1;
        m_discID  = d.myid();
      }
      void visit(BarrelModule& m) {
        ModuleCap* cap = new ModuleCap(m,m_layerID);
        cap->setCategory(MaterialProperties::b_mod);
        cap->setDetName(m_detName);
        //mappaB[&m] = 0;
        //totB ++;
      }
      void visit(EndcapModule& m) {
        ModuleCap* cap = new ModuleCap(m,m_discID);
        cap->setCategory(MaterialProperties::e_mod);
        cap->setDetName(m_detName);
        //mappaE[&m] = 0;
        //totE ++;
      }
    };

    CapsVisitor v;
    tracker.accept(v);
    //std::cout << "Barrel " << v.mappaB.size() << " = " << v.totB << "  ----  Endcap " << v.mappaE.size() << " = " << v.totE << std::endl;
  }
 
  void Materialway::duplicateSections() {

    SectionVector negativeSections;
    
    for (Section* section : sectionsList_) {
      negativeSections.push_back(new Section(*section));
    }
    
    /*
    std::for_each(sectionsList_.begin(), sectionsList_.end(), [negativeSections](Section* section) mutable {
        negativeSections.push_back(new Section(*section));
      });
    */  
    std::for_each(negativeSections.begin(), negativeSections.end(), [](Section* section){
        section->minZ(-1 * section->minZ());
        section->maxZ(-1 * section->maxZ());
        section->inactiveElement()->setZOffset(-1 * section->inactiveElement()->getZOffset() - section->inactiveElement()->getZLength());
      });

    sectionsList_.reserve(sectionsList_.size() + negativeSections.size());
    sectionsList_.insert(sectionsList_.end(), negativeSections.begin(), negativeSections.end());
  }

  void Materialway::populateAllMaterialProperties(Tracker& tracker, WeightDistributionGrid& weightDistribution) {
    //sections
    for(Section* section : sectionsList_) {
      if(section->inactiveElement() != nullptr) {
        //section->inactiveElement()->addLocalMass("Steel", 1000.0*section->inactiveElement()->getZLength());

        section->materialObject().populateMaterialProperties(*section->inactiveElement());
        /*
        double sectionMinZ = undiscretize(section->minZ());
        double sectionMinR = undiscretize(section->minR());
        double sectionMaxZ = undiscretize(section->maxZ());
        double sectionMaxR = undiscretize(section->maxR());
        double sectionLength = undiscretize(section->maxZ() - section->minZ());
        double sectionArea = sectionLength * 2 * M_PI * sectionMinR;
        weightDistribution.addTotalGrams(sectionMinZ, sectionMinR, sectionMaxZ, sectionMaxR, sectionLength, sectionArea, section->materialObject());
        */
      } else {
        logUniqueERROR(inactiveElementError);
      }
    }

    //modules
    class ModuleVisitor : public GeometryVisitor {
    private:
      WeightDistributionGrid& weightDistribution_;
    public:
      ModuleVisitor(WeightDistributionGrid& weightDistribution) :
        weightDistribution_(weightDistribution) {}
      virtual ~ModuleVisitor() {}

      void visit(DetectorModule& module) {
        //ModuleCap* moduleCap = module.getModuleCap();
        //MaterialProperties* materialProperties = ModuleCap;
        //module.materialObject().populateMaterialProperties(*materialProperties);
        module.materialObject().populateMaterialProperties(*module.getModuleCap());

        //weightDistribution_.addTotalGrams(module.minZ(), module.minR(), module.maxZ(), module.maxR(), module.length(), module.area(), module.materialObject());
      }
    };

    ModuleVisitor visitor(weightDistribution);
    tracker.accept(visitor);
  }

  /*
  void Materialway::calculateMaterialValues(Tracker& tracker) {
    //sections
    for(Section* section : sectionsList_) {
      if(section->inactiveElement() != nullptr) {
        section->inactiveElement()->calculateTotalMass();
        section->inactiveElement()->calculateRadiationLength();
        section->inactiveElement()->calculateInteractionLength();
      }
    }

    //modules
    class ModuleVisitor : public GeometryVisitor {
    public:
      ModuleVisitor() {}
      virtual ~ModuleVisitor() {}

      void visit(DetectorModule& module) {
        ModuleCap* moduleCap = module.getModuleCap();
        moduleCap->calculateTotalMass();
        moduleCap->calculateRadiationLength();
        moduleCap->calculateInteractionLength();
      }
    };

    ModuleVisitor visitor;
    tracker.accept(visitor);
  }
  */

  void Materialway::buildInactiveSurface(Tracker& tracker, InactiveSurfaces& inactiveSurface) {
    class SupportVisitor : public GeometryVisitor {
    private:
      InactiveSurfaces& inactiveSurface_;
    public:
      SupportVisitor(InactiveSurfaces& inactiveSurface) : inactiveSurface_(inactiveSurface) {}
    
      void visit (Tracker& tracker) {
        for (auto& supportStructure : tracker.supportStructures()) {
          supportStructure.updateInactiveSurfaces(inactiveSurface_);
        }
      }
      void visit (Barrel& barrel) {
        for (auto& supportStructure : barrel.supportStructures()) {
          supportStructure.updateInactiveSurfaces(inactiveSurface_);
        }
      }
      void visit (Endcap& endcap) {
        for (auto& supportStructure : endcap.supportStructures()) {
          supportStructure.updateInactiveSurfaces(inactiveSurface_);
        }
      }
    };

    SupportVisitor supportVisitor(inactiveSurface);
    tracker.accept(supportVisitor);
    
    
    for(Section* section : sectionsList_) {
      if(section->inactiveElement() != nullptr) {
        /*
        if(section->inactiveElement()->getInteractionLength() < 0){
          std::cout<<"ERRORE INT LEN"<<std::endl;
        }
        */
	inactiveSurface.addBarrelServicePart(*section->inactiveElement());
        if((section->inactiveElement()->localMassCount() == 0)) {
          logWARNING(std::string(Form("Empty inactive element at r=%f, dr=%f, z=%f, dz=%f",
				      section->inactiveElement()->getInnerRadius(),
				      section->inactiveElement()->getRWidth(),
				      section->inactiveElement()->getZOffset(),
				      section->inactiveElement()->getZLength() )));
        }
      } else {
        logUniqueERROR(inactiveElementError);
      }
    }
    /*
    std::vector<InactiveElement>& elements = inactiveSurface.getBarrelServices();
    for (InactiveElement& currElem : elements) {
      if (currElem.getInteractionLength() < 0){
        std::cout<<"ERRORE INT LEN"<<std::endl;
      }
    }
    */
  }

  void Materialway::calculateMaterialValues(InactiveSurfaces& inactiveSurface, Tracker& tracker) {
    //supports
    for (InactiveElement& currElem : inactiveSurface.getSupports()) {
      currElem.calculateTotalMass();
      currElem.calculateRadiationLength();
      currElem.calculateInteractionLength();
    }

    //sections
    for (InactiveElement& currElem : inactiveSurface.getBarrelServices()) {
      currElem.calculateTotalMass();
      currElem.calculateRadiationLength();
      currElem.calculateInteractionLength();
    }

    //modules
    class ModuleVisitor : public GeometryVisitor {
    public:
      ModuleVisitor() {}
      virtual ~ModuleVisitor() {}

      void visit(DetectorModule& module) {
        ModuleCap* moduleCap = module.getModuleCap();
        moduleCap->calculateTotalMass();
        moduleCap->calculateRadiationLength();
        moduleCap->calculateInteractionLength();
      }
    };

    ModuleVisitor visitor;
    tracker.accept(visitor);
  }

  /*
  InactiveElement* Materialway::buildOppositeInactiveElement(InactiveElement* inactiveElement) {
    InactiveElement* newInactiveElement = new InactiveElement();
    newInactiveElement->setVertical(inactiveElement->isVertical());
    newInactiveElement->setFinal(inactiveElement->isFinal());
    newInactiveElement->setZOffset(-1 * inactiveElement->getZOffset());
    newInactiveElement->setZLength(inactiveElement->getZLength());
    newInactiveElement->setInnerRadius(inactiveElement->getInnerRadius());
    newInactiveElement->setRWidth(inactiveElement->getRWidth());
    newInactiveElement->setTotalMass(inactiveElement->getTotalMass());
    newInactiveElement->setLocalMass(inactiveElement->getLocalMass());

    return newInactiveElement;
  }
  */

  //END Materialway

} /* namespace material */


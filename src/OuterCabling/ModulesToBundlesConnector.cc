#include "OuterCabling/ModulesToBundlesConnector.hh"
#include <Tracker.hh>

void ModulesToBundlesConnector::visit(Barrel &b) {
  isBarrel_ = true;
  barrelName_ = b.myid();
}

void ModulesToBundlesConnector::visit(Layer &l) {
  layerNumber_ = l.layerNumber();
  numRods_ = l.numRods();
  totalNumFlatRings_ = l.buildNumModulesFlat() * 2 - 1;
}

void ModulesToBundlesConnector::visit(RodPair &r) {
  bundleType_ = computeBundleType(isBarrel_, barrelName_,
                                  layerNumber_); // Bundle cabling type

  rodPhi_ = r.Phi();
}

void ModulesToBundlesConnector::visit(BarrelModule &m) {
  side_ = (m.uniRef().side > 0.); // geometrical Z-side

  bool isPositiveCablingSide = side_; // By default.
  bool isTilted = false;              // Is the layer tilted ?
  bool isExtraFlatPart =
      false; // Used for Layer 3 flat part only, to add an extra bundle.

  // Here the idea is to compute the cabling side, the isTilted and
  // isExtraFlatPart booleans. 2S BUNDLES
  if (barrelName_ == outer_cabling_tb2s) {

    isPositiveCablingSide = side_;
    isTilted = false;
    isExtraFlatPart = false;
  }

  // PS BUNDLES
  else if (barrelName_ == outer_cabling_tbps) {

    // FLAT PART
    if (!m.isTilted()) {
      double phiSegmentWidth = (2. * M_PI) / numRods_;
      // This is the case where the cabling side can be different from the
      // geometrical side. The full flat rod (both -Z and +Z modules) is
      // assigned to one cabling side only.
      isPositiveCablingSide =
          computeBarrelFlatPartRodCablingSide(rodPhi_, phiSegmentWidth);
      isTilted = false;
      isExtraFlatPart = false;

      // For layer 3, need to add a second bundle for flat part
      if (isPositiveCablingSide &&
          (totalNumFlatRings_ > outer_cabling_maxNumModulesPerBundle) && !side_)
        isExtraFlatPart = true;
      if (!isPositiveCablingSide &&
          (totalNumFlatRings_ > outer_cabling_maxNumModulesPerBundle) && side_)
        isExtraFlatPart = true;
    }

    // TILTED PART
    else if (m.isTilted()) {

      isPositiveCablingSide = side_;
      isTilted = true;
      isExtraFlatPart = false;
    }
  }

  // NOW THAT ALL INFORMATION HAS BEEN GATHERED, COMPUTE PHIPOSITION.
  const PhiPosition &modulePhiPosition =
      PhiPosition(rodPhi_, numRods_, isBarrel_, layerNumber_);

  // BUILD BUNDLE IF NECESSARY, AND CONNECT MODULE TO BUNDLE
  buildBundle(m, bundles_, negBundles_, bundleType_, isBarrel_, barrelName_,
              layerNumber_, modulePhiPosition, isPositiveCablingSide,
              totalNumFlatRings_, isTilted, isExtraFlatPart);
}

void ModulesToBundlesConnector::visit(Endcap &e) {
  isBarrel_ = false;
  endcapName_ = e.myid();
}

void ModulesToBundlesConnector::visit(Disk &d) {
  diskNumber_ = d.diskNumber();
  side_ = d.side(); // geometrical Z-side
}

void ModulesToBundlesConnector::visit(Ring &r) {
  ringNumber_ = r.myid();
  numModulesInRing_ = r.numModules();

  bundleType_ =
      computeBundleType(isBarrel_, endcapName_, diskNumber_, ringNumber_);
}

void ModulesToBundlesConnector::visit(EndcapModule &m) {
  double modPhi = m.center().Phi();

  bool isPositiveCablingSide = side_; // Alyways true in the Endcaps : cabling
                                      // side and geometrical side are the same.

  // NOW THAT ALL INFORMATION HAS BEEN GATHERED, COMPUTE PHIPOSITION.
  const PhiPosition &modulePhiPosition =
      PhiPosition(modPhi, numModulesInRing_, isBarrel_, diskNumber_,
                  endcapName_, bundleType_);

  // BUILD BUNDLE IF NECESSARY, AND CONNECT MODULE TO BUNDLE
  buildBundle(m, bundles_, negBundles_, bundleType_, isBarrel_, endcapName_,
              diskNumber_, modulePhiPosition, isPositiveCablingSide);
}

void ModulesToBundlesConnector::postVisit() {
  // STAGGER MODULES
  staggerModules(bundles_);
  staggerModules(negBundles_);

  // CHECK
  checkModulesToBundlesCabling(bundles_);
  checkModulesToBundlesCabling(negBundles_);

  // ENDCAPS ONLY: ASSIGN THE MFB FANOUTS BRANCHES TO THE MODULES
  connectEndcapModulesToBundlesFanoutBranches(bundles_);
  connectEndcapModulesToBundlesFanoutBranches(negBundles_);
}

/* Compute the cabling side, for flat parts of Barrel rods.
   The full flat rod (both -Z and +Z modules) is assigned to one cabling side
   only. Hence here, the cabling side can be different from the actual
   geometrical side. Modules are alternatively connected to the positive cabling
   side or the negative cabling side, depending on the Phi of the rod.
*/
const bool ModulesToBundlesConnector::computeBarrelFlatPartRodCablingSide(
    const double rodPhi, const double phiSegmentWidth) const {
  // Compute the phi position of the rod (in one cabling side frame of
  // reference).
  const double phiSegmentStartOneCablingSide =
      computePhiSegmentStart(rodPhi, phiSegmentWidth);
  const int phiSegmentRefOneCablingSide = computePhiSegmentRef(
      rodPhi, phiSegmentStartOneCablingSide, phiSegmentWidth);
  // Assign the full rod to the positive cabling side or the negative cabling
  // side alternatively.
  const bool isPositiveCablingSide = ((phiSegmentRefOneCablingSide % 2) == 1);
  return isPositiveCablingSide;
}

/* Compute the bundle cabling type.
 */
const Category ModulesToBundlesConnector::computeBundleType(
    const bool isBarrel, const std::string subDetectorName,
    const int layerDiskNumber, const int ringNumber) const {
  Category bundleType = Category::UNDEFINED;

  // BARREL
  if (isBarrel) {
    // TB2S
    if (subDetectorName == outer_cabling_tb2s) {
      bundleType = Category::SS;
    }
    // TBPS
    else if (subDetectorName == outer_cabling_tbps) {
      bundleType = (layerDiskNumber == 3 ? Category::PS5G : Category::PS10G);
    }
  }

  // ENDCAPS
  else {
    // TEDD_1
    if (subDetectorName == outer_cabling_tedd1) {
      if (ringNumber <= 4)
        bundleType = Category::PS10GA;
      else if (ringNumber >= 5 && ringNumber <= 7)
        bundleType = Category::PS10GB;
      else if (ringNumber >= 8 && ringNumber <= 10)
        bundleType = Category::PS5G;
      else if (ringNumber >= 11)
        bundleType = Category::SS;
    }

    // TEDD_2
    else if (subDetectorName == outer_cabling_tedd2) {
      if (ringNumber <= 3)
        bundleType = Category::UNDEFINED;
      else if (ringNumber >= 4 && ringNumber <= 6)
        bundleType = Category::PS10GB;
      else if (ringNumber >= 7 && ringNumber <= 10)
        bundleType = Category::PS5G;
      else if (ringNumber >= 11)
        bundleType = Category::SS;
    }
  }

  return bundleType;
}

/* Build bundle.
 * The index associated to the bundleType is computed.
 * phiSliceRef: phiSegmentRef in Barrel, phiRegionRef in Endcaps.
 * This allows to compute the bundle Id.
 * Then, the bundle is created, and stored in the bundles_ or negBundles_
 * containers. Lastly, each module is connected to its bundle, and vice-versa.
 */
void ModulesToBundlesConnector::buildBundle(
    DetectorModule &m, std::map<int, OuterBundle *> &bundles,
    std::map<int, OuterBundle *> &negBundles, const Category &bundleType,
    const bool isBarrel, const std::string subDetectorName,
    const int layerDiskNumber, const PhiPosition &modulePhiPosition,
    const bool isPositiveCablingSide, const int totalNumFlatRings,
    const bool isTiltedPart, const bool isExtraFlatPart) {

  // COMPUTE BUNDLE ID
  const int bundleTypeIndex = computeBundleTypeIndex(
      isBarrel, bundleType, totalNumFlatRings, isTiltedPart, isExtraFlatPart);
  const int phiSliceRef = (isBarrel ? modulePhiPosition.phiSegmentRef()
                                    : modulePhiPosition.phiRegionRef());
  const int bundleId =
      computeBundleId(isBarrel, isPositiveCablingSide, layerDiskNumber,
                      phiSliceRef, bundleTypeIndex);

  // COMPUTE STEREO BUNDLE ID (BARREL ONLY, NOT NEEDED FOR THE ENDCAPS SO FAR)
  // The stereo bundle is the bundle located on the other cabling side, by
  // rotation of 180° around CMS_Y.
  const int stereoPhiSliceRef =
      (isBarrel ? modulePhiPosition.stereoPhiSegmentRef()
                : 0); // modulePhiPosition.stereoPhiRegionRef());
  const int stereoBundleId =
      computeStereoBundleId(isBarrel, isPositiveCablingSide, layerDiskNumber,
                            stereoPhiSliceRef, bundleTypeIndex);

  // All Bundles from one cabling side
  const std::map<int, OuterBundle *> &bundlesOneSide =
      (isPositiveCablingSide ? bundles : negBundles);

  // CREATE BUNDLE IF NECESSARY
  OuterBundle *bundle = nullptr;
  auto found = bundlesOneSide.find(bundleId);
  if (found == bundlesOneSide.end()) {
    bundle = createAndStoreBundle(bundles, negBundles, bundleId, stereoBundleId,
                                  bundleType, subDetectorName, layerDiskNumber,
                                  modulePhiPosition, isPositiveCablingSide,
                                  isTiltedPart);
  } else {
    bundle = found->second;
  }

  // CONNECT MODULE TO BUNDLE
  connectModuleToBundle(m, bundle);
}

/* Compute the index associated to each bundle type.
 */
const int ModulesToBundlesConnector::computeBundleTypeIndex(
    const bool isBarrel, const Category &bundleType,
    const int totalNumFlatRings, const bool isTilted,
    const bool isExtraFlatPart) const {
  int bundleTypeIndex;
  // BARREL
  if (isBarrel) {
    if (bundleType == Category::SS)
      bundleTypeIndex = 0;
    else {
      if (!isTilted) {
        if (totalNumFlatRings <= outer_cabling_maxNumModulesPerBundle)
          bundleTypeIndex = 1;
        else {
          if (!isExtraFlatPart)
            bundleTypeIndex = 1;
          else
            bundleTypeIndex = 2;
        }
      } else
        bundleTypeIndex = 0;
    }
  }
  // ENDCAPS
  else {
    if (bundleType == Category::PS10GA)
      bundleTypeIndex = 0;
    else if (bundleType == Category::PS10GB)
      bundleTypeIndex = 1;
    else if (bundleType == Category::PS5G)
      bundleTypeIndex = 2;
    else if (bundleType == Category::SS)
      bundleTypeIndex = 3;
  }
  return bundleTypeIndex;
}

/* Compute the Id associated to each bundle.
 */
const int ModulesToBundlesConnector::computeBundleId(
    const bool isBarrel, const bool isPositiveCablingSide,
    const int layerDiskNumber, const int phiRef,
    const int bundleTypeIndex) const {
  int cablingSideIndex = 0;
  if (isBarrel) {
    cablingSideIndex = (isPositiveCablingSide ? 1 : 3);
  } else {
    cablingSideIndex = (isPositiveCablingSide ? 2 : 4);
  }

  const int bundleId = cablingSideIndex * 10000 + layerDiskNumber * 1000 +
                       phiRef * 10 + bundleTypeIndex;
  return bundleId;
}

/* Compute the Id associated to each stereo bundle.
 * The stereo bundle is the bundle located on the other cabling side, by
 * rotation of 180° around CMS_Y.
 */
const int ModulesToBundlesConnector::computeStereoBundleId(
    const bool isBarrel, const bool isPositiveCablingSide,
    const int layerDiskNumber, const int stereoPhiRef,
    const int bundleTypeIndex) const {
  int stereoBundleId = 0;
  if (isBarrel) {
    const bool stereoSide = !isPositiveCablingSide;
    stereoBundleId = computeBundleId(isBarrel, stereoSide, layerDiskNumber,
                                     stereoPhiRef, bundleTypeIndex);
  }

  return stereoBundleId;
}

/* Create a bundle, if it does not exist yet.
 *  Store it in the bundles_ or negBundles_ containers.
 */
OuterBundle *ModulesToBundlesConnector::createAndStoreBundle(
    std::map<int, OuterBundle *> &bundles,
    std::map<int, OuterBundle *> &negBundles, const int bundleId,
    const int stereoBundleId, const Category &bundleType,
    const std::string subDetectorName, const int layerDiskNumber,
    const PhiPosition &modulePhiPosition, const bool isPositiveCablingSide,
    const bool isTiltedPart) {

  OuterBundle *bundle = new OuterBundle(
      bundleId, stereoBundleId, bundleType, subDetectorName, layerDiskNumber,
      modulePhiPosition, isPositiveCablingSide, isTiltedPart);

  if (isPositiveCablingSide) {
    bundles.insert(std::make_pair(bundleId, bundle));
  } else {
    negBundles.insert(std::make_pair(bundleId, bundle));
  }
  return bundle;
}

/* Connect module to bundle and vice-versa.
 */
void ModulesToBundlesConnector::connectModuleToBundle(
    DetectorModule &m, OuterBundle *bundle) const {
  bundle->addModule(&m);
  m.setBundle(bundle);
}

/* Stagger modules.
   This is a very important step.
   In theory, each module is connected to the bundle corresponding to its
   position in Phi. Though, there can be more modules connected to a bundle than
   possible. If this is the case, one needs to stagger modules : the module
   closest to an adjacent phiRegion is removed and placed in the adjacent
   phiRegion.
*/
void ModulesToBundlesConnector::staggerModules(
    std::map<int, OuterBundle *> &bundles) {

  for (auto &b : bundles) {
    // All this happens in Endcaps only
    if (b.second->subDetectorName() == outer_cabling_tedd1 ||
        b.second->subDetectorName() == outer_cabling_tedd2) {
      const bool isBarrel = false;

      // Too many modules per bundle: staggering needed!
      while (b.second->numModules() > outer_cabling_maxNumModulesPerBundle) {
        const bool isPositiveCablingSide = b.second->isPositiveCablingSide();
        const int diskNumber = b.second->layerDiskNumber();

        const Category &bundleType = b.second->type();
        const int bundleTypeIndex =
            computeBundleTypeIndex(isBarrel, bundleType);

        const PhiPosition &bundlePhiPosition = b.second->phiPosition();

        const double phiRegionStart = bundlePhiPosition.phiRegionStart();

        const double phiRegionWidth = bundlePhiPosition.phiRegionWidth();
        const int numPhiRegions = round(2 * M_PI / phiRegionWidth);

        // Compute the ref of the phiRegion, and of the adjacent phi regions.
        // 'next' or 'previous' refers to:
        // - positiveCablingSide: positive Phi order.
        // - negativeCablingSide: negative Phi order.
        const int phiRegionRef = bundlePhiPosition.phiRegionRef();
        const int nextPhiRegionRef =
            computeNextPhiSliceRef(phiRegionRef, numPhiRegions);
        const int previousPhiRegionRef =
            computePreviousPhiSliceRef(phiRegionRef, numPhiRegions);

        // Compute the associated bundles ids (so that the associated bundles
        // can be accessed).
        // const int bundleId = b.first;
        const int nextBundleId =
            computeBundleId(isBarrel, isPositiveCablingSide, diskNumber,
                            nextPhiRegionRef, bundleTypeIndex);
        const int previousBundleId =
            computeBundleId(isBarrel, isPositiveCablingSide, diskNumber,
                            previousPhiRegionRef, bundleTypeIndex);

        // Distance in Phi from the smallest Phi module to the smallest Phi
        // boundary of the phiRegion.
        const double minPhiBorder =
            fabs(femod((b.second->minPhi() - phiRegionStart), phiRegionWidth));
        // Distance in Phi from the greatest Phi module to the greatest Phi
        // boundary of the phiRegion.
        const double maxPhiBorder =
            fabs(femod((b.second->maxPhi() - phiRegionStart), phiRegionWidth) -
                 phiRegionWidth);

        // Access bundles.
        auto previousBundleSearch = bundles.find(previousBundleId);
        auto nextBundleSearch = bundles.find(nextBundleId);
        if (previousBundleSearch != bundles.end() &&
            nextBundleSearch != bundles.end()) {
          OuterBundle *previousBundle = previousBundleSearch->second;
          OuterBundle *nextBundle = nextBundleSearch->second;

          const int previousBundleNumModules = previousBundle->numModules();
          const int nextBundleNumModules = nextBundle->numModules();
          // A MODULE FROM THE PHI REGION NEEDS TO BE MOVED. LOOK WHERE TO MOVE
          // IT!

          // Cannot assign the extra module : both neighbouring phi regions are
          // full !
          if (previousBundleNumModules >=
                  outer_cabling_maxNumModulesPerBundle &&
              nextBundleNumModules >= outer_cabling_maxNumModulesPerBundle) {
            logERROR(
                any2str("Building cabling map : Staggering modules.") +
                "I am a module in side " + any2str(isPositiveCablingSide) +
                ", disk " + any2str(diskNumber) + ", bundleType " +
                any2str(bundleType) + ", phiRegionRef " +
                any2str(phiRegionRef) + ", phiRegionWidth " +
                any2str(phiRegionWidth) + ". My phiRegion has more than " +
                any2str(outer_cabling_maxNumModulesPerBundle) + " modules." +
                " My 2 neighbouring phiRegions have also more than " +
                any2str(outer_cabling_maxNumModulesPerBundle) + " modules." +
                " Hence, I cannot be staggered to any other phiRegion :/ ");
            break;
          }

          // Assign a module to the next phi region.
          else if (previousBundleNumModules >=
                       outer_cabling_maxNumModulesPerBundle ||
                   maxPhiBorder <= minPhiBorder) {
            logINFO(any2str("Building cabling map : Staggering modules.") +
                    " I am a module in side " + any2str(isPositiveCablingSide) +
                    ", disk " + any2str(diskNumber) + ", bundleType " +
                    any2str(bundleType) + ", phiRegionRef " +
                    any2str(phiRegionRef) + ". maxPhiBorder " +
                    any2str((maxPhiBorder * 180. / M_PI)) + ". My region has " +
                    any2str(b.second->numModules()) +
                    " > maxNumModulesPerBundle = " +
                    any2str(outer_cabling_maxNumModulesPerBundle) +
                    ". I am moved to the next phiRegion, which presently has " +
                    any2str(nextBundleNumModules) + " modules.");

            const auto maxPhiModIt = b.second->maxPhiModule();
            if (maxPhiModIt != b.second->modules().end()) {
              Module *maxPhiMod = *maxPhiModIt;
              maxPhiMod->setBundle(nextBundle);
              nextBundle->moveMaxPhiModuleFromOtherBundle(b.second);
            } else { // Modules not found when trying to access them.
              logERROR(
                  any2str("Building cabling map : Staggering modules.") +
                  "Bundle has more than one module, but can find any module!");
              break;
            }
          }

          // Assign a module to the previous phi region.
          else if (nextBundleNumModules >=
                       outer_cabling_maxNumModulesPerBundle ||
                   minPhiBorder < maxPhiBorder) {
            logINFO(
                any2str("Building cabling map : Staggering modules.") +
                " I am a module in side " + any2str(isPositiveCablingSide) +
                ", disk " + any2str(diskNumber) + ", bundleType " +
                any2str(bundleType) + ", phiRegionRef " +
                any2str(phiRegionRef) + ". minPhiBorder " +
                any2str((minPhiBorder * 180. / M_PI)) + ". My region has " +
                any2str(b.second->numModules()) +
                " > maxNumModulesPerBundle = " +
                any2str(outer_cabling_maxNumModulesPerBundle) +
                ". I am moved to the previous phiRegion, which presently has " +
                any2str(previousBundleNumModules) + " modules.");

            const auto minPhiModIt = b.second->minPhiModule();
            if (minPhiModIt != b.second->modules().end()) {
              Module *minPhiMod = *minPhiModIt;
              minPhiMod->setBundle(previousBundle);
              previousBundle->moveMinPhiModuleFromOtherBundle(b.second);
            } else { // Modules not found when trying to access them.
              logERROR(
                  any2str("Building cabling map : Staggering modules.") +
                  "Bundle has more than one module, but can find any module!");
              break;
            }
          }
        } else { // Bundles not found when trying to access them.
          logERROR(any2str("Building cabling map : Staggering modules.") +
                   "Error building previousBundleId or nextBundleId");
          break;
        }
      }
    }
  }
}

/* Check modules-bundles connections.
 */
void ModulesToBundlesConnector::checkModulesToBundlesCabling(
    const std::map<int, OuterBundle *> &bundles) const {
  for (auto &b : bundles) {

    // CHECK WHETHER THE PHI SLICES REF MAKE SENSE.
    const PhiPosition &bundlePhiPosition = b.second->phiPosition();
    const int phiSegmentRef = bundlePhiPosition.phiSegmentRef();
    const int phiRegionRef = bundlePhiPosition.phiRegionRef();
    const int phiSectorRef = bundlePhiPosition.phiSectorRef();
    if (phiSegmentRef <= -1 || phiRegionRef <= -1 || phiSectorRef <= -1 ||
        phiSectorRef >= outer_cabling_numNonants) {
      logERROR(
          any2str(
              "Building cabling map : a bundle was not correctly created. ") +
          "OuterBundle " + any2str(b.first) +
          ", with bundleType = " + any2str(b.second->type()) +
          ", has phiSegmentRef = " + any2str(phiSegmentRef) +
          ", phiRegionRef = " + any2str(phiRegionRef) +
          ", phiSectorRef = " + any2str(phiSectorRef));
    }

    // CHECK THE NUMBER OF MODULES PER BUNDLE.
    const int bundleNumModules = b.second->numModules();
    if (bundleNumModules > outer_cabling_maxNumModulesPerBundle) {
      logERROR(any2str("Building cabling map : Staggering modules. ") +
               "OuterBundle " + any2str(b.first) + " is connected to " +
               any2str(bundleNumModules) + " modules." +
               " Maximum number of modules per bundle allowed is " +
               any2str(outer_cabling_maxNumModulesPerBundle));
    }
  }
}

/*
  Endcaps only: assign the MFB fanouts branches to the modules.
  The modules are connected to the MFBs through 4 fanouts branches.
  In the Endcaps, this ends up complicated. It was hence requested by the
  Mechanics to have it in at automatic way. This work can obviously only be done
  once the modules are assigned to the MFBs in a final way (hence once the
  modules staggering is already done).

  Concept is the following: for each MFB, one looks at the modules connected to
  it. Standard case: 1 fanout branch per disk surface. Special case: phi sector
  crossing (Y=0): if a MFB has several modules on the same disk (not
  double-disk!) and in different dees, then the fanout branch assignment is per
  dee.

  Assumptions: 4 disk surfaces per dee + dee boundary at (Y=0).
 */
void ModulesToBundlesConnector::connectEndcapModulesToBundlesFanoutBranches(
    std::map<int, OuterBundle *> &bundles) {
  // Loop on all MFBs
  for (auto &b : bundles) {
    // WORK ON A PER MFB BASIS
    const OuterBundle *myBundle = b.second;
    // Endcaps only
    if (myBundle->subDetectorName() == outer_cabling_tedd1 ||
        myBundle->subDetectorName() == outer_cabling_tedd2) {

      // STEP A: FIND, FOR EACH DISK (NOT DOUBLE-DISK!!), HOW THE BUNDLE'S
      // MODULES SHOULD BE SORTED.

      // BY DEFAULT, SORTING IS PER DISK SURFACE.
      bool sortModulesInSmallAbsZDiskPerDee = false;
      bool sortModulesInBigAbsZDiskPerDee = false;

      // Assumes 4 disk surfaces per double-disk
      const int numDiskSurfacesPerDoubleDisk = 4;

      // Get all the modules connected to the MFB
      const std::vector<Module *> &myModules = myBundle->modules();

      // IF THE PHI SECTOR IS CROSSING SEVERAL DEES:
      // WE NEED TO LOOK ON WHICH DEE MODULES ARE ON.
      const int &myPhiSectorRef = myBundle->phiPosition().phiSectorRef();
      if (myPhiSectorRef == outer_cabling_phiNonantCrossingDifferentDeesIndex) {

        // Store, per disk surface, which dee modules are on.
        std::map<int, bool> hasDiskSurfaceAnyModuleInYPosDee;
        std::map<int, bool> hasDiskSurfaceAnyModuleInYNegDee;
        // Initialize: by default, there is no module.
        for (int diskSurfaceIndex = 1;
             diskSurfaceIndex <= numDiskSurfacesPerDoubleDisk;
             diskSurfaceIndex++) {
          hasDiskSurfaceAnyModuleInYPosDee[diskSurfaceIndex] = false;
          hasDiskSurfaceAnyModuleInYNegDee[diskSurfaceIndex] = false;
        }

        // loop on all modules connected to a MFB
        for (const auto &module : myModules) {
          const double modPhiModuloPi =
              femod(module->center().Phi(), M_PI); // modulo PI

          // STORE, PER DISK SURFACE, WHICH DEE MODULES ARE ON.
          if (fabs(modPhiModuloPi) > outer_cabling_roundingTolerance &&
              fabs(modPhiModuloPi - M_PI) >
                  outer_cabling_roundingTolerance // module at (Y=0), if any,
                                                  // does not matter!
          ) {
            const int diskSurfaceIndex = module->diskSurface();
            const double modPhi =
                femod(module->center().Phi(), 2. * M_PI); // modulo 2 PI

            // module is on (Y > 0) dee
            if (modPhi < M_PI) {
              hasDiskSurfaceAnyModuleInYPosDee[diskSurfaceIndex] = true;
            }
            // module is on (Y < 0) dee
            else {
              hasDiskSurfaceAnyModuleInYNegDee[diskSurfaceIndex] = true;
            }
          }
        }

        // IF MODULES ARE ON 2 DEES OF THE SAME DISK SURFACE: CHOOSE DEE
        // SORTING!
        for (int diskSurfaceIndex = 1;
             diskSurfaceIndex <= numDiskSurfacesPerDoubleDisk;
             diskSurfaceIndex++) {
          if (hasDiskSurfaceAnyModuleInYPosDee.at(diskSurfaceIndex) &&
              hasDiskSurfaceAnyModuleInYNegDee.at(diskSurfaceIndex)) {
            // If surface 1 or 2 has modules on both dees
            if (diskSurfaceIndex == 1 || diskSurfaceIndex == 2) {
              // then small |Z| disk should have a sorting per dee!
              sortModulesInSmallAbsZDiskPerDee = true;
            }
            // If surface 3 or 4 has modules on both dees
            else if (diskSurfaceIndex == 3 || diskSurfaceIndex == 4) {
              // then big |Z| disk should have a sorting per dee!
              sortModulesInBigAbsZDiskPerDee = true;
            } else {
              logERROR(any2str(
                  "Unexpected number of disk surfaces per double-disk"));
            }
          }
        }
      }

      // STEP B: ASSIGN THE FANOUT BRANCH INDEX TO THE MFB'S MODULES.
      for (const auto &module : myModules) {
        const int diskSurfaceIndex = module->diskSurface();

        // Find out which sorting mode is relevant
        const bool sortModulesInDiskPerDee =
            ((diskSurfaceIndex == 1 || (diskSurfaceIndex == 2))
                 ? sortModulesInSmallAbsZDiskPerDee
                 : sortModulesInBigAbsZDiskPerDee);
        // Sort modules per disk surface:
        // (+Z) end: fanout branch index = module's disk surface index.
        // (-Z) end: fanout branch index = module's disk surface complementary
        // index. The index itself does not matter anyway (only the grouping
        // does). This indexing scheme, though, allows to show the cabling
        // similitude on both (Z) ends:
        // ((+Z) end, disk surface) <-> ((-Z) end, complementary disk surface)
        // have the exact same cabling. Complementary disk surface: as defined
        // in the code below.
        if (!sortModulesInDiskPerDee) {
          const int isPositiveCablingSide = module->isPositiveCablingSide();
          int endcapFiberFanoutBranchIndex = 0;
          // (+Z) end: fanout branch index = module's disk surface index.
          if (isPositiveCablingSide > 0) {
            endcapFiberFanoutBranchIndex = diskSurfaceIndex;

          }
          // (-Z) end: fanout branch index = module's disk surface complementary
          // index.
          else {
            // This defines the complementary disk surface, to the disk surface
            // with diskSurfaceIndex:
            if (diskSurfaceIndex == 1) {
              endcapFiberFanoutBranchIndex = 2;
            } else if (diskSurfaceIndex == 2) {
              endcapFiberFanoutBranchIndex = 1;
            } else if (diskSurfaceIndex == 3) {
              endcapFiberFanoutBranchIndex = 4;
            } else if (diskSurfaceIndex == 4) {
              endcapFiberFanoutBranchIndex = 3;
            } else {
              logERROR(any2str(
                  "Unexpected number of disk surfaces per double-disk"));
            }
          }
          module->setEndcapFiberFanoutBranch(endcapFiberFanoutBranchIndex);
        }
        // Sort modules per dee:
        // Small |Z| disk, (Y > 0) dee: fanout branch index = 1
        // Small |Z| disk, (Y < 0) dee: fanout branch index = 2
        // Big |Z| disk, (Y > 0) dee: fanout branch index = 3
        // Big |Z| disk, (Y < 0) dee: fanout branch index = 4
        else {
          const double modPhi = femod(module->center().Phi(), 2. * M_PI);
          int endcapFiberFanoutBranchIndex = ((modPhi < M_PI) ? 1 : 2);
          if (diskSurfaceIndex >= 3)
            endcapFiberFanoutBranchIndex += 2;
          module->setEndcapFiberFanoutBranch(endcapFiberFanoutBranchIndex);
        }
      }
    }
  }

  checkEndcapModulesToBundlesFanoutBranchesCabling(bundles);
}

/* TEDD: Check modules - bundles fanouts branches connections.
 */
void ModulesToBundlesConnector::
    checkEndcapModulesToBundlesFanoutBranchesCabling(
        const std::map<int, OuterBundle *> &bundles) const {

  // Loop on all MFBs
  for (auto &b : bundles) {
    // Work on a per MFB basis
    const OuterBundle *myBundle = b.second;
    // Endcaps only
    if (myBundle->subDetectorName() == outer_cabling_tedd1 ||
        myBundle->subDetectorName() == outer_cabling_tedd2) {
      std::set<int> allBranchesIndexes;

      // Loop on all the modules connected to the MFB, ancd check their fanout
      // branch index.
      const std::vector<Module *> &myModules = myBundle->modules();
      for (const auto &module : myModules) {
        const int endcapFiberFanoutBranchIndex =
            module->getEndcapFiberFanoutBranch();
        if (endcapFiberFanoutBranchIndex <= 0 ||
            endcapFiberFanoutBranchIndex >
                outer_cabling_maxNumFanoutBranchesPerEndcapBundle) {
          logERROR(any2str("TEDD: assigning MFB fanouts to modules. ") +
                   "OuterBundle " + any2str(b.first) +
                   ". Found a module connected to fanout index = " +
                   any2str(endcapFiberFanoutBranchIndex) +
                   ", which is not supported.");
        }
        allBranchesIndexes.insert(endcapFiberFanoutBranchIndex);
      }

      // Check that, for each MFB, all fanout branches indexes are met at least
      // once.
      for (int branchIndex = 1;
           branchIndex <= outer_cabling_maxNumFanoutBranchesPerEndcapBundle;
           branchIndex++) {
        if (allBranchesIndexes.find(branchIndex) == allBranchesIndexes.end()) {
          logERROR(any2str("TEDD: assigning MFB fanouts to modules. ") +
                   "OuterBundle " + any2str(b.first) +
                   ". Found 0 module with fanout branch index = " +
                   any2str(branchIndex));
        }
      }
    }
  }
}

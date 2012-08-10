#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>

// My descriptor
#include "tracker.hh"
#include "global_funcs.h"

// ROOT objects
#include "TCanvas.h"
#include "TGeoManager.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TFile.h"
#include "TPolyLine3D.h"
#include "TText.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TStyle.h"

// Stuff to create directories
#include <sys/stat.h>
#include <sys/types.h>

// Date and time
#include <ctime> // for debug

// obsolete
double diffclock(clock_t clock1, clock_t clock2) {
  double diffticks=clock1-clock2;
  double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
  return diffms;
}

// comparators
bool smallerRho(Layer* l1, Layer* l2) { return l1->getMinRho() < l2->getMinRho(); }
bool smallerZ(Layer* l1, Layer* l2) { return l1->getMinZ() < l2->getMinZ(); }

// Endcap Module sorting for vectors
// Returns true if m1 is "lower" than m2
bool moduleSortEndcapStyle(const Module* m1, const Module* m2) {
  XYZVector position[2];
  int radius_mm[2];
  int phi_deg[2];

  // Put position in a comfortable array
  position[0] = m1->getMeanPoint();
  position[1] = m2->getMeanPoint();

  // First sort on radius (mm precision)
  for (int i=0; i<2; ++i) {
    radius_mm[i] = int(position[i].Rho());
  }
  if (radius_mm[0]<radius_mm[1]) return true;
  if (radius_mm[0]>radius_mm[1]) return false;

  // If radius iis equal, then sort on phi (degree precision)
  for (int i=0; i<2; ++i) {
    phi_deg[i] = int(position[i].Phi()/M_PI*180);
  }
  if (phi_deg[0]<phi_deg[1]) return true;
  if (phi_deg[0]>phi_deg[1]) return false;

  // Otherwise use Z (double precision) to sort
  return (position[0].Z()<position[1].Z());
}


using namespace ROOT::Math;

Tracker::~Tracker()  {
  LayerVector::iterator layIt;
  for (layIt=layerSet_.begin(); layIt!=layerSet_.end(); layIt++) {
    if ((*layIt)!=NULL) {
      delete (*layIt);
    }
  }
  layerSet_.clear();
}

Tracker::Tracker() {
  setDefaultParameters();
}

Tracker::Tracker(std::string trackerName) {
  setDefaultParameters();
  setName(trackerName);
}

void Tracker::setDefaultParameters() {
  nMB_ = defaultNMB_;
  bunchSpacingNs_ = defaultBunchSpacingNs_;
  rError_ = defaultRError_;
  zError_ = defaultZError_;
  useIPConstraint_ = defaultUseIPConstraint_;
  smallDelta_ = defaultSmallDelta_;
  bigDelta_ = defaultBigDelta_;
  overlap_ = defaultOverlap_;
  storeDirectory_ = "store";
  summaryDirectory_ = "summaries";
  trackerName_ = "aTracker";
  arguments_= "";
  maxL_ = 0;
  maxR_ = 0;
  phiSegments_ = 4;
  efficiency_ = 1;
  pixelEfficiency_ = 1;
  //colorPicker("ptOut"); // remove these three from here
  //colorPicker("rphi");
  //colorPicker("stereo");
  //colorPicker("ptIn");

  numInvFemtobarns_ = 3000;
  chargeDepletionVoltage_ = 600;
  operatingTemp_ = -20;
  alphaParam_ = 4e-17;
  referenceTemp_ = +20;

  triggerProcessorsPhi_ = 1;
  triggerProcessorsEta_ = 1;
  triggerEtaCut_        = 2.5;
  triggerPtCut_         = 1;

  std::string testMe = "test me!";
}

void Tracker::shapeVolume() {
  // TODO:
  // Should build only the containing volume
}

void Tracker::shapeLayerVolumes() {
  // TODO:
  // Will build one volume per layer
}


void Tracker::createGeometry(bool lite /*= false*/ ) {

  if (!lite) {

    // Define the geometry
    myGeom_ = new TGeoManager("trackerGeometry", "Tracker geometry");
    // Define materials and media
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial* matSi_ = new TGeoMaterial("Si", 26.98, 13, 2.7);
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum", 1, matVacuum);
    myMed_ = new TGeoMedium("Si", 2, matSi_);
    // Define and set the top volume
    myVolume_ = myGeom_->MakeBox("TOP", Vacuum, 270., 270., 120.);
    myGeom_->SetTopVolume(myVolume_);
    // Reset the module counter
    iModule_=0;
    // DO THE THING
    shapeModuleVolumes(false);
    myGeom_->CloseGeometry();
    savingV_.push_back(myGeom_);


  } else {

    // ***************************************
    // *                                     *
    // * Everything here for the line drawer *
    // *                                     *
    // ***************************************

    // TODO: if these objects were already created they must be
    // deleted togherther with their references in the savingV_

    geomLite_ = new TCanvas("geometryLite", "Modules geometry", 800, 800);
    geomLite_->cd();
    shapeModuleVolumes(true);
    savingV_.push_back(geomLite_);

    geomLiteXY_ = new TCanvas("geometryLiteXY", "Modules geometry (XY Section)", 800, 800);
    geomLiteXY_->cd();
    shapeModuleVolumes(true, Layer::XYSection);
    savingV_.push_back(geomLiteXY_);

    geomLiteYZ_ = new TCanvas("geometryLiteYZ", "Modules geometry (YZ Section)", 800, 800);
    geomLiteYZ_->cd();
    shapeModuleVolumes(true, Layer::YZSection|Layer::Forward);
    savingV_.push_back(geomLiteYZ_);

    geomLiteEC_ = new TCanvas("geometryLiteEC", "Modules geometry (Endcap)", 800, 800);
    geomLiteEC_->cd();
    shapeModuleVolumesEndcapSample(true);
    savingV_.push_back(geomLiteEC_);
  }
}

void Tracker::shapeModuleVolumesEndcapSample(bool lite /* = false */) {
  // TODO:
  // Will build all the modules volumes / contours
  ModuleVector::iterator modIt;

  for (modIt=endcapSample_.begin(); modIt!=endcapSample_.end(); modIt++) {
    if (!lite) {
      placeModule(*modIt);
    } else {
      placeModuleLite(*modIt);
    }
  }
}


void Tracker::shapeModuleVolumes(bool lite /* = false */, int section /* = Layer::NoSection*/ ) {
  // TODO:
  // Will build all the modules volumes / contours
  LayerVector::iterator layIt;
  ModuleVector::iterator modIt;
  ModuleVector* aLay;
  bool placeThis;
  int realSection=section & (~Layer::Forward);
  int thisRealSection;


  for (layIt=layerSet_.begin(); layIt!=layerSet_.end(); layIt++) {
    aLay = (*layIt)->getModuleVector();
    for (modIt=aLay->begin(); modIt!=aLay->end(); modIt++) {
      thisRealSection = (*modIt)->getSection() & (~Layer::Forward);
      if ((realSection&thisRealSection)==(realSection)) {
        placeThis=true;
        if ((((section)&(Layer::Forward))!=0) && ( (*modIt)->getMeanPoint().Z()<0  )) placeThis=false;
        if (placeThis) {
          if (!lite) {
            placeModule(*modIt);
          } else {
            placeModuleLite(*modIt);
          }
          //BarrelModule* bm;
          //if ((bm = dynamic_cast<BarrelModule*>(*modIt))) cout << "mod: " << bm->getContainerName() << bm->getLayer() << ", mp: " << bm->getMeanPoint().Z() 
          //                                                   << ", ez: " << bm->getEdgeZ(-1).first << ", " << bm->getEdgeZ(-1).second 
          //                                                   << ", ezs: " << bm->getEdgeZSide(-1).first << ", " << bm->getEdgeZSide(-1).second << endl;
        }
      }
    }
  }
}

void Tracker::placeModule(Module* aModule) {
  std::string moduleName;
  char moduleNumber[100];
  sprintf(moduleNumber, "%d", iModule_++);
  moduleName = "mod_" + std::string(moduleNumber);
  aModule->setId(moduleName);
  aModule->shapeVolume(myVolume_, myMed_, myGeom_);
}

void Tracker::placeModuleLite(Module* aModule) {
  std::string moduleName;
  char moduleNumber[100];
  sprintf(moduleNumber, "%d", iModule_++);
  moduleName = "mod_" + std::string(moduleNumber);
  aModule->setId(moduleName);
  TPolyLine3D* contour = aModule->getContour();
  contour->SetLineWidth(2);
  contour->Draw();

}


//void Tracker::buildBarrel(int nLayer,
//        double minRadius,
//        double maxRadius,
//        int nModules,
//        BarrelModule* sampleModule,
//        int section /* = NoSection */,
//       bool compressed /* = false */) {

//    buildBarrel(nLayer, minRadius, maxRadius, nModules, sampleModule, DEFAULTBARRELNAME, section, compressed);
//    
//}


// All the measures in mm as usual!
LayerVector Tracker::buildBarrel(int nLayer,
                                 double minRadius,
                                 double maxRadius,
                                 double maxZ, // maxZ and nModules should be alternative (i.e if one is set the other should be zero)
                                 int nModules,
                                 BarrelModule* sampleModule,
                                 std::string barrelName,
                                 int section /* = NoSection */,
                                 bool compressed /* = false */,
                                 bool shortBarrel /* = false, used to build mezzanine barrels */,
                                 bool sameRods /* = false, used to build all the layers in the barrel with identical (over-hermetic) rods */) {

  maxR_=(maxRadius>maxR_)?maxRadius:maxR_;

  int push;
  std::map<int, double>::iterator aDirective;
  std::map<int, LayerOption>::iterator anOption;

  LayerVector thisBarrelLayerSet;
  std::ostringstream layerName;
  BarrelLayer* aBarrelLayer;
  std::vector<double> prevDsDistances;

  std::pair<double, double> worstCaseRadii = sameRods ? std::make_pair(minRadius, maxRadius) : std::make_pair(-1.0, -1.0);

  if (maxZ != 0 && nModules != 0) { logWARNING("Barrel " + barrelName + " has both MaxZ and nModules defined. nModules will be ignored."); }

  const BarrelLayer::ModulePlacementStrategy& placeWithMaxZ = BarrelLayer::PlaceWithMaxZ(maxZ);
  const BarrelLayer::ModulePlacementStrategy& placeWithNumModules = BarrelLayer::PlaceWithNumModules(nModules);
  const BarrelLayer::ModulePlacementStrategy& moduleStrategy = maxZ > 0 ? placeWithMaxZ : placeWithNumModules;

  if (sameRods) { logINFO("Building barrel " + barrelName + " with identical rods"); }

  for (int i=0; i<nLayer; i++) {
    double radius = minRadius + (nLayer > 1 ? (maxRadius-minRadius)/double(nLayer-1)*i : 0);

    std::vector<double> geometryDsDistances = getGeometryDsDistances(barrelName, i+1, maxZ > 0 ? REASONABLE_MAX_ROD_MODULES : nModules);
    if (i > 0 && sameRods && !std::equal(prevDsDistances.begin(), prevDsDistances.end(), geometryDsDistances.begin())) {
      logWARNING("Requested identical rods for barrel " + barrelName + " but layer " + any2str(i+1) + " has different dsDistances than previous layers in the Types file. Using dsDistances found for layer 1.");
      geometryDsDistances = prevDsDistances;
    } else {
      prevDsDistances = geometryDsDistances;
    }

    sampleModule->setLayer(i+1);

    if ((i==0)||(i==(nLayer-1))) {
      push = Layer::FIXED;
    } else {
      push = Layer::AUTO;
    }

    if ((i==0)||(i==(nLayer-1))) {
      push = Layer::FIXED;
    } else {
      push = Layer::AUTO;
    }


    aDirective = layerDirectives_.find(i+1);
    if (aDirective!=layerDirectives_.end()) {
      if  ((i==0)||(i==(nLayer-1))) {
        logWARNING("We just read a directive for the first or last layer. This will be ignored");
        /*std::cout << "*******************************" << std::endl;
          std::cout << "*                             *" << std::endl;
          std::cout << "* WARNING:              /\\    *" << std::endl;
          std::cout << "*                      /!!\\   *" << std::endl;
          std::cout << "* We just read a      / !! \\  *" << std::endl;
          std::cout << "* directive for the   ^^^^^^  *" << std::endl;
          std::cout << "* first or last layer...      *" << std::endl;
          std::cout << "*                             *" << std::endl;
          std::cout << "*******************************" << std::endl;*/
      }
      std::ostringstream tempString;
      tempString.str(""); tempString << "Found a directive: " << layerDirectives_[i+1];
      logINFO(tempString.str());
      if (layerDirectives_[i+1]>0) {
        radius = layerDirectives_[i+1];
        push   = Layer::FIXED;
        tempString.str(""); tempString << "Fixing radius of layer " << i+1 << " to " << radius;
        logINFO(tempString.str());
      } else {
        push = int(layerDirectives_[i+1]);
      }
    } else {
      logINFO("Found no directive: auto adjusting layer");
    }

    aBarrelLayer = new BarrelLayer(sampleModule);
    layerName.str("");
    layerName << "L" << std::dec << i+1;
    aBarrelLayer->setName(layerName.str(), i+1);
    aBarrelLayer->setContainerName(barrelName);

    std::ostringstream tempString;
    tempString.str(""); tempString << "Desired radius: " << radius;
    logINFO(tempString.str());

    if (!shortBarrel) { // Standard Barrel
      aBarrelLayer->buildLayer(radius,       // averageRadius
                               worstCaseRadii,
                               getBigDelta(i+1),
                               getSmallDelta(i+1) ,
                               geometryDsDistances,
                               overlap_,     // overlap
                               zError_,      // safetyOrigin
                               moduleStrategy,
                               push,
                               phiSegments_, // modules multiple of ...
                               true,        // false = Strings with opposite parity
                               sampleModule, 
                               section);

      addLayer(aBarrelLayer, barrelName, TypeBarrel);
      thisBarrelLayerSet.push_back(aBarrelLayer);
    } else { // Mezzanine Barrel
      double farthestZ=getMaxBarrelZ(+1); // get the farthest Z reached up to now (WARNING: this implies that mezzanine barrels need to be declared LAST in the config file!)
      aBarrelLayer->buildLayer(radius,       // averageRadius
                               worstCaseRadii,
                               getBigDelta(i+1),
                               getSmallDelta(i+1) ,
                               geometryDsDistances,
                               overlap_,     // overlap
                               zError_,      // safetyOrigin
                               moduleStrategy,    
                               push,
                               phiSegments_, // modules multiple of ...
                               true,        // false = Strings with opposite parity
                               sampleModule, 
                               section, 
                               farthestZ);

      addLayer(aBarrelLayer, barrelName, TypeBarrel);
      thisBarrelLayerSet.push_back(aBarrelLayer);
    }

    anOption = layerOptions_.find(i+1);
    if (anOption!=layerOptions_.end()) {
      BarrelLayer* anotherLayer;
      LayerOption myOption=layerOptions_[i+1];
      if (myOption.first==Layer::Stacked) {
        anotherLayer = new BarrelLayer(*aBarrelLayer);
        anotherLayer->shiftRho(myOption.second);
        addLayer(anotherLayer, barrelName, TypeBarrel);
        thisBarrelLayerSet.push_back(anotherLayer);
      }
    }

  }

  if (compressed) {
    compressBarrelLayers(thisBarrelLayerSet, shortBarrel);
  }


  // Mezzanine barrel needs to be duplicated and reflected
  /* if (shortBarrel) { // CUIDADO: made obsolete by new mezzanine string generation code - TBR
     LayerVector::iterator layIt;
     BarrelLayer* anotherLayer;
     LayerVector justDoneBarrelLayerSet=thisBarrelLayerSet;

     for (layIt = justDoneBarrelLayerSet.begin(); layIt!= justDoneBarrelLayerSet.end(); layIt++) {
     if ( (aBarrelLayer=dynamic_cast<BarrelLayer*>(*layIt)) ) {
     anotherLayer = new BarrelLayer(*aBarrelLayer);
     anotherLayer->reflectZ();
     addLayer(anotherLayer, barrelName, TypeBarrel);
     thisBarrelLayerSet.push_back(anotherLayer);
     }
     }
     }*/

  maxZ=getMaxBarrelZ(+1);
  maxL_=(maxZ>maxL_)?maxZ:maxL_;
  // TODO: update this value if you want an independent compacting of the barrel section

  LpB_.push_back(nLayer);
  rMinpB_.push_back(minRadius);
  rMaxpB_.push_back(maxRadius);
  dZpB_.push_back(maxZ);

  return thisBarrelLayerSet;
}

// Barrel "compactification"
// destZ is 0 by default
void Tracker::compressBarrelLayers(LayerVector aLayerSet, bool oneSided, double destZ ) {
  LayerVector::iterator layIt;
  BarrelLayer* aBarrelLayer;

  double minZm = 0;
  double minZp = 0;
  double aZp;
  double aZm;

  // Take the shortest barrel
  for (layIt = aLayerSet.begin(); layIt!= aLayerSet.end(); layIt++) {
    if ( (aBarrelLayer=dynamic_cast<BarrelLayer*>(*layIt)) ) {
      aZp = aBarrelLayer->getMaxZ(+1);
      aZm = aBarrelLayer->getMaxZ(-1);
      // std::cerr << "it's a barrel layer in the range " << aZm << ".." << aZp; // debug
      if (layIt==aLayerSet.begin()) {
        minZp=aZp;
        minZm=aZm;
      } else {
        if (aZm>minZm) {
          minZm=aZm;
        }
        if (aZp<minZp) {
          minZp=aZp;
        }
      }
    }
    // std::cerr << std::endl; //debug
  }

  // std::cerr << "Shortest layer on minus is " << minZm << std::endl; // debug
  // std::cerr << "Shortest layer on plus  is " << minZp << std::endl; // debug

  double minZt;
  double compactOrigin;
  if (!oneSided) { // Standard barrel compressing
    minZt = (minZp>minZm) ? minZp : minZm;
    compactOrigin=0.;
  } else { // Mezzanine barrel compressing
    // Compact to the value with higher fabs()
    // use the other one as zero reference
    if (fabs(minZp)>fabs(minZm)) {
      minZt = minZp;
      compactOrigin = minZm;
    } else {
      minZt = minZp;
      compactOrigin = minZm;
    }
  }

  if (destZ!=0) minZt=destZ;

  // std::cerr << "compact origin : " << compactOrigin << std::endl; // debug
  // std::cerr << "compact to z   : " << minZt << std::endl; // debug

  for (layIt = aLayerSet.begin(); layIt!= aLayerSet.end(); layIt++) {
    if ( (aBarrelLayer=dynamic_cast<BarrelLayer*>(*layIt)) ) {
      if (!oneSided) { // Normal barrel
        aBarrelLayer->compressToZ(minZt);
      } else { // Mezzanine barrel
        aBarrelLayer->compressExceeding(minZt, compactOrigin);
      }
    } else {
      std::cerr << "ERROR: trying to compact a non-barrel layer" ;
    }
  }
}

void Tracker::alignShortBarrels() {
  if (barrelLayerSet_.size() > 1) {
    bool is_short, is_long;
    LayerVector::iterator iter = barrelLayerSet_.begin();
    LayerVector::iterator guard = barrelLayerSet_.end();
    LayerVector::iterator first_long = iter;
    // find first long layer for start and stop
    is_long = ((*first_long)->getMinZ() < 0) && ((*first_long)->getMaxZ() > 0);
    while (!is_long) {
      first_long++;
      if (first_long == guard) break;
      is_long = ((*first_long)->getMinZ() < 0) && ((*first_long)->getMaxZ() > 0);
    }
    while (iter != guard) {
      is_short = ((*iter)->getMaxZ() < 0) || ((*iter)->getMinZ() > 0);
      if (is_short) {
        bool change;
        LayerVector::iterator start;
        LayerVector::iterator stop;
        LayerVector::iterator cmp;
        start = first_long;
        stop = first_long;
        cmp = first_long;
        if (cmp != guard) cmp++;
        while (cmp != guard) {
          is_long = ((*cmp)->getMinZ() < 0) && ((*cmp)->getMaxZ() > 0);
          if (is_long) {
            if ((*iter)->getMinRho() > (*cmp)->getMinRho()) {
              change = ((*start)->getMinRho() > (*iter)->getMinRho()) || ((*cmp)->getMinRho() > (*start)->getMinRho());
            }
            else change = false;
            if (change) start = cmp;
          }
          cmp++;
        }
        if ((start == first_long) && ((*start)->getMinRho() > (*iter)->getMinRho())) start = guard;
        cmp = stop;
        if (cmp != guard) cmp++;
        while (cmp != guard) {
          is_long = ((*cmp)->getMinZ() < 0) && ((*cmp)->getMaxZ() > 0);
          if (is_long) {
            if ((*iter)->getMinRho() < (*cmp)->getMinRho()) {
              change = ((*stop)->getMinRho() < (*iter)->getMinRho()) || ((*cmp)->getMinRho() < (*stop)->getMinRho());
            }
            else change = false;
            if (change)  stop = cmp;
          }
          cmp++;
        }
        if ((stop == first_long) && ((*stop)->getMinRho() < (*iter)->getMinRho())) stop = guard;
        if ((stop == guard) && (start != guard)) {
          if ((*iter)->getMinZ() > 0) { // right of the origin
            XYZVector dz(0, 0, (*start)->getMaxZ() - (*iter)->getMaxZ());
            (*iter)->translate(dz);
          }
          else { // left of the origin
            XYZVector dz(0, 0, (*start)->getMinZ() - (*iter)->getMinZ());
            (*iter)->translate(dz);
          }
        }
        else if((start == guard) && (stop != guard)) {
          if ((*iter)->getMinZ() > 0) { // right of the origin
            if ((*stop)->getMaxZ() < (*iter)->getMaxZ()) {
              XYZVector dz(0, 0, (*stop)->getMaxZ() - (*iter)->getMaxZ());
              (*iter)->translate(dz);
            }
          }
          else { // left of the origin
            if ((*stop)->getMinZ() > (*iter)->getMinZ()) {
              XYZVector dz(0, 0, (*stop)->getMinZ() - (*iter)->getMinZ());
              (*iter)->translate(dz);
            }
          }
        }
        else if ((start != guard) && (stop != guard)) {
          if ((*iter)->getMinZ() > 0) { // right of the origin
            if ((*start)->getMaxZ() == (*stop)->getMaxZ()) {
              XYZVector dz(0, 0, (*start)->getMaxZ() - (*iter)->getMaxZ());
              (*iter)->translate(dz);
            }
            else {
              double dist1, dist2;
              XYZVector dz;
              dist1 = (*start)->getMaxZ() - (*iter)->getMaxZ();
              dist2 = (*stop)->getMaxZ() - (*iter)->getMaxZ();
              if (fabs(dist1) < fabs(dist2)) dz.SetZ(dist1);
              else dz.SetZ(dist2);
              (*iter)->translate(dz);
            }
          }
          else { // left of the origin
            if ((*start)->getMinZ() == (*stop)->getMinZ()) {
              XYZVector dz(0, 0, (*start)->getMinZ() - (*iter)->getMinZ());
              (*iter)->translate(dz);
            }
            else {
              double dist1, dist2;
              XYZVector dz;
              dist1 = (*start)->getMinZ() - (*iter)->getMinZ();
              dist2 = (*stop)->getMinZ() - (*iter)->getMinZ();
              if (fabs(dist1) < fabs(dist2)) dz.SetZ(dist1);
              else dz.SetZ(dist2);
              (*iter)->translate(dz);
            }
          }
        }
      }
      iter++;
    }
  }
}

void Tracker::sortLayers() {
  std::stable_sort(barrelLayerSet_.begin(), barrelLayerSet_.end(), smallerZ);
  std::stable_sort(barrelLayerSet_.begin(), barrelLayerSet_.end(), smallerRho);
  std::stable_sort(endcapLayerSet_.begin(), endcapLayerSet_.end(), smallerZ);
}

double Tracker::getMaxBarrelZ(int direction) {
  double maxZ = 0;
  double aZ;
  LayerVector::iterator layIt;
  BarrelLayer* aBarrelLayer;

  if (direction==0) {
    std::cerr << "Tracker::getMaxBarrelZ was called with zero direction. Assuming +1" << std::endl;
    direction=+1;
  }
  direction/=int(fabs(direction));

  // Take the shortest barrel
  for (layIt = barrelLayerSet_.begin(); layIt!= barrelLayerSet_.end(); layIt++) {
    if ( (aBarrelLayer=dynamic_cast<BarrelLayer*>(*layIt)) ) {
      aZ = aBarrelLayer->getMaxZ(direction);
      if (layIt==barrelLayerSet_.begin()) {
        maxZ=aZ;
      } else {
        if (direction*aZ>direction*maxZ) {
          maxZ=aZ;
        }
      }
    }
  }

  return maxZ;
}

void Tracker::buildEndcapsAtEta(int nDisks, int nRings, double minZ, double maxZ, double maxEta, double maxRadius,
                                std::map<int, EndcapModule*> sampleModule, std::string endcapName, int diskParity,
                                bool oddSegments, bool alignEdges,
                                int sectioned /* = Layer::NoSection */ ) {

  double minTheta = 2*atan(exp(-1*maxEta));
  double minRadius = minZ * tan(minTheta);

  buildEndcaps(nDisks, nRings, minZ, maxZ, minRadius, maxRadius,
               sampleModule, endcapName, diskParity, oddSegments, alignEdges, sectioned );

}


void Tracker::buildEndcaps(int nDisks, int nRings, double minZ, double maxZ, double minRadius, double maxRadius,
                           std::map<int, EndcapModule*> sampleModule, std::string endcapName, int diskParity,
                           bool oddSegments, bool alignEdges,
                           int sectioned /* = Layer::NoSection */ ) {

  maxR_=(maxRadius>maxR_)?maxRadius:maxR_;
  maxL_=(maxZ>maxL_)?maxZ:maxL_;

  double thisZ;
  double deltaZ;

  // Geometric progression factor
  double alpha = pow(maxZ/minZ, 1/double(nDisks-1));

  EndcapLayer* defaultDisk = new EndcapLayer();
  EndcapLayer* anotherDisk;

  int estimatedNrings = nRings > 0 ? nRings : REASONABLE_MAX_DISK_RINGS;
  std::vector<double> geometryDsDistances(estimatedNrings, 0);
  for (int i = 1; i <= nDisks; i++) {
    std::vector<double> tempDsDistances = getGeometryDsDistances(endcapName, i, estimatedNrings);
    for (int j = 0; j < estimatedNrings; j++) {
      if (tempDsDistances[j] > geometryDsDistances[j]) geometryDsDistances[j] = tempDsDistances[j]; 
    }
  }
  
  if (nRings > 0) {
    if (minRadius > 0) { logWARNING("Endcap " + endcapName + " has both nRings and innerRadius defined. innerRadius will be ignored."); } 
    for (std::map<int, EndcapModule*>::iterator it = sampleModule.begin(); it != sampleModule.end(); ++it) {
      if (it->second->getShape()==Module::Wedge) { logERROR("Option nRings is incompatible with wedge-shaped modules. Endcap " + endcapName + " might be built incorrectly."); break; }
    }
    if (diskParity == -1) { logWARNING("Endcap " + endcapName + " will be built top-to-bottom, but has diskParity = -1. This will lead to non-optimal coverage."); }
    logINFO("Endcap " + endcapName + " will be built top-to-bottom with a fixed number of rings (" + any2str(nRings) + ").");
    defaultDisk->buildSingleDisk( nRings, maxRadius, smallDelta_,
                                  bigDelta_, (minZ+maxZ)/2, overlap_,
                                  zError_+(maxZ-minZ)/2,
                                  geometryDsDistances,
                                  phiSegments_, // Base
                                  oddSegments, alignEdges,
                                  sampleModule,
                                  ringDirectives_,
                                  diskParity,
                                  sectioned );
  } else {
    logINFO("Endcap " + endcapName + " will be built bottom-to-top with a number of rings depending on the innerRadius (" + any2str(minRadius) + ").");
    defaultDisk->buildSingleDisk( minRadius, maxRadius, smallDelta_,
                                  bigDelta_, (minZ+maxZ)/2, overlap_,
                                  zError_+(maxZ-minZ)/2,
                                  geometryDsDistances,
                                  phiSegments_, // Base
                                  oddSegments, alignEdges,
                                  sampleModule,
                                  ringDirectives_,
                                  diskParity,
                                  sectioned );
  }

  std::ostringstream layerName;
  EndcapModule* anEndcapModule;

  for (int iDisk=0; iDisk<nDisks; iDisk++) {
    // Set the disk number in all the modules
    for (ModuleVector::iterator modIt = defaultDisk->getModuleVector()->begin();
         modIt!=defaultDisk->getModuleVector()->end();
         modIt++) {
      if ( (anEndcapModule=dynamic_cast<EndcapModule*>(*modIt)) ) {
        anEndcapModule->setDisk(iDisk+1);
      } else {
        std::cerr << "ERROR IN Tracker::buildEndcaps this shoundn't happen!" << std::endl;
      }
    }

    layerName.str("");
    layerName << "D" << std::dec << iDisk+1;
    thisZ = pow(alpha, iDisk) * minZ;
    deltaZ=-1*(minZ+maxZ)/2+thisZ;
    anotherDisk = new EndcapLayer(*defaultDisk);
    anotherDisk->setName(layerName.str(), iDisk+1);
    anotherDisk->setContainerName(endcapName);
    anotherDisk->translateZ(deltaZ);
    addLayer(anotherDisk, endcapName, TypeEndcap);
    anotherDisk = new EndcapLayer(*anotherDisk);
    anotherDisk->rotateY_PI();
    addLayer(anotherDisk, endcapName, TypeEndcap);
  }

  for (ModuleVector::iterator modIt = defaultDisk->getModuleVector()->begin();
       modIt!=defaultDisk->getModuleVector()->end();
       modIt++) {
    endcapSample_.push_back(*modIt);
  }

  // TODO: decide how to handle this
  // delete defaultDisk;

  DpE_.push_back(nDisks);
  rMinpE_.push_back(minRadius);
  rMaxpE_.push_back(maxRadius);
  dZpE_.push_back((maxZ - minZ) / 2.0);
}

// Function used to remove some endcaps rings
// sectionName: name of the endcap to operate onto
// iDisk: number of the disk onto which operate
// iRing: number of the first ring to remove
// directionOuter: if it's true we must remove all the rings outer than iRing
// while if it's false we remove all the rings inner that iRing
void Tracker::removeDiskRings(std::string sectionName, int iDisk, int iRing, bool directionOuter) {
  LayerVector::iterator layIt;
  ModuleVector::iterator modIt;
  ModuleVector::iterator anotherModIt;
  ModuleVector* aLay;
  Module* aModule;
  EndcapModule* anEndcapModule;

  // Take the vector of layers in the section
  LayerVector myLayers = sectionMap_[sectionName];

  for (layIt=myLayers.begin(); layIt!=myLayers.end(); layIt++) {
    aLay = (*layIt)->getModuleVector();
    for (modIt=aLay->begin(); modIt!=aLay->end(); modIt++) {
      aModule=(*modIt);

      if ( (anEndcapModule=dynamic_cast<EndcapModule*>(aModule)) ) {
        if (directionOuter) {
          if ((anEndcapModule->getDisk()==iDisk)
              && (anEndcapModule->getRing()>=iRing)) {
            delete aModule;
            anotherModIt=modIt-1;
            aLay->erase(modIt);
            modIt=anotherModIt;
            (*layIt)->decreaseModCount(anEndcapModule->getRing() - 1);
          }
        } else {
          if ((anEndcapModule->getDisk()==iDisk)
              && (anEndcapModule->getRing()<=iRing)) {
            delete aModule;
            anotherModIt=modIt-1;
            aLay->erase(modIt);
            modIt=anotherModIt;
            (*layIt)->decreaseModCount(anEndcapModule->getRing() - 1);
          }
        }
      } else if (dynamic_cast<BarrelModule*>(aModule)) {
        // ERROR: this should not happen
        std::cerr << "ERROR: a barrel module was found in section " << sectionName
          << " while we are trying to remove rings from there. It should be an endcap module!" << std::endl;
      }
    }
  }
}

// *******************************
// *                             *
// * Geometry analysis functions *
// *                             *
// *******************************

std::pair<double, double> Tracker::getEtaMinMax() {
  std::pair<double, double> result;
  LayerVector::iterator layIt;
  ModuleVector* moduleV;
  ModuleVector::iterator modIt;

  double theta;
  double minTheta=M_PI+1; // (!) :-)
  double maxTheta=-1;     // idem...

  for (layIt=layerSet_.begin(); layIt!=layerSet_.end(); layIt++) {
    moduleV = (*layIt)->getModuleVector();
    for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
      theta=(*modIt)->getMinTheta();
      if (theta<minTheta) minTheta=theta;
      theta=(*modIt)->getMaxTheta();
      if (theta>maxTheta) maxTheta=theta;
    }
  }

  result.first = -1*log(tan(maxTheta/2.));
  result.second = -1*log(tan(minTheta/2.));

  return result;
}

int Tracker::cutOverEta(double etaCut) {
  int nCut = 0;
  LayerVector::iterator layIt;
  ModuleVector::iterator modIt;

  for (layIt=layerSet_.begin(); layIt!=layerSet_.end(); layIt++) {
    nCut += (*layIt)->cutOverEta(etaCut);
  }

  return nCut;
}


ModuleVector Tracker::trackHit(const XYZVector& origin, const XYZVector& direction, ModuleVector* moduleV) {
  ModuleVector result;
  ModuleVector::iterator modIt;
  double distance;

  for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
    // A module can be hit if it fits the phi (precise) contraints
    // and the eta constaints (taken assuming origin within 5 sigma)
    if ((*modIt)->couldHit(direction.Eta(), direction.Phi())) {
      distance=(*modIt)->trackCross(origin, direction);
      if (distance>0) {
        result.push_back(*modIt);
      }
    }
  }

  return result;
}



// Resets a module type counter
void Tracker::resetTypeCounter(std::map <std::string, int> &modTypes) {
  for (std::map <std::string, int>::iterator it = modTypes.begin();
       it!=modTypes.end(); it++) {
    (*it).second = 0;
  }
}


// Shoots directions with random (flat) phi, random (flat) pseudorapidity
// gives also the direction's eta
std::pair <XYZVector, double > Tracker::shootFixedDirection(double phi, double eta) {
  std::pair <XYZVector, double> result;

  double theta=2*atan(exp(-1*eta));

  // Direction
  result.first  = XYZVector(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
  result.second = eta;
  return result;
}

// Shoots directions with random (flat) phi, random (flat) pseudorapidity
// gives also the direction's eta
std::pair <XYZVector, double > Tracker::shootDirection(double minEta, double spanEta) {
  std::pair <XYZVector, double> result;

  double eta;
  double phi;
  double theta;

  // phi is random [0, 2pi)
  phi = myDice_.Rndm() * 2 * M_PI; // debug

  // eta is random (-4, 4]
  eta = myDice_.Rndm() * spanEta + minEta;
  theta=2*atan(exp(-1*eta));

  // Direction
  result.first  = XYZVector(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
  result.second = eta;
  return result;
}

// Shoots directions with random (flat) phi, random (flat) pseudorapidity
// gives also the direction's eta
std::pair <XYZVector, double > Tracker::shootDirectionFixedPhi(double minEta, double spanEta) {
  std::pair <XYZVector, double> result;

  double eta;
  double theta;

  // eta is random (-4, 4]
  eta = myDice_.Rndm() * spanEta + minEta;
  theta=2*atan(exp(-1*eta));

  // Direction
  result.first  = XYZVector(0, sin(theta), cos(theta));
  result.second = eta;
  return result;
}

// Prints the positions of barrel modules to file or cout
void Tracker::printBarrelModuleZ(ostream& outfile) {
  ModuleVector* myModules;
  ModuleVector::iterator itModule;
  LayerVector* myBarrels;
  LayerVector::iterator itLayer;
  std::pair<int, int> myRZ;
  int myR, myZ;
  XYZVector meanPoint;

  int mmFraction = 1000;

  outfile << "BarrelLayer name, r(mm), z(mm), number of modules" <<std::endl;
  myBarrels = getBarrelLayers();
  for (itLayer = myBarrels->begin();
       itLayer != myBarrels->end();
       itLayer++) {

    std::map< std::pair<int,int>, int > posCount;

    myModules = (*itLayer)->getModuleVector();
    for (itModule = myModules->begin();
         itModule != myModules->end();
         itModule++) {
      meanPoint = (*itModule)->getMeanPoint();
      myR = int(ceil(meanPoint.Rho()*mmFraction-0.5));
      myZ = int(ceil(meanPoint.Z()*mmFraction-0.5));
      myRZ.first=myR;
      myRZ.second=myZ;
      posCount[myRZ]++;
    }

    std::map<std::pair<int,int>, int>::iterator itPos;
    for (itPos = posCount.begin();
         itPos != posCount.end();
         itPos++) {
      // BarrelLayer name
      outfile << (*itLayer)->getContainerName() << "-"
        << (*itLayer)->getName() << ", ";

      outfile << (*itPos).first.first/(double)mmFraction << ", " // r
        << (*itPos).first.second/(double)mmFraction << ", "   // z
        << (*itPos).second // number of modules
        << std::endl;
    }
  }
}

// Prints the positions of endcap modules to file or cout
void Tracker::printEndcapModuleRPhiZ(ostream& outfile) {
  ModuleVector::iterator itModule;
  XYZVector meanPoint;
  double base_inner, base_outer, height;

  EndcapModule* anEC;
  int aRing;
  // Sort the endcap module vector to have it ordered by
  // rings (order by r, then phi, then z)
  std::sort(endcapSample_.begin(), endcapSample_.end(), moduleSortEndcapStyle);

  // Compute the average Z as the mean between max and min Z
  double maxZ=0;
  double minZ=0;
  bool first=true;
  for (ModuleVector::iterator moduleIt=endcapSample_.begin(); moduleIt!=endcapSample_.end(); moduleIt++) {
    if ( (anEC=dynamic_cast<EndcapModule*>(*moduleIt)) ) {
      if (first) {
        first=false;
        minZ=anEC->getMeanPoint().Z();
        maxZ=minZ;
      } else {
        if (anEC->getMeanPoint().Z()<minZ) minZ=anEC->getMeanPoint().Z();
        if (anEC->getMeanPoint().Z()>maxZ) maxZ=anEC->getMeanPoint().Z();
      }
    } else {
      std::cerr << "ERROR: found a non-Endcap module in the map of ring types. This should not happen. Contact the developers." << std::endl;
    }
  }
  double averageZ = (minZ+maxZ)/2.;

  // Look into the endcap sample to get modules' positions 
  outfile << "Ring, r(mm), phi(deg), z(mm), base_inner(mm), base_outer(mm), height(mm)" <<std::endl;
  for (ModuleVector::iterator moduleIt=endcapSample_.begin(); moduleIt!=endcapSample_.end(); moduleIt++) {
    if ( (anEC=dynamic_cast<EndcapModule*>(*moduleIt)) ) {
      aRing=anEC->getRing();
      meanPoint = anEC->getMeanPoint();
      base_inner = anEC->getWidthLo();
      base_outer = anEC->getWidthHi();
      height = anEC->getHeight();

      // Print the data in fixed-precision
      // Limit the precision to one micron for lengths and 1/1000 degree for angles
      outfile << std::fixed;
      outfile << aRing << ", " 
        << std::fixed << std::setprecision(3) << meanPoint.Rho() << ", "
        << std::fixed << std::setprecision(3) << meanPoint.Phi()/M_PI*180. << ", "
        << std::fixed << std::setprecision(3) << meanPoint.Z()-averageZ << ", "
        << std::fixed << std::setprecision(3) << base_inner << ", "
        << std::fixed << std::setprecision(3) << base_outer << ", "
        << std::fixed << std::setprecision(3) << height << std::endl;
    } else {
      std::cerr << "ERROR: found a non-Endcap module in the map of ring types. This should not happen. Contact the developers." << std::endl;
    }
  }
}


void Tracker::printHtmlTableRow(ofstream *output, std::vector<std::string> myRow) {
  std::vector<std::string>::iterator strIt;
  (*output) << "<tr>" << std::endl;
  for (strIt=myRow.begin(); strIt!=myRow.end(); strIt++) {
    (*output) << "<td>" << (*strIt) << "</td> ";
  }
  (*output) << "</tr>" << std::endl;
}

void Tracker::printHtmlTableRow(ofstream *output, std::vector<double> myRow, int coordPrecision /* = 0 */, bool skimZero /*=false*/) {
  std::vector<double>::iterator strIt;
  (*output) << "<tr>" << std::endl;
  for (strIt=myRow.begin(); strIt!=myRow.end(); strIt++) {
    (*output) << "<td>";
    if ((!skimZero)||((*strIt)!=0)) (*output) << std::fixed << std::setprecision(coordPrecision) << (*strIt);
    (*output)<< "</td>";
  }
  (*output) << "</tr>" << std::endl;
}

// TODO: use the boost library to handle directories
// (see http://boost.org/libs/filesystem/doc/index.htm )
void Tracker::createDirectories() {

  if (activeDirectory_=="") {
    activeDirectory_ = summaryDirectory_ + "/" + trackerName_;
    mkdir(summaryDirectory_.c_str(), 0755);
    mkdir(activeDirectory_.c_str(), 0755);
    mkdir(storeDirectory_.c_str(), 0755);
  } else {
    mkdir(activeDirectory_.c_str(), 0755);
    mkdir(storeDirectory_.c_str(), 0755);
  }

}

void Tracker::save() {
  std::vector<TObject* >::iterator itObj;
  std::string fileName;

  fileName = storeDirectory_ + "/" + trackerName_ + ".root";
  createDirectories();
  TFile* myFile = new TFile(fileName.c_str(), "RECREATE", trackerName_.c_str());
  for (itObj = savingV_.begin(); itObj != savingV_.end(); itObj++) {
    (*itObj)->Write();
  }
  myFile->Close();
}


// Enumerate sections by
// axis index normal to draw plane
// (if x=1, y=2, z=3)
void Tracker::drawGrid(double maxL, double maxR, int noAxis/*=1*/, double spacing /*= 100.*/, Option_t* option /*= "same"*/) {
  TPolyLine3D* aLine;
  Color_t gridColor = COLOR_GRID;
  Color_t gridColor_hard = COLOR_HARD_GRID;
  Color_t thisLineColor;

  std::string theOption(option);

  int i;
  int j;
  int k;

  double topMax = (maxL > maxR) ? maxL : maxR;
  topMax = ceil(topMax/spacing)*spacing;

  double aValue[3];
  double minValue[3];
  double maxValue[3];
  double runValue;
  int thisLineStyle;

  i=(noAxis)%3;
  j=(noAxis+1)%3;
  k=(noAxis+2)%3;

  maxL *= 1.1;
  maxR *= 1.1;

  if (noAxis==1) {
    minValue[0]=0;
    maxValue[0]=+maxR;
    minValue[1]=0;
    maxValue[1]=+maxR;
    minValue[2]=0;
    maxValue[2]=+maxL;
  } else {
    minValue[0]=-maxR;
    maxValue[0]=+maxR;
    minValue[1]=-maxR;
    maxValue[1]=+maxR;
    minValue[2]=0;
    maxValue[2]=+maxL;
  }

  aValue[k]=-topMax;
  for(runValue = -topMax; runValue<=topMax; runValue+=spacing) {

    // Special line for axis
    if (runValue==0) {
      thisLineStyle=1;
      thisLineColor=gridColor_hard;
    } else {
      thisLineStyle=2;
      thisLineColor=gridColor;
    }

    // Parallel to j
    if ((runValue<=maxValue[i])&&(runValue>=minValue[i])) {
      aValue[i] = runValue;
      aLine = new TPolyLine3D(2);
      aValue[j] = minValue[j];
      aLine->SetPoint(0, aValue[0], aValue[1], aValue[2]);
      aValue[j] = maxValue[j];
      aLine->SetPoint(1, aValue[0], aValue[1], aValue[2]);
      aLine->SetLineStyle(thisLineStyle);
      aLine->SetLineColor(thisLineColor);
      aLine->Draw(theOption.c_str());
      theOption="same";
    };

    // Parallel to i
    if ((runValue<=maxValue[j])&&(runValue>=minValue[j])) {
      aValue[j] = runValue;
      aLine = new TPolyLine3D(2);
      aValue[i] = minValue[i];
      aLine->SetPoint(0, aValue[0], aValue[1], aValue[2]);
      aValue[i] = maxValue[i];
      aLine->SetPoint(1, aValue[0], aValue[1], aValue[2]);
      aLine->SetLineStyle(thisLineStyle);
      aLine->SetLineColor(thisLineColor);
      aLine->Draw(theOption.c_str());
      theOption="same";
    };

  }


}

// Enumerate sections by
// axis index normal to draw plane
// (if x=1, y=2, z=3)
void Tracker::drawTicks(TView* myView, double maxL, double maxR, int noAxis/*=1*/, double spacing /*= 100.*/, Option_t* option /*= "same"*/) {
  TPolyLine3D* aLine;
  Color_t gridColor_hard = COLOR_HARD_GRID;
  int gridStyle_solid = 1;

  std::string theOption(option);

  double topMax = (maxL > maxR) ? maxL : maxR;
  topMax = ceil(topMax/spacing)*spacing;

  maxL *= 1.1;
  maxR *= 1.1;

  if (noAxis==1) {
    double etaStep=.2;
    double etaMax = 2.1;
    // Add the eta ticks
    double theta;
    double tickLength = 2 * spacing;
    double tickDistance = spacing;
    double startR = maxR + tickDistance;
    double startL = maxL + tickDistance;
    double endR = maxR + tickDistance + tickLength;
    double endL = maxL + tickDistance + tickLength;
    XYZVector startTick;
    XYZVector endTick;
    Double_t pw[3];
    Double_t pn[3];
    TText* aLabel;
    char labelChar[10];
    for (double eta=0; eta<etaMax; eta+=etaStep) {
      aLine = new TPolyLine3D(2);
      theta = 2 * atan(exp(-eta));
      startTick = XYZVector(0, sin(theta), cos(theta));
      startTick *= startR/startTick.Rho();
      endTick = startTick / startTick.Rho() * endR;
      if (startTick.Z()>startL) {
        startTick *= startL/startTick.Z();
        endTick *=  endL/endTick.Z();
      }
      pw[0]=0.;
      pw[1]=endTick.Y();
      pw[2]=endTick.Z();
      myView->WCtoNDC(pw, pn);
      sprintf(labelChar, "%.01f", eta);
      aLabel = new TText(pn[0], pn[1], labelChar);
      aLabel->SetTextSize(aLabel->GetTextSize()*.6);
      aLabel->SetTextAlign(21);
      aLabel->Draw(theOption.c_str());
      theOption="same";
      endTick = (endTick+startTick)/2.;
      aLine->SetPoint(0, 0., startTick.Y(), startTick.Z());
      aLine->SetPoint(1, 0., endTick.Y(), endTick.Z());
      aLine->SetLineStyle(gridStyle_solid);
      aLine->SetLineColor(gridColor_hard);
      aLine->Draw("same");
    }
    aLine = new TPolyLine3D(2);
    theta = 2 * atan(exp(-2.5));
    startTick = XYZVector(0, sin(theta), cos(theta));
    startTick *= startR/startTick.Rho();
    endTick = startTick / startTick.Rho() * endR;
    if (startTick.Z()>startL) {
      startTick *= startL/startTick.Z();
      endTick *=  endL/endTick.Z();
    }
    pw[0]=0.;
    pw[1]=endTick.Y();
    pw[2]=endTick.Z();
    myView->WCtoNDC(pw, pn);
    sprintf(labelChar, "%.01f", 2.5);
    aLabel = new TText(pn[0], pn[1], labelChar);
    aLabel->SetTextSize(aLabel->GetTextSize()*.8);
    aLabel->SetTextAlign(21);
    aLabel->Draw(theOption.c_str());
    theOption="same";
    endTick = (endTick+startTick)/2.;
    aLine->SetPoint(0, 0., 0., 0.);
    aLine->SetPoint(1, 0., endTick.Y(), endTick.Z());
    aLine->SetLineStyle(gridStyle_solid);
    aLine->SetLineColor(gridColor_hard);
    aLine->Draw("same");



    for (double z=0; z<=maxL ; z+=(4*spacing)) {
      aLine = new TPolyLine3D(2);
      startTick = XYZVector(0, 0, z);
      endTick = XYZVector(0, -(tickLength/2), z);
      aLine->SetPoint(0, 0., startTick.Y(), startTick.Z());
      aLine->SetPoint(1, 0., endTick.Y(), endTick.Z());
      pw[0]=0.;
      pw[1]=-tickLength;
      pw[2]=endTick.Z();
      myView->WCtoNDC(pw, pn);
      sprintf(labelChar, "%.0f", z);
      aLabel = new TText(pn[0], pn[1], labelChar);
      aLabel->SetTextSize(aLabel->GetTextSize()*.6);
      aLabel->SetTextAlign(23);
      aLabel->Draw(theOption.c_str());
      theOption="same";
      aLine->SetLineStyle(gridStyle_solid);
      aLine->SetLineColor(gridColor_hard);
      aLine->Draw("same");
    }

    for (double y=0; y<=maxR ; y+=(2*spacing)) {
      aLine = new TPolyLine3D(2);
      startTick = XYZVector(0, y, 0);
      endTick = XYZVector(0, y, -(tickLength/2));
      aLine->SetPoint(0, 0., startTick.Y(), startTick.Z());
      aLine->SetPoint(1, 0., endTick.Y(), endTick.Z());
      pw[0]=0.;
      pw[1]=endTick.Y();
      pw[2]=-tickLength;
      myView->WCtoNDC(pw, pn);
      sprintf(labelChar, "%.0f", y);
      aLabel = new TText(pn[0], pn[1], labelChar);
      aLabel->SetTextSize(aLabel->GetTextSize()*.6);
      aLabel->SetTextAlign(32);
      aLabel->Draw(theOption.c_str());
      theOption="same";
      aLine->SetLineStyle(gridStyle_solid);
      aLine->SetLineColor(gridColor_hard);
      aLine->Draw("same");
    }
  }
}


void Tracker::setModuleTypes(std::string sectionName,
                             std::map<int, int> nStripsAcross,
                             std::map<int, int> nFaces,
                             std::map<int, int> nSegments,
                             std::map<int, std::string> myType,
                             std::map<int, double> dsDistance,
                             std::map<int, int> triggerWindow,
                             std::map<int, double> dsRotation,
                             std::map<int, int> divideBack,
                             std::map<int, double> xResolution,
                             std::map<int, double> yResolution,
                             std::map<std::pair<int, int>, int> nStripsAcrossSecond,
                             std::map<std::pair<int, int>, int> nFacesSecond,
                             std::map<std::pair<int, int>, int> nSegmentsSecond,
                             std::map<std::pair<int, int>, std::string> myTypeSecond,
                             std::map<std::pair<int, int>, double> dsDistanceSecond,
                             std::map<std::pair<int, int>, int> triggerWindowSecond,
                             std::map<std::pair<int, int>, double> dsRotationSecond,
                             std::map<std::pair<int, int>, int> divideBackSecond,
                             std::map<std::pair<int, int>, bool> specialSecond) {

  LayerVector::iterator layIt;
  ModuleVector::iterator modIt;
  ModuleVector* aLay;
  Module* aModule;
  BarrelModule* aBarrelModule;
  EndcapModule* anEndcapModule;

  std::map<int, bool> warningStrips;
  std::map<int, bool> warningFaces;
  std::map<int, bool> warningSegments;
  std::map<int, bool> warningType;
  std::map<int, bool> warningDistance;
  std::map<int, bool> warningRotation;


  int aStripsAcross;
  int aFaces;
  int aSegments;
  std::string aType;
  double aDistance;
  int aTriggerWindow;
  double aRotation;
  int aDivideBack;
  double aXResolution;
  double aYResolution;

  std::pair<int, int> mySpecialIndex;

  std::ostringstream myTag; // This must be set according to the sectionName and index
  int myReadoutType; // (this must be set according to the delcared "type" )
  int myIndex; // This must be set according to the module layer (barrel) or disk (endcap)

  // Take the vector of layers in the section
  LayerVector myLayers = sectionMap_[sectionName];

  for (layIt=myLayers.begin(); layIt!=myLayers.end(); layIt++) {
    aLay = (*layIt)->getModuleVector();
    for (modIt=aLay->begin(); modIt!=aLay->end(); modIt++) {
      aModule=(*modIt);

      // Tag definition and index picking definition
      myTag.str("");
      myTag << sectionName << std::dec ;
      myIndex = -1;
      if ( (aBarrelModule=dynamic_cast<BarrelModule*>(aModule)) ) {
        myTag << "L" << setfill('0') << setw(2) << aBarrelModule->getLayer();
        myIndex = aBarrelModule->getLayer();
        mySpecialIndex.first= myIndex;
        mySpecialIndex.second = aBarrelModule->getRing();
        if (specialSecond[mySpecialIndex]) {
          myTag << "R" << setfill('0') << setw(2) << aBarrelModule->getRing();
        }
      } else if ( (anEndcapModule=dynamic_cast<EndcapModule*>(aModule)) ) {
        myTag << "R" << setfill('0') << setw(2) << anEndcapModule->getRing();
        myIndex = anEndcapModule->getRing();
        // If special rules are applied here, we add the disk id to the tag
        mySpecialIndex.first = myIndex;
        mySpecialIndex.second = anEndcapModule->getDisk();
        if (specialSecond[mySpecialIndex]) {
          myTag << "D" << setfill('0') << setw(2) << anEndcapModule->getDisk();
        }
      } else {
        // This shouldnt happen
        std::cerr << "ERROR! in function Tracker::setModuleTypes() "
          << "I found a module which is not Barrel nor Endcap module. What's this?!?" << std::endl;
        mySpecialIndex.first= -1;
        mySpecialIndex.second = -1;
      }

      aStripsAcross = nStripsAcross[myIndex];
      aFaces = nFaces[myIndex];
      aSegments = nSegments[myIndex];
      aType = myType[myIndex];
      aDistance = dsDistance[myIndex];
      aTriggerWindow = triggerWindow[myIndex];
      aRotation = dsRotation[myIndex];
      aDivideBack = divideBack[myIndex];
      aXResolution = xResolution[myIndex];
      aYResolution = yResolution[myIndex];

      if (specialSecond[mySpecialIndex]) {
        if (nStripsAcrossSecond[mySpecialIndex]!=0) {
          aStripsAcross = nStripsAcrossSecond[mySpecialIndex];
        }
        if (nFacesSecond[mySpecialIndex]!=0) {
          aFaces = nFacesSecond[mySpecialIndex];
        }
        if (nSegmentsSecond[mySpecialIndex]!=0) {
          aSegments = nSegmentsSecond[mySpecialIndex];
        }
        if (myTypeSecond[mySpecialIndex]!="") {
          aType = myTypeSecond[mySpecialIndex];
        }
        if (dsDistanceSecond[mySpecialIndex]!=0) {
          aDistance = dsDistanceSecond[mySpecialIndex];
        }
        if (triggerWindowSecond[mySpecialIndex]!=0) {
          aTriggerWindow = triggerWindowSecond[mySpecialIndex];
        }
        if (dsRotationSecond[mySpecialIndex]!=0) {
          aRotation = dsRotationSecond[mySpecialIndex];
        }
        if (divideBackSecond[mySpecialIndex]!=0) {
          aDivideBack = divideBackSecond[mySpecialIndex];
        }
      }

      // Readout type definition, according to the module type
      if (aType == "pt") {
        myReadoutType = Module::Pt;
      } else if (myType[myIndex] == "pt2") {
        myReadoutType = Module::Pt;
      } else if (myType[myIndex] == "ptIn") {
        myReadoutType = Module::Pt;
      } else if (myType[myIndex] == "ptOut") {
        myReadoutType = Module::Pt;
      } else if (myType[myIndex] == "ptMixed") {
        myReadoutType = Module::Pt;
      } else if (myType[myIndex] == "rphi") {
        myReadoutType = Module::Strip;
      } else if (myType[myIndex] == "stereo") {
        myReadoutType = Module::Strip;
      } else if (myType[myIndex] == "pixel") {
        myReadoutType = Module::Pixel;
      } else {
        myReadoutType = Module::Undefined;
      }


      aModule->setNFaces(aFaces);
      aModule->setNStripsAcross(aStripsAcross);
      aModule->setNSegments(aSegments);
      aModule->setType(aType);
      aModule->setModuleType(&(mapType_[aType]));
      aModule->setStereoDistance(aDistance);
      aModule->setTriggerWindow(aTriggerWindow);
      aModule->setStereoRotation(aRotation);
      aModule->setTag(myTag.str());
      aModule->setContainerName(sectionName);
      // aModule->setColor(colorPicker(aType));
      aModule->setColor(Palette::color(aType));
      aModule->setReadoutType(myReadoutType);

      if (aXResolution) aModule->setResolutionRphi(aXResolution);
      else aModule->setResolutionRphi();
      if (aYResolution) {
        aModule->setResolutionY(aYResolution);
      } else {
        aModule->setResolutionY();
      }

      // Make the back side of the module with longer strips
      // if applicable (and complain otherwise)
      if (aDivideBack>1) {
        if (aFaces!=2) {
          std::cerr << "ERROR: trying to make modules with longer back strips (divideBack = )"
            << aDivideBack << ", but the number of faces is " << aFaces
            << " is not 2" << std::endl << "Ignoring the statement" << std::endl;
        } else {
          if ((aSegments%aDivideBack)!=0) {
            std::cerr << "ERROR: trying to divide the number of segments in the back sensor by "
              << aDivideBack << ", but the number of segments is " << aSegments
              << ", I cannot perform the integer division" << std::endl;
          } else {
            aModule->setNSegments(2, aSegments/aDivideBack);
          }
        }
      }

      // Check if a given module type was not assigned fr a group
      // We do not check the case for special assignment (optional)
      if ((nStripsAcross[myIndex]==0)&&(!warningStrips[myIndex])) {
        std::cerr << "WARNING: undefined or zero nStripsAcross: \"" << nStripsAcross[myIndex] << "\" "
          << "for tracker section " << sectionName << "[" << myIndex << "]" << std::endl;
        warningStrips[myIndex]=true;
      }
      if ((nFaces[myIndex]==0)&&(!warningFaces[myIndex])) {
        std::cerr << "WARNING: undefined or zero nFaces: \"" << nFaces[myIndex] << "\" "
          << "for tracker section " << sectionName << "[" << myIndex << "]" << std::endl;
        warningFaces[myIndex]=true;
      }
      if ((nSegments[myIndex]==0)&&(!warningSegments[myIndex])) {
        std::cerr << "WARNING: undefined or zero nSegments: \"" << nSegments[myIndex] << "\" "
          << "for tracker section " << sectionName << "[" << myIndex << "]" << std::endl;
        warningSegments[myIndex]=true;
      }
      if ((myReadoutType==Module::Undefined)&&(!warningType[myIndex])) {
        std::cerr << "WARNING: undefined or void module type: \"" << myType[myIndex] << "\" "
          << "for tracker section " << sectionName << "[" << myIndex << "]" << std::endl;
        warningType[myIndex]=true;
      }
      if ((dsDistance[myIndex]<0)&&(!warningDistance[myIndex])) {
        std::cerr << "WARNING: negative distance for stereo sensors: \"" << dsDistance[myIndex] << "\" "
          << "for tracker section " << sectionName << "[" << myIndex << "]" << std::endl;
        warningDistance[myIndex]=true;
      }
      if ((dsRotation[myIndex]>360.0)&&(!warningRotation[myIndex])) {
        std::cerr << "WARNING: rotation for stereo sensors is greater than 360*deg: \"" << dsRotation[myIndex] << "\" "
          << "for tracker section " << sectionName << "[" << myIndex << "]" << std::endl;
        warningRotation[myIndex]=true;
      }

    }
  }

}

void Tracker::setGeometryDsDistance(std::string cntName, int firstIndex, int secondIndex, double value) {
  geometryDsDistanceSecond_[cntName][std::make_pair(firstIndex, secondIndex)] = value;
}

void Tracker::setGeometryDsDistance(std::string cntName, int firstIndex, double value) {
  geometryDsDistance_[cntName][firstIndex] = value;
}


std::vector<double> Tracker::getGeometryDsDistances(std::string cntName, int index, int numModules) const { // index will be either a Layer in case of barrel or a Ring in case of endcap
  std::vector<double> dsDistances(numModules, 0);
  if (geometryDsDistance_.count(cntName) && geometryDsDistance_.at(cntName).count(0)) {  // DS DISTANCE OVERRIDE - whatever is written at index 0 takes precedence over any other setting
    dsDistances.assign(numModules, geometryDsDistance_.at(cntName).at(0));
    return dsDistances;
  }
  if (geometryDsDistance_.count(cntName) && geometryDsDistance_.at(cntName).count(index)) {
    dsDistances.assign(numModules, geometryDsDistance_.at(cntName).at(index));
  }
  if (geometryDsDistanceSecond_.count(cntName) > 0) {
    for (std::map<std::pair<int, int>, double>::const_iterator it = geometryDsDistanceSecond_.at(cntName).begin(); it != geometryDsDistanceSecond_.at(cntName).end(); ++it) {
      if (index == it->first.first && it->first.second <= numModules) dsDistances[it->first.second-1] = it->second; // it->first is the coordinate pair (layer, ring) for barrel or (ring, disk) for endcap. numbering starts from 1
    }
  }
  return dsDistances;
}

void Tracker::changeRingModules(std::string diskName, int ringN, std::string newType, Color_t newColor) {
  ModuleVector::iterator modIt;
  ModuleVector* moduleV;
  LayerVector::iterator itLay;
  std::string aTag;
  std::ostringstream myTag;

  for (itLay=endcapLayerSet_.begin(); itLay!=endcapLayerSet_.end(); itLay++) {
    if ((*itLay)->getName()==diskName) {
      moduleV = (*itLay)->getModuleVector();
      for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
        aTag = (*modIt)->getTag();
        myTag.str("");
        myTag<<"R"<<ringN;
        if (myTag.str()==aTag) {
          (*modIt)->setType(newType);
          (*modIt)->setModuleType(&(mapType_[newType]));
          (*modIt)->setColor(newColor);
        }
      }
    }
  }
}

void Tracker::drawSummary(double maxZ, double maxRho, std::string fileName) {
  TCanvas* summaryCanvas;
  TVirtualPad* myPad;
  Int_t irep;

  TCanvas* YZCanvas = new TCanvas("YZCanvas", "YZView Canvas", 800, 800 );
  summaryCanvas = new TCanvas("summaryCanvas", "Summary Canvas", 800, 800);
  summaryCanvas->SetFillColor(COLOR_BACKGROUND);
  summaryCanvas->Divide(2, 2);
  summaryCanvas->SetWindowSize(800, 800);

  for (int i=1; i<=4; i++) { myPad=summaryCanvas->GetPad(i); myPad->SetFillColor(kWhite);  }

  // First pad
  // YZView
  myPad = summaryCanvas->GetPad(1);
  myPad->SetFillColor(COLOR_PLOT_BACKGROUND);
  myPad->cd();
  if (geomLiteYZ_) {
    drawGrid(maxZ, maxRho, ViewSectionYZ);
    geomLiteYZ_->DrawClonePad();
    myPad->GetView()->SetParallel();
    myPad->GetView()->SetRange(0, 0, 0, maxZ, maxZ, maxZ);
    myPad->GetView()->SetView(0 /*long*/, 270/*lat*/, 270/*psi*/, irep);
    drawTicks(myPad->GetView(), maxZ, maxRho, ViewSectionYZ);

    YZCanvas->cd();
    myPad = YZCanvas->GetPad(0);
    myPad->SetFillColor(COLOR_PLOT_BACKGROUND);
    drawGrid(maxZ, maxRho, ViewSectionYZ);
    geomLiteYZ_->DrawClonePad();
    myPad->GetView()->SetParallel();
    myPad->GetView()->SetRange(0, 0, 0, maxZ, maxZ, maxZ);
    myPad->GetView()->SetView(0 /*long*/, 270/*lat*/, 270/*psi*/, irep);
    drawTicks(myPad->GetView(), maxZ, maxRho, ViewSectionYZ);
  }

  // First pad
  // XYView (barrel)
  myPad = summaryCanvas->GetPad(2);
  myPad->cd();
  myPad->SetFillColor(COLOR_PLOT_BACKGROUND);
  if (geomLiteXY_) {
    drawGrid(maxZ, maxRho, ViewSectionXY);
    geomLiteXY_->DrawClonePad();
    myPad->GetView()->SetParallel();
    myPad->GetView()->SetRange(-maxRho, -maxRho, -maxRho, maxRho, maxRho, maxRho);
    myPad->GetView()->SetView(0 /*long*/, 0/*lat*/, 270/*psi*/, irep);
  }

  // Third pad
  // Plots
  myPad = summaryCanvas->GetPad(3);
  myPad->cd();
  myPad->SetFillColor(COLOR_PLOT_BACKGROUND);
  if (etaProfileCanvas_) {
    etaProfileCanvas_->DrawClonePad();
  }

  // Fourth pad
  // XYView (EndCap)
  myPad = summaryCanvas->GetPad(4);
  myPad->cd();
  myPad->SetFillColor(COLOR_PLOT_BACKGROUND);
  if (geomLiteEC_) {
    drawGrid(maxZ, maxRho, ViewSectionXY);
    geomLiteEC_->DrawClonePad();
    myPad->GetView()->SetParallel();
    myPad->GetView()->SetRange(-maxRho, -maxRho, -maxRho, maxRho, maxRho, maxRho);
    myPad->GetView()->SetView(0 /*long*/, 0/*lat*/, 270/*psi*/, irep);
  }


  for (int i=1; i<=4; i++) { myPad=summaryCanvas->GetPad(i); myPad->SetBorderMode(0); }

  summaryCanvas->Modified();

  std::string pngFileName = fileName+".png";
  std::string YZpngFileName = fileName+"YZ.png";
  std::string svgFileName = fileName+".svg";
  std::string YZsvgFileName = fileName+"YZ.svg";
  //std::string gifFileName = fileName+".gif";

  summaryCanvas->SaveAs(pngFileName.c_str());
  summaryCanvas->SaveAs(svgFileName.c_str());
  YZCanvas->SaveAs(YZpngFileName.c_str());
  YZCanvas->SaveAs(YZsvgFileName.c_str());


  if (mapPhiEta_) {
    int prevStat = gStyle->GetOptStat();
    gStyle->SetOptStat(0);
    TCanvas* hitMapCanvas = new TCanvas("hitmapcanvas", "Hit Map", 800, 800);
    hitMapCanvas->cd();
    gStyle->SetPalette(1);
    pngFileName = fileName+"_hitmap.png";                                                                                                                                  
    svgFileName = fileName+"_hitmap.svg";
    hitMapCanvas->SetFillColor(COLOR_PLOT_BACKGROUND);
    hitMapCanvas->SetBorderMode(0);
    hitMapCanvas->SetBorderSize(0);
    mapPhiEta_->Draw("colz");
    hitMapCanvas->Modified();
    hitMapCanvas->SaveAs(pngFileName.c_str());
    hitMapCanvas->SaveAs(svgFileName.c_str());
    gStyle->SetOptStat(prevStat);
  }

  if (etaProfileCanvas_) {
    summaryCanvas = new TCanvas("etaprofilebig", "big etaprofile plot", 1000, 700);
    summaryCanvas->cd();
    pngFileName = fileName+"_nhitplot.png";
    svgFileName = fileName+"_nhitplot.svg";
    etaProfileCanvas_->DrawClonePad();
    etaProfileCanvas_->SetFillColor(COLOR_PLOT_BACKGROUND);
    //     TFrame* myFrame = summaryCanvas->GetFrame();
    //     if (myFrame) {
    //       myFrame->SetFillColor(kYellow-10);
    //     } else {
    //       std::cerr << "myFrame is NULL" << std::endl;
    //     }
    summaryCanvas->SetFillColor(COLOR_BACKGROUND);
    summaryCanvas->SetBorderMode(0);
    summaryCanvas->SetBorderSize(0);
    summaryCanvas->Modified();
    summaryCanvas->SaveAs(pngFileName.c_str());
    summaryCanvas->SaveAs(svgFileName.c_str());
  }


  if (bandWidthCanvas_) {
    pngFileName = fileName+"_bandwidth.png";
    svgFileName = fileName+"_bandwidth.svg";
    bandWidthCanvas_->DrawClonePad();
    bandWidthCanvas_->SaveAs(pngFileName.c_str());
    bandWidthCanvas_->SaveAs(svgFileName.c_str());
  }

  //summaryCanvas->SaveAs(epsFileName.c_str());
  //summaryCanvas->SaveAs(gifFileName.c_str());
}

void Tracker::drawLayout(double maxZ, double maxRho, std::string fileName) {
  TCanvas* layoutCanvas;
  Int_t irep;

  layoutCanvas = new TCanvas("layoutCanvas", "Layout Canvas", 400, 400);
  layoutCanvas->SetFillColor(kWhite);
  layoutCanvas->SetWindowSize(400, 400);

  // YZView only is our layout canvas
  layoutCanvas->cd();
  if (geomLiteYZ_) {
    drawGrid(maxZ, maxRho, ViewSectionYZ);
    geomLiteYZ_->DrawClonePad();
    layoutCanvas->GetView()->SetParallel();
    layoutCanvas->GetView()->SetRange(0, 0, 0, maxZ, maxZ, maxZ);
    layoutCanvas->GetView()->SetView(0 /*long*/, 270/*lat*/, 270/*psi*/, irep);
    drawTicks(layoutCanvas->GetView(), maxZ, maxRho, ViewSectionYZ);
  }

  layoutCanvas->SetBorderMode(0);
  layoutCanvas->Modified();

  layoutCanvas->SaveAs(fileName.c_str());
}

double Tracker::getSmallDelta(const int& index) {
  if (specialSmallDelta_[index]==0) {
    return smallDelta_;
  } else {
    return specialSmallDelta_[index];
  }
}

double Tracker::getBigDelta(const int& index) {
  if (specialBigDelta_[index]==0) {
    return bigDelta_;
  } else {
    return specialBigDelta_[index];
  }
}


// private
/**
 * Creates a module type map
 * It sets a different integer for each one
 * @param tracker the tracker to be analyzed
 * @param moduleTypeCount the map to count the different module types
 * @return the total number of module types
 */
int Tracker::createResetCounters(std::map <std::string, int> &moduleTypeCount) {
  ModuleVector result;
  LayerVector::iterator layIt;
  ModuleVector* moduleV;
  ModuleVector::iterator modIt;

  std::string aType;
  int typeCounter=0;

  for (layIt=layerSet_.begin(); layIt!=layerSet_.end(); layIt++) {
    moduleV = (*layIt)->getModuleVector();
    for (modIt=moduleV->begin(); modIt!=moduleV->end(); modIt++) {
      aType = (*modIt)->getType();
      (*modIt)->resetNHits();
      if (moduleTypeCount.find(aType)==moduleTypeCount.end()) {
        moduleTypeCount[aType]=typeCounter++;
      }
    }
  }

  return(typeCounter);
}


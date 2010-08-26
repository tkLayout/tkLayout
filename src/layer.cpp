#include <iostream>
#include <sstream>
#include <cmath>
#include "layer.hh"
#include "Math/RotationZ.h"
#include "Math/Vector3D.h"

using namespace ROOT::Math;

/******************/
/*                */
/* Generic module */
/*                */
/******************/

Layer::~Layer() {
    ModuleVector::iterator modIt;
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        if ((*modIt)!=NULL) {
            delete (*modIt);
        }
    }
    moduleSet_.clear();
}

Layer::Layer() : MessageLogger("Layer") {
    setDefaultParameters();
}

void Layer::setDefaultParameters() {
    layerName_ = "Layer";
}

void Layer::translate(XYZVector Delta) {
    ModuleVector::iterator modIt;
    
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        (*modIt)->translate(Delta);
    }
}


// TODO: tidy up the next two functions: the Endcap-specific onw should
// just call the generic Layer::rotateY_PI(), and multiply the averageZ_ by -1
void BarrelLayer::rotateY_PI() {
    ModuleVector::iterator modIt;
    
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        (*modIt)->rotateY_PI();
    }
}

void BarrelLayer::reflectZ() {
    ModuleVector::iterator modIt;
    
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        (*modIt)->reflectZ();
    }
}

void BarrelLayer::shiftRho(double Delta) {
    ModuleVector::iterator modIt;
    
    averageRadius_ = InvalidRadius;
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        (*modIt)->shiftRho(Delta);
    }
}

void EndcapLayer::rotateY_PI() {
    ModuleVector::iterator modIt;
    
    averageZ_ *= -1;
    
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        (*modIt)->rotateY_PI();
    }
}


void Layer::shapeVolume(TGeoVolume* container,
        TGeoMedium* medium,
        TGeoManager* geom) {
    // TODO:
    // Should build only the containing volume
}

void Layer::shapeModuleVolumes(TGeoVolume* container,
        TGeoMedium* medium,
        TGeoManager* geom) {
    
    // TODO:
    // Should build only the containing volume
    
}

double Layer::getMaxZ() {
    double maxZ;
    
    ModuleVector::iterator modIt;
    maxZ=(*moduleSet_.begin())->getMaxZ();
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        if ((*modIt)->getMaxZ()>maxZ) maxZ=(*modIt)->getMaxZ();
    }
    
    return maxZ;
};

double Layer::getMinZ() {
    double minZ;
    
    ModuleVector::iterator modIt;
    minZ=(*moduleSet_.begin())->getMinZ();
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        if ((*modIt)->getMinZ()<minZ) minZ=(*modIt)->getMinZ();
    }
    
    return minZ;
};

double Layer::getMaxRho() {
    double maxRho;
    
    ModuleVector::iterator modIt;
    maxRho=(*moduleSet_.begin())->getMaxRho();
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        if ((*modIt)->getMaxRho()>maxRho) maxRho=(*modIt)->getMaxRho();
    }
    
    return maxRho;
};

double Layer::getMinRho() {
    double minRho;
    
    ModuleVector::iterator modIt;
    minRho=(*moduleSet_.begin())->getMinRho();
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        if ((*modIt)->getMinRho()<minRho) minRho=(*modIt)->getMinRho();
    }
    
    return minRho;
};

/******************/
/*                */
/* Barrel layer   */
/*                */
/******************/

BarrelLayer::~BarrelLayer() {
    if (sampleModule_) delete sampleModule_;
}

BarrelLayer::BarrelLayer() {
    sampleModule_ = NULL;
}

BarrelLayer::BarrelLayer(BarrelLayer& inputLayer) {
    ModuleVector::iterator modIt;
    ModuleVector inputModuleV = inputLayer.moduleSet_;
    BarrelModule* aModule;
    BarrelModule* stdModule;
    
    setName(inputLayer.layerName_);
    averageRadius_ = inputLayer.averageRadius_;
    nOfRods_ = inputLayer.getRods();
    nModsOnString_ = inputLayer.getModulesOnRod();
    
    if ( (sampleModule_=dynamic_cast<BarrelModule*>(inputLayer.getSampleModule())) ) {
        sampleModule_ = new BarrelModule(*(inputLayer.getSampleModule()));
    }
    
    for (modIt=inputModuleV.begin(); modIt!=inputModuleV.end(); modIt++) {
        if ( (stdModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
            aModule = new BarrelModule(*stdModule);
            moduleSet_.push_back(aModule);
        } else {
	  addMessage("In BarrelLayer::BarrelLayer(BarrelLayer& inputLayer) I found a non-barrel module in the source", ERROR);
        }
    }
}

BarrelLayer::BarrelLayer(double waferDiameter, double heightOverWidth) {
    sampleModule_ = new BarrelModule(waferDiameter, heightOverWidth);
}

BarrelLayer::BarrelLayer(double heightOverWidth ) {
    sampleModule_ = new BarrelModule(heightOverWidth);
}

BarrelLayer::BarrelLayer(const BarrelModule& mySample) {
    sampleModule_ = new BarrelModule(mySample);
}

BarrelLayer::BarrelLayer(BarrelModule* mySample) {
    sampleModule_ = new BarrelModule(*mySample);
}


// Builds a string and places into a module vector
// Obsolete function to be removed)
int BarrelLayer::buildString(ModuleVector& thisModuleSet,
        double stringAverageRadius,
        double smallDelta, // Half the distance between inner and outer modules
        // In smalldelta it is included the string's parity
        double zOverlap,
        double safetyOrigin,
        double maxZ,
        BarrelModule* sampleModule) {
    
    addMessage("buildString() with maxZ option is an obsolete function: no more mantained, you should not be using this function", ERROR);
    
    int parity;
    int nModules;
    
    // Parity of modules in the string
    // for n=1 smallParity=1
    // for n=2 smallParity=-1
    int smallParity;
    edge myEdge;
    
    if (maxZ>0) {
        parity = +1;
    } else {
        parity = -1;
    }
    
    BarrelModule* aModule;
    BarrelModule* trialModule[3];
    
    aModule = new BarrelModule(*sampleModule);
    
    aModule->translate(XYZVector(0, stringAverageRadius-(parity)*smallDelta, 0));
    
    // The first module is half displeced because the module from the opposite
    // string will also be displaced.
    aModule->setEdgeZ(-1*parity*zOverlap/2., parity);
    thisModuleSet.push_back(aModule);
    
    // Number of modules already in this string
    nModules=1;
    
    edge minSafetyEdge;
    
    while( fabs(aModule->getEdgeZSide(parity, 0).first)<fabs(maxZ) ) {
        // I already take into account the new module
        nModules++;
        
        // The smallParity parameter sets the direction of
        // the smallDelta (rho difference of modules within the same string)
        // smallParity will be -1 the first time here
        // [nModules=2]
        smallParity=(nModules%2)*2-1;
        
        // Two conditions may apply: either we ask for an overlap
        // with respect to the projection from the origin (CASE 0)
        // or we ask for the exact projection from a displaced origin (CASE 1,2)
        
        // CASE 0
        trialModule[0] = new BarrelModule(*aModule);             // Clone the previous
        myEdge = trialModule[0]->getEdgeZSide(parity, zOverlap); // Measure the overlapping point at the far end
        trialModule[0]->setEdgeZ(myEdge.first, parity);          // Move the near end there
        myEdge = trialModule[0]->getEdgeZSide(-1*parity, 0.);     // Pick the near end
        trialModule[0]->projectSideRho(myEdge.second,            // Project the module with no origin displacement
                stringAverageRadius-(smallParity)*(parity)*smallDelta);
        // ... and let's write down were we put it...
        minSafetyEdge=trialModule[0]->getEdgeZSide(-1*parity, 0.);
        // and its index
        minSafetyEdge.second=0;
        
        // CASE 1 and 2
        trialModule[1] = new BarrelModule(*aModule);          // Clone the previous
        myEdge = trialModule[1]->getEdgeZSide(parity, 0.);    // Measure the far end with no overlap
        trialModule[1]->setEdgeZ(myEdge.first, parity);       // Move the near end there
        myEdge = trialModule[0]->getEdgeZSide(-1*parity, 0.);  // Pick the near end
        trialModule[2] = new BarrelModule(*trialModule[1]);   // Clone the module
        // Now we only need to project it twice, with two different origin
        // displacements: in fact we project the module to its destination twice
        // using +- safetyOrigin safety margin along the beam axis
        for (int i=1; i<3; i++) {
            int safetyDirection=i*2-3;
            trialModule[i]->projectSideRho(myEdge.second,
                    stringAverageRadius-(smallParity)*(parity)*smallDelta,
                    safetyOrigin*safetyDirection);
            myEdge = trialModule[i]->getEdgeZSide(-1*parity, 0.);
            
            // If the measured margin is the safest...
            if ((parity*myEdge.first)<(parity*minSafetyEdge.first)) {
                // ..then we record it
                minSafetyEdge.first=myEdge.first;
                minSafetyEdge.second=i;
            }
        }
        
        // Debug output
        // std::cerr << "Safety index: " << minSafetyEdge.second << std::endl; // debug
        
        // Now we just keep the safest, while we delete the others
        for (int i=0; i<3; i++) {
            if (minSafetyEdge.second==i) {
                aModule=trialModule[i];
            } else {
                delete trialModule[i];
            }
        }
        
        // Finally we add this module in the list of built
        thisModuleSet.push_back(aModule);
    }
    
    return nModules;
    
}

// Builds a string and places into a module vector
int BarrelLayer::buildString(ModuleVector& thisModuleSet,
        double stringAverageRadius,
        double smallDelta, // Half the distance between inner and outer modules
        // In smalldelta it is included the string's parity
        double zOverlap,
        double safetyOrigin,
        int nDesiredModules,
        BarrelModule* sampleModule,
        double minZ /* = 0 used to build mezzanine layers */ ) {
    
    int parity;
    int nModules;
    
    // Parity of modules in the string
    // for n=1 smallParity=1
    // for n=2 smallParity=-1
    int smallParity;
    edge myEdge;
    
    if (nDesiredModules>0) {
        parity = +1;
    } else {
        parity = -1;
    }
    
    nDesiredModules /= parity;
    
    BarrelModule* aModule;
    BarrelModule* trialModule[3];
    
    aModule = new BarrelModule(*sampleModule);
    
    aModule->translate(XYZVector(0, stringAverageRadius-(parity)*smallDelta, 0));
    
    // The first module is half displeced because the module from the opposite
    // string will also be displaced.
    if (minZ==0) { // A string of a standard barrel
        aModule->setEdgeZ(-1*parity*zOverlap/2., parity);
    } else { // A special string starting at minZ (mezzanine barrel)
        aModule->setEdgeZ(parity*minZ, parity);
    }
    // std::cout << "pushing back one module out of " << nDesiredModules << " desired modules" << std::endl; // debug
    aModule->setRing(1);
    thisModuleSet.push_back(aModule);
    
    edge minSafetyEdge;
    
    // I start counting from one, as I already created module #0
    for (nModules = 1; nModules<nDesiredModules; nModules++) {
        
        // The smallParity parameter sets the direction of
        // the smallDelta (rho difference of modules within the same string)
        // smallParity will be -1 the first time here
        // [nModules=2]
        smallParity=1-(nModules%2)*2;
        
        // Two conditions may apply: either we ask for an overlap
        // with respect to the projection from the origin (CASE 0)
        // or we ask for the exact projection from a displaced origin (CASE 1,2)
        
        // CASE 0
        trialModule[0] = new BarrelModule(*aModule);             // Clone the previous
        myEdge = trialModule[0]->getEdgeZSide(parity, zOverlap); // Measure the overlapping point at the far end
        trialModule[0]->setEdgeZ(myEdge.first, parity);          // Move the near end there
        myEdge = trialModule[0]->getEdgeZSide(-1*parity, 0.);     // Pick the near end
        trialModule[0]->projectSideRho(myEdge.second,            // Project the module with no origin displacement
                stringAverageRadius-(smallParity)*(parity)*smallDelta);
        // ... and let's write down were we put it...
        minSafetyEdge=trialModule[0]->getEdgeZSide(-1*parity, 0.);
        // and its index
        minSafetyEdge.second=0;
        
        // CASE 1 and 2
        trialModule[1] = new BarrelModule(*aModule);          // Clone the previous
        myEdge = trialModule[1]->getEdgeZSide(parity, 0.);    // Measure the far end with no overlap
        trialModule[1]->setEdgeZ(myEdge.first, parity);       // Move the near end there
        myEdge = trialModule[0]->getEdgeZSide(-1*parity, 0.);  // Pick the near end
        trialModule[2] = new BarrelModule(*trialModule[1]);   // Clone the module
        // Now we only need to project it twice, with two different origin
        // displacements: in fact we project the module to its destination twice
        // using +- safetyOrigin safety margin along the beam axis
        for (int i=1; i<3; i++) {
            int safetyDirection=i*2-3;
            trialModule[i]->projectSideRho(myEdge.second,
                    stringAverageRadius-(smallParity)*(parity)*smallDelta,
                    safetyOrigin*safetyDirection);
            myEdge = trialModule[i]->getEdgeZSide(-1*parity, 0.);
            
            // If the measured margin is the safest...
            if ((parity*myEdge.first)<(parity*minSafetyEdge.first)) {
                // ..then we record it
                minSafetyEdge.first=myEdge.first;
                minSafetyEdge.second=i;
            }
        }
        
        // Debug output
        // std::cerr << "Safety index: " << minSafetyEdge.second << std::endl; // debug
        
        // Now we just keep the safest, while we delete the others
        for (int i=0; i<3; i++) {
            if (minSafetyEdge.second==i) {
                aModule=trialModule[i];
            } else {
                delete trialModule[i];
            }
        }
        
        // std::cout << "pushing back module "<< nModules <<" out of " << nDesiredModules << " desired modules" << std::endl; // debug
        // Finally we add this module in the list of built
        aModule->setRing(nModules+1);
        thisModuleSet.push_back(aModule);
    }
    
    return nModules;
    
}


void BarrelLayer::buildStringPair(ModuleVector& thisModuleSet,
        double stringAverageRadius,
        double smallDelta, // Half the distance between inner and outer modules
        double zOverlap,
        double safetyOrigin,
        double maxZ,
        BarrelModule* sampleModule) {
    int n1; int n2;
    
    n1 = buildString(thisModuleSet, stringAverageRadius, smallDelta, zOverlap, safetyOrigin, maxZ, sampleModule);
    n2 = buildString(thisModuleSet, stringAverageRadius, smallDelta, zOverlap, safetyOrigin, -1*maxZ, sampleModule);
    
    if (n1!=n2) addMessage("Building an asymmetric layer. This should not happen", ERROR);
}

void BarrelLayer::buildStringPair(ModuleVector& thisModuleSet,
        double stringAverageRadius,
        double smallDelta, // Half the distance between inner and outer modules
        double zOverlap,
        double safetyOrigin,
        int nModules,
        BarrelModule* sampleModule) {
    int n1; int n2;
    
    n1 = buildString(thisModuleSet, stringAverageRadius, smallDelta, zOverlap, safetyOrigin, nModules, sampleModule);
    n2 = buildString(thisModuleSet, stringAverageRadius, smallDelta, zOverlap, safetyOrigin, -1*nModules, sampleModule);
    
    if (n1!=n2) addMessage("Building an asymmetric layer. This should not happen", ERROR);
}

// Computes the optimal radius of the inner of two
// layers, which should cover each other
double BarrelLayer::layerRadius(const double& nMod,  // number of strings in a layer
        const double& g,     // Gap between inner and outer
        const double& o,     // Needed overlap in mm
        const double& b      // Module's half width
) {
    double result;
    
    double f = b - o;
    
    double t = tan(2*M_PI/nMod);
    double B = t*g - f - b;
    double C = -1 * (t*f*b+f*g);
    
    // The solution is always the one with plus
    result = -1*B+sqrt(pow(B, 2)-4*t*C);
    result /= 2*t;
    
    return result;
}


// Computes the optimal radius of the inner of two
// layers, which should cover each other, using the desired
// optimization
std::pair<double, int> BarrelLayer::computeRadius(const double& x,    // radius of the inner module
        const double& g,    // Gap between inner and outer
        const double& o,    // Needed overlap in mm
        const double& b,    // Module's half width
        const int& optimal, // wether to shrink or enlarge, fix or auto
        const int& base )   // Base number of modules
{
    
    std::pair<double, int> result;
    
    // Module's half width minus the desired overlap
    double f = b - o;
    
    // (modules=2*nHalf)
    //  double nHalf = 2*M_PI/(2*atan((f*(x+g)+b*x)/(((x+g)*x)-f*b)));
    
    // Tentative number of modules on a base
    // (2-> half shell)
    // (4-> quarter shell)
    double nBase = 2*M_PI/(base*atan((f*(x+g)+b*x)/(((x+g)*x)-f*b)));
    
    // Here the behaviour depends on the pushDirection
    switch (optimal) {
        case SHRINK:
            nBase = floor(nBase);
            result.first = layerRadius(base*nBase, g, o, b);
            break;
        case ENLARGE:
            nBase = ceil(nBase);
            result.first = layerRadius(base*nBase, g, o, b);
            break;
        case FIXED:
            result.first=x;
            nBase = ceil(nBase);
            break;
        case AUTO:
            double nLo, nHi;
            double rLo, rHi;
            
            nLo = floor(nBase);
            nHi = ceil(nBase);
            rLo = layerRadius(base*nLo, g, o, b);
            rHi = layerRadius(base*nHi, g, o, b);
            
            if (fabs(rHi-x)<fabs(rLo-x)) {
                result.first = rHi;
                nBase = nHi;
            } else {
                result.first = rLo;
                nBase = nLo;
            }
            break;
    }
    
    result.second=int(nBase)*base;
    
    return result;
    
}


// Returns the optimal average radius
// along with the total number of strings in a layer
// It only requires the number of strings in a layer to be
// multiple of 2
std::pair<double, int> BarrelLayer::layerPhi(double avRad,          // Desired layer average radius
        double smallDelta,     // Half distance between modules in the same string
        double bigDelta,       // String half gap
        double o,              // overlap
        double b,              // Module half width
        int optim,             // wether to shrink or enlarge, fix or auto
        int base,              // base number of modules (nMod=base*m)
        bool stringSameParity  // Wether the strings should have
// the same or opposite module parity (in/out)
) {
    
    // I chose OppdInner for opposite-parity strings because it's the worse condition
    // (of the smartest configuration).
    
    std::pair<double, int> result;
    
    if (stringSameParity) {
        // SameOuter; (it's always like this)
        result = computeRadius(avRad-bigDelta+smallDelta, 2*bigDelta, o, b, optim, base);
        result.first -= -bigDelta+smallDelta;
        
        // SameInner;
        // result = computeRadius(avRad-bigDelta-smallDelta, 2*bigDelta, o, b, optim, base);
        // result.first -= bigDelta+smallDelta;
    } else {
        // Here the worst condition depends on the parameters
        // OppdInner;
        result = computeRadius(avRad-bigDelta+smallDelta, 2*(bigDelta-smallDelta), o, b, optim, base);
        result.first -= -bigDelta+smallDelta;
        
        // OppdOuter;
        // result = computeRadius(avRad-bigDelta-smallDelta, 2*(bigDelta-smallDelta), o, b, optim, base);
        // result.first -= bigDelta+smallDelta;
    }
    
    return result;
}


void BarrelLayer::buildLayer(double averageRadius,
        double smallDelta,
        double bigDelta,
        double overlap,
        double safetyOrigin,
        int nModules,
        int pushDirection,      // Direction where the average
        // radius is pushed to optimize module usage
        // +1: towards outside, 0 no optimize, -1 towards inside
        int base,               // Modules base number
        bool stringSameParity,  // Wether the strings should have
        // the same or opposite module parity (in/out)
        BarrelModule* sampleModule,
        int sectioned, // = NoSection  (deprecated)
        double minZ // = 0 not zero in case you build a mezzanine
) {
    
    int stringParity;      // Parity of string in the shell (inner, outer)
    int stringSmallParity; // Parity of module in the string
    std::pair <double, int> optimalBarrel;
    double goodRadius;
    int nStrings;
    
    nOfRods_ = 0;
    nModsOnString_ = 0;
    optimalBarrel = layerPhi(averageRadius   , // tentativeX
            smallDelta,
            bigDelta,
            overlap,
            sampleModule->getWidth()/2. , // Module half width
            pushDirection, // wether to shrink or enlarge (>0 means enlarge)
            base,
            stringSameParity
            );
    
    goodRadius = optimalBarrel.first;
    nStrings = optimalBarrel.second;
    averageRadius_ = goodRadius;
    
    tempString.str(""); tempString << "GoodRadius: " << goodRadius
				   << ", nStrings:   " << nStrings;
    addMessage(tempString, INFO);
    
    
    if (nStrings%2!=0) {
      addMessage("WARNING: you just asked for a layer with a number of strings not multiple of 2", WARNING);
    }
    
    double stringPhiShift = 2*M_PI/double(nStrings);
    double rightAngle;
    std::vector<Module*>::iterator itMod;
    
    for (int i=0; i<nStrings; i++) {
        stringParity = (i%2)*2-1;
        if (stringSameParity) {
            stringSmallParity = 1;
        } else {
            stringSmallParity = -1 * stringParity;
        }
        // std::cout << "Building a string" << std::endl; // debug
        
        std::vector<Module*> aString;
        if (minZ==0) { // Build a standard (double) string
            buildStringPair(aString,
                    goodRadius+stringParity*bigDelta,
                    stringSmallParity*smallDelta,
                    overlap,
                    safetyOrigin,
                    nModules,
                    sampleModule);
        } else { // Build a mezzanine (single) string
            buildString(aString,
                    goodRadius+stringParity*bigDelta,
                    stringSmallParity*smallDelta,
                    overlap,
                    safetyOrigin,
                    nModules,
                    sampleModule,	minZ);
        }
        
        if (i==0) {
            //std::cerr << "Computing min edge: "; // debug
            edge phiEdge;
            BarrelModule* stdBarrelMod;
            for (itMod=aString.begin(); itMod!=aString.end(); itMod++) {
                if ( (stdBarrelMod=dynamic_cast<BarrelModule*>(*itMod)) ) {
                    phiEdge=((BarrelModule*)(*itMod))->getEdgePhiSide(-1);
                    if (itMod==aString.begin()) {
                        rightAngle=phiEdge.first;
                        //std::cerr << "first angle is " << rightAngle << " and then: "; // debug
                    } else {
                        //std::cerr << phiEdge.first << " ";// debug
                        if (phiEdge.first<rightAngle) rightAngle=phiEdge.first;
                    }
                } else {
		  // This should never happen
		  addMessage("In building a string: found a non-barrel module", ERROR);
                }
            }
            //std::cerr << "done. The best is: " << rightAngle << std::endl; // debug
        }
        
        bool firstOnes;
        int aSection;
        std::vector<Module* >::iterator secondMod;
        secondMod=aString.begin()++;
        
        for (itMod=aString.begin(); itMod!=aString.end(); itMod++) {
            aSection=NoSection;
            if ((itMod==aString.begin())||(itMod==secondMod)) {
                firstOnes=true;
            } else {
                firstOnes=false;
            }
            
            // TODO: find a way to conjugate the need for a simple geometry
            // for a fast simulation (first string up) and a real geometry with modules
            // not crossing the XZ plane
            //(*itMod)->rotatePhi(i*stringPhiShift-rightAngle);
            (*itMod)->rotatePhi(i*stringPhiShift);
            
            if (i==0) {
                aSection |= YZSection;
            }
            if (firstOnes) {
                aSection |= XYSection;
                // 	std::cout << "it's one of the first: it belongs to section" << XYSection << std::endl; // debug
            }
            
            switch (sectioned) {
                case NoSection:
                    moduleSet_.push_back(*itMod);
                    break;
                case XYSection:
                    if (firstOnes) {
                        moduleSet_.push_back(*itMod);
                    }
                    break;
                case YZSection:
                    if (i==0) {
                        //(*itMod)->rotatePhi(rightAngle);
                        moduleSet_.push_back(*itMod);
                    }
                    break;
            }
            
            (*itMod)->setSection(aSection);
            
            //  if (firstOnes)
            // 	std::cout << "its section is " << (*itMod)->getSection() << std::endl; //debug
        }
    }
    nOfRods_ = nStrings;
    nModsOnString_ = nModules;
}

// Always look for this plot when changing the geometry!
// Execute the program with 2> aaa
// and then "sort -u aaa"
void BarrelLayer::neededModulesPlot(double smallDelta, // Half distance between modules in the same string
        double bigDelta,   // String half gap
        double o,          // overlap
        double b,          // Module half width
        int base
) {
    double minR=400;
    double maxR=1100;
    int steps=50000;
    
    
    //   TProfile* nmod[2];
    
    //   char name[100];
    //   sprintf(name, "NeededModulesParitySame");
    //   nmod[0] = new TProfile(name, name, steps+1, minR, maxR);
    //   nmod[0]->SetLineColor(1);
    
    //   sprintf(name, "NeededModulesParityOppd");
    //   nmod[1] = new TProfile(name, name, steps+1, minR, maxR);
    //   nmod[1]->SetLineColor(2);
    
    std::pair<double, int> SameOuter, SameInner, OppdOuter, OppdInner;
    
    
    int maxnSame, maxnOppd;
    
    for (double x=minR; x<maxR; x+=(maxR-minR)/double(steps)) {
        SameOuter = computeRadius(x-bigDelta+smallDelta, 2*bigDelta, o, b, FIXED, base);
        SameInner = computeRadius(x-bigDelta-smallDelta, 2*bigDelta, o, b, FIXED, base);
        OppdInner = computeRadius(x-bigDelta+smallDelta, 2*(bigDelta-smallDelta), o, b, FIXED, base);
        OppdOuter = computeRadius(x-bigDelta-smallDelta, 2*(bigDelta-smallDelta), o, b, FIXED, base);
        
        maxnSame = (SameOuter.second > SameInner.second) ? SameOuter.second : SameInner.second ;
        maxnOppd = (OppdOuter.second > OppdInner.second) ? OppdOuter.second : OppdInner.second ;
        
        if (SameOuter.second > SameInner.second) {
	  addMessage("SameOuter is bigger than SameInner", DEBUG);
        }
        if (SameOuter.second < SameInner.second) {
	  addMessage("SameInner is bigger than SameOuter", DEBUG);
        }
        
        if (OppdOuter.second > OppdInner.second) {
	  addMessage("OppdOuter is bigger than OppdInner", DEBUG);
        }
        if (OppdOuter.second < OppdInner.second) {
	  addMessage("OppdInner is bigger than OppdOuter", DEBUG);
        }
        
        if (maxnSame>maxnOppd) addMessage("Same is bigger than Oppd", DEBUG);
        if (maxnSame<maxnOppd) addMessage("Oppd is bigger than Same", DEBUG);
    }
    
    return;
}

int BarrelLayer::cutOverEta(double etaCut) {
    int nCut = 0;
    ModuleVector::iterator modIt;
    
    double theta;
    double eta;
    
    //double zmax = Layer::getMaxZ();
    //double zmin = getMinZ();
    
    int minRingDelete=0;
    bool firstCutFound=false;
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        theta=(*modIt)->getMeanTheta();
        eta = -1*log(tan(theta/2.));
        if (fabs(eta)>etaCut) {
            if (!firstCutFound) {
                minRingDelete=(*modIt)->getRing();
                firstCutFound=true;
            } else {
                if ((*modIt)->getRing() < minRingDelete) minRingDelete = (*modIt)->getRing();
            }
        }
    }
    
    if (firstCutFound) {
        // Bookkeeping
        nModsOnString_ = minRingDelete -1;
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); ) {
        if ((*modIt)->getRing()>=minRingDelete) {
            // Erase the useless module!
            delete (*modIt);
            moduleSet_.erase(modIt);
            nCut++;
        } else {
            modIt++;
        }
    }
    }
    return nCut;
}

double BarrelLayer::getMaxZ(int direction) {
    double maxZ;
    //double minZ;
    double aZ;
    ModuleVector::iterator modIt;
    BarrelModule* aBarrelModule;
    bool firstModule=true;
    
    if (direction==0) {
      addMessage("BarrelLayer::getMaxZ was called with direction == 0", ERROR);
      return 0;
    }
    direction/=int(fabs(direction));
    
    maxZ=0;
    
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
            aZ = aBarrelModule->getEdgeZSide(direction).first;
            if (firstModule) {
                firstModule=false;
                maxZ=aZ;
            }
            if (direction*aZ>direction*maxZ) {
                maxZ=aZ;
            }
        }
    }
    
    return maxZ;
}



void BarrelLayer::compressToZ(double newMaxZ) {
    double maxZ;
    double minZ;
    double aZ;
    ModuleVector::iterator modIt;
    BarrelModule* aBarrelModule;
    BarrelModule* maxBarrelModule=NULL;
    BarrelModule* minBarrelModule=NULL;
    bool firstModule=true;
    
    maxZ=0;
    minZ=0;
    
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
            aZ = aBarrelModule->getEdgeZSide(+1).first;
            if ((aZ>maxZ)||firstModule) {
                maxBarrelModule = aBarrelModule;
                maxZ=aZ;
            }
            aZ = aBarrelModule->getEdgeZSide(-1).first;
            if ((aZ<minZ)||firstModule) {
                minBarrelModule = aBarrelModule;
                minZ=aZ;
            }
            if (firstModule) firstModule=false;
        }
    }
    
    maxZ=fabs(maxZ);
    minZ=fabs(minZ);
    
    //   std::cout << "The layer's highest z+ is " << maxZ << std::endl; //debug
    //   std::cout << "The layer's highest z- is " << minZ << std::endl; //debug
    
    // Measure the Delta of the farthest module
    double Deltap;
    double Deltam;
    Deltap = fabs(newMaxZ) - maxZ;
    Deltam = -1*fabs(newMaxZ) + minZ;
    
    // std::cerr << "Delta-: " << Deltam << std::endl; // debug
    // std::cerr << "Delta+: " << Deltap << std::endl; // debug
    // TODO: raise an alarm if Delta > 0
    
    // delta is the ratio between Delta and the module's position
    double deltam=0;
    double deltap=0;
    double myMeanZ;
    XYZVector modShift;
    if (maxBarrelModule!=NULL) {
        deltam=Deltam/((minBarrelModule->getMeanPoint()).Z());
        deltap=Deltap/((maxBarrelModule->getMeanPoint()).Z());
        // std::cerr << "delta-: " << deltam << std::endl; // debug
        // std::cerr << "delta+: " << deltap << std::endl; // debug
        
        for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
            if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
                myMeanZ = ((aBarrelModule->getMeanPoint()).Z());
                if (myMeanZ>0) {
                    modShift=XYZVector(0, 0, deltap*myMeanZ);
                } else {
                    modShift=XYZVector(0, 0, deltam*myMeanZ);
                }
                aBarrelModule->translate(modShift);
            }
        }
    } else {
      addMessage("int BarrelLayer::compactToZ couldn't find the maxZ module", ERROR);
    }
    
    
    // TODO: add exit code
}

// This procedure assumes there is no module with a border < newMinZ
void BarrelLayer::compressExceeding(double newMaxZ, double newMinZ) {
    double maxZ;
    double aZ;
    ModuleVector::iterator modIt;
    BarrelModule* aBarrelModule;
    BarrelModule* maxBarrelModule=NULL;
    bool firstModule=true;
    
    int direction = int(newMaxZ/newMaxZ);
    
    maxZ=0;
    
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
            aZ = aBarrelModule->getEdgeZSide(direction).first;
            if ((aZ>maxZ)||firstModule) {
                maxBarrelModule = aBarrelModule;
                maxZ=aZ;
            }
            if (firstModule) firstModule=false;
        }
    }
    
    // Measure the Delta of the farthest module
    double Delta;
    Delta = newMaxZ - maxZ;
    
    // TODO: if Deltam < 0 raise an error
    
    // delta is the ratio between Delta and the module's loewst Z w.r.t. newMinZ
    double delta=0;
    double myMinZ;
    XYZVector modShift;
    if (maxBarrelModule!=NULL) {
        delta=Delta/((maxBarrelModule->getEdgeZSide(direction*-1)).first-newMinZ);
        
        for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
            if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
                myMinZ = ((aBarrelModule->getEdgeZSide(direction*-1)).first-newMinZ);
                if (myMinZ>0) {
                    modShift=XYZVector(0, 0, delta*myMinZ);
                }
                aBarrelModule->translate(modShift);
            }
        }
    } else {
      addMessage("int BarrelLayer::compressExceeding couldn't find the maxZ module", ERROR);
    }
    
    
    // TODO: add exit code
}

double BarrelLayer::computeAverageRadius() {
    ModuleVector::iterator modIt;
    BarrelModule* aBarrelModule;
    double myRadius=0;
    int nModules=0;
    
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        if ( (aBarrelModule=dynamic_cast<BarrelModule*>(*modIt)) ) {
            myRadius+=aBarrelModule->getMeanPoint().Rho();
            nModules++;
        }
    }
    
    myRadius /= nModules;
    averageRadius_=myRadius;
    return myRadius;
}

/******************/
/*                */
/* Endcap layer   */
/*                */
/******************/

EndcapLayer::~EndcapLayer() {
    if (sampleModule_) delete sampleModule_;
}

EndcapLayer::EndcapLayer() {
    sampleModule_ = NULL;
}

EndcapLayer::EndcapLayer(EndcapLayer& inputLayer) {
    ModuleVector::iterator modIt;
    ModuleVector inputModuleV = inputLayer.moduleSet_;
    EndcapModule* aModule;
    EndcapModule* stdModule;
    
    setName(inputLayer.layerName_);
    averageZ_  = inputLayer.averageZ_;
    nOfRings_ = inputLayer.getRings();
    nModsOnRing_.clear();
    for (unsigned int i = 0; i < inputLayer.getModulesOnRing().size(); i++) {
        nModsOnRing_.push_back(inputLayer.getModulesOnRing().at(i));
    }
    
    if ( (sampleModule_=dynamic_cast<EndcapModule*>(inputLayer.getSampleModule())) ) {
        sampleModule_ = new EndcapModule(*(inputLayer.getSampleModule()));
    }
    
    for (modIt=inputModuleV.begin(); modIt!=inputModuleV.end(); modIt++) {
        if ( (stdModule=dynamic_cast<EndcapModule*>(*modIt)) ) {
            aModule = new EndcapModule(*stdModule);
            moduleSet_.push_back(aModule);
        } else {
	  addMessage("in EndcapLayer::EndcapLayer(EndcapLayer& inputLayer) I found a non-endcap module in the source", ERROR);
        }
    }
}


EndcapLayer::EndcapLayer(const EndcapModule& sample, double alpha, double d) {
    sampleModule_ = new EndcapModule(sample, alpha, d, -1);
}

EndcapLayer::EndcapLayer(double alpha, double d) {
    sampleModule_ = new EndcapModule(alpha, d, -1);
}

EndcapLayer::EndcapLayer(const EndcapModule& mySample) {
    sampleModule_ = new EndcapModule(mySample);
}

//   void rotateY(double alpha);
void EndcapLayer::translateZ(const double& zShift) {
    ModuleVector::iterator modIt;
    
    averageZ_+=zShift;
    
    XYZVector myShift(0, 0, zShift);
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); modIt++) {
        (*modIt)->translate(myShift);
    }
}


#define MAX_LOOPS 10000

// y=r^2/l^2
// R: circle radius
// L: line origin (negative)
// x=sin^2(alpha)
// alpha is the line angle wrt x axis
double EndcapLayer::solvex(double y) {
    return(1/8.*(3*y+2-sqrt(9*pow(y, 2)-4*y+4)));
}

double EndcapLayer::gamma1(double x, double y, double r) {
    double result;
    result = sqrt(1-x) - sqrt(y-x);
    result *= r / sqrt(y);
    return result;
}

double EndcapLayer::gamma2(double x, double y, double r) {
    double result;
    result = sqrt(1-x) + sqrt(y-x);
    result *= r / sqrt(y);
    return result;
}

double EndcapLayer::Area(double x, double y, double r) {
    double result;
    result = sqrt((1-(x/y))*x/y);
    result *= 4*pow(r, 2)*(1-x);
    return result;
}

double EndcapLayer::compute_l(double x, double y, double d) {
    double result = d;
    result /= (1-x)-sqrt((1-x)*(y-x));
    return result;
}

double EndcapLayer::compute_d(double x, double y, double l) {
    double result = l;
    result *= (1-x)-sqrt((1-x)*(y-x));
    return result;
}

void EndcapLayer::buildSingleDisk(double minRadius,
				  double maxRadius,
				  double smallDelta,
				  double bigDelta,
				  double diskZ,
				  double overlap,
				  double zError,
				  int base,
				  EndcapModule* sampleModule,
				  std::map<int, int> ringDirectives,
				  int diskParity, /*=-1*/
				  int sectioned /*=NoSection*/) {
    
    averageZ_=diskZ;
    
    EndcapModule* aRingModule;
    
    EndcapModule* trialModule[3];
    
    std::map<int, int>::iterator aDirective;
    
    
    int ringParity;
    int nearDirection = int(diskZ/fabs(diskZ))*-1;
    int nRing;
    int addModules;
    
    double lastRho = 0;
    double nextRho = minRadius;
    double destZ;
    XYZVector shiftThis;
    edge aSide;
    edge minSafetyEdge;
    std::ostringstream tag;
    
    for (nRing=1; lastRho<maxRadius; nRing++) {
        // Ring parity is 1, -1, 1, -1 for
        //      nRing=    1,  2, 3,  4
        // if diskParity==1 and the opposite, otherwise
        ringParity = diskParity*(((nRing%2)*2)-1);
        
        sampleModule->setRing(nRing);
        
        // Debug
	tempString.str(""); tempString  << "Looking for directives of ring " << nRing;
	addMessage(tempString, INFO);
        
        aDirective = ringDirectives.find(nRing);
        if (aDirective!=ringDirectives.end()) {
            addModules = ringDirectives[nRing];
            addMessage("Found a directive", INFO);
        } else {
            addModules = 0;
            addMessage("Found no directive", INFO);
        }
        
        // ringParity = 1 means the ring is nearer to the interaction point
        lastRho = buildRing(nextRho,
                smallDelta,
                ringParity*bigDelta,
                diskZ,
                overlap,
                base,
                nearDirection,
                sampleModule,
                maxRadius,
                addModules,
                sectioned);
        
	tempString.str(""); tempString << "maxRadius: " << maxRadius;
	addMessage(tempString, INFO);
	tempString.str(""); tempString << "lastrho: " << lastRho;
	addMessage(tempString, INFO);
        aRingModule = new EndcapModule(*sampleModule, 100/lastRho, lastRho, -1);
        aRingModule->setRing(nRing);
        
        
        // Place the mock-up module so as to simulate the worst position for this ring
        shiftThis = XYZVector(0,
                0,
                diskZ + -1*nearDirection*smallDelta + ringParity*nearDirection*bigDelta);
        
        // Compute the destination z
        destZ = diskZ + nearDirection*smallDelta + -1*ringParity*nearDirection*bigDelta;
        
        aRingModule->translate(shiftThis);
	tempString.str(""); tempString << "Pushing from z=" << diskZ + -1*nearDirection*smallDelta + ringParity*nearDirection*bigDelta
				       << " to z=" << destZ;
	addMessage(tempString, INFO);
        
        for (int i=0; i<3; i++) {
            trialModule[i] = new EndcapModule(*aRingModule);
        }
        delete aRingModule;
        
        // Three tests: with a border and projecting form +- zError
        
        int direct;
        
        for (int i=0; i<3; i++) {
            direct=i-1;
            aSide = trialModule[i]->getEdgeRhoSide(-1);
            trialModule[i]->projectSideZ(aSide.second, destZ, direct*zError);
	    tempString.str(""); tempString << "Fake module base: " << trialModule[i]->getEdgeRhoSide(-1).first;
	    addMessage(tempString, INFO);
            
            if (i==1) {
                // trialModule[1] is the one with no zError, but overlap
                // TODO: here I will use the right setedgeside and getedgeside of the
                // endcapmodule, after they're debugged
                trialModule[i]->translate(XYZVector(0., -1*overlap, 0.));
            }
            
        }
        
        minSafetyEdge = trialModule[0]->getEdgeRhoSide(-1);
        minSafetyEdge.second = 0;
        
        // Debug:
        for (int i=0; i<3; i++) {
	  tempString.str("");
	  tempString << "safetyEdge["
		     << i << "]="
		     << trialModule[i]->getEdgeRhoSide(-1).first;
	  addMessage(tempString, INFO);
        }
        for (int i=1; i<3; i++) {
            aSide = trialModule[i]->getEdgeRhoSide(-1);
            if (aSide.first<minSafetyEdge.first) {
                minSafetyEdge.first = aSide.first;
                minSafetyEdge.second = i;
            }
        }
        
        // Now we just keep the safest, while we delete the others
        for (int i=0; i<3; i++) {
            delete trialModule[i];
        }
        
        tempString.str(""); tempString << "Ring computation: using safety rule #" << minSafetyEdge.second;
        nextRho = minSafetyEdge.first;
        tempString << " next radius at rho: " << nextRho;
	addMessage(tempString, INFO);
    }
    nOfRings_ = nRing - 1;
    
}

// Returns a 'standard' module
double EndcapLayer::buildRing(double minRadius,
        double smallDelta,
        double bigDelta,
        double diskZ,
        double overlap,
        int base,
        int nearDirection,
        EndcapModule* sampleModule,
        double maxRadius /* = -1 */,
        int addModules, /* = 0 */
        int sectioned /* = NoSection */) {
    
    
    double alpha; // Module angular aperture
    bool wedges=false;
    
    // Used only for barrel modules
    double modMaxRadius=1;
    
    if ( sampleModule->getShape()==Module::Wedge ) {
        wedges=true;
        double r = sampleModule->getDiameter()/2.;
        
	tempString.str("");
        tempString << "Desired distance: " << minRadius << std::endl;
        tempString << "Desired overlap: " << overlap << std::endl;
        double l=(minRadius-r);
        tempString << "l: " << l << std::endl;
        tempString << "r: " << r << std::endl;
        double y=pow(r/l, 2);
        tempString << "y: " << y << std::endl;
        double x=solvex(y);
        tempString << "x: " << x << std::endl;
        tempString << "gamma1: " << gamma1(x, y, r) << std::endl;
        tempString << "gamma2: " << gamma2(x, y, r) << std::endl;
        
        // Max area
        tempString << "Area : " << Area(x, y, r) << std::endl;
        tempString << "Optimization: " << std::endl;
        
        
        double tempd;
        bool computeFinish=false;
        
        for (int i=0; i < MAX_LOOPS; i++) {
            l=compute_l(x, y, minRadius);
            y=pow(r/l, 2);
            x=solvex(y);
            
            tempd=compute_d(x, y, l);
            
            if (fabs(minRadius-tempd)<1e-15) {
                computeFinish=true;
                break;
            }
            
        }
        
        if (!computeFinish) {
	  tempString << "Computation not finished" << std::endl;
          tempString << "Delta_d is now: " << minRadius-tempd << std::endl;
        }
        tempString << "x: " << x << std::endl;
        tempString << "Area : " << Area(x, y, r) << std::endl;
        
        alpha=asin(sqrt(x))*2;
        
        tempString << "mod aperture: " << alpha*180./M_PI;
	addMessage(tempString, INFO);
    } else if ( sampleModule->getShape()==Module::Rectangular ) { // it's a barrel module
        modMaxRadius=minRadius + sampleModule->getHeight();
        // widthHi or widthLo would give the same number: it's a square module!
        alpha=2*asin(sampleModule->getWidthHi()/2. / modMaxRadius);
    } else { // It's not a barrel or an endcap module
      tempString.str("");
      tempString << "The sample EndcapModule was not Barrel nor Endcap... this should never happen!";
      addMessage(tempString, ERROR);
    }
    
    // The needed overlap becomes an angle delta by
    // checking the unsafest point (r=r_min)
    double delta, effectiveAlpha;
    
    if (wedges) {
        delta = overlap / minRadius;
    } else {
        delta = overlap / modMaxRadius;
    }
    effectiveAlpha = alpha - delta;
    
    
    tempString.str(""); tempString << "overlap delta: " << delta*180./M_PI;
    addMessage(tempString, INFO);
    tempString.str(""); tempString << "mod aperture (effective): " << effectiveAlpha*180./M_PI;
    addMessage(tempString, INFO);
    
    
    double n;
    n = 2*M_PI / effectiveAlpha;
    
    tempString.str(""); tempString << "Number of modules: " << n;
    addMessage(tempString, INFO);
   
    // Optimal number of modules 
    int nOpt;
    
    if (wedges) {
        // We round to the nearest multiple of base
        nOpt = int(floor((n/double(base))+.5)) * base;
	tempString.str(""); tempString << "I would use " << nOpt << " modules (optimized), ";
	nOpt += (addModules*base);
        tempString << "i will use " << nOpt << " modules (user request)";
    } else {
        // For square modules we can only increase the number of sensors to
        // cover the whole area
        nOpt = int(ceil(n/double(base))) * base;
        tempString.str(""); tempString << "I will use " << nOpt << " modules";
    }
    addMessage(tempString, INFO);
    
    double goodAlpha;
    goodAlpha = 2*M_PI/double(nOpt);
    goodAlpha += delta;
    
    //   XYZVector diskShift   = XYZVector(0, 0, diskZ);
    //   XYZVector petalShift  = XYZVector(0, 0, smallDelta);
    //   XYZVector ringShift   = XYZVector(0, 0, bigDelta);
    int ringParity;
    EndcapModule *myModule;
    //BarrelModule *myBarrelModule;
    //Module* myModule;
    
    //   std::cout << "diskShift:  " << diskZ << std::endl
    // 	    << "petalShift: " << smallDelta << std::endl
    // 	    << "ringShift:  " << bigDelta << std::endl;
    
    for (int i=0; i<nOpt; i++) {
        ringParity = ((i%2)*2)-1;
        if (wedges) {
            myModule = new EndcapModule(*sampleModule, goodAlpha, minRadius, maxRadius);
        } else {
            myModule = new EndcapModule(*sampleModule, minRadius);
        }
        myModule->rotatePhi(2.*M_PI*i/double(nOpt));
        XYZVector shift = XYZVector(0, 0, diskZ + nearDirection*ringParity*smallDelta + nearDirection*bigDelta);
        myModule->translate(shift);
        if (i==0) myModule->setSection(YZSection);
        if (sectioned == NoSection) {
            moduleSet_.push_back(myModule);
        } else {
            if (sectioned == YZSection) {
                if (i==0) {
                    moduleSet_.push_back(myModule);
                }
            }
        }
    }
    nModsOnRing_.push_back(nOpt);
    
    double lastRho;
    if (wedges) {
        EndcapModule* aRingModule = new EndcapModule(*sampleModule, goodAlpha, minRadius, maxRadius);
        lastRho = minRadius + aRingModule->getHeight();
	tempString.str(""); tempString << "Actual module area: " << aRingModule->getArea();
	addMessage(tempString, INFO);
        // Just to be sure...!
        if (aRingModule->wasCut()) {
	  tempString.str(""); tempString << "The ring Module was cut, losing: " << aRingModule->getLost();
	  addMessage(tempString, INFO);
	  lastRho=maxRadius;
        }
        delete aRingModule;
    } else {
        lastRho = modMaxRadius;
    }
    
    return lastRho;
}

int EndcapLayer::cutOverEta(double etaCut) {
    int nCut = 0;
    ModuleVector::iterator modIt;
    
    double theta;
    double eta;
    
    for (modIt=moduleSet_.begin(); modIt!=moduleSet_.end(); ) {
        theta=(*modIt)->getMeanTheta();
        eta = -1*log(tan(theta/2.));
        if (fabs(eta)>etaCut) {
            // Bookkeeping
            decreaseModCount((*modIt)->getRing() - 1);
            // Erase the useless module!
            delete (*modIt);
            moduleSet_.erase(modIt);
            nCut++;
        } else {
            modIt++;
        }
    }
    
    return nCut;
}

double EndcapLayer::getMaxModuleThickness() {
    ModuleVector::iterator iter, guard = moduleSet_.end();
    double res = 0.0;
    for (iter = moduleSet_.begin(); iter != guard; iter++) {
        if ((*iter)->getModuleThickness() > res) res = (*iter)->getModuleThickness();
    }
    return res;
}

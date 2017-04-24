// #include <time.h>
// #include <sensor.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>

// #include <TProfile.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
// #include <TCanvas.h>
// #include <TStyle.h>
// #include <TFile.h>
// #include <TRandom3.h>
// #include <stdexcept>
#include "module.hh"
#include "Math/Vector3D.h"
#include "TGeoManager.h"
#include "TFile.h"
#include <TRandom3.h>
#include "TStyle.h"
#include "TPolyLine3D.h"


#include "stylePlot.h"

#define DEBUG_DIAM 14.142135624 // Diameter needed to make the edge = 10mm (1 cm)
#define RANDOM_SEED 0xcaffe


// The complete geometry
TGeoManager  *tkGeom;

// Materials and media
TGeoMaterial *matVacuum;
TGeoMaterial *matSi;
TGeoMedium *Vacuum;
TGeoMedium *Si;

// The top volume
TGeoVolume *top;

// A module counter
int iModule;


void initializeGeometry() {
  // Define the geometry
  tkGeom = new TGeoManager("trackerGeom", "Tracker simple geometry");
  // Define materials and media
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMaterial *matSi = new TGeoMaterial("Si", 26.98,13,2.7);
  TGeoMedium *Vacuum = new TGeoMedium("Vacuum", 1, matVacuum);
  /* TGeoMedium *Si = */ new TGeoMedium("Si", 2, matSi);
  // Define and set the top volume
  top = tkGeom->MakeBox("TOP", Vacuum, 270., 270., 120.);
  tkGeom->SetTopVolume(top);
  // Reset the module counter
  iModule=0;
}

void closeSaveGeometry() {
  //  tkGeom->SetVisLevel(4); 
  tkGeom->CloseGeometry();
  tkGeom->Export("tkGeometry.root");
  delete tkGeom;
  delete matVacuum;
  delete matSi;
  delete Vacuum;
  delete Si;
}

void placeModule(Module* aModule) {
  std::string moduleName;
  char moduleNumber[100];
  sprintf(moduleNumber, "%d", iModule++);
  moduleName = "mod_" + std::string(moduleNumber);
  aModule->setId(moduleName);
  aModule->shapeVolume(top, Si, tkGeom);
}

void placeModuleLite(Module* aModule) {
  std::string moduleName;
  char moduleNumber[100];
  sprintf(moduleNumber, "%d", iModule++);
  moduleName = "mod_" + std::string(moduleNumber);
  aModule->setId(moduleName);
  TPolyLine3D* contour = aModule->getContour();
  contour->Draw();
}

void createGeom(std::vector<std::vector<Module*>* > myGeo) {
  std::vector<std::vector<Module*>* >::iterator layIt;
  std::vector<Module*>::iterator modIt;


  initializeGeometry();
  
  for (layIt=myGeo.begin(); layIt!=myGeo.end(); layIt++) {
    for (modIt=(*layIt)->begin(); modIt!=(*layIt)->end(); modIt++) {
      placeModule(*modIt);
    }
  }
  
  closeSaveGeometry();
}

void createGeomLite(std::vector<std::vector<Module*>* > myGeo) {
  std::vector<std::vector<Module*>* >::iterator layIt;
  std::vector<Module*>::iterator modIt;

  TCanvas *c1 = new TCanvas("mainCanvas");

  for (layIt=myGeo.begin(); layIt!=myGeo.end(); layIt++) {
    for (modIt=(*layIt)->begin(); modIt!=(*layIt)->end(); modIt++) {
      placeModuleLite(*modIt);
    }
  }
  
  c1->SaveAs("geomLite.root");
}




// TEST 1: basic movements
void test1(std::vector<Module*> &trialLayer)
{
  BarrelModule* myModule  = new BarrelModule(DEBUG_DIAM,1);
  BarrelModule* myModule2 = new BarrelModule(DEBUG_DIAM,1);
  trialLayer.push_back(myModule);
  trialLayer.push_back(myModule2);
  
  myModule2->rotatePhi(M_PI/2.*1.5);
  XYZVector shift(0, 6, 0);  
  myModule2->translate(shift);
  myModule2->setEdgeZ(5, +1);
}

// TEST 2: performance
void test2(std::vector<Module*> &trialLayer)
{
  BarrelModule* aModule;
  XYZVector shift(0, 0, 0);
  
  for (int i=0; i<10000; i++) {
    aModule = new BarrelModule(DEBUG_DIAM,1);
    aModule->translate(shift);
    shift += XYZVector(5, 5, 5);
    aModule->rotatePhi(M_PI/90.*i); // 2 degrees
    trialLayer.push_back(aModule);
  }
}

// TEST 3: edge movements
void test3(std::vector<Module*> &trialLayer)
{
  BarrelModule* myModule = new BarrelModule(DEBUG_DIAM, 1);
  trialLayer.push_back(myModule);

  edge myEdge = myModule->getEdgeZSide(+1, 2);
  
  /* XYZVector* pippo = */ myModule->marginBorderSide(10,myEdge.second);
  
  BarrelModule* anotherModule = new BarrelModule(*myModule);
  anotherModule->setEdgeZ(myEdge.first, +1);
  trialLayer.push_back(anotherModule);

  
//   XYZVector shhh(0,0,97.3462);
//   anotherModule.projectSideRho(myEdge.second, 1010, pippo, &shhh);
//   anotherModule.print();
}



// TEST 4: projection movement
void test4(std::vector<Module*> &trialLayer)
{
  BarrelModule* myModule = new BarrelModule(DEBUG_DIAM, 1);
  trialLayer.push_back(myModule);
  myModule = new BarrelModule(DEBUG_DIAM, 1);
  trialLayer.push_back(myModule);
  myModule->translate(XYZVector(0, 4, 0));
  edge myEdge = myModule->getEdgeZSide(+1, 0.);


  BarrelModule* anotherModule = new BarrelModule(*myModule);
  trialLayer.push_back(anotherModule);
  anotherModule->setEdgeZ(myEdge.first, +1);

  for (double destr=5; destr<10; destr+=1.) {
    myModule = new BarrelModule(*anotherModule);
    trialLayer.push_back(myModule);
    myEdge = myModule->getEdgeZSide(-1, 0.);
    myModule->projectSideRho(myEdge.second, destr);
  }


  for (double destr=10; destr<15; destr+=1.) {
    myModule = new BarrelModule(*anotherModule);
    trialLayer.push_back(myModule);
    myEdge = myModule->getEdgeZSide(-1, 0.);
    myModule->projectSideRho(myEdge.second, destr, +5.);
  }

  for (double destr=-10; destr>-15; destr-=1.) {
    myModule = new BarrelModule(*anotherModule);
    trialLayer.push_back(myModule);
    myEdge = myModule->getEdgeZSide(-1, 0.);
    myModule->projectSideRho(myEdge.second, destr, -5.);
  }
  
}


// Test 5: getEdgeZSide, setEdgeZSide
// with the margin option
void test5(std::vector<Module*> &trialLayer) {

  edge myEdge;
  BarrelModule* myModule;

  myModule= new BarrelModule(DEBUG_DIAM, 1);
  trialLayer.push_back(myModule);
  myModule->translate(XYZVector(0, 10, 0));

  // This will place a new module just touching the previous
  myEdge = myModule->getEdgeZSide(+1,0);
  myModule = new BarrelModule(*myModule);
  trialLayer.push_back(myModule);
  myModule->setEdgeZ(myEdge.first, +1);

  
  // This will place a new module penetrating the previous by 1
  myEdge = myModule->getEdgeZSide(+1,1);
  myModule = new BarrelModule(*myModule);
  trialLayer.push_back(myModule);
  myModule->setEdgeZ(myEdge.first, +1);


}

// Test 6: creation of an EndcapModule
void test6(std::vector<Module*> &trialLayer) {

  EndcapModule* myModule;
  int nModules = 20;

  for (int i=0; i<nModules; i++) {
    myModule = new EndcapModule(DEBUG_DIAM, 2.*M_PI/double(nModules), 10);
    myModule->rotatePhi(2.*M_PI*i/double(nModules));
    trialLayer.push_back(myModule);
  }

}


// Passed TESTs: 1, 2, 3, 4, 5, 6

// Verified functions:
// BarrelModule, EndcapModule
// rotatePhi
// translate
// getEdgeZSide (with margin option)
// setEdgeZ
// marginBorderSide
// projectSideRho

void tests(std::vector<std::vector<Module*>* >& myGeom) {

  std::vector<Module*>* trialLayer;

  trialLayer = new std::vector<Module*>;

  myGeom.push_back(trialLayer);

  // test1(*trialLayer); 
  // test2(*trialLayer);
  // test3(*trialLayer);
  // test4(*trialLayer);
  // test5(*trialLayer);
  test6(*trialLayer);

}

int main (int argc, char* argv[]) {
  std::vector<std::vector<Module*>* > myGeom;

  tests(myGeom);

  createGeom(myGeom);
  //createGeomLite(myGeom);
  
  return 0;

}

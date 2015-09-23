#ifndef _PIXELEXTRACTOR_H
#define	_PIXELEXTRACTOR_H

#include <tk2CMSSW_datatypes.h>
#include <tk2CMSSW_strings.h>
#include <set>
#include <cmath>
#include <sstream>
#include <Tracker.h>
#include <MaterialTable.h>
#include <MaterialBudget.h>
#include <MaterialObject.h>
#include "global_funcs.h"
#include "Property.h"
#include "Module.h"
#include "RodPair.h"
#include "Visitable.h"
#include<mainConfigHandler.h>
namespace insur {

  //typedef PtrVector<RodPair> Container;
    //Container& rods;
    typedef PtrVector<Ring> EndcapRings;
    typedef PtrVector<BarrelModule> BmoduleContainer;
    typedef PtrVector<EndcapModule> EmoduleContainer;
    typedef PtrVector<Sensor> Sensors;

  
  class PixelExtractor {
      
      enum volumeType{ hybrid=1,sensor,chip };
  
      class LayerAggregator : public GeometryVisitor { // CUIDADO quick'n'dirty visitor-based adaptor to interface with legacy spaghetti code
        std::vector<Layer*> barrelLayers_;
        std::vector<Disk*> endcapDisks_;
        std::vector<std::vector<ModuleCap> > layerCap_;
        public:
        void visit(Layer& l) { barrelLayers_.push_back(&l); }
        void visit(Disk& d) { endcapDisks_.push_back(&d); }
        std::vector<Layer*>* getBarrelLayers() { return &barrelLayers_; }
        std::vector<Disk*>* getEndcapDisks() { return &endcapDisks_; }
        //Extract barrel modules layer by layer
        std::vector<std::vector<ModuleCap> >& getBarrelCap() {
          if (layerCap_.size()) return layerCap_;
        //std:cout << "Computing my layerCap_ for the first time" << std::endl;
          for (auto it : barrelLayers_ ) {
            class ModuleVisitor : public GeometryVisitor {
            public:
               std::vector<ModuleCap> *myLayerModuleCaps_;
               ModuleVisitor() { myLayerModuleCaps_ = new std::vector<ModuleCap>(); }
               void visit(BarrelModule& bm) { myLayerModuleCaps_->push_back(*bm.getModuleCap()) ; }
            };
            ModuleVisitor mv;
            it->accept(mv);
            layerCap_.push_back(*mv.myLayerModuleCaps_);
            auto lastLayerCap = layerCap_.end();
            lastLayerCap--;
            }
            return layerCap_;
          }
        //Extract endcap modules
        std::vector<std::vector<ModuleCap> >& getEndcapCap() {
          if (layerCap_.size()) return layerCap_;
        //std:cout << "Computing my layerCap_ for the first time" << std::endl;
          for (auto it : endcapDisks_ ) {
            class ModuleVisitor : public GeometryVisitor {
            public:
               std::vector<ModuleCap> *myLayerModuleCaps_;
               ModuleVisitor() { myLayerModuleCaps_ = new std::vector<ModuleCap>(); }
               void visit(EndcapModule& em) { myLayerModuleCaps_->push_back(*em.getModuleCap()) ; }
            };
            ModuleVisitor mv;
            it->accept(mv);
            layerCap_.push_back(*mv.myLayerModuleCaps_);
            auto lastLayerCap = layerCap_.end();
            lastLayerCap--;
            }
            return layerCap_;
          }
      };
      
     public:
       void analyse( MaterialTable& mt,MaterialBudget& mb );

       //Functions Related to Barrel
       void createPixelBarrelActive(Tracker& pix);
       void getPixelBarrelModuleInfo( std::vector<BarrelModule> moduleVec ,int LayerNum,std::string rodName, double rod_Rmin, 
                                      double rod_dR, std::string barrelRmatpathCommon );
       void createPixelBarrelServices(InactiveSurfaces& inatcvser);
       //Functions Related to Endcap
       void createPixelEndcapActive(Tracker& pix);
       void getPixelEndcapModuleInfo( EndcapModule emodule ,int ringNo, int discno,
                                      std::string ringName, double ring_thickness, double rod_dR,
                                      std::string ecapRmatpathCommon );
       void createPixelEndcapServices(InactiveSurfaces& inatcvser);
      
       int getAtomicNumber(double x0, double A);
       void analyseElements(MaterialTable& mattab);
       void analyseCompositeElements( std::string name, double density,MaterialProperties& mp, bool nosensors);
       void addMaterialToVolume( std::string name, Module& module,int targetVoltag,double volume);

       double compositeDensity(ModuleCap& mc, bool nosensors);
       double compositeDensity(InactiveElement& ie);

       void printXml(mainConfigHandler& mainConfiguration, std::string outsubdir);
   private:
       std::string whichShape( ShapeType  s );
       unsigned int numBarrelLayers;
       std::vector<std::string> barrelRmatpath;
       std::vector<std::pair<unsigned int,unsigned int>> discRingpair;
       std::vector<std::string> ecapRmatpath;
       CMSSWBundle cmsswXmlInfo;
  };
}

#endif

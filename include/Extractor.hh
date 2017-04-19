// 
// File:   Extractor.h
// Author: ndemaio
//
// Created on January 31, 2010, 3:33 AM
//

/**
 * @file Extractor.h
 * @brief This is the header file for the class that analyses a full tracker and material budget for export
 */

#ifndef _EXTRACTOR_H
#define	_EXTRACTOR_H

#include <tk2CMSSW_datatypes.hh>
#include <tk2CMSSW_strings.hh>
#include <set>
#include <cmath>
#include <limits.h>
#include <sstream>
#include <Tracker.hh>
#include <MaterialTable.hh>
#include <MaterialBudget.hh>
#include <MaterialObject.hh>

namespace insur {
  /**
   * @class Extractor
   * @brief This class bundles the analysis functions that prepare an existing material budget and table for output to CMSSW XML.
   *
   * The only public function of the class receives the material budget and table that make up the input. The output goes into an instance
   * of a struct - provided as a reference - that bundles a series of vectors of internal datatypes. ons of internal data types. The analysis
   * results will be stored in those, ready to be formatted and written to file.
   */
  class Extractor {
    class LayerAggregator : public GeometryVisitor { // CUIDADO quick'n'dirty visitor-based adaptor to interface with legacy spaghetti code
      std::vector<Layer*> barrelLayers_;
      std::vector<Disk*> endcapLayers_;
      std::vector<std::vector<ModuleCap> > layerCap_;
    public:
      void visit(Layer& l) { barrelLayers_.push_back(&l); }
      void visit(Disk& d) { endcapLayers_.push_back(&d); }
      void postVisit() {
        //std::cout << "I am post-visiting" << std::endl;
        std::sort (barrelLayers_.begin(), barrelLayers_.end(), [ ]( Layer* lhs, Layer* rhs )
                   {
                      return lhs->minR() < rhs->minR();
                   });
        //for (auto it : barrelLayers_ ) {
        //  std::cout << "HERE layer " << it->minR() << std::endl;
        //}
      }
      std::vector<Layer*>* getBarrelLayers() { return &barrelLayers_; }
      std::vector<Disk*>* getEndcapLayers() { return &endcapLayers_; }
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
          //std::cout << "Added a layer with " << mv.myLayerModuleCaps_->size() << " modules" << std::endl;
          layerCap_.push_back(*mv.myLayerModuleCaps_);
          //std::cout << "layerCap_.size()=" << layerCap_.size() << std::endl;
          auto lastLayerCap = layerCap_.end();
          lastLayerCap--;
          //std::cout << "last layer cap has " << lastLayerCap->size() << " modules" << std::endl;
        }
        return layerCap_;
      }
    };
  public:
    void analyse(MaterialTable& mt, MaterialBudget& mb, XmlTags& trackerXmlTags, CMSSWBundle& d, bool wt = false);
  protected:
    void analyseElements(MaterialTable&mattab, std::vector<Element>& elems);
    void analyseBarrelContainer(Tracker& t, XmlTags& trackerXmlTags, std::vector<std::pair<double, double> >& up,
                                std::vector<std::pair<double, double> >& down);
    void analyseEndcapContainer(std::vector<std::vector<ModuleCap> >& ec, Tracker& t, XmlTags& trackerXmlTags, std::vector<std::pair<double, double> >& up,
                                std::vector<std::pair<double, double> >& down);
    void analyseLayers(MaterialTable& mt, Tracker& tr, XmlTags& trackerXmlTags, std::vector<Composite>& c,
                       std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<ShapeOperationInfo>& so, std::vector<PosInfo>& p, 
		       std::vector<AlgoInfo>& a, std::map<std::string,Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri, 
		       bool wt = false);
    void analyseDiscs(MaterialTable& mt, std::vector<std::vector<ModuleCap> >& ec, Tracker& tr, XmlTags& trackerXmlTags, std::vector<Composite>& c,
                      std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s, std::vector<ShapeOperationInfo>& so, std::vector<PosInfo>& p,
		      std::vector<AlgoInfo>& a, std::map<std::string,Rotation>& r, std::vector<SpecParInfo>& t, std::vector<RILengthInfo>& ri,
		      bool wt = false);
    void analyseServices(InactiveSurfaces& is, bool& isPixelTracker, XmlTags& trackerXmlTags,
			 std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
			 std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt = false);
    void analyseSupports(InactiveSurfaces& is, bool& isPixelTracker, XmlTags& trackerXmlTags,
			 std::vector<Composite>& c, std::vector<LogicalInfo>& l, std::vector<ShapeInfo>& s,
                         std::vector<PosInfo>& p, std::vector<SpecParInfo>& t, bool wt = false);
  private:
    Composite createComposite(std::string name, double density, MaterialProperties& mp, bool nosensors = false);
    std::vector<ModuleCap>::iterator findPartnerModule(std::vector<ModuleCap>::iterator i,
                                                       std::vector<ModuleCap>::iterator g, int ponrod, bool find_first = false);
    double findDeltaR(std::vector<Module*>::iterator start, std::vector<Module*>::iterator stop, double middle);
    double findDeltaZ(std::vector<Module*>::iterator start, std::vector<Module*>::iterator stop, double middle);
    int findSpecParIndex(std::vector<SpecParInfo>& specs, std::string label);
    double calculateSensorThickness(ModuleCap& mc, MaterialTable& mt);
    std::string stringParam(std::string name, std::string value);
    std::string numericParam(std::string name, std::string value);
    std::string vectorParam(double x, double y, double z);
    double compositeDensity(ModuleCap& mc, bool nosensors = false);
    double compositeDensity(InactiveElement& ie);
    double fromRim(double r, double w);
    std::pair<double, int> getAZ(double radiationLength, double interactionLength);
    int findClosest_Z(double pA, double radiationLength);
    double compute_X0(int pZ, double pA); 
  };

#if defined(__FLIPSENSORS_IN__) || defined(__FLIPSENSORS_OUT__)
  // For flipping one of sensors
  // This should be moved to tk2CMSSW_strings.h at some point.
  static const std::string rot_sensor_tag = "SensorFlip";
#endif
  class ModuleComplex {
    public :
     ModuleComplex(std::string moduleName, std::string parentName, ModuleCap& modcap);
     ~ModuleComplex();
     void buildSubVolumes();
     void addShapeInfo   (std::vector<ShapeInfo>&   vec);
     void addLogicInfo   (std::vector<LogicalInfo>& vec);
     void addPositionInfo(std::vector<PosInfo>&     vec);
     void addMaterialInfo(std::vector<Composite>&   vec);
     void print() const;

     const double getServiceHybridWidth() const { return serviceHybridWidth; }
     const double getFrontEndHybridWidth() const { return frontEndHybridWidth; }
     const double getExpandedModuleWidth() const { return expandedModWidth; }
     const double getExpandedModuleLength() const { return expandedModLength; }
     const double getExpandedModuleThickness() const { return expandedModThickness; }
     double getRmin() const { return rmin; }
     double getRmax() const { return rmax; }
     double getXmin() const { return xmin; }
     double getXmax() const { return xmax; }
     double getYmin() const { return ymin; }
     double getYmax() const { return ymax; }
     double getZmin() const { return zmin; }
     double getZmax() const { return zmax; }
     double getRminatZmin() const { return rminatzmin; }
     double getRmaxatZmax() const { return rmaxatzmax; }
     double getHybridTotalVolume_mm3() const { return hybridTotalVolume_mm3; }
     void   setHybridTotalVolume_mm3( double v ) { hybridTotalVolume_mm3 = v; }

    private :
      static const double kmm3Tocm3;
      static const int HybridFBLR_0; // Front Back Left Right
      static const int InnerSensor;
      static const int OuterSensor;
      static const int HybridFront;
      static const int HybridBack;
      static const int HybridLeft;
      static const int HybridRight;
      static const int HybridBetween;
      static const int SupportPlate;
      static const int HybridFB;
      static const int HybridLR;
      static const int HybridFBLR_3456; // Front Back Left Right (ailias of HybridFBLR_0)
      static const int PixelModuleNull;
      static const int PixelModuleHybrid;
      static const int PixelModuleSensor;
      static const int PixelModuleChip;    
      
      class Volume {
        public :
          Volume(std::string name, const int type, std::string pname, 
                 double dx,   double dy,   double dz,
                 double posx, double posy, double posz) : fname(name),
                                                          ftype(type),
                                                          fparentname(pname),
                                                          fdx(dx),
                                                          fdy(dy),
                                                          fdz(dz),
                                                          fx(posx),
                                                          fy(posy),
                                                          fz(posz),
                                                          fdxyz(dx*dy*dz*kmm3Tocm3), // mm3 to cm3
                                                          fdensity(-1.),
                                                          fmass(0.){}

          const double      getVolume()    const { return fdxyz; }
          const double      getDx()        const { return fdx;   }
          const double      getDy()        const { return fdy;   }
          const double      getDz()        const { return fdz;   }
          const double      getX()         const { return fx;    }
          const double      getY()         const { return fy;    }
          const double      getZ()         const { return fz;    }
          const std::string getName()      const { return fname; }
          const int         getType()      const { return ftype; }
          const std::string getParentName()const { return fparentname; }

          double getDensity() {
            if ( fdensity < 0. ) fdensity = fmass/fdxyz;
            return fdensity;
          }
          const std::map<std::string, double>& getMaterialList() const { return fmatlist; }
          void   addMaterial( std::string tag, double val ) { fmatlist[tag] += val; }
          void   addMass( double dm ) { fmass += dm; }
          double getMass() const { return fmass; }
          void print() const { 
            if (fdensity<0.) std::cout << "!!! Error : please call setDensity() before print()." << std:: endl; 
            else             std::cout << "   " << fname << " mass =" << fmass 
                                       <<" volume=" << fdxyz << " cm3, density=" << fdensity << " g/cm3" <<std::endl; 
          }
        private :
          std::string  fname;
          const int    ftype;
          std::string  fparentname;
          const double fdx,fdy,fdz;
          const double fx,fy,fz;
          const double fdxyz;
          double       fdensity;
          double       fmass;
          std::map<std::string, double> fmatlist;
      };

      ModuleCap&           modulecap;
      Module&              module;
      std::vector<Volume*> volumes;
      std::string          moduleId;
      std::string          parentId;
      const double         modThickness;
      const double         sensorThickness;
      const double         sensorDistance;
      const double         modWidth;  // Sensor width
      const double         modLength; // Sensor length
      const double         frontEndHybridWidth;
      const double         serviceHybridWidth;
      const double         hybridThickness;
      const double         supportPlateThickness;
      const double         chipThickness;
            double         hybridTotalMass;
            double         hybridTotalVolume_mm3;
            double         hybridFrontAndBackVolume_mm3;
            double         hybridLeftAndRightVolume_mm3;
            double         moduleMassWithoutSensors_expected;
            double         rmin;
            double         rmax;
	    double         xmin;
            double         xmax;
	    double         ymin;
            double         ymax;
            double         zmin;
            double         zmax;
	    double         rminatzmin;
	    double         rmaxatzmax;
	    const double         expandedModWidth;
	    const double         expandedModLength;
	    double               expandedModThickness; 
	    const XYZVector      center; 
	    const XYZVector      normal; 
	    int                  nTypes;
      std::vector<XYZVector> vertex; 
      std::string    prefix_xmlfile;
      const std::string    prefix_material;
  };
}
#endif	/* _EXTRACTOR_H */


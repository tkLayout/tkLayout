#ifndef PTERROR_H
#define PTERROR_H

#include "Module.hh"

#include <MaterialProperties.hh>


class ptError {
  public:
  private:
   
   // Global parameters
   static double IP_length;  // in mm // TODO: setter and getter for this parameter
   static double B;          // in T

   static double minimumPt;
   static double maximumPt;

   // Module parameters
   // (These parameters need get and set) // TODO
   double Module_pitch;    // in mm
   double Module_strip_l;  // in mm
   double Module_z;        // in mm
   double Module_r;        // Module's distance form z axis [mm] 
   double Module_d;        // Module thickness [mm]
   double Module_h;        // Module height [mm]
   double Module_ed;       // Module effective ds distance [mm]
   double Module_tilt;     // Module tilt angle [radians]
   ZCorrelation zCorrelation;
 
   RILength material;      // Material amount prior to this module
   int moduleType;
   int endcapType;

   static constexpr double defaultModulePitch = 90/1000.; // 90 um
   static constexpr double defaultStripLength = 23.2; // mm
   static constexpr double defaultModuleZ = 1346; // mm // Disk 1
   static constexpr double defaultModuleR = 1080; // last layer 
   static constexpr double defaultModuleD = 2; // mm
   static constexpr double defaultModuleHeight = 100;
   static const ZCorrelation defaultZCorrelation = SAMESEGMENT;

   // Internal computers
   void defaultParameters();
   double oneSidedGaussInt(double x);
   double probabilityInside_norm(double cut, double value);

  public:
   // Constructor and destructor
   ptError() {}
   ~ptError() {}

   // Parameter setters
   void setPitch(double newPitch) { Module_pitch = newPitch ; }
   void setStripLength(double newLength) { Module_strip_l = newLength ; }
   void setZ(double newz) { Module_z = newz ; }
   void setR(double newr) { Module_r = newr ; }
   void setDistance(double newDistance) { Module_d = newDistance ; }
   void setHeight(double newHeight) { Module_h = newHeight; }
   void setZCorrelation(ZCorrelation newZCorrelation) { zCorrelation = newZCorrelation; }
   void setModuleType(int newType) { moduleType = newType ; }
   void setEndcapType(int newType) { endcapType = newType ; }
   void setMaterial(RILength newMaterial) { material = newMaterial; }
   void setEffectiveDistance(double newDistance) { Module_ed = newDistance; }
   void setTilt(double newTilt) { Module_tilt = newTilt; }

   // Parameter getters
   double getPitch() const { return Module_pitch ; }
   double getStripLength() const { return Module_strip_l ; }
   double getZ() const { return Module_z ; }
   double getR() const { return Module_r ; }
   double getDistance() const { return Module_d ; }
   double getHeight() const { return Module_h ; }
   ZCorrelation getZCorrelation() const { return zCorrelation; }
   int getModuleType() const { return moduleType ; }
   int getEndcapType() { return endcapType ; }
   RILength getMaterial() { return material ; }

   // Conversion between strips and p
//   double stripsToP(double strips);
//   double pToStrips(double p);
   // Error computation
   double computeErrorBE(double p); // CUIDADO only for testing, called by computeError after previously setting the correct parameters
   double computeError(double p);
   // Efficiency computation
   double probabilityInside(double cut, double value, double value_err);
};

#endif

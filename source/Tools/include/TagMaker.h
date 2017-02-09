#ifndef TAGMAKER_H
#define TAGMAKER_H

#include <string>
#include <sstream>

#include "Visitor.h"
#include "DetectorModule.h"

using std::string;
using std::stringstream;

class TagMaker : public ConstGeometryVisitor {
public:
  string sensorTag, sensorGeoTag, sensorResTag, posTag;

  TagMaker(const DetectorModule& m) { m.accept(*this); }

  void visit(const BarrelModule& m) {
    sensorTag    = makeSensorTag(m);
    sensorGeoTag = makeSensorGeoTag(m);
    sensorResTag = makeSensorResolTag(m);
    posTag       = makePosTag(m);
  }

  void visit(const EndcapModule& m) {
    sensorTag    = makeSensorTag(m);
    sensorGeoTag = makeSensorGeoTag(m);
    sensorResTag = makeSensorResolTag(m);
    posTag       = makePosTag(m);
  }


  static string makeSensorTag(const BarrelModule& m) {
    stringstream ss;
    ss << "Barrel" 
       << "/Tag=" << m.cntName() << "L" << setfill('0') << setw(2) << m.layer() // << "R" << m.ring();
       << "/Width=" << m.maxWidth() 
       << "/Height=" << m.length()
       << "/Thickness=" << m.thickness()
       << "/Type=" << m.moduleType();
    return ss.str();
  }

  static string makeSensorTag(const EndcapModule& m) {
    stringstream ss;
    ss << "Endcap"
       << "/Tag=" << m.cntName() << "R" << setfill('0') << setw(2) << m.ring() << "D" << m.disk() 
       << "/WidthLo=" << m.minWidth()
       << "/WidthHi=" << m.maxWidth()
       << "/Height="  << m.length()
       << "/Thickness=" << m.thickness()
       << "/Type=" << m.moduleType();
    return ss.str();
  } 

  static string makeSensorGeoTag(const BarrelModule& m) {
    stringstream ss;
    ss << "Barrel"
       << "/Width=" << int(m.maxWidth()*1000)/1000.
       << "/Height=" << int(m.length()*1000)/1000.
       << "/Thickness=" << int(m.thickness()*1000)/1000.
       << "/Type=" << m.moduleType()
       << "/StripsAcross=" << m.outerSensor().numStripsAcross()
       << "/Faces=" << m.numSensors();
    for (const auto& s : m.sensors()) ss << "/Segments=" << s.numSegments();

    return ss.str();
  }

  static string makeSensorGeoTag(const EndcapModule& m) {
    stringstream ss;
    ss << "Endcap"
       << "/WidthLo=" << int(m.minWidth()*1000)/1000.
       << "/WidthHi=" << int(m.maxWidth()*1000)/1000.
       << "/Height=" << int(m.length()*1000)/1000.
       << "/Thickness=" << int(m.thickness()*1000)/1000.
       << "/Type=" << m.moduleType()
       << "/StripsAcross=" << m.outerSensor().numStripsAcross()
       << "/Faces=" << m.numSensors();
    for (const auto& s : m.sensors()) ss << "/Segments=" << s.numSegments();

    return ss.str();
  }

  static string makeSensorResolTag(const BarrelModule& m) {
    stringstream ss;
    ss << "Barrel"
       << "/Width=" << int(m.maxWidth()*1000)/1000.
       << "/Height=" << int(m.length()*1000)/1000.
       << "/Thickness=" << int(m.thickness()*1000)/1000.
       << "/Type=" << m.moduleType()
       << "/ResolutionX=" << m.resolutionLocalX()
       << "/ResolutionY=" << m.resolutionLocalY()
       << "/Faces=" << m.numSensors();
    for (const auto& s : m.sensors()) ss << "/Segments=" << s.numSegments();

    return ss.str();
  }

  static string makeSensorResolTag(const EndcapModule& m) {
    stringstream ss;
    ss << "Barrel"
       << "/Width=" << int(m.maxWidth()*1000)/1000.
       << "/Height=" << int(m.length()*1000)/1000.
       << "/Thickness=" << int(m.thickness()*1000)/1000.
       << "/Type=" << m.moduleType()
       << "/ResolutionX=" << m.resolutionLocalX()
       << "/ResolutionY=" << m.resolutionLocalY()
       << "/Faces=" << m.numSensors();
    for (const auto& s : m.sensors()) ss << "/Segments=" << s.numSegments();

    return ss.str();
  }

  static string makePosTag(const BarrelModule& m) {
    stringstream ss;
    ss << m.cntName() << "_L" << setfill('0') << setw(2) << m.layer() /*<< "R" << m.ring() */;
    return ss.str();
  }

  static string makePosTag(const EndcapModule& m) {
    stringstream ss;
    ss << m.cntName() << "_R" << setfill('0') << setw(2) << m.ring() << "D" << m.disk();
    return ss.str();
  }

};


#endif

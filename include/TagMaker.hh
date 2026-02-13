#ifndef TAGMAKER_H
#define TAGMAKER_H

#include <string>
#include <sstream>

#include "Visitor.hh"
#include "Module.hh"

class TagMaker : public ConstGeometryVisitor {
public:
  std::string sensorTag, sensorGeoTag, posTag;

  TagMaker(const Module& m) { m.accept(*this); }

  void visit(const BarrelModule& m) {
    sensorTag = makeSensorTag(m);
    sensorGeoTag = makeSensorGeoTag(m);
    posTag = makePosTag(m);
  }

  void visit(const EndcapModule& m) {
    sensorTag = makeSensorTag(m);
    sensorGeoTag = makeSensorGeoTag(m);
    posTag = makePosTag(m);
  }


  static std::string makeSensorTag(const BarrelModule& m) {
    std::stringstream ss;
    ss << "Barrel" 
       << "/Tag=" << m.subdetectorName() << "L" << std::setfill('0') << std::setw(2) << m.layer() // << "R" << m.ring();
       << "/Width=" << m.maxWidth() 
       << "/Height=" << m.length()
       << "/Thickness=" << m.thickness()
       << "/Type=" << m.moduleType();
    return ss.str();
  }

  static std::string makeSensorTag(const EndcapModule& m) {
    std::stringstream ss;
    ss << "Endcap"
       << "/Tag=" << m.subdetectorName() << "R" << std::setfill('0') << std::setw(2) << m.ring() << "D" << m.disk() 
       << "/WidthLo=" << m.minWidth()
       << "/WidthHi=" << m.maxWidth()
       << "/Height="  << m.length()
       << "/Thickness=" << m.thickness()
       << "/Type=" << m.moduleType();
    return ss.str();
  } 

  static std::string makeSensorGeoTag(const BarrelModule& m) {
    std::stringstream ss;
    ss << "Barrel"
       << "/Width=" << int(m.maxWidth()*1000)/1000.
       << "/Height=" << int(m.length()*1000)/1000.
       << "/Thickness=" << int(m.thickness()*1000)/1000.
       << "/Type=" << m.moduleType()
       << "/StripsAcross=" << m.outerSensor().numStripsAcrossEstimate()
       << "/Faces=" << m.numSensors();
    for (const auto& s : m.sensors()) ss << "/Segments=" << s.numSegmentsEstimate();

    return ss.str();
  }
                          
  static std::string makeSensorGeoTag(const EndcapModule& m) {
    std::stringstream ss;
    ss << "Endcap"
       << "/WidthLo=" << int(m.minWidth()*1000)/1000.
       << "/WidthHi=" << int(m.maxWidth()*1000)/1000.
       << "/Height=" << int(m.length()*1000)/1000.
       << "/Thickness=" << int(m.thickness()*1000)/1000.
       << "/Type=" << m.moduleType()
       << "/StripsAcross=" << m.outerSensor().numStripsAcrossEstimate()
       << "/Faces=" << m.numSensors();
    for (const auto& s : m.sensors()) ss << "/Segments=" << s.numSegmentsEstimate();

    return ss.str();
  }    

  static std::string makePosTag(const BarrelModule& m) {
    std::stringstream ss;
    ss << m.subdetectorName() << "L" << std::setfill('0') << std::setw(2) << m.layer() /*<< "R" << m.ring() */;
    return ss.str();
  }

  static std::string makePosTag(const EndcapModule& m) {
    std::stringstream ss;
    ss << m.subdetectorName() << "R" << std::setfill('0') << std::setw(2) << m.ring();
    //<< "D" << m.disk();
    return ss.str();
  }
};

#endif

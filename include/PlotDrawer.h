#ifndef PLOT_DRAWER_H
#define PLOT_DRAWER_H

#include <utility>
#include <map>
#include <set>
#include <iostream>

#include <TPolyLine.h>
#include <TPaletteAxis.h>
#include <TList.h>
#include <TH2D.h>
#include <TH2C.h>
#include <TText.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

#include <Palette.h>
#include <Module.h>





// ========================================================================================
// Here be STATISTICS
// If possible, additional statistics should be local to the plot drawer instantiation
// ========================================================================================

class NoStat {
  double value_;
public:
  NoStat() : value_(0) {}
  void fill(double value) { value_ = value; }
  double get() const { return value_; }
};

class Average {
  int counts_;
  double total_;
public:
  Average() : counts_(0), total_(0) {}
  void fill(double value) {
    total_ += value;
    counts_+= 1;
  }
  double get() const { return total_/counts_; }
};

class Max {
  double curr_;
public:
  Max() : curr_(0) {}
  void fill(double value) { curr_ = (value > curr_) ? value : curr_; }
  double get() const { return curr_; }
};

class Min {
  double curr_;
public:
  Min() : curr_(0) {}
  void fill(double value) { curr_ = (value < curr_) ? value : curr_; }
  double get() const { return curr_; }
};

class Sum {
  double value_;
public:
  Sum() : value_(0) {}
  void fill(double value) { value_ += value; }
  double get() const { return value_; }
};



// ==============================================================================================
// Here be VALUEGETTERS
// If possible, additional value getters should be local to the plot drawer instantiation
// ==============================================================================================


template<class RetType, RetType (Module::*ModuleMethod)() const>
struct Method {
  double operator()(const Module& m) const { return (double)(m.*ModuleMethod)(); }
};



struct TypeAutoColor { // Auto-assign colors
  std::set<std::string> colorSet_;
  double operator()(const Module& m) { 
    std::pair<std::set<std::string>::iterator, bool> it = colorSet_.insert(m.moduleType());
    return Palette::color(std::distance(colorSet_.begin(), it.first)+1);
  }
};


struct Type { // Module-maintained color
  double operator()(const Module& m) { 
    return Palette::color(m.plotColor());
  }
};


struct CoordZ {
  double operator()(const Module& m) { return m.center().Z(); }
};




// =============================================================================================
// Here be DRAWSTYLES 
// If possible, additional draw styles should be local to the plot drawer instantiation
// =============================================================================================

class DrawerPalette {
  double minValue_, maxValue_;
  TPaletteAxis* framePalette_;
public:
  DrawerPalette() : framePalette_(NULL) {}
  void setMinMaxValues(double minValue, double maxValue) {
    minValue_ = minValue * 1.001 - maxValue * 0.001;
    maxValue_ = maxValue * 1.001 - minValue * 0.001;
  }
  double getMinValue() const { return minValue_; }
  double getMaxValue() const { return maxValue_; }

  void setFramePalette(TPaletteAxis* framePalette) { 
    framePalette_ = framePalette; 
  }

  int getColor(double value) const {
    return framePalette_ ? framePalette_->GetValueColor(value) : value;
  }
};




struct FillStyle {
  template <class StatType> void operator()(TPolyLine& line, StatType& bin, const DrawerPalette& palette) const {
    line.SetFillColor(palette.getColor(bin.get()));
    line.DrawPolyLine(line.GetN(),line.GetX(),line.GetY(),"f");
  }
};

extern int
g;

class ContourStyle {
  const int lineWidth_;
public:
  ContourStyle(int lineWidth = 2)  : lineWidth_(lineWidth) {}
  template <class StatType> void operator()(TPolyLine& line, StatType& bin, const DrawerPalette& palette) const {
    line.SetLineColor(palette.getColor(bin.get()));
    line.SetLineWidth(lineWidth_);
    line.DrawPolyLine(line.GetN(),line.GetX(),line.GetY());
  }
};


// ===============================================================================================
// Here be VALIDATORS 
// If possible, additional module validators should be local to the plot drawer instantiation
// ===============================================================================================



template<const int SubdetType>
struct CheckType {
  bool operator()(const Module& m) const { return m.subdet() & SubdetType; }
};

template<const int PhiIndex>
struct CheckPhiIndex {
  bool operator()(const Module& m) const { return m.posRef().phi == PhiIndex; }
};

// ===============================================================================================
// Here be COORDINATE CLASSES and LINEGETTER
// Careful what you do
// ===============================================================================================



struct Rounder {
  static const int mmFraction = 1000;  // Drawing with 1 micrometer precision
  int round(double x) { return floor(x*mmFraction+0.5); }
};

struct XY : public std::pair<int, int>, private Rounder {
  const bool valid;
 XY(const Module& m) : std::pair<int, int>(round(m.center().X()), round(m.center().Y())), valid(m.center().Z() >= 0) {}
 XY(const XYZVector& v) : std::pair<int, int>(round(v.X()), round(v.Y())), valid(v.Z() >= 0) {}
 XY(const XYZVector& v, const Module& m) : XY(v) {}
  // bool operator<(const XY& other) const { return (x() < other.x()) || (x() == other.x() && y() < other.y()); }
  int x() const { return this->first; }
  int y() const { return this->second; }
};


struct YZ : public std::pair<int, int>, private Rounder {
  const bool valid;
 YZ(const Module& m) : std::pair<int,int>(round(m.center().Z()), round(m.center().Rho())), valid(m.center().Z() >= 0) {}
 YZ(const XYZVector& v) : std::pair<int, int>(round(v.Z()), round(v.Rho())), valid(v.Z() >= 0) {}
 YZ(const XYZVector& v, const Module& m) : valid(v.Z() >= 0) {
    this->first = round(v.Z());

    XYZVector vProjected;
    XYZVector z(0., 0., 1.);

    XYZVector basePolyCenter = m.basePoly().getCenter();
    XYZVector normal = crossProduct(z, basePolyCenter);
    normal = normal.Unit();

    vProjected = v - v.Dot(normal) * normal;
    
    this->second = round(vProjected.Rho());
    /*if (vProjected.Rho() > 1000. && v.Z() > 500. && v.Z() < 650.) {
      //std::cout << "v.Rho() = " << v.Rho() << std::endl;
      std::cout << " basePolyCenter.X() = " << basePolyCenter.X() <<  "basePolyCenter.Y() = " << basePolyCenter.Y() <<  "basePolyCenter.Z() = " << basePolyCenter.Z() << std::endl;
      //std::cout << "v.Dot(normal) = " << v.Dot(normal) << std::endl;
      //std::cout << "vProjected.Rho() = " << vProjected.Rho() << "v.Z() = " << v.Z() << "vProjected.Z() = " << vProjected.Z() << std::endl;
      std::cout << "(m.basePolyCenter().Rho() - m.length() / 2.) = " << (m.basePolyCenter().Rho() - m.length() / 2.) << std::endl;
      }*/
    
  }
  //  bool operator<(const YZ& other) const { return (y() < other.y()) || (y() == other.y() && z() < other.z()); }
  int y() const { return this->second; }
  int z() const { return this->first; }
};

struct YZFull : public YZ {
  const bool valid;
 YZFull(const Module& m) : YZ(m), valid(true) {}
 YZFull(const XYZVector& v) : YZ(v), valid(true) {}
 YZFull(const XYZVector& v, const Module& m) : YZ(v, m), valid(true) {}
};

TPolyLine* drawMod();

template<class CoordType> class LineGetter {
  typedef typename CoordType::first_type CoordTypeX; 
  typedef typename CoordType::first_type CoordTypeY;
  CoordTypeX maxx_, minx_;
  CoordTypeY maxy_, miny_;
public:
  LineGetter() : maxx_(std::numeric_limits<CoordTypeX>::min()), minx_(std::numeric_limits<CoordTypeX>::max()), maxy_(std::numeric_limits<CoordTypeY>::min()), miny_(std::numeric_limits<CoordTypeY>::max()) {}
  CoordTypeX maxx() const { return double(maxx_)/Rounder::mmFraction; }
  CoordTypeX minx() const { return double(minx_)/Rounder::mmFraction; }
  CoordTypeY maxy() const { return double(maxy_)/Rounder::mmFraction; }
  CoordTypeY miny() const { return double(miny_)/Rounder::mmFraction; }
  TPolyLine* operator()(const Module& m) {
    std::set<CoordType> xy; // duplicate detection
    double x[] = {0., 0., 0., 0., 0.}, y[] = {0., 0., 0., 0., 0.};
    int j=0;
    for (int i=0; i<4; i++) {
      CoordType c(m.basePoly().getVertex(i), m);
      if (xy.insert(c).second == true) {
        x[j] = double(c.first)/Rounder::mmFraction;
        y[j++] = double(c.second)/Rounder::mmFraction;
      } 
      maxx_ = MAX(c.first, maxx_);
      minx_ = MIN(c.first, minx_);
      maxy_ = MAX(c.second, maxy_);
      miny_ = MIN(c.second, miny_);
    }
    if (j==4) { // close the poly line in case it's made of 4 distinct points, to get the closing line drawn
      x[j] = x[0]; 
      y[j++] = y[0];
    }
    return !g ? new TPolyLine(j, x, y) : drawMod();
  }
};





// ===============================================================================================
// Here be FRAMEGETTERS & FRAMESTYLERS
// Careful what you do
// ===============================================================================================


class IdMaker {
  static int id;
public:
  int nextId() const { return id++; }
  std::string nextString() const { 
    std::stringstream sid(""); 
    sid << nextId();
    return sid.str(); 
  }
protected:
  static const int nBinsZoom;
};



template<class CoordType> class FrameGetter : private IdMaker {
public:
  TH2C* operator()(double viewportX, double viewportY) const;
};



template<class CoordType> class SummaryFrameStyle {
  void drawEtaTicks(double maxL, double maxR, double tickDistance, double tickLength, double textDistance,
                    Style_t labelFont, Float_t labelSize, double etaStep, double etaMax, double etaLongLine) const;
public:
  void operator()(TH2C& frame, TCanvas& canvas, DrawerPalette&) const;

};

template<class CoordType>
struct HistogramFrameStyle {
  void operator()(TH2C& frame, TCanvas& canvas, DrawerPalette& palette) const {
    frame.Fill(1,1,0);
    frame.SetMaximum(palette.getMaxValue());
    frame.SetMinimum(palette.getMinValue());

    frame.Draw("colz");
    canvas.Update();
    palette.setFramePalette((TPaletteAxis*)frame.GetListOfFunctions()->FindObject("palette"));
    frame.Draw("AXIS colz");
  }
};


// ===============================================================================================
// Here be PLOTDRAWER MAIN CLASS
// Do not touch please
// ===============================================================================================



/// PlotDrawer is the class drawing module plots on user-supplied canvas. Alas, not visitor-enabled, since it was written before the new visitable geometry was conceived
/// Usage:
/// 1) Construct a PlotDrawer object: PlotDrawer<CoordType, ValueGetterType, StatType> drawer(viewportMaxX, viewportMaxY, valueGetter);
///    - CoordType is the coordinate class: XY, YZ or YZFull (for Z- and Z+ sections together)
///    - ValueGetterType obtains a value from modules to decide their color. Any functor or lambda taking a Module& and returning a double can be used here. 
///      Some ValueGetters are pre-defined. In case of user-defined ValueGetters, if possible use lambda or classes local to the instantiation of your PlotDrawer to avoid polluting this header with additional declarations
///    - StatType is the type of statistic to do on the module values obtained with the ValueGetters, in case two modules occupy the same map bin.
///      Default is NoStat, where values overwrite each other and the last one counts. Average, Max, Min, Sum are also available and custom statistics can be defined by the user.
///    - viewportMaxX, viewportMaxY specify the size of the viewport on the canvas. If set to 0, the viewport is automatically calculated to include all the modules. Default is 0 for both.
///    - valueGetter is an instance of a ValueGetterType that can be used by the user to pass a custom valueGetter. Default is valueGetter().
/// 2) Add modules to the internal maps, using either:
///    - void addModulesType(begin, end, moduleTypes); where the last argument moduleTypes can be the constants BARREL, ENDCAP or BARREL | ENDCAP, to restrict to one subdetector or both
///    - void addModules<ModuleValidator>(begin, end, isValid);  where the last argument isValid is of ModuleValidator type. 
///      A ModuleValidator is a lambda or functor taking a const Module& and returning a bool, used to decide whether a module should be included or not in the plot
/// 3) Draw the plot frame: void drawFrame<FrameStyleType>(canvas, frameStyle)
///   - FrameStyleType is the type of frame to draw. The predefined classes are SummaryFrameStyle (which draws eta lines) or HistogramFrameStyle (which draws the legend colour bar)
///   - canvas is the TCanvas to draw on. cd() is called automatically by the PlotDrawer
///   - frameStyle is the instance of a FrameStyleType class, which can be used in case of custom frame styles. Default is FrameStyleType<CoordType>()
/// 4) Draw the modules: void drawModules<DrawStyleType>(canvas, drawStyle)
///   - DrawStyleType is the style of module drawing. The predefined classes are ContourStyle (which only draws the contours of modules) and FillStyle (which draws solid modules).
///   - canvas is the TCanvas to draw on. cd() is called automatically by the PlotDrawer
///   - drawStyle is the instance of a DrawStyleType class, which can be used in case of custom draw styles. Default is DrawStyleType<CoordType>()

template<class CoordType, class ValueGetterType, class StatType = NoStat >
class PlotDrawer {
  typedef typename CoordType::first_type CoordTypeX;
  typedef typename CoordType::second_type CoordTypeY;
  CoordTypeX viewportMaxX_;
  CoordTypeY viewportMaxY_;
  ValueGetterType getValue;
  FrameGetter<CoordType> getFrame;
  LineGetter<CoordType> getLine;

  DrawerPalette palette_;

  std::map<CoordType, StatType*> bins_;
  std::map<CoordType, TPolyLine*> lines_;

public: 
  PlotDrawer(CoordTypeX viewportMaxX = 0, CoordTypeY viewportMaxY = 0, const ValueGetterType& valueGetter = ValueGetterType()) : viewportMaxX_(viewportMaxX), viewportMaxY_(viewportMaxY), getValue(valueGetter) {}

  ~PlotDrawer();

  template<template<class> class FrameStyleType> void drawFrame(TCanvas& canvas, const FrameStyleType<CoordType>& frameStyle = FrameStyleType<CoordType>());
  template<class DrawStyleType> void drawModules(TCanvas& canvas, const DrawStyleType& drawStyle = DrawStyleType());

  void add(const Module& m);
  template<class InputIterator> void addModulesType(InputIterator begin, InputIterator end, int moduleTypes = BARREL | ENDCAP);
  template<class ModuleValidator, class InputIterator> void addModules(InputIterator begin, InputIterator end, const ModuleValidator& isValid = ModuleValidator());

};


template<class CoordType, class ValueGetterType, class StatType> PlotDrawer<CoordType, ValueGetterType, StatType>::~PlotDrawer() {
  for (typename std::map<CoordType, StatType*>::iterator it = bins_.begin(); it != bins_.end(); ++it) {
    delete it->second;
    delete lines_[it->first];
  }
  bins_.clear();
  lines_.clear();
}

template<class CoordType, class ValueGetterType, class StatType>
template<template<class> class FrameStyleType>
void PlotDrawer<CoordType, ValueGetterType, StatType>::drawFrame(TCanvas& canvas, const FrameStyleType<CoordType>& frameStyle) {
  double minValue = std::numeric_limits<double>::max();
  double maxValue = 0;
  for (typename std::map<CoordType, StatType*>::const_iterator it = bins_.begin(); it != bins_.end(); ++it) {
    double value = it->second->get();
    minValue = value < minValue ? value : minValue;
    maxValue = value > maxValue ? value : maxValue; 
  }
  palette_.setMinMaxValues(minValue, maxValue);
  canvas.SetFillColor(kWhite);
  canvas.cd();
  viewportMaxX_ = viewportMaxX_ == 0 ? getLine.maxx()*1.1 : viewportMaxX_;  // in case the viewport coord is 0, auto-viewport mode is used and getLine is queried for the farthest X or Y it has registered
  viewportMaxY_ = viewportMaxY_ == 0 ? getLine.maxy()*1.1 : viewportMaxY_;
  TH2C* frame = getFrame(viewportMaxX_, viewportMaxY_);
  frameStyle(*frame, canvas, palette_);
}



template<class CoordType, class ValueGetterType, class StatType>
template<class DrawStyleType>
void PlotDrawer<CoordType, ValueGetterType, StatType>::drawModules(TCanvas& canvas, const DrawStyleType& drawStyle) {
  canvas.cd();
  for (typename std::map<CoordType, StatType*>::const_iterator it = bins_.begin(); it != bins_.end(); ++it) {
    StatType* bin = it->second;
    TPolyLine* line = lines_.at(it->first);
    drawStyle(*line, *bin, palette_);
  }
}

template<class CoordType, class ValueGetterType, class StatType>
void PlotDrawer<CoordType, ValueGetterType, StatType>::add(const Module& m) {
  CoordType c(m);
  if (!c.valid) return;
  if (bins_[c] == NULL) {
    bins_[c] = new StatType();
    lines_[c] = getLine(m);
  }
  double value = getValue(m);
  bins_[c]->fill(value);
}

template<class CoordType, class ValueGetterType, class StatType>
template<class InputIterator>
void PlotDrawer<CoordType, ValueGetterType, StatType>::addModulesType(InputIterator begin, InputIterator end, int moduleTypes) {
  for (InputIterator it = begin; it != end; ++it) {
    int subDet = (*it)->subdet();
    if (subDet & moduleTypes) {
      /*if ((*it)->center().Rho() == 33) {
	std::cout << "Rho = " << (*it)->center().Rho() << std::endl; 


for (int i=0; i<4; i++) {
        std::cout << "i = " << i << std::endl;
        std::cout << " (*it)->basePoly().getVertex(i).Z() = " <<  (*it)->basePoly().getVertex(i).Z() << std::endl;
/*std::cout << " (*it)->basePoly().getVertex(i).X() = " <<  (*it)->basePoly().getVertex(i).X() << std::endl;
       std::cout << " (*it)->basePoly().getVertex(i).Y() = " <<  (*it)->basePoly().getVertex(i).Y() << std::endl;
        std::cout << " sqrt = " <<  sqrt(pow((*it)->basePoly().getVertex(i).X(),2.) + pow((*it)->basePoly().getVertex(i).Y(),2.)) << std::endl;
        std::cout << " (*it)->basePoly().getVertex(i).Rho() = " <<  (*it)->basePoly().getVertex(i).Rho() << std::endl;
      }
      }*/


      add(**it); 
    }
  }
}

template<class CoordType, class ValueGetterType, class StatType>
template<class ModuleValidator, class InputIterator>
void PlotDrawer<CoordType, ValueGetterType, StatType>::addModules(InputIterator begin, InputIterator end, const ModuleValidator& isValid) {
  for (InputIterator it = begin; it != end; ++it) {
      if (isValid(**it)) add(**it);
  }
}
#endif

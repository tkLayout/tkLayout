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

#include <Palette.hh>
#include <Module.hh>





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
  // XY coordinates of the centre of module m.
 XY(const Module& m) : std::pair<int, int>(round(m.center().X()), round(m.center().Y())), valid(m.center().Z() >= 0) {}
  // XY coordinates of vector v.
 XY(const XYZVector& v) : std::pair<int, int>(round(v.X()), round(v.Y())), valid(v.Z() >= 0) {}
  // XY coordinates of vector v, in the (XY) plane passing by the center of module m.
 XY(const XYZVector& v, const Module& m) : XY(v) {}
  // bool operator<(const XY& other) const { return (x() < other.x()) || (x() == other.x() && y() < other.y()); }
  int x() const { return this->first; }
  int y() const { return this->second; }
};


struct YZ : public std::pair<int, int>, private Rounder {
  const bool valid;
  // RZ coordinates of the centre of module m, in the plane (RZ) defined by ((Z axis), moduleCenter).
 YZ(const Module& m) : std::pair<int,int>(round(m.center().Z()), round(m.center().Rho())), valid(m.center().Z() >= 0) {}
  // RZ coordinates of vector v, in the plane (RZ) defined by ((Z axis), v).
 YZ(const XYZVector& v) : std::pair<int, int>(round(v.Z()), round(v.Rho())), valid(v.Z() >= 0) {}
  // RZ coordinates of vector v, in the plane (RZ) defined by ((Z axis), moduleCenter).
 YZ(const XYZVector& v, const Module& m) : valid(v.Z() >= 0) {
    this->first = round(v.Z());

    // This calculates the projection of vector v into plane ((Z axis), moduleCenter).
    XYZVector vProjected;
    XYZVector z(0., 0., 1.);

    XYZVector basePolyCenter = m.basePoly().getCenter();
    XYZVector normal = crossProduct(z, basePolyCenter); // normal of plane ((Z axis), moduleCenter).
    normal = normal.Unit(); // normalize.

    vProjected = v - v.Dot(normal) * normal; // Calculate projected vector. TO DO : Introduce this as a method !!  
    this->second = round(vProjected.Rho()); // Take the Rho() of the projected vector.   
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
  bool isContour = false;
  LineGetter() : maxx_(std::numeric_limits<CoordTypeX>::min()), minx_(std::numeric_limits<CoordTypeX>::max()), maxy_(std::numeric_limits<CoordTypeY>::min()), miny_(std::numeric_limits<CoordTypeY>::max()) {}
  CoordTypeX maxx() const { return double(maxx_)/Rounder::mmFraction; }
  CoordTypeX minx() const { return double(minx_)/Rounder::mmFraction; }
  CoordTypeY maxy() const { return double(maxy_)/Rounder::mmFraction; }
  CoordTypeY miny() const { return double(miny_)/Rounder::mmFraction; }
  TPolyLine* operator()(const Module& m) {
    if (!isContour) {
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
    } else {
      // well... it's a contour! let's hack a drawing here

      int contourSize = m.decorated().contour().size();
      if (contourSize==0) return nullptr;
     
      // Our local axes in global coordinates
      XYZVector ey = m.basePoly().getVertex(0) - m.basePoly().getVertex(1) ;
      XYZVector ex = m.basePoly().getVertex(2) - m.basePoly().getVertex(1) ;
      XYZVector center = m.center();
      ex = ex / sqrt(ex.Mag2());
      ey = ey / sqrt(ey.Mag2());
      double x[contourSize+1];
      double y[contourSize+1];
      std::set<CoordType> xy; // duplicate detection

      int j=0;
      for (int i=0; i<contourSize; i++) {
	const XYZVector& contourLocal = m.decorated().contour().at(i);
	// std::cerr << contourLocal.X() << "," << contourLocal.Y() << " ";
	XYZVector contourGlobal = ex * contourLocal.X() + ey * contourLocal.Y() + center;
	
	CoordType c(contourGlobal, m);
	if (xy.insert(c).second == true) {
	  x[j] = double(c.first)/Rounder::mmFraction;
	  y[j++] = double(c.second)/Rounder::mmFraction;
	}
	maxx_ = MAX(c.first, maxx_);
	minx_ = MIN(c.first, minx_);
	maxy_ = MAX(c.second, maxy_);
	miny_ = MIN(c.second, miny_);
      }
      x[j] = x[0];
      y[j++] = y[0];

      // std::cerr << std::endl;
      return new TPolyLine(j, x, y);
    }
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
  void drawEtaTicks(double maxL, double maxR,
		    double tickDistanceRRatio, double tickLengthRRatio, double textDistanceRRatio,
		    double tickDistanceLRatio, double tickLengthLRatio, double textDistanceLRatio,
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

    frame.Draw("colz HIST");
    canvas.Update();
    palette.setFramePalette((TPaletteAxis*)frame.GetListOfFunctions()->FindObject("palette"));
    frame.Draw("AXIS HIST colz");
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
/// 5) Draw the modules with outer contour: void drawModuleContours<DrawStyleType>(canvas, drawStyle)
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
  LineGetter<CoordType> getContour;

  DrawerPalette palette_;

  std::map<CoordType, StatType*> bins_;
  std::map<CoordType, TPolyLine*> lines_;
  std::map<CoordType, TPolyLine*> contour_;

public:
   PlotDrawer(CoordTypeX viewportMaxX = 0, CoordTypeY viewportMaxY = 0, const ValueGetterType& valueGetter = ValueGetterType()) : viewportMaxX_(viewportMaxX), viewportMaxY_(viewportMaxY), getValue(valueGetter) { getContour.isContour = true; }

  ~PlotDrawer();

  template<template<class> class FrameStyleType> void drawFrame(TCanvas& canvas, const FrameStyleType<CoordType>& frameStyle = FrameStyleType<CoordType>());
  template<class DrawStyleType> void drawModules(TCanvas& canvas, const DrawStyleType& drawStyle = DrawStyleType());
  template<class DrawStyleType> void drawModuleContours(TCanvas& canvas, const DrawStyleType& drawStyle = DrawStyleType());

  void add(const Module& m);
  template<class InputIterator> void addModulesType(InputIterator begin, InputIterator end, int moduleTypes = BARREL | ENDCAP);
  void addModules(const Visitable& structure);
  void addModulesType(const Visitable& structure, int moduleTypes = BARREL | ENDCAP);
  template<class ModuleValidator> void addModules(const Visitable& structure, const ModuleValidator& isValid = ModuleValidator());
  template<class ModuleValidator, class InputIterator> void addModules(InputIterator begin, InputIterator end, const ModuleValidator& isValid = ModuleValidator());
};


template<class CoordType, class ValueGetterType, class StatType> PlotDrawer<CoordType, ValueGetterType, StatType>::~PlotDrawer() {
  for (typename std::map<CoordType, StatType*>::iterator it = bins_.begin(); it != bins_.end(); ++it) {
    delete it->second;
    delete lines_[it->first];
    if (contour_[it->first]) delete contour_[it->first];
  }
  bins_.clear();
  lines_.clear();
  contour_.clear();
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
template<class DrawStyleType>
  void PlotDrawer<CoordType, ValueGetterType, StatType>::drawModuleContours(TCanvas& canvas, const DrawStyleType& drawStyle) {
  canvas.cd();
  for (typename std::map<CoordType, StatType*>::const_iterator it = bins_.begin(); it != bins_.end(); ++it) {
    StatType* bin = it->second;
    TPolyLine* contour = contour_.at(it->first);
    if (contour) drawStyle(*contour, *bin, palette_);
  }
}

template<class CoordType, class ValueGetterType, class StatType>
void PlotDrawer<CoordType, ValueGetterType, StatType>::add(const Module& m) {
  CoordType c(m);
  if (!c.valid) return;
  if (bins_[c] == NULL) {
    bins_[c] = new StatType();
    lines_[c] = getLine(m);
    contour_[c] = getContour(m);
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
      add(**it);
    }
  }
}

template<class CoordType, class ValueGetterType, class StatType>
  void PlotDrawer<CoordType, ValueGetterType, StatType>::addModules(const Visitable& structure) {
  class ModuleVisitor : public ConstGeometryVisitor {
  private:
    PlotDrawer<CoordType, ValueGetterType, StatType>* pd_;
  public:
    ModuleVisitor(PlotDrawer<CoordType, ValueGetterType, StatType>* pd) { pd_ = pd; };
    void visit(const Module& m) { pd_->add(m); }
  };
  ModuleVisitor v(this);
  structure.accept(v);
}

template<class CoordType, class ValueGetterType, class StatType>
  void PlotDrawer<CoordType, ValueGetterType, StatType>::addModulesType(const Visitable& structure, int moduleTypes) {
  class ModuleVisitor : public ConstGeometryVisitor {
  private:
    PlotDrawer<CoordType, ValueGetterType, StatType>* pd_;
    int moduleTypes_;
  public:
    ModuleVisitor(PlotDrawer<CoordType, ValueGetterType, StatType>* pd, int modType) {
      pd_ = pd;
      moduleTypes_ = modType;
    };
    void visit(const Module& m) {
      int subDet = m.subdet();
      if (subDet & moduleTypes_) pd_->add(m);
    }
  };
  ModuleVisitor v(this, moduleTypes);
  structure.accept(v);
}

template<class CoordType, class ValueGetterType, class StatType>
template<class ModuleValidator>
void PlotDrawer<CoordType, ValueGetterType, StatType>::addModules(const Visitable& structure, const ModuleValidator& isValid) {
  class ModuleVisitor : public ConstGeometryVisitor {
  private:
    PlotDrawer<CoordType, ValueGetterType, StatType>* pd_;
    const ModuleValidator& isValid_;
  public:
    ModuleVisitor(PlotDrawer<CoordType, ValueGetterType, StatType>* pd, const ModuleValidator& isValid) : pd_(pd), isValid_(isValid) {};
    void visit(const Module& m) {
      int subDet = m.subdet();
      if (isValid_(m)) pd_->add(m);
    }
  };
  ModuleVisitor v(this, isValid);
  structure.accept(v);
}

template<class CoordType, class ValueGetterType, class StatType>
template<class ModuleValidator, class InputIterator>
void PlotDrawer<CoordType, ValueGetterType, StatType>::addModules(InputIterator begin, InputIterator end, const ModuleValidator& isValid) {
  for (InputIterator it = begin; it != end; ++it) {
      if (isValid(**it)) add(**it);
  }
}
#endif

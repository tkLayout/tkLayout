#ifndef INCLUDE_PLOT_DRAWER_H
#define INCLUDE_PLOT_DRAWER_H

#include <utility>
#include <map>
#include <set>
#include <iostream>

#include <TPolyLine.h>
#include <TPaletteAxis.h>
#include <TList.h>
#include <TH2C.h>
#include <TText.h>
#include <TLine.h>
#include <TCanvas.h>

#include <Palette.h>
#include <DetectorModule.h>

/*
 * PlotDrawer is the class drawing tracker modules on user-supplied canvas. Alas, not visitor-enabled, since it was written before the new visitable geometry was conceived
 * Usage:
 * 1) Construct a PlotDrawer object: PlotDrawer<CoordType, ValueGetterType, StatType> drawer(viewportMaxX, viewportMaxY, valueGetter);
 *    - CoordType is the coordinate class: XY, RZ or RZFull (for Z- and Z+ sections together)
 *    - ValueGetterType obtains a value from modules to decide their color. Any functor or lambda taking a DetectorModule& and returning a double can be used here.
 *      Some ValueGetters are pre-defined. In case of user-defined ValueGetters, if possible use lambda or classes local to the instantiation of your PlotDrawer to avoid polluting this header with additional declarations
 *    - StatType is the type of statistic to do on the module values obtained with the ValueGetters, in case two modules occupy the same map bin.
 *      Default is NoStat, where values overwrite each other and the last one counts. Average, Max, Min, Sum are also available and custom statistics can be defined by the user.
 *    - viewportMaxX, viewportMaxY specify the size of the viewport on the canvas. If set to 0, the viewport is automatically calculated to include all the modules. Default is 0 for both.
 *    - valueGetter is an instance of a ValueGetterType that can be used by the user to pass a custom valueGetter. Default is valueGetter().
 * 2) Add modules to the internal maps, using either:
 *    - void addModulesType(begin, end, moduleTypes); where the last argument moduleTypes can be the constants BARREL, ENDCAP or BARREL | ENDCAP, to restrict to one subdetector or both
 *    - void addModules<ModuleValidator>(begin, end, isValid);  where the last argument isValid is of ModuleValidator type.
 *      A ModuleValidator is a lambda or functor taking a const DetectorModule& and returning a bool, used to decide whether a module should be included or not in the plot
 * 3) Draw the plot frame: void drawFrame<FrameStyleType>(canvas, frameStyle)
 *   - FrameStyleType is the type of frame to draw. The predefined classes are TicksFrameStyle (which draws eta lines) or HistogramFrameStyle (which draws the legend colour bar)
 *   - canvas is the TCanvas to draw on. cd() is called automatically by the PlotDrawer
 *   - frameStyle is the instance of a FrameStyleType class, which can be used in case of custom frame styles. Default is FrameStyleType<CoordType>()
 * 4) Draw the modules: void drawModules<DrawStyleType>(canvas, drawStyle)
 *   - DrawStyleType is the style of module drawing. The predefined classes are ContourStyle (which only draws the contours of modules) and FillStyle (which draws solid modules).
 *   - canvas is the TCanvas to draw on. cd() is called automatically by the PlotDrawer
 *   - drawStyle is the instance of a DrawStyleType class, which can be used in case of custom draw styles. Default is DrawStyleType<CoordType>()
 */

// ========================================================================================
// Here be STATISTICS
// If possible, additional statistics should be local to the plot drawer instantiation
// ========================================================================================

//! Keep statistics
class NoStat {
  double m_value;
public:
  NoStat() : m_value(0) {}
  void   fill(double value) { m_value = value; }
  double get() const { return m_value; }
};

//! Keep average value
class Average {
  int    m_counts;
  double m_total;
public:
  Average() : m_counts(0), m_total(0) {}
  void   fill(double value) {
    m_total += value;
    m_counts+= 1;
  }
  double get() const { return m_total/m_counts; }
};

// Calculate maximum value
class Max {
  double m_max;
public:
  Max() : m_max(-std::numeric_limits<double>::max()) {}
  void   fill(double value) { m_max = (value > m_max) ? value : m_max; }
  double get() const { return m_max; }
};

//! Calculate minimum value
class Min {
  double m_min;
public:
  Min() : m_min(std::numeric_limits<double>::max()) {}
  void fill(double value) { m_min = (value < m_min) ? value : m_min; }
  double get() const { return m_min; }
};

//! Calculate sum
class Sum {
  double m_sum;
public:
  Sum() : m_sum(0) {}
  void   fill(double value) { m_sum += value; }
  double get() const { return m_sum; }
};

// ==============================================================================================
// Here be VALUEGETTERS
// If possible, additional value getters should be local to the plot drawer instantiation
// ==============================================================================================

//! Access arbitrary module method
template<class RetType, RetType (DetectorModule::*ModuleMethod)() const>
struct Method {
  double operator()(const DetectorModule& m) const { return (double)(m.*ModuleMethod)(); }
};

//! Assign automatically color to a module based on its type
struct TypeAutoColor { // Auto-assign colors
  std::set<std::string> colorSet_;
  double operator()(const DetectorModule& m) {
    std::pair<std::set<std::string>::iterator, bool> it = colorSet_.insert(m.moduleType());
    return Palette::color(std::distance(colorSet_.begin(), it.first)+1);
  }
};

//! Assign color to a module based on plotColor variable defined in configuration file
struct Type { // Module-maintained color
  double operator()(const DetectorModule& m) {
    return Palette::color(m.plotColor());
  }
};

//! Get module central Z position
struct CoordZ {
  double operator()(const DetectorModule& m) { return m.center().Z(); }
};

// =============================================================================================
// Here be DRAWSTYLES 
// If possible, additional draw styles should be local to the plot drawer instantiation
// =============================================================================================

//! Draw style class providing ROOT color palette
class DrawerPalette {
  double m_minValue, m_maxValue;
  TPaletteAxis* m_framePalette;
public:
  DrawerPalette() : m_framePalette(nullptr), m_minValue(0), m_maxValue(0) {}
  void   setMinMaxValues(double minValue, double maxValue) {
    m_minValue = minValue * 1.001 - maxValue * 0.001;
    m_maxValue = maxValue * 1.001 - minValue * 0.001;
  }
  double getMinValue() const { return m_minValue; }
  double getMaxValue() const { return m_maxValue; }
  void   setFramePalette(TPaletteAxis* framePalette) {
    m_framePalette = framePalette;
  }
  int    getColor(double value) const {
    return m_framePalette ? m_framePalette->GetValueColor(value) : value;
  }
};

//! Draws modules as solid lines with defined color
struct FillStyle {
  template <class StatType> void operator()(TPolyLine& line, StatType& bin, const DrawerPalette& palette) const {
    line.SetFillColor(palette.getColor(bin.get()));
    line.DrawPolyLine(line.GetN(),line.GetX(),line.GetY(),"f");
  }
};

//! Draw modules as contour lines with defined color & width of contour line
class ContourStyle {
  const int m_lineWidth;
public:
  ContourStyle(int lineWidth = 2)  : m_lineWidth(lineWidth) {}
  template <class StatType> void operator()(TPolyLine& line, StatType& bin, const DrawerPalette& palette) const {
    line.SetLineColor(palette.getColor(bin.get()));
    line.SetLineWidth(m_lineWidth);
    line.DrawPolyLine(line.GetN(),line.GetX(),line.GetY());
  }
};


// ===============================================================================================
// Here be VALIDATORS 
// If possible, additional module validators should be local to the plot drawer instantiation
// ===============================================================================================

//! Check module type: ModuleSubdetector { BARREL = 1, ENDCAP = 2 };
template<const int SubdetType>
struct CheckType {
  bool operator()(const DetectorModule& m) const { return m.subdet() & SubdetType; }
};

template<const int PhiIndex>
struct CheckPhiIndex {
  bool operator()(const DetectorModule& m) const { return m.posRef().phi == PhiIndex; }
};

// ===============================================================================================
// Here be COORDINATE CLASSES and LINEGETTER
// Careful what you do
// ===============================================================================================

//! Helper struct to inherit round method
struct Rounder {
  static const int mmFraction = 1000;  // Drawing with 1 micrometer precision
  int round(double x) { return floor(x*mmFraction+0.5); }
};

//! Get module XY coordinates rounded with micron precision
struct XY : public std::pair<int, int>, private Rounder {
  const bool valid;
  XY(const DetectorModule& m) : std::pair<int, int>(round(m.center().X()), round(m.center().Y())), valid(m.center().Z() >= 0) {}
  XY(const XYZVector& v)      : std::pair<int, int>(round(v.X()), round(v.Y())), valid(v.Z() >= 0) {}
  // bool operator<(const XY& other) const { return (x() < other.x()) || (x() == other.x() && y() < other.y()); }
  int x() const { return this->first; }
  int y() const { return this->second; }
};

//! Get module RZ (R=>0, Z>=0) coordinates rounded with micron precision
struct RZ : public std::pair<int, int>, private Rounder {
  const bool valid;
  RZ(const DetectorModule& m) : std::pair<int,int>(round(m.center().Z()), round(m.center().Rho())), valid(m.center().Z() >= 0) {}
  RZ(const XYZVector& v)      : std::pair<int,int>(round(v.Z()), round(v.Rho())), valid(v.Z() >= 0) {}
  //  bool operator<(const YZ& other) const { return (y() < other.y()) || (y() == other.y() && z() < other.z()); }
  int r() const { return this->second; }
  int z() const { return this->first; }
};

//! Get module RZ (R=>0, Z full) coordinates rounded with micron precision
struct RZFull : public RZ {
  const bool valid;
  RZFull(const DetectorModule& m) : RZ(m), valid(true) {}
  RZFull(const XYZVector& v) : RZ(v), valid(true) {}
};

//! Based on module coordinate type <XY>, <RZ>, <RZFull> get module parameters: x,y ; r,z(>=0) ; r,z -> module expressed as PolyLine
template<class CoordType> class LineGetter {
  typedef typename CoordType::first_type CoordTypeX;
  typedef typename CoordType::first_type CoordTypeY;
  CoordTypeX m_xMax, m_xMin;
  CoordTypeY m_yMax, m_yMin;
public:
  LineGetter() : m_xMax(std::numeric_limits<CoordTypeX>::min()), m_xMin(std::numeric_limits<CoordTypeX>::max()), m_yMax(std::numeric_limits<CoordTypeY>::min()), m_yMin(std::numeric_limits<CoordTypeY>::max()) {}
  CoordTypeX xMax() const { return double(m_xMax)/Rounder::mmFraction; }
  CoordTypeX xMin() const { return double(m_xMin)/Rounder::mmFraction; }
  CoordTypeY yMax() const { return double(m_yMax)/Rounder::mmFraction; }
  CoordTypeY yMin() const { return double(m_yMin)/Rounder::mmFraction; }

  TPolyLine* operator()(const DetectorModule& m) {
    std::set<CoordType> xy; // duplicate detection
    double x[] = {0., 0., 0., 0., 0.}, y[] = {0., 0., 0., 0., 0.};
    int j=0;
    for (int i=0; i<4; i++) {
      CoordType c(m.basePoly().getVertex(i));
      if (xy.insert(c).second == true) {
        x[j] = double(c.first)/Rounder::mmFraction;
        y[j++] = double(c.second)/Rounder::mmFraction;
      }
      m_xMax = MAX(c.first, m_xMax);
      m_xMin = MIN(c.first, m_xMin);
      m_yMax = MAX(c.second, m_yMax);
      m_yMin = MIN(c.second, m_yMin);
    }
    if (j==4) { // close the poly line in case it's made of 4 distinct points, to get the closing line drawn
      x[j]   = x[0];
      y[j++] = y[0];
    }
    return (new TPolyLine(j, x, y));
  }
};

// ===============================================================================================
// Here be FRAMEGETTERS & FRAMESTYLERS
// Careful what you do
// ===============================================================================================

//! Frame (histogram) base class -> assign unique id to plot drawer & name
class FrameGetterBase {
  static int s_id;
public:
  int nextId() const { return s_id++; }
  std::string nextString() const {
    std::stringstream sid("");
    sid << nextId();
    return sid.str();
  }
protected:
  static const int s_nBinsZoom;
};

//! Frame (histogram & canvas) class for given coordinate type: XY, RZ, RZFull in given range viewportX, viewportY
template<class CoordType> class FrameGetter : private FrameGetterBase {
public:
  TH2C* operator()(double viewportX, double viewportY) const;
};

//! Frame style with ticks to be drawn for given coordinate type: XY, RZ, RZFull -> used to plot eta range
template<class CoordType> class TicksFrameStyle {
  void drawEtaTicks(double maxZ, double maxR, double tickDistance, double tickLength, double textDistance,
                    Style_t labelFont, Float_t labelSize, double etaStep, double etaMax, double etaLongLine) const;
  void drawEtaTicks(double maxZ, double maxR, double tickDistance, double tickScaleFactor, double textScaleFactor,
                    Style_t labelFont, Float_t labelSize, double etaStepShort, double etaStepLong, double etaMax, double etaLongLineI, 
                    double etaLongLineII, bool zPlsMin=false) const;
  
public:
  void operator()(TH2C& frame, TCanvas& canvas, DrawerPalette& palette, bool isPixelType) const;

};

//! Frame style used to draw histograms for given coordinate type: XY, RZ, RZFull
template<class CoordType> struct HistogramFrameStyle {
  void operator()(TH2C& frame, TCanvas& canvas, DrawerPalette& palette, bool isPixelType) const {
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

//! Main plot drawer class:
//! Define coordinate type (coordinate projection) to be used to draw individual modules: XY, RZ (Z>=0), RZFull (full Z axis)
//! Define value getter type -> to define what information about each module will be kept (Type = color, etc.)
template<class CoordType, class ValueGetterType, class StatType = NoStat > class PlotDrawer {

  typedef typename CoordType::first_type  CoordTypeX;
  typedef typename CoordType::second_type CoordTypeY;

  CoordTypeX m_viewportMaxX;
  CoordTypeY m_viewportMaxY;

  ValueGetterType        getValue; //!> Information to be kept about each module
  FrameGetter<CoordType> getFrame; //!> Histogram, into which all modules are drawn
  LineGetter<CoordType>  getLine;  //!> Line, representing given module as an object

  DrawerPalette m_palette;

  std::map<CoordType, StatType*>  m_bins;  // Information kept about each module: color
  std::map<CoordType, TPolyLine*> m_lines; // Graphical line representing the given module

public: 

  //! Constructor
  PlotDrawer(CoordTypeX viewportMaxX = 0, CoordTypeY viewportMaxY = 0, const ValueGetterType& valueGetter = ValueGetterType()) : m_viewportMaxX(viewportMaxX), m_viewportMaxY(viewportMaxY), getValue(valueGetter) {}

  //! Destructor
  ~PlotDrawer() {

    for (typename std::map<CoordType, StatType*>::iterator it = m_bins.begin(); it != m_bins.end(); ++it) {
      delete it->second;
      delete m_lines[it->first];
    }
    m_bins.clear();
    m_lines.clear();
  }

  //! Draw frame, i.e. histogram, using given frame style (a way how canvas is drawn): TicksFrameStyle, HistogramFrameStyle -> the frame style is templated by coordinate type, i.e. canvas in XY, RZ, RZFull projection
  template<template<class> class FrameStyleType> void drawFrame(TCanvas& canvas, bool isPixelType=false, const FrameStyleType<CoordType>& frameStyle = FrameStyleType<CoordType>()) {

    double minValue = std::numeric_limits<double>::max();
    double maxValue = 0;

    // Define color scale to be used to draw modules
    for (typename std::map<CoordType, StatType*>::const_iterator it = m_bins.begin(); it != m_bins.end(); ++it) {
      double value = it->second->get();
      minValue = value < minValue ? value : minValue;
      maxValue = value > maxValue ? value : maxValue;
    }
    m_palette.setMinMaxValues(minValue, maxValue);

    // Set canvas properties
    canvas.SetFillColor(kWhite);
    canvas.cd();
    m_viewportMaxX = m_viewportMaxX == 0 ? getLine.xMax()*1.1 : m_viewportMaxX;  // in case the viewport coord is 0, auto-viewport mode is used and getLine is queried for the farthest X or Y that has been registered
    m_viewportMaxY = m_viewportMaxY == 0 ? getLine.yMax()*1.1 : m_viewportMaxY;

    // Draw histogram containing all modules in given projection: CoordType - XY, RZ, RZFull
    TH2C* frame = getFrame(m_viewportMaxX, m_viewportMaxY);
    frameStyle(*frame, canvas, m_palette, isPixelType);
  }

  //! Draw all modules that have been added to the plot drawer to the given canvas using internal frame, i.e. histogram. Specify draw style to be used: ContourStyle or FillStyle
  template<class DrawStyleType> void drawModules(TCanvas& canvas, const DrawStyleType& drawStyle = DrawStyleType()) {
    canvas.cd();
    for (typename std::map<CoordType, StatType*>::const_iterator it = m_bins.begin(); it != m_bins.end(); ++it) {

      // Get module as info container (first) and graphical object (line)
      StatType* bin   = it->second;
      TPolyLine* line = m_lines.at(it->first);

      // Draw using Contour or FillStyle
      drawStyle(*line, *bin, m_palette);
    }
  }

  //! Add detector module to the plot drawer, to be drawn in the future
  void add(const DetectorModule& m) {

    // Get projection of given module
    CoordType c(m);

    // Is module within projection criteria
    if (!c.valid) return;

    // Set module as a container & graphical object: line
    if (m_bins[c] == nullptr) {
      m_bins[c]  = new StatType();
      m_lines[c] = getLine(m);
    }

    // Set value to be filled, e.g. color
    double value = getValue(m);
    m_bins[c]->fill(value);
  }

  //! Add a set of modules to be drawn (define the set using iterators). In order to specify which module types will be added to plot drawer use BARREL and/or ENDCAP
  template<class InputIterator> bool addModulesType(InputIterator begin, InputIterator end, int moduleTypes = BARREL | ENDCAP) {

    bool foundModule = false;

    // If module of given type, add it to be drawn later
    for (InputIterator it = begin; it != end; ++it) {
      int subDet = (*it)->subdet();
      if (subDet & moduleTypes) {
        foundModule = true;
        add(**it);
      }
    }
    return foundModule;
  }

  //! Add a set of modules to be drawn (define the set using iterators). Define validator (boolean function) to be used if module should be drawn
  template<class ModuleValidator, class InputIterator> void addModules(InputIterator begin, InputIterator end, const ModuleValidator& isValid = ModuleValidator()) {

    // If module valid, add it to be drawn later
    for (InputIterator it = begin; it != end; ++it) {
        if (isValid(**it)) add(**it);
    }
  }


}; // Class

#endif /* INCLUDE_PLOT_DRAWER_H */

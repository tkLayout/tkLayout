#ifndef INCLUDE_TRACKER_H_
#define INCLUDE_TRACKER_H_

// System include files
#include <vector>
#include <string>
#include <memory>

#include <Barrel.h>
#include <DetectorModule.h>
#include <Endcap.h>
#include "Property.h"
#include <SupportStructure.h>
#include "Visitor.h"
#include "Visitable.h"

// Forward
class HierarchicalNameVisitor;
class ModulesSetVisitor;

// Typedefs
typedef PtrVector<SupportStructure> SupportStructures;

/*
 * @class Tracker
 * @details Tracker class holds information about individual detection systems - trackers: pixel tracker, strip tracker,
 * forward tracker, etc... Information read in from an xml tree (boost property tree) passed as an argument to a constructor.
 * The whole tracker structure is built via Tracker build() method. This method builds recursively all individual
 * subdetectors: barrels -> layers -> rods -> modules, endcaps -> disks -> rings -> modules, etc.
 */
class Tracker : public PropertyObject, public Buildable, public Identifiable<string>, Clonable<Tracker>, Visitable {

 public:
  
  //! Constructor - parse geometry config file using boost property tree & read-in Barrel, Endcap & Support nodes
  Tracker(const PropertyTree& treeProperty);

  //! Destructor
  ~Tracker();

  //! Build recursively individual subdetector systems: Barrels, Endcaps
  void build();

  //! Setup: link functions to various tracker related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Return tracker barrels
  const Barrels& barrels() const { return m_barrels; }

  //! Update tracker barrels (with services etc.)
  Barrels& updateBarrels() { return m_barrels; }

  //! Return tracker endcaps
  const Endcaps& endcaps() const { return m_endcaps; }

  //! Return tracker supports
  const SupportStructures& supports() const {return m_supportStructures;}

  //! Return all tracker modules
  const std::vector<DetectorModule*>& modules() const;

  //Modules& modules() { return modulesSetVisitor_.modules(); }

  //! GeometryVisitor pattern -> tracker visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> tracker visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  ReadonlyProperty<double, Computable> minR;        //!< Minimum radius of a tracker
  ReadonlyProperty<double, Computable> maxR;        //!< Maximum radius of a tracker
  ReadonlyProperty<double, Computable> maxZ;        //!< Maximum Z position of a tracker
  ReadonlyProperty<double, Computable> minEta;      //!< Minimum tracker eta
  ReadonlyProperty<double, Computable> maxEta;      //!< Maximum tracker eta
  ReadonlyProperty<double, Default>    etaCut;      //!< Eta cut @ which automated building of a tracker is stopped
  ReadonlyProperty<bool,   Default>    isPixelType; //!< Is tracker built of pixels
  ReadonlyProperty<bool,   Default>    servicesForcedUp;
  ReadonlyProperty<bool,   Default>    skipAllServices;
  ReadonlyProperty<bool,   Default>    skipAllSupports;

 private:

  //! Copy constructor
  Tracker(const Tracker&) = default;

  Barrels           m_barrels;              //!< Barrel components of tracker
  Endcaps           m_endcaps;              //!< Endcap components of tracker
  SupportStructures m_supportStructures;    //!< Tracker supports

  ModulesSetVisitor*       m_modulesSetVisitor; //!< Visitor pattern setting all modules assigned to the tracker
  HierarchicalNameVisitor* m_cntNameVisitor;    //!< Visitor pattern getting names of all tracker components

  PropertyNode<string>       m_barrelNode;  //!< Property tree nodes for barrel (included geometry config file)
  PropertyNode<string>       m_endcapNode;  //!< Property tree nodes for endcap (included geometry config file)
  PropertyNodeUnique<string> m_supportNode; //!< Property tree nodes for support (included geometry config file)

  MultiProperty<set<string>, ','> m_containsOnly;

}; // Class

/*
 * Helper class: Module visitor -> get all tracker modules
 */
class ModulesSetVisitor : public GeometryVisitor {

 public:

  // Typedef
  typedef std::vector<DetectorModule*> Modules;

  // Destructor
  virtual ~ModulesSetVisitor() {};

  // Visit pattern - main method
  void visit(DetectorModule& m) override { m_modules.push_back(&m); }

  // Return modules
  Modules& modules()             { return m_modules; }
  const Modules& modules() const { return m_modules; }

  // Iterators
  Modules::iterator begin() { return m_modules.begin(); }
  Modules::iterator end()   { return m_modules.end(); }
  Modules::const_iterator begin() const { return m_modules.begin(); }
  Modules::const_iterator end()   const { return m_modules.end(); }

 private:
  Modules m_modules;
}; // Helper class

/*
 * Helper class: Name visitor -> get names of all tracker components
 */
class HierarchicalNameVisitor : public GeometryVisitor {

 public:

  virtual ~HierarchicalNameVisitor() {}
  void visit(Barrel& b);
  void visit(Endcap& e);
  void visit(Layer& l);
  void visit(Disk& d);
  void visit(RodPair& r);
  void visit(Ring& r);
  void visit(DetectorModule& m);
  void visit(BarrelModule& m);
  void visit(EndcapModule& m);

 private:

  int cntId = 0;
  string cnt;
  int c1, c2;
}; // Helper class

#endif /* INCLUDE_TRACKER_H_ */

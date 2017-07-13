#ifndef CABLING_CONSTANTS_HH
#define CABLING_CONSTANTS_HH

#include <math.h>
#include <string>
#include "global_funcs.hh"


// DESIGN
// Maximum number of modules per bundle
static const int cabling_maxNumModulesPerBundle = 12;
// Maximum number of bundles per cable
static const int cabling_maxNumBundlesPerCable = 6;


// PHI SLICES
// Size of the phi slices used
static const double cabling_nonantWidth = 2. * M_PI / 9.;                          // 40°
static const double cabling_semiNonantWidth = 2. * M_PI / 18.;                     // 20°
static const double cabling_endcapStripStripPhiRegionWidth = 2. * M_PI / 27.;      // 13.333°

// Offset are sometimes used to set the phi slices.
// This has been tried to be reduced to the bare minimum: only 2 hardcoded constants :)
static const double cabling_tedd1StripStripPhiRegionStart = -0.55 * M_PI / 180.;
static const double cabling_tedd2StripStripPhiRegionStart = -0.001 * M_PI / 180.;


// ROUNDING
static const double cabling_roundingTolerance = 1.E-4;


// GEOMETRY NAMES
// It is unavailable by design to use the geometry blocks names for the cabling map.
// At least, they are placed here to ease renaming if needed.
static const std::string cabling_tbps = "TBPS";
static const std::string cabling_tb2s = "TB2S";
static const std::string cabling_tedd1 = "TEDD_1";
static const std::string cabling_tedd2 = "TEDD_2";


// NEGATIVE CABLING SIDE
static const std::string cabling_negativePrefix = "neg";


// CABLING TYPE
// PS or SS stands for the module type.
// 10G or 5G stands for the speed in the optical fibers (Gb/s).
// A or B is an extra distinction between groups of modules, used in TEDD.
// Modules associated to different cabling types should not be mixed with each other.
// Modules with the same cabling type have to be conected to the same DTC type (PS10G, PS5G, or SS).
enum Category { UNDEFINED, PS10G, PS10GA, PS10GB, PS5G, SS };


#endif  // CABLING_CONSTANTS_HH

#ifndef CABLING_CONSTANTS_HH
#define CABLING_CONSTANTS_HH

#include <math.h>
#include <string>
#include "global_funcs.hh"


static const int cabling_maxNumModulesPerBundle = 12;
static const int cabling_maxNumBundlesPerCable = 6;


static const double cabling_nonantWidth = 2. * M_PI / 9.;
static const double cabling_semiNonantWidth = 2. * M_PI / 18.;
static const double cabling_endcapStripStripPhiRegionWidth = 2. * M_PI / 27.;

static const double cabling_tedd1StripStripPhiRegionStart = -0.55 * M_PI / 180.;
static const double cabling_tedd2StripStripPhiRegionStart = -0.001 * M_PI / 180.;

static const double cabling_roundingTolerance = 1.E-4;


static const std::string cabling_tbps = "TBPS";
static const std::string cabling_tb2s = "TB2S";
static const std::string cabling_tedd1 = "TEDD_1";
static const std::string cabling_tedd2 = "TEDD_2";

enum Category { UNDEFINED, PS10G, PS5G, PS5GA, PS5GB, SS };



#endif  // CABLING_CONSTANTS_HH

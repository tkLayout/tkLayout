#ifndef INNER_CABLING_CONSTANTS_HH
#define INNER_CABLING_CONSTANTS_HH

#include <math.h>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include "global_funcs.hh"
#include "OuterCabling/outer_cabling_constants.hh"


// DESIGN
// Maximum number of modules per serial power chain
static const int inner_cabling_maxNumModulesPerPowerChain = 12;
// Maximum number of ELinks per GBT
static const int inner_cabling_maxNumELinksPerGBT = 7;
// Maximum number of modules per GBT
static const int inner_cabling_maxNumModulesPerGBT = 6;
// Maximum number of GBTs per bundle
static const int inner_cabling_maxNumGBTsPerBundle = 12;
// Maximum number of bundles per cable
static const int inner_cabling_maxNumBundlesPerCable = 6;

// GEOMETRY NAMES
// The cabling map is geometry-dependent, hence the subdetector names are placed here.
static const std::string inner_cabling_tbpx = "TBPX";
static const std::string inner_cabling_tfpx = "TFPX";
static const std::string inner_cabling_tepx = "TEPX";

// READOUT 1.28 GB/s E-LINKS
// Each value's index corresponds to the specific layer (TBPX) or disk (TFPX/TEPX) number.
static const std::unordered_map<std::string_view, std::unordered_map<int, int>> inner_cabling_numELinksPerModule = {
  {inner_cabling_tbpx, {{1, 6}, {2, 2}, {3, 2}, {4, 1}}},
  {inner_cabling_tfpx, {{1, 3}, {2, 3}, {3, 2}, {4, 2}}},
  {inner_cabling_tepx, {{1, 4}, {2, 3}, {3, 2}, {4, 2}, {5, 1}}}
};

// MAX NUMBER OF POWER CHAINS PER FIBER BUNDLE, in TBPX
// Each key corresponds to the specific layer number.
static const std::unordered_map<int, int> maxPowerChainsPerBundleTBPX = {
  {1, 1}, {2, 3}, {3, 3}, {4, 4}
};

// POWER CHAIN TYPE
// 4 ampere or 8 ampere.
enum PowerChainType { IUNDEFINED, I4A, I8A };


// ROUNDING
static const double inner_cabling_roundingTolerance = outer_cabling_roundingTolerance;

// CMSSW IDS
// Total number of DTCs in OT
static const int outer_cabling_totalNumDTCs = 216; // This could be obtained by:
// outer_cabling_numNonants * outer_cabling_maxNumDTCsPerNonantPerZEnd,
// but ugly to introduce dependency of IT cabling map on OT, because the maps are actually fully independent.


#endif  // INNER_CABLING_CONSTANTS_HH

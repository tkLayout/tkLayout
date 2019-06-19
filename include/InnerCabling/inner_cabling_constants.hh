#ifndef INNER_CABLING_CONSTANTS_HH
#define INNER_CABLING_CONSTANTS_HH

#include <math.h>
#include <string>
#include "global_funcs.hh"
#include "OuterCabling/outer_cabling_constants.hh"


// DESIGN
// Maximum number of modules per serial power chain
static const int inner_cabling_maxNumModulesPerPowerChain = 10;
// Maximum number of ELinks per GBT
static const int inner_cabling_maxNumELinksPerGBT = 7;
// Maximum number of GBTs per bundle
static const int inner_cabling_maxNumGBTsPerBundle = 12;
// Maximum number of bundles per cable
static const int inner_cabling_maxNumBundlesPerCable = 6;


// READOUT 1.28 GB/s E-LINKS
static const int inner_cabling_numELinksPerModuleBarrelLayer1 = 6;
static const int inner_cabling_numELinksPerModuleBarrelLayer2 = 2;
static const int inner_cabling_numELinksPerModuleBarrelLayer3 = 2;
static const int inner_cabling_numELinksPerModuleBarrelLayer4 = 1;

static const int inner_cabling_numELinksPerModuleForwardRing1 = 3;
static const int inner_cabling_numELinksPerModuleForwardRing2 = 2;
static const int inner_cabling_numELinksPerModuleForwardRing3 = 2;
static const int inner_cabling_numELinksPerModuleForwardRing4 = 1;

static const int inner_cabling_numELinksPerModuleEndcapRing1 = 3;
static const int inner_cabling_numELinksPerModuleEndcapRing2 = 2;
static const int inner_cabling_numELinksPerModuleEndcapRing3 = 1;
static const int inner_cabling_numELinksPerModuleEndcapRing4 = 1;
static const int inner_cabling_numELinksPerModuleEndcapRing5 = 1;


// MAX NUMBER OF POWER CHAINS PER FIBER BUNDLE, in TBPX
// This would be uselessly complicated to make this automatic.
static const int maxNumPowerChainsPerBundleBarrelLayer1 = 1;
static const int maxNumPowerChainsPerBundleBarrelLayer2 = 3;
static const int maxNumPowerChainsPerBundleBarrelLayer3 = 3;
static const int maxNumPowerChainsPerBundleBarrelLayer4 = 4;


// POWER CHAIN TYPE
// 4 ampere or 8 ampere.
enum PowerChainType { IUNDEFINED, I4A, I8A };


// ROUNDING
static const double inner_cabling_roundingTolerance = outer_cabling_roundingTolerance;


// GEOMETRY NAMES
// The cabling map is geometry-dependent, hence the subdetector names are placed here.
static const std::string inner_cabling_tbpx = "PXB";
static const std::string inner_cabling_tfpx = "FPIX_1";
static const std::string inner_cabling_tepx = "FPIX_2";


#endif  // INNER_CABLING_CONSTANTS_HH

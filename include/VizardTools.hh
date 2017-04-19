#ifndef VIZARDTOOLS_HH
#define VIZARDTOOLS_HH

#include <string>
#include <cmath>
#include <vector>
#include <Palette.hh>

/*
 * Assorted messages that may pop up
 */
static const std::string msg_uninitialised = "Vizard::buildVisualization(am, is) needs to be called first to build the visual geometry objects.";
static const std::string root_wrong = "Something went wrong creating output file. Existing geometry was not written to file.";
static const std::string graph_wrong = "File stream reported error state: neighbour graph not written to file.";
static const std::string exc_badalloc_graph = "Error: caught bad_alloc exception in Vizard::writeNeighbourGraph(). ";
static const std::string graph_nowrite = "Neighbour graph was not written to file.";

//clearStart="<tt>";
//clearEnd="</tt>";

// gStyle stuff
static const int style_grid = 3;
static const int materialNBins = 300;

// Colors for plot background and such
static const int color_plot_background       = kWhite;
static const int color_pad_background        = kGray;
static const int color_grid                  = kGreen-10;
static const int color_hard_grid             = kGray;
static const std::vector<std::string> color_names = {"Black","BrightBlue","Red","BrightGreen","Yellow","Pink","Aqua","Green","Blue"};


// Pads to plot the tracker ortho views
static const unsigned int padYZ = 1;
static const unsigned int padXY = 2;
static const unsigned int padProfile = 3;
static const unsigned int padEC = 4;

// Formatting parameters
static const int coordPrecision = 3;
static const int zOverlapPrecision = 5;
static const int anglePrecision = 1;
static const int areaPrecision = 1;
static const int occupancyPrecision = 1;
static const int rphiResolutionPrecision = 0;
static const int rphiResolutionRmsePrecision = 1;
static const int pitchPrecision = 0;
static const int stripLengthPrecision = 1;
static const int millionChannelPrecision = 2;
static const int totalPowerPrecision = 2;
static const int modulePowerPrecision = 0;
static const int powerPrecision = 1;
static const int powerPerUnitPrecision = 2;
static const int minimumBiasPrecision = 0;
static const int luminosityPrecision = 0;
static const int alphaParamPrecision = 4;
static const int weightPrecision = 0;

static const int tableResolutionPrecisionHigh = 5;
static const int tableResolutionPrecisionStd  = 3;


// Some strings for the html formatting
static const std::string subStart   = "<sub>";
static const std::string subEnd     = "</sub>";
static const std::string superStart = "<sup>";
static const std::string superEnd   = "</sup>";
static const std::string smallStart = "<small>";
static const std::string smallEnd   = "</small>";
static const std::string emphStart  = "<b>";
static const std::string emphEnd    = "</b>";
static const std::string muLetter   = "&mu;";
static const std::string sigmaLetter= "&sigma;";
static const std::string etaLetter  = "&eta;";
static const std::string phiLetter  = "&phi;";
static const std::string thetaLetter= "&theta;";
static const std::string deltaLetter= "&delta;";
static const std::string tauLetter  = "&tau;";

#endif

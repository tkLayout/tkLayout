/**
  * @file gui_constants.h
  * @author Nicoletta De Maio
  * @brief Messages, string and integer constants and custom data types are bundled here.
  */

#include <vector>
#include <utility>
#include <qstring.h>
#include <qpixmap.h>

/**
 * Messages that may pop up inside the status bar or on stdout
 */
static const QString msgSimulationExit = " application exited with status ";
static const QString msgCriticalErrorConfigFile = "Critical error handling config files. Aborting...";
static const QString msgErrValueLoad ="Error loading values: ";
static const QString msgErrSysCall = "Error during system() call.";
static const QString msgTemplateError = "Error choosing geometry template: ";
static const QString msgParamReset = "Parameters reset.";
static const QString msgErrParamTableAccess = "Unable to access parameter table: ";
static const QString msgErrParamCacheAccess = "Unable to access parameter backup table: ";
static const QString msgDefaultOverride ="Default settings saved." ;
static const QString msgSettingsRestored = "Settings file restored.";
static const QString msgErrSettingsRestore = "No file to restore from found.";
static const QString msgParamsRead = "Parameters read from file.";
static const QString msgParamsWritten = "Parameters written to file.";
static const QString msgErrReadFile = "Error opening read file.";
static const QString msgErrWriteFile = "Error opening write file.";
static const QString msgErrConfigFileParse = "Error parsing configuration file.";
static const QString msgResultsSaved = "Results saved.";
static const QString msgDirNotFound = "Error: directory does not exist.";
static const QString msgValidationError = "Validation error: check cost and power consumption input fields.";
static const QString msgSpinValidationError = "Error: values typed into the spinners must be integers.";
static const QString msgSpinValidationFuzzy = "Error: integer typed into spinner must be in range.";
static const QString msgErrValidationStrange = "A strange validation error happened...";

/**
 * String constants used throughout the GUI
 */
static const QString cAlphanumStartFilter = "[a-zA-Z0-9]*";
static const QString cGuiExtension = "/gui";
static const QString cResExtension = "/res";
static const QString cSettingsExtension = "/settings";
static const QString cSummaryExtension = "/summaries";
static const QString cRootDirExtension = "/store";
static const QString cRootFileExt = ".root";
static const QString cTmpDir = "/temp";
static const QString cCommand = "TrackerGeom2 ";
static const QString cRadioButtonBase = "radioButton";
static const QString cDescriptionName = "desc.txt";
static const QString cImageName = "layout.png";
static const QString cConfigFileName = "geometry.cfg";
static const QString cDefaultConfig = "defaultgeometry.cfg";
static const QString cDefaultSettings = "defaultsettings.cfg";
static const QString cSettingsBackup = "settings.cfg";
static const QString cSummaryIndex = "index.html";
static const QString cTempConfig = "tmpc.cfg";
static const QString cTempSettings = "tmps.cfg";
static const QString cDefaultTrackerName = "aTracker";

/**
 * Integer constants used throughout the GUI
 */
static const int cLayerChipModulus = 256;
static const int cRingChipModulus = 256;
static const int cDecimals = 15;
static const int cNonnegativeNumbers = 0;
static const int cPositiveNumbers = 1;
static const int cMaxChipsInSpinner = 8;

/**
 * Enumeration of the available module types: single sided, double sided and pt
 */
enum moduletype { none = -1, rphi = 0, stereo = 1, pt = 2 };

/**
 * Encapsulating struct for information about a pre-packaged geometry.
 * @param layoutDescription A short text description of the geometry; read from a file containing plaintext or simple HTML
 * @param layoutImage A diagram showing a cross-section of the detector geometry; read from a .png file
 * @param configFile The name of the original config file that comes with the geometry
 */
typedef struct geominfo {
    QString layoutDescription;
    QPixmap layoutImage;
    QString configFile;
};

/**
 * Encapsulating struct for parameteres in a pre-packaged geometry that may be modified by the user
 * @param trackerName A name for the optimisation experiment...
 * @param nlayers The number of layers, listed per barrel; displayed only, not customisable by the user
 * @param ndiscs The number of discs, listed per endcap; displayed only, not customisable by the user
 * @param nrings The number of rings in the selected geometry; displayed only, not customisable by the user
 * @param barrelnames A list of unique names to identify each barrel
 * @param endcapnames A list of unique names to identify each endcap
 * @param nchipslayer The number of chips across a module in each layer, listed in one vector per barrel
 * @param nchipsring The number of chips across a module in each ring, listed in one vector per endcap
 * @param nsegmentslayer The number of segments along a module in each layer, listed in one vector per barrel
 * @param nsegmentsring The number of segments along a module in each ring, listed in one vector per endcap
 * @param mtypeslayers The module type (rphi, stereo, pt or none) of each layer, listed in one vector per barrel
 * @param mtypesrings The module type (rphi, stereo, pt or none) of each ring, listed in one vector per endcap
 * @param costpersqcm What it costs, per square cm, to produce one single sided module; currency is left deliberately vague...
 * @param ptcostpersqcm What it costs, per square cm, to produce one pt module; see above for currency...
 * @param powerperchannel The current required by a single sided channel to function; in mA
 * @param ptpowerperchannel The current required by a pt channel to function; in mA
 */
typedef struct paramaggreg {
    QString trackerName;
    std::vector<int> nlayers;
    std::vector<int> ndiscs;
    std::vector<int> nrings;
    std::vector<QString> barrelnames;
    std::vector<QString> endcapnames;
    std::vector<std::vector<int> > nchipslayer;
    std::vector<std::vector<std::vector<int> > > nchipsring;
    std::vector<std::vector<int> > nsegmentslayer;
    std::vector<std::vector<std::vector<int> > > nsegmentsring;
    std::vector<std::vector<moduletype> > mtypeslayers;
    std::vector<std::vector<std::vector<moduletype> > > mtypesrings;
    double costpersqcm;
    double ptcostpersqcm;
    double powerperchannel;
    double ptpowerperchannel;
};

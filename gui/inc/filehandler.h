#include <stdexcept>
#include <iostream>
#include <qtextstream.h>
#include <qfile.h>
#include <qdir.h>
#include "gui_constants.h"

static const QString msgNoTracker = " No tracker definition found.";
static const QString msgNoBarrel = " No barrel definition found.";
static const QString msgNoEndcap = " No endcap definition found.";
static const QString msgErrTrackerBlock = " Error parsing tracker block.";
static const QString msgErrBarrelBlock = " Error parsing barrel block.";
static const QString msgErrEndcapBlock = " Error parsing endcap block.";
static const QString msgErrModuleDressing = " Error parsing additional module parameters.";
static const QString msgErrUnexpectedEndOfInput = " Unexpected end of input file.";
static const QString msgCleanupDirContents = "Cleanup: could not remove directory contents.";
static const QString msgCleanupNoDir = "Cleanup: requested result directory doesn't exist.";

static const QString trackerblock = "Tracker";
static const QString barrelblock = "Barrel";
static const QString endcapblock = "Endcap";
static const QString barreltypeblock = "BarrelType";
static const QString endcaptypeblock = "EndcapType";
static const QString outputblock = "Output";

static const QString scost = "stripCost";
static const QString pcost = "ptCost";
static const QString spow = "stripPower";
static const QString ppow = "ptPower";
static const QString layers = "nLayers";
static const QString chips = "nStripsAcross";
static const QString segs = "nSegments";
static const QString sides = "nSides";
static const QString type = "type";
static const QString outpath = "Path";

static const QString trphi = "rphi";
static const QString tstereo = "stereo";
static const QString tpt = "pt";

static const QString sob = "{";
static const QString eob = "}";
static const QString sep = "=";
static const QString soi = "[";
static const QString eoi = "]";
static const QString eol = ";";
static const QString pad = "  ";

class FileHandler {
public:
    void readConfigurationFromFile(QFile& readFile, paramaggreg& paramrow);
    void dressGeometry(QFile& readFile, paramaggreg& paramrow);
    void configureTracker(QFile& geometryFile, const paramaggreg& paramrow);
    void writeSettingsToFile(QFile& writeFile, const paramaggreg& paramrow, const QString& outputPath);
    void copyTextFile(QFile& inFile, QFile& outFile);
    void removeOutputDir(const QString outDir);
    void removeTmpConfigFile(const QString& fileName);
protected:
    bool cleanOutDirectory(QDir& workingDir);
private:
    void parseConfigFile(const QStringList& lineList, paramaggreg& paramrow);
    void parseSettingsFile(const QStringList& lineList, paramaggreg& paramrow);
    QString assembleSettingsFile(const paramaggreg& paramrow, const QString& relativeOutputPath);
    QStringList::const_iterator parseTrackerBlock(QStringList::const_iterator& iter,
						  const QStringList::const_iterator& end, paramaggreg& paramrow);
    QStringList::const_iterator parseBarrelBlock(QStringList::const_iterator& iter,
						 const QStringList::const_iterator& end, paramaggreg& paramrow);
    QStringList::const_iterator parseEndcapBlock(QStringList::const_iterator& iter, paramaggreg& paramrow);
    QStringList::const_iterator parseBarrelTypeBlock(QStringList::const_iterator& iter,
						     const QStringList::const_iterator& end, paramaggreg& paramrow);
    QStringList::const_iterator parseEndcapTypeBlock(QStringList::const_iterator& iter,
						     const QStringList::const_iterator& end, paramaggreg& paramrow);
    QStringList::iterator assembleTrackerBlock(QStringList::iterator& iter, const paramaggreg& paramrow);
    void appendBarrelTypeBlocks(QStringList& fileContents, const paramaggreg& paramrow);
    void appendEndcapTypeBlocks(QStringList& fileContents, const paramaggreg& paramrow);
    void appendOutputBlock(QStringList& fileContents, const QString& outputPath);
    int parseIndex(QString line);
    moduletype assignModuleType(QString toconvert);
};

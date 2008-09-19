#include "inc/filehandler.h"

/**
 * This function attempts to read a set of parameters for the parameter page from a file.
 * It will throw a runtime exception if it cannot open the file.
 * @param readFile A reference to the input file
 * @param paramrow The destination struct that absorbs the information found in the file
 */
void FileHandler::readConfigurationFromFile(QFile& readFile, paramaggreg& paramrow)
{
    if (readFile.open(IO_ReadOnly)) {
	QTextStream instream(&readFile);
	QString contents = instream.read();
	if (contents && !contents.isEmpty()) {
	    QStringList contentlist = QStringList::split("\n", contents);
	    QStringList::iterator iter = contentlist.begin();
	    while (iter != contentlist.end()) {
		(*iter) = (*iter).simplifyWhiteSpace();
		iter++;
	    }
	    parseConfigFile(contentlist, paramrow);
	}
	else {
	    readFile.close();
	    throw std::runtime_error(msgErrConfigFileParse);
	}
	readFile.close();
    }
    else throw std::runtime_error(msgErrReadFile);
}

void FileHandler::dressGeometry(QFile& readFile, paramaggreg& paramrow)
{
    if (readFile.open(IO_ReadOnly)) {
	QTextStream instream(&readFile);
	QString contents = instream.read();
	if (contents && !contents.isEmpty()) {
	    QStringList contentlist = QStringList::split("\n", contents);
	    QStringList::iterator iter = contentlist.begin();
	    while (iter != contentlist.end()) {
		(*iter) = (*iter).simplifyWhiteSpace();
		iter++;
	    }
	    parseSettingsFile(contentlist, paramrow);
	}
	else {
	    readFile.close();
	    throw std::runtime_error(msgErrConfigFileParse);
	}
	readFile.close();
    }
    else throw std::runtime_error(msgErrReadFile);
}

void FileHandler::configureTracker(QFile& geometryFile, const paramaggreg& paramrow)
{
    QString contents;
    if (geometryFile.open(IO_ReadOnly)) {
	QTextStream instream(&geometryFile);
	contents = instream.read();
	QStringList fileContents = QStringList::split("\n", contents);
	QStringList::iterator iter = fileContents.begin();
	while (iter != fileContents.end()) {
	    bool lineofinterest = FALSE;
	    if ((*iter).startsWith(trackerblock)) {
		iter = assembleTrackerBlock(iter, paramrow);
		lineofinterest = TRUE;
	    }
	    if (!lineofinterest) iter++;
	}
	contents = fileContents.join("\n");
	geometryFile.close();
    }
    else throw std::runtime_error(msgErrReadFile);
    if (geometryFile.open(IO_WriteOnly)) {
	QTextStream outstream(&geometryFile);
	outstream << contents;
	geometryFile.close();
    }
    else throw std::runtime_error(msgErrWriteFile);
}

/**
 * This function attempts to write the contents of the parameter page to a file (which is overwritten).
 * If an error occurs while opening or creating the output file, it throws a runtime exception.
 * @param writeFile A reference to the destination file
 * @param paramrow The struct containing the information about the selected geometry's current parameter values
 */
void FileHandler::writeSettingsToFile(QFile& writeFile, const paramaggreg& paramrow, const QString& outputPath)
{
    if (writeFile.open(IO_WriteOnly)) {
	QString contents;
	contents = assembleSettingsFile(paramrow, outputPath);
	QTextStream outstream(&writeFile);
	outstream << contents;
	writeFile.close();
    }
    else throw std::runtime_error(msgErrWriteFile);
}

/**
 * The process of copying the contents of one file to another (which is overwritten) is bundled here for convenience.
 * The function also displays a message in the status bar confirming the process.
 * If there are errors opening either file, a runtime exception is thrown.
 * @param inFile The data source
 * @param outFile The data destination
 * @param msg A message that will be displayed in the status bar
 */
void FileHandler::copyTextFile(QFile& inFile, QFile& outFile)
{
    if (inFile.open(IO_ReadOnly)) {
	if (outFile.open(IO_WriteOnly)) {
	    QTextStream fromStream(&inFile);
	    QTextStream toStream(&outFile);
	    toStream << fromStream.read();
	    outFile.close();
	    inFile.close();
	}
	else {
	    inFile.close();
	    throw std::runtime_error(msgErrWriteFile);
	}
    }
    else throw std::runtime_error(msgErrReadFile);
}

/**
 * Removal of the directory that contains the HTML summary happens in here.
 * If the directory is not empty, the files it contains are deleted first.
 */
void FileHandler::removeOutputDir(const QString outDir)
{
    QDir resultdir(outDir);
    if (resultdir.exists()) {
	if (cleanOutDirectory(resultdir)) resultdir.rmdir(resultdir.canonicalPath());
	else std::cout << msgCleanupDirContents << std::endl;
    }
    else std::cout << msgCleanupNoDir << std::endl;
}

void FileHandler::removeTmpConfigFile(const QString& fileName)
{
    QFile tmpFile(fileName);
    if (tmpFile.exists()) tmpFile.remove();
    tmpFile.setName(fileName + "~");
    if (tmpFile.exists()) tmpFile.remove();
}

/**
 * This function removes every file - but nothing else - from the given directory.
 * Note that it does not recurse through subfolders;
 * this is adequate for the summaries as they are now but may have to be adapted later.
 * @param workingDir The directory that is to be emptied
 * @return A boolean value indicating success or failure
 */
bool FileHandler::cleanOutDirectory(QDir& workingDir)
{
    bool success = TRUE;
    QStringList contents = workingDir.entryList("*", QDir::Files | QDir::Hidden);
    while (!contents.empty()) {
	if (!workingDir.remove(contents.first())) {
	    success = FALSE;
	    break;
	}
	contents.pop_front();
    }
    return success;
}

void FileHandler::parseConfigFile(const QStringList& lineList, paramaggreg& paramrow)
{
    bool tb = FALSE, bb = FALSE, eb = FALSE;
    QStringList::const_iterator iter = lineList.begin();
    while (iter != lineList.end()) {
	bool lineofinterest = FALSE;
	if ((*iter).startsWith(trackerblock)) {
	    lineofinterest = TRUE;
	    iter = parseTrackerBlock(iter, lineList.end(), paramrow);
	    tb = TRUE;
	}
	if ((*iter).startsWith(barrelblock) && !(*iter).startsWith(barreltypeblock)) {
	    lineofinterest = TRUE;
	    iter = parseBarrelBlock(iter, lineList.end(), paramrow);
	    bb = TRUE;
	}
	if ((*iter).startsWith(endcapblock) && !(*iter).startsWith(endcaptypeblock)) {
	    lineofinterest = TRUE;
	    iter = parseEndcapBlock(iter, paramrow);
	    eb = TRUE;
	}
	if (!lineofinterest) iter++;
    }
    if (!tb) throw std::runtime_error(msgErrConfigFileParse + msgNoTracker);
    if (!bb) throw std::runtime_error(msgErrConfigFileParse + msgNoBarrel);
    if (!eb) throw std::runtime_error(msgErrConfigFileParse + msgNoEndcap);
}

void FileHandler::parseSettingsFile(const QStringList& lineList, paramaggreg& paramrow)
{
    QStringList::const_iterator iter = lineList.begin();
    while (iter != lineList.end()) {
	bool lineofinterest = FALSE;
	if ((*iter).startsWith(barreltypeblock)) {
	    lineofinterest = TRUE;
	    iter = parseBarrelTypeBlock(iter, lineList.end(), paramrow);
	}
	if ((*iter).startsWith(endcaptypeblock)) {
	    lineofinterest = TRUE;
	    iter = parseEndcapTypeBlock(iter, lineList.end(), paramrow);
	}
	if (!lineofinterest) iter++;
    }
}

QString FileHandler::assembleSettingsFile(const paramaggreg& paramrow, const QString& relativeOutputPath)
{
    QStringList fileContents;
    appendBarrelTypeBlocks(fileContents, paramrow);
    appendEndcapTypeBlocks(fileContents, paramrow);
    appendOutputBlock(fileContents, relativeOutputPath);
    return fileContents.join("\n");
}

QStringList::const_iterator FileHandler::parseTrackerBlock(QStringList::const_iterator& iter,
							   const QStringList::const_iterator& end,paramaggreg& paramrow)
{
    bool sc = FALSE, pc = FALSE, sp = FALSE, pp = FALSE;
    QStringList::const_iterator itrator = iter;
    QString line = *itrator;
    QStringList words = QStringList::split(" ", line);
    if ((++words.begin()) == words.end()) throw std::runtime_error(msgErrConfigFileParse + msgErrTrackerBlock);
    if ((*(++words.begin())).startsWith(sob)) paramrow.trackerName = cDefaultTrackerName;
    else paramrow.trackerName = *(++words.begin());
    itrator++;
    while (*itrator != eob) {
	if (itrator == end) throw std::runtime_error(msgErrConfigFileParse + msgErrUnexpectedEndOfInput);
	line = *itrator;
	if (line.startsWith(scost)) {
	    int pos = line.find(sep) + sep.length();
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.costpersqcm = toconvert.stripWhiteSpace().toDouble();
	    sc = TRUE;
	}
	if (line.startsWith(pcost)) {
	    int pos = line.find(sep) + sep.length();
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.ptcostpersqcm = toconvert.stripWhiteSpace().toDouble();
	    pc = TRUE;
	}
	if (line.startsWith(spow)) {
	    int pos = line.find(sep) + sep.length();
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.powerperchannel = toconvert.stripWhiteSpace().toDouble();
	    sp = TRUE;
	}
	if (line.startsWith(ppow)) {
	    int pos = line.find(sep) + sep.length();
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.ptpowerperchannel = toconvert.stripWhiteSpace().toDouble();
	    pp = TRUE;
	}
	itrator++;
    }
    if (!sc) paramrow.costpersqcm = 0.0;
    if (!pc) paramrow.ptcostpersqcm = 0.0;
    if (!sp) paramrow.powerperchannel = 0.0;
    if (!pp) paramrow.ptpowerperchannel = 0.0;
    return itrator;
}

QStringList::const_iterator FileHandler::parseBarrelBlock(QStringList::const_iterator& iter,
							  const QStringList::const_iterator& end,paramaggreg& paramrow)
{
    bool nl = FALSE;
    QStringList::const_iterator itrator = iter;
    QString line = *itrator;
    QStringList words = QStringList::split(" ", line);
    if ((++words.begin()) == words.end() || (*(++words.begin())).startsWith(sob))
	throw std::runtime_error(msgErrConfigFileParse + msgErrBarrelBlock);
    paramrow.barrelnames.push_back(*(++words.begin()));
    itrator++;
    while (*itrator != eob) {
	if (itrator == end) throw std::runtime_error(msgErrConfigFileParse + msgErrUnexpectedEndOfInput);
	line = *itrator;
	if (line.startsWith(layers)) {
	    int pos = line.find(sep) + sep.length();
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.nlayers.push_back(toconvert.stripWhiteSpace().toInt());
	    std::vector<int> chips;
	    std::vector<int> segs;
	    std::vector<moduletype> types;
	    for (int i = 0; i < paramrow.nlayers.back(); i++) {
		chips.push_back(cLayerChipModulus);
		segs.push_back(1);
		types.push_back(none);
	    }
	    paramrow.nchipslayer.push_back(chips);
	    paramrow.nsegmentslayer.push_back(segs);
	    paramrow.mtypeslayers.push_back(types);
	    nl = TRUE;
	    break;
	}
	itrator++;
    }
    if (!nl) throw std::runtime_error(msgErrConfigFileParse + msgErrBarrelBlock);
    return itrator;
}

QStringList::const_iterator FileHandler::parseEndcapBlock(QStringList::const_iterator& iter, paramaggreg& paramrow)
{
    QStringList::const_iterator itrator = iter;
    QString line = *itrator;
    QStringList words = QStringList::split(" ", line);
    if ((++words.begin()) == words.end() || (*(++words.begin())).startsWith(sob))
	throw std::runtime_error(msgErrConfigFileParse + msgErrEndcapBlock);
    paramrow.endcapnames.push_back(*(++words.begin()));
    paramrow.nrings.push_back(0);
    itrator++;
    return itrator;
}

QStringList::const_iterator FileHandler::parseBarrelTypeBlock(QStringList::const_iterator& iter,
							      const QStringList::const_iterator& end,paramaggreg& paramrow)
{
    int pos;
    QStringList::const_iterator itrator = iter;
    QString line = *itrator;
    QStringList words = QStringList::split(" ", line);
    uint index = 0, idx;
    while (index < paramrow.barrelnames.size()) {
	if (paramrow.barrelnames.at(index) == *(++words.begin())) break;
	index++;
    }
    itrator++;
    while (*itrator != eob) {
	if (itrator == end) throw std::runtime_error(msgErrConfigFileParse + msgErrUnexpectedEndOfInput);
	line = *itrator;
	if (line.startsWith(chips)) {
	    idx = parseIndex(*itrator);
	    pos = line.find(sep) + 1;
	    if (idx >= paramrow.nchipslayer.size())  paramrow.nchipslayer.at(index).resize(idx + 1, cLayerChipModulus);
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.nchipslayer.at(index).at(idx) = toconvert.stripWhiteSpace().toInt();
	}
	if (line.startsWith(segs)) {
	    idx = parseIndex(*itrator);
	    pos = line.find(sep) + 1;
	    if (idx >= paramrow.nsegmentslayer.size())  paramrow.nsegmentslayer.at(index).resize(idx + 1, 1);
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.nsegmentslayer.at(index).at(idx) = toconvert.stripWhiteSpace().toInt();
	}
	if (line.startsWith(type)) {
	    idx = parseIndex(*itrator);
	    pos = line.find(sep) + 1;
	    if (idx >= paramrow.mtypeslayers.size())  paramrow.mtypeslayers.at(index).resize(idx + 1, rphi);
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    toconvert = toconvert.stripWhiteSpace();
	    paramrow.mtypeslayers.at(index).at(idx) = assignModuleType(toconvert);
	}
	itrator++;
    }
    return itrator;
}

QStringList::const_iterator FileHandler::parseEndcapTypeBlock(QStringList::const_iterator& iter,
							      const QStringList::const_iterator& end,paramaggreg& paramrow)
{
    int pos;
    QStringList::const_iterator itrator = iter;
    QString line = *itrator;
    QStringList words = QStringList::split(" ", line);
    uint index = 0, idx;
    while (index < paramrow.endcapnames.size()) {
	if (paramrow.endcapnames.at(index) == *(++words.begin())) break;
	index++;
    }
    if (index >= paramrow.nchipsring.size()) paramrow.nchipsring.resize(index + 1);
    if (index >= paramrow.nsegmentsring.size()) paramrow.nsegmentsring.resize(index + 1);
    if (index >= paramrow.mtypesrings.size()) paramrow.mtypesrings.resize(index + 1);
    itrator++;
    while (*itrator != eob) {
	if (itrator == end) throw std::runtime_error(msgErrConfigFileParse + msgErrUnexpectedEndOfInput);
	line = *itrator;
	if (line.startsWith(chips)) {
	    idx = parseIndex(*itrator);
	    pos = line.find(sep) + 1;
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    if (idx >= paramrow.nchipsring.at(index).size())  paramrow.nchipsring.at(index).resize(idx + 1, cRingChipModulus);
	    if (idx >= (uint)paramrow.nrings.at(index)) paramrow.nrings.at(index) = idx + 1;
	    paramrow.nchipsring.at(index).at(idx) = toconvert.stripWhiteSpace().toInt();
	}
	if (line.startsWith(segs)) {
	    idx = parseIndex(*itrator);
	    pos = line.find(sep) + 1;
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    if (idx >= paramrow.nsegmentsring.at(index).size())  paramrow.nsegmentsring.at(index).resize(idx + 1, 1);
	    if (idx > (uint)paramrow.nrings.at(index)) paramrow.nrings.at(index) = idx;
	    paramrow.nsegmentsring.at(index).at(idx) = toconvert.stripWhiteSpace().toInt();
	}
	if (line.startsWith(type)) {
	    idx = parseIndex(*itrator);
	    pos = line.find(sep) + 1;
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    toconvert = toconvert.stripWhiteSpace();
	    if (idx >= paramrow.mtypesrings.at(index).size())  paramrow.mtypesrings.at(index).resize(idx + 1, none);
	    if (idx > (uint)paramrow.nrings.at(index)) paramrow.nrings.at(index) = idx;
	    paramrow.mtypesrings.at(index).at(idx) = assignModuleType(toconvert);
	}
	itrator++;
    }
    return itrator;
}

QStringList::iterator FileHandler::assembleTrackerBlock(QStringList::iterator& iter, const paramaggreg& paramrow)
{
    int pos, len;
    QString value;
    QStringList::iterator itrator = iter;
    QString line = *itrator;
    pos = 0;
    len = line.length();
    value = trackerblock + " " + paramrow.trackerName + " " + sob;
    line.replace(pos, len, value);
    *itrator = line;
    itrator ++;
    while (*itrator != eob) {
	line = (*itrator).simplifyWhiteSpace();
	if (line.startsWith(scost)) {
	    pos = (*itrator).find(sep) + sep.length();
	    len = (*itrator).find(eol) - pos;
	    value = " " + QString::number(paramrow.costpersqcm);
	    *itrator = (*itrator).replace(pos, len, value);
	}
	if (line.startsWith(pcost)) {
	    pos = (*itrator).find(sep) + sep.length();
	    len = (*itrator).find(eol) - pos;
	    value = " " + QString::number(paramrow.ptcostpersqcm);
	    *itrator = (*itrator).replace(pos, len, value);
	}
	if (line.startsWith(spow)) {
	    pos = (*itrator).find(sep) + sep.length();
	    len = (*itrator).find(eol) - pos;
	    value = " " +QString::number(paramrow.powerperchannel);
	    *itrator = (*itrator).replace(pos, len, value);
	}
	if (line.startsWith(ppow)) {
	    pos = (*itrator).find(sep) + sep.length();
	    len = (*itrator).find(eol) - pos;
	    value = " " + QString::number(paramrow.ptpowerperchannel);
	    *itrator = (*itrator).replace(pos, len, value);
	}
	itrator++;
    }
    return itrator;
}

void FileHandler::appendBarrelTypeBlocks(QStringList& fileContents, const paramaggreg& paramrow)
{
    uint index;
    int idx;
    for (index = 0; index < paramrow.barrelnames.size(); index++) {
	QString line = barreltypeblock + " " + paramrow.barrelnames.at(index) + " " + sob;
	fileContents.append(line);
	for (idx = 0; idx < paramrow.nlayers.at(index); idx++) {
	    line = pad + chips + soi + QString::number(idx + 1) + eoi + " " + sep + " ";
	    line = line + QString::number(paramrow.nchipslayer.at(index).at(idx)) + eol;
	    fileContents.append(line);
	    line = "";
	    if (paramrow.mtypeslayers.at(index).at(idx) == rphi) line = "1" + eol;
            if (paramrow.mtypeslayers.at(index).at(idx) == stereo ||
            					paramrow.mtypeslayers.at(index).at(idx) == pt) line = "2" + eol;
            if (line.length() > 0) {
                line = pad + sides + soi + QString::number(idx + 1) + eoi + " " + sep + " " + line;
                fileContents.append(line);
            }
	    line = pad + segs + soi + QString::number(idx + 1) + eoi + " " + sep + " ";
	    line = line + QString::number(paramrow.nsegmentslayer.at(index).at(idx)) + eol;
	    fileContents.append(line);
            line = "";
	    if (paramrow.mtypeslayers.at(index).at(idx) == rphi) line = trphi + eol;
	    if (paramrow.mtypeslayers.at(index).at(idx) == stereo) line = tstereo + eol;
	    if (paramrow.mtypeslayers.at(index).at(idx) == pt) line = tpt + eol;
	    if (line.length() > 0) {
		line = pad + type + soi + QString::number(idx + 1) + eoi + " " + sep + " " + line;
		fileContents.append(line);
	    }
	}
	line = eob;
	fileContents.append(line);
    }
}

void FileHandler::appendEndcapTypeBlocks(QStringList& fileContents, const paramaggreg& paramrow)
{
    uint index;
    int idx;
    for (index = 0; index < paramrow.endcapnames.size(); index++) {
        if (paramrow.nrings.at(index) > 0) {
	QString line = endcaptypeblock + " " + paramrow.endcapnames.at(index) + " " + sob;
	fileContents.append(line);
	for (idx = 0; idx < paramrow.nrings.at(index); idx++) {
	    line = pad + chips + soi + QString::number(idx + 1) + eoi + " " + sep + " ";
	    line = line + QString::number(paramrow.nchipsring.at(index).at(idx)) + eol;
	    fileContents.append(line);
            line = "";
	    if (paramrow.mtypesrings.at(index).at(idx) == rphi) line = "1" + eol;
            if (paramrow.mtypesrings.at(index).at(idx) == stereo ||
            					paramrow.mtypesrings.at(index).at(idx) == pt) line = "2" + eol;
            if (line.length() > 0) {
                line = pad + sides + soi + QString::number(idx + 1) + eoi + " " + sep + " " + line;
                fileContents.append(line);
            }
	    line = pad + segs + soi + QString::number(idx + 1) + eoi + " " + sep + " ";
	    line = line + QString::number(paramrow.nsegmentsring.at(index).at(idx)) + eol;
	    fileContents.append(line);
            line = "";
	    if (paramrow.mtypesrings.at(index).at(idx) == rphi) line = trphi + eol;
	    if (paramrow.mtypesrings.at(index).at(idx) == stereo) line = tstereo + eol;
	    if (paramrow.mtypesrings.at(index).at(idx) == pt) line = tpt + eol;
	    if (line.length() != 0) {
		line = pad + type + soi + QString::number(idx + 1) + eoi + " " + sep + " " + line;
		fileContents.append(line);
	    }
	}
	line = eob;
	fileContents.append(line);
        }
    }
}

void FileHandler::appendOutputBlock(QStringList& fileContents, const QString& relativeOutputPath)
{
    if (relativeOutputPath != "") {
	QString line = outputblock + " " + sob;
	fileContents.append(line);
	line = pad + outpath + " " + sep + " " + relativeOutputPath + eol;
	fileContents.append(line);
	line = eob;
	fileContents.append(line);
    }
}

int FileHandler::parseIndex(QString line)
{
    int idx, start, stop;	    
    start = line.find(soi);
    stop = line.find(eoi);
    if (start > stop) throw std::runtime_error(msgErrConfigFileParse + msgErrModuleDressing);
    idx = line.mid(start + 1, stop - start - 1).toInt();
    return idx - 1;
}

moduletype FileHandler::assignModuleType(QString toconvert)
{
    if (toconvert == trphi) {
	return rphi;
    }
    if (toconvert == tstereo) {
	return stereo;
    }
    if (toconvert == tpt) {
	return pt;
    }
    return none;
}

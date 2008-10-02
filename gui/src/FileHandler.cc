/**
  * @file filehandler.cc
  * @author Nicoletta De Maio
  * @brief This is the implementation of the FileHandler helper class.
  */

#include "inc/filehandler.h"

/**
 * This function attempts to read a set of parameters for the parameter page from a geometry configuration file.
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

/**
 * This function attempts to read a set of module options for the parameter page from a settings file.
 * It will throw a runtime exception if it cannot open the file.
 * @param readFile A reference to the input file
 * @param paramrow The destination struct that absorbs the information found in the file
 */
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

/**
 * This function attempts to write the geometry-related contents of the parameter page to a file (which is modified).
 * If an error occurs while opening or creating the output file, it throws a runtime exception.
 * @param writeFile A reference to the destination file
 * @param paramrow The struct containing the information about the selected geometry
 */
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
 * This function attempts to write module options from the parameter page to a file (which is overwritten).
 * If an error occurs while opening or creating the output file, it throws a runtime exception.
 * @param writeFile A reference to the destination file
 * @param paramrow The struct containing the information about the selected geometry
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
 * If there are errors opening either file, a runtime exception is thrown.
 * @param inFile The data source
 * @param outFile The data destination
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
  * This is an overloaded function that does the same as <i>copyTextFile(QFile& inFile, QFile& outFile)</i>.
  * @param inName The source filename
  * @param outFile The destination file
  */
void FileHandler::copyTextFile(QString& inName, QFile& outFile)
{
    QFile inFile(inName);
    copyTextFile(inFile, outFile);
}

void FileHandler::copyDataFile(QFile& inFile, QFile& outFile)
{
    if (inFile.open(IO_ReadOnly)) {
	if (outFile.open(IO_WriteOnly)) {
	    QDataStream fromStream(&inFile);
	    QDataStream toStream(&outFile);
	    char byte;
	    while (!fromStream.atEnd()) {
		fromStream.readRawBytes(&byte, 1);
		toStream.writeRawBytes(&byte, 1);
	    }
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
 * Removal of the temporary directory that contains the HTML summary happens in here.
 * If the directory is not empty, the files it contains are deleted first.
 */
void FileHandler::removeOutputDir(const QString outDir)
{
    QDir resultdir(outDir);
    if (resultdir.exists()) {
	if (cleanOutDirectory(resultdir)) resultdir.rmdir(resultdir.canonicalPath());
	else std::cout << msgCleanupDirContents << std::endl;
    }
}

/**
 * If the file with the name <i>fileName</i> exists, it is removed in here.
 */
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

/**
 * This function looks for the block headers in the geometry configuration file and calls the appropriate parser for the block.
 * Unless a tracker and at least one barrel and one endcap definition are found, it throws an exception.
 * @param lineList The contents of the config file, broken up into individual lines
 * @param paramrow The destination data struct
 */
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
	    iter = parseEndcapBlock(iter, lineList.end(), paramrow);
	    eb = TRUE;
	}
	if (!lineofinterest) iter++;
    }
    if (!tb) throw std::runtime_error(msgErrConfigFileParse + msgNoTracker);
    if (!bb) throw std::runtime_error(msgErrConfigFileParse + msgNoBarrel);
    if (!eb) throw std::runtime_error(msgErrConfigFileParse + msgNoEndcap);
}

/**
 * This function looks for the block headers in a configuration file that dresses existng modules.
*  If it finds one, it calls the appropriate parser for the block.
 * @param lineList The contents of the config file, broken up into individual lines
 * @param paramrow The destination data struct
 */
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

/**
 * The information on how to dress the modules is read from the data struct and assembled here.
 * The output directory where the simulation results are going to be stored is appended last in a separate block.
 * @param paramrow The data source for the module configurations
 * @param relativeOutputPath The path to the designated output directory, relative to where the application is running from
 * @return The contents of the settings file in one long QString
 */
QString FileHandler::assembleSettingsFile(const paramaggreg& paramrow, const QString& relativeOutputPath)
{
    QStringList fileContents;
    appendBarrelTypeBlocks(fileContents, paramrow);
    appendEndcapTypeBlocks(fileContents, paramrow);
    appendOutputBlock(fileContents, relativeOutputPath);
    return fileContents.join("\n");
}

/**
 * This parser function reads parameters of interest from a tracker block into a <i>paramaggreg</i> data struct.
 * In case of missing data, it uses semi-sensible default values.
 * If the beginning of the tracker block is found to be malformed or the iteration encounters an unexpected end of file,
 * an exception is thrown.
 * @param iter An iterator to the current line in the configuration file
 * @param end The iterator pointing to the end of the list of lines
 * @param paramrow The destination data struct
 * @return An iterator pointing to the line starting with the closing brace of the tracker block
 */
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
    while (!(*itrator).startsWith(eob)) {
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

/**
 * This parser function reads parameters of interest from a barrel block into a <i>paramaggreg</i> data struct.
 * In case of missing data, it uses semi-sensible default values.
 * If the beginning of the block is found to be malformed or the iteration encounters an unexpected end of file,
 * an exception is thrown. The same thing happens if it cannot find any information about the number of layers.
 * @param iter An iterator to the current line in the configuration file
 * @param end The iterator pointing to the end of the list of lines
 * @param paramrow The destination data struct
 * @return An iterator pointing to the line starting with the closing brace of the barrel block
 */
QStringList::const_iterator FileHandler::parseBarrelBlock(QStringList::const_iterator& iter,
							  const QStringList::const_iterator& end, paramaggreg& paramrow)
{
    bool nl = FALSE;
    QStringList::const_iterator itrator = iter;
    QString line = *itrator;
    QStringList words = QStringList::split(" ", line);
    if ((++words.begin()) == words.end() || (*(++words.begin())).startsWith(sob))
	throw std::runtime_error(msgErrConfigFileParse + msgErrBarrelBlock);
    paramrow.barrelnames.push_back(*(++words.begin()));
    itrator++;
    while (!(*itrator).startsWith(eob)) {
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

/**
 * This parser function reads the endcap name from the first line of an edcap block into a <i>paramaggreg</i> data struct.
 * It also initialises the corresponding vector entry for the number of rings, but none of the other vectors (the GUI handles those).
 * If the beginning of the block is found to be malformed or the iteration encounters an unexpected end of file,
 * an exception is thrown.
 * @param iter An iterator to the current line in the configuration file
 * @param paramrow The destination data struct
 * @return An iterator pointing to the line after the one starting with the endcap block definition
 */
QStringList::const_iterator FileHandler::parseEndcapBlock(QStringList::const_iterator& iter,
							  const QStringList::const_iterator& end, paramaggreg& paramrow)
{
    bool nl = FALSE;
    QStringList::const_iterator itrator = iter;
    QString line = *itrator;
    QStringList words = QStringList::split(" ", line);
    if ((++words.begin()) == words.end() || (*(++words.begin())).startsWith(sob))
	throw std::runtime_error(msgErrConfigFileParse + msgErrEndcapBlock);
    paramrow.endcapnames.push_back(*(++words.begin()));
    paramrow.nrings.push_back(0);
    itrator++;
    while (!(*itrator).startsWith(eob)) {
	if (itrator == end) throw std::runtime_error(msgErrConfigFileParse + msgErrUnexpectedEndOfInput);
	line = *itrator;
	if (line.startsWith(discs)) {
	    int pos = line.find(sep) + sep.length();
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.ndiscs.push_back(toconvert.stripWhiteSpace().toInt());
	    std::vector<std::vector<int> > ringschips;
	    std::vector<std::vector<int> > ringssegs;
	    std::vector<std::vector<moduletype> > ringstypes;
	    for (int i = 0; i < paramrow.ndiscs.back(); i++) {
		paramrow.nchipsring.push_back(ringschips);
		paramrow.nsegmentsring.push_back(ringssegs);
		paramrow.mtypesrings.push_back(ringstypes);
	    }
	    nl = TRUE;
	    break;
	}
	itrator++;
    }
    if (!nl) throw std::runtime_error(msgErrConfigFileParse + msgErrBarrelBlock);
    return itrator;
}

/**
 * This parser function reads parameters of interest from a block containing module options per layer.
 * The extracted information is added to a <i>paramaggreg</i> data struct.
 * If the iteration encounters an unexpected end of file, an exception is thrown.
 * @param iter An iterator to the current line in the configuration file
 * @param end The iterator pointing to the end of the list of lines
 * @param paramrow The destination data struct
 * @return An iterator pointing to the line starting with the closing brace of the block
 */
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
    while (!(*itrator).startsWith(eob)) {
	if (itrator == end) throw std::runtime_error(msgErrConfigFileParse + msgErrUnexpectedEndOfInput);
	line = *itrator;
	if (line.startsWith(chips)) {
	    idx = parse1DIndex(*itrator);
	    pos = line.find(sep) + 1;
	    if (idx >= paramrow.nchipslayer.size())  paramrow.nchipslayer.at(index).resize(idx + 1, cLayerChipModulus);
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.nchipslayer.at(index).at(idx) = toconvert.stripWhiteSpace().toInt();
	}
	if (line.startsWith(segs)) {
	    idx = parse1DIndex(*itrator);
	    pos = line.find(sep) + 1;
	    if (idx >= paramrow.nsegmentslayer.size())  paramrow.nsegmentslayer.at(index).resize(idx + 1, 1);
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    paramrow.nsegmentslayer.at(index).at(idx) = toconvert.stripWhiteSpace().toInt();
	}
	if (line.startsWith(type)) {
	    idx = parse1DIndex(*itrator);
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

/**
 * This parser function reads parameters of interest from a block containing module options per ring.
 * The extracted information is added to a <i>paramaggreg</i> data struct, resizing the vectors as necessary.
 * If the iteration encounters an unexpected end of file, an exception is thrown.
 * @param iter An iterator to the current line in the configuration file
 * @param end The iterator pointing to the end of the list of lines
 * @param paramrow The destination data struct
 * @return An iterator pointing to the line starting with the closing brace of the block
 */
QStringList::const_iterator FileHandler::parseEndcapTypeBlock(QStringList::const_iterator& iter,
							      const QStringList::const_iterator& end,paramaggreg& paramrow)
{
    int pos;
    QStringList::const_iterator itrator = iter;
    QString line = *itrator;
    QStringList words = QStringList::split(" ", line);
    uint index = 0;
    std::pair<int,int> idx;
    while (index < paramrow.endcapnames.size()) {
	if (paramrow.endcapnames.at(index) == *(++words.begin())) break;
	index++;
    }
    for (uint i = 0; i < paramrow.ndiscs.size(); i++) {
	if (index >= paramrow.nchipsring.at(i).size()) paramrow.nchipsring.at(i).resize(index + 1);
	if (index >= paramrow.nsegmentsring.at(i).size()) paramrow.nsegmentsring.at(i).resize(index + 1);
	if (index >= paramrow.mtypesrings.at(i).size()) paramrow.mtypesrings.at(i).resize(index + 1);
    }
    itrator++;
    while (!(*itrator).startsWith(eob)) {
	if (itrator == end) throw std::runtime_error(msgErrConfigFileParse + msgErrUnexpectedEndOfInput);
	line = *itrator;
	if (line.startsWith(chips)) {
	    idx = parse2DIndex(*itrator);
	    pos = line.find(sep) + 1;
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    if ((uint)idx.second >= paramrow.nchipsring.at(index).size())
		throw std::runtime_error(msgErrConfigFileParse + msgErrEndcapTypeBlock);
	    if ((uint)idx.first >= paramrow.nchipsring.at(index).at(idx.second).size())
		paramrow.nchipsring.at(index).at(idx.second).resize(idx.first + 1, cRingChipModulus);
	    if (idx.first >= paramrow.nrings.at(index)) paramrow.nrings.at(index) = idx.first + 1;
	    paramrow.nchipsring.at(index).at(idx.second).at(idx.first) = toconvert.stripWhiteSpace().toInt();
	}
	if (line.startsWith(segs)) {
	    idx = parse2DIndex(*itrator);
	    pos = line.find(sep) + 1;
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    if ((uint)idx.second >= paramrow.nsegmentsring.at(index).size())
		throw std::runtime_error(msgErrConfigFileParse + msgErrEndcapTypeBlock);
	    if ((uint)idx.first >= paramrow.nsegmentsring.at(index).at(idx.second).size())
		paramrow.nsegmentsring.at(index).at(idx.second).resize(idx.first + 1, 1);
	    if (idx.first >= paramrow.nrings.at(index)) paramrow.nrings.at(index) = idx.first + 1;
	    paramrow.nsegmentsring.at(index).at(idx.second).at(idx.first) = toconvert.stripWhiteSpace().toInt();
	}
	if (line.startsWith(type)) {
	    idx = parse2DIndex(*itrator);
	    pos = line.find(sep) + 1;
	    QString toconvert = line.remove(0, pos);
	    toconvert.truncate(toconvert.find(eol));
	    toconvert = toconvert.stripWhiteSpace();
	    if ((uint)idx.second >= paramrow.mtypesrings.at(index).size())
		throw std::runtime_error(msgErrConfigFileParse + msgErrEndcapTypeBlock);
	    if ((uint)idx.first >= paramrow.mtypesrings.at(index).at(idx.second).size())
		paramrow.mtypesrings.at(index).at(idx.second).resize(idx.first + 1, none);
	    if (idx.first >= paramrow.nrings.at(index)) paramrow.nrings.at(index) = idx.first + 1;
	    paramrow.mtypesrings.at(index).at(idx.second).at(idx.first) = assignModuleType(toconvert);
	}
	itrator++;
    }
    for (int i = 0; i < paramrow.ndiscs.at(index); i++) {
	if ((uint)paramrow.nrings.at(index) >= paramrow.nchipsring.at(i).at(index).size())
	    paramrow.nchipsring.at(index).at(i).resize(paramrow.nrings.at(index), cRingChipModulus);
	if ((uint)paramrow.nrings.at(index) >= paramrow.nsegmentsring.at(i).at(index).size())
	    paramrow.nsegmentsring.at(index).at(i).resize(paramrow.nrings.at(index), 1);
	if ((uint)paramrow.nrings.at(index) >= paramrow.mtypesrings.at(i).at(index).size())
	    paramrow.mtypesrings.at(index).at(i).resize(paramrow.nrings.at(index), none);
    }
    return itrator;
}

/**
 * This function replaces name, cost parameters and power parameters in an existing tracker block with
 * the information given by the user in a <i>paramaggreg</i> data struct.
 * @param iter An iterator to the current line in the configuration file
 * @param paramrow The source data struct
 * @return An iterator pointing to the line starting with the closing brace of the block
 */
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
    while (!(*itrator).startsWith(eob)) {
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

/**
 * This function appends a barrel type block to a module settings file.
 * @param fileContents The list of config file lines that the block is to be appended to
 * @param paramrow The source data struct
 */
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

/**
 * This function appends an endcap type block to a module settings file.
 * @param fileContents The list of config file lines that the block is to be appended to
 * @param paramrow The source data struct
 */
void FileHandler::appendEndcapTypeBlocks(QStringList& fileContents, const paramaggreg& paramrow)
{
    uint index;
    int idx, id;
    for (index = 0; index < paramrow.endcapnames.size(); index++) {
	if (paramrow.nrings.at(index) > 0) {
	    QString line = endcaptypeblock + " " + paramrow.endcapnames.at(index) + " " + sob;
	    fileContents.append(line);
	    for (idx = 0; idx < paramrow.nrings.at(index); idx++) {
		for (id = 0; id < paramrow.ndiscs.at(index); id++) {
		    line = pad + chips + soi + QString::number(idx + 1) + "," + QString::number(id + 1) + eoi + " ";
		    line = line + sep + " " + QString::number(paramrow.nchipsring.at(index).at(id).at(idx)) + eol;
		    fileContents.append(line);
		    line = "";
		    if (paramrow.mtypesrings.at(index).at(id).at(idx) == rphi) line = "1" + eol;
		    if (paramrow.mtypesrings.at(index).at(id).at(idx) == stereo ||
			paramrow.mtypesrings.at(index).at(id).at(idx) == pt) line = "2" + eol;
		    if (line.length() > 0) {
			line = pad + sides + soi + QString::number(idx + 1) + "," + QString::number(id + 1);
			line = line + eoi + " " + sep + " " + line;
			fileContents.append(line);
		    }
		    line = pad + segs + soi + QString::number(idx + 1) + "," + QString::number(id + 1) + eoi + " ";
		    line = line + sep + " " + QString::number(paramrow.nsegmentsring.at(index).at(id).at(idx)) + eol;
		    fileContents.append(line);
		    line = "";
		    if (paramrow.mtypesrings.at(index).at(id).at(idx) == rphi) line = trphi + eol;
		    if (paramrow.mtypesrings.at(index).at(id).at(idx) == stereo) line = tstereo + eol;
		    if (paramrow.mtypesrings.at(index).at(id).at(idx) == pt) line = tpt + eol;
		    if (line.length() != 0) {
			line = pad + type + soi + QString::number(idx + 1) + "," + QString::number(id + 1);
			line = line + eoi + " " + sep + " " + line;
			fileContents.append(line);
		    }
		}
	    }
	    line = eob;
	    fileContents.append(line);
	}
    }
}

/**
 * This function appends an output block to a module settings file.
 * @param fileContents The list of config file lines that the block is to be appended to
 * @param paramrow The destination path, relative to where the application is running from, in QString form
 */
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

/**
 * Parsing of a layer or ring index from an indexed line in a module dressing block is bundled in here.
 * @param line The line within the config file containing the parameter of interest
 * @return The index as it will be used by the internal data structures (i.e. the converted number - 1)
 */
int FileHandler::parse1DIndex(QString line)
{
    int idx, start, stop;	    
    start = line.find(soi);
    stop = line.find(eoi);
    if (start > stop) throw std::runtime_error(msgErrConfigFileParse + msgErrModuleDressing);
    idx = line.mid(start + 1, stop - start - 1).toInt();
    return idx - 1;
}

std::pair<int,int> FileHandler::parse2DIndex(QString line)
{
    int start, stop;
    std::pair<int,int> rd;
    start = line.find(soi);
    if (start < 0) throw std::runtime_error(msgErrConfigFileParse + msgErrModuleDressing);
    stop = line.find(",", start);
    if (stop < 0) throw std::runtime_error(msgErrConfigFileParse + msgErrModuleDressing);
    rd.first = line.mid(start + 1, stop - start - 1).toInt() - 1;
    start = stop;
    stop = line.find(eoi, start);
    if (stop < 0) throw std::runtime_error(msgErrConfigFileParse + msgErrModuleDressing);
    rd.second = line.mid(start + 1, stop - start - 1).toInt() - 1;
    return rd;
}

/**
 * This little convenience function matches a string from an appropriate line in the config file to its module type.
 * @param toconvert The string value from the config file
 * @return One of the values listed in <i>enum moduletype</i>
 */
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

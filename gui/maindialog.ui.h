/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you want to add, delete, or rename functions or slots, use
** Qt Designer to update this file, preserving your code.
**
** You should not define a constructor or destructor in this file.
** Instead, write your code in functions called init() and destroy().
** These will automatically be called by the form's constructor and
** destructor.
*****************************************************************************/

/**
 * @file maindialog.ui.h
 * @author Nicoletta De Maio
 * @brief This class creates the tkgeometry GUI and handles events related to it. It starts the simulation in the background and displays the results. It also gives an option to save those parameters that can be modified by the user in a file. 
 * Non-widget class variables:
 * - Paths (string values encapsulated in QString instances)
 *   - @param basePath The path to the tkgeometry installation
 *   - @param guiExtension The location of the GUI executable
 *   - @param resExtension The base directory for geometry resources
 *   - @param settingsExtension The location of the parameter setting files within a geometry resource
 *   - @param summaryExtension The base directory for the HTML summaries
 *   - @param outDirExtension The output directory for the chosen geometry and parameters; constructed at runtime
 *   - @param cmdLineStub A buffer for the argument that is passed to <i>signal()</i>; assebled from various parameters
 * - Data collections (vectors of <i>geominfo</i> or <i>paramaggreg</i>)
 *   - @param geometryTable A vector of <i>geominfo</i> for immutable information about the available geometries
 *   - @param parameterTable A vector of <i>paramaggreg</i> for parameters that can be changed by the user
 *   - @param widgetCache A vector of <i>paramaggreg</i> that backs up defaults for parameters exposed to the user
 * - Pointer to popup
 *   - @param *settingsPopup A popup menu containing the options to save and retrieve parameter settings from file
 */

/**
 * Messages that may pop up inside the status bar or on stdout
 */
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
static const QString msgErrInputFileParse = "Error parsing input file.";
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
static const QString cCommand = "/TrackerGeom ";
static const QString cRadioButtonBase = "radioButton";
static const QString cDescriptionName = "desc.txt";
static const QString cImageName = "layout.png";
static const QString cCmdLineName = "cmdline";
static const QString cDefaultSettings = "defaultsettings";
static const QString cSummaryIndex = "index.html";
static const QString cDefaultTrackerName = "aTracker";

/**
 * Integer constants used throughout the GUI
 */
static const int cDefaultChipsPerSide = 128;
static const int cDefaultSegsPerRing = 128;
static const int cChipModulus = 2;
static const int cDecimals = 15;
static const int cNonnegativeNumbers = 0;
static const int cPositiveNumbers = 1;

/**
 * Enumeration of the available module types: single sided, double sided and pt
 */
enum moduletype { ss, ds, pt };

/**
 * Encapsulating struct for information about a pre-packaged geometry.
 * @param layoutName A string value giving the name that the geometry should be listed as
 * @param layoutDescription A short text description of the geometry; read from a file containing plaintext or simple HTML
 * @param layoutImage A diagram showing a cross-section of the detector geometry; read from a .png file
 * @param cmdLineParams Those parameters not exposed to the user, aggregated in a single string value
 */
typedef struct geominfo {
    QString layoutName;
    QString layoutDescription;
    QPixmap layoutImage;
    QString cmdLineParams;
};

/**
 * Encapsulating struct for parameteres in a pre-packaged geometry that may be modified by the user
 * @param trackerName A name for the optimisation experiment...
 * @param nlayers The number of layers in the selected geometry; displayed only, not customisable by the user
 * @param nrings The number of rings in the selected geometry; displayed only, not customisable by the user
 * @param nchipsperside The number of chips on a module in each layer (per side, not in total), listed in a vector
 * @param nsegmentsperring The number of segments on a module in each ring, listed in a vector
 * @param mtypeslayers The module type (single sided, double sided or pt) of each layer, listed in a vector
 * @param mtypesrings The module type (single sided, double sided or pt) of each ring, listed in a vector
 * @param costpersqcm What it costs, per square cm, to produce one single sided module; currency is left deliberately vague...
 * @param ptcostpersqcm What it cost, per square cm, to produce one pt module; see above for currency...
 * @param powerperchannel The current required by a single sided channel to function; in mA
 * @param ptpowerperchannel The current required by a pt channel to function; in mA
 */
typedef struct paramaggreg {
    QString trackerName;
    int nlayers;
    int nrings;
    std::vector<int> nchipsperside;
    std::vector<int> nsegmentsperring;
    std::vector<moduletype> mtypeslayers;
    std::vector<moduletype> mtypesrings;
    double costpersqcm;
    double ptcostpersqcm;
    double powerperchannel;
    double ptpowerperchannel;
};

/**
 * This is the initialisation function that serves as a programmer-defined addition to the constructor.
 * (The constructor itself is private and auto-generated.)
 * It initialises some of the widgets programmatically, reads the information about the available geometries from the
 * resource directory and caches it in the global table vectors <i>geometryTable , parameterTable</i> and <i>widgetCache</i>.
 * If a resource directory is found to contain nonsensical parameter files it is skipped.
 */
void MainDialog::init()
{
    /// Add the menu items to the settings popup menu and connect the widget to the main form
    settingsPopup = new QPopupMenu(settingsButton, "Settings");
    settingsPopup->insertItem("Default Values", 0, 0);
    settingsPopup->insertItem("Load Settings...", 1, 1);
    settingsPopup->insertItem("Save Settings...", 2, 2);
    settingsPopup->insertItem("Save As Default", 3, 3);
    settingsPopup->insertItem("Restore Defaults", 4, 4);
    connect(settingsPopup, SIGNAL(activated(int)), this, SLOT(settingsChanged(int)));
    settingsButton->setPopup(settingsPopup);

    /// Add validators to the spinners and those text input fields that only accept numbers
    QIntValidator *intval = new QIntValidator(chipsPerSideSpinner);
    intval->setBottom(cPositiveNumbers);
    chipsPerSideSpinner->setValidator(intval);
    
    intval = new QIntValidator(segmentsPerRingSpinner);
    intval->setBottom(cPositiveNumbers);
    segmentsPerRingSpinner->setValidator(intval);
    
    QDoubleValidator *doubleval = new QDoubleValidator(costPerSqCmEdit);
    doubleval->setBottom(cNonnegativeNumbers);
    doubleval->setDecimals(cDecimals);
    costPerSqCmEdit->setValidator(doubleval);
    
    doubleval = new QDoubleValidator(costPtPerSqCmEdit);
    doubleval->setBottom(cNonnegativeNumbers);
    doubleval->setDecimals(cDecimals);
    costPtPerSqCmEdit->setValidator(doubleval);
    
    doubleval = new QDoubleValidator(powerEdit);
    doubleval->setBottom(cNonnegativeNumbers);
    doubleval->setDecimals(cDecimals);
    powerEdit->setValidator(doubleval);
    
    doubleval = new QDoubleValidator(ptPowerEdit);
    doubleval->setBottom(cNonnegativeNumbers);
    doubleval->setDecimals(cDecimals);
    ptPowerEdit->setValidator(doubleval);

    /// Set those paths that will remain constant at runtime
    guiExtension = cGuiExtension;
    resExtension = cResExtension;
    settingsExtension = cSettingsExtension;
    summaryExtension = cSummaryExtension;

    /// Create directory access to read the geometry resources
    QDir workingDir;
    workingDir.setSorting(QDir::Name | QDir::IgnoreCase);
    workingDir.setFilter(QDir::Dirs);
    workingDir.setNameFilter(cAlphanumStartFilter);

    /// Extract the base path for the application
    basePath = workingDir.canonicalPath();
    if (basePath.endsWith(guiExtension)) basePath.truncate(basePath.length() - guiExtension.length());

    /// Go to the resources base directory and prepare the radio button group
    workingDir.setPath(basePath + resExtension);
    geometryPicker->setColumnLayout(1, Qt::Vertical);

    /// Work your way through the contents of the resource directory (set to only include subfolders); ignore nonsensical entries
    QStringList::const_iterator iter;
    for (iter = workingDir.entryList().constBegin(); iter != workingDir.entryList().constEnd(); iter++) {
	geominfo geomrow;
	geomrow.layoutName = *iter;

	/// Read the text description of the geometry from file
	QFile workingFile(workingDir.canonicalPath() + "/" + *iter + "/" + cDescriptionName);
	if (workingFile.exists() && workingFile.open(IO_ReadOnly)) {
	    QTextStream instream(&workingFile);
	    geomrow.layoutDescription = instream.read();
	    workingFile.close();
	}
	else geomrow.layoutDescription = "";

	/// Read the geometry diagram from file
	workingFile.setName(workingDir.canonicalPath() + "/" + *iter + "/" + cImageName);
	if (!geomrow.layoutImage.load(workingFile.name())) geomrow.layoutImage.resize(1,1);

	/// Read the immutable command line parameters from file, then add the collected information to <i>geometryTable</i>
	workingFile.setName(workingDir.canonicalPath() + "/" + *iter + "/" + cCmdLineName);
	if (workingFile.exists() && workingFile.open(IO_ReadOnly)) {
	    QTextStream instream(&workingFile);
	    geomrow.cmdLineParams = instream.readLine();
	    workingFile.close();
	}
	else continue;

	/// Read the cusomisable parameters from their settings file and add them to <i>parameterTable</i> and <i>widgetCache</i>
	paramaggreg paramrow;
	workingFile.setName(workingDir.canonicalPath() + "/" + *iter + settingsExtension + "/" + cDefaultSettings);
	if (workingFile.exists()) {
	    try {
		readParametersFromFile(workingFile, paramrow);
	    }
	    catch (std::runtime_error re) {
		std::cout << re.what() << std::endl;
	    }
	}
	else continue;
	
	/// Add the assembled info structs to their vector containers
	geometryTable.push_back(geomrow);
	parameterTable.push_back(paramrow);
	widgetCache.push_back(paramrow);

	/// Add a radio button for the geometry to the widget so that the user can select it later
	QString name = cRadioButtonBase + "_";
	name += *iter;
	QRadioButton *geometrySelection = new QRadioButton(*iter, geometryPicker, name.ascii());
	geometryPicker->insert(geometrySelection, workingDir.entryList().findIndex(*iter));
    }
}

/**
 * Since the destructor is as private and auto-generated as the constructor, cleanup of programmer-defined stuff happens in here.
 * This consists of deleting the settings popup menu (which was created with <i>new</i>) and of removing the output folder.
 */
void MainDialog::destroy()
{
    if (settingsPopup) delete settingsPopup;
    removeOutputDir();
}

/**
 * Removal of the directory that contains the HTML summary happens in here.
 * If the directory is not empty, the files it contains are deleted first.
 */
void MainDialog::removeOutputDir()
{
    QDir resultdir(basePath + summaryExtension + outDirExtension);
    if (resultdir.exists()) {
	if (cleanOutDirectory(resultdir)) resultdir.rmdir(basePath + summaryExtension + outDirExtension);
	else std::cout << "Cleanup: could not remove directory contents." << std::endl;
    }
    else std::cout << "Cleanup: requested result directory doesn't exist." << std::endl;
}

/**
 * This function removes every file - but nothing else - from the given directory.
 * Note that it does not recurse through subfolders;
 * this is adequate for the summaries as they are now but may have to be adapted later.
 * @param workingDir The directory that is to be emptied
 * @return A boolean value indicating success or failure
 */
bool MainDialog::cleanOutDirectory(QDir& workingDir)
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
 * This is the navigation function that moves focus from the parameter page back to the geometry selection page.
 * It is connected to the <i>Back</i> button on the parameter page.
 * It resets the content of the combo boxes to none so that they may be filled again for the next selected geometry.
 */
void MainDialog::settingsToGeometry()
{
    layerSelection->clear();
    ringSelection->clear();
    backOnePage();
}

/**
 * This is the navigation function that moves focus from the results page back to the parameter page.
 * It is connected to the <i>Back</i> button on the results page.
 * It eliminates the output files from the previous simulation since they are no longer needed.
 * (The simulation will have to be re-run in order to get back to the results page.)
 */
void MainDialog::resultsToSettings()
{
    removeOutputDir();
    backOnePage();
}

/**
 * The command that tells the widget stack to raise the previous page is encapsulated in this function for convenience.
 */
void MainDialog::backOnePage()
{
    mainWidgetStack->raiseWidget(mainWidgetStack->id(mainWidgetStack->visibleWidget()) - 1);
}

/**
 * This is the navigation function that moves focus from the geometry selection page to the parameter page.
 * It is connected to the <i>Next</i> button on the geometry selection page.
 * Before raising the parameter page on the widget stack, it initialises it based on the selected geometry.
 */
void MainDialog::nextPage()
{
    statusBar->clear();
    try {
	for (int i = 0; i < parameterTable.at(geometryPicker->selectedId()).nlayers; i++ ) {
	    layerSelection->insertItem(QString::number(i + 1), i);
	}
	for (int i = 0; i < parameterTable.at(geometryPicker->selectedId()).nrings; i++) {
	    ringSelection->insertItem(QString::number(i + 1), i);
	}
	valuesToWidgets(parameterTable.at(geometryPicker->selectedId()));
    }
    catch (std::out_of_range oor) {
	QString statusText(msgErrValueLoad);
	statusText += oor.what();
	statusBar->setText(statusText);
    }
    mainWidgetStack->raiseWidget(mainWidgetStack->id(mainWidgetStack->visibleWidget()) + 1);
}

/**
 * This is the navigation function that moves focus from the results page back to the geometry selection page.
 * It is connected to the <i>New</i> button on the results page.
 * It resets the contents of the combo boxes to none and removes the summary files in preparation of a new simulation run.
 * Then it raises the first widget on the stack.
 */
void MainDialog::backToStart()
{
    removeOutputDir();
    layerSelection->clear();
    ringSelection->clear();
    clearParameters();
    mainWidgetStack->raiseWidget(0);
}


// TODO: finish this
void MainDialog::go()
{
    if (validateInput()) {
	QString command = basePath + cCommand + " ";
	// collect parameters from parameterTable, appending to cmdLineStub
	// set outDirExtension
	outDirExtension = "/testdir"; // temporary hack...
	// if save folder exists, confirm overwrite
	// complete command line string
	command += cmdLineStub;
	try  {
	    int simulstatus;	
	    simulstatus = simulate(command);
	    // write evaluation of simulstatus to stdout?
	    // update summaryPage in background: read index file in basePath + summaryExtension + outDirExtension + "/index.html"
	    summaryTextEdit->mimeSourceFactory()->setFilePath(basePath + summaryExtension + outDirExtension);
	    QFile workingFile(basePath + summaryExtension + outDirExtension + "/" + cSummaryIndex);
	    if (workingFile.exists() && workingFile.open(IO_ReadOnly)) {
		QTextStream instream(&workingFile);
		summaryTextEdit->setText(instream.read());
	    }
	    else summaryTextEdit->clear();
	    mainWidgetStack->raiseWidget(mainWidgetStack->id(mainWidgetStack->visibleWidget()) + 1);
	}
	catch (std::runtime_error re) {
	    // evaluate exception "something went wrong during simulation"
	    // statusBar->setText(/*error message*/);
	}
    }
}
// End of TODO

/**
 * Control is transferred to the <i>TrackerGeom</i> background application in here.
 * Once that exits, the function throws an exception if the attempt to run <i>TrackerGeom</i> was unsuccessful.
 * It returns normally otherwise.
 * @param command The command line string that will be passed through to the <i>system()</i> call
 * @return The exit code of the background application
 */
int MainDialog::simulate( const QString& command)
{
    int ret = system(command.ascii());
    if (ret < 0) throw std::runtime_error(msgErrSysCall);
    return ret;
}

/**
 * This is the event handler for the radio button group on the geometry selection page.
 * It reacts by displaying the available information about the chosen geometry in the diagram and description areas.
 * It also adds the immutable parts of the parameter list to <i>cmdLineStub</i> and enables the <i>Next</i> navigation button.
 * @param radiobuttonid The index of the selected radio button; at the same time, the required index for <i>geometryTable</i>
 */
void MainDialog::geometryPicked( int radiobuttonid )
{
    try {
	geometryLayoutBox->setPixmap(geometryTable.at(radiobuttonid).layoutImage);
	geometryInfoBox->setText(geometryTable.at(radiobuttonid).layoutDescription);
	cmdLineStub = geometryTable.at(radiobuttonid).cmdLineParams;
	nextButton->setEnabled(TRUE);
    }
    catch (std::out_of_range oor) {
	std::cout << msgTemplateError << oor.what() << std::endl;
    }
}

/**
 * This is the event handler for the <i>Settings...</i> button on the parameter page.
 * It reacts by opening the popup menu that contains the options for saving and loading settings files.
 */
void MainDialog::settingsDialog()
{
    settingsButton->openPopup();
}

/**
 * This is the event handler for the <i>Default Values</i> option in the settings popup menu.
 * It reacts by resetting the widgets on the parameter page to their initial defaults (stored in <i>widgetCache</i>).
 * It also reloads the information for the selected geometry on the geometry selection page and resets <i>cmdLineStub</i>.
 */
void MainDialog::clearParameters()
{
    try {
	geometryLayoutBox->setPixmap(geometryTable.at(geometryPicker->selectedId()).layoutImage);
	geometryInfoBox->setText(geometryTable.at(geometryPicker->selectedId()).layoutDescription);
	cmdLineStub = geometryTable.at(geometryPicker->selectedId()).cmdLineParams;
	defaultsFromCache(parameterTable.at(geometryPicker->selectedId()), geometryPicker->selectedId());
	statusBar->setText(msgParamReset);
    }
    catch (std::out_of_range oor) {
	std::cout << oor.what() << std::endl;
	cmdLineStub = " ";
	geometryLayoutBox->clear();
	geometryInfoBox->clear();
	if (geometryPicker->selected()) {
	    static_cast<QRadioButton *>(geometryPicker->selected())->setChecked(FALSE);
	}
	nextButton->setEnabled(FALSE);
	backOnePage();
    }
}

/**
 * This is the event handler for the <i>Load settings...</i> option in the settings popup menu.
 * It reacts by opening a file dialog and attempting to read the contents of the chosen settings file.
 */
void MainDialog::loadSettingsDialog()
{
    QString settingsPath;
    settingsPath = buildSettingsPath(settingsPath);
    QString fromFile = QFileDialog::getOpenFileName(settingsPath, QString::null, this, "load settings",
						    QString("Load Settings from File"));
    if (!fromFile.isEmpty()) {
	try {
	    QFile workingFile(fromFile);
	    readParametersFromFile(workingFile, parameterTable.at(geometryPicker->selectedId()));
	}
	catch (std::runtime_error re) {
	    statusBar->setText(re.what());
	}
	catch (std::out_of_range oor) {
	    QString statustext(msgErrParamTableAccess);
	    statustext += oor.what();
	    statusBar->setText(statustext);
	}
    }
}

/**
 * This is the event handler for the <i>Save settings...</i> option in the settings popup menu.
 * It reacts by opening a file dialog and writing the contents of the parameter page to a file with the chosen name.
 */
void MainDialog::saveSettingsDialog()
{
    if (validateInput()) {
	QString settingsPath;
	settingsPath = buildSettingsPath(settingsPath);
	QString toFile = QFileDialog::getSaveFileName(settingsPath, QString::null, this, "save settings",
						      QString("Save Settings to File"));
	if (!toFile.isEmpty()) {
	    try {
		QFile workingFile(toFile);
		writeParametersToFile(workingFile, parameterTable.at(geometryPicker->selectedId()));
	    }
	    catch (std::runtime_error re) {
		statusBar->setText(re.what());
	    }
	    catch (std::out_of_range oor) {
		QString statustext(msgErrParamTableAccess);
		statustext += oor.what();
		statusBar->setText(statustext);
	    }
	}
    }
}

/**
 * This is the event handler for the <i>Save As Default</i> option in the settings popup menu.
 * It reacts by first backing up the original <i>defaultsettings</i> in a hidden file if necessary.
 * It then replaces <i>defaultsettings</i> by the contents of the parameter page.
 */
void MainDialog::overwriteDefaultSettings()
{
    if (validateInput()) {
	try {
	    QString settingsPath;
	    settingsPath = buildSettingsPath(settingsPath);
	    QDir settingsDir(settingsPath);
	    QFile workingFile(settingsDir.canonicalPath() + "/" + cDefaultSettings);
	    QFile backupFile(settingsDir.canonicalPath() + "/." + cDefaultSettings);
	    defaultsToCache(parameterTable.at(geometryPicker->selectedId()), geometryPicker->selectedId());
	    writeParametersToFile(workingFile, parameterTable.at(geometryPicker->selectedId()));
	    if (!settingsDir.exists("." + cDefaultSettings, FALSE)) {
		copyTextFile(workingFile, backupFile, QString(msgDefaultOverride));
	    }
	}
	catch (std::runtime_error re) {
	    statusBar->setText(re.what());
	}
	catch (std::out_of_range oor) {
	    QString statustext(msgErrParamTableAccess);
	    statustext += oor.what();
	    statusBar->setText(statustext);
	}
    }
}

/**
 * This is the event handler for the <i>Restore Default</i> option in the settings popup menu.
 * It reacts by looking up the hidden file <i>.defaultsettings</i> and copying its contents to <i>defaultsettings</i>.
 * If <i>.defaultsettings</i> is not found, it throws a runtime exception.
 */
void MainDialog::restoreDefaultSettings()
{
    try {
	QString settingsPath;
	settingsPath = buildSettingsPath(settingsPath);
	QDir settingsDir(settingsPath);
	if (settingsDir.exists("." + cDefaultSettings, FALSE)) {
	    QFile backupFile(settingsDir.canonicalPath() + "/." + cDefaultSettings);
	    QFile workingFile(settingsDir.canonicalPath() + "/" + cDefaultSettings);
	    copyTextFile(backupFile, workingFile, msgSettingsRestored);
	}
	else throw std::runtime_error(msgErrSettingsRestore);
    }
    catch (std::runtime_error re) {
	statusBar->setText(re.what());
    }
}

/**
 * The process of copying the contents of one file to another (which is overwritten) is bundled here for convenience.
 * The function also displays a message in the status bar confirming the process.
 * If there are errors opening either file, a runtime exception is thrown.
 * @param inFile
 * @param outFile
 * @param msg
 */
void MainDialog::copyTextFile(QFile& inFile, QFile& outFile, QString msg)
{
    if (inFile.open(IO_ReadOnly)) {
	if (outFile.open(IO_WriteOnly)) {
	    QTextStream fromStream(&inFile);
	    QTextStream toStream(&outFile);
	    toStream << fromStream.read();
	    outFile.close();
	    inFile.close();
	    statusBar->setText(msg);
	}
	else {
	    inFile.close();
	    throw std::runtime_error(msgErrWriteFile);
	}
    }
    else throw std::runtime_error(msgErrReadFile);
}

/**
 * This is the event handler that deals with activation of one of the items in the settings popup menu.
 * @param button The index of the selected menu item
 */
void MainDialog::settingsChanged(int button)
{
    statusBar->clear();
    switch (button) {
    case 0 : clearParameters();
	break;
    case 1 : loadSettingsDialog();
	break;
    case 2 : saveSettingsDialog();
	break;
    case 3 : overwriteDefaultSettings();
	break;
    case 4 : restoreDefaultSettings();
	break;
    }
}

/**
 * Building the path to the settings directory from internal variables and the name of the current geometry is bundled in here.
 * @param sPath A reference to the string variable that will hold the assembled path
 * @return The modified input parameter <i>sPath</i>
 */
QString& MainDialog::buildSettingsPath(QString& sPath)
{
    sPath = basePath + resExtension +"/";
    sPath += QString::number(geometryPicker->selectedId() + 1);
    sPath += settingsExtension;
    return sPath;
}

/**
 * This function attempts to read a set of parameters for the parameter page from a file.
 * It will throw a runtime exception if it cannot read the file
 * @param readFile A reference to the input file
 * @param paramrow The destination struct that absorbs the information found in the file
 */
void MainDialog::readParametersFromFile(QFile& readFile, paramaggreg& paramrow)
{
    if (readFile.open(IO_ReadOnly)) {
	QTextStream instream(&readFile);
	QString contents = instream.read();
	if (contents && !contents.isEmpty()) {
	    contents = contents.stripWhiteSpace().simplifyWhiteSpace();
	    QStringList contentlist = QStringList::split(" ", contents);
	    QStringList::const_iterator iter = contentlist.begin();
	    paramrow.trackerName = *iter;
	    iter++;
	    if (iter != contentlist.end()) { paramrow.nlayers = (*iter).toInt(); iter++; }
	    int i = 0;
	    while (i < paramrow.nlayers && iter != contentlist.end()) {
		paramrow.nchipsperside.push_back((*iter).toInt());
		iter++;
		i++;
	    }
	    i = 0;
	    while (i < paramrow.nlayers && iter != contentlist.end()) {
		switch ((*iter).toInt()) {
		case 0 : paramrow.mtypeslayers.push_back(ss);
		    break;
		case 1 : paramrow.mtypeslayers.push_back(ds);
		    break;
		case 2 : paramrow.mtypeslayers.push_back(pt);
		    break;
		default: paramrow.mtypeslayers.push_back(ss);
		}
		iter++;
		i++;
	    }
	    if (iter != contentlist.end()) { paramrow.nrings = (*iter).toInt(); iter++; }
	    i = 0;
	    while (i < paramrow.nrings && iter != contentlist.end()) {
		paramrow.nsegmentsperring.push_back((*iter).toInt());
		iter++;
		i++;
	    }
	    i = 0;
	    while (i < paramrow.nrings && iter != contentlist.end()) {
		switch ((*iter).toInt()) {
		case 0 : paramrow.mtypesrings.push_back(ss);
		    break;
		case 1 : paramrow.mtypesrings.push_back(ds);
		    break;
		case 2 : paramrow.mtypesrings.push_back(pt);
		    break;
		default: paramrow.mtypesrings.push_back(ss);
		}
		iter++;
		i++;
	    }
	    if (iter != contentlist.end()) { paramrow.costpersqcm = (*iter).toDouble(); iter++; }
	    if (iter != contentlist.end()) { paramrow.ptcostpersqcm = (*iter).toDouble(); iter++; }
	    if (iter != contentlist.end()) { paramrow.powerperchannel = (*iter).toDouble(); iter++; }
	    if (iter != contentlist.end()) paramrow.ptpowerperchannel = (*iter).toDouble();
	    statusBar->setText(msgParamsRead);
	}
	else statusBar->setText(msgErrInputFileParse);
	readFile.close();
    }
    else throw std::runtime_error(msgErrReadFile);
}

/**
 * This function attempts to write the contents of the parameter page to a file (which is overwritten).
 * If an error occurs while opening or creating the output file, it throws a runtime exception.
 * @param writeFile A reference to the destination file
 * @param paramrow The struct containing the information about the selected geometry's current parameter values
 */
void MainDialog::writeParametersToFile(QFile& writeFile, const paramaggreg& paramrow)
{
    if (writeFile.open(IO_WriteOnly)) {
	QString contents = paramrow.trackerName + " " + QString::number(paramrow.nlayers);
	for (int i = 0; i < paramrow.nlayers; i++) {
	    contents = contents + " " + QString::number(paramrow.nchipsperside.at(i));
	}
	for (int i = 0; i < paramrow.nlayers; i++) {
	     contents = contents + " " + QString::number(paramrow.mtypeslayers.at(i));
	}
	contents = contents + " " + QString::number(paramrow.nrings);
	for (int i = 0; i < paramrow.nrings; i++) {
	    contents = contents + " " + QString::number(paramrow.nsegmentsperring.at(i));
	}
	for (int i = 0; i < paramrow.nrings; i++) {
	    contents = contents + " " + QString::number(paramrow.mtypesrings.at(i));
	}
	contents = contents + " " + QString::number(paramrow.costpersqcm);
		   contents = contents + " " + QString::number(paramrow.ptcostpersqcm);
		   contents = contents + " " + QString::number(paramrow.powerperchannel);
		   contents = contents + " " + QString::number(paramrow.ptpowerperchannel);
	QTextStream outstream(&writeFile);
	outstream << contents;
	writeFile.close();
	statusBar->setText(msgParamsWritten);
    }
    else throw std::runtime_error(msgErrWriteFile);
}

/**
 * Copying the parameter cache entry at a given position to an instance of <i>paramaggreg</i> is bundled here for convenience.
 * @param paramrow A reference to the destination object
 * @param pos The index pointing to the entry in <i>widgetCache</i> that contains the requested information
 */
void MainDialog::defaultsFromCache(paramaggreg& paramrow, int pos)
{
    try {
	paramrow.trackerName = widgetCache.at(pos).trackerName;
	paramrow.nlayers = widgetCache.at(pos).nlayers;
	paramrow.nrings = widgetCache.at(pos).nrings;
	paramrow.nchipsperside = widgetCache.at(pos).nchipsperside;
	paramrow.nsegmentsperring = widgetCache.at(pos).nsegmentsperring;
	paramrow.mtypeslayers = widgetCache.at(pos).mtypeslayers;
	paramrow.mtypesrings = widgetCache.at(pos).mtypesrings;
	paramrow.costpersqcm = widgetCache.at(pos).costpersqcm;
	paramrow.ptcostpersqcm = widgetCache.at(pos).ptcostpersqcm;
	paramrow.powerperchannel = widgetCache.at(pos).powerperchannel;
	paramrow.ptpowerperchannel = widgetCache.at(pos).ptpowerperchannel;
    }
    catch (std::out_of_range oor) {
	std::cout << msgErrParamCacheAccess << oor.what() << std::endl;
    }
}

/**
 * Copying the contents of an instance of <i>paramggreg</i> to an entry in the parameter cache is bundled here for convenience.
 * @param paramrow A reference to the information source
 * @param pos The index pointing to the entry in <i>widgetCache</i> that will receive the data
 */
void MainDialog::defaultsToCache(const paramaggreg& paramrow, int pos)
{
    try {
	widgetCache.at(pos).trackerName = paramrow.trackerName;
	widgetCache.at(pos).nlayers = paramrow.nlayers;
	widgetCache.at(pos).nrings = paramrow.nrings;
	widgetCache.at(pos).nchipsperside = paramrow.nchipsperside;
	widgetCache.at(pos).nsegmentsperring = paramrow.nsegmentsperring;
	widgetCache.at(pos).mtypeslayers = paramrow.mtypeslayers;
	widgetCache.at(pos).mtypesrings = paramrow.mtypesrings;
	widgetCache.at(pos).costpersqcm = paramrow.costpersqcm;
	widgetCache.at(pos).ptcostpersqcm = paramrow.ptcostpersqcm;
	widgetCache.at(pos).powerperchannel = paramrow.powerperchannel;
	widgetCache.at(pos).ptpowerperchannel = paramrow.ptpowerperchannel;
    }
    catch (std::out_of_range oor) {
	std::cout << msgErrParamCacheAccess << oor.what() << std::endl;
    }
}

/**
 * Copying the contents of an instance of <i>paramaggreg</i> to the parameter page widgets is bundled here for convenience.
 * @param paramrow A reference to the information source
 */
void MainDialog::valuesToWidgets(const paramaggreg& paramrow)
{
    trackerNameLineEdit->setText(paramrow.trackerName);
    QString layersrings = QString::number(paramrow.nlayers);
    layersrings = "<b>" + layersrings + "</b>";
    mLabel->setText(layersrings);
    layersrings = QString::number(paramrow.nrings);
    layersrings = "<b>" + layersrings + "</b>";
    nLabel->setText(layersrings);
    if (layerSelection->count() > 0) layerSelected(0);
    if (ringSelection->count() > 0) ringSelected(0);
    costPerSqCmEdit->setText(QString::number(paramrow.costpersqcm));
    costPtPerSqCmEdit->setText(QString::number(paramrow.ptcostpersqcm));
    powerEdit->setText(QString::number(paramrow.powerperchannel));
    ptPowerEdit->setText(QString::number(paramrow.ptpowerperchannel));
}

/**
 * This is the event handler that translates a selection in the ring module type listbox to an entry in <i>parameterTable</i>.
 * @param index The index of the selected list item
 */
void MainDialog::ringTypeSelected(int index)
{
    try {
	switch (index) {
	case 0 : parameterTable.at(geometryPicker->selectedId())
		    .mtypesrings.at(ringSelection->currentItem()) = ss;
	    break;
	case 1 : parameterTable.at(geometryPicker->selectedId())
		    .mtypesrings.at(ringSelection->currentItem()) = ds;
	    break;
	case 2 : parameterTable.at(geometryPicker->selectedId())
		    .mtypesrings.at(ringSelection->currentItem()) = pt;
	    break;
	}
    }
    catch (std::out_of_range oor) {
	QString statustext(msgErrParamTableAccess);
	statustext += oor.what();
	statusBar->setText(statustext);
    }
}

/**
 * This is the event handler that translates a selection in the module type listbox to an entry in <i>parameterTable</i>.
 * @param index The index of the selected list item
 */
void MainDialog::layerTypeSelected(int index)
{
    try {
	switch (index) {
	case 0 : parameterTable.at(geometryPicker->selectedId())
		    .mtypeslayers.at(layerSelection->currentItem()) = ss;
	    break;
	case 1 : parameterTable.at(geometryPicker->selectedId())
	    .mtypeslayers.at(layerSelection->currentItem()) = ds;
	    break;
	case 2 : parameterTable.at(geometryPicker->selectedId())
	    .mtypeslayers.at(layerSelection->currentItem()) = pt;
	    break;
	}
    }
    catch (std::out_of_range oor) {
	QString statustext(msgErrParamTableAccess);
	statustext += oor.what();
	statusBar->setText(statustext);
    }
}

/**
 * This is the event handler that updates the connected spinner and listbox when a layer is selected on the parameter page.
 * @param index The index of the selected combobox item
 */
void MainDialog::layerSelected( int index)
{
    try {
	QString totalchips = QString::number(parameterTable.at(geometryPicker->selectedId())
					     .nchipsperside.at(index));
	chipsPerSideSpinner->setValue(totalchips.toInt() / cChipModulus);
	totalchips = "<b>" + totalchips + "</b>";
	totalChipsLabel->setText(totalchips);
	int idx;
	switch (parameterTable.at(geometryPicker->selectedId()).mtypeslayers.at(index)) {
	case ss : idx = 0;
	    break;
	case ds : idx = 1;
	    break;
	case pt : idx = 2;
	    break;
	default : idx = 0;
	}
	layerTypeListBox->setCurrentItem(idx);
    }
    catch (std::out_of_range oor) {
	QString statustext(msgErrParamTableAccess);
	statustext += oor.what();
	statusBar->setText(statustext);
    }
}

/**
 * This is the event handler that updates the connected spinner and listbox when a ring is selected on the parameter page.
 * @param index The index of the selected combobox item
 */
void MainDialog::ringSelected( int index)
{
    try {
	segmentsPerRingSpinner->setValue(
		parameterTable.at(geometryPicker->selectedId()).nsegmentsperring.at(index));
	int idx;
	switch (parameterTable.at(geometryPicker->selectedId()).mtypesrings.at(index)) {
	case ss : idx = 0;
	    break;
	case ds : idx = 1;
	    break;
	case pt : idx = 2;
	    break;
	default : idx = 0;
	}
	ringTypeListBox->setCurrentItem(idx);
    }
    catch (std::out_of_range oor) {
	QString statustext(msgErrParamTableAccess);
	statustext += oor.what();
	statusBar->setText(statustext);
    }    
}

/**
 * This is the event handler that validates the input when the value of the <i>chips per side</i> spinner changes.
 * It also updates the value of the label that displays the final number of chips per side, and the entry in <i>parameterTable</i>.
 * @param value The new value that was entered into the spinner
 */
void MainDialog::chipsPerSideChanged(int value)
{
    try {
	int pos = 0;
	QString inttostring = chipsPerSideSpinner->text();
	switch (chipsPerSideSpinner->validator()->validate(inttostring, pos)) {
	case QValidator::Acceptable : parameterTable.at(geometryPicker->selectedId())
		    .nchipsperside.at(layerSelection->currentItem()) = value * cChipModulus;
	    {
		QString totalchips = QString::number(parameterTable.at(geometryPicker->selectedId())
						     .nchipsperside.at(layerSelection->currentItem()));
		totalchips = "<b>" + totalchips + "</b>";
		totalChipsLabel->setText(totalchips);
	    }
	    break;
	case QValidator::Intermediate : statusBar->setText(msgSpinValidationFuzzy);
	    break;
	case QValidator::Invalid : statusBar->setText(msgSpinValidationError);
	    break;
	default : statusBar->setText(msgErrValidationStrange);
    }
    }
    catch (std::out_of_range oor) {
	QString statustext(msgErrParamTableAccess);
	statustext += oor.what();
	statusBar->setText(statustext);
    }   
}

/**
 * This is the event handler that validates the input when the value of the <i>segments per ring</i> spinner changes.
 * It also updates the corresponding entry in <i>parameterTable</i>.
 * @param value The new value that was entered into the spinner
 */
void MainDialog::segmentsPerRingChanged(int value)
{
    try {
	int pos = 0;
	QString inttostring = segmentsPerRingSpinner->text();
	switch (segmentsPerRingSpinner->validator()->validate(inttostring, pos)) {
	case QValidator::Acceptable : parameterTable.at(geometryPicker->selectedId())
		    .nsegmentsperring.at(layerSelection->currentItem()) = value;
	    break;
	case QValidator::Intermediate : statusBar->setText(msgSpinValidationFuzzy);
	    break;
	case QValidator::Invalid : statusBar->setText(msgSpinValidationError);
	    break;
	default : statusBar->setText(msgErrValidationStrange);
              }
    }
    catch (std::out_of_range oor) {
	QString statustext(msgErrParamTableAccess);
	statustext += oor.what();
	statusBar->setText(statustext);
    }      
}

/**
 * This function validates all input that may not be in the correct format when the <i>Go</i> button is clicked.
 * If it doesn't accept the input, it displays a message in the status bar to give a clue about what may have happened.
 * @return success or failure; a result of <i>Intermediate</i> is considered a failure
 */
bool MainDialog::validateInput()
{
    try {
	int pos = 0;
	QString inttostring = chipsPerSideSpinner->text();
	switch (chipsPerSideSpinner->validator()->validate(inttostring, pos)) {
	    case QValidator::Acceptable : break;
	    case QValidator::Intermediate : statusBar->setText(msgSpinValidationFuzzy);
		return FALSE;
	    case QValidator::Invalid : statusBar->setText(msgSpinValidationError);
		return FALSE;
	    default : statusBar->setText(msgErrValidationStrange);
		return FALSE;
	}
	inttostring = segmentsPerRingSpinner->text();
	switch (segmentsPerRingSpinner->validator()->validate(inttostring, pos)) {
	    case QValidator::Acceptable : break;
	    case QValidator::Intermediate : statusBar->setText(msgSpinValidationFuzzy);
		return FALSE;
	    case QValidator::Invalid : statusBar->setText(msgSpinValidationError);
		return FALSE;
	    default : statusBar->setText(msgErrValidationStrange);
		return FALSE;
	}
	if (!costPerSqCmEdit->hasAcceptableInput()) {
	    statusBar->setText(msgValidationError);
	    return FALSE;
	}
	if (!costPtPerSqCmEdit->hasAcceptableInput()) {
	    statusBar->setText(msgValidationError);
	    return FALSE;
	}
	if (!powerEdit->hasAcceptableInput()) {
	    statusBar->setText(msgValidationError);
	    return FALSE;
	}
	if (!ptPowerEdit->hasAcceptableInput()) {
	    statusBar->setText(msgValidationError);
	    return FALSE;
	}
	return TRUE;
    }
    catch (std::out_of_range oor) {
	QString statustext(msgErrParamTableAccess);
	statustext += oor.what();
	statusBar->setText(statustext);
	return FALSE;
    }   
}

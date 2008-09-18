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
 * This is the initialisation function that serves as a programmer-defined addition to the constructor.
 * (The constructor itself is private and auto-generated.)
 * It initialises some of the widgets programmatically, reads the information about the available geometries from the
 * resource directory and caches it in the global table vectors <i>geometryTable , parameterTable</i> and <i>widgetCache</i>.
 * If a resource directory is found to contain nonsensical parameter files it is skipped.
 */
void MainDialog::init()
{
    /// Instantiate the file handler helper class
    fh = new FileHandler();
    
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
    QIntValidator *intval = new QIntValidator(layerChipsSpinner);
    intval->setBottom(cPositiveNumbers);
    intval->setTop(cMaxChipsInSpinner);
    layerChipsSpinner->setValidator(intval);
    
    intval = new QIntValidator(layerSegmentsSpinner);
    intval->setBottom(cPositiveNumbers);
    layerSegmentsSpinner->setValidator(intval);
    
    intval = new QIntValidator(ringChipsSpinner);
    intval->setBottom(cPositiveNumbers);
    intval->setTop(cMaxChipsInSpinner);
    ringChipsSpinner->setValidator(intval);
    
    intval = new QIntValidator(ringSegmentsSpinner);
    intval->setBottom(cPositiveNumbers);
    ringSegmentsSpinner->setValidator(intval);
    
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
    outDirExtension = cTmpDir;

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

	/// Set the location of the initial configuration file
	workingFile.setName(workingDir.canonicalPath() + "/" + *iter + "/" + cConfigFileName);
	if (workingFile.exists()) {
	    geomrow.configFile = workingFile.name();
	}
	else continue;

	/// Read the cusomisable parameters from the config file
	paramaggreg paramrow;
	workingFile.setName(workingDir.canonicalPath() + "/" + *iter + cSettingsExtension + "/" + cDefaultSettings);
	try {
	    if (workingFile.exists()) {
		fh->readParametersFromFile(workingFile, paramrow);
		QFileInfo fi1(geomrow.configFile);
		QFileInfo fi2(workingFile);
		if (fi1.size() == fi2.size()) {
		    fh->writeParametersToFile(workingFile, paramrow, summaryExtension + outDirExtension);
		}
	    }
	    else continue;
	}
	catch (std::runtime_error re) {
	    std::cout << re.what() << std::endl;
	    continue;
	}
	
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
    fh->removeOutputDir(basePath + summaryExtension + outDirExtension);
    fh->removeTmpConfigFile(tmpConfig.name());
    delete fh;
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
    fh->removeTmpConfigFile(tmpConfig.name());
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
    fh->removeOutputDir(basePath + summaryExtension + outDirExtension);
    fh->removeTmpConfigFile(tmpConfig.name());
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
	tmpConfig.setName(basePath + "/" + cTempConfig);
	QString configpath = basePath + resExtension + "/" + geometryPicker->selected()->text();
	configpath = configpath + settingsExtension + "/";
	QFile config(configpath + cDefaultSettings);
	fh->copyTextFile(config, tmpConfig);
	for (uint i = 0; i < parameterTable.at(geometryPicker->selectedId()).barrelnames.size(); i++) {
	    barrelSelection->insertItem(parameterTable.at(geometryPicker->selectedId()).barrelnames.at(i), i);
	}
	for (uint i = 0; i < parameterTable.at(geometryPicker->selectedId()).endcapnames.size(); i++) {
	    endcapSelection->insertItem(parameterTable.at(geometryPicker->selectedId()).endcapnames.at(i), i);
	}
	valuesToWidgets(parameterTable.at(geometryPicker->selectedId()));
	barrelSelected(0);
	endcapSelected(0);
    }
    catch (std::out_of_range oor) {
	QString statusText(msgErrValueLoad);
	statusText += oor.what();
	statusBar->setText(statusText);
    }
    catch (std::runtime_error re) {
	std::cout << re.what() << std::endl;
	std::cout << msgCriticalErrorConfigFile << std::endl;
	exit(-1);
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
    fh->removeOutputDir(basePath + summaryExtension + outDirExtension);
    layerSelection->clear();
    ringSelection->clear();
    removeTmpConfigFile(tmpConfig.name());
    clearParameters();
    mainWidgetStack->raiseWidget(0);
}


void MainDialog::go()
{
    if (validateInput()) {
	QString command = basePath + cCommand + " ";
	fh->writeParametersToFile(tmpConfig, parameterTable.at(geometryPicker->selectedId()),
				  summaryExtension + outDirExtension);
	cmdLineStub = tmpConfig.name();
	command += cmdLineStub;
	try  {
	    int simulstatus;	
	    simulstatus = simulate(command);
	    std::cout << msgSimulationExit << simulstatus << "." << std::endl;
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
	    statusBar->setText(re.what());
	}
    }
}

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
	cmdLineStub = "";
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
	    fh->readParametersFromFile(workingFile, parameterTable.at(geometryPicker->selectedId()));
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
		fh->writeParametersToFile(workingFile, parameterTable.at(geometryPicker->selectedId()),
					  summaryExtension + outDirExtension);
		statusBar->setText(msgParamsWritten);
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
	    fh->writeParametersToFile(workingFile, parameterTable.at(geometryPicker->selectedId()),
				      summaryExtension + outDirExtension);
	    statusBar->setText(msgParamsWritten);
	    if (!settingsDir.exists("." + cDefaultSettings, FALSE)) {
		fh->copyTextFile(workingFile, backupFile);
		statusBar->setText(msgDefaultOverride);
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
	    fh->copyTextFile(backupFile, workingFile);
	    statusBar->setText(msgSettingsRestored);
	}
	else throw std::runtime_error(msgErrSettingsRestore);
    }
    catch (std::runtime_error re) {
	statusBar->setText(re.what());
    }
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
	paramrow.barrelnames = widgetCache.at(pos).barrelnames;
	paramrow.endcapnames = widgetCache.at(pos).endcapnames;
	paramrow.nchipslayer = widgetCache.at(pos).nchipslayer;
	paramrow.nchipsring = widgetCache.at(pos).nchipsring;
	paramrow.nsegmentslayer = widgetCache.at(pos).nsegmentslayer;
	paramrow.nsegmentsring = widgetCache.at(pos).nsegmentsring;
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
	widgetCache.at(pos).barrelnames = paramrow.barrelnames;
	widgetCache.at(pos).endcapnames = paramrow.endcapnames;
	widgetCache.at(pos).nchipslayer = paramrow.nchipslayer;
	widgetCache.at(pos).nchipsring = paramrow.nchipsring;
	widgetCache.at(pos).nsegmentslayer = paramrow.nsegmentslayer;
	widgetCache.at(pos).nsegmentsring = paramrow.nsegmentsring;
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
    if (barrelSelection->count() > 0) barrelSelected(0);
    if (endcapSelection->count() > 0) endcapSelected(0);
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
		    .mtypesrings.at(endcapSelection->currentItem()).at(ringSelection->currentItem()) = rphi;
	    break;
	case 1 : parameterTable.at(geometryPicker->selectedId())
		    .mtypesrings.at(endcapSelection->currentItem()).at(ringSelection->currentItem()) = stereo;
	    break;
	case 2 : parameterTable.at(geometryPicker->selectedId())
		    .mtypesrings.at(endcapSelection->currentItem()).at(ringSelection->currentItem()) = pt;
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
		    .mtypeslayers.at(barrelSelection->currentItem()).at(layerSelection->currentItem()) = rphi;
	    break;
	case 1 : parameterTable.at(geometryPicker->selectedId())
	    .mtypeslayers.at(barrelSelection->currentItem()).at(layerSelection->currentItem()) = stereo;
	    break;
	case 2 : parameterTable.at(geometryPicker->selectedId())
	    .mtypeslayers.at(barrelSelection->currentItem()).at(layerSelection->currentItem()) = pt;
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
	QString total = QString::number(parameterTable.at(geometryPicker->selectedId())
					.nchipslayer.at(barrelSelection->currentItem()).at(index));
	layerChipsSpinner->setValue(total.toInt() / cLayerChipModulus);
	total = "<b>" + total + "</b>";
	layerTotalChipsLabel->setText(total);
	layerSegmentsSpinner->setValue(parameterTable.at(geometryPicker->selectedId())
				     .nsegmentslayer.at(barrelSelection->currentItem()).at(index));
	int idx;
	switch (parameterTable.at(geometryPicker->selectedId())
		.mtypeslayers.at(barrelSelection->currentItem()).at(index)) {
	case rphi : idx = 0;
	    break;
	case stereo : idx = 1;
	    break;
	case pt : idx = 2;
	    break;
	default : idx = -1;
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
	QString total = QString::number(parameterTable.at(geometryPicker->selectedId())
					.nchipsring.at(endcapSelection->currentItem()).at(index));
	ringChipsSpinner->setValue(total.toInt() / cRingChipModulus);
	total = "<b>" + total + "</b>";
	ringTotalChipsLabel->setText(total);
	ringSegmentsSpinner->setValue(parameterTable.at(geometryPicker->selectedId())
				      .nsegmentsring.at(endcapSelection->currentItem()).at(index));
	int idx;
	switch (parameterTable.at(geometryPicker->selectedId())
		.mtypesrings.at(endcapSelection->currentItem()).at(index)) {
	case rphi : idx = 0;
	    break;
	case stereo : idx = 1;
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

void MainDialog::barrelSelected(int index)
{
    try {
	layerSelection->clear();
	for (int i = 0; i < parameterTable.at(geometryPicker->selectedId()).nlayers.at(index); i++) {
	    layerSelection->insertItem(QString::number(i + 1), i);
	}
	layerSelection->setCurrentItem(0);
	layerSelected(0);
    }
    catch (std::out_of_range oor) {
	QString statustext(msgErrParamTableAccess);
	statustext += oor.what();
	statusBar->setText(statustext);
    }
}

void MainDialog::endcapSelected(int index)
{
    try {
	ringSelection->clear();
	if (parameterTable.at(geometryPicker->selectedId()).nrings.at(index) > 0) {
	    for (int i = 0; i < parameterTable.at(geometryPicker->selectedId()).nrings.at(index); i++) {
		ringSelection->insertItem(QString::number(i + 1), i);
	    }
	    ringSelection->setCurrentItem(0);
	    ringSelected(0);
	}
	else {
	    ringChipsSpinner->setEnabled(FALSE);
	    ringSegmentsSpinner->setEnabled(FALSE);
	    ringTypeListBox->setEnabled(FALSE);
	}
    }
    catch (std::out_of_range oor) {
	QString statustext(msgErrParamTableAccess);
	statustext += oor.what();
	statusBar->setText(statustext);
    }
}

void MainDialog::addRing()
{
    ringSelection->insertItem(QString::number(ringSelection->count() + 1), ringSelection->count());
    parameterTable.at(geometryPicker->selectedId()).nrings.at(endcapSelection->currentItem())++;
    parameterTable.at(geometryPicker->selectedId()).nchipsring.at(endcapSelection->currentItem())
	    .push_back(cRingChipModulus);
    parameterTable.at(geometryPicker->selectedId()).nsegmentsring.at(endcapSelection->currentItem())
	    .push_back(ringSegmentsSpinner->minValue());
    parameterTable.at(geometryPicker->selectedId()).mtypesrings.at(endcapSelection->currentItem())
	    .push_back(rphi);
    widgetCache.at(geometryPicker->selectedId()).nrings.at(endcapSelection->currentItem())++;
    widgetCache.at(geometryPicker->selectedId()).nchipsring.at(endcapSelection->currentItem())
	    .push_back(cRingChipModulus);
    widgetCache.at(geometryPicker->selectedId()).nsegmentsring.at(endcapSelection->currentItem())
	    .push_back(ringSegmentsSpinner->minValue());
    widgetCache.at(geometryPicker->selectedId()).mtypesrings.at(endcapSelection->currentItem())
	    .push_back(rphi);
    if (ringSelection->count() == 1) {
	ringChipsSpinner->setEnabled(TRUE);
	ringSegmentsSpinner->setEnabled(TRUE);
	ringTypeListBox->setEnabled(TRUE);
    }
    ringSelection->setCurrentItem(ringSelection->count() - 1);
    ringSelected(ringSelection->currentItem());
}

void MainDialog::removeRing()
{
    if (ringSelection->count() > 0) {
	ringSelection->removeItem(ringSelection->count() - 1);
	parameterTable.at(geometryPicker->selectedId()).nrings.at(endcapSelection->currentItem())--;
	parameterTable.at(geometryPicker->selectedId()).nchipsring
		.at(endcapSelection->currentItem()).pop_back();
	parameterTable.at(geometryPicker->selectedId()).nsegmentsring
		.at(endcapSelection->currentItem()).pop_back();
	parameterTable.at(geometryPicker->selectedId()).mtypesrings
		.at(endcapSelection->currentItem()).pop_back();
	widgetCache.at(geometryPicker->selectedId()).nrings.at(endcapSelection->currentItem())--;
	widgetCache.at(geometryPicker->selectedId()).nchipsring
		.at(endcapSelection->currentItem()).pop_back();
	widgetCache.at(geometryPicker->selectedId()).nsegmentsring
		.at(endcapSelection->currentItem()).pop_back();
	widgetCache.at(geometryPicker->selectedId()).mtypesrings
		.at(endcapSelection->currentItem()).pop_back();
	if (ringSelection->count() < 1) {
	    ringChipsSpinner->setValue(ringChipsSpinner->minValue());
	    ringChipsSpinner->setEnabled(FALSE);
	    ringSegmentsSpinner->setValue(ringSegmentsSpinner->minValue());
	    ringSegmentsSpinner->setEnabled(FALSE);
	    ringTypeListBox->setCurrentItem(none);
	    ringTypeListBox->setEnabled(FALSE);
	}
	else {
	    ringSelection->setCurrentItem(ringSelection->count() -1);
	    ringSelected(ringSelection->currentItem());
	}
    }
}

/**
 * This is the event handler that validates the input when the value of the <i>chips per side</i> spinner changes.
 * It also updates the value of the label that displays the final number of chips per side, and the entry in <i>parameterTable</i>.
 * @param value The new value that was entered into the spinner
 */
void MainDialog::layerChipsAcrossChanged(int value)
{
    try {
	int pos = 0;
	QString inttostring = layerChipsSpinner->text();
	switch (layerChipsSpinner->validator()->validate(inttostring, pos)) {
	case QValidator::Acceptable : parameterTable.at(geometryPicker->selectedId())
		    .nchipslayer.at(barrelSelection->currentItem())
		    .at(layerSelection->currentItem()) = value * cLayerChipModulus;
	    {
		QString totalchips = QString::number(parameterTable.at(geometryPicker->selectedId())
						     .nchipslayer.at(barrelSelection->currentItem())
						     .at(layerSelection->currentItem()));
		totalchips = "<b>" + totalchips + "</b>";
		layerTotalChipsLabel->setText(totalchips);
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

void MainDialog::layerSegmentsAlongChanged(int value)
{
    try {
	int pos = 0;
	QString inttostring = layerSegmentsSpinner->text();
	switch (layerSegmentsSpinner->validator()->validate(inttostring, pos)) {
	case QValidator::Acceptable : parameterTable.at(geometryPicker->selectedId()).nsegmentslayer
		    .at(barrelSelection->currentItem()).at(layerSelection->currentItem()) = value;
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

void MainDialog::ringChipsAcrossChanged(int value)
{
    try {
	int pos = 0;
	QString inttostring = ringChipsSpinner->text();
	switch (ringChipsSpinner->validator()->validate(inttostring, pos)) {
	case QValidator::Acceptable : parameterTable.at(geometryPicker->selectedId())
		    .nchipsring.at(endcapSelection->currentItem())
		    .at(ringSelection->currentItem()) = value * cRingChipModulus;
	    {
		QString totalchips = QString::number(parameterTable.at(geometryPicker->selectedId())
						     .nchipsring.at(endcapSelection->currentItem())
						     .at(ringSelection->currentItem()));
		totalchips = "<b>" + totalchips + "</b>";
		ringTotalChipsLabel->setText(totalchips);
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
void MainDialog::ringSegmentsAlongChanged(int value)
{
    try {
	int pos = 0;
	QString inttostring = ringSegmentsSpinner->text();
	switch (ringSegmentsSpinner->validator()->validate(inttostring, pos)) {
	case QValidator::Acceptable : parameterTable.at(geometryPicker->selectedId()).nsegmentsring
		    .at(endcapSelection->currentItem()).at(ringSelection->currentItem()) = value;
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
	QString inttostring = layerChipsSpinner->text();
	switch (layerChipsSpinner->validator()->validate(inttostring, pos)) {
	    case QValidator::Acceptable : break;
	    case QValidator::Intermediate : statusBar->setText(msgSpinValidationFuzzy);
		return FALSE;
	    case QValidator::Invalid : statusBar->setText(msgSpinValidationError);
		return FALSE;
	    default : statusBar->setText(msgErrValidationStrange);
		return FALSE;
	}
	inttostring = layerSegmentsSpinner->text();
	switch (layerSegmentsSpinner->validator()->validate(inttostring, pos)) {
	    case QValidator::Acceptable : break;
	    case QValidator::Intermediate : statusBar->setText(msgSpinValidationFuzzy);
		return FALSE;
	    case QValidator::Invalid : statusBar->setText(msgSpinValidationError);
		return FALSE;
	    default : statusBar->setText(msgErrValidationStrange);
		return FALSE;
	}
	 inttostring = ringChipsSpinner->text();
	switch (ringChipsSpinner->validator()->validate(inttostring, pos)) {
	    case QValidator::Acceptable : break;
	    case QValidator::Intermediate : statusBar->setText(msgSpinValidationFuzzy);
		return FALSE;
	    case QValidator::Invalid : statusBar->setText(msgSpinValidationError);
		return FALSE;
	    default : statusBar->setText(msgErrValidationStrange);
		return FALSE;
	}
	inttostring = ringSegmentsSpinner->text();
	switch (ringSegmentsSpinner->validator()->validate(inttostring, pos)) {
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
	if (trackerNameLineEdit->text().length() == 0) {
	    trackerNameLineEdit->setText(cDefaultTrackerName);
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

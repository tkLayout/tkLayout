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
 * @file resultspopup.ui.h
 * @brief This file describes the popup that shows the results of <i>tkgeomgui</i> after a simulation has completed. It gives the option of saving those results in a directory chosen by the user
 */

/**
 * @class resultspopup
 * @brief This form creates a non-modal dialog popup that displays the simulation results.
 */

/**
 * A filter for the filetypes that need to be taken from the original results directory to the user-defined one
 */
static const QString saveFilter = "*.html *.png";

/**
 * This is the initialisation function that serves as a programmer-defined addition to the constructor.
 * (The constructor itself is private and auto-generated.)
 * It creates an instance of the <i>FileHandler</i> helper class that it needs to save the displayed results.
 */
void ResultsPopup::init()
{
    fh = new FileHandler();
}

/**
 * Since the destructor is as private and auto-generated as the constructor, cleanup of the <i>FileHandler</i> helper
 * class and removal of the results files that are no longer needed happens in here.
 */
void ResultsPopup::destroy()
{
    fh->removeOutputDir(resultsPath);
    QString rootfile = basePath + cRootDirExtension + "/" + trackerName + cRootFileExt;
    fh->removeTmpConfigFile(rootfile);
    delete fh;
}

/**
 * This is the event handler for the <i>Save...</i> button inside the popup.
 * It reacts by opening a file dialog that returns the name of the directory where the results files will be saved.
 * If that name is not null, it copies those results files that match the filter string <i>saveFileter</i> there.
 */
void ResultsPopup::saveResults()
{
    QString savePath = QFileDialog::getExistingDirectory(QString::null, this, "saveres", "Save in...");
    if (!savePath.isEmpty()) {
	try {
	    QDir resDir(resultsPath);
	    resDir.setNameFilter(saveFilter);
	    resDir.setFilter(QDir::Files | QDir::NoSymLinks);
	    QDir saveDir(savePath);
	    if (!saveDir.exists()) saveDir.mkdir(savePath);
	    QStringList fileNames = resDir.entryList();
	    for (QStringList::const_iterator iter = fileNames.begin(); iter != fileNames.end(); iter++) {
		QFile infile(resDir.canonicalPath() + "/" + *iter);
		QFile outfile(saveDir.canonicalPath() + "/" + *iter);
		fh->copyDataFile(infile, outfile);
	    }
	    statusBar->setText(msgResultsSaved);
	}
	catch (std::runtime_error re) {
	    statusBar->setText(re.what());
	}
    }
}

/**
 * This is a setter for the popup's internal copy of the base path.
 * @param path The designated base path as an instance of QString
 */
void ResultsPopup::setBasePath(QString path)
{
    basePath = path;
}

/**
 * This is a getter for the popup's internal copy of the base path.
 * @return A reference to the internal copy of the base path as an instance of QString
 */
QString& ResultsPopup::getResultsPath()
{
    return resultsPath;
}

/**
 * This is a setter for the popup's internal copy of the temporary results path.
 * @param path The designated temporary results path as an instance of QString
 */
void ResultsPopup::setResultsPath(QString path)
{
    resultsPath = path;
}

/**
 * This is a setter for the popup's internal copy of the tracker name.
 * It is used to find and delete the <i>.root</i> file output when the popup is closed.
 * @param path The tracker name as an instance of QString
 */
void ResultsPopup::setTrackerName(QString name)
{
    trackerName = name;
}

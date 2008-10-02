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
static const QString saveFilter = "*.html *.png";

void ResultsPopup::init()
{
    fh = new FileHandler();
}

void ResultsPopup::destroy()
{
    fh->removeOutputDir(resultsPath);
    QString rootfile = basePath + cRootDirExtension + "/" + trackerName + cRootFileExt;
    fh->removeTmpConfigFile(rootfile);
    delete fh;
}

void ResultsPopup::saveResults()
{
    QString savePath = QFileDialog::getExistingDirectory(QString::null, this, "saveres", "Save in...");
    if (!savePath.isEmpty()) {
	try {
	    QDir resDir(resultsPath);
	    resDir.setNameFilter(saveFilter);
	    resDir.setFilter(QDir::Files | QDir::NoSymLinks);
	    QDir saveDir(savePath);
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

void ResultsPopup::setBasePath(QString path)
{
    basePath = path;
}

QString& ResultsPopup::getResultsPath()
{
    return resultsPath;
}

void ResultsPopup::setResultsPath(QString path)
{
    resultsPath = path;
}

void ResultsPopup::setTrackerName(QString name)
{
    trackerName = name;
}

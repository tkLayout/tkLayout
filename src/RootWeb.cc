// Project includes
#include <RootWeb.hh>
#include "global_funcs.hh"

// standard includes
#include <limits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

// ROOT includes
#include <TCanvas.h>
#include <TError.h>
#include <time.h>
#include <TView.h>
#include <TFile.h>
#include <vector>
#include <TColor.h>
#include <TROOT.h>

// boost includes
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>



int RootWImage::imageCounter_ = 0;
std::map <std::string, int> RootWImage::imageNameCounter_;

using namespace RootWeb;

//*******************************************//
// RootWeb                                   //
//*******************************************//
std::string RootWeb::cleanUpObjectName(const std::string& source) {
  std::string dest = "";
  const std::string allowedChars  = "abcdefghijklmnopqrstuvwxyz"
                                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                    "0123456789";
  const std::string convertedChars = "._-() "; // including space here
  auto N=source.size();
  auto i=N;
  bool pleaseConvert = false;
  for (i=0; i<N; ++i) {
    if (allowedChars.find(source[i])<allowedChars.size()) {
      dest+=source[i];
      pleaseConvert = true;
    }
    else if (convertedChars.find(source[i])<allowedChars.size() && pleaseConvert) {
      dest+="_";
      // avoid converting many charachters in a row
      pleaseConvert = false;
    }
  }
  return dest;
}


//*******************************************//
// RootWTable                                //
//*******************************************//
RootWTable::RootWTable() {
  serialRow_ = 0;
  serialCol_ = 0;
}

ostream& RootWTable::dump(ostream& output) {
  RootWTable& myRootWTable = (*this);

  int minRow = tableContent_.begin()->first.first; int maxRow=0;
  int minCol = tableContent_.begin()->first.second; int maxCol=0;
  bool firstNotFound = true;

  pair<int, int> myIndex;
  rootWTableContent::iterator tableContentIt;
  rootWTableContent& myTableContent = myRootWTable.tableContent_;

  //std::cerr << "Table: "; // debug
  for (tableContentIt = tableContent_.begin();
       tableContentIt != tableContent_.end(); ++tableContentIt) {
    myIndex = (*tableContentIt).first;
    //std::cerr << "(" << myIndex.first << ", " << myIndex.second << ")" << std::endl; // debug
    if (firstNotFound) {
      firstNotFound = false;
      minRow = myIndex.first;
      maxRow = myIndex.first;
      minCol = myIndex.second;
      maxCol = myIndex.second;
    } else {
      minRow = ( myIndex.first < minRow ) ? myIndex.first : minRow;
      maxRow = ( myIndex.first > maxRow ) ? myIndex.first : maxRow;
      minCol = ( myIndex.second < minCol ) ? myIndex.second : minCol;
      maxCol = ( myIndex.second > maxCol ) ? myIndex.second : maxCol;
    }
  }
  //std::cerr << "Size:" << tableContent_.size() << "; Rows: " << minRow << "-" << maxRow <<"; Cols: " << minCol << "-" << maxCol << endl; // debug
  if (firstNotFound) return output;
  std::string myColorCode;
  std::string myCellCode;
  int myColorIndex;

  output << "<table>";
  for (int iRow = minRow; iRow<=maxRow; ++iRow) {
    output << "<tr>";
    for (int iCol = minCol; iCol<=maxCol; ++iCol) {
      if ((iRow==minRow)&&(iRow==0)) myCellCode = "th";
      else myCellCode = "td";
      output << "<" << myCellCode;
      myColorIndex = tableContentColor_[make_pair(iRow, iCol)];
      if (myColorIndex!=0) {
	      myColorCode=gROOT->GetColor(myColorIndex)->AsHexString();
	      output << " style=\"color:" << myColorCode << ";\" " << std::endl;
      }
      output << ">";
      output << myTableContent[make_pair(iRow, iCol)];
      output << "</" << myCellCode << ">" << " ";
    }
    output << "</tr>";
  }
  output << "</table>" << endl;

  return output;
}

void RootWTable::setColor(int row, int column, int newColor) {
  tableContentColor_[make_pair(row, column)] = newColor;
}

void RootWTable::setContent(int row, int column, string content) {
  // std::cerr << "setContent("<<row<<", "<<column<<", "<<content<<")"<<endl; // debug
  tableContent_[make_pair(row, column)] = content;
}

void RootWTable::setContent(int row, int column, int number) {
  // std::cerr << "setContent("<<row<<", "<<column<<", "<<number<<")"<<endl; // debug
  stringstream myNum_;
  myNum_.clear();
  myNum_ << dec << number;
  tableContent_[make_pair(row, column)] = myNum_.str();
}

void RootWTable::setContent(int row, int column, double number, int precision) {
  // std::cerr << "setContent("<<row<<", "<<column<<", "<<number<<")"<<endl; // debug
  stringstream myNum_;
  myNum_.clear();
  myNum_ << dec << fixed << setprecision(precision) << number;
  tableContent_[make_pair(row, column)] = myNum_.str();
}

pair<int, int> RootWTable::addContent(string myContent) {
  setContent(serialRow_, serialCol_++, myContent);
  return make_pair(serialRow_, serialCol_);
}

pair<int, int> RootWTable::addContent(int number) {
  stringstream myNum_;
  myNum_.clear();
  myNum_ << dec << number;
  return(addContent(myNum_.str()));
}

pair<int, int> RootWTable::addContent(double number, int precision) {
  stringstream myNum_;
  myNum_.clear();
  myNum_ << dec << fixed << setprecision(precision) << number;
  return(addContent(myNum_.str()));
}

pair<int, int> RootWTable::newLine() {
  serialRow_++;
  serialCol_=0;
  return make_pair(serialRow_, serialCol_);
}

//*******************************************//
// RootWImage                                //
//*******************************************//

RootWImage::RootWImage() {
  imageCounter_++;
  myCanvas_ = NULL;
  zoomedWidth_ = 0; zoomedHeight_ = 0;
  relativeHtmlDirectory_ = "";
  targetDirectory_ = "";
  comment_ = "";
  name_ = "img";
  allowedExtensions_ = DEFAULTALLOWEDEXTENSIONS;
  setDefaultExtensions();
}

RootWImage::RootWImage(TCanvas* myCanvas, int witdh, int height) {
  imageCounter_++;
  myCanvas_ = NULL;
  setCanvas(myCanvas);
  setZoomedSize(witdh, height);
  relativeHtmlDirectory_ = "";
  targetDirectory_ = "";
  comment_ = "";
  name_ = "img";
  allowedExtensions_ = DEFAULTALLOWEDEXTENSIONS;
  setDefaultExtensions();
}

RootWImage::RootWImage(TCanvas* myCanvas, int witdh, int height, string relativehtmlDirectory) {
  imageCounter_++;
  myCanvas_ = NULL;
  setCanvas(myCanvas);
  setZoomedSize(witdh, height);
  setRelativeHtmlDirectory(relativehtmlDirectory);
  targetDirectory_ = "";
  comment_ = "";
  name_ = "img";
  allowedExtensions_ = DEFAULTALLOWEDEXTENSIONS;
  setDefaultExtensions();
}

RootWImage::RootWImage(TCanvas& myCanvas, int witdh, int height) {
  imageCounter_++;
  myCanvas_ = NULL;
  setCanvas(myCanvas);
  setZoomedSize(witdh, height);
  relativeHtmlDirectory_ = "";
  targetDirectory_ = "";
  comment_ = "";
  name_ = "img";
  allowedExtensions_ = DEFAULTALLOWEDEXTENSIONS;
  setDefaultExtensions();
}

RootWImage::RootWImage(TCanvas& myCanvas, int witdh, int height, string relativehtmlDirectory) {
  imageCounter_++;
  myCanvas_ = NULL;
  setCanvas(myCanvas);
  setZoomedSize(witdh, height);
  setRelativeHtmlDirectory(relativehtmlDirectory);
  targetDirectory_ = "";
  comment_ = "";
  name_ = "img";
  allowedExtensions_ = DEFAULTALLOWEDEXTENSIONS;
  setDefaultExtensions();
}

RootWImage::~RootWImage() {
  if (myCanvas_) delete myCanvas_;
}

void RootWImage::setDefaultExtensions() {
   //addExtension("C");
   addExtension("pdf");
   addExtension("root");
}

void RootWImage::setComment(string newComment) {
  comment_ = newComment;
}

void RootWImage::setName(string newName) {
  name_ = newName;
}

std::string RootWImage::getName() {
  return name_ ;
}

void RootWImage::setCanvas(TCanvas* myCanvas) {
  if (myCanvas_) delete myCanvas_;
  myCanvas_ = (TCanvas*)myCanvas->DrawClone();
  std::ostringstream canvasName("");
  canvasName << "canvas" << setfill('0') << setw(3) << imageCounter_;
  myCanvas_->SetName(canvasName.str().c_str());
  TView* myView = myCanvas->GetView();
  if (myView) {
    TView* newView = myCanvas_->GetView();
    if (newView) {
      double min[3], max[3];
      Int_t irep;
      newView->SetView(myView->GetLongitude(),
		       myView->GetLatitude(),
		       myView->GetPsi(), irep);
      myView->GetRange(min, max);
      newView->SetRange(min, max);
    }
  }
}

void RootWImage::setCanvas(TCanvas& myCanvas) {
  setCanvas(&myCanvas);
}

void RootWImage::setZoomedSize(int witdh, int height) {
  zoomedWidth_ = witdh;
  zoomedHeight_ = height;
}

void RootWImage::setRelativeHtmlDirectory(string newDirectory) {
  relativeHtmlDirectory_ = newDirectory;
}

void RootWImage::setTargetDirectory(string newDirectory) {
  targetDirectory_ = newDirectory;
}

string RootWImage::saveFiles(int smallWidth, int smallHeight, int largeWidth, int largeHeight) {
  setZoomedSize(largeWidth, largeHeight);
  return saveFiles(smallWidth, smallHeight);
}

RootWImageSize RootWImage::makeSizeCode(int sw, int sh, int lw, int lh) {
  stringstream aString;
  aString << sw <<"x" << sh << "-" << lw <<"x" << lh;
  return aString.str();
}

string RootWImage::saveFiles(int smallWidth, int smallHeight) {
  RootWImageSize myImageSize;
  myImageSize = makeSizeCode(smallWidth,smallHeight,zoomedWidth_,zoomedHeight_);
  std::ostringstream tmpCanvasName("");
  tmpCanvasName << name_ << setfill('0') << setw(3) << imageNameCounter_[name_]++;
  string canvasName = tmpCanvasName.str();
  myCanvas_->SetName(canvasName.c_str());
  if (fileSaved_[myImageSize]) {
    return myText_[myImageSize];
  }

  stringstream thisText;
  thisText.clear();
  if (relativeHtmlDirectory_ == "") relativeHtmlDirectory_ = ".";
  if (targetDirectory_ == "" ) targetDirectory_ = ".";

  if (!myCanvas_) return "";

  //string canvasName = myCanvas_->GetName();
  string smallCanvasFileName = canvasName + "-" + myImageSize + "_small.png";
  string largeCanvasFileBaseName = canvasName;

  double wScale = double(smallWidth) / double(zoomedWidth_);
  double hScale = double(smallHeight) / double(zoomedHeight_);
  double myScale;

  myScale =  (hScale > wScale) ? hScale : wScale;
  int imgW = int (zoomedWidth_ * myScale);
  int imgH = int (zoomedHeight_ * myScale);
  int canW = int(imgW * thumb_compression_);
  int canH = int(imgH * thumb_compression_);

  string smallCanvasCompleteFileName = targetDirectory_+"/"+smallCanvasFileName;
  string largeCanvasCompleteFileName = targetDirectory_+"/"+largeCanvasFileBaseName + ".png";

  gErrorIgnoreLevel = 1500;
  myCanvas_->cd();
  myCanvas_->SetCanvasSize(canW, canH);
  // std::cerr << "Saving " << largeCanvasCompleteFileName << std::endl; // debug
  myCanvas_->Print(smallCanvasCompleteFileName.c_str());
  myCanvas_->SetCanvasSize(zoomedWidth_, zoomedHeight_);
  myCanvas_->Print(largeCanvasCompleteFileName.c_str());

  string fileTypeList;
  for (vector<string>::iterator it=fileTypeV_.begin(); it!=fileTypeV_.end(); ++it) {
    largeCanvasCompleteFileName = targetDirectory_+"/"+largeCanvasFileBaseName + "." + (*it);
    myCanvas_->Print(largeCanvasCompleteFileName.c_str());
    if (it!=fileTypeV_.begin()) fileTypeList+="|";
    fileTypeList+=(*it);
  }

  thisText << "<img class='wleftzoom' width='"<< imgW<<"' height='"<<imgH<<"' src='"<< relativeHtmlDirectory_ << "/" << smallCanvasFileName <<"' alt='"<< comment_ <<"' onclick=\"popupDiv('" << relativeHtmlDirectory_ << "/" << largeCanvasFileBaseName << "', "<< zoomedWidth_<<", "<< zoomedHeight_ << " , '"<< comment_ << "','"<<fileTypeList << "');\" />";

  myText_[myImageSize] = thisText.str();
  fileSaved_[myImageSize] = true;
  lastSize_ = myImageSize;
  return thisText.str();
}

ostream& RootWImage::dump(ostream& output) {
  output << myText_[lastSize_];
  return output;
}

bool RootWImage::addExtension(string myExtension) {
  // First check if the extension has not a "|" inside, which
  // would invalidate the search
  if (myExtension.find("|")!=string::npos) {
    cerr << "In bool RootWImage::addExtension(string myExtension) extension contains '|', which is not allowed" << endl;
    return false;
  }
  bool result = false;

  // Then create the real search pattern
  string mySearchExtension="|"+myExtension+"|";

  // Last: check if the extension is in the list of allowed extensions
  if (allowedExtensions_.find(mySearchExtension)!=string::npos) {
    // Try to find if the extension is already registered
    vector<string>::iterator anExtension = find(fileTypeV_.begin(), fileTypeV_.end(), myExtension);
    if (anExtension == fileTypeV_.end()) {
      fileTypeV_.push_back(myExtension);
      result=true;
    }
  } else {
    cerr << "Extension '" << myExtension << "' is not allowed for canvas saving: new file type not added." << endl;
  }

  return result;
}

void RootWImage::saveSummary(std::string baseName, TFile* myTargetFile) {
  if (!myCanvas_) return;
  baseName += name_;
  saveSummaryLoop(myCanvas_, baseName, myTargetFile);
}

void RootWImage::saveSummaryLoop(TPad* basePad, std::string baseName, TFile* myTargetFile) {
  TList* aList;
  TObject* anObject;
  TPad* myPad;
  std::string myClass;
  std::string myName;

  TNamed* aNamed;

  // TSystemFile* aFile;
  // string aFileName;
  // string aFileNameTail;
  // TFile* myRootFile;

  aList = basePad->GetListOfPrimitives();
  for (int i=0; i<aList->GetEntries(); ++i) {
    anObject = aList->At(i);
    myClass = anObject->ClassName();
    if (myClass=="TPad") { // Go one step inside
      myPad = (TPad*) anObject;
      saveSummaryLoop(myPad, baseName, myTargetFile);
    } else if (
	       (myClass=="TProfile") ||
	       (myClass=="TGraph") ||
	       (myClass=="TGraphErrors") ||
	       (myClass=="TH1I") ||
	       (myClass=="TH1F") ||
	       (myClass=="TH1D") ||
	       (myClass=="TH2C") ||
	       (myClass=="TH2I") ||
	       (myClass=="TH2F") ||
	       (myClass=="TH2D") ||
	       (myClass=="THStack")
	       ) {
      aNamed = (TNamed*) anObject;
      myTargetFile->cd();
      myName = Form("%s.%s", baseName.c_str(), aNamed->GetName());
      myName = RootWeb::cleanUpObjectName(myName);
      aNamed->SetName(myName.c_str());
      aNamed->Write();
    } else if (
	       (myClass=="TEllipse") ||
	       (myClass=="TFrame") ||
	       (myClass=="TLatex") ||
	       (myClass=="TLegend") ||
	       (myClass=="TLine") ||
	       (myClass=="TPaveText") ||
	       (myClass=="TPolyLine") ||
	       (myClass=="TText") ||
	       (myClass=="TArc") ||
	       (myClass=="TBox")
	       ) {
    } else {
      std::cerr << Form("Unhandled class %s", myClass.c_str()) << std::endl;
    }
  }
}

//*******************************************//
// RootWContent                              //
//*******************************************//

RootWContent::RootWContent() {
  title_ = "Untitled";
  targetDirectory_ = "";
  visible_ = true;
}

RootWContent::RootWContent(string title, bool visible /*=true*/) {
  title_ = title;
  targetDirectory_ = "";
  visible_ = visible;
}

void RootWContent::setTargetDirectory(string newTargetDirectory) {
  targetDirectory_ = newTargetDirectory;
}

RootWContent::~RootWContent() {
  RootWItem* myItem;

  while(!itemList_.empty()) {
    myItem = (*itemList_.begin());
    if (myItem!=NULL) {
      //delete myItem;
      myItem = NULL;
    }
    itemList_.erase(itemList_.begin());
  }
}

void RootWContent::addItem(RootWItem* newItem) {
  if (newItem) {
    itemList_.push_back(newItem);
  }
}

void RootWContent::setTitle(string newTitle) {
  title_ = newTitle;
}

void RootWContent::addParagraph(string parText) {
  RootWText* newText = new RootWText;
  //newText->addText("<p>");
  newText->addText(parText);
  newText->addText("<br/>\n");
  //newText->addText("</p>\n");


  addItem(newText);
}

RootWText& RootWContent::addText() {
  RootWText* newText = new RootWText();
  addItem(newText);
  return (*newText);
}

RootWText& RootWContent::addText(string newText) {
  RootWText* newTextItem = new RootWText(newText);
  addItem(newTextItem);
  return (*newTextItem);
}

RootWInfo& RootWContent::addInfo() {
  RootWInfo* newInfo = new RootWInfo();
  addItem(newInfo);
  return (*newInfo);
}

RootWInfo& RootWContent::addInfo(string description)  {
  RootWInfo* newInfo = new RootWInfo(description);
  addItem(newInfo);
  return (*newInfo);
}

RootWInfo& RootWContent::addInfo(string description, string value)  {
  RootWInfo* newInfo = new RootWInfo(description, value);
  addItem(newInfo);
  return (*newInfo);
}

RootWTable& RootWContent::addTable() {
  RootWTable* newTable = new RootWTable();
  addItem(newTable);
  return (*newTable);
}

RootWImage& RootWContent::addImage() {
  RootWImage* newImage = new RootWImage();
  addItem(newImage);
  return (*newImage);
}

RootWImage& RootWContent::addImage(TCanvas* myCanvas, int witdh, int height) {
  RootWImage* newImage = new RootWImage(myCanvas, witdh, height);
  addItem(newImage);
  return (*newImage);
}

RootWImage& RootWContent::addImage(TCanvas* myCanvas, int witdh, int height, string relativeHtmlDirectory) {
  RootWImage* newImage = new RootWImage(myCanvas, witdh, height, relativeHtmlDirectory);
  addItem(newImage);
  return (*newImage);
}

RootWImage& RootWContent::addImage(TCanvas& myCanvas, int witdh, int height) {
  RootWImage* newImage = new RootWImage(myCanvas, witdh, height);
  addItem(newImage);
  return (*newImage);
}

RootWImage& RootWContent::addImage(TCanvas& myCanvas, int witdh, int height, string relativeHtmlDirectory) {
  RootWImage* newImage = new RootWImage(myCanvas, witdh, height, relativeHtmlDirectory);
  addItem(newImage);
  return (*newImage);
}

RootWTextFile& RootWContent::addTextFile() {
  RootWTextFile* newTextFile = new RootWTextFile();
  addItem(newTextFile);
  return (*newTextFile);
}

RootWTextFile& RootWContent::addTextFile(string newFileName) {
  RootWTextFile* newTextFile = new RootWTextFile(newFileName);
  addItem(newTextFile);
  return (*newTextFile);
}

RootWTextFile& RootWContent::addTextFile(string newFileName, string newDescription) {
  RootWTextFile* newTextFile = new RootWTextFile(newFileName, newDescription);
  addItem(newTextFile);
  return (*newTextFile);
}

RootWBinaryFile& RootWContent::addBinaryFile() {
  RootWBinaryFile* newBinaryFile = new RootWBinaryFile();
  addItem(newBinaryFile);
  return (*newBinaryFile);
}

RootWBinaryFile& RootWContent::addBinaryFile(string newFileName) {
  RootWBinaryFile* newBinaryFile = new RootWBinaryFile(newFileName);
  addItem(newBinaryFile);
  return (*newBinaryFile);
}

RootWBinaryFile& RootWContent::addBinaryFile(string newFileName, string newDescription) {
  RootWBinaryFile* newBinaryFile = new RootWBinaryFile(newFileName, newDescription);
  addItem(newBinaryFile);
  return (*newBinaryFile);
}

RootWBinaryFile& RootWContent::addBinaryFile(string newFileName, string newDescription, string newOriginalFile) {
  RootWBinaryFile* newBinaryFile = new RootWBinaryFile(newFileName, newDescription, newOriginalFile);
  addItem(newBinaryFile);
  return (*newBinaryFile);
}


ostream& RootWContent::dump(ostream& output) {
  vector<RootWItem*>::iterator it;
  RootWItem* myItem;
  RootWImage* myImage;
  RootWFile* myFile;
  TFile* summaryFile = getSummaryFile();
  std::string baseName="";
  if (page_ != nullptr) baseName = page_->getTitle() + ".";

  //std::cerr << "Content: " << title_ <<endl; //debug
  output << "<h2 class=\"hidingTitle\">"<<title_<<"</h2>" << endl;
  if (visible_) {
    output << "<div class=\"hideable\"> ";
  } else {
    output << "<div class=\"hidden\"> ";
  }
  for (it = itemList_.begin(); it != itemList_.end(); ++it) {
    myItem =(*it);
    if (myItem->isImage()) {
      if ( (myImage=dynamic_cast<RootWImage*>(myItem)) ) {
        myImage->setTargetDirectory(targetDirectory_);
	if (summaryFile) myImage->saveSummary(baseName, summaryFile);
        myImage->saveFiles(THUMBSMALLSIZE, THUMBSMALLSIZE);
      } else {
        cout << "WARNING: this should never happen. contact the author immediately!" << endl;
      }
    }
    if (myItem->isFile()) {
      if ( (myFile=dynamic_cast<RootWFile*>(myItem))) { // CUIDADO This is not good OOP for God's sake!!!
        myFile->setTargetDirectory(targetDirectory_);
      } else if (RootWFileList* myFileList=dynamic_cast<RootWFileList*>(myItem)) {
        myFileList->setTargetDirectory(targetDirectory_);
      } else {
        cout << "WARNING: this should never happen. contact the author immediately!" << endl;
      }
    }
    myItem->dump(output);
  }
  output << "<div class=\"clearer\">&nbsp;</div>";
  output << "</div>" << endl;
  return output;
}

void RootWContent::setPage(RootWPage* newPage) {
  page_ = newPage;
}

TFile* RootWContent::getSummaryFile() {
  if (page_==nullptr) {
    return nullptr;
  } else {
    return page_->getSummaryFile();
  }
}

//*******************************************//
// RootWPage                                 //
//*******************************************//

RootWPage::RootWPage() {
  address_ = "";
  title_ = "Untitled";
  site_ = NULL;
  targetDirectory_ = "";
  alert_ = 0;
}

RootWPage::RootWPage(string title) {
  setTitle(title);
  site_ = NULL;
  targetDirectory_ = "";
  alert_ = 0;
}

RootWPage::~RootWPage() {
}

void RootWPage::setTargetDirectory(string newTargetDirectory) {
  targetDirectory_ = newTargetDirectory;
}

void RootWPage::setTitle(string newTitle) {
  title_ = newTitle;
  if (address_=="") address_ = title_+".html";
}

string RootWPage::getTitle() {
  return title_;
}

void RootWPage::setAddress(string newAddress) {
  address_ = newAddress;
}

string RootWPage::getAddress() {
  return address_;
}

void RootWPage::setSite(RootWSite* newSite) {
  site_ = newSite;
}

void RootWPage::addContent(RootWContent* newContent) {
  if (newContent) {
    contentList_.push_back(newContent);
  }
}

RootWContent& RootWPage::addContent(string title, bool visible /*=true*/) {
  RootWContent* newContent = new RootWContent(title, visible);
  addContent(newContent);
  return (*newContent);
}

void RootWPage::setAlert(double alert) {
  if (alert<0) alert_=0;
  else if (alert>1) alert_=1;
  else alert_=alert;
}

double RootWPage::getAlert() {
  return alert_;
}

ostream& RootWPage::dump(ostream& output) {
  bool fakeSite=false;
  if (site_==NULL) {
    fakeSite = true;
    site_ = new RootWSite();
    site_->addPage(this);
  }

  // Header standard part, with the menu list, etc...
  // up to opening the primaryContent div
  site_->dumpHeader(output, this);

  vector<RootWContent*>::iterator it;
  RootWContent* myContent;
  for (it = contentList_.begin(); it != contentList_.end(); ++it) {
    myContent =(*it);
    myContent->setTargetDirectory(targetDirectory_);
    myContent->setPage(this);
    myContent->dump(output);
  }

  // Header standard part, with copyright, etc...
  // starting from closing the primaryContent div
  site_->dumpFooter(output);

  if (fakeSite) {
    delete site_;
    site_=NULL;
  }

  return output;

}


void RootWPage::setRelevance(int newRelevance) {
  relevance = newRelevance;
}

int RootWPage::getRelevance() {
  return relevance;
}

TFile* RootWPage::getSummaryFile() {
  if (site_) {
    return site_->getSummaryFile();
  } else {
    return nullptr;
  }
}


//*******************************************//
// RootWSite                                 //
//*******************************************//

RootWSite::RootWSite() {
  title_ = "Untitled";
  comment_ = "";
  commentLink_ = "";
  programName_ = DEFAULTPROGRAMNAME;
  programSite_ = DEFAULTPROGRAMSITE;
  revision_="";
  targetDirectory_ = ".";
  summaryFile_ = nullptr;
  createSummaryFile_ = true;
  summaryFileName_ = "summary.root";
}

RootWSite::RootWSite(string title) {
  title_ = title;
  comment_ = "";
  commentLink_ = "";
  programName_ = DEFAULTPROGRAMNAME;
  programSite_ = DEFAULTPROGRAMSITE;
  revision_="";
  targetDirectory_ = ".";
  summaryFile_ = nullptr;
  createSummaryFile_ = true;
  summaryFileName_ = "summary.root";
}

RootWSite::RootWSite(string title, string comment) {
  title_ = title;
  comment_ = comment;
  commentLink_ = "";
  programName_ = DEFAULTPROGRAMNAME;
  programSite_ = DEFAULTPROGRAMSITE;
  revision_="";
  targetDirectory_ = ".";
  summaryFile_ = nullptr;
  createSummaryFile_ = true;
  summaryFileName_ = "summary.root";
}

RootWSite::~RootWSite() {
}

void RootWSite::setTitle(string newTitle) {
  title_ = newTitle;
}

void RootWSite::setComment(string newComment) {
  comment_ = newComment;
}

void RootWSite::setCommentLink(string newCommentLink) {
  commentLink_ = newCommentLink;
}

string RootWSite::getTitle() {
  return title_;
}

string RootWSite::getComment() {
  return comment_;
}

string RootWSite::getCommentLink() {
  return commentLink_;
}

string RootWSite::getRevision() {
  return revision_;
}

void RootWSite::setRevision(string newRevision) {
  revision_ = newRevision;
}

void RootWSite::addPage(RootWPage* newPage, int relevance /* = least_relevant */) {
  newPage->setRelevance(relevance);
  if (relevance == least_relevant) {
    pageList_.push_back(newPage);
  } else if (relevance == most_relevant) {
    pageList_.insert(pageList_.begin(), newPage);
  } else {
    // Go through the vector to find the first element which has a
    // lower relevance than 'relevance'
    vector<RootWPage*>::iterator it = pageList_.begin();
    vector<RootWPage*>::iterator beforeThis = pageList_.end();
    RootWPage* aPage;
    for (; it!=pageList_.end(); ++it) {
      aPage = *it;
      if (aPage->getRelevance()<relevance) {
        beforeThis = it;
        break;
      }
    }
    pageList_.insert(beforeThis, newPage);
  }
  newPage->setSite(this);
}

RootWPage& RootWSite::addPage(string newTitle, int relevance /* = least_relevant */ ) {
  RootWPage* newPage = new RootWPage(newTitle);
  addPage(newPage, relevance);
  return (*newPage);
}

void RootWSite::addAuthor(string newAuthor) {
  authorList_.push_back(newAuthor);
}

ostream& RootWSite::dumpFooter(ostream& output) {
  output  << "          </div>" << endl
	  << "        </div>" << endl
	  << "        <div class=\"clear\"></div>" << endl
	  << "      </div>" << endl
	  << "      <div id=\"footer\">" << endl;

  // Add the list of authors if any
  bool firstAuthorFound=false;
  for (vector<string>::iterator it= authorList_.begin();
       it != authorList_.end(); ++it) {
    if (!firstAuthorFound) {
      firstAuthorFound = true;
      output << "        <p>&copy; ";
      output << (*it);
    } else {
      output << ", " << (*it);
    }
  }
  if (firstAuthorFound) output <<"</p>" << endl;

  time_t rawtime;
  time ( &rawtime );
  output << "        <p>Page created on "<< asctime(gmtime ( &rawtime )) << " GMT</p>" << endl
	 << "        <p>by <a href=\""<< programSite_ <<"\">"<<programName_<<"</a>";
  if (revision_!="")
    output << " revision " << revision_;
  output << "</p>" << endl
	 << "      </div>" << endl
	 << "    </div>" << endl
	 << "  </body>" << endl
	 << "</html>" << endl;
  return output;
}

ostream& RootWSite::dumpHeader(ostream& output, RootWPage* thisPage) {
  output << "<html xmlns=\"http://www.w3.org/1999/xhtml\">" << endl
	 << "  <head>" << endl
	 << "    <meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />" << endl
	 << "    <title>"<< title_ << " - " << thisPage->getTitle() <<"</title>" << endl
         << "    <meta http-equiv=\"cache-control\" content=\"max-age=0\" />" << endl
         << "    <meta http-equiv=\"cache-control\" content=\"no-cache\" />" << endl
         << "    <meta http-equiv=\"expires\" content=\"0\" />" << endl
         << "    <meta http-equiv=\"expires\" content=\"Tue, 01 Jan 1980 1:00:00 GMT\" />" << endl
         << "    <meta http-equiv=\"pragma\" content=\"no-cache\" />" << endl
	 << "    <meta name=\"keywords\" content=\"CERN CMS tracker upgrade\" />" << endl
	 << "    <meta name=\"description\" content=\"CMS Tracker upgrade summary page\" />" << endl
	 << "    <link href=\"../style/default.css\" rel=\"stylesheet\" type=\"text/css\" />" << endl
	 << "    <link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"../style/images/favicon.ico\">" << endl
	 << " </head>" << endl
	 << "  <script type=\"text/javascript\" src=\"../style/functions.js\"></script>" << endl
	 << "  <body onload=\"preparePageGoodies();\">" << endl
	 << "    <div id=\"outer\">" << endl
	 << "      <div id=\"header\">" << endl
	 << "        <h1>" << endl
	 << "    <a href=\"index.html\">" << title_ << "</a>" << endl
	 << "        </h1>" << endl;
  if (commentLink_=="") {
    output << "    <h2>"<<comment_<<"</h2>" << endl;
  } else {
    output << "    <h2>"
	   << "<a href=\"" << commentLink_ << "\">"
	   << comment_
	   << "</a>"
	   << "</h2>" << endl;
  }
  output << "      </div>" << endl
	 << "      <div id=\"menu\">" << endl
	 << "        <ul>" << endl;
  for (vector<RootWPage*>::iterator it = pageList_.begin();
       it !=pageList_.end(); ++it) {
    output << "          <li";
    if ((*it)==(thisPage)) output << " class=\"current\"";
    output << ">";
    output  << "    <a href=\"" << (*it)->getAddress() << "\" ";
    if ((*it)->getAlert()>0) {
      ostringstream anAlarmColor;
      anAlarmColor.str("");
      anAlarmColor << "rgb(255, "
		   << int(255 * (*it)->getAlert())
		   << ", 0)";
      output << " style = 'background:"
		   << anAlarmColor.str()
		   << ";' ";
    }
    output << ">";
    output << (*it)->getTitle();
    output <<"</a>" << "          </li>" << endl;
  }
  output << "        </ul>" << endl
	 << "      </div>" << endl
	 << "      <div id=\"content\">" << endl
	 << "        <div id=\"primaryContentContainer\">" << endl
	 << "          <div id=\"primaryContent\">" << endl;

  return output;

}

bool RootWSite::makeSite(bool verbose) {
  ofstream myPageFile;
  RootWPage* myPage;
  string myPageFileName;
  //string targetStyleDirectory = targetDirectory_ + "/style";


  // Check if the directory already exists
  if (boost::filesystem::exists( targetDirectory_ )) {
    if (! boost::filesystem::is_directory(targetDirectory_) ) {
      cerr << "A file named " << targetDirectory_ << " already exists. Cannot create target directory" << endl;
      return false;
    }
  } else {
    // If the directory does not exist, then try to create it
    if (!boost::filesystem::create_directory(targetDirectory_)) {
      cerr << "Couldn't create directory " << targetDirectory_ << endl;
      return false;
    }
  }

  vector<RootWPage*>::iterator it;
  if (createSummaryFile_) {
    summaryFile_ = new TFile(Form("%s/%s",
				  targetDirectory_.c_str(),
				  summaryFileName_.c_str()), "RECREATE");
  } else summaryFile_ = nullptr;
  for (it=pageList_.begin(); it!=pageList_.end(); it++) {
    myPage = (*it);
    if (verbose) std::cout << " " << myPage->getTitle() << std::flush;
    myPageFileName = targetDirectory_+"/"+myPage->getAddress();
    myPageFile.open(myPageFileName.c_str(), ios::out);
    myPage->setTargetDirectory(targetDirectory_);
    myPage->dump(myPageFile);
    myPageFile.close();
    if (it==pageList_.begin()) {
      namespace fs = boost::filesystem;
      fs::path const from_path = myPageFileName;
      fs::path const to_path = targetDirectory_+"/index.html";
      fs::copy_file(from_path, to_path, fs::copy_option::overwrite_if_exists);
    }
  }
  if (verbose) std::cout << " ";
  if (summaryFile_) summaryFile_->Close();

  return true;
}

TFile* RootWSite::getSummaryFile() {
  return summaryFile_;
}

void RootWSite::setSummaryFile(bool doSummary) {
  createSummaryFile_ = doSummary;
}

void RootWSite::setSummaryFileName(std::string newName) {
  summaryFileName_ = newName;
  createSummaryFile_ = true;
}


//*******************************************//
// RootWItemCollection                       //
//*******************************************//


RootWItem* RootWItemCollection::getItem(string itemName) {
  RootWItem* result;
  result=itemCollection_[itemName];
  if (result!=NULL) result->taken=true;
  return result;
}

void RootWItemCollection::addItem(RootWItem* anItem, string itemName) {
  RootWItem* findItem;
  findItem=itemCollection_[itemName];
  if (findItem!=NULL) {
    cerr << "Warning in RootWItemCollection::addItem(): item named " << itemName << " already exists! Overwriting it." << endl;
  }
  itemCollection_[itemName]=anItem;
}

vector<RootWItem*> RootWItemCollection::getOtherItems() {
  vector<RootWItem*> result;
  map<string, RootWItem*>::iterator itItem=itemCollection_.begin();
  map<string, RootWItem*>::iterator endItem=itemCollection_.end();
  for (; itItem!=endItem; ++itItem) {
    if (!(*itItem).second->taken) result.push_back((*itItem).second);
  }
  return result;
}

//*******************************************//
// RootWTextFile                             //
//*******************************************//

ostream& RootWTextFile::dump(ostream& output) {

  if (fileName_=="") {
    cerr << "Warning: RootWTextFile::dump() was called without prior setting the destination file name" << endl;
    return output;
  }

  std::ofstream outputFile;
  string destinationFileName = targetDirectory_ +"/" + fileName_;
  outputFile.open(destinationFileName.c_str());
  // TODO: add a check here if the file was correctly opened
  // outputFile << myText_.str();
  // myText_ >> outputFile;
  outputFile << myText_.str() << endl;
  outputFile.close();

  output << "<b>" << description_ << ":</b> <a href=\""
         << fileName_ << "\">"
         << fileName_ << "</a></tt><br/>";
  return output;
};


//*******************************************//
// RootWBinaryFile                           //
//*******************************************//

ostream& RootWBinaryFile::dump(ostream& output) {
  if ((originalFileName_=="")&&(!noCopy_)) {
    cerr << "Warning: RootWBinaryFile::dump() was called without prior setting the original file name" << endl;
    return output;
  }
  if (fileName_=="") {
    cerr << "Warning: RootWBinaryFile::dump() was called without prior setting the destination file name" << endl;
    return output;
  }

  string destinationFileName = targetDirectory_ +"/" + fileName_;

  if ((!noCopy_)&&(boost::filesystem::exists(originalFileName_) && originalFileName_ != destinationFileName)) { // CUIDADO: naive control on copy on itself. it only matches the strings, not taking into account relative paths and symlinks
    try {
      if (boost::filesystem::exists(destinationFileName))
        boost::filesystem::remove(destinationFileName);
      boost::filesystem::copy_file(originalFileName_, destinationFileName);
    } catch (boost::filesystem::filesystem_error e) {
      cerr << e.what() << endl;
      return output;
    }
  }

  output << "<b>" << description_ << ":</b> <a href=\""
         << fileName_ << "\">"
         << fileName_ << "</a></tt><br/>";
  return output;
};

//*******************************************//
// RootWBinaryFileList                       //
//*******************************************//

ostream& RootWBinaryFileList::dump(ostream& output) {
  if (originalFileNames_.empty()) {
    cerr << "Warning: RootWBinaryFileList::dump() was called without prior setting the original file names" << endl;
    return output;
  }
  if (fileNames_.empty()) {
    cerr << "Warning: RootWBinaryFileList::dump() was called without prior setting the destination file names" << endl;
    return output;
  }
  if (fileNames_.size() != originalFileNames_.size()) {
    cerr << "Warning: RootWBinaryFileList::dump() was called with original and destination file name lists differing in length" << endl;
    return output;
  }

  auto it = originalFileNames_.begin();
  for (auto fn : fileNames_) {
    //auto path = split(fn, "/");
    string destinationFileName = targetDirectory_;
    //std::for_each(path.begin(), path.end()-1, [&destinationFileName](string p) {
    //  destinationFileName += "/" + p;
    //  if (!boost::filesystem::exists(destinationFileName)) boost::filesystem::create_directory(destinationFileName);
    //});
    destinationFileName += "/" + fn; //path.back();
    try {
      if (boost::filesystem::exists(destinationFileName)) boost::filesystem::remove(destinationFileName);
      boost::filesystem::copy_file(*it++, destinationFileName);
    }
    catch (boost::filesystem::filesystem_error e) {
      cerr << e.what() << endl;
      return output;
    }
  }

  std::vector<std::string> cleanedUpFileNames;
  std::transform(originalFileNames_.begin(), originalFileNames_.end(), std::back_inserter(cleanedUpFileNames), [](const std::string& s) {
      auto pos = s.find("stdinclude");
      if (s.find("xml") != string::npos) pos = s.rfind("/") + 1;
      return pos != string::npos ? s.substr(pos) : s;
    });
  output << "<b>" << description_ << ":</b>";
  auto dfn = fileNames_.begin();
  for (auto cfn : cleanedUpFileNames) {
    output << (cleanedUpFileNames.size() > 1 ? "<br>&nbsp&nbsp&nbsp&nbsp" : "") << " <a href=\"" << *(dfn++) << "\">" << cfn << "</a></tt>" << " ";
  }
  output << "<br/>";
  return output;
}


//*******************************************//
// RootWInfo                                 //
//*******************************************//

string RootWInfo::setValue(string newText) {
  value_ = newText;
  return newText;
}

string RootWInfo::setValue(int number) {
  stringstream myNum_;
  myNum_.clear();
  myNum_ << dec << number;
  return(setValue(myNum_.str()));
}

string RootWInfo::setValue(double number, int precision) {
  stringstream myNum_;
  myNum_.clear();
  myNum_ << dec << fixed << setprecision(precision) << number;
  return(setValue(myNum_.str()));
}

string RootWInfo::setValueSci(double number, int precision) {
  stringstream myNum_;
  myNum_.clear();
  int nearestLog = floor(log10(number));
  double mantissa = number/pow(10, nearestLog);
  myNum_.unsetf(ios_base::floatfield);
  myNum_.precision(5);
  myNum_ << dec << mantissa;
  myNum_ << "&times;10<sup>" << nearestLog << "</sup>";
  return(setValue(myNum_.str()));
}

string RootWInfo::appendValue(string aString) {
  stringstream myNum_;
  return setValue(value_ + aString);
}

ostream& RootWInfo::dump(ostream& output) {
  std::ofstream outputFile;

  output << "<b>" << description_ << ":</b> "
         << value_ << "</a></tt><br/>";
  return output;
}

//*******************************************//
// RootWGraphViz                             //
//*******************************************//

ostream& RootWGraphViz::dump(ostream& output) {
  string myDescription = description_;
  description_ += " (GraphViz)";
  RootWTextFile::dump(output);
  int result = system("which dot > /dev/null");
  if (result==0) {
    string svgFileName = fileName_+".svg";
    string command = "dot -Tsvg " + targetDirectory_ +"/" + fileName_ + " > " + targetDirectory_ + "/" + svgFileName;
    result = system(command.c_str());
    if (result==0) {
      RootWBinaryFile* mySvg = new RootWBinaryFile(svgFileName, myDescription+" (SVG)" );
      mySvg->setNoCopy(true);
      mySvg->dump(output);
    }
  }
  return output;
}

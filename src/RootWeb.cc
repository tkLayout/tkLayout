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
RootWTable::RootWTable(bool isPlacedBelow /*= false*/) {
  isPlacedBelow_ = isPlacedBelow;
  serialRow_ = 0;
  serialCol_ = 0;
  maxRow_ = 0;
  maxCol_ = 0;
}

std::ostream& RootWTable::dump(std::ostream& output) {
  RootWTable& myRootWTable = (*this);

  int minRow = tableContent_.begin()->first.first; int maxRow=0;
  int minCol = tableContent_.begin()->first.second; int maxCol=0;
  bool firstNotFound = true;

  std::pair<int, int> myIndex;
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

      if ((iRow==minRow) && (iRow==0)) myCellCode = "th";
      else myCellCode = "td";

      output << "<" << myCellCode;
      myColorIndex = tableContentColor_[std::make_pair(iRow, iCol)];
      if (myColorIndex!=0) {
	      myColorCode=gROOT->GetColor(myColorIndex)->AsHexString();
	      output << " style=\"color:" << myColorCode << ";\" " << std::endl;
      }
      output << ">";

      // FINT OUT WHETHER THE CELL SHOULD BE BOLD
      const std::pair<int, int> myCellPos = std::make_pair(iRow, iCol);
      const auto& isBoldIt = tableContentBold_.find(myCellPos);
      const bool isBoldCell = (isBoldIt != tableContentBold_.end() ? isBoldIt->second : false);
      if (isBoldCell) output << "<strong>";

      // ADD CONTENT
      output << myTableContent[std::make_pair(iRow, iCol)];

      // BOLD CELL CASE
      if (isBoldCell) output << "</strong>";

      output << "</" << myCellCode << ">" << " ";
    }
    output << "</tr>";
  }
  output << "</table>" << std::endl;

  return output;
}

void RootWTable::setColor(const int row, const int column, const int color) {
  if (color != kBlack) {
    tableContentColor_[std::make_pair(row, column)] = color;
  }
}

void RootWTable::setBold(const int row, const int column, const bool isBold) {
  tableContentBold_[std::make_pair(row, column)] = isBold;
}

void RootWTable::setContent(int row, int column, std::string content, const bool isBold, const int color) {
  tableContent_[std::make_pair(row, column)] = content;
  maxRow_ = MAX(maxRow_, row);
  maxCol_ = MAX(maxCol_, column);
  setBold(row, column, isBold);
  setColor(row, column, color);
}

void RootWTable::setContent(int row, int column, int number, const bool isBold, const int color) {
  std::stringstream myNum_;
  myNum_.clear();
  myNum_ << std::dec << number;
  tableContent_[std::make_pair(row, column)] = myNum_.str();
  maxRow_ = MAX(maxRow_, row);
  maxCol_ = MAX(maxCol_, column);
  setBold(row, column, isBold);
  setColor(row, column, color);
}

void RootWTable::setContent(int row, int column, double number, int precision, const bool isBold, const int color) {
  std::stringstream myNum_;
  myNum_.clear();
  myNum_ << std::dec << std::fixed << std::setprecision(precision) << number;
  tableContent_[std::make_pair(row, column)] = myNum_.str();
  maxRow_ = MAX(maxRow_, row);
  maxCol_ = MAX(maxCol_, column);
  setBold(row, column, isBold);
  setColor(row, column, color);
}

std::pair<int, int> RootWTable::addContent(std::string myContent) {
  setContent(serialRow_, serialCol_++, myContent);
  return std::make_pair(serialRow_, serialCol_);
}

std::pair<int, int> RootWTable::addContent(int number) {
  std::stringstream myNum_;
  myNum_.clear();
  myNum_ << std::dec << number;
  return(addContent(myNum_.str()));
}

std::pair<int, int> RootWTable::addContent(double number, int precision) {
  std::stringstream myNum_;
  myNum_.clear();
  myNum_ << std::dec << std::fixed << std::setprecision(precision) << number;
  return(addContent(myNum_.str()));
}

std::pair<int, int> RootWTable::newLine() {
  serialRow_++;
  serialCol_=0;
  return std::make_pair(serialRow_, serialCol_);
}

//*******************************************//
// RootWImage                                //
//*******************************************//

RootWImage::RootWImage() {
  imageCounter_++;
  myCanvas_ = nullptr;
  zoomedWidth_ = 0; zoomedHeight_ = 0;
  relativeHtmlDirectory_ = "";
  targetDirectory_ = "";
  comment_ = "";
  name_ = "img";
  allowedExtensions_ = DEFAULTALLOWEDEXTENSIONS;
  setDefaultExtensions();
}

RootWImage::RootWImage(std::unique_ptr<TCanvas> myCanvas, int witdh, int height) {
  imageCounter_++;
  myCanvas_ = nullptr;
  setCanvas(std::move(myCanvas));
  setZoomedSize(witdh, height);
  relativeHtmlDirectory_ = "";
  targetDirectory_ = "";
  comment_ = "";
  name_ = "img";
  allowedExtensions_ = DEFAULTALLOWEDEXTENSIONS;
  setDefaultExtensions();
}

RootWImage::RootWImage(std::unique_ptr<TCanvas> myCanvas, int witdh, int height, std::string relativehtmlDirectory) {
  imageCounter_++;
  myCanvas_ = nullptr;
  setCanvas(std::move(myCanvas));
  setZoomedSize(witdh, height);
  setRelativeHtmlDirectory(relativehtmlDirectory);
  targetDirectory_ = "";
  comment_ = "";
  name_ = "img";
  allowedExtensions_ = DEFAULTALLOWEDEXTENSIONS;
  setDefaultExtensions();
}

void RootWImage::setDefaultExtensions() {
   //addExtension("C");
   addExtension("pdf");
   addExtension("root");
}

void RootWImage::setComment(std::string newComment) {
  comment_ = newComment;
}

void RootWImage::setName(std::string newName) {
  name_ = newName;
}

std::string RootWImage::getName() {
  return name_ ;
}

void RootWImage::setCanvas(std::unique_ptr<TCanvas> myCanvas) {
  myCanvas_.reset(myCanvas.release()); // KEY POINT: instead of using TCanvas::DrawClone(), just transfer TCanvas object ownership from Vizard to RootWeb!
  std::ostringstream canvasName("");
  canvasName << "canvas" << std::setfill('0') << std::setw(3) << imageCounter_;
  myCanvas_->SetName(canvasName.str().c_str());
}

void RootWImage::setZoomedSize(int witdh, int height) {
  zoomedWidth_ = witdh;
  zoomedHeight_ = height;
}

void RootWImage::setRelativeHtmlDirectory(std::string newDirectory) {
  relativeHtmlDirectory_ = newDirectory;
}

void RootWImage::setTargetDirectory(std::string newDirectory) {
  targetDirectory_ = newDirectory;
}

std::string RootWImage::saveFiles(int smallWidth, int smallHeight, int largeWidth, int largeHeight) {
  setZoomedSize(largeWidth, largeHeight);
  return saveFiles(smallWidth, smallHeight);
}

RootWImageSize RootWImage::makeSizeCode(int sw, int sh, int lw, int lh) {
  std::stringstream aString;
  aString << sw <<"x" << sh << "-" << lw <<"x" << lh;
  return aString.str();
}

std::string RootWImage::saveFiles(int smallWidth, int smallHeight) {
  RootWImageSize myImageSize;
  myImageSize = makeSizeCode(smallWidth,smallHeight,zoomedWidth_,zoomedHeight_);
  std::ostringstream tmpCanvasName("");
  tmpCanvasName << name_ << std::setfill('0') << std::setw(3) << imageNameCounter_[name_]++;
  std::string canvasName = tmpCanvasName.str();
  myCanvas_->SetName(canvasName.c_str());
  if (fileSaved_[myImageSize]) {
    return myText_[myImageSize];
  }

  std::stringstream thisText;
  thisText.clear();
  if (relativeHtmlDirectory_ == "") relativeHtmlDirectory_ = ".";
  if (targetDirectory_ == "" ) targetDirectory_ = ".";

  if (!myCanvas_) return "";

  //string canvasName = myCanvas_->GetName();
  std::string smallCanvasFileName = canvasName + "-" + myImageSize + "_small.png";
  std::string largeCanvasFileBaseName = canvasName;

  double wScale = double(smallWidth) / double(zoomedWidth_);
  double hScale = double(smallHeight) / double(zoomedHeight_);
  double myScale;

  myScale =  (hScale > wScale) ? hScale : wScale;
  int imgW = int (zoomedWidth_ * myScale);
  int imgH = int (zoomedHeight_ * myScale);
  int canW = int(imgW * thumb_compression_);
  int canH = int(imgH * thumb_compression_);

  std::string smallCanvasCompleteFileName = targetDirectory_+"/"+smallCanvasFileName;
  std::string largeCanvasCompleteFileName = targetDirectory_+"/"+largeCanvasFileBaseName + ".png";

  gErrorIgnoreLevel = 1500;
  myCanvas_->cd();
  myCanvas_->SetCanvasSize(canW, canH);
  // std::cerr << "Saving " << largeCanvasCompleteFileName << std::endl; // debug
  myCanvas_->Print(smallCanvasCompleteFileName.c_str());
  myCanvas_->SetCanvasSize(zoomedWidth_, zoomedHeight_);
  myCanvas_->Print(largeCanvasCompleteFileName.c_str());

  std::string fileTypeList;
  for (std::vector<std::string>::iterator it=fileTypeV_.begin(); it!=fileTypeV_.end(); ++it) {
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

std::ostream& RootWImage::dump(std::ostream& output) {
  output << myText_[lastSize_];
  return output;
}

bool RootWImage::addExtension(std::string myExtension) {
  // First check if the extension has not a "|" inside, which
  // would invalidate the search
  if (myExtension.find("|")!=std::string::npos) {
    std::cerr << "In bool RootWImage::addExtension(std::string myExtension) extension contains '|', which is not allowed" << std::endl;
    return false;
  }
  bool result = false;

  // Then create the real search pattern
  std::string mySearchExtension="|"+myExtension+"|";

  // Last: check if the extension is in the list of allowed extensions
  if (allowedExtensions_.find(mySearchExtension)!=std::string::npos) {
    // Try to find if the extension is already registered
    std::vector<std::string>::iterator anExtension = find(fileTypeV_.begin(), fileTypeV_.end(), myExtension);
    if (anExtension == fileTypeV_.end()) {
      fileTypeV_.push_back(myExtension);
      result=true;
    }
  } else {
    std::cerr << "Extension '" << myExtension << "' is not allowed for canvas saving: new file type not added." << std::endl;
  }

  return result;
}

void RootWImage::saveSummary(std::string baseName, TFile* myTargetFile) {
  if (!myCanvas_) return;
  baseName += name_;
  saveSummaryLoop(myCanvas_.get(), baseName, myTargetFile);
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
	       (myClass=="TArrow") ||
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

RootWContent::RootWContent(std::string title, bool visible /*=true*/) {
  title_ = title;
  targetDirectory_ = "";
  visible_ = visible;
}

void RootWContent::setTargetDirectory(std::string newTargetDirectory) {
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

void RootWContent::addItem(std::unique_ptr<RootWItem> newItem) {
  addItem(newItem.release());
}

void RootWContent::setTitle(std::string newTitle) {
  title_ = newTitle;
}

void RootWContent::addParagraph(std::string parText) {
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

RootWText& RootWContent::addText(std::string newText) {
  RootWText* newTextItem = new RootWText(newText);
  addItem(newTextItem);
  return (*newTextItem);
}

RootWInfo& RootWContent::addInfo() {
  RootWInfo* newInfo = new RootWInfo();
  addItem(newInfo);
  return (*newInfo);
}

RootWInfo& RootWContent::addInfo(std::string description)  {
  RootWInfo* newInfo = new RootWInfo(description);
  addItem(newInfo);
  return (*newInfo);
}

RootWInfo& RootWContent::addInfo(std::string description, std::string value)  {
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

RootWImage& RootWContent::addImage(std::unique_ptr<TCanvas> myCanvas, int witdh, int height) {
  RootWImage* newImage = new RootWImage(std::move(myCanvas), witdh, height);
  addItem(newImage);
  return (*newImage);
}

RootWImage& RootWContent::addImage(std::unique_ptr<TCanvas> myCanvas, int witdh, int height, std::string relativeHtmlDirectory) {
  RootWImage* newImage = new RootWImage(std::move(myCanvas), witdh, height, relativeHtmlDirectory);
  addItem(newImage);
  return (*newImage);
}

RootWTextFile& RootWContent::addTextFile() {
  RootWTextFile* newTextFile = new RootWTextFile();
  addItem(newTextFile);
  return (*newTextFile);
}

RootWTextFile& RootWContent::addTextFile(std::string newFileName) {
  RootWTextFile* newTextFile = new RootWTextFile(newFileName);
  addItem(newTextFile);
  return (*newTextFile);
}

RootWTextFile& RootWContent::addTextFile(std::string newFileName, std::string newDescription) {
  RootWTextFile* newTextFile = new RootWTextFile(newFileName, newDescription);
  addItem(newTextFile);
  return (*newTextFile);
}

RootWBinaryFile& RootWContent::addBinaryFile() {
  RootWBinaryFile* newBinaryFile = new RootWBinaryFile();
  addItem(newBinaryFile);
  return (*newBinaryFile);
}

RootWBinaryFile& RootWContent::addBinaryFile(std::string newFileName) {
  RootWBinaryFile* newBinaryFile = new RootWBinaryFile(newFileName);
  addItem(newBinaryFile);
  return (*newBinaryFile);
}

RootWBinaryFile& RootWContent::addBinaryFile(std::string newFileName, std::string newDescription) {
  RootWBinaryFile* newBinaryFile = new RootWBinaryFile(newFileName, newDescription);
  addItem(newBinaryFile);
  return (*newBinaryFile);
}

RootWBinaryFile& RootWContent::addBinaryFile(std::string newFileName, std::string newDescription, std::string newOriginalFile) {
  RootWBinaryFile* newBinaryFile = new RootWBinaryFile(newFileName, newDescription, newOriginalFile);
  addItem(newBinaryFile);
  return (*newBinaryFile);
}


std::ostream& RootWContent::dump(std::ostream& output) {
  std::vector<RootWItem*>::iterator it;
  RootWItem* myItem;
  RootWImage* myImage;
  RootWFile* myFile;
  TFile* summaryFile = getSummaryFile();
  std::string baseName="";
  if (page_ != nullptr) baseName = page_->getTitle() + ".";

  //std::cerr << "Content: " << title_ <<endl; //debug
  output << "<h2 class=\"hidingTitle\">"<<title_<<"</h2>" << std::endl;
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
        std::cout << "WARNING: this should never happen. contact the author immediately!" << std::endl;
      }
    }
    if (myItem->isFile()) {
      if ( (myFile=dynamic_cast<RootWFile*>(myItem))) { // CUIDADO This is not good OOP for God's sake!!!
        myFile->setTargetDirectory(targetDirectory_);
      } else if (RootWFileList* myFileList=dynamic_cast<RootWFileList*>(myItem)) {
        myFileList->setTargetDirectory(targetDirectory_);
      } else {
        std::cout << "WARNING: this should never happen. contact the author immediately!" << std::endl;
      }
    }
    if (myItem->isPlacedBelow()) { output << "<div class=\"clearer\"> <br> <br> &nbsp;</div>"; }
    myItem->dump(output);
  }
  output << "<div class=\"clearer\">&nbsp;</div>";
  output << "</div>" << std::endl;
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

RootWPage::RootWPage(std::string title) {
  setTitle(title);
  site_ = NULL;
  targetDirectory_ = "";
  alert_ = 0;
}

RootWPage::~RootWPage() {
}

void RootWPage::setTargetDirectory(std::string newTargetDirectory) {
  targetDirectory_ = newTargetDirectory;
}

void RootWPage::setTitle(std::string newTitle) {
  title_ = newTitle;
  if (address_=="") address_ = title_+".html";
}

std::string RootWPage::getTitle() {
  return title_;
}

void RootWPage::setAddress(std::string newAddress) {
  address_ = newAddress;
}

std::string RootWPage::getAddress() {
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

void RootWPage::addContent(std::unique_ptr<RootWContent> newContent) {
  addContent(newContent.release());
}

RootWContent& RootWPage::addContent(std::string title, bool visible /*=true*/) {
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

std::ostream& RootWPage::dump(std::ostream& output) {
  bool fakeSite=false;
  if (site_==NULL) {
    fakeSite = true;
    site_ = new RootWSite();
    site_->addPage(this);
  }

  // Header standard part, with the menu list, etc...
  // up to opening the primaryContent div
  site_->dumpHeader(output, this);

  std::vector<RootWContent*>::iterator it;
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
  toolkitName_ = toolkit_name;
  toolkitGithub_ = toolkit_github;
  toolkitContributors_ = toolkit_contributors;
  revision_="";
  targetDirectory_ = ".";
  summaryFile_ = nullptr;
  createSummaryFile_ = true;
  summaryFileName_ = "summary.root";
}

RootWSite::RootWSite(std::string title) {
  title_ = title;
  comment_ = "";
  commentLink_ = "";
  toolkitName_ = toolkit_name;
  toolkitGithub_ = toolkit_github;
  toolkitContributors_ = toolkit_contributors;
  revision_="";
  targetDirectory_ = ".";
  summaryFile_ = nullptr;
  createSummaryFile_ = true;
  summaryFileName_ = "summary.root";
}

RootWSite::RootWSite(std::string title, std::string comment) {
  title_ = title;
  comment_ = comment;
  commentLink_ = "";
  toolkitName_ = toolkit_name;
  toolkitGithub_ = toolkit_github;
  toolkitContributors_ = toolkit_contributors;
  revision_="";
  targetDirectory_ = ".";
  summaryFile_ = nullptr;
  createSummaryFile_ = true;
  summaryFileName_ = "summary.root";
}

RootWSite::~RootWSite() {
}

void RootWSite::setTitle(std::string newTitle) {
  title_ = newTitle;
}

void RootWSite::setComment(std::string newComment) {
  comment_ = newComment;
}

void RootWSite::setCommentLink(std::string newCommentLink) {
  commentLink_ = newCommentLink;
}

std::string RootWSite::getTitle() {
  return title_;
}

std::string RootWSite::getComment() {
  return comment_;
}

std::string RootWSite::getCommentLink() {
  return commentLink_;
}

std::string RootWSite::getRevision() {
  return revision_;
}

void RootWSite::setRevision(std::string newRevision) {
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
    std::vector<RootWPage*>::iterator it = pageList_.begin();
    std::vector<RootWPage*>::iterator beforeThis = pageList_.end();
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

RootWPage& RootWSite::addPage(std::string newTitle, int relevance /* = least_relevant */ ) {
  RootWPage* newPage = new RootWPage(newTitle);
  addPage(newPage, relevance);
  return (*newPage);
}

std::ostream& RootWSite::dumpFooter(std::ostream& output) {
  output  << "          </div>" << std::endl
	  << "        </div>" << std::endl
	  << "        <div class=\"clear\"></div>" << std::endl
	  << "      </div>" << std::endl
	  << "      <div id=\"footer\">" << std::endl;

  output << "<p> &copy; tkLayout developers: " 
	 << "<a href=\"" << toolkitContributors_ << "\">" << toolkitContributors_ << "</a>"
	 << "</p>" 
	 << std::endl;


  time_t rawtime;
  time ( &rawtime );
  output << "        <p>Page created on "<< asctime(gmtime ( &rawtime )) << " GMT</p>" << std::endl
	 << "        <p>by <a href=\""<< toolkitGithub_ <<"\">"<< toolkitName_ <<"</a>";
  if (revision_!="")
    output << " revision " << revision_;
  output << "</p>" << std::endl
	 << "      </div>" << std::endl
	 << "    </div>" << std::endl
	 << "  </body>" << std::endl
	 << "</html>" << std::endl;
  return output;
}

std::ostream& RootWSite::dumpHeader(std::ostream& output, RootWPage* thisPage) {
  output << "<html xmlns=\"http://www.w3.org/1999/xhtml\">" << std::endl
	 << "  <head>" << std::endl
	 << "    <meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />" << std::endl
	 << "    <title>"<< title_ << " - " << thisPage->getTitle() <<"</title>" << std::endl
         << "    <meta http-equiv=\"cache-control\" content=\"max-age=0\" />" << std::endl
         << "    <meta http-equiv=\"cache-control\" content=\"no-cache\" />" << std::endl
         << "    <meta http-equiv=\"expires\" content=\"0\" />" << std::endl
         << "    <meta http-equiv=\"expires\" content=\"Tue, 01 Jan 1980 1:00:00 GMT\" />" << std::endl
         << "    <meta http-equiv=\"pragma\" content=\"no-cache\" />" << std::endl
	 << "    <meta name=\"keywords\" content=\"CERN CMS tracker upgrade\" />" << std::endl
	 << "    <meta name=\"description\" content=\"CMS Tracker upgrade summary page\" />" << std::endl
	 << "    <link href=\"../style/default.css\" rel=\"stylesheet\" type=\"text/css\" />" << std::endl
	 << "    <link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"../style/images/favicon.ico\">" << std::endl
	 << " </head>" << std::endl
	 << "  <script type=\"text/javascript\" src=\"../style/functions.js\"></script>" << std::endl
	 << "  <body onload=\"preparePageGoodies();\">" << std::endl
	 << "    <div id=\"outer\">" << std::endl
	 << "      <div id=\"header\">" << std::endl
	 << "        <h1>" << std::endl
	 << "    <a href=\"index.html\">" << title_ << "</a>" << std::endl
	 << "        </h1>" << std::endl;
  if (commentLink_=="") {
    output << "    <h2>"<<comment_<<"</h2>" << std::endl;
  } else {
    output << "    <h2>"
	   << "<a href=\"" << commentLink_ << "\">"
	   << comment_
	   << "</a>"
	   << "</h2>" << std::endl;
  }
  output << "      </div>" << std::endl
	 << "      <div id=\"menu\">" << std::endl
	 << "        <ul>" << std::endl;
  for (std::vector<RootWPage*>::iterator it = pageList_.begin();
       it !=pageList_.end(); ++it) {
    output << "          <li";
    if ((*it)==(thisPage)) output << " class=\"current\"";
    output << ">";
    output  << "    <a href=\"" << (*it)->getAddress() << "\" ";
    if ((*it)->getAlert()>0) {
      std::ostringstream anAlarmColor;
      anAlarmColor.str("");
      anAlarmColor << "rgb(255, "
		   << int(255 * (1-(*it)->getAlert())) // the bigger the alert int, the more red we want to be!
	                                               // rgb(255,255,0) is yellow, rgb(255,0,0) is red
		   << ", 0)";
      output << " style = 'background:"
		   << anAlarmColor.str()
		   << ";' ";
    }
    output << ">";
    output << (*it)->getTitle();
    output <<"</a>" << "          </li>" << std::endl;
  }
  output << "        </ul>" << std::endl
	 << "      </div>" << std::endl
	 << "      <div id=\"content\">" << std::endl
	 << "        <div id=\"primaryContentContainer\">" << std::endl
	 << "          <div id=\"primaryContent\">" << std::endl;

  return output;

}

bool RootWSite::prepareTargetDirectory() {
  // Check if the directory already exists
  if (boost::filesystem::exists( targetDirectory_ )) {
    if (! boost::filesystem::is_directory(targetDirectory_) ) {
      std::cerr << "A file named " << targetDirectory_ << " already exists. Cannot create target directory" << std::endl;
      return false;
    }
  } else {
    // If the directory does not exist, then try to create it
    if (!boost::filesystem::create_directory(targetDirectory_)) {
      std::cerr << "Couldn't create directory " << targetDirectory_ << std::endl;
      return false;
    }
  }
  return true;
}

bool RootWSite::makeSite(bool verbose) {
  std::ofstream myPageFile;
  RootWPage* myPage;
  std::string myPageFileName;

  prepareTargetDirectory();

  std::vector<RootWPage*>::iterator it;
  if (createSummaryFile_) {
    summaryFile_.reset(new TFile(Form("%s/%s",
				  targetDirectory_.c_str(),
				      summaryFileName_.c_str()), "RECREATE"));
  } else summaryFile_ = nullptr;
  for (it=pageList_.begin(); it!=pageList_.end(); it++) {
    myPage = (*it);
    if (verbose) std::cout << " " << myPage->getTitle() << std::flush;
    myPageFileName = targetDirectory_+"/"+myPage->getAddress();
    myPageFile.open(myPageFileName.c_str(), std::ios::out);
    myPage->setTargetDirectory(targetDirectory_);
    myPage->dump(myPageFile);
    myPageFile.close();
    if (it==pageList_.begin()) {
      namespace fs = boost::filesystem;
      fs::path const from_path = myPageFileName;
      fs::path const to_path = targetDirectory_+"/index.html";
#if BOOST_VERSION >= 106000
      fs::copy_file(from_path, to_path, fs::copy_options::overwrite_existing);
#else
      fs::copy_file(from_path, to_path, fs::copy_option::overwrite_if_exists);
#endif
    }
  }
  if (verbose) std::cout << " ";
  if (summaryFile_) summaryFile_->Close();

  return true;
}

TFile* RootWSite::getSummaryFile() {
  return summaryFile_.get();
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


RootWItem* RootWItemCollection::getItem(std::string itemName) {
  RootWItem* result;
  result=itemCollection_[itemName];
  if (result!=NULL) result->taken=true;
  return result;
}

void RootWItemCollection::addItem(RootWItem* anItem, std::string itemName) {
  RootWItem* findItem;
  findItem=itemCollection_[itemName];
  if (findItem!=NULL) {
    std::cerr << "Warning in RootWItemCollection::addItem(): item named " << itemName << " already exists! Overwriting it." << std::endl;
  }
  itemCollection_[itemName]=anItem;
}

std::vector<RootWItem*> RootWItemCollection::getOtherItems() {
  std::vector<RootWItem*> result;
  std::map<std::string, RootWItem*>::iterator itItem=itemCollection_.begin();
  std::map<std::string, RootWItem*>::iterator endItem=itemCollection_.end();
  for (; itItem!=endItem; ++itItem) {
    if (!(*itItem).second->taken) result.push_back((*itItem).second);
  }
  return result;
}

//*******************************************//
// RootWTextFile                             //
//*******************************************//

std::ostream& RootWTextFile::dump(std::ostream& output) {

  if (fileName_=="") {
    std::cerr << "Warning: RootWTextFile::dump() was called without prior setting the destination file name" << std::endl;
    return output;
  }

  std::ofstream outputFile;
  std::string destinationFileName = targetDirectory_ +"/" + fileName_;
  outputFile.open(destinationFileName.c_str());
  // TODO: add a check here if the file was correctly opened
  // outputFile << myText_.str();
  // myText_ >> outputFile;
  outputFile << myText_.str() << std::endl;
  outputFile.close();

  output << "<b>" << description_ << ":</b> <a href=\""
         << fileName_ << "\">"
         << fileName_ << "</a></tt><br/>";
  return output;
};


//*******************************************//
// RootWBinaryFile                           //
//*******************************************//

std::ostream& RootWBinaryFile::dump(std::ostream& output) {
  if ((originalFileName_=="")&&(!noCopy_)) {
    std::cerr << "Warning: RootWBinaryFile::dump() was called without prior setting the original file name" << std::endl;
    return output;
  }
  if (fileName_=="") {
    std::cerr << "Warning: RootWBinaryFile::dump() was called without prior setting the destination file name" << std::endl;
    return output;
  }

  std::string destinationFileName = targetDirectory_ +"/" + fileName_;

  if ((!noCopy_)&&(boost::filesystem::exists(originalFileName_) && originalFileName_ != destinationFileName)) { // CUIDADO: naive control on copy on itself. it only matches the strings, not taking into account relative paths and symlinks
    try {
      if (boost::filesystem::exists(destinationFileName))
        boost::filesystem::remove(destinationFileName);
      boost::filesystem::copy_file(originalFileName_, destinationFileName);
    } catch (boost::filesystem::filesystem_error& e) {
      std::cerr << e.what() << std::endl;
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

std::ostream& RootWBinaryFileList::dump(std::ostream& output) {
  if (originalFileNames_.empty()) {
    std::cerr << "Warning: RootWBinaryFileList::dump() was called without prior setting the original file names" << std::endl;
    return output;
  }
  if (fileNames_.empty()) {
    std::cerr << "Warning: RootWBinaryFileList::dump() was called without prior setting the destination file names" << std::endl;
    return output;
  }
  if (fileNames_.size() != originalFileNames_.size()) {
    std::cerr << "Warning: RootWBinaryFileList::dump() was called with original and destination file name lists differing in length" << std::endl;
    return output;
  }

  auto it = originalFileNames_.begin();
  for (auto fn : fileNames_) {
    std::string destinationFileName = targetDirectory_;
    destinationFileName += "/" + fn; //path.back();
    try {
      if (boost::filesystem::exists(destinationFileName)) boost::filesystem::remove(destinationFileName);
      boost::filesystem::copy_file(*it++, destinationFileName);
    }
    catch (boost::filesystem::filesystem_error& e) {
      std::cerr << e.what() << std::endl;
      return output;
    }
  }

  std::vector<std::string> cleanedUpFileNames;
  std::transform(originalFileNames_.begin(), originalFileNames_.end(), std::back_inserter(cleanedUpFileNames), [](const std::string& s) {
      auto pos = s.find("stdinclude");
      if (s.find("xml") != std::string::npos) pos = s.rfind("/") + 1;
      return pos != std::string::npos ? s.substr(pos) : s;
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

std::string RootWInfo::setValue(std::string newText) {
  value_ = newText;
  return newText;
}

std::string RootWInfo::setValue(int number) {
  std::stringstream myNum_;
  myNum_.clear();
  myNum_ << std::dec << number;
  return(setValue(myNum_.str()));
}

std::string RootWInfo::setValue(double number, int precision) {
  std::stringstream myNum_;
  myNum_.clear();
  myNum_ << std::dec << std::fixed << std::setprecision(precision) << number;
  return(setValue(myNum_.str()));
}

std::string RootWInfo::setValueSci(double number, int precision) {
  std::stringstream myNum_;
  myNum_.clear();
  int nearestLog = floor(log10(number));
  double mantissa = number/pow(10, nearestLog);
  myNum_.unsetf(std::ios_base::floatfield);
  myNum_.precision(5);
  myNum_ << std::dec << mantissa;
  myNum_ << "&times;10<sup>" << nearestLog << "</sup>";
  return(setValue(myNum_.str()));
}

std::string RootWInfo::appendValue(std::string aString) {
  std::stringstream myNum_;
  return setValue(value_ + aString);
}

std::ostream& RootWInfo::dump(std::ostream& output) {
  std::ofstream outputFile;

  output << "<b>" << description_ << ":</b> "
         << value_ << "</a></tt><br/>";
  return output;
}

//*******************************************//
// RootWGraphViz                             //
//*******************************************//

std::ostream& RootWGraphViz::dump(std::ostream& output) {
  std::string myDescription = description_;
  description_ += " (GraphViz)";
  RootWTextFile::dump(output);
  int result = system("which dot > /dev/null");
  if (result==0) {
    std::string svgFileName = fileName_+".svg";
    std::string command = "dot -Tsvg " + targetDirectory_ +"/" + fileName_ + " > " + targetDirectory_ + "/" + svgFileName;
    result = system(command.c_str());
    if (result==0) {
      RootWBinaryFile* mySvg = new RootWBinaryFile(svgFileName, myDescription+" (SVG)" );
      mySvg->setNoCopy(true);
      mySvg->dump(output);
    }
  }
  return output;
}

/*
 * RootWContent.cc
 *
 *  Created on: 20. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */
#include "RootWContent.h"

#include "RootWBinaryFile.h"
#include "RootWBinaryFileList.h"
#include "RootWGraphVizFile.h"
#include "RootWInfo.h"
#include "RootWImage.h"
#include "RootWItem.h"
#include "RootWTable.h"
#include "RootWText.h"
#include "RootWTextFile.h"
#include "TCanvas.h"

//
// Default constructor
//
RootWContent::RootWContent() {

  m_title = "Untitled";
  m_targetDirectory = "";
  m_visible = true;
}

//
// Constructor defining content title and visibility
//
RootWContent::RootWContent(std::string title, bool visible /*=true*/) {
  m_title = title;
  m_targetDirectory = "";
  m_visible = visible;
}

//
// Default desctructor
//
RootWContent::~RootWContent() {

  m_itemList.clear();
}

//
// Method adding any item derived from RootWItem: Text, Info, Image, ... - authorship of pointer is handed over to this class
//
void RootWContent::addItem(std::unique_ptr<RootWImage> newImage)
{
  std::unique_ptr<RootWItem> newItem(std::move(newImage));
  m_itemList.push_back(std::move(newItem));
}
void RootWContent::addItem(std::unique_ptr<RootWInfo> newInfo)
{
  std::unique_ptr<RootWItem> newItem(std::move(newInfo));
  m_itemList.push_back(std::move(newItem));
}
void RootWContent::addItem(std::unique_ptr<RootWTable> newTable)
{
  std::unique_ptr<RootWItem> newItem(std::move(newTable));
  m_itemList.push_back(std::move(newItem));
}
void RootWContent::addItem(std::unique_ptr<RootWText> newText)
{
  std::unique_ptr<RootWItem> newItem(std::move(newText));
  m_itemList.push_back(std::move(newItem));
}
void RootWContent::addItem(std::unique_ptr<RootWTextFile> newTextFile)
{
  std::unique_ptr<RootWItem> newItem(std::move(newTextFile));
  m_itemList.push_back(std::move(newItem));
}
void RootWContent::addItem(std::unique_ptr<RootWBinaryFile> newBinaryFile)
{
  std::unique_ptr<RootWItem> newItem(std::move(newBinaryFile));
  m_itemList.push_back(std::move(newItem));
}
void RootWContent::addItem(std::unique_ptr<RootWBinaryFileList> newBinaryFileList)
{
  std::unique_ptr<RootWItem> newItem(std::move(newBinaryFileList));
  m_itemList.push_back(std::move(newItem));
}
void RootWContent::addItem(std::unique_ptr<RootWGraphVizFile> newGraphVizFile)
{
  std::unique_ptr<RootWItem> newItem(std::move(newGraphVizFile));
  m_itemList.push_back(std::move(newItem));
}

//
// Add RootWText in form of paragraph - define text to be given on the output
//
void RootWContent::addParagraph(std::string parText) {

  std::unique_ptr<RootWText> newText(new RootWText);
  //newText->addText("<p>");
  newText->addText(parText);
  newText->addText("<br/>\n");
  //newText->addText("</p>\n");

  std::unique_ptr<RootWItem> newItem(std::move(newText));
  m_itemList.push_back(std::move(newItem));
}

//
// Add RootWText - content can be changed using the reference
//
RootWText& RootWContent::addText() {

  std::unique_ptr<RootWItem> newText(new RootWText());
  m_itemList.push_back(std::move(newText));

  return (*(dynamic_cast<RootWText*>(m_itemList.back().get())));
}
RootWText& RootWContent::addText(std::string newText) {

  std::unique_ptr<RootWItem> newTextItem(new RootWText(newText));
  m_itemList.push_back(std::move(newTextItem));
  return (*(dynamic_cast<RootWText*>(m_itemList.back().get())));
}

//
// Add RootWInfo - content can be changed using the reference
//
RootWInfo& RootWContent::addInfo() {

  std::unique_ptr<RootWItem> newInfo(new RootWInfo());
  m_itemList.push_back(std::move(newInfo));
  return (*(dynamic_cast<RootWInfo*>(m_itemList.back().get())));
}
RootWInfo& RootWContent::addInfo(std::string description)  {

  std::unique_ptr<RootWItem> newInfo(new RootWInfo(description));
  m_itemList.push_back(std::move(newInfo));
  return (*(dynamic_cast<RootWInfo*>(m_itemList.back().get())));
}
RootWInfo& RootWContent::addInfo(std::string description, std::string value)  {

  std::unique_ptr<RootWItem> newInfo(new RootWInfo(description, value));
  m_itemList.push_back(std::move(newInfo));
  return (*(dynamic_cast<RootWInfo*>(m_itemList.back().get())));
}

//
// Add RootWTable - content can be changed using the reference
//
RootWTable& RootWContent::addTable() {

  std::unique_ptr<RootWItem> newTable(new RootWTable());
  m_itemList.push_back(std::move(newTable));
  return (*(dynamic_cast<RootWTable*>(m_itemList.back().get())));
}

//
// Add RootWImage - content can be changed using the reference
//
RootWImage& RootWContent::addImage(const TCanvas& myCanvas) {

  int width = myCanvas.GetWindowWidth();
  int height= myCanvas.GetWindowHeight();
  std::unique_ptr<RootWItem> newImage(new RootWImage(myCanvas, width, height));
  m_itemList.push_back(std::move(newImage));
  return (*(dynamic_cast<RootWImage*>(m_itemList.back().get())));
}
RootWImage& RootWContent::addImage(const TCanvas& myCanvas, int width, int height) {

  std::unique_ptr<RootWItem> newImage(new RootWImage(myCanvas, width, height));
  m_itemList.push_back(std::move(newImage));
  return (*(dynamic_cast<RootWImage*>(m_itemList.back().get())));
}
RootWImage& RootWContent::addImage(const TCanvas& myCanvas, int width, int height, std::string relativeHtmlDirectory) {

  std::unique_ptr<RootWItem> newImage(new RootWImage(myCanvas, width, height, relativeHtmlDirectory));
  m_itemList.push_back(std::move(newImage));
  return (*(dynamic_cast<RootWImage*>(m_itemList.back().get())));
}

//
// Find image & return canvas assigned to it
//
const TCanvas& RootWContent::findImage(std::string imageName) const {

  const TCanvas* canvas = nullptr;
  for (const auto& iContent : m_itemList) {

    const RootWImage* image = dynamic_cast<const RootWImage*>(iContent.get());

    if (image!=nullptr) {
      if (image->getName()==imageName) {
        canvas = &(image->getCanvas());
        break;
      }
    }
  }

  // Check that TCanvas defined for searched image
  if (canvas) return *canvas;
  else throw std::invalid_argument( "RootWContent::findImage - TCanvas not defined (null reference), check!!!" );
}

//
// Add RootWTextFile - content can be changed using the reference
//
RootWTextFile& RootWContent::addTextFile(std::string newFileName) {

  std::unique_ptr<RootWItem> newTextFile(new RootWTextFile(newFileName));
  m_itemList.push_back(std::move(newTextFile));
  return (*(dynamic_cast<RootWTextFile*>(m_itemList.back().get())));
}
RootWTextFile& RootWContent::addTextFile(std::string newFileName, std::string newDescription) {

  std::unique_ptr<RootWItem> newTextFile(new RootWTextFile(newFileName, newDescription));
  m_itemList.push_back(std::move(newTextFile));
  return (*(dynamic_cast<RootWTextFile*>(m_itemList.back().get())));
}

//
// Add RootWBinaryFile - content can be changed using the reference
//
RootWBinaryFile& RootWContent::addBinaryFile(std::string newFileName) {

  std::unique_ptr<RootWItem> newBinaryFile(new RootWBinaryFile(newFileName));
  m_itemList.push_back(std::move(newBinaryFile));
  return (*(dynamic_cast<RootWBinaryFile*>(m_itemList.back().get())));
}
RootWBinaryFile& RootWContent::addBinaryFile(std::string newFileName, std::string newDescription) {

  std::unique_ptr<RootWItem> newBinaryFile(new RootWBinaryFile(newFileName, newDescription));
  m_itemList.push_back(std::move(newBinaryFile));
  return (*(dynamic_cast<RootWBinaryFile*>(m_itemList.back().get())));
}
RootWBinaryFile& RootWContent::addBinaryFile(std::string newFileName, std::string newDescription, std::string newOriginalFile) {

  std::unique_ptr<RootWItem> newBinaryFile(new RootWBinaryFile(newFileName, newDescription, newOriginalFile));
  m_itemList.push_back(std::move(newBinaryFile));
  return (*(dynamic_cast<RootWBinaryFile*>(m_itemList.back().get())));
}

//
// Dump method - dump the content block in html format on a web-page
//
std::ostream& RootWContent::dump(std::ostream& output) {

  RootWImage*    myImage    = nullptr;
  RootWFile*     myFile     = nullptr;
  RootWFileList* myFileList = nullptr;

  //std::cerr << "Content: " << m_title <<std::endl;; //debug
  output << "      <h2 class=\"hidingTitle\">"<<m_title<<"</h2>" << std::endl;;
  if (m_visible) {
    output << "      <div class=\"hideable\"> ";
  } else {
    output << "      <div class=\"hidden\"> ";
  }
  for (auto& myItem : m_itemList) {

    if ((myImage=dynamic_cast<RootWImage*>(myItem.get()))!=nullptr) {
      myImage->setTargetDirectory(m_targetDirectory);
      myImage->saveFiles(c_thumbSmallSize, c_thumbSmallSize);
    }
    if ((myFile=dynamic_cast<RootWFile*>(myItem.get()))!=nullptr) {
      myFile->setTargetDirectory(m_targetDirectory);
    }
    if ((myFileList=dynamic_cast<RootWFileList*>(myItem.get()))!=nullptr) {
      myFileList->setTargetDirectory(m_targetDirectory);
    }
    myItem->dump(output);
  }
  output << "      <div class=\"clearer\">&nbsp;</div>";
  output << "      </div>" << std::endl;;
  return output;
}


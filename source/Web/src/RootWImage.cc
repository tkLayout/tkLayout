/*
 * RootWImage.cc
 *
 *  Created on: 21. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.cc to standalone file
 */
#include "RootWImage.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

#include "MessageLogger.h"
#include "TCanvas.h"
#include "TView.h"
#include "TError.h"

//
// Define static variables
//
int                         RootWImage::s_imageCounter = 0;
std::map <std::string, int> RootWImage::s_imageNameCounter;

//
// Default constructor - canvas is cloned and dynamically allocated
//
RootWImage::RootWImage(const TCanvas& myCanvas) {
  s_imageCounter++;
  setZoomedSize(myCanvas.GetWindowWidth(), myCanvas.GetWindowHeight());
  setCanvas(myCanvas);
  m_relativeHtmlDirectory = "";
  m_targetDirectory = "";
  m_comment = "";
  m_name = "img";
  m_allowedExtensions = c_defaultAllowedExtensions;
  setDefaultExtensions();
}

//
// Constructor - clone canvas and define image width & height
//
RootWImage::RootWImage(const TCanvas& myCanvas, int witdh, int height) {
  s_imageCounter++;
  setCanvas(myCanvas);
  setZoomedSize(witdh, height);
  m_relativeHtmlDirectory = "";
  m_targetDirectory = "";
  m_comment = "";
  m_name = "img";
  m_allowedExtensions = c_defaultAllowedExtensions;
  setDefaultExtensions();
}

//
// Constructor - clone canvas and define image width & height, define relativeHtml directory wrt target directory
//
RootWImage::RootWImage(const TCanvas& myCanvas, int witdh, int height, std::string relativehtmlDirectory) {
  s_imageCounter++;
  setCanvas(myCanvas);
  setZoomedSize(witdh, height);
  setRelativeHtmlDirectory(relativehtmlDirectory);
  m_targetDirectory = "";
  m_comment = "";
  m_name = "img";
  m_allowedExtensions = c_defaultAllowedExtensions;
  setDefaultExtensions();
}

//
// Default destructor
//
RootWImage::~RootWImage() {}

//
// Private method - set default supported format: png (automatically) + those defined here
//
void RootWImage::setDefaultExtensions() {
   //addExtension("C");
   addExtension("pdf");
   addExtension("root");
}

//
// Return reference to the saved TCanvas to get access to plots etc. saved in the canvas
//
const TCanvas& RootWImage::getCanvas() const
{
  const TCanvas* canvas = nullptr;

  // Check that TCanvas defined for this image
  if (m_canvas.get()) canvas = m_canvas.get();
  else throw std::invalid_argument( "RootWImage::getCanvas - TCanvas not defined (null reference), check!!!" );;

  return *canvas;
}

//
// Private helper method cloning the given canvas and setting its properties
//
void RootWImage::setCanvas(const TCanvas& myCanvas) {

  // Create exact copy of TCanvas & allocate space for it -> assign to member variable
  m_canvas.reset((TCanvas*)myCanvas.DrawClone());

  // Set unique name
  std::ostringstream canvasName;
  canvasName << "canvas" << std::setfill('0') << std::setw(3) << s_imageCounter;
  m_canvas->SetName(canvasName.str().c_str());

  TView* myView = m_canvas->GetView();
  if (myView) {
    TView* newView = m_canvas->GetView();
    if (newView) {
      double min[3], max[3];
      Int_t irep;
      newView->SetView(myView->GetLongitude(),myView->GetLatitude(),myView->GetPsi(), irep);
      myView->GetRange(min, max);
      newView->SetRange(min, max);
    }
  }
}

//
// Method preceding dump to save files on disk and prepare images
//
std::string RootWImage::saveFiles(int smallWidth, int smallHeight, int largeWidth, int largeHeight) {
  setZoomedSize(largeWidth, largeHeight);
  return saveFiles(smallWidth, smallHeight);
}

std::string RootWImage::saveFiles(int smallWidth, int smallHeight) {

  // Suppress print info coming from ROOT
  gErrorIgnoreLevel = kWarning;

  auto myImageSize = makeSizeCode(smallWidth,smallHeight,m_zoomedWidth,m_zoomedHeight);

  std::ostringstream tmpCanvasName("");
  tmpCanvasName << m_name << std::setfill('0') << std::setw(3) << s_imageNameCounter[m_name]++;

  std::string canvasName = tmpCanvasName.str();
  m_canvas->SetName(canvasName.c_str());
  if (m_fileSaved[myImageSize]) return m_text[myImageSize];

  std::stringstream thisText;
  thisText.clear();
  if (m_relativeHtmlDirectory == "") m_relativeHtmlDirectory = ".";
  if (m_targetDirectory       == "" ) m_targetDirectory = ".";

  if (!m_canvas.get()) return "";

  //std::string canvasName = m_canvas->GetName();
  std::string smallCanvasFileName     = canvasName + "-" + myImageSize + "_small.png";
  std::string largeCanvasFileBaseName = canvasName;

  double wScale = double(smallWidth) / double(m_zoomedWidth);
  double hScale = double(smallHeight) / double(m_zoomedHeight);
  double myScale;

  myScale =  (hScale > wScale) ? hScale : wScale;
  int imgW = int (m_zoomedWidth * myScale);
  int imgH = int (m_zoomedHeight * myScale);
  int canW = int(imgW * c_thumbCompression);
  int canH = int(imgH * c_thumbCompression);

  std::string smallCanvasCompleteFileName = m_targetDirectory+"/"+smallCanvasFileName;
  std::string largeCanvasCompleteFileName = m_targetDirectory+"/"+largeCanvasFileBaseName + ".png";

  m_canvas->cd();
  m_canvas->SetCanvasSize(canW, canH);
  // std::cerr << "Saving " << largeCanvasCompleteFileName << std::endl; // debug
  m_canvas->Print(smallCanvasCompleteFileName.c_str());
  m_canvas->SetCanvasSize(m_zoomedWidth, m_zoomedHeight);
  m_canvas->Print(largeCanvasCompleteFileName.c_str());

  std::string fileTypeList;
  for (auto it=m_fileTypeV.begin(); it!=m_fileTypeV.end(); ++it) {
    largeCanvasCompleteFileName = m_targetDirectory+"/"+largeCanvasFileBaseName + "." + (*it);
    m_canvas->Print(largeCanvasCompleteFileName.c_str());
    if (it!=m_fileTypeV.begin()) fileTypeList+="|";
    fileTypeList+=(*it);
  }

  thisText << "<img class='wleftzoom' width='"<< imgW<<"' height='"<<imgH<<"' src='"
           << m_relativeHtmlDirectory << "/" << smallCanvasFileName <<"' alt='"<< m_comment
           <<"' onclick=\"popupDiv('" << m_relativeHtmlDirectory << "/" << largeCanvasFileBaseName
           << "', "<< m_zoomedWidth<<", "<< m_zoomedHeight << " , '"<< m_comment << "','"<<fileTypeList << "');\" />";

  m_text[myImageSize] = thisText.str();
  m_fileSaved[myImageSize] = true;
  m_lastSize = myImageSize;
  return thisText.str();
}

//
// Dump method - printing image in the html format with description etc. (first call of saveFiles is needed)
//
std::ostream& RootWImage::dump(std::ostream& output) {
  output << m_text[m_lastSize];
  return output;
}

RootWImageSize RootWImage::makeSizeCode(int sw, int sh, int lw, int lh) {
  std::stringstream aString;
  aString << sw <<"x" << sh << "-" << lw <<"x" << lh;
  return aString.str();
}

//
// Image will be saved in default formats (setDefaultExtensions) + all formats defined by this method (needs
// to be a supported format, cf. c_defaultAllowedExtensions)
//
bool RootWImage::addExtension(std::string myExtension) {
  // First check if the extension has not a "|" inside, which
  // would invalidate the search
  if (myExtension.find("|")!=std::string::npos) {
    logERROR("In bool RootWImage::addExtension(std::string myExtension) extension contains '|', which is not allowed");
    return false;
  }
  bool result = false;

  // Then create the real search pattern
  std::string mySearchExtension="|"+myExtension+"|";

  // Last: check if the extension is in the list of allowed extensions
  if (m_allowedExtensions.find(mySearchExtension)!=std::string::npos) {

    // Try to find if the extension is already registered
    auto anExtension = find(m_fileTypeV.begin(), m_fileTypeV.end(), myExtension);
    if (anExtension == m_fileTypeV.end()) {
      m_fileTypeV.push_back(myExtension);
      result=true;
    }
  } else {
    std::ostringstream message;
    message << "Extension '" << myExtension << "' is not allowed for canvas saving: new file type not added." << endl;
    logERROR(message);
  }

  return result;
}




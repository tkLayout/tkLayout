/*
 * RootWContent.h
 *
 *  Created on: 20. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */

#ifndef INCLUDE_ROOTWCONTENT_H_
#define INCLUDE_ROOTWCONTENT_H_

#include <memory>
#include <string>
#include <ostream>
#include <vector>

// Forward declaration
class TCanvas;
class RootWBinaryFile;
class RootWBinaryFileList;
class RootWGraphVizFile;
class RootWImage;
class RootWInfo;
class RootWItem;
class RootWTable;
class RootWText;
class RootWTextFile;

/*
 * @class RootWContent
 * Individual html blocks are put into RootWPage as instances of RootWContent class.
 * The following blocks are supported: RooWText (text), RootWImage (image), RoowInfo
 * (info container), RootWTable (table), RootWTextFile (link to text file), RootWBinaryFile
 * (link to binary file) or any other item derived from RootWItem. If addItem() method is used
 * to fill an item to the content, a unique_ptr needs to be used on input to handle correctly the
 * pointer ownership (use std::move() method).
 */
class RootWContent {

 public:

  //! Default constructor
  RootWContent();
  //! Constructor defining content title and visibility
  RootWContent(std::string title, bool visible=true);
  //! Default desctructor
  ~RootWContent();

  //! Dump method - dump the content block in html format on a web-page
  std::ostream& dump(std::ostream& output);

  //! Method adding any item derived from RootWItem: Text, Info, Image, ... - authorship
  //! of pointer is handed over to this class (unique_ptrs used)
  void addItem(std::unique_ptr<RootWImage> newImage);
  void addItem(std::unique_ptr<RootWInfo> newInfo);
  void addItem(std::unique_ptr<RootWTable> newTable);
  void addItem(std::unique_ptr<RootWText> newText);
  void addItem(std::unique_ptr<RootWTextFile> newTextFile);
  void addItem(std::unique_ptr<RootWBinaryFile> newBinaryFile);
  void addItem(std::unique_ptr<RootWBinaryFileList> newBinaryFileList);
  void addItem(std::unique_ptr<RootWGraphVizFile> newGraphVizFile);

  //! Add RootWText in form of paragraph - define text to be given on the output
  void addParagraph(std::string parText) ;

  //! Add RootWText - content can be changed using the reference
  RootWText& addText();
  RootWText& addText(std::string newText);

  //! Add RootWInfo - content can be changed using the reference
  RootWInfo& addInfo();
  RootWInfo& addInfo(std::string description);
  RootWInfo& addInfo(std::string description, std::string value);

  //! Add RootWTable - content can be changed using the reference
  RootWTable& addTable();

  //! Add RootWImage - content can be changed using the reference
  RootWImage& addImage(const TCanvas& myCanvas);
  RootWImage& addImage(const TCanvas& myCanvas, int width, int height);
  RootWImage& addImage(const TCanvas& myCanvas, int width, int height, std::string relativeHtmlDirectory);

  //! Add RootWTextFile - content can be changed using the reference
  RootWTextFile& addTextFile(std::string newFileName);
  RootWTextFile& addTextFile(std::string newFileName, std::string newDescription);

  //! Add RootWBinaryFile - content can be changed using the reference
  RootWBinaryFile& addBinaryFile(std::string newFileName);
  RootWBinaryFile& addBinaryFile(std::string newFileName, std::string newDescription);
  RootWBinaryFile& addBinaryFile(std::string newFileName, std::string newDescription, std::string newOriginalFile);

  //! Find image & return canvas assigned to it
  const TCanvas& findImage(std::string imageName) const;

  // Setter methods
  void setTargetDirectory(std::string newTargetDirectory) {m_targetDirectory = newTargetDirectory; }
  void setTitle(std::string newTitle)                     {m_title = newTitle;}

  // Getter methods
  std::string getTitle() const {return m_title;}

 private:

  std::vector<std::unique_ptr<RootWItem>> m_itemList; //! List of assigned web items (pointers owned by this class)

  bool        m_visible;         //!< Is object rolled-out on the web page (true) or one needs to click on it to roll it out (false)
  std::string m_title;           //!< Content title displayed on the web page
  std::string m_targetDirectory; //!< Target directory - set the same as target directory of the web-page

  const int c_thumbSmallSize = 120;

}; // Class

#endif /* INCLUDE_ROOTWCONTENT_H_ */

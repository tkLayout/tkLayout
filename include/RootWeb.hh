#ifndef _ROOTWEB_HH_
#define _ROOTWEB_HH_

#define BOOST_NO_CXX11_SCOPED_ENUMS
#define BOOST_NO_CXX11_SCOPED_ENUM

// standard includes
#include <map>
#include <sstream>
#include <string>
#include <list>
#include <limits>
// ROOT includes
#include <TCanvas.h>
#include <TFile.h>
#include <vector>
#include <memory>

#define THUMBSMALLSIZE 120
// The following is a list of allowed file etensions for TCanvas::SaveAs
// It should be separated, start and end with '|'
#define DEFAULTALLOWEDEXTENSIONS "|C|png|gif|svg|root|eps|pdf|ps|"

namespace RootWeb {
  static const std::string toolkit_name = "tkLayout";
  static const std::string toolkit_github = "https://github.com/tkLayout/tkLayout";
  static const std::string toolkit_contributors = "https://github.com/tkLayout/tkLayout/graphs/contributors";
  std::string cleanUpObjectName(const std::string&);
  static const int least_relevant = std::numeric_limits<int>::min();
  static const int most_relevant = std::numeric_limits<int>::max();
};

using namespace RootWeb;

class RootWItem {
public:
  virtual ~RootWItem() {};
  RootWItem() {taken=false; isPlacedBelow_=false;};
  virtual bool isTable() {return false;};
  virtual bool isImage() {return false;};
  virtual bool isText() {return false;};
  virtual bool isFile() {return false;};
  virtual std::ostream& dump(std::ostream& output) {return output;};
  bool taken;
  virtual bool isPlacedBelow() { return isPlacedBelow_; }
protected:
  bool isPlacedBelow_;
};

class RootWText: public RootWItem {
public:
  RootWText() {myText_.clear();};
  RootWText(std::string newText) {myText_.str(newText); };
  ~RootWText() {};
  void addText(std::string newText) { myText_ << newText; };
  std::ostream& dump(std::ostream& output) {output << myText_.str(); return output;};
  bool isText() {return true;};
protected:
  std::stringstream myText_;
};

class RootWInfo: public RootWItem {
public:
  RootWInfo() {description_="missing_description";value_="missing_value";};
  RootWInfo(std::string description) {description_=description; value_="missing_value";};
  RootWInfo(std::string description, std::string value) {description_=description;value_=value;};
  ~RootWInfo() {};
  std::string setDescription(std::string newText) { description_ = newText; return newText; };
  std::string setValue(std::string newText);
  std::string setValue(int number);
  std::string setValue(double number, int precision);
  std::string setValueSci(double number, int precision);
  std::string appendValue(std::string);
  std::string addValueText(std::string newText) {value_+=newText; return value_;};
  std::ostream& dump(std::ostream& output);
protected:
  std::string description_;
  std::string value_;
};

typedef std::map<std::pair<int,int>, std::string> rootWTableContent;
typedef std::map<std::pair<int,int>, int> rootWTableContentColor;
typedef std::map<std::pair<int,int>, bool> rootWTableContentBold;
class RootWTable : public RootWItem {
public:
  ~RootWTable() {};
  RootWTable(bool isPlacedBelow = false);
  void setContent(int row, int column, std::string content, const bool isBold = false, const int color = kBlack);
  void setContent(int row, int column, int number, const bool isBold = false, const int color = kBlack);
  void setContent(int row, int column, double number, int precision, const bool isBold = false, const int color = kBlack);
  void setContent(const rootWTableContent& newContent) { tableContent_ = newContent; };
  void setColor(const int row,const  int column, const int newColor);
  void setBold(const int row, const int column, const bool isBold);
  std::ostream& dump(std::ostream& output);
  std::pair<int, int> addContent(std::string content);
  std::pair<int, int> addContent(int number);
  std::pair<int, int> addContent(double number, int precision);
  std::pair<int, int> newLine();
  bool isTable() {return true;};
  const int maxRow() const { return maxRow_; }
  const int maxCol() const { return maxCol_; }
private:
  rootWTableContent tableContent_;
  rootWTableContentColor tableContentColor_;
  rootWTableContentBold tableContentBold_;
  int serialRow_, serialCol_;
  int maxRow_, maxCol_;
};

typedef std::string RootWImageSize;

class RootWImage : public RootWItem {
public:
  RootWImage();
  // TODO: the methods with TCanvas* (pointer) should be made obsolete
  RootWImage(std::unique_ptr<TCanvas> myCanvas, int witdh, int height);
  RootWImage(std::unique_ptr<TCanvas> myCanvas, int witdh, int height, std::string relativeHtmlDirectory); // TODO: is this used for real?
  void setCanvas(std::unique_ptr<TCanvas> myCanvas);
  void setComment(std::string newComment);
  void setName(std::string newName);
  std::string getName();
  void setZoomedSize(int witdh, int height);
  void setRelativeHtmlDirectory(std::string newDirectory);
  void setTargetDirectory(std::string newDirectory);
  std::string saveFiles(int smallWidth, int smallHeight);
  std::string saveFiles(int smallWidth, int smallHeight, int largeWidth, int largeHeight);
  std::ostream& dump(std::ostream& output);
  bool isImage() {return true;};
  bool addExtension(std::string newExt);
  void saveSummary(std::string baseName, TFile* myTargetFile);
private:
  std::unique_ptr<TCanvas> myCanvas_ = nullptr;
  int zoomedWidth_;
  int zoomedHeight_;
  std::string relativeHtmlDirectory_;
  std::string targetDirectory_;
  std::map<RootWImageSize, std::string> myText_;
  std::map<RootWImageSize, bool> fileSaved_;
  RootWImageSize lastSize_;
  std::string comment_;
  std::string name_;
  RootWImageSize makeSizeCode(int sw, int sh, int lw, int lh);
  std::vector<std::string> fileTypeV_;
  void saveSummaryLoop(TPad* basePad, std::string baseName, TFile* myTargetFile);

  static constexpr double thumb_compression_ = 2.;
  std::string allowedExtensions_; // Will be initialized in the constructor
  void setDefaultExtensions();

  static int imageCounter_;
  static std::map<std::string, int> imageNameCounter_;
};

class RootWFile : public RootWItem {
protected:
  std::string fileName_; // The destination file name
  std::string description_;
  std::string targetDirectory_;
public:
  RootWFile() {fileName_="aFile.txt";};
  ~RootWFile() {};
  RootWFile(std::string newFileName) { setFileName(newFileName); setDescription(""); };
  RootWFile(std::string newFileName, std::string newDescription) {setFileName(newFileName); setDescription(newDescription); };
  void setFileName(std::string newFileName) { fileName_ = newFileName;};
  void setDescription(std::string newDescription) { description_=newDescription; };
  std::string getFileName() { return fileName_ ; };
  std::string getDescription() { return description_ ; };
  void setTargetDirectory(std::string newTargetDirectory) {targetDirectory_ = newTargetDirectory; };
  bool isFile() {return true;};
};

class RootWFileList : public RootWItem {
protected:
  std::list<std::string> fileNames_;
  std::string description_;
  std::string targetDirectory_;
public:
  RootWFileList() {};
  template<class I> RootWFileList(I begin, I end) { addFileNames(begin, end); }
  template<class I> RootWFileList(I begin, I end, const std::string& description) : description_(description) { addFileNames(begin, end); }
  void addFileName(std::string newFileName) { fileNames_.push_back(newFileName); }
  template<class I> void addFileNames(I begin, I end) { fileNames_.insert(fileNames_.end(), begin, end); }
  void setDescription(std::string newDescription) { description_ = newDescription; }
  const std::list<std::string>& getFileNames() const { return fileNames_; }
  std::string getDescription() { return description_; }
  void setTargetDirectory(std::string newTargetDirectory) {targetDirectory_ = newTargetDirectory; }
  bool isFile() {return true;}
};

class RootWTextFile : public RootWFile {
private:
  std::stringstream myText_;
public:
  RootWTextFile() {};
  ~RootWTextFile() {};
  RootWTextFile(std::string newFileName) {setFileName(newFileName); setDescription(""); };
  RootWTextFile(std::string newFileName, std::string newDescription) {setFileName(newFileName); setDescription(newDescription); };
  void addText(std::string newText) {myText_ << newText ; };
  void addText(std::stringstream& newText) {myText_ << newText.str() ; };
  std::ostream& dump(std::ostream& output);
};

class RootWBinaryFile : public RootWFile {
private:
  std::string originalFileName_;
  bool noCopy_ = false;
public:
  RootWBinaryFile() {};
  ~RootWBinaryFile() {};
  RootWBinaryFile(std::string newFileName) {setFileName(newFileName); setDescription(""); setOriginalFile(""); };
  RootWBinaryFile(std::string newFileName, std::string newDescription) {setFileName(newFileName); setDescription(newDescription); setOriginalFile("");};
  RootWBinaryFile(std::string newFileName, std::string newDescription, std::string newOriginalFile) {
    setFileName(newFileName); setDescription(newDescription); setOriginalFile(newOriginalFile);};
  void setOriginalFile(std::string newFile) {originalFileName_ = newFile ; };
  std::ostream& dump(std::ostream& output);
  void setNoCopy(bool newNoCopy) { noCopy_ = newNoCopy ; }
};

class RootWBinaryFileList : public RootWFileList {
private:
  std::list<std::string> originalFileNames_;
public:
  RootWBinaryFileList() {};
  ~RootWBinaryFileList() {};
  template<class I> RootWBinaryFileList(I begin, I end) : RootWFileList(begin, end) {}
  template<class I> RootWBinaryFileList(I begin, I end, std::string newDescription) : RootWFileList(begin, end, newDescription) {}
  template<class I, class J> RootWBinaryFileList(I beginDestNames, I endDestNames, std::string newDescription, J beginOrigNames, J endOrigNames) : RootWFileList(beginDestNames, endDestNames, newDescription) {
    setOriginalFiles(beginOrigNames, endOrigNames);
  }
  template<class I> void setOriginalFiles(I begin, I end) { originalFileNames_.insert(originalFileNames_.end(), begin, end); }
  std::ostream& dump(std::ostream& output);
};

class RootWPage;

class RootWContent {
public:
  ~RootWContent();
  RootWContent();
  RootWContent(std::string title, bool visible=true);
  void setTargetDirectory(std::string newTargetDirectory);
  void addParagraph(std::string parText) ;
  void setTitle(std::string newTitle) ;
  void addItem(RootWItem* newItem);
  void addItem(std::unique_ptr<RootWItem> newItem);
  std::ostream& dump(std::ostream& output);
  RootWText& addText();
  RootWText& addText(std::string newText);
  RootWInfo& addInfo();
  RootWInfo& addInfo(std::string description);
  RootWInfo& addInfo(std::string description, std::string value);
  RootWTable& addTable();
  RootWImage& addImage();
  RootWImage& addImage(std::unique_ptr<TCanvas> myCanvas, int witdh, int height);
  RootWImage& addImage(std::unique_ptr<TCanvas> myCanvas, int witdh, int height, std::string relativeHtmlDirectory); // TODO: is this used for real?
  RootWTextFile& addTextFile();
  RootWTextFile& addTextFile(std::string newFileName);
  RootWTextFile& addTextFile(std::string newFileName, std::string newDescription);
  RootWBinaryFile& addBinaryFile();
  RootWBinaryFile& addBinaryFile(std::string newFileName);
  RootWBinaryFile& addBinaryFile(std::string newFileName, std::string newDescription);
  RootWBinaryFile& addBinaryFile(std::string newFileName, std::string newDescription, std::string newOriginalFile);
  void setPage(RootWPage*);
  TFile* getSummaryFile();
private:
  bool visible_;
  std::string title_;
  RootWPage* page_;
  //void setDefaultParameters();
  std::vector<RootWItem*> itemList_;
  std::string targetDirectory_;
};

class RootWPage;

class RootWSite {
protected:
private:
  std::vector<RootWPage*> pageList_;
  std::string title_;
  std::string comment_;
  std::string commentLink_;
  std::string toolkitName_;
  std::string toolkitGithub_;
  std::string toolkitContributors_;
  std::string revision_;
  std::string targetDirectory_;
  std::unique_ptr<TFile> summaryFile_;
  bool createSummaryFile_;
  std::string summaryFileName_;
public:
  ~RootWSite();
  RootWSite();
  RootWSite(std::string title);
  RootWSite(std::string title, std::string comment);
  void setTitle(std::string newTitle);
  void setComment(std::string newComment);
  void setCommentLink(std::string newCommentLink);
  std::string getTitle();
  std::string getComment();
  std::string getCommentLink();
  std::string getRevision();
  void setRevision (std::string newRevision);
  std::ostream& dumpHeader(std::ostream& output, RootWPage* thisPage);
  std::ostream& dumpFooter(std::ostream& output);
  RootWPage& addPage(std::string title, int relevance = least_relevant);
  void addPage(RootWPage* newPage, int relevance = least_relevant);
  void setTargetDirectory(std::string newTargetDirectory) {targetDirectory_ = newTargetDirectory; };
  bool prepareTargetDirectory();
  bool makeSite(bool verbose);
  void setSummaryFile(bool);
  TFile* getSummaryFile();
  void setSummaryFileName(std::string);
};

class RootWPage {
private:
  std::string title_;
  std::string address_;
  std::vector<RootWContent*> contentList_;
  RootWSite* site_;
  std::string targetDirectory_;
  double alert_;
  int relevance;
public:
  ~RootWPage();
  RootWPage();
  RootWPage(std::string title);
  void setTargetDirectory(std::string newTargetDirectory);
  void setTitle(std::string newTitle);
  std::string getTitle();
  void setAddress(std::string newAddress);
  std::string getAddress();
  void setSite(RootWSite* newSite);
  std::ostream& dump(std::ostream& output);
  void addContent(RootWContent* newContent);
  RootWContent& addContent(std::string title, bool visible=true);
  void addContent(std::unique_ptr<RootWContent> newContent);
  void setAlert(double alert);
  double getAlert();
  void setRelevance(int newRelevance);
  int getRelevance();
  TFile* getSummaryFile();
};

class RootWItemCollection {
private:
  std::map<std::string, RootWItem*> itemCollection_;
public:
  RootWItemCollection() {};
  ~RootWItemCollection() {};
  RootWItem* getItem(std::string itemName);
  void  addItem(RootWItem* anItem, std::string itemName);
  std::vector<RootWItem*> getOtherItems();
};

class RootWGraphViz : public RootWTextFile {
public:
  RootWGraphViz(std::string a, std::string b) : RootWTextFile(a, b) {} ;
  std::ostream& dump(std::ostream& output);
};

#endif

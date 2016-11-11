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

using namespace std;

#define THUMBSMALLSIZE 120
#define DEFAULTPROGRAMNAME "tkLayout"
#define DEFAULTPROGRAMSITE "https://github.com/tkLayout/tkLayout"
// The following is a list of allowed file etensions for TCanvas::SaveAs
// It should be separated, start and end with '|'
#define DEFAULTALLOWEDEXTENSIONS "|C|png|gif|svg|root|eps|pdf|ps|"

namespace RootWeb {
  std::string cleanUpObjectName(const std::string&);
  static const int least_relevant = std::numeric_limits<int>::min();
  static const int most_relevant = std::numeric_limits<int>::max();
};

using namespace RootWeb;

class RootWItem {
public:
  virtual ~RootWItem() {};
  RootWItem() {taken=false;};
  virtual bool isTable() {return false;};
  virtual bool isImage() {return false;};
  virtual bool isText() {return false;};
  virtual bool isFile() {return false;};
  virtual ostream& dump(ostream& output) {return output;};
  bool taken;
};

class RootWText: public RootWItem {
public:
  RootWText() {myText_.clear();};
  RootWText(string newText) {myText_.str(newText); };
  ~RootWText() {};
  void addText(string newText) { myText_ << newText; };
  ostream& dump(ostream& output) {output << myText_.str(); return output;};
  bool isText() {return true;};
protected:
  stringstream myText_;
};

class RootWInfo: public RootWItem {
public:
  RootWInfo() {description_="missing_description";value_="missing_value";};
  RootWInfo(string description) {description_=description; value_="missing_value";};
  RootWInfo(string description, string value) {description_=description;value_=value;};
  ~RootWInfo() {};
  string setDescription(string newText) { description_ = newText; return newText; };
  string setValue(string newText);
  string setValue(int number);
  string setValue(double number, int precision);
  string setValueSci(double number, int precision);
  string appendValue(string);
  string addValueText(string newText) {value_+=newText; return value_;};
  ostream& dump(ostream& output);
protected:
  string description_;
  string value_;
};

typedef std::map<pair<int,int>, string> rootWTableContent;
typedef std::map<pair<int,int>, int> rootWTableContentColor;
class RootWTable : public RootWItem {
public:
  ~RootWTable() {};
  RootWTable();
  void setContent(int row, int column, string content);
  void setContent(int row, int column, int number);
  void setContent(int row, int column, double number, int precision);
  void setContent(const rootWTableContent& newContent) { tableContent_ = newContent; };
  void setColor(int row, int column, int newColor);
  ostream& dump(ostream& output);
  pair<int, int> addContent(string content);
  pair<int, int> addContent(int number);
  pair<int, int> addContent(double number, int precision);
  pair<int, int> newLine();
  bool isTable() {return true;};
private:
  rootWTableContent tableContent_;
  rootWTableContentColor tableContentColor_;
  int serialRow_, serialCol_;
};

typedef string RootWImageSize;

class RootWImage : public RootWItem {
public:
  ~RootWImage();
  RootWImage();
  // TODO: the methods with TCanvas* (pointer) should be made obsolete
  RootWImage(TCanvas* myCanvas, int witdh, int height);
  RootWImage(TCanvas* myCanvas, int witdh, int height, string relativeHtmlDirectory); // TODO: is this used for real?
  RootWImage(TCanvas& myCanvas, int witdh, int height);
  RootWImage(TCanvas& myCanvas, int witdh, int height, string relativeHtmlDirectory); // TODO: is this used for real?
  void setCanvas(TCanvas* myCanvas);
  void setCanvas(TCanvas& myCanvas);
  void setComment(string newComment);
  void setName(string newName);
  std::string getName();
  void setZoomedSize(int witdh, int height);
  void setRelativeHtmlDirectory(string newDirectory);
  void setTargetDirectory(string newDirectory);
  string saveFiles(int smallWidth, int smallHeight);
  string saveFiles(int smallWidth, int smallHeight, int largeWidth, int largeHeight);
  ostream& dump(ostream& output);
  bool isImage() {return true;};
  bool addExtension(string newExt);
  void saveSummary(std::string baseName, TFile* myTargetFile);
private:
  TCanvas* myCanvas_;
  int zoomedWidth_;
  int zoomedHeight_;
  string relativeHtmlDirectory_;
  string targetDirectory_;
  map<RootWImageSize, string> myText_;
  map<RootWImageSize, bool> fileSaved_;
  RootWImageSize lastSize_;
  string comment_;
  string name_;
  RootWImageSize makeSizeCode(int sw, int sh, int lw, int lh);
  vector<string> fileTypeV_;
  void saveSummaryLoop(TPad* basePad, std::string baseName, TFile* myTargetFile);

  static constexpr double thumb_compression_ = 2.;
  string allowedExtensions_; // Will be initialized in the constructor
  void setDefaultExtensions();

  static int imageCounter_;
  static std::map<std::string, int> imageNameCounter_;
};

class RootWFile : public RootWItem {
protected:
  string fileName_; // The destination file name
  string description_;
  string targetDirectory_;
public:
  RootWFile() {fileName_="aFile.txt";};
  ~RootWFile() {};
  RootWFile(string newFileName) { setFileName(newFileName); setDescription(""); };
  RootWFile(string newFileName, string newDescription) {setFileName(newFileName); setDescription(newDescription); };
  void setFileName(string newFileName) { fileName_ = newFileName;};
  void setDescription(string newDescription) { description_=newDescription; };
  string getFileName() { return fileName_ ; };
  string getDescription() { return description_ ; };
  void setTargetDirectory(string newTargetDirectory) {targetDirectory_ = newTargetDirectory; };
  bool isFile() {return true;};
};

class RootWFileList : public RootWItem {
protected:
  std::list<string> fileNames_;
  string description_;
  string targetDirectory_;
public:
  RootWFileList() {};
  template<class I> RootWFileList(I begin, I end) { addFileNames(begin, end); }
  template<class I> RootWFileList(I begin, I end, const string& description) : description_(description) { addFileNames(begin, end); }
  void addFileName(string newFileName) { fileNames_.push_back(newFileName); }
  template<class I> void addFileNames(I begin, I end) { fileNames_.insert(fileNames_.end(), begin, end); }
  void setDescription(string newDescription) { description_ = newDescription; }
  const std::list<string>& getFileNames() const { return fileNames_; }
  string getDescription() { return description_; }
  void setTargetDirectory(string newTargetDirectory) {targetDirectory_ = newTargetDirectory; }
  bool isFile() {return true;}
};

class RootWTextFile : public RootWFile {
private:
  stringstream myText_;
public:
  RootWTextFile() {};
  ~RootWTextFile() {};
  RootWTextFile(string newFileName) {setFileName(newFileName); setDescription(""); };
  RootWTextFile(string newFileName, string newDescription) {setFileName(newFileName); setDescription(newDescription); };
  void addText(string newText) {myText_ << newText ; };
  void addText(stringstream& newText) {myText_ << newText.str() ; };
  ostream& dump(ostream& output);
};

class RootWBinaryFile : public RootWFile {
private:
  string originalFileName_;
  bool noCopy_ = false;
public:
  RootWBinaryFile() {};
  ~RootWBinaryFile() {};
  RootWBinaryFile(string newFileName) {setFileName(newFileName); setDescription(""); setOriginalFile(""); };
  RootWBinaryFile(string newFileName, string newDescription) {setFileName(newFileName); setDescription(newDescription); setOriginalFile("");};
  RootWBinaryFile(string newFileName, string newDescription, string newOriginalFile) {
    setFileName(newFileName); setDescription(newDescription); setOriginalFile(newOriginalFile);};
  void setOriginalFile(string newFile) {originalFileName_ = newFile ; };
  ostream& dump(ostream& output);
  void setNoCopy(bool newNoCopy) { noCopy_ = newNoCopy ; }
};

class RootWBinaryFileList : public RootWFileList {
private:
  std::list<string> originalFileNames_;
public:
  RootWBinaryFileList() {};
  ~RootWBinaryFileList() {};
  template<class I> RootWBinaryFileList(I begin, I end) : RootWFileList(begin, end) {}
  template<class I> RootWBinaryFileList(I begin, I end, string newDescription) : RootWFileList(begin, end, newDescription) {}
  template<class I, class J> RootWBinaryFileList(I beginDestNames, I endDestNames, string newDescription, J beginOrigNames, J endOrigNames) : RootWFileList(beginDestNames, endDestNames, newDescription) {
    setOriginalFiles(beginOrigNames, endOrigNames);
  }
  template<class I> void setOriginalFiles(I begin, I end) { originalFileNames_.insert(originalFileNames_.end(), begin, end); }
  ostream& dump(ostream& output);
};

class RootWPage;

class RootWContent {
public:
  ~RootWContent();
  RootWContent();
  RootWContent(string title, bool visible=true);
  void setTargetDirectory(string newTargetDirectory);
  void addParagraph(string parText) ;
  void setTitle(string newTitle) ;
  void addItem(RootWItem* newItem);
  ostream& dump(ostream& output);
  RootWText& addText();
  RootWText& addText(string newText);
  RootWInfo& addInfo();
  RootWInfo& addInfo(string description);
  RootWInfo& addInfo(string description, string value);
  RootWTable& addTable();
  RootWImage& addImage();
  RootWImage& addImage(TCanvas* myCanvas, int witdh, int height);
  RootWImage& addImage(TCanvas* myCanvas, int witdh, int height, string relativeHtmlDirectory); // TODO: is this used for real?
  RootWImage& addImage(TCanvas& myCanvas, int witdh, int height);
  RootWImage& addImage(TCanvas& myCanvas, int witdh, int height, string relativeHtmlDirectory); // TODO: is this used for real?
  RootWTextFile& addTextFile();
  RootWTextFile& addTextFile(string newFileName);
  RootWTextFile& addTextFile(string newFileName, string newDescription);
  RootWBinaryFile& addBinaryFile();
  RootWBinaryFile& addBinaryFile(string newFileName);
  RootWBinaryFile& addBinaryFile(string newFileName, string newDescription);
  RootWBinaryFile& addBinaryFile(string newFileName, string newDescription, string newOriginalFile);
  void setPage(RootWPage*);
  TFile* getSummaryFile();
private:
  bool visible_;
  string title_;
  RootWPage* page_;
  //void setDefaultParameters();
  vector<RootWItem*> itemList_;
  string targetDirectory_;
};

class RootWPage;

class RootWSite {
protected:
private:
  vector<RootWPage*> pageList_;
  string title_;
  string comment_;
  string commentLink_;
  vector<string> authorList_;
  string programName_;
  string programSite_;
  string revision_;
  string targetDirectory_;
  TFile* summaryFile_;
  bool createSummaryFile_;
  string summaryFileName_;
public:
  ~RootWSite();
  RootWSite();
  RootWSite(string title);
  RootWSite(string title, string comment);
  void setTitle(string newTitle);
  void setComment(string newComment);
  void setCommentLink(string newCommentLink);
  string getTitle();
  string getComment();
  string getCommentLink();
  string getRevision();
  void setRevision (string newRevision);
  ostream& dumpHeader(ostream& output, RootWPage* thisPage);
  ostream& dumpFooter(ostream& output);
  RootWPage& addPage(string title, int relevance = least_relevant);
  void addPage(RootWPage* newPage, int relevance = least_relevant);
  void addAuthor(string newAuthor);
  void setTargetDirectory(string newTargetDirectory) {targetDirectory_ = newTargetDirectory; };
  //void setStyleDirectory(string newStyleDirectory) {styleDirectory_ = newStyleDirectory; } ;
  bool makeSite(bool verbose);
  void setSummaryFile(bool);
  TFile* getSummaryFile();
  void setSummaryFileName(std::string);
};

class RootWPage {
private:
  string title_;
  string address_;
  vector<RootWContent*> contentList_;
  RootWSite* site_;
  string targetDirectory_;
  double alert_;
  int relevance;
public:
  ~RootWPage();
  RootWPage();
  RootWPage(string title);
  void setTargetDirectory(string newTargetDirectory);
  void setTitle(string newTitle);
  string getTitle();
  void setAddress(string newAddress);
  string getAddress();
  void setSite(RootWSite* newSite);
  ostream& dump(ostream& output);
  void addContent(RootWContent* newContent);
  RootWContent& addContent(string title, bool visible=true);
  void setAlert(double alert);
  double getAlert();
  void setRelevance(int newRelevance);
  int getRelevance();
  TFile* getSummaryFile();
};

class RootWItemCollection {
private:
  map<string, RootWItem*> itemCollection_;
public:
  RootWItemCollection() {};
  ~RootWItemCollection() {};
  RootWItem* getItem(string itemName);
  void  addItem(RootWItem* anItem, string itemName);
  vector<RootWItem*> getOtherItems();
};

class RootWGraphViz : public RootWTextFile {
public:
  RootWGraphViz(string a, string b) : RootWTextFile(a, b) {} ;
  ostream& dump(ostream& output);
};




#endif

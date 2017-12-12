#include <rootweb.hh>
#include <iostream>
#include <sstream>
#include <fstream>

#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <boost/filesystem/operations.hpp>

//#define STYLEDIRECTORY "TKG_STYLE_DIR"
//#define TARGETDIRECTORY "TKG_TARGET_DIR"

using namespace std;
using namespace boost::filesystem;

#include <mainConfigHandler.h>

int main(int argc, char** argv) {
  RootWItemCollection myCollection;
  
  RootWTable* myTable = new RootWTable();
  RootWTable* myTable2 = new RootWTable();
  myCollection.addItem(myTable, "table1");
  myCollection.addItem(myTable2, "table2");
  
  for (int x=0; x<10; x++) {
    for(int y=0; y<10; y++) {
      if ((x+y)%3) {
	myTable->setContent(x, y, x);
	myTable2->addContent(x*100);
      } else {
	myTable2->addContent("");
      }
    }
    myTable2->newLine();
  }
  
  TH1F* pippo = new TH1F("pippo", "Test plot here...", 100, 0, 10);
  for (int i = 0; i < 1E4; i++ ) {
    pippo->Fill(int((i*i/3.241234235562433))%10);
  }
  TCanvas* myCanvas = new TCanvas();
  myCanvas->cd();
  pippo->Draw();
  
  // This is to test the binary file linking:
  RootWBinaryFile *myBinaryFile = new RootWBinaryFile("testFile.png", "A test binary file", "/tmp/testme.png");
  myCanvas->SaveAs("/tmp/testme.png");
  myCollection.addItem(myBinaryFile, "Abinaryfile");

  RootWImage* myPlot = new RootWImage(myCanvas, 1024, 600);
  myPlot->setComment("This is the best plot ever seen...");
  //myPlot->saveFiles(200, 200);

  RootWContent* myContent = new RootWContent("Look at this");
  RootWContent* myContent2 = new RootWContent("Pictures here");
  stringstream testOutput;
  testOutput << "Hello, Europe!" << endl;
  testOutput << "How are you doing?" << endl;

  RootWTextFile* myTextFile = new RootWTextFile("test.txt");
  myTextFile->setDescription("Just another test obect (a file)");
  myTextFile->addText(testOutput);
  

  myContent->addParagraph("Hello, world!");
  myContent->addItem(myTextFile);
  myContent->addParagraph("Can you read this?");
  
  vector<RootWItem*> itemV = myCollection.getOtherItems();
  for (vector<RootWItem*>::iterator it = itemV.begin(); it!=itemV.end(); ++it) {
    myContent->addItem(*it);
  }
  myContent2->addItem(myPlot);
  myContent2->addItem(myPlot);

  RootWSite* mySite = new RootWSite();
  mySite->setTitle("Standard tracker");
  mySite->setComment("layout summary");
  mySite->addAuthor("Nicoletta De Maio");
  mySite->addAuthor("Stefano Mersi");
    

  RootWPage* myPage = new RootWPage("First");
  myPage->setAddress("index.html");
  myPage->addContent(myContent2);
  myPage->addContent(myContent);
 
  RootWPage* myPage2 = new RootWPage("SeconD");
  myPage2->setAddress("two.html");
  myPage2->addContent(myContent2);

  mySite->addPage(myPage);
  mySite->addPage(myPage2);
  //myPage->dump(cout);
  //myPage2->dump(cout);

  string targetDirectory;

  /* if you want to use the environment variables: 
  styleDirectory = string(getenv(STYLEDIRECTORY));
  if (!exists(styleDirectory)) {
    cerr << "Warning: cannot find the style directory " << styleDirectory 
         << " defined by the environment variable " << STYLEDIRECTORY << endl
         << "Using current directory as style directory (probably you will see the site without formatting)" << endl;
    styleDirectory=".";
  }

  targetDirectory = string(getenv(TARGETDIRECTORY));
  if (targetDirectory=="") {
    cerr << "Warning: cannot find the definition of target directory "
         << "by the environment variable " << STYLEDIRECTORY << endl
         << "Using current directory as target directory" << endl;
    targetDirectory=".";
  }
  if (!boost::filesystem::exists( targetDirectory )) {
    if (!boost::filesystem::create_directory(targetDirectory)) {
      cerr << "Couldn't create directory " << targetDirectory << endl;
      return 0;
    }
  } */

  mainConfigHandler myConfigHandler;
  if (myConfigHandler.getConfiguration(targetDirectory)) {
    targetDirectory+="/trackerName"; 
    mySite->setTargetDirectory(targetDirectory);
    mySite->makeSite();
  }

  //myContent.dump(cout);
  //myContent2.dump(cout);


  return 0;
}


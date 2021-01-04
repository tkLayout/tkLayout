/*
 * RootWSite.cc
 *
 *  Created on: 19. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.cc to standalone file
 */

#include "RootWSite.h"

#include <iomanip>
#include <fstream>

#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>

#include "global_constants.h"
#include "MainConfigHandler.h"
#include "MessageLogger.h"
#include "RootWPage.h"

// Used namespaces
using std::ostream;
using std::string;
using std::vector;

// Defined constants
const string RootWSite::c_toolkit_name          = "tkLayout";
const string RootWSite::c_toolkit_github        = "https://github.com/tkLayout/tkLayout/tree/masterLite";
const string RootWSite::c_toolkit_developers    = "https://github.com/tkLayout/tkLayout/graphs/contributors";

//
// Default constructor
//
RootWSite::RootWSite() {
  m_title                = "Untitled";
  m_comment              = "";
  m_commentLink          = "";
  m_toolkitName          = c_toolkit_name;
  m_toolkitGithub        = c_toolkit_github;
  m_toolkitDevelopers    = c_toolkit_developers;
  m_revision             = "";
  m_targetDirectory      = ".";
  //styleDirectory_ = ".";
}

//
// Constructor defining web-site title
//
RootWSite::RootWSite(string title) {
  m_title                = title;
  m_comment              = "";
  m_commentLink          = "";
  m_toolkitName          = c_toolkit_name;
  m_toolkitGithub        = c_toolkit_github;
  m_toolkitDevelopers    = c_toolkit_developers;
  m_revision             = "";
  m_targetDirectory      = ".";
  //styleDirectory_ = ".";
}

//
// Constructor defining web-site title and comment
//
RootWSite::RootWSite(string title, string comment) {
  m_title                = title;
  m_comment              = comment;
  m_commentLink          = "";
  m_toolkitName          = c_toolkit_name;
  m_toolkitGithub        = c_toolkit_github;
  m_toolkitDevelopers    = c_toolkit_developers;
  m_revision             = "";
  m_targetDirectory      = ".";
  //styleDirectory_ = ".";
}

//
// Default destructor
//
RootWSite::~RootWSite() {
}

//
// Create new web page & return reference to it
//
RootWPage& RootWSite::addPage(string newTitle, int relevance /* = c_leastRelevant */ ) {

  // Create pointer & get raw pointer to be returned
  std::unique_ptr<RootWPage> newPage(new RootWPage(newTitle, relevance, *this));
  auto newPagePtr = newPage.get();

  if (relevance == c_leastRelevant) {
    m_pageList.push_back(std::move(newPage));
  }
  else {
    // Go through the vector to find the first element which has a
    // lower relevance than 'relevance'
    auto beforeThis = m_pageList.end();

    for (auto itPage=m_pageList.begin(); itPage!=m_pageList.end(); ++itPage) {

      const RootWPage& aPage = **itPage;
      if (aPage.getRelevance()<relevance) {
        beforeThis = itPage;
        break;
      }
    }
    m_pageList.insert(beforeThis, std::move(newPage));
  }

  return (*newPagePtr);
}

//
// Main method creating the web-site (calling internally: dumpHeader(), dumpFooter() & dump() for each web page)
//
bool RootWSite::makeSite(bool verbose) {

  std::ofstream myPageFile;
  string        myPageFileName;
  //string targetStyleDirectory = m_targetDirectory + "/style";

  // Check if the directory already exists
  std::ostringstream message;
  if (boost::filesystem::exists( m_targetDirectory )) {
    if (! boost::filesystem::is_directory(m_targetDirectory) ) {

      message << "A file named " << m_targetDirectory << " already exists. Cannot create target directory";
      logERROR(message.str());
      return false;
    }
  } else {
    // If the directory does not exist, then try to create it
    if (!boost::filesystem::create_directory(m_targetDirectory)) {

      message << "Couldn't create directory " << m_targetDirectory;
      logERROR(message.str());
      return false;
    }
  }

  // Recreate the style symlink
  //if (boost::filesystem::exists( targetStyleDirectory )) {
  //boost::filesystem::remove_all( targetStyleDirectory );
  //}
  //boost::filesystem::create_symlink(styleDirectory_, targetStyleDirectory);

  for (auto itPage=m_pageList.begin(); itPage!=m_pageList.end(); itPage++) {

    RootWPage& myPage = **itPage;

    // Print name
    if (verbose) {
      if (itPage==m_pageList.begin()) std::cout << std::endl;
      std::cout << " " << myPage.getTitle() << std::endl;
    }

    myPageFileName = m_targetDirectory+"/"+myPage.getAddress();
    myPageFile.open(myPageFileName.c_str(), std::ios::out);
    myPage.setTargetDirectory(m_targetDirectory);
    myPage.dump(myPageFile);
    myPageFile.close();
  }

  return true;
}

//
// Dump method - printing header in the html format using the general web-site info
//
ostream& RootWSite::dumpHeader(ostream& output, const RootWPage& thisPage) const {

  time_t rawtime;
  time ( &rawtime );
  char timeBuffer [80];
  strftime (timeBuffer,80,"%a, %b %d %G %T",gmtime(&rawtime));

  output << "<html xmlns=\"http://www.w3.org/1999/xhtml\">" << std::endl
         << " <head>" << std::endl
         << "  <meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />" << std::endl
         << "  <title>"<< m_title << " - " << thisPage.getTitle() <<"</title>" << std::endl
         << "  <meta http-equiv=\"cache-control\" content=\"max-age=0\" />" << std::endl
         << "  <meta http-equiv=\"cache-control\" content=\"no-cache\" />" << std::endl
         << "  <meta http-equiv=\"expires\" content=\"0\" />" << std::endl
         << "  <meta http-equiv=\"expires\" content=\""<< timeBuffer<< " GMT\" />" << std::endl
         << "  <meta http-equiv=\"pragma\" content=\"no-cache\" />" << std::endl
         << "  <meta name=\"keywords\" content=\"CERN tracker design\" />" << std::endl
         << "  <meta name=\"description\" content=\" Tracker design summary page\" />" << std::endl
         << "  <link href=\"../style/default.css\" rel=\"stylesheet\" type=\"text/css\" />" << std::endl
         << "  <link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"../style/images/favicon.ico\">" << std::endl
         << " </head>" << std::endl
         << " <script type=\"text/javascript\" src=\"../style/functions.js\"></script>" << std::endl
         << " <body onload=\"preparePageGoodies();\">" << std::endl
         << "  <div id=\"outer\">" << std::endl
         << "   <div id=\"header\">" << std::endl
         << "    <h1>" << std::endl
         << "     <a href=\"index.html\">" << m_title << "</a>" << std::endl
         << "    </h1>" << std::endl;

  if (m_commentLink=="") {
    output << "    <h2>"<<m_comment<<"</h2>" << std::endl;
  }
  else {
    output << "    <h2>"
           << "<a href=\"" << m_commentLink << "\">"
           << m_comment
           << "</a>"
           << "</h2>" << std::endl;
  }
  output << "   </div>" << std::endl
         << "   <div id=\"menu\">" << std::endl
         << "    <ul>" << std::endl;

  for (auto& iPage : m_pageList) {

    output << "     <li";
    if (iPage.get()==&thisPage) output << " class=\"current\"";
    output << "> ";
    output << "<a href=\"" << iPage->getAddress() << "\"";
    if (iPage->getAlert()>0) {
      std::ostringstream anAlarmColor;
      anAlarmColor.str("");
      anAlarmColor << "rgb(255, "
                   << int(255 * iPage->getAlert())
                   << ", 0)";
      output << " style = 'background:"
             << anAlarmColor.str()
             << ";' ";
    }
    output << ">";
    output << iPage->getTitle();
    output <<"</a>" << " </li>" << std::endl;
  }
  output << "    </ul>" << std::endl
         << "   </div>" << std::endl
         << "   <div id=\"content\">" << std::endl
         << "    <div id=\"primaryContentContainer\">" << std::endl
         << "     <div id=\"primaryContent\">" << std::endl;

  return output;
}

//
//! Dump method - printing footer in the html format using the general web-site info
//
ostream& RootWSite::dumpFooter(ostream& output) const {
  output << "     </div>" << std::endl
         << "    </div>" << std::endl
         << "    <div class=\"clear\"></div>" << std::endl
         << "   </div>" << std::endl
         << "   <div id=\"footer\">" << std::endl;

  output << "<p> &copy; tkLayout developers: " 
	 << "<a href=\"" << m_toolkitDevelopers << "\">" << m_toolkitDevelopers << "</a>"
	 << "</p>" 
	 << std::endl;

  time_t rawtime;
  time ( &rawtime );
  char timeBuffer [80];
  strftime (timeBuffer,80,"%a, %b %d %G, %T",gmtime(&rawtime));

  output << "    <p>Page created on "<< timeBuffer << " GMT</p>" << std::endl
         << "    <p>By <a href=\""<< m_toolkitGithub <<"\">"<<m_toolkitName<<"</a>";

  if (m_revision!="") output << ", revision " << m_revision;

  output << "    </p>" << std::endl
         << "   </div>" << std::endl
         << "  </div>" << std::endl
         << " </body>" << std::endl
         << "</html>" << std::endl;

  return output;
}

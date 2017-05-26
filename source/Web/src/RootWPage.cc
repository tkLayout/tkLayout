/*
 * RootWPage.cc
 *
 *  Created on: 19. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.cc to standalone file
 */
#include "RootWPage.h"

#include "RootWContent.h"
#include "RootWSite.h"

//
// Default constructor defining name of web page, its priority when displayed on the web-site & reference to the main web-site
//
RootWPage::RootWPage(std::string title, int relevance, const RootWSite& webSite) :
 m_title(title),
 m_site(webSite),
 m_targetDirectory(webSite.getTargetDirectory()),
 m_alert(0),
 m_address(std::string(title+".html")),
 m_relevance(relevance)
{}

//
// Default destructor
//
RootWPage::~RootWPage()
{}

//
// Add new content holding different items: text, image, info, table, link to text file or link to binary file,
// reference to the created content container is returned.
//
RootWContent& RootWPage::addContent(std::string title, bool visible)
{
  std::unique_ptr<RootWContent> newContent(new RootWContent(title, visible));
  m_contentList.push_back(std::move(newContent));
  return (*(m_contentList.back().get()));
}

//
// Dump method - printing tab in the html format using the general web-site info
//
std::ostream& RootWPage::dump(std::ostream& output) {

  // Header standard part, with the menu list, etc...
  // up to opening the primaryContent div
  m_site.dumpHeader(output, *this);

  for (auto& iContent : m_contentList) {

    iContent->setTargetDirectory(m_targetDirectory);
    iContent->dump(output);
  }

  // Header standard part, with copyright, etc...
  // starting from closing the primaryContent div
  m_site.dumpFooter(output);

  return output;

}

//
// Find certain content on the web page using title as an identifier
//
const RootWContent& RootWPage::findContent(std::string title) const {

  const RootWContent* content = nullptr;

  for (auto& iContent : m_contentList) {
    if (iContent->getTitle()==title) {
      content = iContent.get();
    }
  }

  // Check that RootWContent of given name exists
  if (content) return *content;
  else throw std::invalid_argument( "RootWPage::findContent - RootWContent not defined (null reference), check!!!" );
}

//
// Set web-page title and html file name if not defined
//
void RootWPage::setTitle(std::string newTitle) {
  m_title = newTitle;
  if (m_address=="") m_address = m_title+".html";
}

//
// Set alert
//
void RootWPage::setAlert(double alert) {
  if      (alert<0) m_alert=0;
  else if (alert>1) m_alert=1;
  else              m_alert=alert;
}

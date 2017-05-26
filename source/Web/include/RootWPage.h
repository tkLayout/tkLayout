/*
 * RootWPage.h
 *
 *  Created on: 19. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */

#ifndef INCLUDE_ROOTWPAGE_H_
#define INCLUDE_ROOTWPAGE_H_

#include <memory>
#include <vector>
#include <string>
#include <ostream>

// Forward declaration
class RootWContent;
class RootWSite;

/*
 * @class RootWPage
 * Web container that is printed out as a html tab on the main web-site (RootWSite). A user forms
 * its content by adding various blocks, so-called RootWContents - using addContent() method. The
 * RootWContent can hold different types: text, image, info, table, link to text file or link to
 * binary file.
 */
class RootWPage {

 public:

  //! Default constructor defining name of web page, its priority when displayed on the web-site & reference to the main web-site
  RootWPage(std::string title, int relevance, const RootWSite& webSite);
  //! Default destructor
  ~RootWPage();

  //! Add new content holding different items: text, image, info, table, link to text file or link to binary file,
  //! reference to the created content container is returned.
  RootWContent& addContent(std::string title, bool visible=true);

  //! Dump method - printing tab in the html format using the general web-site info
  std::ostream& dump(std::ostream& output);

  //! Find certain content on the web page using title as an identifier
  const RootWContent& findContent(std::string title) const;

  // Setter methods
  void setTargetDirectory(std::string newTargetDirectory) {m_targetDirectory = newTargetDirectory;}
  void setTitle(std::string newTitle);
  void setAddress(std::string newAddress)                 {m_address = newAddress;}
  void setAlert(double alert);
  void setRelevance(int newRelevance)                     {m_relevance = newRelevance;}

  // Getter methods
  std::string getTitle()   const {return m_title;}
  std::string getAddress() const {return m_address;}
  double getAlert()        const {return m_alert;}
  int getRelevance()       const {return m_relevance;}

 private:

  std::vector<std::unique_ptr<RootWContent>> m_contentList; //! List of assigned web-blocks (pointers owned by this class)

  std::string m_title;          //!< Web page name
  std::string m_address;        //!< Html file name - if not defined defined as title.html
  std::string m_targetDirectory;//!< Target directory - set the same as target directory of the main web-site
  double      m_alert;
  int         m_relevance;      //!< Web page priority according to which the web pages are displayed on the main web-site

  const RootWSite& m_site;

}; // Class

#endif /* INCLUDE_ROOTWPAGE_H_ */

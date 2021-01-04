/*
 * RootWSite.h
 *
 *  Created on: 19. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */

#ifndef INCLUDE_ROOTWSITE_H_
#define INCLUDE_ROOTWSITE_H_

#include <memory>
#include <vector>
#include <string>
#include <ostream>

// Forward declaration
class RootWPage;

/*
 * @class RootWSite
 * Web interface container (web-site) collecting all pieces of information (analysis results) and printing them
 * in the html format to the target directory. The web-site consists of individual web-pages (web tabs), a header
 * and footer. Individual web pages are added to the web-site container using addPage() method. The web-site is
 * printed out by calling makeSite(), which practically means an internal call of dumpHeader(), dumpFooter() and
 * dump() method for each web page.
 */
class RootWSite {

 public:

  //! Default constructor
  RootWSite();
  //! Constructor defining web-site title
  RootWSite(std::string title);
  //! Constructor defining web-site title and comment
  RootWSite(std::string title, std::string comment);
  //! Default destructor
  ~RootWSite();

  //! Create new web page & return reference to it
  RootWPage& addPage(std::string title, int relevance = c_leastRelevant);

  //! Main method creating the web-site (calling internally: dumpHeader(), dumpFooter() & dump() for each web page)
  bool makeSite(bool verbose);

  //! Dump method - printing header in the html format using the general web-site info
  std::ostream& dumpHeader(std::ostream& output, const RootWPage& thisPage) const;

  //! Dump method - printing footer in the html format using the general web-site info
  std::ostream& dumpFooter(std::ostream& output) const;

  // Setter methods
  inline void setTitle(std::string newTitle)             {m_title = newTitle;}
  inline void setComment(std::string newComment)         {m_comment = newComment;}
  inline void setCommentLink(std::string newCommentLink) {m_commentLink = newCommentLink; }
  inline void setRevision (std::string newRevision)      {m_revision = newRevision;}

  //! Set target directory, where all html generated results will be saved
  inline void setTargetDirectory(std::string newTargetDirectory) {m_targetDirectory = newTargetDirectory; };

  //void setStyleDirectory(std::string newStyleDirectory) {m_styleDirectory = newStyleDirectory; } ;

  // Getter methods
  inline std::string getTitle()          const {return m_title;}
  inline std::string getComment()        const {return m_comment;}
  inline std::string getCommentLink()    const {return m_commentLink;}
  inline std::string getRevision()       const {return m_revision;}
  inline std::string getTargetDirectory()const {return m_targetDirectory;}

 private:

  std::vector<std::unique_ptr<RootWPage>>  m_pageList; //!< List of all web pages assigned to the web site (pointers owned by this class)
  std::string              m_title;                    //!< Web site title
  std::string              m_comment;                  //!< Comment to the web-site name (e.g. All layouts)
  std::string              m_commentLink;              //!< Link connecting the comment with other web-site (e.g. ../ to get to the uppermost directory with all layouts)
  std::string              m_toolkitName;              //!< Program name, i.e. tkLayout
  std::string              m_toolkitGithub;            //!< Github address, where the program is located
  std::string              m_toolkitDevelopers;        //!< Core toolkit developers
  std::string              m_revision;                 //!< Footer prints out current software revision
  std::string              m_targetDirectory;          //!< Directory, where all the html results are saved
  std::string              m_styleDirectory;

  static const int         c_leastRelevant= -1000;     //!< Lowest priority to order web pages on the web site
  static const std::string c_toolkit_name;             //!< Default program name, i.e. tkLayout
  static const std::string c_toolkit_github;           //!< Default github address, where the program is located
  static const std::string c_toolkit_developers;       //!< Core toolkit developers

}; // Class

#endif /* INCLUDE_ROOTWSITE_H_ */

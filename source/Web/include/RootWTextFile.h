/*
 * RootWTextFile.h
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */
#ifndef INCLUDE_ROOTWTEXTFILE_H_
#define INCLUDE_ROOTWTEXTFILE_H_

#include <string>
#include <sstream>
#include <ostream>

#include "RootWFile.h"

/*
 * @class RootWTextFile
 * Create a text, flush its content to a text file, copy the text file to the target directory &
 * create a file reference on the web page. All done by internall dump() method.
 */
class RootWTextFile : public RootWFile {

 public:

  //! Default constructor setting the text file name
  RootWTextFile(std::string newFileName);
  //! Default constructor setting the text file name & file description on the web site
  RootWTextFile(std::string newFileName, std::string newDescription);
  //! Default destructor
  ~RootWTextFile();

  //! Add text to the text file
  void addText(std::string newText);
  void addText(std::stringstream& newText);

  //! Dump method - create text, print its content to text file & create a link on the web page to this file
  std::ostream& dump(std::ostream& output);

 private:

  std::stringstream m_text;

}; // Class

#endif /* SOURCE_WEB_INCLUDE_ROOTWTEXTFILE_H_ */

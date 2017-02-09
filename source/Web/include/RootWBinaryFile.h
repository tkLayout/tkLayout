/*
 * RootWBinaryFile.h
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */

#ifndef INCLUDE_ROOTWBINARYFILE_H_
#define INCLUDE_ROOTWBINARYFILE_H_

#define BOOST_NO_CXX11_SCOPED_ENUMS
#define BOOST_NO_CXX11_SCOPED_ENUM

#include <ostream>
#include <string>

#include "RootWFile.h"

/*
 * @class RootWBinaryFile
 * Copy a defined binary file to a target directory (if not dontCopyFile variable specified)
 * and then print a reference to the file with some text description (in html format) on
 * the required web page. Use dump() method to make it.
 */
class RootWBinaryFile : public RootWFile {

 public:

  //! Default constructor setting the binary file name
  RootWBinaryFile(std::string newFileName);
  //! Default constructor setting the binarry file name & file description on the web site
  RootWBinaryFile(std::string newFileName, std::string newDescription);
  //! Default constructor setting the binarry file name, file description on the web site & original file from which the binary is cloned
  RootWBinaryFile(std::string newFileName, std::string newDescription, std::string newOriginalFile);
  //! Default destructor
  ~RootWBinaryFile();

  void setOriginalFile(std::string newFile) {m_originalFileName = newFile ; };
  void setDontCopyFile(bool noCopy) { m_dontCopyFile = noCopy;}

  //! Dump method - create copy of original binary file & create a link on the web page to this file
  std::ostream& dump(std::ostream& output);

 private:
  std::string m_originalFileName;
  bool        m_dontCopyFile = false; //! Don't copy the given file to the target directory? (Implicitely copied)

}; // Class

#endif /* INCLUDE_ROOTWBINARYFILE_H_ */

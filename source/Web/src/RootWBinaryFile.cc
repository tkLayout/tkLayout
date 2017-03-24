/*
 * RootWBinaryFile.cc
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.cc to standalone file
 */
#include "RootWBinaryFile.h"

#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>

#include "MessageLogger.h"

//
// Default constructor setting the binary file name
//
RootWBinaryFile::RootWBinaryFile(std::string newFileName) :
 RootWFile(newFileName)
{}

//
// Default constructor setting the binarry file name & file description on the web site
//
RootWBinaryFile::RootWBinaryFile(std::string newFileName, std::string newDescription) :
 RootWFile(newFileName, newDescription)
{}

//
// Default constructor setting the binarry file name, file description on the web site & original file from which the binary is cloned
//
RootWBinaryFile::RootWBinaryFile(std::string newFileName, std::string newDescription, std::string newOriginalFile) :
 RootWFile(newFileName, newDescription),
 m_originalFileName(newOriginalFile)
{}

//
// Default destructor
//
RootWBinaryFile::~RootWBinaryFile() {};

//
// Dump method - create copy of original binary file (if dontCopyFile not specified) & create a link on the web page to this file
//
ostream& RootWBinaryFile::dump(ostream& output) {

  if ((m_originalFileName=="") && (!m_dontCopyFile)) {
    logERROR("RootWBinaryFile::dump() was called without prior setting the original file name");
    return output;
  }
  if (m_fileName=="") {
    logERROR("RootWBinaryFile::dump() was called without prior setting the destination file name");
    return output;
  }

  string destinationFileName = m_targetDirectory +"/" + m_fileName;

  if ((!m_dontCopyFile)&&(boost::filesystem::exists(m_originalFileName) && m_originalFileName != destinationFileName)) { // CUIDADO: naive control on copy on itself. it only matches the strings, not taking into account relative paths and symlinks

    try {
      if (boost::filesystem::exists(destinationFileName)) boost::filesystem::remove(destinationFileName);
        boost::filesystem::copy_file(m_originalFileName, destinationFileName);
    }
    catch (boost::filesystem::filesystem_error& e) {
      std::ostringstream message;
      message << "RootWBinaryFile::dump problem: " << e.what() << endl;
      logERROR(message);
      return output;
    }
  }

  output << "<b>" << m_description << ":</b> <a href=\""
         << m_fileName << "\">"
         << m_fileName << "</a></tt><br/>";
  return output;
};


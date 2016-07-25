/*
 * RootWBinaryFileList.cc
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */
#include "RootWBinaryFileList.h"

#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>

#include "MessageLogger.h"

//
// Default constructor
//
RootWBinaryFileList::RootWBinaryFileList()
{}

//
// Default destructor
//
RootWBinaryFileList::~RootWBinaryFileList()
{}

//
// Dump method - create copy of original binary files & create links on the web page to these files
//
std::ostream& RootWBinaryFileList::dump(std::ostream& output) {

  if (m_originalFileNames.empty()) {
    logERROR("RootWBinaryFileList::dump() was called without prior setting the original file names");
    return output;
  }
  if (m_fileNames.empty()) {
    logERROR("RootWBinaryFileList::dump() was called without prior setting the destination file names");
    return output;
  }
  if (m_fileNames.size() != m_originalFileNames.size()) {
    logERROR("RootWBinaryFileList::dump() was called with original and destination file name lists differing in length");
    return output;
  }

  auto it = m_originalFileNames.begin();
  for (auto fn : m_fileNames) {
    //auto path = split(fn, "/");
    string destinationFileName = m_targetDirectory;
    //std::for_each(path.begin(), path.end()-1, [&destinationFileName](string p) {
    //  destinationFileName += "/" + p;
    //  if (!boost::filesystem::exists(destinationFileName)) boost::filesystem::create_directory(destinationFileName);
    //});
    destinationFileName += "/" + fn; //path.back();
    if (boost::filesystem::exists(*it) && *it != destinationFileName) { // CUIDADO: naive control on copy on itself.
      try {
        if (boost::filesystem::exists(destinationFileName)) boost::filesystem::remove(destinationFileName);
        boost::filesystem::copy_file(*it++, destinationFileName);
      }
      catch (boost::filesystem::filesystem_error& e) {
        std::ostringstream message;
        message << "RootWBinaryFileList::dump problem: " << e.what() << endl;
        logERROR(message);cerr << e.what() << endl;
        return output;
      }
    }
  }
  std::vector<std::string> cleanedUpFileNames;
  std::transform(m_originalFileNames.begin(), m_originalFileNames.end(), std::back_inserter(cleanedUpFileNames), [](const std::string& s) {
    auto pos = s.find("stdinclude");
    return pos != string::npos ? s.substr(pos) : s;
  });
  output << "<b>" << m_description << ":</b>";
  auto dfn = m_fileNames.begin();
  for (auto cfn : cleanedUpFileNames) {
    output << (cleanedUpFileNames.size() > 1 ? "<br>&nbsp&nbsp&nbsp&nbsp" : "") << " <a href=\"" << *(dfn++) << "\">" << cfn << "</a></tt>" << " ";
  }
  output << "<br/>";
  return output;
}

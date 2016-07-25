/*
 * RootWTextFile.cc
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.cc to standalone file
 */
#include "RootWTextFile.h"

#include <fstream>

#include "MessageLogger.h"

//
// Default constructor setting the text file name
//
RootWTextFile::RootWTextFile(std::string newFileName) :
 RootWFile(newFileName)
{}

//
// Default constructor setting the text file name & file description on the web site
//
RootWTextFile::RootWTextFile(std::string newFileName, std::string newDescription) :
 RootWFile(newFileName, newDescription)
{}

//
// Default destructor
//
RootWTextFile::~RootWTextFile()
{}

//
// Add text to the text file
//
void RootWTextFile::addText(std::string newText)
{
  m_text << newText ;
}

void RootWTextFile::addText(std::stringstream& newText)
{
  m_text << newText.str();
}

//
// Dump method - create text, print its content to text file & create a link on the web page to this file
//
std::ostream& RootWTextFile::dump(std::ostream& output) {

  if (m_fileName=="") {
    logERROR("RootWTextFile::dump() was called without prior setting the destination file name");
    return output;
  }

  std::ofstream outputFile;
  string destinationFileName = m_targetDirectory +"/" + m_fileName;

  outputFile.open(destinationFileName.c_str());
  if (outputFile.is_open()) {

    outputFile << m_text.str() << endl;
    outputFile.close();

    output << "<b>" << m_description << ":</b> <a href=\""
           << m_fileName << "\">"
           << m_fileName << "</a></tt><br/>";
  }
  else {
    logERROR("RootWTextFile::dump() - text file not correcty opened!!!");
  }

  return output;
};

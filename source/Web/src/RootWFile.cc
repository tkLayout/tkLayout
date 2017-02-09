/*
 * RootWFile.cc
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */
#include "RootWFile.h"

//
// Default constructor setting the file name
//
RootWFile::RootWFile(std::string newFileName)
{
  m_fileName    = newFileName;
  m_description = "";
}

//
// Default constructor setting the file name & file description on the web site
//
RootWFile::RootWFile(std::string newFileName, std::string newDescription)
{
  m_fileName    = newFileName;
  m_description = newDescription;
}

//
// Default destructor
//
RootWFile::~RootWFile()
{}

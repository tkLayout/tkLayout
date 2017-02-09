/*
 * RootWFileList.cc
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.cc to standalone file
 */
#include "RootWFileList.h"

//
// Default constructor
//
RootWFileList::RootWFileList()
{}

//
// Add file to the list
//
void RootWFileList::addFileName(std::string newFileName)
{
  m_fileNames.push_back(newFileName);
}


/*
 * RootWFileList.h
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */

#ifndef INCLUDE_ROOTWFILELIST_H_
#define INCLUDE_ROOTWFILELIST_H_

#include <list>
#include <string>

#include "RootWItem.h"

/*
 * @class RootWFileList
 * An abstract class for list of files (binary, text) to be linked on the web page.
 */
class RootWFileList : public RootWItem {

 public:

  //! Default constructor
  RootWFileList();
  template<class I> RootWFileList(I begin, I end) {addFileNames(begin, end);}
  template<class I> RootWFileList(I begin, I end, std::string description) : m_description(description) {addFileNames(begin, end);}

  //! Add file to the list
  void addFileName(std::string newFileName);

  //! Add files to the list using iterators
  template<class I> void addFileNames(I begin, I end) {m_fileNames.insert(m_fileNames.end(), begin, end);}

  // Setter methods
  void setDescription(std::string newDescription)         {m_description = newDescription; }
  void setTargetDirectory(std::string newTargetDirectory) {m_targetDirectory = newTargetDirectory; }

  // Getter methods
  const std::list<std::string>& getFileNames() const { return m_fileNames; }
  std::string getDescription()                 const { return m_description; }

 protected:

  std::list<std::string> m_fileNames;
  std::string            m_description;
  std::string            m_targetDirectory;

}; // Class

#endif /* INCLUDE_ROOTWFILELIST_H_ */

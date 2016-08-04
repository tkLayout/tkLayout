/*
 * RootWFile.h
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */
#ifndef INCLUDE_ROOTWFILE_H_
#define INCLUDE_ROOTWFILE_H_

#include <string>

#include "RootWItem.h"

/*
 * @class RootWFile
 * An abstract class for file content (binary, text) to be linked on the web page.
 */
class RootWFile : public RootWItem {

 public:

  //! Default constructor setting the file name
  RootWFile(std::string newFileName);
  //! Default constructor setting the file name & file description on the web site
  RootWFile(std::string newFileName, std::string newDescription);
  //! Default destructor
  ~RootWFile();

  // Setter methods
  void setFileName(std::string newFileName)               { m_fileName = newFileName;};
  void setDescription(std::string newDescription)         { m_description=newDescription; };
  void setTargetDirectory(std::string newTargetDirectory) { m_targetDirectory = newTargetDirectory; };

  // Getter methods
  std::string getFileName()    const { return m_fileName; };
  std::string getDescription() const { return m_description; };

 protected:
   std::string m_fileName; // The destination file name
   std::string m_description;
   std::string m_targetDirectory;

}; // Class

#endif /* SOURCE_WEB_INCLUDE_ROOTWFILE_H_ */

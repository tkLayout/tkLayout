/*
 * RootWBinaryFileList.h
 *
 *  Created on: 22. 7. 2016
 *      Author: Drasal (CERN) - Extracted from rootweb.h to standalone file
 */
#ifndef INCLUDE_ROOTWBINARYFILELIST_H_
#define INCLUDE_ROOTWBINARYFILELIST_H_

#define BOOST_NO_CXX11_SCOPED_ENUMS
#define BOOST_NO_CXX11_SCOPED_ENUM

#include <ostream>
#include <string>

#include "RootWFileList.h"

/*
 * @class RootWBinaryFileList
 * Copy list of binary files to a target directory & print their names and description in
 * the html format on the web page using the dump() method. Don't forget to set the list of
 * original binary files to be copied to the target directory and the names of copied files to
 * be referenced on the web page.
 */
class RootWBinaryFileList : public RootWFileList {

 public:

  //! Default constructor
  RootWBinaryFileList();
  //! Constructor defining list of file names using iterators
  template<class I> RootWBinaryFileList(I begin, I end) : RootWFileList(begin, end) {};
  template<class I> RootWBinaryFileList(I begin, I end, std::string newDescription) : RootWFileList(begin, end, newDescription) {};
  template<class I, class J> RootWBinaryFileList(I beginDestNames, I endDestNames, std::string newDescription, J beginOrigNames, J endOrigNames) :
   RootWFileList(beginDestNames, endDestNames, newDescription) {setOriginalFiles(beginOrigNames, endOrigNames);};
  // Default destructor
  ~RootWBinaryFileList();

  //! Define names of original binary files, to be copied to the target directory & linked on the web page
  template<class I> void setOriginalFiles(I begin, I end) {m_originalFileNames.insert(m_originalFileNames.end(), begin, end);};

  //! Dump method - create copy of original binary files & create links on the web page to these files
  std::ostream& dump(std::ostream& output);

 private:
  std::list<std::string> m_originalFileNames;
};

#endif /* INCLUDE_ROOTWBINARYFILELIST_H_ */

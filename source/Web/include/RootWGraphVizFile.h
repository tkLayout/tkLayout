/*
 * RootWGraphVizFile.h
 *
 *  Created on: 5. 10. 2016
 */
#ifndef INCLUDE_ROOTWGRAPHVIZFILE_H_
#define INCLUDE_ROOTWGRAPHVIZFILE_H_

#include <string>
#include <ostream>

#include "RootWTextFile.h"

/*
 * @class RootWGraphViz
 * Extension to the RootWTextFile class, writing out the description file for graphviz dot program
 * into the text file, which is linked then in the target directory. If graphviz dot program exists
 * on the system, it will also call this program and directly create the output: graph. This graph
 * show complex relations as decribed by the input file.
 */
class RootWGraphVizFile : public RootWTextFile {

public:

  //! Constructor calling RootWTextFile constructor
  RootWGraphVizFile(std::string fileName, std::string description);

  //! Dump method - main method building graph based on the given file and linking it on the web
  std::ostream& dump(std::ostream& output);
};

#endif /* INCLUDE_ROOTWGRAPHVIZFILE_H_ */

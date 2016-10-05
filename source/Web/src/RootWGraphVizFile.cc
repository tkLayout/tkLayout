/*
 * RootWGraphVizFile.cc
 *
 *  Created on: 5. 10. 2016
 */
#include "RootWGraphVizFile.h"

#include "MessageLogger.h"
#include "RootWBinaryFile.h"

RootWGraphVizFile::RootWGraphVizFile(std::string fileName, std::string description) :
  RootWTextFile(fileName, description)
{}

std::ostream& RootWGraphVizFile::dump(std::ostream& output) {

  std::string myDescription = m_description;
  m_description += " (GraphViz)";

  RootWTextFile::dump(output);

  // Check if graphviz dot program exists on the system
  int result = system("which dot > /dev/null");
  if (result==0) {

    std::string svgFileName = m_fileName+".svg";
    std::string command = "dot -Tsvg " + m_targetDirectory +"/" + m_fileName + " > " + m_targetDirectory + "/" + svgFileName;

    // Call graphviz dot program to create a graph (using graph components description given by the provided text file)
    result = system(command.c_str());

    if (result==0) {
      RootWBinaryFile* mySvg = new RootWBinaryFile(svgFileName, myDescription+" (SVG)" );
      mySvg->setDontCopyFile(true);
      mySvg->dump(output);
    }
    else {
      logWARNING("RootWGraphVizFile couldn't produce the GraphViz plot, graphviz program 'dot' missing on the system");
    }
  }
  return output;
}


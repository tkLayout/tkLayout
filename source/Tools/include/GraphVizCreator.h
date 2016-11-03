#ifndef INCLUDE_GRAPHVIZCREATOR_H_
#define INCLUDE_GRAPHVIZCREATOR_H_

#include <map>
#include <string>
#include <vector>

/*
 * @class GraphVizCreator
 * @brief Using graphviz library one can illustrate a graph description file, created by this class & illustrating
 * recursive @include complexity & relations between individual @include files (used for geometry configuration files)
 */
class GraphVizCreator {

public:

  int         getFileId(std::string istreamid);
  std::string getNodeName(int nodeId);
  void        addGraphLink(int parentNodeId, int childNodeId);
  void        clearGraphLinks(int parentNodeId);
  std::string createGraphVizFile();
  void        addNodeRename(std::string nodeName, std::string symbolic);
  void        addNodeUrl(std::string nodeName, std::string url);
  void        setNodeLocal(std::string nodeName, bool isLocal);
  void        prepareNodeOutput(std::string absoluteFileName, std::string relativeFileName, bool webOutput);
  void        createEncodedFileList(std::vector<std::string>& originalFileNames, std::vector<std::string>& destinationFileNames);

private:

  std::string m_graphVizConnections;

  std::map<std::string, int>         m_graphVizFilesToNodes;
  std::map<std::string, std::string> m_nodeRenameMap;
  std::map<std::string, std::string> m_nodeDestinationFilename;
  std::map<std::string, std::string> m_nodeUrl;
  std::map<int, std::vector<int> >   m_nodeLinks;
  std::map<std::string, bool>        m_nodeLocal;

  std::string destinationEncode(const std::string& originalFileName);

  const std::string c_indent = "    ";
}; // Class

#endif /*INCLUDE_GRAPHVIZCREATOR_H_ */

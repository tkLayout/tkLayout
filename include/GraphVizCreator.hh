#ifndef GRAPHVIZCREATOR_H
#define GRAPHVIZCREATOR_H

#include <map>
#include <string>
#include <list>

class GraphVizCreator {
public:
  int getFileId(std::string istreamid);
  std::string getNodeName(int nodeId);
  void addGraphLink(int parentNodeId, int childNodeId);
  std::string createGraphVizFile();
  void addNodeRename(std::string nodeName, std::string symbolic);
  void addNodeUrl(std::string nodeName, std::string url);
private:
  std::map<std::string, int> graphVizFilesToNodes_;
  std::string graphVizConnections_;
  std::map<std::string, std::string> nodeRenameMap_;
  std::map<std::string, std::string> nodeUrl_;
  const std::string indent_ = "    ";
};

#endif

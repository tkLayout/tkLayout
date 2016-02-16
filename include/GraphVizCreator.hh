#ifndef GRAPHVIZCREATOR_H
#define GRAPHVIZCREATOR_H

#include <map>
#include <string>
#include <vector>

class GraphVizCreator {
public:
  int getFileId(std::string istreamid);
  std::string getNodeName(int nodeId);
  void addGraphLink(int parentNodeId, int childNodeId);
  void clearGraphLinks(int parentNodeId);
  std::string createGraphVizFile();
  void addNodeRename(std::string nodeName, std::string symbolic);
  void addNodeUrl(std::string nodeName, std::string url);
  void setNodeLocal(std::string nodeName, bool isLocal);
  void prepareNodeOutput(std::string absoluteFileName, std::string relativeFileName, bool webOutput);
  void createEncodedFileList(std::vector<std::string>& originalFileNames,
			     std::vector<std::string>& destinationFileNames);
private:
  std::map<std::string, int> graphVizFilesToNodes_;
  std::string graphVizConnections_;
  std::map<std::string, std::string> nodeRenameMap_;
  std::map<std::string, std::string> nodeDestinationFilename_;
  std::map<std::string, std::string> nodeUrl_;
  std::map<int, std::vector<int> > nodeLinks_;
  std::map<std::string, bool> nodeLocal_;
  std::string destinationEncode(const std::string& originalFileName);

  const std::string indent_ = "    ";
};

#endif

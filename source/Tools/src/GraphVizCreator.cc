#include <GraphVizCreator.h>

// Boost include files
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>

int GraphVizCreator::getFileId(std::string istreamid) {
  if (m_graphVizFilesToNodes[istreamid]==0) m_graphVizFilesToNodes[istreamid]=m_graphVizFilesToNodes.size();
  return m_graphVizFilesToNodes[istreamid];
}

std::string GraphVizCreator::getNodeName(int nodeId) {
  return "node_"+std::to_string(nodeId);
}

void GraphVizCreator::addGraphLink(int parentNodeId, int childNodeId) {
  std::vector<int>& linkedNodes = m_nodeLinks[parentNodeId];
  linkedNodes.push_back(childNodeId);
}

void GraphVizCreator::clearGraphLinks(int parentNodeId) {
  std::vector<int>& linkedNodes = m_nodeLinks[parentNodeId];
  linkedNodes.clear();
}

void GraphVizCreator::addNodeRename(std::string nodeName, std::string symbolic) {
  m_nodeRenameMap[nodeName]=symbolic;
}

void GraphVizCreator::setNodeLocal(std::string nodeName, bool isLocal) {
  m_nodeLocal[nodeName]=isLocal;
}

void GraphVizCreator::addNodeUrl(std::string nodeName, std::string url) {
  m_nodeUrl[nodeName]=url;
}

std::string GraphVizCreator::createGraphVizFile() {

  std::string result;
  result += "digraph aGraph {\n";
  result += c_indent + "overlap=false;\n";
  result += c_indent + "rankdir=LR;\n";
  result += c_indent + "graph [ranksep=2, nodesep=0.1];\n";
  std::string URLString;

  for (const auto& it : m_graphVizFilesToNodes) {
    if (m_nodeUrl[it.first]!="") {
      URLString=" URL=\""+m_nodeUrl[it.first]+"\"";
    } else {
      URLString="";
    }
    std::string addBox = "";
    if (!m_nodeLocal[it.first]) addBox = " shape=box";
    if (m_nodeRenameMap[it.first]=="") {
      result += c_indent + ""+getNodeName(it.second)+"[label=\""+it.first+"\""+URLString+" "+addBox+"];\n";
    } else {
      result += c_indent + ""+getNodeName(it.second)+"[label=\""+m_nodeRenameMap[it.first]+"\""+addBox+URLString+"];\n";
    }
  }

  result += "\n";
  for (const auto& it : m_nodeLinks) {
    const std::vector<int>& linkedNodes = it.second;
    for (const auto& childIt : linkedNodes) {
      m_graphVizConnections += c_indent + getNodeName(it.first)+":e -> "+getNodeName(childIt)+":w;\n";
    }
  }

  result += m_graphVizConnections;
  result += "}\n";

  return result;
}

void GraphVizCreator::prepareNodeOutput(std::string absoluteFileName, std::string relativeFileName, bool webOutput) {
  bool isLocal = m_nodeLocal[absoluteFileName];

  // Name the node correctly on the output graph
  if (isLocal==false) addNodeRename(absoluteFileName, "STD/"+relativeFileName); // @include-std 
  else addNodeRename(absoluteFileName, relativeFileName); // @include

  // Prepare the correct link
  if (webOutput) addNodeUrl(absoluteFileName, destinationEncode(absoluteFileName));
  else addNodeUrl(absoluteFileName, "file://"+absoluteFileName);

}

std::string GraphVizCreator::destinationEncode(const std::string& originalFileName) {
  boost::filesystem::path p(originalFileName);
  return std::to_string(getFileId(originalFileName))+"-"+p.filename().string();
}

void GraphVizCreator::createEncodedFileList(std::vector<std::string>& originalFileNames, std::vector<std::string>& destinationFileNames) {

  originalFileNames.clear();
  destinationFileNames.clear();

  for (auto& aFileToNode : m_graphVizFilesToNodes ) {
    const std::string& origFileName = aFileToNode.first;
    std::string destFileName = destinationEncode(aFileToNode.first);
    originalFileNames.push_back(origFileName);
    destinationFileNames.push_back(destFileName);
    m_nodeDestinationFilename[origFileName]=destFileName;
  }
}

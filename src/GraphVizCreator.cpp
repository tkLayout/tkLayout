#include <GraphVizCreator.hh>

int GraphVizCreator::getFileId(std::string istreamid) {
  if (graphVizFilesToNodes_[istreamid]==0) graphVizFilesToNodes_[istreamid]=graphVizFilesToNodes_.size();
  return graphVizFilesToNodes_[istreamid];
}

std::string GraphVizCreator::getNodeName(int nodeId) {
  return "node_"+std::to_string(nodeId);
}

void GraphVizCreator::addGraphLink(int parentNodeId, int childNodeId) {
  graphVizConnections_ += indent_ + getNodeName(parentNodeId)+":e -> "+getNodeName(childNodeId)+":w;\n";
}

void GraphVizCreator::addNodeRename(std::string nodeName, std::string symbolic) {
  nodeRenameMap_[nodeName]=symbolic;
}

void GraphVizCreator::addNodeUrl(std::string nodeName, std::string url) {
  nodeUrl_[nodeName]=url;
}

std::string GraphVizCreator::createGraphVizFile() {
  std::string result;
  result += "digraph aGraph {\n";
  result += indent_ + "overlap=false;\n";
  result += indent_ + "rankdir=LR;\n";
  result += indent_ + "graph [ranksep=2, nodesep=0.1];\n";
  std::string URLString;
  for (const auto& it : graphVizFilesToNodes_) {
    if (nodeUrl_[it.first]!="") {
      URLString=" URL=\""+nodeUrl_[it.first]+"\"";
    } else {
      URLString="";
    }
    if (nodeRenameMap_[it.first]=="") {
      result += indent_ + ""+getNodeName(it.second)+"[label=\""+it.first+"\""+URLString+"];\n";
    } else {
      result += indent_ + ""+getNodeName(it.second)+"[label=\""+nodeRenameMap_[it.first]+"\" shape=box"+URLString+"];\n";
    }
  }
  result += "\n";
  result += graphVizConnections_;
  result += "}\n";

  return result;
}


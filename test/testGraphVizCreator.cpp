#include <GraphVizCreator.hh>
#include <iostream>

const int node_pappa = 0;
const int node_pappa_bis = 1;
const int node_cane = 2;
const int node_cacca = 3;
const int node_altro = 4;

int main(int argc, char* argv[]) {
  GraphVizCreator GVC;
  int fileId[10];

  fileId[node_pappa] = GVC.getFileId("pappa");
  fileId[node_pappa_bis] = GVC.getFileId("pappa");
  fileId[node_cane] = GVC.getFileId("cane");
  fileId[node_cacca] = GVC.getFileId("cacca");
  fileId[node_altro] = GVC.getFileId("altro");

  GVC.addGraphLink(fileId[node_pappa], fileId[node_cane]);
  GVC.addGraphLink(fileId[node_pappa_bis], fileId[node_cane]);
  GVC.addGraphLink(fileId[node_cane], fileId[node_cacca]);
  GVC.addGraphLink(fileId[node_cane], fileId[node_altro]);
  GVC.addGraphLink(fileId[node_cane], fileId[node_altro]);

  std::cout << GVC.createGraphVizFile() ;

  return 0;

}

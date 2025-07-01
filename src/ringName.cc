#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>

void syntax(char* name) { printf("Syntax: %s [TEDD|TEPX|TFPX|TBPS|TB2S|TBPX]\n", name); }

void generateTepxRing(int ring, int nMods) {
  printf("  Ring %d {\n", ring);
  int nModsHalf = nMods/2;
  for (int i=0; i<nMods; ++i) printf("    Module %d { label \"%d%02d\" }\n", i, ring, (i%nModsHalf)+1);
  printf("  }\n");
}

void generateTepx() {
  int ringSize[] = {20, 28, 36, 44, 48};
  int nRings = sizeof(ringSize) / sizeof(ringSize[0]);
  for (int i=0; i<nRings; ++i) generateTepxRing(i+1, ringSize[i]);
}

void generateTeddRing(int ring, int nMods) {
  printf("  Ring %d {\n", ring);
  int nModsHalf = nMods/2;
  for (int i=0; i<nMods; ++i) printf("    Module %d { label \"%d-%d\" }\n", i, ring, (i%nModsHalf)+1);
  printf("  }\n");
}

void generateTedd() {
  int ringSize[] = {20, 24, 24, 28, 32, 32, 36, 40, 40, 44, 52, 60, 64, 72, 76};
  int nRings = sizeof(ringSize) / sizeof(ringSize[0]);
  for (int i=0; i<nRings; ++i) generateTeddRing(i+1, ringSize[i]);
}

int main(int argc, char* argv[]) {
  if ((argc!=2)) { syntax(argv[0]); return(-1); };

  std::string subDet = argv[1];
  if (subDet == "TEDD") generateTedd();
  else if (subDet == "TEPX") generateTepx();
  else {
    std::cerr << "Unhandled subdetector " << subDet << std::endl;
    return -1;
  }

  return 0;
}

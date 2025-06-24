#include <cstdio>
#include <cstdlib>

void syntax(char* name) { printf("Syntax: %s nModulesInRing\n", name); }
int main(int argc, char* argv[]) {
  if ((argc!=2)) { syntax(argv[0]); return(-1); };

  int nMods = atoi(argv[1]);
  double nMods_d = nMods;
  for (int i =0; i<nMods; ++i) {
    printf ("      Module %d {\n        manualPhiCenterDeg %.10f\n        yawAngleFromConfig 0\n      }\n", i, ((i/nMods_d+.5) - int(i/nMods_d+.5) - 0.5)*360 );
  }

  return 0;
}

#include <cmath>
#include <cstdio>
#include <cstdlib>
void syntax(char* name) { printf("Syntax: %s z1 z2 n delta [outerRadius]\n", name); }
int main(int argc, char* argv[]) {
  if ((argc!=5)&&(argc!=6)) { syntax(argv[0]); return(-1); };
  double z1, z2; int n, deltaDisk;
  z1=atof(argv[1]);
  z2=atof(argv[2]);
  n=atoi(argv[3]);
  deltaDisk=atoi(argv[4]);
  double r = pow(z2/z1, 1/(double(n)-1));
  printf("Max/min hit r on disk = %f\n", r);
  if (argc==6) printf("Innermost track r = %f\n\n\n", atof(argv[5])/r);
  for (int i=0; i<n; ++i) printf("    Disk %d { placeZ %.2f }\n", i+1, z1*pow(r,i));
  for (int i=0; i<n; ++i) printf("    Disk %d { destination FPIX%d }\n", i+1, i+1+deltaDisk);
  for (int i=0; i<n; ++i) printf("Station {\n  stationName FPIX%d\n  type second\n"
                                 "  minZ %.2f\n  maxZ %.2f\n"
                                 "  @include-std CMS_Phase2/Pixel/Conversions/TWP_to_GBT\n}\n", i+1+deltaDisk, z1*pow(r,i)+20, z1*pow(r,i)+70);
  return 0;
}

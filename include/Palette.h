#include <TColor.h>
#include <TROOT.h>

class Palette {
 public:
  static void prepare(unsigned int nColors, double phase=0, double luminosity=0.75, double saturation=0.8);
  static unsigned int color(unsigned int colorIndex);
 private:
  static unsigned int myColorBase;
  static unsigned int myColors;
};

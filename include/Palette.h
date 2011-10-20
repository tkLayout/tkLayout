#ifndef palette_h
#define palette_h

#include <TColor.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <map>

class Palette {
 public:
  static Color_t color(const std::string& objectName);
  static Color_t color(const unsigned int& colorIndex);
  static const Color_t color_invalid_module = kGray + 1;
 private:
  static std::map<std::string, int> colorPickMap;
  static Color_t color_int(const unsigned int& plotIndex);
  static bool initialized;
  static void initializeMe();
  // static void skipColors(unsigned int colorIndex);
  //private:
  //static unsigned int myColorBase;
  //static unsigned int myColors;
};

#endif

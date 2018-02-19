#ifndef palette_h
#define palette_h

#include <TColor.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <map>

#include <global_funcs.hh>

class Palette {
 public:
  static Color_t color(const std::string& objectName);
  static Color_t color(const unsigned int& colorIndex, bool isTransparent = false);
  static Color_t colorDTC(const int& colorIndex, bool isTransparent = false);
  static Color_t colorChannel(const int& colorIndex, bool isTransparentActivated = false);
  static Color_t colorScrabble(const int& colorIndex, bool isTransparent = false);
  static const Color_t color_invalid_module = kGray + 1;
 private:
  static std::map<std::string, int> colorPickMap;
  static Color_t color_int(const unsigned int& plotIndex, bool isTransparent = false);
  static bool initialized;
  static void initializeMe();
  static Int_t GetColorTransparent(Int_t colorIndex, Float_t ratio);
  // static void skipColors(unsigned int colorIndex);
  //private:
  //static unsigned int myColorBase;
  //static unsigned int myColors;
};

#endif

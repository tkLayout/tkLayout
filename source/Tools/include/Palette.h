#ifndef palette_h
#define palette_h

#include <TColor.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <map>

/*
 * @class Color Palette
 * @brief Helper static class providing color from ROOT color palette based on index or color name
 */
class Palette {

 public:

  //! Based on module predefined type return standard color
  static       Color_t color(const std::string& objectName);

  //! Based on index return predefined color, if index higher than color range, modulo range operator used
  static       Color_t color(const unsigned int& colorIndex);

  static const Color_t color_invalid_module = kGray + 1;

 private:

  static std::map<std::string, int> colorPickMap;

  //! Internal definition of color map
  static Color_t color_int(const unsigned int& plotIndex);
  static bool initialized;
  static void initializeMe();

}; // Class

#endif

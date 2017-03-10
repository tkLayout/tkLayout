#ifndef INLCUDE_PALETTE_H_
#define INCLUDE_PALETTE_H_

#include <TColor.h>
#include <iostream>
#include <string>
#include <map>

/*
 * @class Color Palette
 * @brief Helper static class providing color from ROOT color palette based on index or color name. In addition defines color constants etc.
 */
class Palette {

 public:

  //! Based on module predefined type return standard color
  static       Color_t color(const std::string& objectName);

  //! Based on index return predefined color, if index higher than color range, modulo range operator used
  static       Color_t color(const short& colorIndex);

  //! Based on index return predefined color for momenta
  static       Color_t colorMomenta(const short& colorIndex);

  //! Based on index return names for predefined momenta colors
  static       std::string colorMomentaNames(const short& colorIndex);

  //! Set one of Root predefined palettes for drawing TH2D etc -> default RainBow (i.e. 55)
  static       void setRootPalette(short palette=55);

  // Colors definitions
  static const Color_t color_invalid_module = kGray + 1;

  static const Color_t color_plot_background= kWhite;
  static const Color_t color_pad_background = kGray;
  static const Color_t color_grid           = kGreen-10;
  static const Color_t color_hard_grid      = kGray;

 private:

  static std::map<std::string, short> colorPickMap;

  //! Internal definition of color map: index <-> color relation
  static Color_t color_int(const short& plotIndex);

  static bool    initialized;    //!< Has been initialized
  static void    initializeMe(); //!< Standard initialization method

}; // Class

#endif /* INCLUDE_PALETTE_H_*/

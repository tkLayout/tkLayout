#include <Palette.h>
#include <TStyle.h>

bool Palette::initialized = false;  
std::map<std::string, short> Palette::colorPickMap;

void Palette::initializeMe() {
  Palette::colorPickMap["pt2S"]   = 2;
  Palette::colorPickMap["rphi"]   = 4;
  Palette::colorPickMap["stereo"] = 6;
  Palette::colorPickMap["ptIn"]   = 9;
  Palette::colorPickMap["ptPS"]   = 1;
  Palette::colorPickMap["pixel"]  = 9;
  initialized = true;
}

//
// Based on module predefined type return standard color
//
Color_t Palette::color(const std::string& type) {

  if (!initialized) initializeMe();
  if (type=="") return color_invalid_module;
  if (colorPickMap[type]==0) {
    // New type! I'll pick a new color
    int iColor=0;
    bool found=true;
    while (found) {
      ++iColor;
      found=false;
      for (auto it=colorPickMap.begin(); it!=colorPickMap.end(); ++it) {
        if (it->first==iColor) {
          found = true;
          break;
        }
      }
    }
    colorPickMap[type]=Palette::color(iColor);
  }
  return color_int(colorPickMap[type]);
}

//
// Based on index return predefined color, if index higher than color range, modulo range operator used
//
Color_t Palette::color(const short& plotIndex) {
  if (!initialized) initializeMe();
  return color_int(plotIndex);
}

//
// Based on index return predefined color for momenta
//
Color_t Palette::colorMomenta(const short& colorIndex)
{
  if      (colorIndex==0) return kBlack;
  else if (colorIndex==1) return kBlue+3;
  else if (colorIndex==2) return kBlue;
  else if (colorIndex==3) return kRed;
  else if (colorIndex==4) return kGreen;
  else if (colorIndex==5) return kMagenta;
  else if (colorIndex==6) return kGreen+4;
  else                    return kBlack;
}

//
// Based on index return names for predefined momenta colors
std::string Palette::colorMomentaNames(const short& colorIndex) {

  if      (colorIndex==0) return std::string("Black");
  else if (colorIndex==1) return std::string("DarkBlue");
  else if (colorIndex==2) return std::string("Blue");
  else if (colorIndex==3) return std::string("Red");
  else if (colorIndex==4) return std::string("Green");
  else if (colorIndex==5) return std::string("Magenta");
  else if (colorIndex==6) return std::string("DarkGreen");
  else                    return std::string("Black");
}

//
// Set one of predefined Root palette -> default kRainBow (i.e. 55)
//
void Palette::setRootPalette(short palette) {

  gStyle->SetPalette(palette);
}

//
// Internal definition of color map
//
Color_t Palette::color_int(const short& plotIndex) {
  std::string colorCode;
  
  if (plotIndex==0) colorCode = "#000000";
  else {
    int nColor=(plotIndex-1) % 12;
    switch (nColor) {
    case 0 :
      colorCode="#004586"; // Blue
      break;
    case 1 :
      colorCode="#FF420E"; // Dark Orange
      break;
    case 2 :
      colorCode="#FFD320"; // Yellow
      break;
    case 3 :
      colorCode="#579D1C"; // Green
      break;
    case 4 :
      colorCode="#7E0021"; // brown-red
      break;
    case 5 :
      colorCode="#83CAFF"; // Light blue
      break;
    case 6 :
      colorCode="#314004"; // Dark green
      break;
    case 7 :
      colorCode="#AECF00"; // Light green
      break;
    case 8 :
      colorCode="#4B1F6F"; // Violet
      break;
    case 9 :
      colorCode="#FF950E"; // Ligth orange
      break;
    case 10 :
      colorCode="#C5000B"; // Light brown-red
      break;
    case 11 :
      colorCode="#0084D1"; // Azure blue
      break;
    default :
      std::cerr << "ERROR: in Vizard::getNiceColor() n%12 is not an int between 0 and 11! This should not happen." << std::endl;
      colorCode="#000000";
      break;
    }
  }
  
  return TColor::GetColor(colorCode.c_str());
}


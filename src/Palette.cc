#include <Palette.hh>

bool Palette::initialized = false;  
std::map<std::string, int> Palette::colorPickMap;


void Palette::initializeMe() {
  Palette::colorPickMap["pt2S"] = 2;
  Palette::colorPickMap["rphi"] = 4;
  Palette::colorPickMap["stereo"] = 6;
  Palette::colorPickMap["ptIn"] = 9;
  Palette::colorPickMap["ptPS"] = 1;
  Palette::colorPickMap["pixel"] = 9;
  initialized = true;
}

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
      for (std::map<std::string, int>::iterator it=colorPickMap.begin(); it!=colorPickMap.end(); ++it) {
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

Color_t Palette::color(const unsigned int& plotIndex) {
  if (!initialized) initializeMe();
  return color_int(plotIndex);
}

Color_t Palette::color_int(const unsigned int& plotIndex) {
  std::string colorCode;
  
  if (plotIndex==0) colorCode = "#000000";
  else {
    int nColor=(plotIndex-1) % 12;
    switch (nColor) {
    case 0 :
      colorCode="#004586";
      break;
    case 1 :
      colorCode="#FF420E";
      break;
    case 2 :
      colorCode="#FFD320";
      break;
    case 3 :
      colorCode="#579D1C";
      break;
    case 4 :
      colorCode="#7E0021";
      break;
    case 5 :
      colorCode="#83CAFF";
      break;
    case 6 :
      colorCode="#314004";
      break;
    case 7 :
      colorCode="#AECF00";
      break;
    case 8 :
      colorCode="#4B1F6F";
      break;
    case 9 :
      colorCode="#FF950E";
      break;
    case 10 :
      colorCode="#C5000B";
      break;
    case 11 :
      colorCode="#0084D1";
      break;
    default :
      std::cerr << "ERROR: in Vizard::getNiceColor() n%12 is not an int between 0 and 11! This should not happen." << std::endl;
      colorCode="#000000";
      break;
    }
  }
  
  return TColor::GetColor(colorCode.c_str());
}


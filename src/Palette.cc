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

Color_t Palette::color(const unsigned int& plotIndex, bool isTransparent) {
  if (!initialized) initializeMe();
  return color_int(plotIndex, isTransparent);
}

Color_t Palette::color_int(const unsigned int& plotIndex, bool isTransparent) {
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
  
  short paletteIndex = TColor::GetColor(colorCode.c_str());
  if (isTransparent) paletteIndex = Palette::GetColorTransparent(paletteIndex, 0.2);

  return paletteIndex;
}


/*
  This allows to have one different color for each of the 12 DTC slots.
  A shift in the color scheme is also done for each of the 9 possible phi sectors.
 */
Color_t Palette::colorDTC(const int& colorIndex, bool isTransparent) {
  //TColor::CreateColorWheel();
  //return gROOT->GetColor(paletteIndex);

  const int zone = femod(colorIndex % 12, 12);  // unit digit (in a numbering of base 12)
  const int phiSector = (colorIndex - 1) / 12;  // dizain digit (in a numbering of base 12)
  
  short paletteIndex;
  if (colorIndex == 0) paletteIndex = 1;

  else {
    switch (zone) {
    case 0 :
      paletteIndex= kYellow ;
      break;
    case 1 :
      paletteIndex= kOrange;
      break;
    case 2 :
      paletteIndex= kRed;
      break;
    case 3 :
      paletteIndex=kPink;
      break;
    case 4 :
      paletteIndex=kMagenta;
      break;
    case 5 :
      paletteIndex=kViolet;
      break;
    case 6 :
      paletteIndex=kBlue;
      break;
    case 7 :
      paletteIndex=kAzure;
      break;
    case 8 :
      paletteIndex=kCyan;
      break;
    case 9 :
      paletteIndex=kTeal;
      break;
    case 10 :
      paletteIndex=kGreen;
      break;
    case 11 :
      paletteIndex=kSpring;
      break;
    default :
      std::cerr << "ERROR: modulo 12" << std::endl;
      paletteIndex=kWhite;
      break;
    }

    paletteIndex -= (colorIndex % 10);  // should be -= phiSector, but decision was made to keep things like this, since color scheme cannot be perfectly unique anyway.
    if (isTransparent) paletteIndex = Palette::GetColorTransparent(paletteIndex, 0.2);
  }
 
  return paletteIndex;
}


/* This allows to have one color for each of the 12 services channels.
   If 12 is added, the color is set to transparent (if transparent colors allowed by isTransparentActivated).
 */
Color_t Palette::colorChannel(const int& colorIndex, bool isTransparentActivated) {

  const int zone = femod(colorIndex % 12, 12);  // unit digit (in a numbering of base 12)
  const int shift = (colorIndex - 1) / 12;      // dizain digit (in a numbering of base 12)
  
  short paletteIndex;

  if (colorIndex == 0) paletteIndex = 1;

  else {
    switch (zone) {
    case 1 :
      paletteIndex= kYellow;
      break;
    case 2 :
      paletteIndex= kOrange - 3;
      break;
    case 3 :
      paletteIndex= kOrange + 3;
      break;
    case 4 :
      paletteIndex=kRed;
      break;
    case 5 :
      paletteIndex=kGray + 1;
      break;
    case 6 :
      paletteIndex=kMagenta;
      break;
    case 7 :
      paletteIndex=kViolet - 6;
      break;
    case 8 :
      paletteIndex=kBlue + 1;
      break;
    case 9 :
      paletteIndex=kAzure + 1;
      break;
    case 10 :
      paletteIndex=kCyan;
      break;
    case 11 :
      paletteIndex=kGreen + 2;
      break;
    case 0 :
      paletteIndex=kSpring;
      break;
    default :
      std::cerr << "ERROR: modulo 12" << std::endl;
      paletteIndex=kWhite;
      break;
    }

    if (isTransparentActivated) {
      const bool isTransparent = (shift >= 1); // set transparent if 12 has been added
      if (isTransparent) paletteIndex = Palette::GetColorTransparent(paletteIndex, 0.1);
    }
  }
 
  return paletteIndex;
}


// TO DO : Why the hell is TColor::GetColorTransparent not recognized as a method of TColor ?? 
// Temporary : use this instead.  
Int_t Palette::GetColorTransparent(Int_t colorIndex, Float_t ratio) {
  if (colorIndex < 0) return -1;

  TColor* color = gROOT->GetColor(colorIndex);
  if (color) {
    TColor* transColor = new TColor(gROOT->GetListOfColors()->GetLast()+1,
				    color->GetRed(), color->GetGreen(), color->GetBlue());
    transColor->SetAlpha(ratio);
    transColor->SetName(Form("%s_transparent", color->GetName()));
    return transColor->GetNumber();
  } else {
    std::cerr << "TColor::GetColorTransparent : color with undefined index" << std::endl;
    return -1;
  }
}

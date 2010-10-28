#include <Palette.h>

unsigned int Palette::myColorBase=10000;
unsigned int Palette::myColors=0;

// Skips a number of colors on the palette
// @param nColors the number of colors to be skipped
void Palette::skipColors(unsigned int nColors) {
  myColorBase+=nColors;
}

// Sets the colors of my palette
// between myColorBase and myColorBase+nColors
// @param nColors number of colors to be booked
// @param phase the absolute phase of the color multiplet
//        in the hue circle, starting from 210
// @param luminosity the colors' luminosity (defaults to 0.75)
// @param saturation the colors' saturation (defaults to 0.8)
void Palette::prepare(unsigned int nColors,
		      double phase /*=0*/,
		      double luminosity /*=0.75*/,
		      double saturation /*=0.8*/) {
  myColorBase+=myColors;
  myColors = nColors;
  
  Float_t r, g, b;
  Float_t h, l, s;
  //TColor* aColor;
  
  l=luminosity;  // luminosity
  s=saturation;  // saturation

  // nColors different colors
  for (unsigned int iColor=0; iColor<nColors; ++iColor) {
    h = double(iColor)/double(nColors)*360.; 
    //h += phase+210;
    h -= int(h/360)*360.;
    TColor::HLStoRGB(h, l, s, r, g, b);
    TColor* myColor = gROOT->GetColor(iColor+myColorBase);
    if (myColor) myColor->SetRGB(r, g, b);
    else new TColor(iColor+myColorBase, r, g, b);
  }
}

unsigned int Palette::color(unsigned int colorIndex) {
  if (myColors==0) return 0;
  if (colorIndex<0) return 0;
  colorIndex = (colorIndex % myColors);
  if (colorIndex>=myColors) return 0;
  return myColorBase+colorIndex;
}

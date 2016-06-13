#include "PlotDrawer.h"
#include "TLatex.h"
#include <global_constants.h>
#include <iomanip>

int FrameGetterBase::s_id = 0;
const int FrameGetterBase::s_nBinsZoom = 1000;

//
// Create frame, i.e. histogram containing all modules in RZ projection (R>=0, Z -0+)
//
template<> TH2C* FrameGetter<RZFull>::operator()(double viewportZ, double viewportR) const {
  std::string name = std::string("frameRZFull") + nextString();
  TH2C* frame = new TH2C(name.c_str(), ";z [mm];r [mm]", s_nBinsZoom, -viewportZ, viewportZ, s_nBinsZoom, 0, viewportR);
  frame->GetXaxis()->SetTitleOffset(1.3);
  frame->SetStats(kFALSE);
  return frame;
}

//
// Create frame, i.e. histogram containing all modules in RZ projection (R>=0, Z>=0)
//
template<> TH2C* FrameGetter<RZ>::operator()(double viewportZ, double viewportR) const {
  std::string name = std::string("frameRZ") + nextString();
  TH2C* frame = new TH2C(name.c_str(), ";z [mm];r [mm]", s_nBinsZoom, 0, viewportZ, s_nBinsZoom, 0, viewportR);
  frame->GetXaxis()->SetTitleOffset(1.3);
  frame->SetStats(kFALSE);
  return frame;
}

//
// Create frame, i.e. histogram containing all modules in XY projection
//
template<> TH2C* FrameGetter<XY>::operator()(double viewportX, double viewportY) const {
  std::string name = std::string("frameYZ") + nextString();
  TH2C* frame = new TH2C(name.c_str(), ";x [mm];y [mm]", s_nBinsZoom, -viewportX, viewportX, s_nBinsZoom, -viewportY, viewportY);
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->SetStats(kFALSE);
  return frame;
}

template<class CoordType> void TicksFrameStyle<CoordType>::drawEtaTicks(double maxL, double maxR, double tickDistance, double tickLength, double textDistance,
                                                                       Style_t labelFont, Float_t labelSize, double etaStep, double etaMax, double etaLongLine) const {
  // Add the eta ticks
  double theta;
  double startR = maxR + tickDistance;
  double startL = maxL + tickDistance;
  double endR = maxR + tickLength + tickDistance;
  double endL = maxL + tickLength + tickDistance;
  TLine* aTick = new TLine();
  double textX = 0, textY = 0;

  TText* aLabel;
  char labelChar[10];
  double eta;
  
  // Use epsilon to avoid rounding errors
  etaMax = etaMax+vis_step_eta_epsilon;

  double thetaLimit = atan(startR/startL);
  std::vector<double> etaSteps;
  for (eta=0; eta<=etaMax; eta+=etaStep) {
 
    etaSteps.push_back(eta);
}    

  if (etaLongLine>0) etaSteps.push_back(etaLongLine);

  for (std::vector<double>::iterator it = etaSteps.begin(); it!=etaSteps.end(); ++it) {
    eta=*it;
    theta = 2 * atan(exp(-eta));
    
    aTick->SetLineColor(kBlack);
    aTick->SetLineWidth(1);
    if (theta>thetaLimit) {
      aTick->DrawLine(startR / tan(theta), startR, endR / tan(theta),  endR);
    } else {
      aTick->DrawLine(startL, startL * tan(theta), endL, endL * tan(theta));
    }
    if ((eta==etaLongLine)&&(etaLongLine!=0)) {
      aTick->SetLineColor(kGray);
      if (theta>thetaLimit) aTick->DrawLine(0., 0., endR / tan(theta), endR);
      else aTick->DrawLine(0., 0., endL, endL * tan(theta));
      labelSize*=2;
    }

    if (labelSize!=0) {
      textX = (endR+textDistance) / tan(theta);
      textY = (endL + textDistance) * tan(theta);
      if (textX>endL+textDistance) textX = endL+textDistance;
      if (textY>endR+textDistance) textY = endR+textDistance;

      sprintf(labelChar, "%.01f", eta);
      aLabel = new TText(textX, textY, labelChar);
      aLabel->SetTextAlign(21);
      aLabel->SetTextSize(labelSize);
      aLabel->SetTextFont(labelFont);
      aLabel->Draw();
    }
  }

  if (labelSize>0) {
    textY-=3*tickLength;
    aLabel = new TLatex(textX, textY, "#eta");
    aLabel->SetTextFont(labelFont);
    aLabel->SetTextSize(labelSize);
    aLabel->Draw();
  }
}

template<class CoordType> void TicksFrameStyle<CoordType>::drawEtaTicks(double maxL, double maxR, double tickDistance, double scaleFactorX, double scaleFactorY,
                                                                       Style_t labelFont, Float_t stdLabelSize, double etaStepShort, double etaStepLong, double etaMax, double etaLongLineI, double etaLongLineII,
                                                                       bool zPlsMin) const {
  // Add the eta ticks
  double  theta;
  double  thetaSign;
  double  tickLengthY   = maxR/scaleFactorY;
  double  textDistanceY = maxR/scaleFactorY;
  double  tickLengthX   = maxL/scaleFactorX;
  double  textDistanceX = maxL/scaleFactorX;
  
  double  startR       = maxR + tickDistance;
  double  startL       = maxL + tickDistance;
  double  endR         = maxR + tickLengthY + tickDistance;
  double  endL         = maxL + tickLengthX + tickDistance;
  TLine*  aTick        = new TLine();
  
  double  textX = 0, textY = 0;

  TText*  aLabel;
  int     labelColor = kBlack;
  Float_t labelSize  = stdLabelSize;
  int     labelAlign = 21;
  char    labelChar[10];
  double  eta;
  
  // Use epsilon to avoid rounding errors when drawing eta ticks
  double etaMaxSafe = etaMax+vis_step_eta_epsilon;

  double thetaLimit = atan(startR/startL);
  std::vector<double> etaSteps;
  std::vector<double> etaStepsFull;

  // Start with short eta step, continue with long eta step when drawing eta ticks 
  for (eta=0; eta<=etaMaxSafe+etaStepShort; eta+=etaStepShort) {
 
    etaSteps.push_back(eta);
  }
  for (eta=etaMax+etaStepLong; eta<etaLongLineI; eta+=etaStepLong) {
    
    etaSteps.push_back(eta);
  }    
  for (eta=etaLongLineI; eta<etaLongLineII; eta+=2*etaStepLong) {
    
    etaSteps.push_back(eta);
  }   

  // Add final eta tick
  etaSteps.push_back(etaLongLineII);

  // Create ticks for both sides
  if (zPlsMin) etaStepsFull.resize(2*etaSteps.size());
  else         etaStepsFull.resize(etaSteps.size());

  int iStep=0;
  if (zPlsMin)  for (auto it = etaSteps.end(); it!=etaSteps.begin(); it--, iStep++) etaStepsFull[iStep] = -1 * *it;
  for (auto it = etaSteps.begin(); it!=etaSteps.end(); it++, iStep++) etaStepsFull[iStep] = *it;

  for (auto it = etaStepsFull.begin(); it!=etaStepsFull.end(); ++it) {
    eta       =*it;
    thetaSign = eta>=0 ? +1 : -1;
    theta     = 2 * atan(exp(-fabs(eta)));
    
    aTick->SetLineColor(kBlack);
    aTick->SetLineWidth(1);
    if (theta>thetaLimit) {
      if (zPlsMin && eta==0) aTick->DrawLine(0, 0, endR / tan(theta)*thetaSign,  endR);
      else                   aTick->DrawLine(startR / tan(theta)*thetaSign, startR, endR / tan(theta)*thetaSign,  endR);
    } else {
      aTick->DrawLine(startL, startL * tan(theta), endL, endL * tan(theta));
    }
    if ((eta==etaLongLineI)&&(etaLongLineI!=0)) {
      aTick->SetLineColor(kBlue); // blue
      if (theta>thetaLimit) aTick->DrawLine(0., 0., endR / tan(theta)*thetaSign, endR);
      else aTick->DrawLine(0., 0., endL*thetaSign, endL * tan(theta));
      labelSize  = 1.2*stdLabelSize;
      labelAlign = 21;
      labelColor = kBlue;
    }
    else if ((eta==etaLongLineII)&&(etaLongLineII!=0)) {
      aTick->SetLineColor(kGreen); // green
      if (theta>thetaLimit) aTick->DrawLine(0., 0., endR / tan(theta)*thetaSign, endR);
      else aTick->DrawLine(0., 0., endL*thetaSign, endL * tan(theta));
      labelSize  = 1.5*stdLabelSize;
      labelAlign = 23; 
      labelColor = kGreen;
    }
    else {
      labelSize  = stdLabelSize;
      labelAlign = 21;
      labelColor = kBlack; 
     }

    if (theta<thetaLimit && thetaSign<0) labelSize = 0; // Avoid collision of ticks with axis description

    if (labelSize!=0) {
      textX = (endR +textDistanceY) / tan(theta);
      textY = (endL +textDistanceX) * tan(theta);
      if (textX>endL+textDistanceX) textX = endL+textDistanceX;
      if (textY>endR+textDistanceY) textY = endR+textDistanceY;
      textX *= thetaSign;

      sprintf(labelChar, "%.01f", eta);

      aLabel = new TText(textX, textY, labelChar);
      aLabel->SetTextAlign(labelAlign);
      aLabel->SetTextSize(labelSize);
      aLabel->SetTextFont(labelFont);
      aLabel->SetTextColor(labelColor);
      aLabel->Draw();
    }
  }

  if (labelSize>0) {
    //textY-=3*tickLengthY;
    textY = -2*tickLengthY;
    aLabel = new TLatex(textX, textY, "#eta");
    aLabel->SetTextFont(labelFont);
    aLabel->SetTextSize(1.5*stdLabelSize);
    aLabel->Draw();
  }
}

//
// Draw using style when RZ projection required
//
template<> void TicksFrameStyle<RZ>::operator()(TH2C& frame, TCanvas& canvas, DrawerPalette& palette, bool isPixelType) const {
  frame.Draw();

  //drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 50, 50, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), step_eta_normal, trk_eta_coverage, max_eta_coverage);
  if (isPixelType) drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 60, 24, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), vis_step_eta_long, vis_step_eta_long, vis_short_eta_coverage, vis_trk_eta_coverage, geom_max_eta_coverage);
  else             drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 60, 24, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), vis_step_eta_short, vis_step_eta_long, vis_short_eta_coverage, vis_trk_eta_coverage, geom_max_eta_coverage);
}

//
// Draw using style when RZFull projection required
//
template<> void TicksFrameStyle<RZFull>::operator()(TH2C& frame, TCanvas& canvas, DrawerPalette& palette, bool isPixelType) const {
  frame.Draw();

  //drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 50, 50, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), step_eta_normal, trk_eta_coverage, max_eta_coverage);
  if (isPixelType) drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 60, 24, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), vis_step_eta_long, vis_step_eta_long, vis_short_eta_coverage, vis_trk_eta_coverage, geom_max_eta_coverage, true);
  else             drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 60, 24, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), vis_step_eta_short, vis_step_eta_long, vis_short_eta_coverage, vis_trk_eta_coverage, geom_max_eta_coverage, true);
}

//
// Draw using style when XY projection required
//
template<> void TicksFrameStyle<XY>::operator()(TH2C& frame, TCanvas& canvas, DrawerPalette& palette, bool isPixelType) const {
  frame.Draw();    
}

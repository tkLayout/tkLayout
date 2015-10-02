#include <PlotDrawer.h>
#include <global_constants.h>
#include <iomanip>

int IdMaker::id = 0;
const int IdMaker::nBinsZoom = 1000;

template<> TH2C* FrameGetter<YZFull>::operator()(double viewportX, double viewportY) const {
  std::string name = std::string("frameYZ") + nextString();
  TH2C* frame = new TH2C(name.c_str(), ";z [mm];r [mm]", nBinsZoom, -viewportX, viewportX, nBinsZoom, 0, viewportY);
  frame->GetXaxis()->SetTitleOffset(1.3);
  return frame;
}



template<> TH2C* FrameGetter<YZ>::operator()(double viewportX, double viewportY) const {
  std::string name = std::string("frameYZ") + nextString();
  TH2C* frame = new TH2C(name.c_str(), ";z [mm];r [mm]", nBinsZoom, 0, viewportX, nBinsZoom, 0, viewportY);
  frame->GetXaxis()->SetTitleOffset(1.3);
  //    frame->GetXaxis()->SetTickLength(-0.03);
  //    frame->GetXaxis()->SetLabelOffset(0.03);
  //    frame->GetXaxis()->SetNdivisions(10);

  //    frame->GetYaxis()->SetTitleOffset(1.2);
  //    frame->GetYaxis()->SetTickLength(-0.015);
  //    frame->GetYaxis()->SetLabelOffset(0.015);
  //    frame->GetYaxis()->SetNdivisions(10);

  return frame;
}

int g
;

template<> TH2C* FrameGetter<XY>::operator()(double viewportX, double viewportY) const {
  std::string name = std::string("frameYZ") + nextString();
  TH2C* frame = new TH2C(name.c_str(), ";x [mm];y [mm]", nBinsZoom, -viewportX, viewportX, nBinsZoom, -viewportY, viewportY);
  frame->GetYaxis()->SetTitleOffset(1.3);
  return frame;
}

TPolyLine* drawMod() {
  double x[] = { 131., 31., 31., 131., 131., 101. };
  double y[] = { 132., 132., 32., 32., 81., 81. };
  return new TPolyLine(6, x, y);
}

template<class CoordType> void SummaryFrameStyle<CoordType>::drawEtaTicks(double maxL, double maxR, double tickDistance, double tickLength, double textDistance,
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
  etaMax = etaMax+insur::vis_step_eta_epsilon;

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

template<class CoordType> void SummaryFrameStyle<CoordType>::drawEtaTicks(double maxL, double maxR, double tickDistance, double scaleFactorX, double scaleFactorY,
                                                                          Style_t labelFont, Float_t stdLabelSize, double etaStepShort, double etaStepLong, double etaMax, double etaLongLineI, double etaLongLineII) const {
  // Add the eta ticks
  double  theta;
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
  double etaMaxSafe = etaMax+insur::vis_step_eta_epsilon;

  double thetaLimit = atan(startR/startL);
  std::vector<double> etaSteps;

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
    if ((eta==etaLongLineI)&&(etaLongLineI!=0)) {
      aTick->SetLineColor(4); // blue
      if (theta>thetaLimit) aTick->DrawLine(0., 0., endR / tan(theta), endR);
      else aTick->DrawLine(0., 0., endL, endL * tan(theta));
      labelSize  = 1.2*stdLabelSize;
      labelAlign = 21;
      labelColor = 4; 
    }
    else if ((eta==etaLongLineII)&&(etaLongLineII!=0)) {
      aTick->SetLineColor(3); // green
      if (theta>thetaLimit) aTick->DrawLine(0., 0., endR / tan(theta), endR);
      else aTick->DrawLine(0., 0., endL, endL * tan(theta));
      labelSize  = 1.5*stdLabelSize;
      labelAlign = 23; 
      labelColor = 3; 
    }
    else {
      labelSize  = stdLabelSize;
      labelAlign = 21;
      labelColor = kBlack; 
     }

    if (labelSize!=0) {
      textX = (endR +textDistanceY) / tan(theta);
      textY = (endL +textDistanceX) * tan(theta);
      if (textX>endL+textDistanceX) textX = endL+textDistanceX;
      if (textY>endR+textDistanceY) textY = endR+textDistanceY;

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


template<> void SummaryFrameStyle<YZ>::operator()(TH2C& frame, TCanvas& canvas, DrawerPalette& palette, bool isPixelType) const {
  frame.Draw();

  //drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 50, 50, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), insur::step_eta_normal, insur::trk_eta_coverage, insur::max_eta_coverage);
  if (isPixelType) drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 60, 24, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), insur::vis_step_eta_long, insur::vis_step_eta_long, insur::vis_short_eta_coverage, insur::vis_trk_eta_coverage, insur::geom_max_eta_coverage);
  else             drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 60, 24, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), insur::vis_step_eta_short, insur::vis_step_eta_long, insur::vis_short_eta_coverage, insur::vis_trk_eta_coverage, insur::geom_max_eta_coverage);
}



template<> void SummaryFrameStyle<XY>::operator()(TH2C& frame, TCanvas& canvas, DrawerPalette& palette, bool isPixelType) const {
  frame.Draw();    
}





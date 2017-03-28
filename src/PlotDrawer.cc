#include <PlotDrawer.hh>

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

template<class CoordType> void SummaryFrameStyle<CoordType>::drawEtaTicks(double maxL, double maxR, 
									  double tickDistanceRRatio, double tickLengthRRatio, double textDistanceRRatio,
									  double tickDistanceLRatio, double tickLengthLRatio, double textDistanceLRatio, 
									  Style_t labelFont, Float_t labelSize, 
									  double etaStep, double etaMax, double etaLongLine) const {
  // Add the eta ticks
  double theta;
  double tickDistanceR = tickDistanceRRatio * maxR;
  double tickDistanceL = tickDistanceLRatio * maxL;
  double tickLengthR = tickLengthRRatio * maxR;
  double tickLengthL = tickLengthLRatio * maxL;
  double textDistanceR = textDistanceRRatio * maxR;
  double textDistanceL = textDistanceLRatio * maxL;

  double startR = maxR + tickDistanceR;
  double startL = maxL + tickDistanceL;
  double endR = maxR + tickLengthR + tickDistanceR;
  double endL = maxL + tickLengthL + tickDistanceL;
  TLine* aTick = new TLine();
  double textX = 0, textY = 0;

  TText* aLabel;
  char labelChar[10];
  double eta;

  double thetaLimit = atan(startR/startL);
  std::vector<double> etaSteps;
  for (eta=0; eta<=etaMax; eta+=etaStep) etaSteps.push_back(eta);
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
      textX = (endR + textDistanceR) / tan(theta);
      textY = (endL + textDistanceL) * tan(theta);
      if (textX > (endL + textDistanceL)) textX = endL + textDistanceL;
      if (textY > (endR + textDistanceR)) textY = endR + textDistanceR;

      sprintf(labelChar, "%.01f", eta);
      aLabel = new TText(textX, textY, labelChar);
      aLabel->SetTextAlign(21);
      aLabel->SetTextSize(labelSize);
      aLabel->SetTextFont(labelFont);
      aLabel->Draw();
    }
  }

  if (labelSize>0) {
    textY -= 3*tickLengthR;
    aLabel = new TLatex(textX, textY, "#eta");
    aLabel->SetTextFont(labelFont);
    aLabel->SetTextSize(labelSize);
    aLabel->Draw();
  }
}


template<> void SummaryFrameStyle<YZ>::operator()(TH2C& frame, TCanvas&, DrawerPalette&) const {
  frame.Draw();
  drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 0.04, 0.04, 0, 0.017, 0.017, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), 0.2, 3.4, 4);
}

template<> void SummaryFrameStyle<XY>::operator()(TH2C& frame, TCanvas&, DrawerPalette&) const {
  frame.Draw();    
}





#include <PlotDrawer.h>

int IdMaker::id = 0;


template<> TH2D* FrameGetter<YZ>::operator()(double width, double height) const {
    std::string name = std::string("frameYZ") + nextString();
    TH2D* frame = new TH2D(name.c_str(), ";z [mm];r [mm]", 1, 0, width, 1, 0, height);
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

template<> TH2D* FrameGetter<XY>::operator()(double width, double height) const {
    std::string name = std::string("frameYZ") + nextString();
    TH2D* frame = new TH2D(name.c_str(), ";x [mm];y [mm]", 1, -width/2, width/2, 1, -height/2, height/2);
    frame->GetYaxis()->SetTitleOffset(1.3);
    return frame;
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
    double textX, textY;

    TText* aLabel;
    char labelChar[10];
    double eta;

    double thetaLimit = atan(startR/startL);
    std::vector<double> etaSteps;
    for (eta=0; eta<=etaMax; eta+=etaStep) etaSteps.push_back(eta);
    if (etaLongLine>0) etaSteps.push_back(2.5);

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


template<> void SummaryFrameStyle<YZ>::operator()(TH2D& frame, TCanvas&, DrawerPalette&) const {
    frame.Draw();
    drawEtaTicks(frame.GetXaxis()->GetXmax(), frame.GetYaxis()->GetXmax(), 0, 50, 50, frame.GetXaxis()->GetLabelFont(), frame.GetXaxis()->GetLabelSize(), 0.2, 2.2, 2.5);
}

template<> void SummaryFrameStyle<XY>::operator()(TH2D& frame, TCanvas&, DrawerPalette&) const {
    frame.Draw();    
}





void setCleanStyle() {
  gStyle->SetNdivisions(510, "x");
  gStyle->SetNdivisions(510, "y");
  gStyle->SetNdivisions(510, "z");
  gStyle->SetAxisColor(1, "x");
  gStyle->SetAxisColor(1, "y");
  gStyle->SetAxisColor(1, "z");
  gStyle->SetLabelColor(1, "x");
  gStyle->SetLabelColor(1, "y");
  gStyle->SetLabelColor(1, "z");
  gStyle->SetLabelFont(62, "x");
  gStyle->SetLabelFont(62, "y");
  gStyle->SetLabelFont(62, "z");
  gStyle->SetLabelOffset(0.005, "x");
  gStyle->SetLabelOffset(0.005, "y");
  gStyle->SetLabelOffset(0.005, "z");
  gStyle->SetLabelSize(0.04, "x");
  gStyle->SetLabelSize(0.04, "y");
  gStyle->SetLabelSize(0.04, "z");
  gStyle->SetTickLength(0.03, "x");
  gStyle->SetTickLength(0.03, "y");
  gStyle->SetTickLength(0.03, "z");
  gStyle->SetTitleOffset(1, "x");
  gStyle->SetTitleOffset(1, "y");
  gStyle->SetTitleOffset(1, "z");
  gStyle->SetTitleSize(0.04, "x");
  gStyle->SetTitleSize(0.04, "y");
  gStyle->SetTitleSize(0.04, "z");
  gStyle->SetTitleColor(1, "x");
  gStyle->SetTitleColor(1, "y");
  gStyle->SetTitleColor(1, "z");
  gStyle->SetTitleFont(62, "x");
  gStyle->SetTitleFont(62, "y");
  gStyle->SetTitleFont(62, "z");
  gStyle->SetBarWidth(1);
  gStyle->SetBarOffset(0);
  gStyle->SetDrawBorder(0);
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptFit(0);
  gStyle->SetNumberContours(20);
  gStyle->GetAttDate()->SetTextFont(62);
  gStyle->GetAttDate()->SetTextSize(0.025);
  gStyle->GetAttDate()->SetTextAngle(0);
  gStyle->GetAttDate()->SetTextAlign(11);
  gStyle->GetAttDate()->SetTextColor(1);
  gStyle->SetDateX(0.01);
  gStyle->SetDateY(0.01);
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(0.5);
  gStyle->SetFuncColor(1);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(3);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetLegendBorderSize(4);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameFillStyle(1001);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameBorderMode(1);
  gStyle->SetHistFillColor(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistFillStyle(1001);
  gStyle->SetHistLineStyle(1);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistMinimumZero(kFALSE);
  gStyle->SetCanvasPreferGL(kFALSE);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderSize(2);
  gStyle->SetCanvasBorderMode(1);
  gStyle->SetCanvasDefH(500);
  gStyle->SetCanvasDefW(700);
  gStyle->SetCanvasDefX(10);
  gStyle->SetCanvasDefY(10);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderSize(2);
  gStyle->SetPadBorderMode(1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadGridX(kFALSE);
  gStyle->SetPadGridY(kFALSE);
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetPaperSize(20, 26);
  gStyle->SetScreenFactor(0.88);
  gStyle->SetStatColor(19);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatBorderSize(2);
  gStyle->SetStatFont(62);
  gStyle->SetStatFontSize(0);
  gStyle->SetStatStyle(1001);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.995);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.16);
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTitleAlign(13);
  gStyle->SetTitleFillColor(19);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleBorderSize(2);
  gStyle->SetTitleFont(62);
  gStyle->SetTitleFontSize(0);
  gStyle->SetTitleStyle(1001);
  gStyle->SetTitleX(0.01);
  gStyle->SetTitleY(0.995);
  gStyle->SetTitleW(0);
  gStyle->SetTitleH(0);
  gStyle->SetLegoInnerR(0.5);

  Int_t fPaletteColor[50] = {19, 18, 17, 16, 15, 14, 13, 12, 11, 
			     20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
			     30, 8, 31, 32, 33, 34, 35, 36, 37, 38, 
			     39, 40, 9, 41, 42, 43, 44, 45, 47, 48, 
			     49, 46, 50, 2, 7, 6, 5, 4, 3, 2, 1};
  gStyle->SetPalette(50, fPaletteColor);

  TString fLineStyleArrayTmp[30] = {"", "  ", " 12 12", " 4 8", 
				    " 12 16 4 16", " 20 12 4 12", " 20 12 4 12 4 12 4 12", " 20 20", " 20 12 4 12 4 12", 
				    " 80 20", " 80 40 4 40", "  ", "  ", "  ", 
				    "  ", "  ", "  ", "  ", "  ", 
				    "  ", "  ", "  ", "  ", "  ", 
				    "  ", "  ", "  ", "  ", "  ", "  "};
  for (Int_t i=0; i<30; i++)
    gStyle->SetLineStyleString(i, fLineStyleArrayTmp[i]);

  gStyle->SetHeaderPS("");
  gStyle->SetTitlePS("");
  gStyle->SetFitFormat("5.4g");
  gStyle->SetPaintTextFormat("g");
  gStyle->SetLineScalePS(3);
  gStyle->SetColorModelPS(0);
  gStyle->SetTimeOffset(788918400);

  gStyle->SetLineColor(1);
  gStyle->SetLineStyle(1);
  gStyle->SetLineWidth(1);
  gStyle->SetFillColor(0);
  gStyle->SetFillStyle(1001);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetMarkerStyle(1);
  gStyle->SetTextAlign(11);
  gStyle->SetTextAngle(0);
  gStyle->SetTextColor(1);
  gStyle->SetTextFont(62);
  gStyle->SetTextSize(0.05);
}

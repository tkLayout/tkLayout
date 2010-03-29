#!/bin/bash

FirstPlot=$1
SecondPlot=$2

function usage() {
  myOwnName=`basename $0`
  echo Syntax:  $myOwnName FirstPlot SecondPlot
  echo The plot name comprises the base directory
  echo Example: $myOwnName  ~/Desktop/matsum/Outer_ecsq ./myDir/matsum/2pt_onlyrest
  exit 0
}

if [ "$FirstPlot" == "" ] || [ "$SecondPlot" == "" ] ; then
  usage
else
  FirstFile="$FirstPlot.global.C"
  SecondFile="$SecondPlot.global.C"
  FirstPlot=`basename "$FirstPlot"`
  SecondPlot=`basename "$SecondPlot"`
  echo looking for $FirstFile and $SecondFile
  [ -f "$FirstFile" ] || { echo $FirstFile not found; exit; }
  [ -f "$SecondFile" ] || { echo $SecondFile not found; exit; }
fi

MyScript="
gROOT->ProcessLine(\".x $FirstFile\");
TH1D* FirstPlot = rglobal;
matbudgetcanvas->Close();
gROOT->ProcessLine(\".x $SecondFile\");
TH1D* SecondPlot = rglobal;
matbudgetcanvas->Close();

FirstPlot->SetFillColor(kGray);
FirstPlot->SetFillStyle(1001);
FirstPlot->GetListOfFunctions()->RemoveAll();
SecondPlot->SetFillColor(kRed);
SecondPlot->SetFillStyle(3001);
SecondPlot->GetListOfFunctions()->RemoveAll();

gStyle->SetOptStat(0);
gStyle->SetCanvasColor(kWhite);
c1 = new TCanvas();

FirstPlot->Draw();
SecondPlot->Draw(\"same\");
c1->SaveAs(\"$FirstPlot-vs-$SecondPlot.png\");
printf(\"Plot file '%s' should have been created.\n\", \"$FirstPlot-vs-$SecondPlot.png\");
"

echo "$MyScript" | root -l -b


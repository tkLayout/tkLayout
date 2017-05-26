#include <string>
#include <iostream>
#include <algorithm>
#include "tclap/CmdLine.h"

#include "module.hh"
#include "layer.hh"
#include "tracker.hh"

#include "TStyle.h"
#include "stylePlot.h"

using namespace TCLAP;
using namespace std;

std::map<int,int> ringDirective;
std::map<int,int>::iterator ringDirectiveIt;

std::map<int,double> layerDirective;
std::map<int,double>::iterator layerDirectiveIt;

bool doSimulate;
bool doPlaceOnly;

std::string tkName;
double barrelToEndcap;
double zError;
double bigDelta;
double smallDelta;
double overlap;
double etaCut;

int nBarrelLayers;
int nBarrelModules;
double barrelRhoIn;
double barrelRhoOut;

int nEndcapDisks;
double endcapRhoIn;
double endcapRhoOut;
double endcapMaxZ;
int endcapDiskParity;

double summaryPtCost;
double summaryPtPower;
double summaryStripCost;
double summaryStripPower;

bool fullGeo;

bool parseValues(int argc, char** argv)
{
  // Wrap everything in a try block.  Do this every time, 
  // because exceptions will be thrown for problems. 
  try {  

    // Define the command line object.
    CmdLine cmd("Tests a tracker configuration", '=', "2.5");

    // All-tracker options
    ValueArg<string> nameArg("n","trackerName","Name of the tracker",false,"aTestTracker","string");
    cmd.add( nameArg );
    ValueArg<double> smallDeltaArg("","trackerSmallDelta","Distance between modules on same layer surface",false,2.,"double");
    cmd.add(smallDeltaArg);
    ValueArg<double> bigDeltaArg("","trackerBigDelta","Distance between layer surfaces",false, 12., "double");
    cmd.add(bigDeltaArg);
    ValueArg<double> barrelToEndcapArg("","trackerBarrelToEndcap","Distance between barrel and endcap",false, 200., "double");
    cmd.add(barrelToEndcapArg);
    ValueArg<double> zErrorArg("","trackerZError","The sigma of the z coordinate of interaction point",false, 70., "double");
    cmd.add(zErrorArg);
    ValueArg<double> overlapArg("","trackerOverlap","Desired overlap of modules",false, 1., "double");
    cmd.add(overlapArg);
    ValueArg<double> etaCutArg("","trackerEtaCut","Cut modules with centre over this eta",false, 2.5, "double");
    cmd.add(etaCutArg);


    // Barrel
    ValueArg<int> nBarrelArg("","barrelLayers","Number of barrel layers",false,6,"int");
    cmd.add( nBarrelArg );
    ValueArg<int> nBarrelModulesArg("","barrelModules","Number modules building one half string",false,12,"int");
    cmd.add( nBarrelModulesArg );
    ValueArg<double> barrelRhoInArg("","barrelInnerRadius","Inner radius of barrel",false,340.,"double");
    cmd.add( barrelRhoInArg );
    ValueArg<double> barrelRhoOutArg("","barrelOuterRadius","Outer radius of barrel",false,1080.,"double");
    cmd.add( barrelRhoOutArg );
//     ValueArg<int> barrelModulus("","barrelModulus","String number in a layer will be multiple of this",false,4,"int");
//     cmd.add( barrelModulus );

    // Endcap
    ValueArg<int> nEndcapArg("","endcapDisks","Number of endcap disks per side",false,6,"int");
    cmd.add( nEndcapArg );
    ValueArg<double> endcapRhoInArg("","endcapInnerRadius","Inner radius of disks",false,323.,"double");
    cmd.add( endcapRhoInArg );
    ValueArg<double> endcapRhoOutArg("","endcapOuterRadius","Outer radius of disks",false,1095.,"double");
    cmd.add( endcapRhoOutArg );
    ValueArg<double> endcapMaxZArg("", "endcapMaxZ", "Highest Z available for the endcaps", false, 2650., "double");
    ValueArg<int> endcapDiskParityArg("", "endcapDiskParity", "Endcap disk parity (if the inner ring is next or far from the interaction point). +1 and -1 accepted", false, +1, "int");
    cmd.add( endcapDiskParityArg );
//     ValueArg<int> endcapModulus("","endcapModulus","Modules in a disk will be multiple of this",false,4,"int");
//     cmd.add( endcapModulus );

    // Directives
    MultiArg<string> directiveRing("", "directiveRing", "Asks for a ring increase/decrease modules (e.g. : 3-2 will reduce by 2*endcapModulus the modules in ring 3)", false, "string");
    cmd.add( directiveRing );
    MultiArg<string> directiveLayer("", "directiveLayer", "Asks for a barrel layer radius. Syntax: layer#/directive (directive can be F(ix), S(hrink), E(nlarge), A(uto)or a radius)", false, "string");
    cmd.add( directiveLayer );
    
    // Commands
    SwitchArg goSwitch("s","simulate","Does the actual job", false);
    cmd.add( goSwitch );

    SwitchArg placeSwitch("p","place","Stops after having placed modules", false);
    cmd.add( placeSwitch );

    SwitchArg fullSwitch("f","fullGeometry","Saves also the solid geometry", false);
    cmd.add( fullSwitch );

    // Summary Parameters
    ValueArg<double> summaryPtCostArg("", "ptCost", "Cost of pt silicon strip detector in CHF/cm^2", false, 200., "double");
    ValueArg<double> summaryStripCostArg("", "stripCost", "Cost of silicon strip detector in CHF/cm^2", false, 40., "double");
    ValueArg<double> summaryPtPowerArg("", "ptPower", "Power consumption of pt silicon strip detector in mW/Channel", false, 0.3, "double");
    ValueArg<double> summaryStripPowerArg("", "stripPower", "Power consumption of silicon strip detector in mW/Channel", false, 0.7, "double");
    cmd.add ( summaryPtCostArg );
    cmd.add ( summaryStripCostArg );
    cmd.add ( summaryPtPowerArg );
    cmd.add ( summaryStripPowerArg );

    // Parse the args.
    cmd.parse( argc, argv );

    // Tracker parameters
    tkName = nameArg.getValue();
    barrelToEndcap = barrelToEndcapArg.getValue();
    zError = zErrorArg.getValue();
    bigDelta = bigDeltaArg.getValue();
    smallDelta = smallDeltaArg.getValue();
    overlap = overlapArg.getValue();
    etaCut = etaCutArg.getValue();

    // Parameters of barrel
    nBarrelLayers  = nBarrelArg.getValue();
    nBarrelModules = nBarrelModulesArg.getValue();
    barrelRhoIn    = barrelRhoInArg.getValue();
    barrelRhoOut   = barrelRhoOutArg.getValue();

    
    // Parameters of endcap
    nEndcapDisks     = nEndcapArg.getValue();
    endcapRhoIn      = endcapRhoInArg.getValue();
    endcapRhoOut     = endcapRhoOutArg.getValue();
    endcapMaxZ       = endcapMaxZArg.getValue();
    endcapDiskParity = ( endcapDiskParityArg.getValue() > 0 ) ? +1 : -1 ;

    // Summary parameters
    summaryPtCost     = summaryPtCostArg.getValue();;
    summaryPtPower    = summaryPtPowerArg.getValue();;
    summaryStripCost  = summaryStripCostArg.getValue();;
    summaryStripPower = summaryStripPowerArg.getValue();;


    std::ostringstream aNamePart;

    aNamePart.str("");
    aNamePart << "_" << nBarrelLayers << "L" << nBarrelModules;
    aNamePart << "_" << nEndcapDisks << "D" << ( (endcapDiskParity > 0) ? "+" : "-" );

    tkName += aNamePart.str();

    // Directive for rings
    std::vector<std::string> v = directiveRing.getValue();
    int a, b;
    for ( int z = 0; static_cast<unsigned int>(z) < v.size(); z++ ) {
      if (sscanf(v[z].c_str(), "%d%d", &a, &b)==2) {
	if (b!=0) ringDirective[a]=b;
      }
    }
    for (ringDirectiveIt = ringDirective.begin();
	 ringDirectiveIt!=ringDirective.end();
	 ringDirectiveIt++) {
      aNamePart.str("");
      aNamePart << "_R" << std::noshowpos << (*ringDirectiveIt).first
		<< std::showpos ;
      aNamePart << (*ringDirectiveIt).second;
      tkName += aNamePart.str();
    }

    // Directive for Layers
    v = directiveLayer.getValue();
    char charBuf[100];
    std::string aString;
    double aVal;
    for ( int z = 0; static_cast<unsigned int>(z) < v.size(); z++ ) {
      if (sscanf(v[z].c_str(), "%d/%s", &a, charBuf)==2) {
	aString = charBuf;
	if (aString=="F") {
	  layerDirective[a]=Layer::FIXED;
	}
	if (aString=="S") {
	  layerDirective[a]=Layer::SHRINK;
	}
	if (aString=="E") {
	  layerDirective[a]=Layer::ENLARGE;
	}
	if (aString=="A") {
	  layerDirective[a]=Layer::AUTO;
	}
	aVal = atof(aString.c_str());
	if (aVal!=0) {
	  layerDirective[a]=aVal;
	}
      }
    }
    for (layerDirectiveIt = layerDirective.begin();
	 layerDirectiveIt!=layerDirective.end();
	 layerDirectiveIt++) {
      aNamePart.str("");
      aNamePart << "_L" << std::noshowpos << (*layerDirectiveIt).first << std::showpos ;
      if ((*layerDirectiveIt).second==Layer::FIXED)  aNamePart << "F";
      if ((*layerDirectiveIt).second==Layer::SHRINK) aNamePart << "S";
      if ((*layerDirectiveIt).second==Layer::ENLARGE) aNamePart << "E";
      if ((*layerDirectiveIt).second==Layer::AUTO) aNamePart << "A";
      if ((*layerDirectiveIt).second>0) aNamePart << std::showpos << (*layerDirectiveIt).second;
      tkName += aNamePart.str();
    }
    
    doSimulate = goSwitch.getValue();
    doPlaceOnly = placeSwitch.getValue();
    fullGeo = fullSwitch.getValue();

    std::cout << "Tracker name is: " << tkName << std::endl;

    
    
  } catch (ArgException &e)  // catch any exceptions
    {
      cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
      return false;
    }
  
  return true;
  
}


int main (int argc, char* argv[]) {
  setCleanStyle();
  std::string arguments="";
  std::string anArgument;
  for (int i=1; i<argc; i++) {
    anArgument=argv[i];
    if (i!=1) arguments+=" ";
    arguments+=anArgument;
  }
  parseValues(argc, argv);

  if (!doSimulate) {
    std::cout << "Fake mode: if you're happy with this add -s" << std::endl;
    return 0;
  }

  int mySection = Layer::NoSection;

  // The tracker object (with a nice name)
  Tracker* myTracker = new Tracker(tkName);
  myTracker->setArguments(arguments);
  myTracker->setZError(zError);
  myTracker->setBigDelta(bigDelta);
  myTracker->setSmallDelta(smallDelta);
  myTracker->setOverlap(overlap);
  myTracker->setRingDirectives(ringDirective);
  myTracker->setLayerDirectives(layerDirective);
  
  // Give the summary parameters
  myTracker->setCost(Module::Pt, summaryPtCost);
  myTracker->setCost(Module::Strip, summaryStripCost);
  myTracker->setPower(Module::Pt, summaryPtPower*1e-3);
  myTracker->setPower(Module::Strip, summaryStripPower*1e-3);

  // Build the barrel with square modules
  BarrelModule* sampleBarrelModule = new BarrelModule(1.);   // Square modules of kind rphi
  myTracker->buildBarrel(nBarrelLayers,
			 barrelRhoIn,
			 barrelRhoOut,
			 nBarrelModules,
			 sampleBarrelModule, mySection, true); // Actually build compressed layer

  delete sampleBarrelModule; // Dispose of the sample module


  std::cout << "Tracker stays in the z range "
	    << myTracker->getMaxBarrelZ(-1) << ".." <<  myTracker->getMaxBarrelZ(+1)
	    << std::endl;


  // The same old sample module
  EndcapModule* sampleModule = new EndcapModule();
  myTracker->buildEndcaps(nEndcapDisks,     // nDisks (per side)
			  myTracker->getMaxBarrelZ(+1)+barrelToEndcap,  // minZ
			  endcapMaxZ,
			  endcapRhoIn,
			  endcapRhoOut,
			  sampleModule,
			  endcapDiskParity, 
			  mySection );


  // Color modules, set nStrips, segments, ecc ecc
  myTracker->setModuleTypes();


  
  // A first simple analysis
  std::pair<double, double> minMaxEta;
  minMaxEta = myTracker->getEtaMinMax();
  std::cout << "Eta coverage of the tracker (prior to module purging): " << std::endl
	    << "etaMin: " << minMaxEta.first << std::endl
	    << "etaMax: " << minMaxEta.second << std::endl;
  myTracker->cutOverEta(etaCut);
  minMaxEta = myTracker->getEtaMinMax();
  std::cout << "Eta coverage of the tracker (after module purging at eta " << etaCut << " ): " << std::endl
	    << "etaMin: " << minMaxEta.first << std::endl
	    << "etaMax: " << minMaxEta.second << std::endl;


  if (doPlaceOnly) {
    myTracker->changeRingModules("Disk4", 4, "ciao", kRed);
    std::cout << "Modules placed: now I exit, as you asked" << std::endl;
    return 0;
  }

  // Prepare the geometry objects
  if (fullGeo) {
    myTracker->createGeometry();
  }
  myTracker->createGeometry(true);

  // Optical transmission
  myTracker->computeBandwidth();

  // Analysis
  myTracker->analyze(2000, Layer::YZSection);

  // Summary and save
  myTracker->writeSummary();
  myTracker->save();
  return 0;

}

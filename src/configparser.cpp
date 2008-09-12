#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>
#include "configparser.hh"
#include "module.hh"
#include "layer.hh"
#include "tracker.hh"

using namespace std;

configParser::configParser() {
  myTracker_= NULL;
}

configParser::~configParser() {
}

string configParser::getTill(istream &inStream, char delimiter, bool singleWord,bool allowNothing=false) {
  string result="";

  if (getline(inStream, result, delimiter)) {
    inStream.get();
    if (singleWord) {
      string dummyString;
      istringstream tempStream(result.c_str());
      tempStream >> result;
      if (tempStream >> dummyString) {
	result="";
	cerr << "ERROR: unexpected extra word found: \""<<dummyString<<"\"";
      }
    }
  } else {
    result="";
    if (!allowNothing) {
      cerr << "ERROR: expecting " << delimiter << " and found end-of-file" << endl;
    }
  }

  return result;
}

bool configParser::parseParameter(string &parameterName, string &parameterValue, istream& inStream) {
  string str;
  string myLine;

  myLine=getTill(inStream, ';', false, true);
  
  if (myLine!="") {
    istringstream myLineStr(myLine);
    parameterName=getTill(myLineStr, '=', true);
    if (parameterName!="") {
      if (myLineStr >> parameterValue) {
	string dummy;
	if (myLineStr >> dummy) {
	  cerr << "WARNING: ignoring extra parameter value " << dummy << endl;
	}
	return true;
      }
    }
  } 
  parameterName="";
  parameterValue="";
  return false;
}

bool configParser::parseTracker(string myName, istream& inStream) {
  string parameterName;
  string parameterValue;
  double doubleValue;


  // Tracker is a singleton. Declare just one, please
  if (myTracker_) {
    cout << "Error: tracker is not NULL when declaring tracker " << myName << endl;
    cout << "Double declaration of a Tracker object?" << endl;
    throw parsingException();
  }

  myTracker_ = new Tracker(myName);

  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      if (parameterName=="zError") {
	doubleValue=atof(parameterValue.c_str());
	myTracker_->setZError(doubleValue);
      } else if (parameterName=="smallDelta") {
	doubleValue=atof(parameterValue.c_str());
	myTracker_->setSmallDelta(doubleValue);
      } else if (parameterName=="bigDelta") {
	doubleValue=atof(parameterValue.c_str());
	myTracker_->setBigDelta(doubleValue);
      } else if (parameterName=="overlap") {
	doubleValue=atof(parameterValue.c_str());
	myTracker_->setOverlap(doubleValue);

      } else if (parameterName=="etaCut") {
	doubleValue=atof(parameterValue.c_str());
	myTracker_->setEtaCut(doubleValue);
      } else if (parameterName=="ptCost") {
	doubleValue=atof(parameterValue.c_str());
	myTracker_->setCost(Module::Pt, doubleValue);
      } else if (parameterName=="stripCost") {
	doubleValue=atof(parameterValue.c_str());
	myTracker_->setCost(Module::Strip, doubleValue);
      } else if (parameterName=="ptPower") {
	doubleValue=atof(parameterValue.c_str());
	myTracker_->setPower(Module::Pt, doubleValue*1e-3);
      } else if (parameterName=="stripPower") {
	doubleValue=atof(parameterValue.c_str());
	myTracker_->setPower(Module::Strip, doubleValue*1e-3);
      } else {
	cout << "Unknown parameter name: " << parameterName << endl;
	throw parsingException();
      }
      cout << "\t" << parameterName << " = " << doubleValue << ";" << endl; // debug
    }
  }

  // All the tracker parameters are set now
}

bool configParser::parseBarrel(string myName, istream& inStream) {
  string parameterName;
  string parameterValue;

  int nBarrelLayers = 0;
  double barrelRhoIn = 0;
  double barrelRhoOut = 0;
  int nBarrelModules = 0;
  string aDirective = "";
  BarrelModule* sampleBarrelModule = NULL;

  // Directives (this are communicated to the Tracker object)
  std::map<int,double> layerDirective;

  // Tracker shoujld be already there
  if (!myTracker_) {
    cout << "Error: tracker is NULL when declaring barrel " << myName << endl;
    cout << "Missing declaration of a Tracker object?" << endl;
    throw parsingException();
  }


  // TODO: place the object name in a meaningful place

  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      if (parameterName=="nLayers") {
	nBarrelLayers=atoi(parameterValue.c_str());
      } else if (parameterName=="nModules") {
	nBarrelModules=atoi(parameterValue.c_str());
      } else if (parameterName=="innerRadius") {
	barrelRhoIn=atof(parameterValue.c_str());
      } else if (parameterName=="outerRadius") {
	barrelRhoOut=atof(parameterValue.c_str());
      } else if (parameterName=="directive") {
	char charBuf[100];
	std::string aString;
	double aVal;
	int layerNum;
	bool gotIt=false;
	if (sscanf(parameterValue.c_str(), "%d/%s", &layerNum, charBuf)==2) {
	  aString = charBuf;
	  if (aString=="F") {
	    layerDirective[layerNum]=Layer::FIXED;
	    gotIt=true;
	  } else if (aString=="S") {
	    layerDirective[layerNum]=Layer::SHRINK;
	    gotIt=true;
	  } else if (aString=="E") {
	    layerDirective[layerNum]=Layer::ENLARGE;
	    gotIt=true;
	  } else if (aString=="A") {
	    layerDirective[layerNum]=Layer::AUTO;
	    gotIt=true;
	  }
	  if (!gotIt) {
	    aVal = atof(aString.c_str());
	    if (aVal>0) {
	      layerDirective[layerNum]=aVal;
	      gotIt=true;
	    }
	  }
	  if (!gotIt) {
	    cout << "Parsing layer directive for barrel " << myName
		 << ": unknown/nonsense directive \"" << aString << "\"" << endl;
	    throw parsingException();
	  }
	} else {
	  cout << "Parsing barrel " << myName << endl
	       << "Wrong syntax for a layer directive: \"" << parameterValue
	       << "\" should be layer/command" << endl;
	  throw parsingException();	  
	}
      } else {
	cout << "Unknown parameter \"" << parameterName << "\"" << endl;
	throw parsingException();
      }
      cout << "\t" << parameterName << " = " << parameterValue << ";" << endl; // debug
    }
  }

  
  // Actually creating the barrel if all the mandatory parameters were set
  if ( (nBarrelLayers != 0) &&
       (barrelRhoIn != 0) &&
       (barrelRhoOut != 0) &&
       (nBarrelModules != 0) ) {
    sampleBarrelModule = new BarrelModule(1.);   // Square modules of kind rphi

    // Important: if no directive was given, the following line will clear
    // possible previous directives coming from a different barrel
    myTracker_->setLayerDirectives(layerDirective);

    myTracker_->buildBarrel(nBarrelLayers,
			   barrelRhoIn,
			   barrelRhoOut,
			   nBarrelModules,
			   sampleBarrelModule, Layer::NoSection, true); // Actually build a compressed barrel
    delete sampleBarrelModule; // Dispose of the sample module
  } else {
    cout << "Missing mandatory parameter for barrel " << myName << endl;
    throw parsingException();
  }

}

bool configParser::parseEndcap(string myName, istream &inStream) {
  string parameterName;
  string parameterValue;

  int nDisks = 0;       // nDisks (per side)
  double barrelToEndcap = 0;  // gap between barrel and endcap
  double minZ = 0;            // explicit minimum z
  double maxZ = 0;            // maximum z
  double rhoIn = 0;
  double rhoOut = 0;
  int diskParity = 0; 

  Module* sampleModule = NULL;

  // Directives (this are communicated to the Tracker object)
  std::map<int,int> ringDirective;

  // Tracker should be already there
  if (!myTracker_) {
    cout << "Error: tracker is NULL when declaring endcap " << myName << endl;
    cout << "Missing declaration of a Tracker object?" << endl;
    throw parsingException();
  }

  // TODO: place the object name in a meaningful place

  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      if (parameterName=="nDisks") {
	nDisks=atoi(parameterValue.c_str());
      } else if (parameterName=="innerRadius") {
	rhoIn=atof(parameterValue.c_str());
      } else if (parameterName=="outerRadius") {
	rhoOut=atof(parameterValue.c_str());
      } else if (parameterName=="barrelGap") {
	barrelToEndcap=atof(parameterValue.c_str());
      } else if (parameterName=="minimumZ") {
	minZ=atof(parameterValue.c_str());	
      } else if (parameterName=="maximumZ") {
	maxZ=atof(parameterValue.c_str());
      } else if (parameterName=="diskParity") {
	diskParity=atoi(parameterValue.c_str()); 
	if (diskParity>0) diskParity=1; else diskParity=-1;
      } else if (parameterName=="directive") {
	int ringNum; int increment;
	if (sscanf(parameterValue.c_str(), "%d%d", &ringNum, &increment)==2) {
	  if (increment!=0) ringDirective[ringNum]=increment;
	} else {
	  cout << "Parsing endcap " << myName << endl
	       << "Wrong syntax for a ring directive: \"" << parameterValue
	       << "\" should be ring+increment or ring-decrement )" << endl;
	  throw parsingException();	  
	}
      } else if (parameterName=="removeRings") {
	// TODO: add this
	cout << "WARNING: removeRings not implemented" << endl;
      } else {
	cout << "Unknown parameter \"" << parameterName << "\"" << endl;
	throw parsingException();
      }
      cout << "\t" << parameterName << " = " << parameterValue << ";" << endl; // debug
    }
  }

  // Chose the proper endcap minimum Z according to the parsed parameters
  if ((minZ!=0)&&(barrelToEndcap!=0)) {
    cout << "Parsing endcap " << myName << endl
	 << " I got both minimumZ = " << minZ << " and barrelGap = " << barrelToEndcap << endl
	 << "This is inconsistent: please choose one of the two methods to place the endcap" << endl;
    throw parsingException();
  }
  
  if (minZ==0) {
    minZ=myTracker_->getMaxBarrelZ(+1)+barrelToEndcap;
  }

  // Actually creating the endcap if all the mandatory parameters were set
  if ( (nDisks != 0) &&
       (rhoIn != 0) &&
       (rhoOut != 0) &&
       (minZ != 0) &&
       (maxZ != 0) &&
       (diskParity != 0)) {
    // The same old sample module
    Module* sampleModule = new Module();
    // Important: if no directive was given, the following line will clear
    // possible previous directives coming from a different endcap
    myTracker_->setRingDirectives(ringDirective);

    myTracker_->buildEndcaps(nDisks,     // nDisks (per side)
			    minZ,  // minZ
			    maxZ,
			    rhoIn,
			    rhoOut,
			    sampleModule,
			    diskParity, 
			    Layer::NoSection);

    delete sampleModule; // Dispose of the sample module
  } else {
    cout << "Missing mandatory parameter for endcap " << myName << endl;
    throw parsingException();
  }
}

bool configParser::parseType(string myType) {
  string str;
  string typeConfig;
  
  if (myType=="Tracker") {
    cout << "Reading tracker main parameters and name: ";
    str=getTill(configFile_, '{', true);
    if (str!="") { 
      cout << str << endl;
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
      	istringstream typeStream(typeConfig);
      	parseTracker(str, typeStream);
      }
    }
  } else if (myType=="Barrel") {
    cout << "CREATING BARREL:\t";
    str=getTill(configFile_, '{', true);
    if (str!="") { 
      cout << str << endl;
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
	istringstream typeStream(typeConfig);
	parseBarrel(str, typeStream);
      }
    }
  } else if (myType=="Endcap") {
    cout << "CREATING ENDCAP:\t";
    str=getTill(configFile_, '{', true);
    if (str!="") { 
      cout << str << endl;
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
	istringstream typeStream(typeConfig);
	parseEndcap(str, typeStream);
      }
    }
  } else {
    cerr << "Error: unknown piece of tracker " << myType;
    return false;
  }

  return true;
}

bool configParser::parseFile(string configFileName) {
  string str;

  if (rawConfigFile_.is_open()) {
    cerr << "config file is already open" << endl;
    return false;
  }

  rawConfigFile_.open(configFileName.c_str());
  if (rawConfigFile_.is_open()) {

    // Reset the configFile_ object
    configFile_.~stringstream();                      
    new ( (void*) &configFile_) std::stringstream();

    // Skim comments delimited by // or # and put the skimmed file into configFile_
    string::size_type locComm1; // Location of commenting substring 1
    string::size_type locComm2; // Location of commenting substring 2

    while (getline(rawConfigFile_,str)) {
      locComm1=str.find("//", 0);
      locComm2=str.find("#", 0);
      if ((locComm1==string::npos)&&(locComm2==string::npos)) {
	configFile_ << str << endl;
      } else {
	configFile_ << str.substr(0, (locComm1<locComm2 ? locComm1 : locComm2)) << endl;
      }
    }
    
    try {
      while(configFile_ >> str) {
	if (!parseType(str)) break;
      }
    } catch (exception& e) {
      cerr << e.what() << endl;
    }
    
    rawConfigFile_.close();
  } else {
    cerr << "Error: could not open config file " << configFileName << endl;
    return false;
  }

  // Eta cut and other post-operations
  std::pair<double, double> minMaxEta;
  minMaxEta = myTracker_->getEtaMinMax();
  std::cout << "Eta coverage of the tracker (prior to module purging): " << std::endl
	    << "etaMin: " << minMaxEta.first << std::endl
	    << "etaMax: " << minMaxEta.second << std::endl;
  myTracker_->cutOverEta(myTracker_->getEtaCut());
  minMaxEta = myTracker_->getEtaMinMax();
  std::cout << "Eta coverage of the tracker (after module purging at eta "
	    << myTracker_->getEtaCut() << " ): " << std::endl
	    << "etaMin: " << minMaxEta.first << std::endl
	    << "etaMax: " << minMaxEta.second << std::endl;

  
  return true;
}

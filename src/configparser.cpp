#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>
#include <cctype>
#include "configparser.hh"
#include "module.hh"
#include "layer.hh"
#include "tracker.hh"
#include <stdlib.h>

using namespace std;

configParser::configParser() {
  myTracker_= NULL;
}

configParser::~configParser() {

}

string configParser::getTill(istream &inStream, char delimiter, bool singleWord, bool allowNothing=false) {
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

// Parse a single parameter of the config file, dividing what's before the '=' to what's after
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
          std::ostringstream tempString;
          tempString.str(""); tempString << "Ignoring extra parameter value " << dummy;
          logWARNING(tempString);
        }
        return true;
      }
    }
  }
  parameterName="";
  parameterValue="";
  return false;
}

// Parses the parameter of the Tracker type in the config file.
bool configParser::parseTracker(string myName, istream& inStream) {
  string parameterName;
  string parameterValue;
  int intValue;
  double doubleValue;
  string parameterNameCopy;
  bool correctlyBroken;
  string stringIndex;

  // Tracker is a singleton. Declare just one, please
  if (myTracker_) {
    cerr << "ERROR: tracker is not NULL when declaring tracker " << myName << endl;
    cerr << "Double declaration of a Tracker object?" << endl;
    throw parsingException();
  }

  myTracker_ = new Tracker(myName);

  myTracker_->setGeometryDsDistances(geometryDsDistance_, geometryDsDistanceSecond_);

  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      parameterNameCopy = parameterName;
      correctlyBroken = breakParameterName(parameterNameCopy, stringIndex);
      if (parameterName=="nMB") { // Number of minimum bias events per BX
        doubleValue=atof(parameterValue.c_str());
        myTracker_->setNMB(doubleValue);
      } else if (parameterName=="rError") {
        doubleValue=atof(parameterValue.c_str());
        myTracker_->setRError(doubleValue);
      } else if (parameterName=="zError") {
        logWARNING("zError is deprecated, please use zErrorConstruction or zErrorCollider");
        doubleValue=atof(parameterValue.c_str());
        myTracker_->setZErrorConstruction(doubleValue);
        myTracker_->setZErrorCollider(doubleValue);
      } else if (parameterName=="zErrorConstruction") {
        doubleValue=atof(parameterValue.c_str());
        myTracker_->setZErrorConstruction(doubleValue);
      } else if (parameterName=="zErrorCollider") {
        doubleValue=atof(parameterValue.c_str());
        myTracker_->setZErrorCollider(doubleValue);
      } else if (parameterName=="efficiency") {
        doubleValue=atof(parameterValue.c_str());
        myTracker_->setEfficiency(doubleValue);
      } else if (parameterName=="pixelEfficiency") {
        doubleValue=atof(parameterValue.c_str());
        myTracker_->setPixelEfficiency(doubleValue);
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
      } else if (parameterName=="useIPConstraint") {
        intValue=atoi(parameterValue.c_str());
        if (intValue==1) myTracker_->setUseIPConstraint(true);
        else if (intValue==0) myTracker_->setUseIPConstraint(false);
        else {
          std::cerr << "ERROR: useIPConstraint can be 0 or 1" << std::endl;
          throw parsingException();
        }
      } else if (parameterName=="numInvFemtobarns") {
        doubleValue = atof(parameterValue.c_str());
        myTracker_->setNumInvFemtobarns(doubleValue);
      } else if (parameterName=="operatingTemp") {
        doubleValue = atof(parameterValue.c_str());
        myTracker_->setOperatingTemp(doubleValue);
      } else if (parameterName=="referenceTemp") {
        doubleValue = atof(parameterValue.c_str());
        myTracker_->setReferenceTemp(doubleValue);
      } else if (parameterName=="chargeDepletionVoltage") {
        intValue = atoi(parameterValue.c_str());
        myTracker_->setChargeDepletionVoltage(intValue);
      } else if (parameterName=="alphaParam") {
        doubleValue = atof(parameterValue.c_str());
        myTracker_->setAlphaParam(doubleValue);
      } else if (parameterName=="bunchSpacingNs") {
        doubleValue = atof(parameterValue.c_str());
        myTracker_->setBunchSpacingNs(doubleValue);
      } else if (parameterName=="triggerProcessorsEta") {
        intValue = atoi(parameterValue.c_str());
        myTracker_->setTriggerProcessorsEta(intValue);
      } else if (parameterName=="triggerProcessorsPhi") {
        intValue = atoi(parameterValue.c_str());
        myTracker_->setTriggerProcessorsPhi(intValue);
      } else if (parameterName=="triggerEtaCut") {
        doubleValue = atof(parameterValue.c_str());
        myTracker_->setTriggerEtaCut(doubleValue);
      } else if (parameterName=="triggerPtCut") {
        doubleValue = atof(parameterValue.c_str());
        myTracker_->setTriggerPtCut(doubleValue);        
      } else if (correctlyBroken) { // Per module type parameters
        if (parameterNameCopy == "triggerErrorIncreaseX") {
          doubleValue = atof(parameterValue.c_str());
          myTracker_->setTriggerErrorX(stringIndex, doubleValue);
        } else if (parameterNameCopy == "triggerErrorIncreaseY") {
          doubleValue = atof(parameterValue.c_str());
          myTracker_->setTriggerErrorY(stringIndex, doubleValue);
        } else if (parameterNameCopy == "opticalPowerFixed") {
          // Input is in mW, while we always store in SI units internally
          doubleValue = atof(parameterValue.c_str()) * 1e-3;
          myTracker_->setPowerFixed(stringIndex, ModuleType::OpticalPower, doubleValue);
        } else if (parameterNameCopy == "chipPowerFixed") {
          // Input is in mW, while we always store in SI units internally
          doubleValue = atof(parameterValue.c_str()) * 1e-3;
          myTracker_->setPowerFixed(stringIndex, ModuleType::ChipPower, doubleValue);
        } else if (parameterNameCopy == "opticalPower") {
          // Input is in mW, while we always store in SI units internally
          doubleValue = atof(parameterValue.c_str()) * 1e-3;
          myTracker_->setPower(stringIndex, ModuleType::OpticalPower, doubleValue);
        } else if (parameterNameCopy == "chipPower") {
          // Input is in mW, while we always store in SI units internally
          doubleValue = atof(parameterValue.c_str()) * 1e-3;
          myTracker_->setPower(stringIndex, ModuleType::ChipPower, doubleValue);
        } else if (parameterNameCopy == "sparsifiedHeaderBits") {
          intValue = atoi(parameterValue.c_str());
          myTracker_->setSparsifiedHeaderBits(stringIndex, intValue);
        } else if (parameterNameCopy == "sparsifiedPayloadBits") {
          intValue = atoi(parameterValue.c_str());
          myTracker_->setSparsifiedPayloadBits(stringIndex, intValue);
        } else if (parameterNameCopy == "triggerDataHeaderBits") { 
          intValue = atoi(parameterValue.c_str());
          myTracker_->setTriggerDataHeaderBits(stringIndex, intValue);
        } else if (parameterNameCopy == "triggerDataPayloadBits") {
          intValue = atoi(parameterValue.c_str());
          myTracker_->setTriggerDataPayloadBits(stringIndex, intValue);
        } else if (parameterNameCopy == "sensorThickness") {
          doubleValue = atof(parameterValue.c_str());
          myTracker_->setSensorThickness(stringIndex, doubleValue);
        } else if (parameterNameCopy == "inefficiencyType") {
          if (parameterValue == "edgeonly" || parameterValue == "EdgeOnly" || parameterValue == "eo") myTracker_->setInefficiencyType(stringIndex, ptError::EdgeOnly);
          else if (parameterValue == "stripwise" || parameterValue == "StripWise" || parameterValue == "sw") myTracker_->setInefficiencyType(stringIndex, ptError::StripWise);
          else {
            cerr << "ERROR: Cannot parse value: " << parameterValue << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: Unknown parameter name: " << parameterNameCopy << endl;
          throw parsingException();
        }
      } else {
        cerr << "ERROR: Unknown parameter name: " << parameterName << endl;
        throw parsingException();
      }
      // cout << "\t" << parameterName << " = " << doubleValue << ";" << endl; // debug
    }
  }

  // All the tracker parameters are set now
  return true;

}


// Parses the part of config files between '{' and '}', creating a Barrel
bool configParser::parseBarrel(string myName, istream& inStream) {
  string parameterName;
  string parameterValue;
  string parameterNameCopy;
  int mainIndex, secondaryIndex;
  bool correctlyBroken;

  int readoutMode = Module::Binary;
  int nBarrelLayers = 0;
  //double minZ=0;
  bool sameRods = false;
  bool compress = true;
  bool shortBarrel = false;
  double maxZ=0;
  double barrelRhoIn = 0;
  double barrelRhoOut = 0;
  int nBarrelModules = 0;
  double dsDistanceOverride = -1;
  string aDirective = "";
  BarrelModule* sampleBarrelModule = NULL;
  double aspectRatio = 1.;
  bool aspectRatioManual = false;
  int phiSegments = 4;
  std::pair<double, double> size; size.first=0; size.second=0; // width, length

  // Fetch the generic Delta of the tracker
  double genericSmallDelta = myTracker_->getSmallDelta();
  double genericBigDelta = myTracker_->getBigDelta();

  // Directives (this are communicated to the Tracker object)
  std::map<int, double> layerDirectives;
  std::map<int, LayerOption> layerOptions;


  // Tracker shoujld be already there
  if (!myTracker_) {
    cerr << "ERROR: tracker is NULL when declaring barrel " << myName << endl;
    cerr << "Missing declaration of a Tracker object?" << endl;
    throw parsingException();
  }


  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      parameterNameCopy = parameterName;
      correctlyBroken = breakParameterName(parameterNameCopy, mainIndex, secondaryIndex);
      if (parameterName=="nLayers") {
        nBarrelLayers=atoi(parameterValue.c_str());
      } else if (parameterName=="readout") { // TODO: move this to the types
        if (parameterValue=="cluster") readoutMode=Module::Cluster;
        else if (parameterValue=="binary") readoutMode=Module::Binary;
        else {
          cerr << "Error: unknown readout mode: " << parameterValue << endl;
          cerr << "Allowed values: cluster, binary" << endl;
          throw parsingException();
        }
      } else if (parameterName=="smallDelta") {
        myTracker_->setSmallDelta(atof(parameterValue.c_str()));
      } else if (parameterName=="bigDelta") {
        myTracker_->setBigDelta(atof(parameterValue.c_str()));
      } else if ((correctlyBroken)&&((parameterNameCopy=="smallDelta")||(parameterNameCopy=="bigDelta"))) {
        double deltaValue =  atof(parameterValue.c_str());
        if ((mainIndex==0)||(secondaryIndex!=0)||(deltaValue==0)) {
          cerr << "Error: parameter setting for " << parameterNameCopy << " must have the following structure: "
            << "parameter[ring] = value , with non-zero ring index and non-zero value" << endl;
          throw parsingException();
        } else {
          if (parameterNameCopy=="smallDelta") {
            // put smallDelta [layer] as a special option
            myTracker_->setSpecialSmallDelta(mainIndex, deltaValue);
          } else if (parameterNameCopy=="bigDelta") {
            // put bigDelta [layer] as a special option
            myTracker_->setSpecialBigDelta(mainIndex, deltaValue);
          } else {
            cerr << "BAD error: this should never happen. contact the developers..." << endl;
            throw parsingException();
          }
        }
      } else if (parameterName=="phiSegments") {
        phiSegments=atoi(parameterValue.c_str());
      } else if (parameterName=="shortBarrel") {
        shortBarrel=(parameterValue=="true");
      } else if (parameterName=="sameRods") {
        sameRods=(parameterValue=="true");
      } else if (parameterName=="compress") {
        compress=(parameterValue=="true");
      } else if (parameterName=="minimumZ") {  // deprecated. mantained for backwards compatibility
        shortBarrel=(atof(parameterValue.c_str())>0.);
      } else if (parameterName=="maximumZ") {
        maxZ=atof(parameterValue.c_str());
      } else if (parameterName=="aspectRatio") {
        aspectRatio=atof(parameterValue.c_str());
        aspectRatioManual = true;
        if (aspectRatio<=0) {
          cerr << "ERROR: Parsing barrel \"" << myName
            << "\": wrong aspect ratio (height/width): " << parameterValue
                                                            << " should be a positive number" << endl;
          throw parsingException();
        }
      } else if (parameterName=="size") {
        double widthValue, lengthValue;
        if (sscanf(parameterValue.c_str(), "%lfx%lf", &widthValue, &lengthValue)==2) {
          if ((widthValue>0)&&(lengthValue>0)) {
            size.first=widthValue;
            size.second=lengthValue;
          } else {
            cerr << "ERROR: parsing the module size for barrel " << myName
              << ": \"" << parameterValue.c_str() << "\"I got a negative width or length." << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: Parsing size of modules for barrel \"" << myName
            << "\": unknown/nonsense value \"" << parameterValue << "\". Should be 30x50 if the module is 30mm wide and 50mm long." << endl;
          throw parsingException();
        }
      } else if (parameterName=="nModules") {
        nBarrelModules=atoi(parameterValue.c_str());
      } else if (parameterName=="innerRadius") {
        barrelRhoIn=atof(parameterValue.c_str());
      } else if (parameterName=="outerRadius") {
        barrelRhoOut=atof(parameterValue.c_str());
      } else if (parameterName=="dsDistanceOverride") {
        dsDistanceOverride=atof(parameterValue.c_str());
      } else if (parameterName=="option") {
        char charBuf[100];
        std::string aString;
        double aVal;
        int layerNum;
        bool gotIt=false;
        if (sscanf(parameterValue.c_str(), "%d/%s", &layerNum, charBuf)==2) {
          aString = charBuf;
          // Stacked option, example option = 3/Stacked-30
          // meaning 3rd layer becomes stacked with the inner one at -30 mm w.r.t. the outer
          // Warning: the inner one would be the 3rd layer and the outer one would become the 4th
          // In order to have 2 consecutive stacked layers you should say: 3/Stacked-30 5/Stacked-30
          if (aString.find("Stacked")==0) {
            aVal=atof(aString.substr(7, aString.length()).c_str());
            if (aVal<0) {
              layerOptions[layerNum].first=Layer::Stacked;
              layerOptions[layerNum].second=aVal;
              gotIt=true;
              std::ostringstream tempString;
              tempString.str(""); tempString << "Option: stacked for layer " << layerNum << " at " << aVal;
              logDEBUG(tempString);
            } else {
              cerr << "ERROR: Wrong stack distance: (" << aVal << ") must be negative." << std::endl;
              throw parsingException();
            }
          }
          if (!gotIt) {
            cerr << "ERROR: Parsing layer option for barrel \"" << myName
              << "\': unknown/nonsense option \"" << aString << "\"" << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: Parsing barrel \"" << myName
            << "\": Wrong syntax for a layer option: \"" << parameterValue
            << "\" should be layer/option" << endl;
          throw parsingException();
        }
      } else if (parameterName=="directive") {
        char charBuf[100];
        std::string aString;
        double aVal;
        int layerNum;
        bool gotIt=false;
        if (sscanf(parameterValue.c_str(), "%d/%s", &layerNum, charBuf)==2) {
          aString = charBuf;
          if (aString=="F") {
            layerDirectives[layerNum]=Layer::FIXED;
            gotIt=true;
          } else if (aString=="S") {
            layerDirectives[layerNum]=Layer::SHRINK;
            gotIt=true;
          } else if (aString=="E") {
            layerDirectives[layerNum]=Layer::ENLARGE;
            gotIt=true;
          } else if (aString=="A") {
            layerDirectives[layerNum]=Layer::AUTO;
            gotIt=true;
          }
          if (!gotIt) {
            aVal = atof(aString.c_str());
            if (aVal>0) {
              layerDirectives[layerNum]=aVal;
              gotIt=true;
            }
          }
          if (!gotIt) {
            cerr << "ERROR: Parsing layer directive for barrel \"" << myName
              << "\": unknown/nonsense directive \"" << aString << "\"" << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: Parsing barrel \"" << myName << endl
            << "\": wrong syntax for a layer directive: \"" << parameterValue
            << "\" should be layer/command" << endl;
          throw parsingException();
        }
      } else if (parameterName=="tiltedLayerSpecFile") {
        if (!parseTilted(parameterValue, myName)) {
          cerr << "ERROR: Failure while parsing tilted module parms" << endl;
          throw parsingException();
        }
      } else {
        cerr << "ERROR: Unknown parameter \"" << parameterName << "\"" << endl;
        throw parsingException();
      }
      // cout << "\t" << parameterName << " = " << parameterValue << ";" << endl; // debug
    }
  }


  // Actually creating the barrel if all the mandatory parameters were set
  if (((nBarrelLayers > 1) && (barrelRhoIn != 0) && (barrelRhoOut != 0) && ((nBarrelModules != 0) || (maxZ != 0))  ) || 
      ((nBarrelLayers == 1)&& (barrelRhoIn != 0) && ((nBarrelModules != 0) || (maxZ != 0)) )) { // single-layer barrel disregard the outerRadius parameter and place the layer at the innerRadius


    if ((size.first==0)||(size.second==0)) {
      // cout << "mersidebug: normal module" << endl; // debug
      sampleBarrelModule = new BarrelModule(aspectRatio);   // Square modules of kind rphi
    } else {
      // cout << "mersidebug: normal special size module" << endl; //debug
      if (aspectRatioManual) { // You set the size manually and also specify the aspect ratio!
        // cout << "mersidebug: special size, but also manual aspect ratio" << endl; // debug
        cerr << "ERROR: Parsing barrel \"" << myName << "\" found both the size and the aspect ratio specified!" << endl
          << "Please remove one of the two settings" << endl;
        throw parsingException();
      } else {
        double waferDiameter = pow((pow(size.first, 2) + pow(size.second, 2)), 0.5);
        aspectRatio = size.second/size.first; // heigth/width
        sampleBarrelModule = new BarrelModule(waferDiameter, aspectRatio);
        // cout << "mersidebug: sampleBarrelModule = new BarrelModule(" << waferDiameter << ", " << aspectRatio << ");" <<  endl; // debug
      }
    }

    if (dsDistanceOverride > -1.) myTracker_->setGeometryDsDistance(myName, 0, dsDistanceOverride); // Set the override for dsDistance, so that the values from the Types file are ignored
    // Important: if no directive was given, the following line will clear
    // possible previous directives coming from a different barrel
    myTracker_->setLayerDirectives(layerDirectives);
    myTracker_->setLayerOptions(layerOptions);
    myTracker_->setPhiSegments(phiSegments);

    sampleBarrelModule->setReadoutMode(readoutMode); 
    sampleBarrelModule->setResolutionRphi();
    sampleBarrelModule->setResolutionY();

    if (tiltedBarrelSpecs_[myName].size() == 0) {

      LayerVector myBarrelLayers = myTracker_->buildBarrel(nBarrelLayers,
                                                           barrelRhoIn,
                                                           barrelRhoOut,
                                                           maxZ,
                                                           nBarrelModules,
                                                           sampleBarrelModule,
                                                           myName,
                                                           Layer::NoSection,
                                                           compress,
                                                           shortBarrel, 
                                                           sameRods); // Actually build a compressed barrel (mezzanine or normal)

      delete sampleBarrelModule; // Dispose of the sample module
      if ((maxZ!=0)&&compress)
        myTracker_->compressBarrelLayers(myBarrelLayers, shortBarrel, maxZ);
    } else {
      logDEBUG("Building tilted module barrel " + myName + ", composed of " + any2str(tiltedBarrelSpecs_[myName].size()) + " layers");
      myTracker_->buildTiltedBarrel(myName, tiltedBarrelSpecs_[myName], sampleBarrelModule);

    }

  } else {
    cerr << "ERROR: Missing mandatory parameter for barrel " << myName << endl;
    throw parsingException();
  }

  // Set back the generic small and big deltas
  myTracker_->setSmallDelta(genericSmallDelta);
  myTracker_->setBigDelta(genericBigDelta);
  myTracker_->resetSpecialDeltas();

  return true;
}

// Parses the part of config files between '{' and '}', creating an Endcap
bool configParser::parseEndcap(string myName, istream &inStream) {
  string parameterName;
  string parameterValue;
  string parameterNameCopy;
  int mainIndex, secondaryIndex;
  bool correctlyBroken;

  int readoutMode = Module::Binary;
  int nDisks = 0;       // nDisks (per side)
  int nRings = 0;
  double barrelToEndcap = 0;  // gap between barrel and endcap
  double minZ = 0;            // explicit minimum z
  double maxZ = 0;            // maximum z
  double rhoIn = 0;
  double rhoOut = 0;
  double innerEta = 0;
  double dsDistanceOverride = -1;
  int diskParity = 0;
  std::map<int, int> shapeType;
  std::map<int, bool> explicitShapeType;
  std::map<int, double> aspectRatio;
  std::map<int, std::pair<double, double> > size;
  std::map<int, bool> explicitAspectRatio;
  std::map<int, bool> explicitSize;
  std::map<int, bool> specialRing;
  int phiSegments = 4;
  bool alignEdges = true;
  bool oddSegments = false;

  // Ring 0 represents the default
  shapeType[0] = Module::Wedge;
  explicitShapeType[0] = false;
  aspectRatio[0] = 1.;
  explicitAspectRatio[0] = false;
  explicitSize[0] = false;
  specialRing[0] = true;

  // Fetch the generic Delta of the tracker
  double genericSmallDelta = myTracker_->getSmallDelta();
  double genericBigDelta = myTracker_->getBigDelta();

  map<pair<int, int>, bool> mapDiskRingRemoveToOuter;
  map<pair<int, int>, bool>::iterator mapDiskRingRemoveToOuterIt;

  // Directives (this are communicated to the Tracker object)
  std::map<int, int> ringDirective;
  std::map<int, double> ringGaps;

  // Tracker should be already there
  if (!myTracker_) {
    cerr << "ERROR: tracker is NULL when declaring endcap \"" << myName
      << "\": missing declaration of a Tracker object?" << endl;
    throw parsingException();
  }

  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      parameterNameCopy = parameterName;
      correctlyBroken = breakParameterName(parameterNameCopy, mainIndex, secondaryIndex);
      if (correctlyBroken) {
        if ((mainIndex==0)||(secondaryIndex!=0)) {
          cerr << "ERROR: parameter setting for " << parameterNameCopy << " must have the following structure: "
            << "parameter[ring] = value , with greater-than-zero ring index" << endl;
          throw parsingException();
        }
        if (parameterNameCopy=="smallDelta") {
          double deltaValue =  atof(parameterValue.c_str());
          // Put smallDelta [layer] as a special option
          myTracker_->setSpecialSmallDelta(mainIndex, deltaValue);
        } else if (parameterNameCopy=="bigDelta") {
          double deltaValue =  atof(parameterValue.c_str());
          // Put bigDelta [layer] as a special option
          myTracker_->setSpecialBigDelta(mainIndex, deltaValue);
        } else if (parameterNameCopy=="aspectRatio") {
          double myValue = atof(parameterValue.c_str());
          if (myValue<=0) {
            cerr << "ERROR: Parsing endcap \"" << myName
              << "\": wrong aspect ratio (height/width): " << parameterValue
                                                              << " should be a positive number" << endl;
            throw parsingException();
          }
          aspectRatio[mainIndex] = myValue;
          explicitAspectRatio[mainIndex] = true;
          specialRing[mainIndex] = true;
        } else if (parameterNameCopy=="size") {
          double widthValue, lengthValue;
          if (sscanf(parameterValue.c_str(), "%lfx%lf", &widthValue, &lengthValue)==2) {
            if ((widthValue>0)&&(lengthValue>0)) {
              size[mainIndex].first = widthValue;
              size[mainIndex].second = lengthValue;
              explicitSize[mainIndex] = true;
              explicitShapeType[mainIndex] = true;
              shapeType[mainIndex] = Module::Rectangular;
              specialRing[mainIndex] = true;
            } else {
              cerr << "ERROR: parsing the module size for endcap " << myName
                << ": \"" << parameterValue.c_str() << "\"I got a negative width or length." << endl;
              throw parsingException();
            }
          } else {
            cerr << "ERROR: Parsing size of modules for endcap \"" << myName
              << "\": unknown/nonsense value \"" << parameterValue << "\". Should be 30x50 if the module is 30mm wide and 50mm long." << endl;
            throw parsingException();
          }
        } else if (parameterNameCopy=="shape") {
          if (parameterValue=="wedge") {
            shapeType[mainIndex]=Module::Wedge;
            explicitShapeType[mainIndex] = true;
            specialRing[mainIndex] = true;
          } else if (parameterValue=="rectangular") {
            shapeType[mainIndex]=Module::Rectangular;
            explicitShapeType[mainIndex] = true;
            specialRing[mainIndex] = true;  
          } else {
            cerr << "ERROR: Parsing endcap \"" << myName
              << "\": wrong syntax for a shape: \"" << parameterValue
              << "\" should be \"rectangular\" or \"wedge\"" << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: parameter " << parameterNameCopy << " is not allowed to have"
            << "parameter[ring] = value" << endl;
          throw parsingException();
        }
      } else if (parameterName=="nDisks") {
        nDisks=atoi(parameterValue.c_str());
      } else if (parameterName=="readout") { // TODO: move this to the types
        if (parameterValue=="cluster") readoutMode=Module::Cluster;
        else if (parameterValue=="binary") readoutMode=Module::Binary;
        else {
          cerr << "Error: unknown readout mode: " << parameterValue << endl;
          cerr << "Allowed values: cluster, binary" << endl;
          throw parsingException();
        }
      } else if (parameterName=="smallDelta") {
        myTracker_->setSmallDelta(atof(parameterValue.c_str()));
      } else if (parameterName=="bigDelta") {
        myTracker_->setBigDelta(atof(parameterValue.c_str()));
      } else if (parameterName=="phiSegments") {
        phiSegments=atoi(parameterValue.c_str());
      } else if (parameterName=="alignEdges") {
        if (parameterValue=="true") alignEdges = true;
        else if (parameterValue=="false") alignEdges = false;
        else {
          cerr << "Error in alignEdges. Boolean value (true|false) expected. Value '"
            << parameterValue << "' was found" << std::endl;
          throw parsingException();
        }
      } else if (parameterName=="oddSegments") {
        if (parameterValue=="true") oddSegments = true;
        else if (parameterValue=="false") oddSegments = false;
        else {
          cerr << "Error in oddSegments. Boolean value (true|false) expected. Value '"
            << parameterValue << "' was found" << std::endl;
          throw parsingException();
        }
      } else if (parameterName=="innerEta") {
        innerEta=atof(parameterValue.c_str());
      } else if (parameterName=="innerRadius") {
        rhoIn=atof(parameterValue.c_str());
      } else if (parameterName=="outerRadius") {
        rhoOut=atof(parameterValue.c_str());
      } else if (parameterName=="nRings") {
        nRings=atoi(parameterValue.c_str());
      } else if (parameterName=="barrelGap") {
        barrelToEndcap=atof(parameterValue.c_str());
      } else if (parameterName=="minimumZ") {
        minZ=atof(parameterValue.c_str());
      } else if (parameterName=="maximumZ") {
        maxZ=atof(parameterValue.c_str());
      } else if (parameterName=="dsDistanceOverride") {
        dsDistanceOverride = atof(parameterValue.c_str());
      } else if (parameterName=="diskParity") {
        diskParity=atoi(parameterValue.c_str());
        if (diskParity>0) diskParity=1; else diskParity=-1;
      } else if (parameterName=="aspectRatio") { // Can also be ring-dependent
        aspectRatio[0]=atof(parameterValue.c_str());
        explicitAspectRatio[0] = true;
        if (aspectRatio[0]<=0) {
          cerr << "ERROR: Parsing endcap \"" << myName
            << "\": wrong aspect ratio (height/width): " << parameterValue
                                                            << " should be a positive number" << endl;
          throw parsingException();
        }
      } else if (parameterName=="size") { // Can also be ring-dependent
        double widthValue, lengthValue;
        if (sscanf(parameterValue.c_str(), "%lfx%lf", &widthValue, &lengthValue)==2) {
          if ((widthValue>0)&&(lengthValue>0)) {
            size[0].first = widthValue;
            size[0].second = lengthValue;
            explicitSize[0] = true;
            explicitShapeType[0] = true;
            shapeType[0] = Module::Rectangular;
          } else {
            cerr << "ERROR: parsing the module size for endcap " << myName
              << ": \"" << parameterValue.c_str() << "\"I got a negative width or length." << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: Parsing size of modules for endcap \"" << myName
            << "\": unknown/nonsense value \"" << parameterValue << "\". Should be 30x50 if the module is 30mm wide and 50mm long." << endl;
          throw parsingException();
        }
      } else if (parameterName=="shape") { // Can also be ring-dependent
        bool syntaxOk=true;
        if (parameterValue=="wedge") {
          shapeType[0]=Module::Wedge;
          explicitShapeType[0] = true;
        } else if (parameterValue=="rectangular") {
          shapeType[0]=Module::Rectangular;
          explicitShapeType[0] = true;
        } else {
          syntaxOk=false;
        }
        if (!syntaxOk) {
          cerr << "ERROR: Parsing endcap \"" << myName
            << "\": wrong syntax for a shape: \"" << parameterValue
            << "\" should be \"rectangular\" or \"wedge\"" << endl;
          throw parsingException();
        }
      } else if (parameterName=="directive") {
        int ringNum; int increment;
        if (sscanf(parameterValue.c_str(), "%d%d", &ringNum, &increment)==2) {
          if (increment!=0) ringDirective[ringNum]=increment;
        } else {
          cerr << "ERROR: Parsing endcap \"" << myName
            << "\": wrong syntax for a ring directive: \"" << parameterValue
            << "\" should be ring+increment or ring-decrement )" << endl;
          throw parsingException();
        }
      } else if (parameterName=="ringGap") {
        int ringNum; float gap;
        if (sscanf(parameterValue.c_str(), "%d%f", &ringNum, &gap)==2) {
          if (gap!=0) ringGaps[ringNum]=gap;
        } else {
          cerr << "ERROR: Parsing endcap \"" << myName
               << "\": wrong syntax for a ring gap: \"" << parameterValue
               << "\" should be ring+gap or ring-gap )" << endl;
          throw parsingException();
        }
      } else if (parameterName=="removeRings") {
        // Example of the syntax "D1R4+"
        int diskNum, ringNum;
        char plusMinus;
        bool parseError=false;
        pair <int, int> diskRing;
        if (sscanf(parameterValue.c_str(), "D%dR%d%c", &diskNum, &ringNum, &plusMinus)==3) {
          if ((plusMinus=='+')||(plusMinus=='-')) {
            diskRing.first=diskNum;
            diskRing.second=ringNum;
            bool outer;
            outer=(plusMinus=='+');
            mapDiskRingRemoveToOuter[diskRing]=outer;
          } else {
            parseError=true;
          }
        } else {
          parseError=true;
        }

        if (parseError) {
          cerr << "ERROR: Parsing endcap \"" << myName
            << "\": wrong syntax for a ring directive: \"" << parameterValue
            << "\" should be like D1R4+ (to remove rings 4 or more from disk 1)" << endl;
          throw parsingException();
        }
      } else {
        cerr << "ERROR: Unknown parameter \"" << parameterName << "\"" << endl;
        throw parsingException();
      }
      // cout << "\t" << parameterName << " = " << parameterValue << ";" << endl; // debug
    }
  }

  // Check the consistency in oddSegments and phiSegments
  if (((phiSegments%2)!=0)&&(oddSegments)) {
    cerr << "ERROR: inconsistent requirements in endcap \"" << myName
      << "\": oddSegments was set to true and the number of segments was"
      << " set to " << oddSegments << ". This makes it possible to have an even number"
      << " of modules in a ring, so the first and last one could clash." << endl;
    throw parsingException();
  }

  // Chose the proper endcap minimum Z according to the parsed parameters
  if ((innerEta!=0)&&(rhoIn!=0)) {
    cerr << "ERROR: Parsing endcap \"" << myName
      << "\": I got both innerEta = " << innerEta << " and innerRadius = " << rhoIn << endl
      << "This is inconsistent: please choose one of the two methods to place the endcap" << endl;
    throw parsingException();
  }

  // Chose the proper endcap minimum Z according to the parsed parameters
  if ((minZ!=0)&&(barrelToEndcap!=0)) {
    cerr << "Parsing endcap \"" << myName
      << "\": I got both minimumZ = " << minZ << " and barrelGap = " << barrelToEndcap << endl
      << "This is inconsistent: please choose one of the two methods to place the endcap" << endl;
    throw parsingException();
  }

  if (minZ==0) {
    minZ=myTracker_->getMaxBarrelZ(+1)+barrelToEndcap;
  }


  for (int iRing=0; iRing<1; ++iRing) {
    // Check consistency in module shape type assigned
    if (explicitShapeType[iRing]) {
      if (shapeType[iRing]==Module::Wedge) { 
        // It's wedge-shaped
        if ((explicitAspectRatio[iRing]) || (explicitSize[iRing])) {
          // And the aspect ratio or size was explicitly given!
          // This is an error...
          cerr << "Parsing endcap \"" << myName << "\"";
          if (iRing!=0) cerr << " (ring " << iRing << ")";
          cerr << ": I see module shape 'wedge' and aspect ratio or size assigned. This is inconsistent. I quit." << endl;
          throw parsingException();
        }
      }
    }
    // TODO: check also that aspectratio and size are not given at the same time
  }


  if (dsDistanceOverride > -1.) myTracker_->setGeometryDsDistance(myName, 0, dsDistanceOverride); // Set the override for dsDistance, so that the values from the Types file are ignored

  // Actually creating the endcap if all the mandatory parameters were set
  if ( (nDisks != 0) &&
       ((rhoIn != 0)||(innerEta!=0)||(nRings!=0)) &&
       (rhoOut != 0) &&
       (minZ != 0) &&
       (maxZ != 0) &&
       (diskParity != 0)) {

    std::map<int, EndcapModule*> sampleModule;


    std::string mess;
    for (std::map<int, bool>::iterator it = specialRing.begin(); it != specialRing.end(); ++it) {

      const int& iRing = it->first;
      EndcapModule* aModule;

      // The same old sample module
      if (shapeType[iRing]==Module::Wedge) {
        aModule = new EndcapModule(Module::Wedge);
        mess = "wedge";
      } else if (shapeType[iRing]==Module::Rectangular) {
        mess = "rect_";
        if (explicitSize[iRing]) {
          mess += "size";
          // TODO: create a sample module with the right size here
          // instead of the really big one
          double waferDiameter = pow((pow(size[iRing].first, 2) + pow(size[iRing].second, 2)), 0.5);
          double aspectRatio = size[iRing].second / size[iRing].first; // height / width
          aModule = new EndcapModule(waferDiameter, aspectRatio);
        } else {
          mess += "square";
          aModule = new EndcapModule(aspectRatio[iRing]);
        }
      } else {
        cerr << "ERROR: an unknown module shape type was generated inside the configuration parser."
          << "this should never happen!" << std::endl;
        throw parsingException();
      }

      aModule->setReadoutMode(readoutMode);
      aModule->setResolutionRphi();
      aModule->setResolutionY();
      sampleModule[iRing] = aModule;

      // std::cerr << "Module in ring " << iRing << " is " << mess << std::endl; // debug
    }

    // Important: if no directive was given, the following line will clear
    // possible previous directives coming from a different endcap
    myTracker_->setRingDirectives(ringDirective);
    myTracker_->setPhiSegments(phiSegments);
    myTracker_->setRingGaps(ringGaps);


    if (rhoIn!=0 || nRings != 0) {
      myTracker_->buildEndcaps(nDisks,     // nDisks (per side)
                               nRings,
                               minZ,  // minZ
                               maxZ,
                               rhoIn,
                               rhoOut,
                               sampleModule, // TODO: give the full map of sample modules
                               myName,
                               diskParity,
                               oddSegments, alignEdges,
                               Layer::NoSection);
    } else {
      myTracker_->buildEndcapsAtEta(nDisks,     // nDisks (per side)
                                    nRings,
                                    minZ,  // minZ
                                    maxZ,
                                    innerEta,
                                    rhoOut,
                                    sampleModule, // TODO: give the full map of sample modules
                                    myName,
                                    diskParity,
                                    oddSegments, alignEdges,
                                    Layer::NoSection);
    }

    for (std::map<int, EndcapModule*>::iterator it = sampleModule.begin();
         it != sampleModule.end(); ++it) {
      delete it->second;
    }
    sampleModule.clear();

    // Remove the deleted rings

    for (mapDiskRingRemoveToOuterIt=mapDiskRingRemoveToOuter.begin();
         mapDiskRingRemoveToOuterIt!=mapDiskRingRemoveToOuter.end();
         mapDiskRingRemoveToOuterIt++) {
      int iDisk = ((*mapDiskRingRemoveToOuterIt).first).first;
      int iRing = ((*mapDiskRingRemoveToOuterIt).first).second;
      bool directionOuter = (*mapDiskRingRemoveToOuterIt).second;

      // cerr << "myTracker_->removeDiskRings("<<myName<<","<<iDisk<<", "<<iRing<<", "<<directionOuter<<");" << endl; // debug
      myTracker_->removeDiskRings(myName, iDisk, iRing, directionOuter);
    }

  } else {
    cerr << "ERROR: Missing mandatory parameter for endcap " << myName << endl ;
    cerr << "Mandatory parameters are: nDisks " << endl
      << "                          [ innerRadius | innerEta | nRings ]" << endl
      << "                          outerRadius" << endl
      << "                          [ minimumZ | barrelGap ]" << endl
      << "                          diskParity" << endl;
    throw parsingException();
  }

  // Set back the generic small and big deltas
  myTracker_->setSmallDelta(genericSmallDelta);
  myTracker_->setBigDelta(genericBigDelta);
  myTracker_->resetSpecialDeltas();

  return true;
}

bool configParser::parsePixels(string myName, istream &inStream) {
  string parameterName, parameterValue;
  double waferDiameter, aspectRatio;

  double etaCut = 0.0;
  double rInBar = 0.0;
  double rInEnd = 0.0;
  double rOutBar = 0.0;
  double rOutEnd = 0.0;
  double barrelToEndcap = 0.0;
  double maxZ = 0.0;
  int nLayers = 0;
  int nModules = 0;
  int nDisks = 0;
  int nRings = 0;
  int phiSegments = 4;
  int diskParity = 0;
  double smallDelta = 0.0;
  double bigDelta = 3;
  double dsDistanceOverride = -1;
  bool sameRods = false;
  pair<double, double> bmsize, emsize; // width, length
  map<int, double> layerDirectives;
  map<int, LayerOption> layerOptions;
  map<int, int> ringDirectives;

  // Tracker is a singleton. Declare just one, please
  if (myTracker_) {
    cerr << "ERROR: tracker is not NULL when declaring tracker " << myName << endl;
    cerr << "Double declaration of a Tracker object?" << endl;
    throw parsingException();
  }
  myTracker_ = new Tracker(myName);
  myTracker_->setSmallDelta(smallDelta);
  myTracker_->setBigDelta(bigDelta);

  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      if (parameterName == "etaCut") {
        etaCut = atof(parameterValue.c_str());
      } else if (parameterName == "nLayers") {
        nLayers = atoi(parameterValue.c_str());
      } else if (parameterName == "nModules") {
        nModules = atoi(parameterValue.c_str());
      } else if (parameterName == "nDisks") {
        nDisks = atoi(parameterValue.c_str());
      } else if (parameterName == "nRings") {
        nRings = atoi(parameterValue.c_str());
      } else if (parameterName == "innerRadiusBarrel") {
        rInBar = atof(parameterValue.c_str());
      } else if (parameterName == "innerRadiusEndcap") {
        rInEnd = atof(parameterValue.c_str());
      } else if (parameterName == "outerRadiusBarrel") {
        rOutBar = atof(parameterValue.c_str());
      } else if (parameterName == "outerRadiusEndcap") {
        rOutEnd = atof(parameterValue.c_str());
      } else if (parameterName == "phiSegments") {
        phiSegments = atoi(parameterValue.c_str());
      } else if (parameterName == "barrelGap") {
        barrelToEndcap = atof(parameterValue.c_str());
      } else if (parameterName == "sameRods") {
        sameRods = (parameterValue == "true");
      } else if (parameterName == "dsDistanceOverride") {
        dsDistanceOverride = atof(parameterValue.c_str());
      } else if (parameterName == "maximumZ") {
        maxZ = atof(parameterValue.c_str());
      } else if (parameterName == "diskParity") {
        diskParity = atoi(parameterValue.c_str());
        (diskParity > 0) ? diskParity = 1 : diskParity = -1;
      } else if (parameterName == "bmSize") {
        double widthValue, lengthValue;
        if (sscanf(parameterValue.c_str(), "%lfx%lf", &widthValue, &lengthValue)==2) {
          if ((widthValue>0)&&(lengthValue>0)) {
            bmsize.first=widthValue;
            bmsize.second=lengthValue;
          } else {
            cerr << "ERROR: parsing the module size for pixel detector " << myName
              << ": \"" << parameterValue.c_str() << "\"I got a negative width or length." << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: Parsing size of barrel modules for pixel detector \"" << myName
            << "\": unknown/nonsense value \"" << parameterValue << "\". Should be 30x50 if the module is 30mm wide and 50mm long." << endl;
          throw parsingException();
        }
      } else if (parameterName == "emSize") {
        double widthValue, lengthValue;
        if (sscanf(parameterValue.c_str(), "%lfx%lf", &widthValue, &lengthValue)==2) {
          if ((widthValue>0)&&(lengthValue>0)) {
            emsize.first=widthValue;
            emsize.second=lengthValue;
          } else {
            cerr << "ERROR: parsing the endcap module size for pixel detector " << myName
              << ": \"" << parameterValue.c_str() << "\"I got a negative width or length." << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: Parsing size of endcap modules for pixel detector \"" << myName
            << "\": unknown/nonsense value \"" << parameterValue << "\". Should be 30x50 if the module is 30mm wide and 50mm long." << endl;
          throw parsingException();
        }
      } else if ( parameterName == "option") {
        char charBuf[100];
        std::string aString;
        double aVal;
        int layerNum;
        bool gotIt=false;
        if (sscanf(parameterValue.c_str(), "%d/%s", &layerNum, charBuf)==2) {
          aString = charBuf;
          // Stacked option, example option = 3/Stacked-30
          // meaning 3rd layer becomes stacked with the inner one at -30 mm w.r.t. the outer
          // Warning: the inner one would be the 3rd layer and the outer one would become the 4th
          // In order to have 2 consecutive stacked layers you should say: 3/Stacked-30 5/Stacked-30
          if (aString.find("Stacked")==0) {
            aVal=atof(aString.substr(7, aString.length()).c_str());
            if (aVal<0) {
              layerOptions[layerNum].first=Layer::Stacked;
              layerOptions[layerNum].second=aVal;
              gotIt=true;
              std::ostringstream tempString;
              tempString.str(""); tempString << "Option: stacked for layer " << layerNum << " at " << aVal;
              logDEBUG(tempString);
            } else {
              cerr << "ERROR: Wrong stack distance: (" << aVal << ") must be negative." << std::endl;
              throw parsingException();
            }
          }
          if (!gotIt) {
            cerr << "ERROR: Parsing layer option for barrel \"" << myName
              << "\': unknown/nonsense option \"" << aString << "\"" << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: Parsing barrel \"" << myName
            << "\": Wrong syntax for a layer option: \"" << parameterValue
            << "\" should be layer/option" << endl;
          throw parsingException();
        }
      } else if ( parameterName == "directive") {
        if (parameterValue.find('/') != string::npos) { // barrel directive
          char charBuf[100];
          std::string aString;
          double aVal;
          int layerNum;
          bool gotIt=false;
          if (sscanf(parameterValue.c_str(), "%d/%s", &layerNum, charBuf)==2) {
            aString = charBuf;
            if (aString=="F") {
              layerDirectives[layerNum]=Layer::FIXED;
              gotIt=true;
            } else if (aString=="S") {
              layerDirectives[layerNum]=Layer::SHRINK;
              gotIt=true;
            } else if (aString=="E") {
              layerDirectives[layerNum]=Layer::ENLARGE;
              gotIt=true;
            } else if (aString=="A") {
              layerDirectives[layerNum]=Layer::AUTO;
              gotIt=true;
            }
            if (!gotIt) {
              aVal = atof(aString.c_str());
              if (aVal>0) {
                layerDirectives[layerNum]=aVal;
                gotIt=true;
              }
            }
            if (!gotIt) {
              cerr << "ERROR: Parsing layer directive for pixel barrel in \"" << myName
                << "\": unknown/nonsense directive \"" << aString << "\"" << endl;
              throw parsingException();
            }
          } else {
            cerr << "ERROR: Parsing pixel barrel in \"" << myName << endl
              << "\": wrong syntax for a layer directive: \"" << parameterValue
              << "\" should be layer/command" << endl;
            throw parsingException();
          }
        } else if (parameterValue.find('+') != string::npos) { // endcap directive
          int ringNum; int increment;
          if (sscanf(parameterValue.c_str(), "%d%d", &ringNum, &increment)==2) {
            if (increment!=0) ringDirectives[ringNum]=increment;
          } else {
            cerr << "ERROR: Parsing pixel endcap in \"" << myName
              << "\": wrong syntax for a ring directive: \"" << parameterValue
              << "\" should be ring+increment or ring-decrement )" << endl;
            throw parsingException();
          }
        } else {
          cerr << "ERROR: Unable to parse directive for pixel detector \"" << myName
            << "\": unknown directive \"" << parameterValue << "\"" << endl;
          throw parsingException();
        }
      } else {
        cerr << "ERROR: Unknown parameter name: " << parameterName << endl;
        throw parsingException();
      }
    }
  }
  if (dsDistanceOverride > -1.) myTracker_->setGeometryDsDistance(myName, 0, dsDistanceOverride);
  myTracker_->setEtaCut(etaCut);
  myTracker_->setPhiSegments(phiSegments);
  myTracker_->setLayerDirectives(layerDirectives);
  myTracker_->setLayerOptions(layerOptions);
  myTracker_->setRingDirectives(ringDirectives);

  waferDiameter = pow((pow(bmsize.first, 2) + pow(bmsize.second, 2)), 0.5);
  aspectRatio = bmsize.second / bmsize.first; // height / width
  BarrelModule* sampleBarrelModule = new BarrelModule(waferDiameter, aspectRatio);
  sampleBarrelModule->setResolutionRphi();
  sampleBarrelModule->setResolutionY();

  myTracker_->buildBarrel(nLayers, rInBar, rOutBar, 0, nModules, sampleBarrelModule, myName, Layer::NoSection, true, false, sameRods);
  delete sampleBarrelModule;

  if (nDisks > 0) {
    waferDiameter = pow((pow(emsize.first, 2) + pow(emsize.second, 2)), 0.5);
    aspectRatio = emsize.second / emsize.first; // height / width
    std::map<int, EndcapModule*> sampleEndcapModule;
    sampleEndcapModule[0] = new EndcapModule(waferDiameter, aspectRatio);
    sampleEndcapModule[0]->setResolutionRphi();
    sampleEndcapModule[0]->setResolutionY();
    myTracker_->buildEndcaps(nDisks, nRings, myTracker_->getMaxBarrelZ(+1) + barrelToEndcap,
                             maxZ, rInEnd, rOutEnd, sampleEndcapModule, myName, diskParity, false, false, Layer::NoSection);
    delete sampleEndcapModule[0];
  }

  return true;
}

bool configParser::breakParameterName(string& parameterName, int& ringIndex, int& diskIndex) {
  bool result = false;
  string tempParamStr;
  string tempParamInd;
  string dummy;
  int nIndexes;
  istringstream parameterStream(parameterName);
  diskIndex=0;

  getline(parameterStream, tempParamStr, '[');
  getline(parameterStream, tempParamInd, ']');
  getline(parameterStream, dummy);

  // We check that we found parameter[index] and nothing else
  if ((tempParamStr!="")&&(tempParamInd!="")&&(dummy=="")&&(parameterStream.eof())) {
    nIndexes=sscanf(tempParamInd.c_str(), "%d,%d", &ringIndex, &diskIndex);
    if (nIndexes==2) {
      result=true;
    } else if (nIndexes==1) {
      result=true;
      diskIndex=0;
    }
    parameterName=tempParamStr;
  }

  return result;
}

bool configParser::breakParameterName(string& parameterName, string& stringIndex) {
  bool result = false;
  string tempParamStr;
  string tempParamInd;
  string dummy;
  istringstream parameterStream(parameterName);
  stringIndex = "";

  getline(parameterStream, tempParamStr, '[');
  getline(parameterStream, tempParamInd, ']');
  getline(parameterStream, dummy);

  // We check that we found parameter[index] and nothing else
  if ((tempParamStr!="")&&(tempParamInd!="")&&(dummy=="")&&(parameterStream.eof())) {
    parameterName=tempParamStr;
    stringIndex = tempParamInd;
    result = true;
  }

  return result;
}


// The following two functions are temporary: by now we have no difference in syntax
// between the barrel and the endcap
bool configParser::parseBarrelType(string myName, istream& inStream) {
  return parseAnyType(myName, inStream);
}

bool configParser::parseEndcapType(string myName, istream& inStream) {
  return parseAnyType(myName, inStream);
}

bool configParser::parsePixelType(istream& inStream) {
  string parameterName, parameterValue;

  string type;
  int nSides = 0;
  int nStripsAcross = 0;
  int nSegments = 0;
  double dsDistance = 0.0;
  double dsRotation = 0.0;
  double resolutionRphi = -1.0;
  double resolutionY = -1.0;

  // Tracker should be already there
  if (!myTracker_) {
    cerr << "ERROR: pixel detector is NULL when trying to assign types to its modules." << endl;
    return false;
  }

  // parse parameters
  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      if (parameterName == "type") {
        type = parameterValue;
      } else if (parameterName == "nSides") {
        nSides = atoi(parameterValue.c_str());
      } else if (parameterName == "nStripsAcross") {
        nStripsAcross = atoi(parameterValue.c_str());
      } else if (parameterName == "nSegments") {
        nSegments = atoi(parameterValue.c_str());
      } else if (parameterName == "dsDistance") {
        dsDistance = atof(parameterValue.c_str());
      } else if (parameterName == "dsRotation") {
        dsRotation = atof(parameterValue.c_str());
      } else if (parameterName == "resolutionRphi") {
        resolutionRphi = atof(parameterValue.c_str());
      } else if (parameterName == "resolutionY") {
        resolutionY = atof(parameterValue.c_str());
      } else {
        cerr << "ERROR: Unknown parameter in PixelType block: " << parameterName << ". Ignoring it." << endl;
      }
    }
  }

  // dress barrel
  vector<Module*>::iterator miter, mguard;
  vector<Layer*>::iterator liter, lguard = myTracker_->getBarrelLayers()->end();
  for (liter = myTracker_->getBarrelLayers()->begin(); liter != lguard; liter++) {
    mguard = (*liter)->getModuleVector()->end();
    for (miter = (*liter)->getModuleVector()->begin(); miter != mguard; miter++) {
      (*miter)->setType(type);
      (*miter)->setNFaces(nSides);
      (*miter)->setNStripsAcross(nStripsAcross);
      (*miter)->setNSegments(nSegments);
      (*miter)->setStereoDistance(dsDistance);
      (*miter)->setStereoRotation(dsRotation);
      (*miter)->setResolutionRphi(resolutionRphi);
      (*miter)->setResolutionY(resolutionY);
    }
  }

  // dress endcap
  lguard = myTracker_->getEndcapLayers()->end();
  for (liter = myTracker_->getEndcapLayers()->begin(); liter != lguard; liter++) {
    mguard = (*liter)->getModuleVector()->end();
    for (miter = (*liter)->getModuleVector()->begin(); miter != mguard; miter++) {
      (*miter)->setType(type);
      (*miter)->setNFaces(nSides);
      (*miter)->setNStripsAcross(nStripsAcross);
      (*miter)->setNSegments(nSegments);
      (*miter)->setStereoDistance(dsDistance);
      (*miter)->setStereoRotation(dsRotation);
      (*miter)->setResolutionRphi(resolutionRphi);
      (*miter)->setResolutionY(resolutionY);
    }
  }

  return true;
}

// TODO: put the code into parseBarrelType and parseEndcapType when they will actually
// contain different parameters
bool configParser::parseAnyType(string myName, istream& inStream) {
  string parameterName;
  string parameterValue;
  int mainIndex; // Used for the ring and layer
  int secondaryIndex; // Used for the disk

  map<int, int> nStripsAcross;
  map<int, int> nROCRows;
  map<int, int> nSides;
  map<int, int> nSegments;
  map<int, string> type;
  map<int, double> dsDistance;
  map<int, int> triggerWindow;
  map<int, double> dsRotation;
  map<int, int> divideBack;
  map<int, double> xResolution;
  map<int, double> yResolution;

  pair<int, int> specialIndex; // used to indicate ring,disk

  map<pair<int, int>, int> nStripsAcrossSecond;
  map<pair<int, int>, int> nROCRowsSecond;
  map<pair<int, int>, int> nSidesSecond;
  map<pair<int, int>, int> nSegmentsSecond;
  map<pair<int, int>, string> typeSecond;
  map<pair<int, int>, double> dsDistanceSecond;
  map<pair<int, int>, int> triggerWindowSecond;
  map<pair<int, int>, double> dsRotationSecond;
  map<pair<int, int>, int> divideBackSecond;
  map<pair<int, int>, bool> specialSecond;

  // Tracker should be already there
  if (!myTracker_) {
    cerr << "ERROR: tracker is NULL when trying to assign types to a barrel " << myName << endl;
    return false;
  }

  // TODO: decide whether to use nStripAcross or nStripsAcross
  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      if (!breakParameterName(parameterName, mainIndex, secondaryIndex)) {
        cerr << "ERROR: parameter name " << parameterName << " must have the following structure: "
          << "parameter[ring] or parameter[layer] or parameter[ring,disk] e.g. nModules[2]" << endl;
        throw parsingException();
      } else {
        if (secondaryIndex==0) { // Standard assignment per ring or per layer
          if (parameterName=="nStripsAcross") {
            nStripsAcross[mainIndex]=atoi(parameterValue.c_str());
          } else if (parameterName=="nROCRows") {
            nROCRows[mainIndex]=atoi(parameterValue.c_str());
          } else if (parameterName=="nSides") {
            nSides[mainIndex]=atoi(parameterValue.c_str());
          } else if (parameterName=="nSegments") {
            nSegments[mainIndex]=atoi(parameterValue.c_str());
          } else if (parameterName=="type") {
            type[mainIndex]=parameterValue.c_str();
          } else if (parameterName == "dsDistance") {
            dsDistance[mainIndex]=atof(parameterValue.c_str());
          } else if (parameterName == "triggerWindow") {
            triggerWindow[mainIndex]=atoi(parameterValue.c_str());
          } else if (parameterName == "dsRotation") {
            dsRotation[mainIndex]=atof(parameterValue.c_str());
          } else if (parameterName == "divideBack") {
            divideBack[mainIndex]=atoi(parameterValue.c_str());
          } else if (parameterName == "xResolution") {
            xResolution[mainIndex]=atof(parameterValue.c_str());
          } else if (parameterName == "yResolution") {
            yResolution[mainIndex]=atof(parameterValue.c_str());
          }
          // cout << "\t" << parameterName << "[" << mainIndex << "] = " << parameterValue << ";" << endl; // debug
        } else { // Special assignment per disk/ring
          specialIndex.first = mainIndex;
          specialIndex.second = secondaryIndex;
          bool isSpecial = false;
          if (parameterName=="nStripsAcross") {
            if (atoi(parameterValue.c_str())!=nStripsAcross[mainIndex]) {
              nStripsAcrossSecond[specialIndex]=atoi(parameterValue.c_str());
              // It is "special" only if it differs from the default values
              specialSecond[specialIndex]=true;
              isSpecial = true;
            }
          } else if (parameterName=="nROCRows") {
            if (atoi(parameterValue.c_str())!=nROCRows[mainIndex]) {
              nROCRowsSecond[specialIndex]=atoi(parameterValue.c_str());
              specialSecond[specialIndex]=true;
              isSpecial = true;
            }
          } else if (parameterName=="nSides") {
            if (atoi(parameterValue.c_str())!=nSides[mainIndex]) {
              nSidesSecond[specialIndex]=atoi(parameterValue.c_str());
              specialSecond[specialIndex]=true;
              isSpecial = true;
            }
          } else if (parameterName=="nSegments") {
            if (atoi(parameterValue.c_str())!=nSegments[mainIndex]) {
              nSegmentsSecond[specialIndex]=atoi(parameterValue.c_str());
              specialSecond[specialIndex]=true;
              isSpecial = true;
            }
          } else if (parameterName=="type") {
            if (parameterValue!=type[mainIndex]) {
              typeSecond[specialIndex]=parameterValue.c_str();
              specialSecond[specialIndex]=true;
              isSpecial = true;
            }
          } else if (parameterName=="dsDistance") {
            if (atof(parameterValue.c_str())!=dsDistance[mainIndex]) { // TODO: BUG! if you put atoi it will likely crash in the logger
              dsDistanceSecond[specialIndex]=atof(parameterValue.c_str());
              specialSecond[specialIndex]=true;
              isSpecial = true;
            }
          } else if (parameterName=="triggerWindow") {
            if (atoi(parameterValue.c_str())!=triggerWindow[mainIndex]) {
              triggerWindowSecond[specialIndex]=atoi(parameterValue.c_str());
              specialSecond[specialIndex]=true;
              isSpecial = true;
            }
          } else if (parameterName=="dsRotation") {
            if (atof(parameterValue.c_str())!=dsRotation[mainIndex]) {
              dsRotationSecond[specialIndex]=atof(parameterValue.c_str());
              specialSecond[specialIndex]=true;
              isSpecial = true;
            }
          } else if (parameterName=="divideBack") {
            if (atoi(parameterValue.c_str())!=divideBack[mainIndex]) {
              divideBackSecond[specialIndex]=atoi(parameterValue.c_str());
              specialSecond[specialIndex]=true;
              isSpecial = true;
            }
          }
          if (!isSpecial) {
            std::ostringstream tempString;
            tempString.str(""); tempString << "The special parameter "
              << parameterName << "[" << mainIndex << "," << secondaryIndex << "] is setting the same "
              << "values as the default parameter "
              << parameterName << "[" << mainIndex << "]. Ignoring it.";
            logWARNING(tempString);
          } // else { // debug
          // cout << "\t" << parameterName << "[" << mainIndex << ","<<secondaryIndex<<"] = " << parameterValue << ";" << endl; // debug
          // }
        }
      }
    }
  }

  myTracker_->setModuleTypes(myName,
                             nStripsAcross, nROCRows, nSides, nSegments, type, dsDistance, triggerWindow, dsRotation, divideBack,
                             xResolution, yResolution,
                             nStripsAcrossSecond, nROCRowsSecond, nSidesSecond, nSegmentsSecond, typeSecond,
                             dsDistanceSecond, triggerWindowSecond, dsRotationSecond, divideBackSecond, specialSecond);

  return true;

}

// Output stuff
bool configParser::parseOutput(istream& inStream) {
  string parameterName;
  string parameterValue;

  string outPath="";

  // Tracker should be already there
  if (!myTracker_) {
    cerr << "ERROR: tracker is NULL when trying to assign output options" << endl;
    return false;
  }

  // TODO: decide whether to use nStripAcross or nStripsAcross
  while (!inStream.eof()) {
    while (parseParameter(parameterName, parameterValue, inStream)) {
      if (parameterName=="Path") {
        outPath=parameterValue;
      } else {
        cerr << "ERROR: While parting output parameters, I got an unrecognized parameter: " << parameterName << endl;
        throw parsingException();
      }
    }
  }

  if (outPath!="") myTracker_->setActiveDirectory(outPath);

  return true;

}

bool configParser::parseSupportParameters(std::istream& inStream, std::list<std::pair<int, double> >& plist) {
  std::string name, value;
  int index, dummy;
  double mid_z;
  while (!inStream.eof()) {
    while (parseParameter(name, value, inStream)) {
      std::pair<int, double> p;
      if (name.compare("midZ") == 0) index = 0;
      else {
        if (name.compare("midZ") > 0) {
          if (!breakParameterName(name, index, dummy)) return false;
        }
        else return false;
      }
      p.first = index;
      mid_z = atof(value.c_str());
      p.second = mid_z;
      if (mid_z >= 0) plist.push_back(p);
    }
  }
  return true;
}

// Takes the type and the name of the object, creates a strstream of
// everything between '{' and '}' and passes it to the parser corresponding
// to the declared Type. (parseTracker, parseBarrel, parseEndcap)
bool configParser::parseObjectType(string myType) {
  string str;
  string typeConfig;

  if (myType=="Tracker") {
    // cout << "Reading tracker main parameters and name: "; // debug
    str=getTill(configFile_, '{', true);
    if (str!="") {
      // cout << str << endl; // debug
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
        istringstream typeStream(typeConfig);
        parseTracker(str, typeStream);
      }
    }
  } else if (myType=="Barrel") {
    startTaskClock("Creating barrel");
    str=getTill(configFile_, '{', true);
    if (str!="") {
      addTaskInfo(str);
      //cout << str << endl;
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
        istringstream typeStream(typeConfig);
        parseBarrel(str, typeStream);
      }
    }
    stopTaskClock();
  } else if (myType=="Endcap") {
    startTaskClock("Creating endcap");
    str=getTill(configFile_, '{', true);
    if (str!="") {
      addTaskInfo(str);
      //cout << str << endl;
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
        istringstream typeStream(typeConfig);
        parseEndcap(str, typeStream);
      }
    }
    stopTaskClock();
  } else if (myType=="Support") {
    getTill(configFile_, '}', false, true);
  } else if (myType=="Pixels") {
    getTill(configFile_, '}', false, true);
  } else {
    ostringstream tempSS; tempSS << "Unknown piece of tracker " << myType;
    logERROR(tempSS);
    return false;
  }

  return true;
}

// Takes the type and the name of the object to be dressed, creates a strstream of
// everything between '{' and '}' and passes it to the parser corresponding
// to the declared Type. (parseBarrelType, parseEndcapType)
bool configParser::parseDressType(string myType) {
  string str;
  string typeConfig;

  if (myType=="BarrelType") {
    startTaskClock("Assigning module types to barrel");
    str=getTill(configFile_, '{', true);
    if (str!="") {
      //cout << str << endl;
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
        istringstream typeStream(typeConfig);
        parseBarrelType(str, typeStream);
      }
    }
    stopTaskClock();
  } else if (myType=="EndcapType") {
    startTaskClock("Assigning module types to endcap");
    str=getTill(configFile_, '{', true);
    if (str!="") {
      typeConfig=getTill(configFile_, '}', false);
      if (typeConfig!="") {
        istringstream typeStream(typeConfig);
        parseEndcapType(str, typeStream);
      }
    }
    stopTaskClock();
  } else if (myType=="Output") {
    getTill(configFile_, '{', false, true);
    typeConfig = getTill(configFile_, '}', false, true);
    if (typeConfig!="") {
      istringstream typeStream(typeConfig);
      parseOutput(typeStream);
    }
  } else if (myType=="PixelType") {
    getTill(configFile_, '}', false, true);
  } else {
    std::ostringstream tempSS; tempSS << "Unknown module type assignment keyword: " << myType;
    logERROR(tempSS);
    return false;
  }

  return true;
}


// Main config parser function. It opens the file, if possible
// then it skims the comments away, creating the readable stringstream
// and passes this to the type parser (parseObjectType).
// When this function is over a full tracker object should be there
// unless an exception was raised during the parsing (in which case
// the user will be informed about details through the stderr)
Tracker* configParser::parseFile(string configFileName, string typesFileName) {
  string str;
  Tracker* result = NULL;

  if (!typesFileName.empty()) peekTypes(typesFileName);

  if (rawConfigFile_.is_open()) {
    logERROR("Tracker config file is already open");
    return result;
  }

  rawConfigFile_.~ifstream();
  new ( (void*) &rawConfigFile_) std::ifstream();

  rawConfigFile_.open(configFileName.c_str());
  if (rawConfigFile_.is_open()) {

    // Reset the configFile_ object
    configFile_.~stringstream();
    new ( (void*) &configFile_) std::stringstream();

    // Skim comments delimited by // or # and put the skimmed file into configFile_
    string::size_type locComm1; // Location of commenting substring 1
    string::size_type locComm2; // Location of commenting substring 2

    while (getline(rawConfigFile_, str)) {
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
        if (!parseObjectType(str)) break;
      }
    } catch (exception& e) {
      ostringstream tempSS; tempSS << e.what() << endl;
      logERROR(tempSS);
      rawConfigFile_.close();
      if (myTracker_) delete myTracker_; myTracker_ = NULL;
      return NULL;
    }

    rawConfigFile_.close();
  } else {
    cerr << "ERROR: could not open tracker config file " << configFileName << endl;
    if (myTracker_) delete myTracker_; myTracker_ = NULL;
    return NULL;
  }
  //myTracker_->alignShortBarrels();
  myTracker_->sortLayers();

  // Eta cut and other post-operations
  std::pair<double, double> minMaxEta;
  minMaxEta = myTracker_->getEtaMinMax();
  std::ostringstream tempString;
  tempString.str("");
  tempString << "Eta coverage (min, max) of the tracker (prior to module purging): ("
      << minMaxEta.first << ", " << minMaxEta.second << ")";
  logINFO(tempString);
  myTracker_->cutOverEta(myTracker_->getEtaCut());
  minMaxEta = myTracker_->getEtaMinMax();
  tempString.str("");
  tempString << "Eta coverage (min, max) of the tracker (after module purging at eta "
    << myTracker_->getEtaCut() << "): ("
        << minMaxEta.first << ", " << minMaxEta.second << ")";

  logINFO(tempString);


  result = myTracker_;
  myTracker_ = NULL;
  return result;
}

Tracker* configParser::parsePixelsFromFile(string configFileName) {
  string str;
  Tracker* result = NULL;

  if (rawConfigFile_.is_open()) {
    cerr << "ERROR: Pixel config file is already open" << endl;
    return result;
  }

  rawConfigFile_.~ifstream();
  new ( (void*) &rawConfigFile_) std::ifstream();

  rawConfigFile_.open(configFileName.c_str());
  if (rawConfigFile_.is_open()) {

    // Reset the configFile_ object
    configFile_.~stringstream();
    new ( (void*) &configFile_) std::stringstream();

    // Skim comments delimited by // or # and put the skimmed file into configFile_
    string::size_type locComm1; // Location of commenting substring 1
    string::size_type locComm2; // Location of commenting substring 2

    while (getline(rawConfigFile_, str)) {
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
        if (str == "Pixels") {
          getTill(configFile_, '{', true);
          std::string pixelParams = getTill(configFile_, '}', false, true);
          if (pixelParams != "") {
            istringstream paramStream(pixelParams);
            parsePixels("pixelDetector", paramStream);
          }
        }
      }
    } catch (exception& e) {
      cerr << e.what() << endl;
      rawConfigFile_.close();
      if (myTracker_) delete myTracker_; myTracker_ = NULL;
      return NULL;
    }

    rawConfigFile_.close();
  } else {
    cerr << "ERROR: could not open tracker config file " << configFileName << endl;
    if (myTracker_) delete myTracker_; myTracker_ = NULL;
    return NULL;
  }
  // Stop here if there was no "Pixels" block to parse
  if (myTracker_ == NULL) return NULL;
  myTracker_->sortLayers();
  // Eta cut and other post-operations
  std::pair<double, double> minMaxEta;
  minMaxEta = myTracker_->getEtaMinMax();
  std::ostringstream tempString;
  tempString.str("");
  tempString << "Eta coverage (min, max) of the pixel detector (prior to module purging): ("
      << minMaxEta.first << ", " << minMaxEta.second << ")";
  logINFO(tempString);
  myTracker_->cutOverEta(myTracker_->getEtaCut());
  minMaxEta = myTracker_->getEtaMinMax();
  tempString.str("");
  tempString << "Eta coverage (min, max) of the pixel detector (after module purging at eta "
    << myTracker_->getEtaCut() << "): ("
        << minMaxEta.first << ", " << minMaxEta.second << ")";
  logINFO(tempString);
  result = myTracker_;
  myTracker_ = NULL;
  return result;
}

// "Dressing" config parser function. It opens the file, if possible
// then it skims the comments away, creating the readable stringstream
// and passes this to the dressing type parser (parseDressType).
// This function takes a tracker pointer, works on it according to the
// config file by "dressing" it. Returns a boolean value depending whether there
// were errors
// The user will be informed about details through the stderr
bool configParser::dressTracker(Tracker* aTracker, string configFileName) {
  string str;
  myTracker_=aTracker;

  if (rawConfigFile_.is_open()) {
    cerr << "ERROR: Module type config file is already open" << endl;
    myTracker_=NULL; return false;
  }

  rawConfigFile_.~ifstream();
  new ( (void*) &rawConfigFile_) std::ifstream();

  rawConfigFile_.open(configFileName.c_str());
  if (rawConfigFile_.is_open()) {

    // Reset the configFile_ object
    configFile_.~stringstream();
    new ( (void*) &configFile_) std::stringstream();

    // Skim comments delimited by // or # and put the skimmed file into configFile_
    string::size_type locComm1; // Location of commenting substring 1
    string::size_type locComm2; // Location of commenting substring 2

    while (getline(rawConfigFile_, str)) {
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
        if (!parseDressType(str)) break;
      }
    } catch (exception& e) {
      cerr << e.what() << endl;
      rawConfigFile_.close();
      myTracker_=NULL; return false;
    }

    rawConfigFile_.close();
  } else {
    cerr << "ERROR: could not open module type config file " << configFileName << endl;
    myTracker_=NULL; return false;
  }

  myTracker_=NULL; return true;
}

bool configParser::dressPixels(Tracker* aTracker, string configFileName) {
  string str, typeConfig;
  myTracker_=aTracker;

  if (rawConfigFile_.is_open()) {
    cerr << "ERROR: Module type config file is already open" << endl;
    myTracker_=NULL; return false;
  }

  rawConfigFile_.~ifstream();
  new ( (void*) &rawConfigFile_) std::ifstream();

  rawConfigFile_.open(configFileName.c_str());
  if (rawConfigFile_.is_open()) {

    // Reset the configFile_ object
    configFile_.~stringstream();
    new ( (void*) &configFile_) std::stringstream();

    // Skim comments delimited by // or # and put the skimmed file into configFile_
    string::size_type locComm1; // Location of commenting substring 1
    string::size_type locComm2; // Location of commenting substring 2

    while (getline(rawConfigFile_, str)) {
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
        if (str == "PixelType") {
          startTaskClock("Assigning module types to pixels");
          getTill(configFile_, '{', false, true);
          typeConfig = getTill(configFile_, '}', false, true);
          if (typeConfig != "") {
            istringstream typeStream(typeConfig);
            parsePixelType(typeStream);
          }
          stopTaskClock();  
        }
      }
    } catch (exception& e) {
      cerr << e.what() << endl;
      rawConfigFile_.close();
      myTracker_=NULL; return false;
    }

    rawConfigFile_.close();
  } else {
    cerr << "ERROR: could not open module type config file " << configFileName << endl;
    myTracker_=NULL; return false;
  }

  myTracker_=NULL; return true;
}

std::list<std::pair<int, double> >* configParser::parseSupportsFromFile(string fileName) {
  std::list<std::pair<int, double> >* result = NULL;
  std::ifstream infilestream;
  infilestream.open(fileName.c_str());
  if (infilestream.is_open()) {
    result = new std::list<std::pair<int, double> >;
    std::stringstream filecontents;
    std::string lineorword;
    std::string::size_type locComm1;
    std::string::size_type locComm2;
    while (std::getline(infilestream, lineorword)) {
      locComm1 = lineorword.find("//", 0);
      locComm2 = lineorword.find("#", 0);
      if ((locComm1 == string::npos) && (locComm2 == string::npos)) {
        filecontents << lineorword << endl;
      }
      else {
        filecontents << lineorword.substr(0, (locComm1 < locComm2 ? locComm1 : locComm2)) << endl;
      }
    }
    try {
      while (filecontents >> lineorword) {
        if (lineorword == "Support") {
          getTill(filecontents, '{', true);
          std::string configparams = getTill(filecontents, '}', false, true);
          if (configparams != "") {
            std::istringstream paramstream(configparams);
            if (!parseSupportParameters(paramstream, *result)) {
              if (result) delete result;
              return NULL;
            }
          }
        }
      }
      return result;
    }
    catch (std::exception& e) {
      cerr << e.what() << endl;
      if (result) delete result;
      return NULL;
    }
  }
  if (result) delete result;
  return NULL;
}


bool configParser::irradiateTracker(Tracker* tracker, string fileName) {
  // cout << "Trying to open irradiation map file " << fileName << endl;
  std::ifstream filein(fileName.c_str());
  if (!filein.is_open()) { 
    logERROR("Failed opening irradiation map file!");
    return false;
  }
  std::string line;
  while(std::getline(filein, line)) {
    if (line.find_first_of("#//;")==0 || line=="") continue;
    std::stringstream ss(line);
    double z, r = -1.0, fluence = -1.0, error = -1.0; // set to -1.0 to check for parsing errors
    ss >> z;
    ss >> r;
    ss >> fluence;
    ss >> error;
    if (r < 0.0 || fluence < 0.0) {
      ostringstream tempSS;
      tempSS << "Error while parsing irradiation map line: " << z << " ," << r << " ," << fluence << " ," << error;
      logERROR(tempSS);
    }
    tracker->getIrradiationMap()[make_pair(int(z/2.5),int(r/2.5))] = fluence;
  }
  return true;
}


bool configParser::peekTypes(string configFileName) {

  std::ifstream filein(configFileName.c_str());
  if (!filein.is_open()) {
    logWARNING("Couldn't opening types file for geometry peek. Skipping.");
    return false;
  }

  std::string line;
  std::string keyword, blockName; 
  for(int lineCount = 1; std::getline(filein, line); lineCount++) {
    line.erase(0, line.find_first_not_of(" \t"));
    line.erase(line.find_last_not_of(" \t")+1); // trim line
    if (line.find("//") != string::npos) line.erase(line.find("//")); // remove comments
    if (line.empty() || isspace(line.at(0))) continue; // skip empty lines

    std::istringstream parser(line);
    if (blockName.empty() && (line.find("BarrelType")==0 || line.find("EndcapType")==0)) { // enter block
      parser >> keyword;
      parser >> blockName;
      parser >> keyword;
      if (blockName == "{" || blockName.empty() || keyword != "{") { logERROR("Syntax error on line " + any2str(lineCount) + ": " + line); return false; }
      continue;
    }

    if (line.find("}") != string::npos) blockName.clear();

    if (!blockName.empty()) { // in a block
      std::string parameterName, parameterValue;
      if (!parseParameter(parameterName, parameterValue, parser)) {logERROR("Syntax error on line " + any2str(lineCount) + ": " + line); return false; }
      if (parameterName.find("dsDistance")!=string::npos && !parameterValue.empty()) {
        int firstIndex, secondIndex;
        if (!breakParameterName(parameterName, firstIndex, secondIndex)) { logERROR("Syntax error on line " + any2str(lineCount) + ": " + line); return false; }
        if (secondIndex == 0) geometryDsDistance_[blockName][firstIndex] = str2any<double>(parameterValue);
        else geometryDsDistanceSecond_[blockName][std::make_pair(firstIndex, secondIndex)] = str2any<double>(parameterValue);
      }
    }

  }
  filein.close();
  return true;
}


bool configParser::parseTilted(const std::string& fileName, const std::string& barrelName) {
  std::ifstream filein(fileName.c_str());
  if (!filein.is_open()) {
    logERROR("Cannot open file '"+ fileName +"'.");
    return false;
  }

  TiltedLayerSpecs tiltlay;
  tiltlay.numRods = 0;
  
  using std::string;
  string line;
  for (int lineCount = 1; std::getline(filein, line); lineCount++) {
    line = trim(line);
    vector<string> tokens = split(line, ",");
    if (tokens.size() == 7) { // parsing 2 parm sets per row (inner and outer rod) + num rods in phi
      vector<double> dtokens;
      std::transform(tokens.begin(), tokens.end(), std::back_inserter(dtokens), str2any<double>);
      TiltedModuleSpecs m1  = { dtokens[0], dtokens[1], dtokens[2] };
      TiltedModuleSpecs m2 = { dtokens[3], dtokens[4], dtokens[5] };
      if (m1.valid() && m2.valid()) {
        tiltlay.innerRod.push_back(m1);
        tiltlay.outerRod.push_back(m2);
        tiltlay.numRods = dtokens[6];
      }
    }
  }

  if (!tiltlay.valid()) { 
    logERROR("Failure while parsing spec file " + fileName + ". The resulting layer is invalid: numRods = " + any2str(tiltlay.numRods) + ", numMods/rod = " + any2str(tiltlay.innerRod.size()));
    return false;
  }

  tiltedBarrelSpecs_[barrelName].push_back(tiltlay);

  logDEBUG("Parsed tilted module layer " + any2str(tiltedBarrelSpecs_.size()) + " of barrel " + barrelName + ", composed of " + any2str(tiltlay.numRods) + " rods, with " + any2str(tiltlay.innerRod.size()) + " modules each.");


  return true;
}

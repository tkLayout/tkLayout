#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>
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

// Parses the parameter of the Tracker type in the config file.
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
    return true;
    
}


// Parses the part of config files between '{' and '}', creating a Barrel
bool configParser::parseBarrel(string myName, istream& inStream) {
    string parameterName;
    string parameterValue;
    
    int nBarrelLayers = 0;
    double minZ=0;
    double barrelRhoIn = 0;
    double barrelRhoOut = 0;
    int nBarrelModules = 0;
    string aDirective = "";
    BarrelModule* sampleBarrelModule = NULL;
    double aspectRatio = 1.;
    
    // Directives (this are communicated to the Tracker object)
    std::map<int, double> layerDirectives;
    std::map<int, LayerOption> layerOptions;
    
    
    // Tracker shoujld be already there
    if (!myTracker_) {
        cout << "Error: tracker is NULL when declaring barrel " << myName << endl;
        cout << "Missing declaration of a Tracker object?" << endl;
        throw parsingException();
    }
    
    
    while (!inStream.eof()) {
        while (parseParameter(parameterName, parameterValue, inStream)) {
            if (parameterName=="nLayers") {
                nBarrelLayers=atoi(parameterValue.c_str());
            } else if (parameterName=="minimumZ") {
                minZ=atof(parameterValue.c_str());
            } else if (parameterName=="aspectRatio") {
                aspectRatio=atof(parameterValue.c_str());
                if (aspectRatio<=0) {
		  cout << "Parsing barrel " << myName << endl
		       << "Wrong aspect ratio (height/width): " << parameterValue
		       << " should be a positive number" << endl;
		  throw parsingException();
		}
            } else if (parameterName=="nModules") {
                nBarrelModules=atoi(parameterValue.c_str());
            } else if (parameterName=="innerRadius") {
                barrelRhoIn=atof(parameterValue.c_str());
            } else if (parameterName=="outerRadius") {
                barrelRhoOut=atof(parameterValue.c_str());
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
                            std::cout << "Option: stacked for layer " << layerNum << " at " << aVal << std::endl; // debug
                        } else {
                            cout << "Wrong stack distance: (" << aVal << ") must be negative." << std::endl;
                            throw parsingException();
                        }
                    }
                    if (!gotIt) {
                        cout << "Parsing layer option for barrel " << myName
                                << ": unknown/nonsense option \"" << aString << "\"" << endl;
                        throw parsingException();
                    }
                } else {
                    cout << "Parsing barrel " << myName << endl
                            << "Wrong syntax for a layer option: \"" << parameterValue
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
        
        sampleBarrelModule = new BarrelModule(aspectRatio);   // Square modules of kind rphi
        
        // Important: if no directive was given, the following line will clear
        // possible previous directives coming from a different barrel
        myTracker_->setLayerDirectives(layerDirectives);
        myTracker_->setLayerOptions(layerOptions);
        
        myTracker_->buildBarrel(nBarrelLayers,
                barrelRhoIn,
                barrelRhoOut,
                nBarrelModules,
                sampleBarrelModule,
                myName,
                Layer::NoSection,
                true,
                minZ); // Actually build a compressed barrel (mezzanine or normal)
        
        delete sampleBarrelModule; // Dispose of the sample module
    } else {
        cout << "Missing mandatory parameter for barrel " << myName << endl;
        throw parsingException();
    }
    
    return true;
}

// Parses the part of config files between '{' and '}', creating an Endcap
bool configParser::parseEndcap(string myName, istream &inStream) {
    string parameterName;
    string parameterValue;
    
    int nDisks = 0;       // nDisks (per side)
    double barrelToEndcap = 0;  // gap between barrel and endcap
    double minZ = 0;            // explicit minimum z
    double maxZ = 0;            // maximum z
    double rhoIn = 0;
    double rhoOut = 0;
    double innerEta = 0;
    int diskParity = 0;
    int shapeType = Module::Wedge;
    double aspectRatio = 1.;
    
    map<pair<int, int>, bool> mapDiskRingRemoveToOuter;
    map<pair<int, int>, bool>::iterator mapDiskRingRemoveToOuterIt;
    
    // Directives (this are communicated to the Tracker object)
    std::map<int, int> ringDirective;
    
    // Tracker should be already there
    if (!myTracker_) {
        cout << "Error: tracker is NULL when declaring endcap " << myName << endl;
        cout << "Missing declaration of a Tracker object?" << endl;
        throw parsingException();
    }
    
    while (!inStream.eof()) {
        while (parseParameter(parameterName, parameterValue, inStream)) {
            if (parameterName=="nDisks") {
                nDisks=atoi(parameterValue.c_str());
            } else if (parameterName=="innerEta") {
                innerEta=atof(parameterValue.c_str());
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
            } else if (parameterName=="aspectRatio") {
                aspectRatio=atof(parameterValue.c_str());
                if (aspectRatio<=0) {
		  cout << "Parsing endcap " << myName << endl
		       << "Wrong aspect ratio (height/width): " << parameterValue
		       << " should be a positive number" << endl;
		  throw parsingException();
		}
	    } else if (parameterName=="shape") {
	      bool syntaxOk=true;
	      if (parameterValue=="wedge") {
		shapeType=Module::Wedge;
	      } else if (parameterValue=="rectangular") {
		shapeType=Module::Rectangular;
	      } else {
		syntaxOk=false;
	      }
	      if (!syntaxOk) {
		cout << "Parsing endcap " << myName << endl
		     << "Wrong syntax for a shape: \"" << parameterValue
		     << "\" should be \"rectangular\" or \"wedge\"" << endl;
		throw parsingException();
	      }
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
                    cout << "Parsing endcap " << myName << endl
                            << "Wrong syntax for a ring directive: \"" << parameterValue
                            << "\" should be like D1R4+ (to remove rings 4 or more from disk 1)" << endl;
                    throw parsingException();
                }
            } else {
                cout << "Unknown parameter \"" << parameterName << "\"" << endl;
                throw parsingException();
            }
            cout << "\t" << parameterName << " = " << parameterValue << ";" << endl; // debug
        }
    }
    
    // Chose the proper endcap minimum Z according to the parsed parameters
    if ((innerEta!=0)&&(rhoIn!=0)) {
        cout << "Parsing endcap " << myName << endl
                << " I got both innerEta = " << innerEta << " and innerRadius = " << rhoIn << endl
                << "This is inconsistent: please choose one of the two methods to place the endcap" << endl;
        throw parsingException();
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
    ((rhoIn != 0)||(innerEta!=0)) &&
    (rhoOut != 0) &&
    (minZ != 0) &&
    (maxZ != 0) &&
    (diskParity != 0)) {
      // The same old sample module
      EndcapModule* sampleModule;
      if (shapeType==Module::Wedge) {
	sampleModule = new EndcapModule(Module::Wedge);
      } else if (shapeType==Module::Rectangular) {
	sampleModule = new EndcapModule(aspectRatio);
      } else {
	std::cout << "ERROR: an unknown module shape type was generated inside the configuration parser" 
		  << std::endl << "this should never happen!" << std::endl;
	throw parsingException();
      }
        // Important: if no directive was given, the following line will clear
        // possible previous directives coming from a different endcap
        myTracker_->setRingDirectives(ringDirective);
        
        if (rhoIn!=0) {
            myTracker_->buildEndcaps(nDisks,     // nDisks (per side)
                    minZ,  // minZ
                    maxZ,
                    rhoIn,
                    rhoOut,
                    sampleModule,
                    myName,
                    diskParity,
                    Layer::NoSection);
        } else {
            myTracker_->buildEndcapsAtEta(nDisks,     // nDisks (per side)
                    minZ,  // minZ
                    maxZ,
                    innerEta,
                    rhoOut,
                    sampleModule,
                    myName,
                    diskParity,
                    Layer::NoSection);
        }
        
        delete sampleModule; // Dispose of the sample module
        
        // Remove the deleted rings
        
        for (mapDiskRingRemoveToOuterIt=mapDiskRingRemoveToOuter.begin();
        mapDiskRingRemoveToOuterIt!=mapDiskRingRemoveToOuter.end();
        mapDiskRingRemoveToOuterIt++) {
            int iDisk = ((*mapDiskRingRemoveToOuterIt).first).first;
            int iRing = ((*mapDiskRingRemoveToOuterIt).first).second;
            bool directionOuter = (*mapDiskRingRemoveToOuterIt).second;
            
            cout << "myTracker_->removeDiskRings("<<myName<<","<<iDisk<<", "<<iRing<<", "<<directionOuter<<");" << endl; // debug
            myTracker_->removeDiskRings(myName, iDisk, iRing, directionOuter);
        }
        
    } else {
        cout << "Missing mandatory parameter for endcap " << myName << endl;
        throw parsingException();
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

// The following two functions are temporary: by now we have no difference in syntax
// between the barrel and the endcap
bool configParser::parseBarrelType(string myName, istream& inStream) {
    return parseAnyType(myName, inStream);
}

bool configParser::parseEndcapType(string myName, istream& inStream) {
    return parseAnyType(myName, inStream);
}

// TODO: put the code into parseBarrelType and parseEndcapType when they will actually
// contain different parameters
bool configParser::parseAnyType(string myName, istream& inStream) {
    string parameterName;
    string parameterValue;
    int mainIndex; // Used for the ring and layer
    int secondaryIndex; // Used for the disk
    
    map<int, int> nStripsAcross;
    map<int, int> nSides;
    map<int, int> nSegments;
    map<int, string> type;
    
    pair<int, int> specialIndex; // used to indicate ring,disk
    
    map<pair<int, int>, int> nStripsAcrossSecond;
    map<pair<int, int>, int> nSidesSecond;
    map<pair<int, int>, int> nSegmentsSecond;
    map<pair<int, int>, string> typeSecond;
    map<pair<int, int>, bool> specialSecond;
    
    // Tracker should be already there
    if (!myTracker_) {
        cout << "Error: tracker is NULL when trying to assign types to a barrel " << myName << endl;
        return false;
    }
    
    // TODO: decide whether to use nStripAcross or nStripsAcross
    while (!inStream.eof()) {
        while (parseParameter(parameterName, parameterValue, inStream)) {
            if (!breakParameterName(parameterName, mainIndex, secondaryIndex)) {
                cerr << "Error: parameter name " << parameterName << " must have the following structure: "
                        << "parameter[ring] or parameter[layer] or parameter[ring,disk] e.g. nModules[2]" << endl;
                throw parsingException();
            } else {
                if (secondaryIndex==0) { // Standard assignment per ring or per layer
                    if (parameterName=="nStripsAcross") {
                        nStripsAcross[mainIndex]=atoi(parameterValue.c_str());
                    } else if (parameterName=="nStripAcross") {
                        nStripsAcross[mainIndex]=atoi(parameterValue.c_str());
                    } else if (parameterName=="nSides") {
                        nSides[mainIndex]=atoi(parameterValue.c_str());
                    } else if (parameterName=="nSegments") {
                        nSegments[mainIndex]=atoi(parameterValue.c_str());
                    } else if (parameterName=="type") {
                        type[mainIndex]=parameterValue.c_str();
                    }
                    cout << "\t" << parameterName << "[" << mainIndex << "] = " << parameterValue << ";" << endl; // debug
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
                    }
                    if (!isSpecial) {
                        cerr << "WARNING: the special parameter "
                                << parameterName << "[" << mainIndex << "," << secondaryIndex << "] is setting the same "
                                << "values as the default parameter "
                                << parameterName << "[" << mainIndex << "]. Ignoring it." << endl;
                    } else {
                        cout << "\t" << parameterName << "[" << mainIndex << ","<<secondaryIndex<<"] = " << parameterValue << ";" << endl; // debug
                    }
                }
            }
        }
    }
    
    
    myTracker_->setModuleTypes(myName,
            nStripsAcross, nSides, nSegments, type,
            nStripsAcrossSecond, nSidesSecond, nSegmentsSecond, typeSecond, specialSecond);
    
    return true;
    
}

// Output stuff
bool configParser::parseOutput(istream& inStream) {
    string parameterName;
    string parameterValue;
    
    string outPath="";
    
    // Tracker should be already there
    if (!myTracker_) {
        cout << "Error: tracker is NULL when trying to assign output options" << endl;
        return false;
    }
    
    // TODO: decide whether to use nStripAcross or nStripsAcross
    while (!inStream.eof()) {
        while (parseParameter(parameterName, parameterValue, inStream)) {
            if (parameterName=="Path") {
                outPath=parameterValue;
            } else {
                cerr << "While parting output parameters, I got an unrecognized parameter: " << parameterName << endl;
                throw parsingException();
            }
        }
    }
    
    if (outPath!="") myTracker_->setActiveDirectory(outPath);
    
    return true;
    
}

bool configParser::parseSupportParameters(istream& inStream, list<double>& plist) {
    string name, value;
    double mid_z;
    while (!inStream.eof()) {
        while (parseParameter(name, value, inStream)) {
            if (name == "midZ") {
                mid_z = atof(value.c_str());
                if (mid_z >= 0) plist.push_back(mid_z);
            }
            else {
                return false;
            }
        }
    }
    list<double>::iterator iter = plist.begin();
    while (iter != plist.end()) {
        if (*iter < 0) iter = plist.erase(iter);
        else iter++;
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
    } else if (myType=="Support") {
        getTill(configFile_, '}', false, true);
    } else {
        cerr << "Error: unknown piece of tracker " << myType;
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
        cout << "ASSIGNING MODULE TYPES TO BARREL:\t";
        str=getTill(configFile_, '{', true);
        if (str!="") {
            cout << str << endl;
            typeConfig=getTill(configFile_, '}', false);
            if (typeConfig!="") {
                istringstream typeStream(typeConfig);
                parseBarrelType(str, typeStream);
            }
        }
    } else if (myType=="EndcapType") {
        cout << "ASSIGNING MODULE TYPES TO ENDCAP:\t";
        str=getTill(configFile_, '{', true);
        if (str!="") {
            cout << str << endl;
            typeConfig=getTill(configFile_, '}', false);
            if (typeConfig!="") {
                istringstream typeStream(typeConfig);
                parseEndcapType(str, typeStream);
            }
        }
    } else if (myType=="Output") {
        getTill(configFile_, '{', false, true);
        typeConfig = getTill(configFile_, '}', false, true);
        if (typeConfig!="") {
            istringstream typeStream(typeConfig);
            parseOutput(typeStream);
        }
    } else {
        cerr << "Error: unknown module type assignment keyword: " << myType;
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
Tracker* configParser::parseFile(string configFileName) {
    string str;
    Tracker* result = NULL;
    
    if (rawConfigFile_.is_open()) {
        cerr << "Tracker config file is already open" << endl;
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
            cerr << e.what() << endl;
            rawConfigFile_.close();
            if (myTracker_) delete myTracker_; myTracker_ = NULL;
            return NULL;
        }
        
        rawConfigFile_.close();
    } else {
        cerr << "Error: could not open tracker config file " << configFileName << endl;
        if (myTracker_) delete myTracker_; myTracker_ = NULL;
        return NULL;
    }
    myTracker_->alignShortBarrels();
    myTracker_->sortLayers();
    
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
        cerr << "Module type config file is already open" << endl;
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
        cerr << "Error: could not open module type config file " << configFileName << endl;
        myTracker_=NULL; return false;
    }
    
    myTracker_=NULL; return true;
}

list<double>* configParser::parseSupportsFromFile(string fileName) {
    list<double>* result = NULL;
    ifstream infilestream;
    infilestream.open(fileName.c_str());
    if (infilestream.is_open()) {
        result = new list<double>;
        stringstream filecontents;
        string lineorword;
        string::size_type locComm1;
        string::size_type locComm2;
        while (getline(infilestream, lineorword)) {
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
                    string configparams = getTill(filecontents, '}', false, true);
                    if (configparams != " ") {
                        istringstream paramstream(configparams);
                        if (!parseSupportParameters(paramstream, *result)) {
                            if (result) delete result;
                            return NULL;
                        }
                    }
                }
            }
            return result;
        }
        catch (exception& e) {
            cerr << e.what() << endl;
            if (result) delete result;
            return NULL;
        }
    }
    if (result) delete result;
    return NULL;
}

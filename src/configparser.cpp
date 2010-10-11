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

configParser::configParser() : MessageLogger("configParser") {
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
                    tempString.str(""); tempString << "Ignoring extra parameter value " << dummy;
                    addMessage(tempString, WARNING);
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
        cerr << "ERROR: tracker is not NULL when declaring tracker " << myName << endl;
        cerr << "Double declaration of a Tracker object?" << endl;
        throw parsingException();
    }
    
    myTracker_ = new Tracker(myName);
    
    while (!inStream.eof()) {
        while (parseParameter(parameterName, parameterValue, inStream)) {
            if (parameterName=="nMB") { // Number of minimum bias events per BX
                doubleValue=atof(parameterValue.c_str());
                myTracker_->setNMB(doubleValue);
            } else if (parameterName=="zError") {
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
    
    int nBarrelLayers = 0;
    double minZ=0;
    double barrelRhoIn = 0;
    double barrelRhoOut = 0;
    int nBarrelModules = 0;
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
            } else if (parameterName=="minimumZ") {
                minZ=atof(parameterValue.c_str());
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
                            tempString.str(""); tempString << "Option: stacked for layer " << layerNum << " at " << aVal;
                            addMessage(tempString, DEBUG);
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
            } else {
                cerr << "ERROR: Unknown parameter \"" << parameterName << "\"" << endl;
                throw parsingException();
            }
            // cout << "\t" << parameterName << " = " << parameterValue << ";" << endl; // debug
        }
    }
    
    
    // Actually creating the barrel if all the mandatory parameters were set
    if ( (nBarrelLayers != 0) &&
    (barrelRhoIn != 0) &&
    (barrelRhoOut != 0) &&
    (nBarrelModules != 0) ) {
        
        
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
        
        // Important: if no directive was given, the following line will clear
        // possible previous directives coming from a different barrel
        myTracker_->setLayerDirectives(layerDirectives);
        myTracker_->setLayerOptions(layerOptions);
        myTracker_->setPhiSegments(phiSegments);
        
        sampleBarrelModule->setResolutionRphi();
        sampleBarrelModule->setResolutionY();
        
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
    
    int nDisks = 0;       // nDisks (per side)
    double barrelToEndcap = 0;  // gap between barrel and endcap
    double minZ = 0;            // explicit minimum z
    double maxZ = 0;            // maximum z
    double rhoIn = 0;
    double rhoOut = 0;
    double innerEta = 0;
    int diskParity = 0;
    int shapeType = Module::Wedge;
    bool explicitShapeType = false;
    double aspectRatio = 1.;
    bool explicitAspectRatio = false;
    int phiSegments = 4;
    
    // Fetch the generic Delta of the tracker
    double genericSmallDelta = myTracker_->getSmallDelta();
    double genericBigDelta = myTracker_->getBigDelta();
    
    map<pair<int, int>, bool> mapDiskRingRemoveToOuter;
    map<pair<int, int>, bool>::iterator mapDiskRingRemoveToOuterIt;
    
    // Directives (this are communicated to the Tracker object)
    std::map<int, int> ringDirective;
    
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
            if (parameterName=="nDisks") {
                nDisks=atoi(parameterValue.c_str());
            } else if (parameterName=="smallDelta") {
                myTracker_->setSmallDelta(atof(parameterValue.c_str()));
            } else if (parameterName=="bigDelta") {
                myTracker_->setBigDelta(atof(parameterValue.c_str()));
            } else if ((correctlyBroken)&&((parameterNameCopy=="smallDelta")||(parameterNameCopy=="bigDelta"))) {
                double deltaValue =  atof(parameterValue.c_str());
                if ((mainIndex==0)||(secondaryIndex!=0)||(deltaValue==0)) {
                    cerr << "ERROR: parameter setting for " << parameterNameCopy << " must have the following structure: "
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
                    // TODO: fix this
                    cerr << "ERROR: per-ring small and big deltas are not yet imlpemented" << endl;
                    throw parsingException();
                }
            } else if (parameterName=="phiSegments") {
                phiSegments=atoi(parameterValue.c_str());
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
                explicitAspectRatio = true;
                if (aspectRatio<=0) {
                    cerr << "ERROR: Parsing endcap \"" << myName
                            << "\": wrong aspect ratio (height/width): " << parameterValue
                            << " should be a positive number" << endl;
                    throw parsingException();
                }
            } else if (parameterName=="shape") {
                bool syntaxOk=true;
                if (parameterValue=="wedge") {
                    shapeType=Module::Wedge;
                    explicitShapeType = true;
                } else if (parameterValue=="rectangular") {
                    shapeType=Module::Rectangular;
                    explicitShapeType = true;
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
    
    // Check consistency in module shape type assigned
    if (explicitShapeType) {
        if ((shapeType==Module::Wedge)&&(explicitAspectRatio)) {
            cerr << "Parsing endcap \"" << myName
                    << "\": I see module shape 'wedge' ans aspect ratio assigned. This is inconsistent. I quit." << endl;
            throw parsingException();
        }
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
            cerr << "ERROR: an unknown module shape type was generated inside the configuration parser."
                    << "this should never happen!" << std::endl;
            throw parsingException();
        }
        // Important: if no directive was given, the following line will clear
        // possible previous directives coming from a different endcap
        myTracker_->setRingDirectives(ringDirective);
        myTracker_->setPhiSegments(phiSegments);
        
        sampleModule->setResolutionRphi();
        sampleModule->setResolutionY();
        
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
            
            // cerr << "myTracker_->removeDiskRings("<<myName<<","<<iDisk<<", "<<iRing<<", "<<directionOuter<<");" << endl; // debug
            myTracker_->removeDiskRings(myName, iDisk, iRing, directionOuter);
        }
        
    } else {
        cerr << "ERROR: Missing mandatory parameter for endcap " << myName << endl ;
        cerr << "Mandatory parameters are: nDisks" << endl
                << "                          [ innerRadius | innerEta ]" << endl
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
    double rIn = 0.0;
    double rOut = 0.0;
    double barrelToEndcap = 0.0;
    double maxZ = 0.0;
    int nLayers = 0;
    int nModules = 0;
    int nDisks = 0;
    int phiSegments = 4;
    int diskParity = 0;
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
            } else if (parameterName == "innerRadius") {
                rIn = atof(parameterValue.c_str());
            } else if (parameterName == "outerRadius") {
                rOut = atof(parameterValue.c_str());
            } else if (parameterName == "phiSegments") {
                phiSegments = atoi(parameterValue.c_str());
            } else if (parameterName == "barrelGap") {
                barrelToEndcap = atof(parameterValue.c_str());
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
            } else if ( "option") {
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
                            tempString.str(""); tempString << "Option: stacked for layer " << layerNum << " at " << aVal;
                            addMessage(tempString, DEBUG);
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
            } else if ( "directive") {
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
    myTracker_->buildBarrel(nLayers, rIn, rOut, nModules, sampleBarrelModule, myName, Layer::NoSection, true, 0.0);
    delete sampleBarrelModule;
    
    if (nDisks > 0) {
        waferDiameter = pow((pow(emsize.first, 2) + pow(emsize.second, 2)), 0.5);
        aspectRatio = emsize.second / emsize.first; // height / width
        EndcapModule* sampleEndcapModule = new EndcapModule(waferDiameter, aspectRatio);
        sampleEndcapModule->setResolutionRphi();
        sampleEndcapModule->setResolutionY();
        myTracker_->buildEndcaps(nDisks, myTracker_->getMaxBarrelZ(+1) + barrelToEndcap,
                maxZ, rIn, rOut, sampleEndcapModule, myName, diskParity, Layer::NoSection);
        delete sampleEndcapModule;
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
    map<int, int> nSides;
    map<int, int> nSegments;
    map<int, string> type;
    map<int, double> dsDistance;
    map<int, double> dsRotation;
    
    pair<int, int> specialIndex; // used to indicate ring,disk
    
    map<pair<int, int>, int> nStripsAcrossSecond;
    map<pair<int, int>, int> nSidesSecond;
    map<pair<int, int>, int> nSegmentsSecond;
    map<pair<int, int>, string> typeSecond;
    map<pair<int, int>, double> dsDistanceSecond;
    map<pair<int, int>, double> dsRotationSecond;
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
                    } else if (parameterName=="nStripAcross") {
                        nStripsAcross[mainIndex]=atoi(parameterValue.c_str());
                    } else if (parameterName=="nSides") {
                        nSides[mainIndex]=atoi(parameterValue.c_str());
                    } else if (parameterName=="nSegments") {
                        nSegments[mainIndex]=atoi(parameterValue.c_str());
                    } else if (parameterName=="type") {
                        type[mainIndex]=parameterValue.c_str();
                    } else if (parameterName == "dsDistance") {
                        dsDistance[mainIndex]=atof(parameterValue.c_str());
                    } else if (parameterName == "dsRotation") {
                        dsRotation[mainIndex]=atof(parameterValue.c_str());
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
                        if (atoi(parameterValue.c_str())!=dsDistance[mainIndex]) {
                            dsDistanceSecond[specialIndex]=atof(parameterValue.c_str());
                            specialSecond[specialIndex]=true;
                            isSpecial = true;
                        }
                    } else if (parameterName=="dsRotation") {
                        if (atoi(parameterValue.c_str())!=dsRotation[mainIndex]) {
                            dsRotationSecond[specialIndex]=atof(parameterValue.c_str());
                            specialSecond[specialIndex]=true;
                            isSpecial = true;
                        }
                    }
                    if (!isSpecial) {
                        tempString.str(""); tempString << "The special parameter "
                                << parameterName << "[" << mainIndex << "," << secondaryIndex << "] is setting the same "
                                << "values as the default parameter "
                                << parameterName << "[" << mainIndex << "]. Ignoring it.";
                        addMessage(tempString, WARNING);
                    } // else {
                    // cout << "\t" << parameterName << "[" << mainIndex << ","<<secondaryIndex<<"] = " << parameterValue << ";" << endl; // debug
                    // }
                }
            }
        }
    }
    
    myTracker_->setModuleTypes(myName,
            nStripsAcross, nSides, nSegments, type, dsDistance, dsRotation,
            nStripsAcrossSecond, nSidesSecond, nSegmentsSecond, typeSecond,
            dsDistanceSecond, dsRotationSecond, specialSecond);
    
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
        cout << "Creating barrel ";
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
        cout << "Creating endcap ";
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
    } else if (myType=="Pixels") {
        getTill(configFile_, '}', false, true);
    } else {
        cerr << "ERROR: unknown piece of tracker " << myType;
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
        cout << "Assigning module types to barrel ";
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
        cout << "Assigning module types to endcap ";
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
    } else if (myType=="PixelType") {
        getTill(configFile_, '}', false, true);
    } else {
        cerr << "ERROR: unknown module type assignment keyword: " << myType;
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
        cerr << "ERROR: Tracker config file is already open" << endl;
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
        cerr << "ERROR: could not open tracker config file " << configFileName << endl;
        if (myTracker_) delete myTracker_; myTracker_ = NULL;
        return NULL;
    }
    myTracker_->alignShortBarrels();
    myTracker_->sortLayers();
    
    // Eta cut and other post-operations
    std::pair<double, double> minMaxEta;
    minMaxEta = myTracker_->getEtaMinMax();
    tempString.str("");
    tempString << "Eta coverage (min, max) of the tracker (prior to module purging): ("
            << minMaxEta.first << ", " << minMaxEta.second << ")";
    addMessage(tempString, INFO);
    myTracker_->cutOverEta(myTracker_->getEtaCut());
    minMaxEta = myTracker_->getEtaMinMax();
    tempString.str("");
    tempString << "Eta coverage (min, max) of the tracker (after module purging at eta "
    << myTracker_->getEtaCut() << "): ("
    << minMaxEta.first << ", " << minMaxEta.second << ")";
    addMessage(tempString, INFO);
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
    tempString.str("");
    tempString << "Eta coverage (min, max) of the pixel detector (prior to module purging): ("
            << minMaxEta.first << ", " << minMaxEta.second << ")";
    addMessage(tempString, INFO);
    myTracker_->cutOverEta(myTracker_->getEtaCut());
    minMaxEta = myTracker_->getEtaMinMax();
    tempString.str("");
    tempString << "Eta coverage (min, max) of the pixel detector (after module purging at eta "
    << myTracker_->getEtaCut() << "): ("
    << minMaxEta.first << ", " << minMaxEta.second << ")";
    addMessage(tempString, INFO);
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
                    cout << "Assigning module types to pixels " << endl;
                    getTill(configFile_, '{', false, true);
                    typeConfig = getTill(configFile_, '}', false, true);
                    if (typeConfig != "") {
                        istringstream typeStream(typeConfig);
                        parsePixelType(typeStream);
                    }
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

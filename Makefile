# Define the REVISIONNUMBER variable to have it
# apearing in the root web page output
# -DREVISIONNUMBER=555
REVISION=$(which svnversion > /dev/null && svnversion)
DEFINES=`./getVersionDefine`

ROOTFLAGS=`root-config --cflags`
ROOTLIBDIR=`root-config --libdir`
ROOTLIBFLAGS=`root-config --libs`
BOOSTLIBFLAGS=-lboost_filesystem -lboost_regex
GEOMLIBFLAG=-lGeom
GLIBFLAGS=`root-config --glibs`
INCLUDEFLAGS=-Iinclude/
SRCDIR=src
INCDIR=include
LIBDIR=lib

COMP=g++ -Wall $(INCLUDEFLAGS) $(DEFINES)

all: tkgeometry exocom general elements ushers dressers viz naly squid tkLayout testObj tkmaterial rootwebTest 
	@echo "Full build successful."

#TKGEOMETRY
tkgeometry: $(LIBDIR)/configparser.o $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o
	@echo "Built target 'tkgeometry'."

$(LIBDIR)/configparser.o:	$(SRCDIR)/configparser.cpp $(INCDIR)/configparser.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/configparser.o $(SRCDIR)/configparser.cpp

$(LIBDIR)/module.o: $(SRCDIR)/module.cpp $(INCDIR)/module.hh
	@echo "Building target module.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/module.o $(SRCDIR)/module.cpp 
	@echo "Built target module.o"

$(LIBDIR)/layer.o: $(SRCDIR)/layer.cpp $(INCDIR)/layer.hh
	@echo "Building target layer.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/layer.o $(SRCDIR)/layer.cpp 
	@echo "Built target layer.o"

$(LIBDIR)/tracker.o: $(SRCDIR)/tracker.cpp $(INCDIR)/tracker.hh
	@echo "Building target tracker.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/tracker.o $(SRCDIR)/tracker.cpp 	
	@echo "Built target tracker.o"

#EXOCOM
exocom:  $(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o $(LIBDIR)/XMLWriter.o
	@echo "Built target 'exocom'."

$(LIBDIR)/MatParser.o: $(SRCDIR)/MatParser.cc $(INCDIR)/MatParser.h
	@echo "Building target MatParser.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/MatParser.o $(SRCDIR)/MatParser.cc
	@echo "Built target MatParser.o"
	
$(LIBDIR)/Extractor.o: $(SRCDIR)/Extractor.cc $(INCDIR)/Extractor.h
	@echo "Building target Extractor.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Extractor.o $(SRCDIR)/Extractor.cc
	@echo "Built target Extractor.o"
	
$(LIBDIR)/XMLWriter.o: $(SRCDIR)/XMLWriter.cc $(INCDIR)/XMLWriter.h
	@echo "Building target XMLWriter.o..."
	$(COMP) -c -o $(LIBDIR)/XMLWriter.o $(SRCDIR)/XMLWriter.cc
	@echo "Built target XMLWriter.o"

#GENERAL
general: $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialProperties.o \
	$(LIBDIR)/InactiveSurfaces.o
	@echo "Built target 'general'."

$(LIBDIR)/MaterialBudget.o: $(SRCDIR)/MaterialBudget.cc $(INCDIR)/MaterialBudget.h
	@echo "Building target MaterialBudget.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/MaterialBudget.o $(SRCDIR)/MaterialBudget.cc
	@echo "Built target MaterialBudget.o"

$(LIBDIR)/MaterialTable.o: $(SRCDIR)/MaterialTable.cc $(INCDIR)/MaterialTable.h
	@echo "Building target MaterialTable.o..."
	$(COMP) -c -o $(LIBDIR)/MaterialTable.o $(SRCDIR)/MaterialTable.cc
	@echo "Built target MaterialTable.o"

$(LIBDIR)/MaterialProperties.o: $(SRCDIR)/MaterialProperties.cc $(INCDIR)/MaterialProperties.h
	@echo "Building target MaterialProperties.o..."
	$(COMP) -c -o $(LIBDIR)/MaterialProperties.o $(SRCDIR)/MaterialProperties.cc
	@echo "Built target MaterialProperties.o"

$(LIBDIR)/InactiveSurfaces.o: $(SRCDIR)/InactiveSurfaces.cc $(INCDIR)/InactiveSurfaces.h
	@echo "Building target InactiveSurfaces.o..."
	$(COMP) -c -o $(LIBDIR)/InactiveSurfaces.o $(SRCDIR)/InactiveSurfaces.cc
	@echo "Built target InactiveSurfaces.o"
	
#ELEMENTS
elements: $(LIBDIR)/ModuleCap.o $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o $(LIBDIR)/InactiveTube.o
	@echo "Built target 'elements'."

 $(LIBDIR)/ModuleCap.o: $(SRCDIR)/ModuleCap.cc $(INCDIR)/ModuleCap.h
	@echo "Building target ModuleCap.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/ModuleCap.o $(SRCDIR)/ModuleCap.cc
	@echo "Built target ModuleCap.o"

$(LIBDIR)/InactiveElement.o: $(SRCDIR)/InactiveElement.cc $(INCDIR)/InactiveElement.h
	@echo "Building target InactiveElement.o..."
	$(COMP) -c -o $(LIBDIR)/InactiveElement.o $(SRCDIR)/InactiveElement.cc
	@echo "Built target InactiveElement.o"

$(LIBDIR)/InactiveRing.o: $(SRCDIR)/InactiveRing.cc $(INCDIR)/InactiveRing.h
	@echo "Building target InactiveRing.o..."
	$(COMP) -c -o $(LIBDIR)/InactiveRing.o $(SRCDIR)/InactiveRing.cc
	@echo "Built target InactiveRing.o."

$(LIBDIR)/InactiveTube.o: $(SRCDIR)/InactiveTube.cc $(INCDIR)/InactiveTube.h
	@echo "Building target InactiveTube.o..."
	$(COMP) -c -o $(LIBDIR)/InactiveTube.o $(SRCDIR)/InactiveTube.cc
	@echo "Built target InactiveTube.o."

#USHERS
ushers: $(LIBDIR)/Usher.o
	@echo "Built target 'ushers'."

$(LIBDIR)/Usher.o: $(SRCDIR)/Usher.cc $(INCDIR)/Usher.h
	@echo "Building target Usher.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Usher.o $(SRCDIR)/Usher.cc
	@echo "Built target Usher.o"

#DRESSERS
dressers: $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o
	@echo "Built target 'dressers'."

$(LIBDIR)/MatCalc.o: $(SRCDIR)/MatCalc.cc $(INCDIR)/MatCalc.h
	@echo "Building target MatCalc.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/MatCalc.o $(SRCDIR)/MatCalc.cc
	@echo "Built target MatlCalc.o"

$(LIBDIR)/MatCalcDummy.o: $(SRCDIR)/MatCalcDummy.cc $(INCDIR)/MatCalcDummy.h
	@echo "Building target MatCalcDummy.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/MatCalcDummy.o $(SRCDIR)/MatCalcDummy.cc
	@echo "Built target MatCalcDummy.o"

#VISUALISATION
viz: $(LIBDIR)/Vizard.o $(LIBDIR)/tk2CMSSW.o
	@echo "Built target 'viz'."
	
$(LIBDIR)/Vizard.o: $(SRCDIR)/Vizard.cc $(INCDIR)/Vizard.h
	@echo "Building target Vizard.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Vizard.o $(SRCDIR)/Vizard.cc
	@echo "Built target Vizard.o"

$(LIBDIR)/tk2CMSSW.o: $(SRCDIR)/tk2CMSSW.cc $(INCDIR)/tk2CMSSW.h
	@echo "Building target tk2CMSSW.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/tk2CMSSW.o $(SRCDIR)/tk2CMSSW.cc
	@echo "Built target tk2CMSSW.o"

#ANALYSYS
naly: $(LIBDIR)/Analyzer.o
	@echo "Built target 'naly'."

$(LIBDIR)/Analyzer.o: $(SRCDIR)/Analyzer.cc $(INCDIR)/Analyzer.h
	@echo "Building target Analyzer.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Analyzer.o $(SRCDIR)/Analyzer.cc
	@echo "Built target Analyzer.o"

#SQUID
squid: $(LIBDIR)/Squid.o
	@echo "Built target 'squid'."
	
$(LIBDIR)/Squid.o: $(SRCDIR)/Squid.cc $(INCDIR)/Squid.h
	@echo "Building target Squid.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Squid.o $(SRCDIR)/Squid.cc
	@echo "Built target Squid.o"

#ROOTWEB
$(LIBDIR)/rootweb.o:	src/rootweb.cpp include/rootweb.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/rootweb.o src/rootweb.cpp

# Main (program) configuration file handler
$(LIBDIR)/mainConfigHandler.o: src/mainConfigHandler.cpp include/mainConfigHandler.h
	$(COMP) -c -o $(LIBDIR)/mainConfigHandler.o src/mainConfigHandler.cpp

#FINAL
tkLayout: TrackerGeom2.cpp $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o $(LIBDIR)/configparser.o
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o \
	$(LIBDIR)/configparser.o TrackerGeom2.cpp $(ROOTLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o tkLayout

tkmaterial: $(LIBDIR)/tkmaterial.o $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o $(LIBDIR)/configparser.o \
        $(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o $(LIBDIR)/XMLWriter.o $(LIBDIR)/MaterialTable.o \
        $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialProperties.o $(LIBDIR)/ModuleCap.o $(LIBDIR)/InactiveSurfaces.o  \
        $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o $(LIBDIR)/InactiveTube.o $(LIBDIR)/Usher.o \
        $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o $(LIBDIR)/Vizard.o $(LIBDIR)/tk2CMSSW.o \
        $(LIBDIR)/Analyzer.o $(LIBDIR)/Squid.o
	$(COMP) $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o $(LIBDIR)/configparser.o \
	$(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o $(LIBDIR)/XMLWriter.o $(LIBDIR)/MaterialTable.o \
	$(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialProperties.o $(LIBDIR)/ModuleCap.o  $(LIBDIR)/InactiveSurfaces.o  \
	$(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o $(LIBDIR)/InactiveTube.o $(LIBDIR)/Usher.o \
	$(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o $(LIBDIR)/Vizard.o $(LIBDIR)/tk2CMSSW.o \
	$(LIBDIR)/Analyzer.o $(LIBDIR)/Squid.o $(LIBDIR)/tkmaterial.o \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o tkmaterial

$(LIBDIR)/tkmaterial.o: tkmaterial.cpp
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/tkmaterial.o tkmaterial.cpp

testObj: testObjects.cpp $(LIBDIR)/module.o $(LIBDIR)/layer.o
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/module.o $(LIBDIR)/layer.o testObjects.cpp \
	$(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o testObj

rootwebTest: rootwebTest.cpp $(LIBDIR)/mainConfigHandler.o $(LIBDIR)/rootweb.o 
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/mainConfigHandler.o $(LIBDIR)/rootweb.o rootwebTest.cpp $(ROOTLIBFLAGS) $(BOOSTLIBFLAGS) -o rootwebTest

#CLEANUP
cleantkgeometry:
	-rm -f $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o $(LIBDIR)/configparser.o

cleanexocom:
	-rm -f $(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o $(LIBDIR)/XMLWriter.o

cleangeneral:
	-rm -f $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialProperties.o \
	$(LIBDIR)/InactiveSurfaces.o

cleanelements:
	-rm -f $(LIBDIR)/ModuleCap.o $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o $(LIBDIR)/InactiveTube.o

cleanushers:
	-rm -f $(LIBDIR)/Usher.o

cleandressers:
	-rm -f $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o 
	
cleanviz:
	-rm -r $(LIBDIR)/Vizard.o $(LIBDIR)/tk2CMSSW.o
	
cleannaly:
	-rm -r $(LIBDIR)/Analyzer.o

cleanrootweb:
	rm -f $(LIBDIR)/rootweb.o 
	
cleantkmaine:
	-rm -f $(LIBDIR)/Squid.o $(LIBDIR)/tkmaterial.o ./tkmaterial ./tkLayout ./testObj ./rootwebTest
	
clean: cleantkgeometry cleangeneral cleanelements cleanushers cleandressers cleanviz cleannaly cleanrootweb cleantkmaine





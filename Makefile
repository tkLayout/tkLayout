# Define the REVISIONNUMBER variable to have it
# apearing in the root web page output
# -DREVISIONNUMBER=555
REVISION=$(which svnversion > /dev/null && svnversion)
DEFINES=`./getVersionDefine`
DEFINES+=-DDEBUG_PERFORMANCE

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
BINDIR=bin
TESTDIR=test
DOCDIR=doc
DOXYDIR=doc/doxygen

COMP=g++ -Wall $(INCLUDEFLAGS) $(DEFINES)

bin: tkmaterial tklayout setup
	@echo "Executable built."

all: hit tkgeometry exocom general elements ushers dressers viz naly squid testObjects tkmaterial tklayout rootwebTest 
	@echo "Full build successful."

install:
	./install.sh

#TRACKS
hit: $(LIBDIR)/hit.o
	@echo "Built target 'hit'."

$(LIBDIR)/hit.o: $(SRCDIR)/hit.cpp $(INCDIR)/hit.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/hit.o $(SRCDIR)/hit.cpp

#TKGEOMETRY
tkgeometry: $(LIBDIR)/configparser.o $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o  \
	$(LIBDIR)/messageLogger.o $(LIBDIR)/mainConfigHandler.o
	@echo "Built target 'tkgeometry'."

$(LIBDIR)/configparser.o: $(SRCDIR)/configparser.cpp $(INCDIR)/configparser.hh
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

$(LIBDIR)/mainConfigHandler.o: src/mainConfigHandler.cpp $(INCDIR)/mainConfigHandler.h
	$(COMP) -c -o $(LIBDIR)/mainConfigHandler.o src/mainConfigHandler.cpp

$(LIBDIR)/messageLogger.o: src/messageLogger.cpp $(INCDIR)/messageLogger.h
	$(COMP) -c -o $(LIBDIR)/messageLogger.o src/messageLogger.cpp

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

#ROOT-related stuff
$(LIBDIR)/rootweb.o:	src/rootweb.cpp $(INCDIR)/rootweb.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/rootweb.o src/rootweb.cpp

#$(LIBDIR)/rootutils.o: src/rootutils.cpp $(INCDIR)/rootutils.h
#	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/rootutils.o src/rootutils.cpp

# UTILITIES
setup: $(BINDIR)/setup.bin
	@echo "Building setup..."

$(BINDIR)/setup.bin: setup.cpp $(LIBDIR)/mainConfigHandler.o
	$(COMP) $(BOOSTLIBFLAGS) $(LIBDIR)/mainConfigHandler.o setup.cpp -o $(BINDIR)/setup.bin

#FINAL
#(obsolete)
tkLayout: $(BINDIR)/tkLayout
	@echo "Building tkLayout..."

#(obsolete)
$(BINDIR)/tkLayout: TrackerGeom2.cpp tkgeometry
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o  $(LIBDIR)/messageLogger.o \
	$(LIBDIR)/configparser.o TrackerGeom2.cpp $(ROOTLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $(BINDIR)/tkLayout

tkmaterial: $(BINDIR)/tkmaterial
	@echo "Building tkmaterial..."

$(BINDIR)/tkmaterial: $(LIBDIR)/tkmaterial.o $(LIBDIR)/hit.o $(LIBDIR)/module.o $(LIBDIR)/layer.o \
	$(LIBDIR)/tracker.o  $(LIBDIR)/messageLogger.o $(LIBDIR)/configparser.o $(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o \
	$(LIBDIR)/XMLWriter.o $(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialProperties.o \
	$(LIBDIR)/ModuleCap.o  $(LIBDIR)/InactiveSurfaces.o  $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o \
	$(LIBDIR)/InactiveTube.o $(LIBDIR)/Usher.o $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o \
	$(LIBDIR)/Vizard.o $(LIBDIR)/tk2CMSSW.o $(LIBDIR)/Analyzer.o $(LIBDIR)/Squid.o $(LIBDIR)/mainConfigHandler.o \
	$(LIBDIR)/rootweb.o
	$(COMP) $(LIBDIR)/hit.o $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o  $(LIBDIR)/messageLogger.o \
	$(LIBDIR)/configparser.o $(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o $(LIBDIR)/XMLWriter.o \
	$(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialProperties.o $(LIBDIR)/ModuleCap.o \
	$(LIBDIR)/InactiveSurfaces.o $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o $(LIBDIR)/InactiveTube.o \
	$(LIBDIR)/Usher.o $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o $(LIBDIR)/Vizard.o \
	$(LIBDIR)/tk2CMSSW.o $(LIBDIR)/Analyzer.o $(LIBDIR)/Squid.o $(LIBDIR)/mainConfigHandler.o \
	$(LIBDIR)/rootweb.o $(LIBDIR)/tkmaterial.o \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $(BINDIR)/tkmaterial

$(LIBDIR)/tkmaterial.o: tkmaterial.cpp
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/tkmaterial.o tkmaterial.cpp

tklayout: $(BINDIR)/tklayout
	@echo "Building tklayout..."

$(BINDIR)/tklayout: $(LIBDIR)/tklayout.o $(LIBDIR)/hit.o $(LIBDIR)/module.o $(LIBDIR)/layer.o \
	$(LIBDIR)/tracker.o $(LIBDIR)/configparser.o $(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o \
	$(LIBDIR)/XMLWriter.o $(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialProperties.o \
	$(LIBDIR)/ModuleCap.o  $(LIBDIR)/InactiveSurfaces.o  $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o \
	$(LIBDIR)/InactiveTube.o $(LIBDIR)/Usher.o $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o \
	$(LIBDIR)/Vizard.o $(LIBDIR)/tk2CMSSW.o $(LIBDIR)/Analyzer.o $(LIBDIR)/Squid.o $(LIBDIR)/rootweb.o $(LIBDIR)/mainConfigHandler.o \
	$(LIBDIR)/messageLogger.o
	$(COMP) $(LIBDIR)/hit.o $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o \
	$(LIBDIR)/configparser.o $(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o $(LIBDIR)/XMLWriter.o \
	$(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialProperties.o $(LIBDIR)/ModuleCap.o \
	$(LIBDIR)/InactiveSurfaces.o $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o $(LIBDIR)/InactiveTube.o \
	$(LIBDIR)/Usher.o $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o $(LIBDIR)/Vizard.o \
	$(LIBDIR)/tk2CMSSW.o $(LIBDIR)/Analyzer.o $(LIBDIR)/Squid.o $(LIBDIR)/rootweb.o $(LIBDIR)/mainConfigHandler.o \
	$(LIBDIR)/messageLogger.o \
	$(LIBDIR)/tklayout.o \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $(BINDIR)/tklayout

$(LIBDIR)/tklayout.o: tklayout.cpp
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/tklayout.o tklayout.cpp

testObjects: $(TESTDIR)/testObjects
$(TESTDIR)/testObjects: $(TESTDIR)/testObjects.cpp $(LIBDIR)/module.o $(LIBDIR)/layer.o
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/module.o $(LIBDIR)/layer.o $(TESTDIR)/testObjects.cpp \
	$(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o $(TESTDIR)/testObjects

rootwebTest: $(TESTDIR)/rootwebTest
$(TESTDIR)/rootwebTest: $(TESTDIR)/rootwebTest.cpp $(LIBDIR)/mainConfigHandler.o $(LIBDIR)/rootweb.o 
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/mainConfigHandler.o $(LIBDIR)/rootweb.o $(TESTDIR)/rootwebTest.cpp $(ROOTLIBFLAGS) $(BOOSTLIBFLAGS) -o $(TESTDIR)/rootwebTest

#CLEANUP
cleanhit:
	@rm -f $(LIBDIR)/hit.o

cleantkgeometry:
	@rm -f $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o  $(LIBDIR)/messageLogger.o $(LIBDIR)/configparser.o $(LIBDIR)/mainConfigHandler.o

cleanexocom:
	@rm -f $(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o $(LIBDIR)/XMLWriter.o

cleangeneral:
	@rm -f $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialProperties.o \
	$(LIBDIR)/InactiveSurfaces.o

cleanelements:
	@rm -f $(LIBDIR)/ModuleCap.o $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o $(LIBDIR)/InactiveTube.o

cleanushers:
	@rm -f $(LIBDIR)/Usher.o

cleandressers:
	@rm -f $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o 
	
cleanviz:
	@rm -f $(LIBDIR)/Vizard.o $(LIBDIR)/tk2CMSSW.o
	
cleannaly:
	@rm -f $(LIBDIR)/Analyzer.o

cleanrootweb:
	@rm -f $(LIBDIR)/rootweb.o 

#cleanutils:
#	@rm -f $(LIBDIR)/rootutils.o
	
cleantkmaine:
	@rm -f $(LIBDIR)/Squid.o $(LIBDIR)/tkmaterial.o $(BINDIR)/tkmaterial $(LIBDIR)/tklayout.o $(BINDIR)/tklayout $(BINDIR)/tkLayout $(TESTDIR)/testObjects $(TESTDIR)/rootwebTest $(BINDIR)/setup.bin
	
clean: cleanhit cleanexocom cleantkgeometry cleangeneral cleanelements cleanushers cleandressers cleanviz cleannaly cleanrootweb cleantkmaine 
	
doc: doxydoc
	
doxydoc:
	rm -rf $(DOXYDIR)/html
	doxygen $(DOCDIR)/tkdoc.doxy
	@echo "Created API documentation."

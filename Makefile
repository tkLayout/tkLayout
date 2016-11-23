# Define the REVISIONNUMBER variable to have it
# apearing in the root web page output
SVNREVISIONDEFINE=$(shell ./getRevisionDefine)
ROOTFLAGS=`root-config --cflags`
ROOTLIBDIR=`root-config --libdir`
ROOTLIBFLAGS=`root-config --libs`
ROOTLIBFLAGS+=-lHistPainter
ifneq ($(strip $(BOOST_INCLUDE)),)
INCLUDEFLAGS=-I$(BOOST_INCLUDE)
endif
ifneq ($(strip $(BOOST_LIB)),)
BOOSTLIBFLAGS=-L$(BOOST_LIB)
endif
BOOSTLIBFLAGS+=-L$(BOOST_LIB) -lboost_system$(BOOST_SUFFIX) -lboost_filesystem$(BOOST_SUFFIX) -lboost_program_options$(BOOST_SUFFIX)
GEOMLIBFLAG=-lGeom
GLIBFLAGS=`root-config --glibs`
INCLUDEFLAGS+=-Iinclude/
SRCDIR=src
INCDIR=include
LIBDIR=lib
BINDIR=bin
TESTDIR=test
DOCDIR=doc
DOXYDIR=doc/doxygen
#COMPILERFLAGS+=-Wall
COMPILERFLAGS+=-std=c++11
#COMPILERFLAGS+=-ggdb
COMPILERFLAGS+=-g
COMPILERFLAGS+=-fpermissive
COMPILERFLAGS+=-lstdc++
COMPILERFLAGS+=-fmax-errors=2
#COMPILERFLAGS+=-pg
#COMPILERFLAGS+=-Werror
#COMPILERFLAGS+=-O5
LINKERFLAGS+=-Wl,--copy-dt-needed-entries
#LINKERFLAGS+=-pg

OUT_DIR+=$(LIBDIR)
OUT_DIR+=$(BINDIR)
MKDIR_P = mkdir -p
.PHONY: directories

CXX := g++ $(JUST_DO_IT)
# COLORGCC
PATH := $(addprefix .:, $(PATH))
HASCOLOR = $(shell if test `which colorgcc 2> /dev/null`; then echo true; else echo false; fi)
ifneq ($(HASCOLOR),true)
HASCOLOR = $(shell if test -e colorgcc; then echo true; else echo false; fi)
endif

ifeq ($(HASCOLOR),true)
ifneq ($(COLOR), false)
CXX := colorgcc
endif
endif

COMP=$(CXX) $(COMPILERFLAGS) $(INCLUDEFLAGS)

LINK=$(CXX) $(LINKERFLAGS)

# All objects to compile
OBJS+=GraphVizCreator
OBJS+=PtError
OBJS+=Hit
OBJS+=CoordinateOperations
OBJS+=global_funcs
OBJS+=Property
OBJS+=Sensor
OBJS+=GeometricModule
OBJS+=DetectorModule
OBJS+=RodPair
OBJS+=Layer
OBJS+=Barrel
OBJS+=Disk
OBJS+=Endcap
OBJS+=Tracker
OBJS+=Ring
OBJS+=SimParms
OBJS+=Bag
OBJS+=SummaryTable
OBJS+=AnalyzerVisitor
OBJS+=PtErrorAdapter
OBJS+=PlotDrawer
OBJS+=Polygon3d
OBJS+=messageLogger
OBJS+=MatParser
OBJS+=Extractor
OBJS+=XMLWriter
OBJS+=IrradiationMap
OBJS+=IrradiationMapsManager
OBJS+=MaterialBudget
OBJS+=MaterialTable
OBJS+=MaterialProperties
OBJS+=InactiveSurfaces
OBJS+=ModuleCap
OBJS+=InactiveElement
OBJS+=InactiveRing
OBJS+=InactiveTube
OBJS+=Materialway
OBJS+=MaterialTab
OBJS+=MaterialObject
OBJS+=WeightDistributionGrid
OBJS+=ConversionStation
OBJS+=SupportStructure
OBJS+=Usher
OBJS+=MatCalc
OBJS+=MatCalcDummy
OBJS+=Vizard
OBJS+=tk2CMSSW
OBJS+=Analyzer
OBJS+=mainConfigHandler
OBJS+=Squid
OBJS+=rootweb
OBJS+=Palette
OBJS+=StopWatch

ANALYZERVISITORS+=Bandwidth
ANALYZERVISITORS+=IrradiationPower
ANALYZERVISITORS+=MaterialBillAnalyzer
ANALYZERVISITORS+=TriggerDistanceTuningPlots
ANALYZERVISITORS+=TriggerFrequency
ANALYZERVISITORS+=TriggerProcessorBandwidth

EXES+=tklayout
EXES+=setup
EXES+=diskPlace

OBJECTFILES=$(addsuffix .o,$(addprefix ${LIBDIR}/,${OBJS}))
ANALYZERVISITORFILES=$(addsuffix .o,$(addprefix ${LIBDIR}/AnalyzerVisitors/,${ANALYZERVISITORS}))
EXEFILES=$(addprefix ${BINDIR}/,${EXES})
EXELIBFILES=$(addsuffix .o,$(addprefix ${LIBDIR}/,${EXES}))

all: directories all_exes
	@echo "Full build successful."

# Executables can be compiled just by calling their name
$(EXES) : % : $(BINDIR)/%

directories: ${OUT_DIR}
${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

install:
	./install.sh

all_exes: ${EXEFILES}

clean_exes:
	rm -f ${EXEFILES}

all_objects: ${OBJECTFILES}

clean_objects:
	rm -f ${OBJECTFILES}

all_analyzerVisitors: ${ANALYZERVISITORFILES}

clean_analyzerVisitors:
	rm -f ${ANALYZERVISITORFILES}

# General rule to build objects
$(LIBDIR)/%.o: $(SRCDIR)/%.cpp $(INCDIR)/%.h
	$(COMP) $(ROOTFLAGS) -o $@ -c $<
	@echo "Object $@ built"

# Visitors in their own directory
$(LIBDIR)/AnalyzerVisitors/%.o: $(SRCDIR)/AnalyzerVisitors/%.cpp $(INCDIR)/AnalyzerVisitors/%.h
	@echo "Building AnalyzerVisitor $@..."
	mkdir -p $(LIBDIR)/AnalyzerVisitors
	$(COMP) $(ROOTFLAGS) -c -o $@ $<
	@echo "Built AnalyzerVisitor $@"

# Special rule for objects with a "main"
$(EXELIBFILES): $(LIBDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMP) $(ROOTFLAGS) -o $@ -c $<
	@echo "Main object $@ built"

# General rule to build executables
$(BINDIR)/%: $(LIBDIR)/%.o $(OBJECTFILES) $(ANALYZERVISITORFILES) getRevisionDefine
	# Let's make the revision object first
	$(COMP) $(SVNREVISIONDEFINE) -c $(SRCDIR)/SvnRevision.cpp -o $(LIBDIR)/SvnRevision.o
	# Now we just ahve to link standard objects, revision and main object
	$(LINK) $< $(OBJECTFILES) $(ANALYZERVISITORFILES) $(LIBDIR)/SvnRevision.o \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $@
	@echo "Executable $@ built"

# Speacial rules for light-weight executables
$(BINDIR)/diskPlace: $(SRCDIR)/diskPlace.cpp
	$(COMP) $(SRCDIR)/diskPlace.cpp -lm -o $(BINDIR)/diskPlace


# Other stuff for testing


testObjects: $(TESTDIR)/testObjects
$(TESTDIR)/testObjects: $(TESTDIR)/testObjects.cpp $(LIBDIR)/module.o $(LIBDIR)/layer.o
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/messageLogger.o $(TESTDIR)/testObjects.cpp \
	$(LIBDIR)/ptError.o $(LIBDIR)/moduleType.o \
	$(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o $(TESTDIR)/testObjects

testGraphVizCreator: $(TESTDIR)/testGraphVizCreator
$(TESTDIR)/testGraphVizCreator: $(TESTDIR)/testGraphVizCreator.cpp $(LIBDIR)/GraphVizCreator.o
	g++ $(COMPILERFLAGS) $(INCLUDEFLAGS) $(LIBDIR)/GraphVizCreator.o $(TESTDIR)/testGraphVizCreator.cpp -o $(TESTDIR)/testGraphVizCreator

rootwebTest: $(TESTDIR)/rootwebTest
$(TESTDIR)/rootwebTest: $(TESTDIR)/rootwebTest.cpp $(LIBDIR)/mainConfigHandler.o $(LIBDIR)/rootweb.o
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/mainConfigHandler.o $(LIBDIR)/rootweb.o $(TESTDIR)/rootwebTest.cpp $(ROOTLIBFLAGS) $(BOOSTLIBFLAGS) -o $(TESTDIR)/rootwebTest

test: $(TESTDIR)/ModuleTest

$(TESTDIR)/%: $(SRCDIR)/Tests/%.cpp $(INCDIR)/Tests/%.h
	@echo "Building target $@..."
	$(COMP) $(ROOTFLAGS) $(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAGS) -o $@ $< $(SRCDIR)/DetectorModule.cpp $(SRCDIR)/GeometricModule.cpp $(SRCDIR)/Sensor.cpp $(SRCDIR)/global_funcs.cpp $(SRCDIR)/Polygon3d.cpp
	@echo "Built target $@"

clean: clean_exes clean_objects clean_analyzerVisitors

doc: doxydoc

doxydoc:
	rm -rf $(DOXYDIR)/html
	doxygen $(DOCDIR)/tkdoc.doxy
	@echo "Created API documentation."

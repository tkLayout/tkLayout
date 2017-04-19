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
OUT_DIR+=$(LIBDIR)/AnalyzerVisitors
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
OBJS+=Analyzer
OBJS+=AnalyzerTools
OBJS+=AnalyzerVisitor
OBJS+=Bag
OBJS+=Barrel
OBJS+=ConversionStation
OBJS+=CoordinateOperations
OBJS+=DetectorModule
OBJS+=Disk
OBJS+=Endcap
OBJS+=Extractor
OBJS+=GeometricModule
OBJS+=global_funcs
OBJS+=GraphVizCreator
OBJS+=Histo
OBJS+=Hit
OBJS+=HitNew
OBJS+=InactiveElement
OBJS+=InactiveRing
OBJS+=InactiveSurfaces
OBJS+=InactiveTube
OBJS+=IrradiationMap
OBJS+=IrradiationMapsManager
OBJS+=Layer
OBJS+=MainConfigHandler
OBJS+=MatCalc
OBJS+=MatCalcDummy
OBJS+=MaterialBudget
OBJS+=MaterialObject
OBJS+=MaterialProperties
OBJS+=MaterialTab
OBJS+=MaterialTable
OBJS+=Materialway
OBJS+=MatParser
OBJS+=MessageLogger
OBJS+=ModuleCap
OBJS+=Module
OBJS+=Palette
OBJS+=PlotDrawer
OBJS+=Polygon3d
OBJS+=Property
OBJS+=PtError
OBJS+=PtErrorAdapter
OBJS+=ReportIrradiation
OBJS+=ReportModuleCount
OBJS+=Ring
OBJS+=RodPair
OBJS+=RootWeb
OBJS+=Sensor
OBJS+=SimParms
OBJS+=Squid
OBJS+=StopWatch
OBJS+=SummaryTable
OBJS+=SupportStructure
OBJS+=tk2CMSSW
OBJS+=Tracker
OBJS+=TrackNew
OBJS+=Usher
OBJS+=Vizard
OBJS+=WeightDistributionGrid
OBJS+=XMLWriter

ANALYZERVISITORS+=Bandwidth
ANALYZERVISITORS+=GeometricInfo
ANALYZERVISITORS+=IrradiationPower
ANALYZERVISITORS+=MaterialBillAnalyzer
ANALYZERVISITORS+=TriggerDistanceTuningPlots
ANALYZERVISITORS+=TriggerFrequency
ANALYZERVISITORS+=TriggerProcessorBandwidth
ANALYZERVISITORS+=ModuleCount

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
	@rm -f ${EXEFILES} ${EXELIBFILES}

all_objects: ${OBJECTFILES}

clean_objects:
	@rm -f ${OBJECTFILES} lib/SvnRevision.o

all_analyzerVisitors: ${ANALYZERVISITORFILES}

clean_analyzerVisitors:
	@rm -f ${ANALYZERVISITORFILES}

# General rule to build objects
$(LIBDIR)/%.o: $(SRCDIR)/%.cc $(INCDIR)/%.hh
	$(COMP) $(ROOTFLAGS) -o $@ -c $<
	@echo "Object $@ built"

# Visitors in their own directory
$(LIBDIR)/AnalyzerVisitors/%.o: $(SRCDIR)/AnalyzerVisitors/%.cc $(INCDIR)/AnalyzerVisitors/%.hh
	@echo "Building AnalyzerVisitor $@..."
	mkdir -p $(LIBDIR)/AnalyzerVisitors
	$(COMP) $(ROOTFLAGS) -c -o $@ $<
	@echo "Built AnalyzerVisitor $@"

# Special rule for objects with a "main"
$(EXELIBFILES): $(LIBDIR)/%.o: $(SRCDIR)/%.cc
	$(COMP) $(ROOTFLAGS) -o $@ -c $<
	@echo "Main object $@ built"

# General rule to build executables
$(BINDIR)/%: $(LIBDIR)/%.o $(OBJECTFILES) $(ANALYZERVISITORFILES) getRevisionDefine
	# Let's make the revision object first
	$(COMP) $(SVNREVISIONDEFINE) -c $(SRCDIR)/SvnRevision.cc -o $(LIBDIR)/SvnRevision.o
	# Now we just have to link standard objects, revision and main object
	$(LINK) $< $(OBJECTFILES) $(ANALYZERVISITORFILES) $(LIBDIR)/SvnRevision.o \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $@
	@echo "Executable $@ built"

# Speacial rules for light-weight executables
$(BINDIR)/diskPlace: $(SRCDIR)/diskPlace.cc
	$(COMP) $(SRCDIR)/diskPlace.cc -lm -o $(BINDIR)/diskPlace

$(BINDIR)/setup: $(LIBDIR)/MainConfigHandler.o $(LIBDIR)/global_funcs.o $(LIBDIR)/GraphVizCreator.o $(SRCDIR)/setup.cc
	$(COMP) $(LINKERFLAGS) $(LIBDIR)/MainConfigHandler.o $(LIBDIR)/global_funcs.o $(LIBDIR)/GraphVizCreator.o $(SRCDIR)/setup.cc \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $(BINDIR)/setup

# Clean and documentation
clean: clean_exes clean_objects clean_analyzerVisitors

doc: doxydoc

doxydoc:
	rm -rf $(DOXYDIR)/html
	doxygen $(DOCDIR)/tkdoc.doxy
	@echo "Created API documentation."






# Other stuff for testing
testGraphVizCreator: $(TESTDIR)/testGraphVizCreator
$(TESTDIR)/testGraphVizCreator: $(TESTDIR)/testGraphVizCreator.cc $(LIBDIR)/GraphVizCreator.o
	g++ $(COMPILERFLAGS) $(INCLUDEFLAGS) $(LIBDIR)/GraphVizCreator.o $(TESTDIR)/testGraphVizCreator.cc -o $(TESTDIR)/testGraphVizCreator

rootwebTest: $(TESTDIR)/rootwebTest
$(TESTDIR)/rootwebTest: $(TESTDIR)/rootwebTest.cc $(LIBDIR)/MainConfigHandler.o $(LIBDIR)/rootweb.o
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/MainConfigHandler.o $(LIBDIR)/rootweb.o $(TESTDIR)/rootwebTest.cc $(ROOTLIBFLAGS) $(BOOSTLIBFLAGS) -o $(TESTDIR)/rootwebTest

test: $(TESTDIR)/ModuleTest

$(TESTDIR)/%: $(SRCDIR)/Tests/%.cc $(INCDIR)/Tests/%.hh
	@echo "Building target $@..."
	$(COMP) $(ROOTFLAGS) $(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAGS) -o $@ $< $(SRCDIR)/DetectorModule.cc $(SRCDIR)/GeometricModule.cc $(SRCDIR)/Sensor.cc $(SRCDIR)/global_funcs.cc $(SRCDIR)/Polygon3d.cc
	@echo "Built target $@"



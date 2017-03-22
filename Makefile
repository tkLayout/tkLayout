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

all: directories tklayout setup diskplace
	@echo "Full build successful."

directories: ${OUT_DIR}
${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

install:
	./install.sh

# Pt computation
$(LIBDIR)/GraphVizCreator.o: $(SRCDIR)/GraphVizCreator.cpp $(INCDIR)/GraphVizCreator.hh
	$(COMP) -c -o $(LIBDIR)/GraphVizCreator.o $(SRCDIR)/GraphVizCreator.cpp

$(LIBDIR)/ptError.o: $(SRCDIR)/ptError.cpp $(INCDIR)/ptError.h
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/ptError.o $(SRCDIR)/ptError.cpp

$(BINDIR)/tunePtParam: $(SRCDIR)/tunePtParam.cpp $(LIBDIR)/ptError.o
	$(LINK) $(LIBDIR)/ptError.o $(SRCDIR)/tunePtParam.cpp \
	$(ROOTLIBFLAGS) $(ROOTFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $(BINDIR)/tunePtParam

#TRACKS
hit: $(LIBDIR)/hit.o
	@echo "Built target 'hit'."

$(LIBDIR)/CoordinateOperations.o: $(SRCDIR)/CoordinateOperations.cpp $(INCDIR)/CoordinateOperations.h
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/CoordinateOperations.o $(SRCDIR)/CoordinateOperations.cpp

$(LIBDIR)/hit.o: $(SRCDIR)/hit.cpp $(INCDIR)/hit.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/hit.o $(SRCDIR)/hit.cpp

$(LIBDIR)/global_funcs.o: $(SRCDIR)/global_funcs.cpp $(INCDIR)/global_funcs.h
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/global_funcs.o $(SRCDIR)/global_funcs.cpp
	
$(LIBDIR)/HitNew.o: $(SRCDIR)/HitNew.cc $(INCDIR)/HitNew.h
	@echo "Building target HitNew.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/HitNew.o $(SRCDIR)/HitNew.cc 
	@echo "Built target HitNew.o"
	
$(LIBDIR)/TrackNew.o: $(SRCDIR)/TrackNew.cc $(INCDIR)/TrackNew.h
	@echo "Building target TrackNew.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/TrackNew.o $(SRCDIR)/TrackNew.cc 
	@echo "Built target TrackNew.o"

$(LIBDIR)/Property.o: $(SRCDIR)/Property.cpp $(INCDIR)/Property.h
	@echo "Building target Property.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Property.o $(SRCDIR)/Property.cpp 
	@echo "Built target Property.o"

$(LIBDIR)/Sensor.o: $(SRCDIR)/Sensor.cpp $(INCDIR)/Sensor.h
	@echo "Building target Sensor.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Sensor.o $(SRCDIR)/Sensor.cpp 
	@echo "Built target Sensor.o"

$(LIBDIR)/GeometricModule.o: $(SRCDIR)/GeometricModule.cpp $(INCDIR)/GeometricModule.h
	@echo "Building target GeometricModule.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/GeometricModule.o $(SRCDIR)/GeometricModule.cpp 
	@echo "Built target GeometricModule.o"

$(LIBDIR)/DetectorModule.o: $(SRCDIR)/DetectorModule.cpp $(INCDIR)/DetectorModule.h
	@echo "Building target DetectorModule.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/DetectorModule.o $(SRCDIR)/DetectorModule.cpp 
	@echo "Built target DetectorModule.o"

$(LIBDIR)/RodPair.o: $(SRCDIR)/RodPair.cpp $(INCDIR)/RodPair.h
	@echo "Building target RodPair.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/RodPair.o $(SRCDIR)/RodPair.cpp 
	@echo "Built target RodPair.o"

$(LIBDIR)/Layer.o: $(SRCDIR)/Layer.cpp $(INCDIR)/Layer.h
	@echo "Building target Layer.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Layer.o $(SRCDIR)/Layer.cpp 
	@echo "Built target Layer.o"

$(LIBDIR)/Barrel.o: $(SRCDIR)/Barrel.cpp $(INCDIR)/Barrel.h
	@echo "Building target Barrel.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Barrel.o $(SRCDIR)/Barrel.cpp 
	@echo "Built target Barrel.o"

$(LIBDIR)/Ring.o: $(SRCDIR)/Ring.cpp $(INCDIR)/Ring.h
	@echo "Building target Ring.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Ring.o $(SRCDIR)/Ring.cpp 
	@echo "Built target Ring.o"

$(LIBDIR)/Disk.o: $(SRCDIR)/Disk.cpp $(INCDIR)/Disk.h
	@echo "Building target Disk.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Disk.o $(SRCDIR)/Disk.cpp 
	@echo "Built target Disk.o"

$(LIBDIR)/Endcap.o: $(SRCDIR)/Endcap.cpp $(INCDIR)/Endcap.h
	@echo "Building target Endcap.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Endcap.o $(SRCDIR)/Endcap.cpp 
	@echo "Built target Endcap.o"

$(LIBDIR)/Tracker.o: $(SRCDIR)/Tracker.cpp $(INCDIR)/Tracker.h
	@echo "Building target Tracker.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Tracker.o $(SRCDIR)/Tracker.cpp 
	@echo "Built target Tracker.o"

$(LIBDIR)/SimParms.o: $(SRCDIR)/SimParms.cpp $(INCDIR)/SimParms.h
	@echo "Building target SimParms.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/SimParms.o $(SRCDIR)/SimParms.cpp 
	@echo "Built target SimParms.o"

$(LIBDIR)/Bag.o: $(SRCDIR)/Bag.cpp $(INCDIR)/Bag.h
	@echo "Building target Bag.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Bag.o $(SRCDIR)/Bag.cpp 
	@echo "Built target Bag.o"

$(LIBDIR)/SummaryTable.o: $(SRCDIR)/SummaryTable.cpp $(INCDIR)/SummaryTable.h
	@echo "Building target SummaryTable.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/SummaryTable.o $(SRCDIR)/SummaryTable.cpp 
	@echo "Built target SummaryTable.o"

$(LIBDIR)/AnalyzerVisitors/%.o: $(SRCDIR)/AnalyzerVisitors/%.cpp $(INCDIR)/AnalyzerVisitors/%.h
	@echo "Building target $@..."
	mkdir -p $(LIBDIR)/AnalyzerVisitors
	$(COMP) $(ROOTFLAGS) -c -o $@ $< 
	@echo "Built target $@"

$(LIBDIR)/AnalyzerVisitor.o: $(SRCDIR)/AnalyzerVisitor.cpp $(INCDIR)/AnalyzerVisitor.h
	@echo "Building target AnalyzerVisitor.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/AnalyzerVisitor.o $(SRCDIR)/AnalyzerVisitor.cpp 
	@echo "Built target AnalyzerVisitor.o"

$(LIBDIR)/PtErrorAdapter.o: $(SRCDIR)/PtErrorAdapter.cpp $(INCDIR)/PtErrorAdapter.h
	@echo "Building target PtErrorAdapter.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/PtErrorAdapter.o $(SRCDIR)/PtErrorAdapter.cpp 
	@echo "Built target PtErrorAdapter.o"

$(LIBDIR)/ReportModuleCount.o: $(SRCDIR)/ReportModuleCount.cpp $(INCDIR)/ReportModuleCount.hh
	@echo "Building target ReportModuleCount.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/ReportModuleCount.o $(SRCDIR)/ReportModuleCount.cpp
	@echo "Built target ReportModuleCount.o"

$(LIBDIR)/PlotDrawer.o: $(SRCDIR)/PlotDrawer.cpp $(INCDIR)/PlotDrawer.h
	@echo "Building target PlotDrawer.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/PlotDrawer.o $(SRCDIR)/PlotDrawer.cpp
	@echo "Built target PlotDrawer.o"

$(LIBDIR)/Polygon3d.o: $(SRCDIR)/Polygon3d.cpp $(INCDIR)/Polygon3d.h
	@echo "Building target Polygon3d.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Polygon3d.o $(SRCDIR)/Polygon3d.cpp
	@echo "Built target Polygon3d.o"

$(LIBDIR)/TrackShooter.o: $(SRCDIR)/TrackShooter.cpp $(INCDIR)/TrackShooter.h
	@echo "Building target TrackShooter.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/TrackShooter.o $(SRCDIR)/TrackShooter.cpp
	@echo "Built target TrackShooter.o"

$(LIBDIR)/messageLogger.o: $(SRCDIR)/messageLogger.cpp $(INCDIR)/messageLogger.h
	$(COMP) -c -o $(LIBDIR)/messageLogger.o $(SRCDIR)/messageLogger.cpp

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

$(LIBDIR)/IrradiationMap.o: $(SRCDIR)/IrradiationMap.cpp $(INCDIR)/IrradiationMap.h
	@echo "Building target IrradiationMap.o..."
	$(COMP) -c -o $(LIBDIR)/IrradiationMap.o $(SRCDIR)/IrradiationMap.cpp
	@echo "Built target IrradiationMap.o"

$(LIBDIR)/IrradiationMapsManager.o: $(SRCDIR)/IrradiationMapsManager.cpp $(INCDIR)/IrradiationMapsManager.h
	@echo "Building target IrradiationMap.o..."
	$(COMP) -c -o $(LIBDIR)/IrradiationMapsManager.o $(SRCDIR)/IrradiationMapsManager.cpp
	@echo "Built target IrradiationMap.o"

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

#MATERIALWAY
materialways: $(LIBDIR)/Materialway.o $(LIBDIR)/MaterialTab.o $(LIBDIR)/WeightDistributionGrid.o $(LIBDIR)/MaterialObject.o $(LIBDIR)/ConversionStation.o $(LIBDIR)/SupportStructure.o
	@echo "Built target 'Materialways'."

$(LIBDIR)/Materialway.o: $(SRCDIR)/Materialway.cpp $(INCDIR)/Materialway.h
	@echo "Building target Materialway.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Materialway.o $(SRCDIR)/Materialway.cpp
	@echo "Built target Materialway.o"

$(LIBDIR)/MaterialTab.o: $(SRCDIR)/MaterialTab.cpp $(INCDIR)/MaterialTab.h
	@echo "Building target MaterialTab.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/MaterialTab.o $(SRCDIR)/MaterialTab.cpp
	@echo "Built target MaterialTab.o"

$(LIBDIR)/MaterialObject.o: $(SRCDIR)/MaterialObject.cpp $(INCDIR)/MaterialObject.h
	@echo "Building target MaterialObject.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/MaterialObject.o $(SRCDIR)/MaterialObject.cpp
	@echo "Built target MaterialObject.o"

$(LIBDIR)/WeightDistributionGrid.o: $(SRCDIR)/WeightDistributionGrid.cpp $(INCDIR)/WeightDistributionGrid.h
	@echo "Building target WeightDistributionGrid.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/WeightDistributionGrid.o $(SRCDIR)/WeightDistributionGrid.cpp
	@echo "Built target WeightDistributionGrid.o"

$(LIBDIR)/ConversionStation.o: $(SRCDIR)/ConversionStation.cpp $(INCDIR)/ConversionStation.h
	@echo "Building target ConversionStation.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/ConversionStation.o $(SRCDIR)/ConversionStation.cpp
	@echo "Built target ConversionStation.o"

$(LIBDIR)/SupportStructure.o: $(SRCDIR)/SupportStructure.cpp $(INCDIR)/SupportStructure.h
	@echo "Building target SupportStructure.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/SupportStructure.o $(SRCDIR)/SupportStructure.cpp
	@echo "Built target SupportStructure.o"

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

$(LIBDIR)/Analyzer.o: $(SRCDIR)/Analyzer.cpp $(INCDIR)/Analyzer.h
	@echo "Building target Analyzer.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Analyzer.o $(SRCDIR)/Analyzer.cpp
	@echo "Built target Analyzer.o"

#SQUID
squid: $(LIBDIR)/Squid.o
	@echo "Built target 'squid'."

$(LIBDIR)/mainConfigHandler.o: $(SRCDIR)/mainConfigHandler.cpp $(INCDIR)/mainConfigHandler.h
	@echo "Building target mainConfigHandler.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/mainConfigHandler.o $(SRCDIR)/mainConfigHandler.cpp
	@echo "Built target mainConfigHandler.o"

$(LIBDIR)/Squid.o: $(SRCDIR)/Squid.cc $(INCDIR)/Squid.h
	@echo "Building target Squid.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Squid.o $(SRCDIR)/Squid.cc
	@echo "Built target Squid.o"

#ROOT-related stuff
$(LIBDIR)/rootweb.o: $(SRCDIR)/rootweb.cpp $(INCDIR)/rootweb.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/rootweb.o $(SRCDIR)/rootweb.cpp

$(LIBDIR)/Palette.o: $(SRCDIR)/Palette.cc  $(INCDIR)/Palette.h
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Palette.o $(SRCDIR)/Palette.cc

# Helper objects
$(LIBDIR)/StopWatch.o: $(SRCDIR)/StopWatch.cpp $(INCDIR)/StopWatch.h
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/StopWatch.o $(SRCDIR)/StopWatch.cpp

#$(LIBDIR)/rootutils.o: $(SRCDIR)/rootutils.cpp $(INCDIR)/rootutils.h
#	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/rootutils.o $(SRCDIR)/rootutils.cpp

# UTILITIES
setup: $(BINDIR)/setup.bin
	@echo "setup built"

$(BINDIR)/setup.bin: $(LIBDIR)/mainConfigHandler.o $(LIBDIR)/global_funcs.o $(LIBDIR)/GraphVizCreator.o $(SRCDIR)/setup.cpp
	$(COMP) $(LINKERFLAGS) $(LIBDIR)/mainConfigHandler.o $(LIBDIR)/global_funcs.o $(LIBDIR)/GraphVizCreator.o $(SRCDIR)/setup.cpp \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $(BINDIR)/setup.bin

houghtrack: $(BINDIR)/houghtrack
	@echo "houghtrack built"

$(LIBDIR)/Histo.o: $(SRCDIR)/Histo.cpp
	@echo "Building target Histo.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Histo.o $(SRCDIR)/Histo.cpp
	@echo "Built target Histo.o"

$(BINDIR)/houghtrack: $(LIBDIR)/TrackShooter.o $(LIBDIR)/module.o $(LIBDIR)/moduleType.o $(LIBDIR)/global_funcs.o $(LIBDIR)/ptError.o $(LIBDIR)/Histo.o $(SRCDIR)/HoughTrack.cpp $(INCDIR)/HoughTrack.h
	$(COMP) $(LINKERFLAGS) $(ROOTFLAGS) $(LIBDIR)/TrackShooter.o $(LIBDIR)/module.o $(LIBDIR)/moduleType.o $(LIBDIR)/global_funcs.o $(LIBDIR)/ptError.o $(LIBDIR)/Histo.o $(SRCDIR)/HoughTrack.cpp \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $(BINDIR)/houghtrack


#FINAL
tklayout: $(BINDIR)/tklayout
	@echo "tklayout built"

diskplace: $(BINDIR)/diskplace
	@echo "diskplace built"

delphize: $(BINDIR)/delphize
	@echo "delphize built"

tunePtParam: $(BINDIR)/tunePtParam
	@echo "tunePtParam built"

$(BINDIR)/diskplace: $(SRCDIR)/diskPlace.cpp
	$(COMP) $(SRCDIR)/diskPlace.cpp -lm -o $(BINDIR)/diskplace

$(BINDIR)/tklayout: $(LIBDIR)/tklayout.o $(LIBDIR)/CoordinateOperations.o $(LIBDIR)/hit.o $(LIBDIR)/HitNew.o $(LIBDIR)/TrackNew.o $(LIBDIR)/global_funcs.o $(LIBDIR)/Polygon3d.o \
	$(LIBDIR)/Property.o \
	$(LIBDIR)/Sensor.o $(LIBDIR)/GeometricModule.o $(LIBDIR)/DetectorModule.o $(LIBDIR)/RodPair.o $(LIBDIR)/Layer.o $(LIBDIR)/Barrel.o $(LIBDIR)/Ring.o $(LIBDIR)/Disk.o $(LIBDIR)/Endcap.o $(LIBDIR)/Tracker.o $(LIBDIR)/SimParms.o \
  $(LIBDIR)/AnalyzerVisitors/MaterialBillAnalyzer.o \
	$(LIBDIR)/AnalyzerVisitors/TriggerFrequency.o $(LIBDIR)/AnalyzerVisitors/Bandwidth.o $(LIBDIR)/AnalyzerVisitors/IrradiationPower.o $(LIBDIR)/AnalyzerVisitors/TriggerProcessorBandwidth.o $(LIBDIR)/AnalyzerVisitors/TriggerDistanceTuningPlots.o \
	$(LIBDIR)/AnalyzerVisitor.o $(LIBDIR)/Bag.o $(LIBDIR)/SummaryTable.o $(LIBDIR)/PtErrorAdapter.o $(LIBDIR)/Analyzer.o $(LIBDIR)/ptError.o \
  $(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o \
	$(LIBDIR)/XMLWriter.o $(LIBDIR)/IrradiationMap.o $(LIBDIR)/IrradiationMapsManager.o $(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialProperties.o \
	$(LIBDIR)/ModuleCap.o  $(LIBDIR)/InactiveSurfaces.o  $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o \
	$(LIBDIR)/InactiveTube.o $(LIBDIR)/Usher.o $(LIBDIR)/Materialway.o $(LIBDIR)/MaterialTab.o $(LIBDIR)/WeightDistributionGrid.o $(LIBDIR)/MaterialObject.o $(LIBDIR)/ConversionStation.o $(LIBDIR)/SupportStructure.o $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o $(LIBDIR)/PlotDrawer.o $(LIBDIR)/ReportModuleCount.o \
	$(LIBDIR)/Vizard.o $(LIBDIR)/tk2CMSSW.o $(LIBDIR)/Squid.o $(LIBDIR)/rootweb.o $(LIBDIR)/mainConfigHandler.o \
	$(LIBDIR)/messageLogger.o $(LIBDIR)/Palette.o $(LIBDIR)/StopWatch.o $(LIBDIR)/GraphVizCreator.o getRevisionDefine
	#
	# Let's make the revision object first
	$(COMP) $(SVNREVISIONDEFINE) -c $(SRCDIR)/SvnRevision.cpp -o $(LIBDIR)/SvnRevision.o
	#
	# And compile the executable by linking the revision too
	$(LINK)	$(LIBDIR)/CoordinateOperations.o $(LIBDIR)/hit.o $(LIBDIR)/HitNew.o $(LIBDIR)/TrackNew.o $(LIBDIR)/global_funcs.o $(LIBDIR)/Polygon3d.o \
	$(LIBDIR)/Property.o \
	$(LIBDIR)/Sensor.o $(LIBDIR)/GeometricModule.o $(LIBDIR)/DetectorModule.o $(LIBDIR)/RodPair.o $(LIBDIR)/Layer.o $(LIBDIR)/Barrel.o $(LIBDIR)/Ring.o $(LIBDIR)/Disk.o $(LIBDIR)/Endcap.o $(LIBDIR)/Tracker.o $(LIBDIR)/SimParms.o \
  $(LIBDIR)/AnalyzerVisitors/MaterialBillAnalyzer.o \
	$(LIBDIR)/AnalyzerVisitors/TriggerFrequency.o $(LIBDIR)/AnalyzerVisitors/Bandwidth.o $(LIBDIR)/AnalyzerVisitors/IrradiationPower.o $(LIBDIR)/AnalyzerVisitors/TriggerProcessorBandwidth.o $(LIBDIR)/AnalyzerVisitors/TriggerDistanceTuningPlots.o \
	$(LIBDIR)/AnalyzerVisitor.o $(LIBDIR)/Bag.o $(LIBDIR)/SummaryTable.o $(LIBDIR)/PtErrorAdapter.o $(LIBDIR)/Analyzer.o $(LIBDIR)/ptError.o \
	$(LIBDIR)/MatParser.o $(LIBDIR)/Extractor.o \
	$(LIBDIR)/XMLWriter.o $(LIBDIR)/IrradiationMap.o $(LIBDIR)/IrradiationMapsManager.o $(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialProperties.o \
	$(LIBDIR)/ModuleCap.o $(LIBDIR)/InactiveSurfaces.o $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o \
	$(LIBDIR)/InactiveTube.o $(LIBDIR)/Usher.o $(LIBDIR)/Materialway.o $(LIBDIR)/MaterialTab.o $(LIBDIR)/WeightDistributionGrid.o $(LIBDIR)/MaterialObject.o $(LIBDIR)/ConversionStation.o $(LIBDIR)/SupportStructure.o $(LIBDIR)/MatCalc.o $(LIBDIR)/MatCalcDummy.o $(LIBDIR)/PlotDrawer.o $(LIBDIR)/ReportModuleCount.o \
	$(LIBDIR)/Vizard.o $(LIBDIR)/tk2CMSSW.o $(LIBDIR)/Squid.o $(LIBDIR)/rootweb.o $(LIBDIR)/mainConfigHandler.o \
	$(LIBDIR)/messageLogger.o $(LIBDIR)/Palette.o $(LIBDIR)/StopWatch.o $(LIBDIR)/GraphVizCreator.o \
	$(LIBDIR)/SvnRevision.o \
	$(LIBDIR)/tklayout.o \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $(BINDIR)/tklayout

$(BINDIR)/delphize: $(LIBDIR)/delphize.o
	# And compile the executable by linking the revision too
	$(LINK)	$(LIBDIR)/delphize.o \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o $(BINDIR)/delphize

$(LIBDIR)/tklayout.o: $(SRCDIR)/tklayout.cpp
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/tklayout.o $(SRCDIR)/tklayout.cpp

$(LIBDIR)/delphize.o: $(SRCDIR)/delphize.cpp
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/delphize.o $(SRCDIR)/delphize.cpp

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


#CLEANUP
cleanall:
	@rm -rf $(LIBDIR)/*

clean: cleanall

doc: doxydoc

doxydoc:
	rm -rf $(DOXYDIR)/html
	doxygen $(DOCDIR)/tkdoc.doxy
	@echo "Created API documentation."

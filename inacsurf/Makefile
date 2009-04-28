
ROOTFLAGS=`root-config --cflags`
ROOTLIBDIR=`root-config --libdir`
ROOTLIBFLAGS=`root-config --libs`
BOOSTLIBFLAGS=-lboost_filesystem
GEOMLIBFLAG=-lGeom
GLIBFLAGS=`root-config --glibs`
INCLUDEFLAGS=-Iinc -I../include/
TKSRCDIR=../src
TKINCDIR=../include
SRCDIR=src
INCDIR=inc
LIBDIR=lib

COMP=g++ -Wall $(INCLUDEFLAGS)

all: tkgeometry general elements ushers dressers viz test demo
	@echo "Full build successful."
#TODO: define targets for all: tkgeometry, insur, densities, demo

#TKGEOMETRY
tkgeometry: $(LIBDIR)/configparser.o $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o
	@echo "Built target 'tkgeometry'."

$(LIBDIR)/configparser.o:	$(TKSRCDIR)/configparser.cpp $(TKINCDIR)/configparser.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/configparser.o $(TKSRCDIR)/configparser.cpp

$(LIBDIR)/module.o: $(TKSRCDIR)/module.cpp $(TKINCDIR)/module.hh
	@echo "Building target module.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/module.o $(TKSRCDIR)/module.cpp 
	@echo "Built target module.o"

$(LIBDIR)/layer.o: $(TKSRCDIR)/layer.cpp $(TKINCDIR)/layer.hh
	@echo "Building target layer.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/layer.o $(TKSRCDIR)/layer.cpp 
	@echo "Built target layer.o"

$(LIBDIR)/tracker.o: $(TKSRCDIR)/tracker.cpp $(TKINCDIR)/tracker.hh
	@echo "Building target tracker.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/tracker.o $(TKSRCDIR)/tracker.cpp 	
	@echo "Built target tracker.o"

#GENERAL
general: $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialProperties.o \
	$(LIBDIR)/InactiveSurfaces.o $(LIBDIR)/TrackerActions.o $(LIBDIR)/MatParser.o
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

$(LIBDIR)/TrackerActions.o: $(SRCDIR)/TrackerActions.cc $(INCDIR)/TrackerActions.h
	@echo "Building target TrackerActions.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/TrackerActions.o $(SRCDIR)/TrackerActions.cc
	@echo "Built target TrackerActions.o"

$(LIBDIR)/MatParser.o: $(SRCDIR)/MatParser.cc $(INCDIR)/MatParser.h
	@echo "Building target MatParser.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/MatParser.o $(SRCDIR)/MatParser.cc
	@echo "Built target MatParser.o"
	
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
dressers: $(LIBDIR)/MatCalc.o
	@echo "Built target 'dressers'."

$(LIBDIR)/MatCalc.o: $(SRCDIR)/MatCalc.cc $(INCDIR)/MatCalc.h
	@echo "Building target MatCalc.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/MatCalc.o $(SRCDIR)/MatCalc.cc
	@echo "Built target MatlCalc.o"

#VISUALISATION
viz: $(LIBDIR)/Vizard.o
	@echo "Built target 'viz'."
	
$(LIBDIR)/Vizard.o: $(SRCDIR)/Vizard.cc $(INCDIR)/Vizard.h
	@echo "Building target Vizard.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/Vizard.o $(SRCDIR)/Vizard.cc
	@echo "Built target Vizard.o"

#TESTS
test: $(LIBDIR)/MatCalcDummy.o
	@echo "Built target 'test'."

$(LIBDIR)/MatCalcDummy.o: $(SRCDIR)/MatCalcDummy.cc $(INCDIR)/MatCalcDummy.h
	@echo "Building target MatCalcDummy.o..."
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/MatCalcDummy.o $(SRCDIR)/MatCalcDummy.cc
	@echo "Built target MatCalcDummy.o"

#FINAL
demo: $(LIBDIR)/demo.o
	$(COMP) $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o \
	$(LIBDIR)/configparser.o $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MatParser.o \
	$(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialProperties.o $(LIBDIR)/TrackerActions.o \
	$(LIBDIR)/InactiveSurfaces.o  $(LIBDIR)/ModuleCap.o $(LIBDIR)/InactiveElement.o \
	$(LIBDIR)/InactiveRing.o $(LIBDIR)/InactiveTube.o $(LIBDIR)/Usher.o \
	$(LIBDIR)/MatCalc.o $(LIBDIR)/Vizard.o $(LIBDIR)/MatCalcDummy.o \
	$(LIBDIR)/demo.o \
	$(ROOTLIBFLAGS) $(GLIBFLAGS) $(BOOSTLIBFLAGS) $(GEOMLIBFLAG) \
	-o demo

$(LIBDIR)/demo.o: demo.cpp
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/demo.o demo.cpp

cleantkgeometry:
	-rm -f include/tkgeometry/*~ include/*~ *~ $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o $(LIBDIR)/configparser.o

cleangeneral:
	-rm -f $(LIBDIR)/MaterialBudget.o $(LIBDIR)/MaterialTable.o $(LIBDIR)/MaterialProperties.o \
	$(LIBDIR)/TrackerActions.o $(LIBDIR)/MatParser.o $(LIBDIR)/InactiveSurfaces.o

cleanelements:
	-rm -f $(LIBDIR)/ModuleCap.o $(LIBDIR)/InactiveElement.o $(LIBDIR)/InactiveRing.o $(LIBDIR)/InactiveTube.o

cleanushers:
	-rm -f $(LIBDIR)/Usher.o

cleandressers:
	-rm -f $(LIBDIR)/MatCalc.o
	
cleanviz:
	-rm -r $(LIBDIR)/Vizard.o

cleandemo:
	-rm -f $(LIBDIR)/demo.o ./demo

clean: cleantkgeometry cleangeneral cleanelements cleanushers cleandressers cleanviz cleandemo
	
#TODO: add final build targets to list
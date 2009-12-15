
ROOTFLAGS=`root-config --cflags`
ROOTLIBDIR=`root-config --libdir`
ROOTLIBFLAGS=`root-config --libs`
BOOSTLIBFLAGS=-lboost_filesystem
GEOMLIBFLAG=-lGeom
INCLUDEFLAGS=-Iinclude

LIBDIR=lib

COMP=g++ -Wall $(INCLUDEFLAGS)

all: tkLayout testObj Gui 

Gui:
	cd gui && qmake tkgeomgui.pro && make && cd ..

$(LIBDIR)/configparser.o:	src/configparser.cpp include/configparser.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/configparser.o src/configparser.cpp

$(LIBDIR)/module.o: src/module.cpp include/module.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/module.o src/module.cpp 

$(LIBDIR)/layer.o: src/layer.cpp include/layer.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/layer.o src/layer.cpp 

$(LIBDIR)/tracker.o: src/tracker.cpp include/tracker.hh
	$(COMP) $(ROOTFLAGS) -c -o $(LIBDIR)/tracker.o src/tracker.cpp 	

TrackerGeom: TrackerGeom.cpp $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o TrackerGeom.cpp $(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o TrackerGeom

tkLayout: TrackerGeom2.cpp $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o $(LIBDIR)/configparser.o
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o $(LIBDIR)/configparser.o TrackerGeom2.cpp $(BOOSTLIBFLAGS) $(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o tkLayout

testObj: testObjects.cpp $(LIBDIR)/module.o $(LIBDIR)/layer.o
	$(COMP) $(ROOTFLAGS) $(LIBDIR)/module.o $(LIBDIR)/layer.o testObjects.cpp $(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o testObj

clean:
	rm -f include/*~ *~ $(LIBDIR)/module.o $(LIBDIR)/layer.o $(LIBDIR)/tracker.o $(LIBDIR)/configparser.o TrackerGeom tkLayout tkGeometry.root testObj cmsTest

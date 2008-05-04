
ROOTFLAGS=`root-config --cflags`
ROOTLIBDIR=`root-config --libdir`
ROOTLIBFLAGS=`root-config --libs`
GEOMLIBFLAG=-lGeom
INCLUDEFLAGS=-Iinclude

COMP=g++ $(INCLUDEFLAGS)

all: TrackerGeom test

module.o: src/module.cpp include/module.hh
	$(COMP) $(ROOTFLAGS) -c -o module.o src/module.cpp 

layer.o: src/layer.cpp include/layer.hh
	$(COMP) $(ROOTFLAGS) -c -o layer.o src/layer.cpp 

tracker.o: src/tracker.cpp include/tracker.hh
	$(COMP) $(ROOTFLAGS) -c -o tracker.o src/tracker.cpp 	

TrackerGeom: TrackerGeom.cpp module.o layer.o tracker.o
	$(COMP) $(ROOTFLAGS) module.o layer.o tracker.o TrackerGeom.cpp $(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o TrackerGeom

test: testObjects.cpp module.o
	$(COMP) $(ROOTFLAGS) module.o layer.o testObjects.cpp $(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o test

clean:
	rm -f include/*~ *~ module.o layer.o tracker.o TrackerGeom tkGeometry.root test cmsTest

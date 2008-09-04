
ROOTFLAGS=`root-config --cflags`
ROOTLIBDIR=`root-config --libdir`
ROOTLIBFLAGS=`root-config --libs`
GEOMLIBFLAG=-lGeom
INCLUDEFLAGS=-Iinclude

COMP=g++ $(INCLUDEFLAGS)

all: TrackerGeom testObj testConfig

configparser.o:	src/configparser.cpp include/configparser.hh
	$(COMP) $(ROOTFLAGS) -c -o configparser.o src/configparser.cpp

testConfig: configparser.o testConfig.cpp
	$(COMP) configparser.o testConfig.cpp -o testConfig

module.o: src/module.cpp include/module.hh
	$(COMP) $(ROOTFLAGS) -c -o module.o src/module.cpp 

layer.o: src/layer.cpp include/layer.hh
	$(COMP) $(ROOTFLAGS) -c -o layer.o src/layer.cpp 

tracker.o: src/tracker.cpp include/tracker.hh
	$(COMP) $(ROOTFLAGS) -c -o tracker.o src/tracker.cpp 	

TrackerGeom: TrackerGeom.cpp module.o layer.o tracker.o
	$(COMP) $(ROOTFLAGS) module.o layer.o tracker.o TrackerGeom.cpp $(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o TrackerGeom

testObj: testObjects.cpp module.o
	$(COMP) $(ROOTFLAGS) module.o layer.o testObjects.cpp $(ROOTLIBFLAGS) $(GEOMLIBFLAG) -o testObj

clean:
	rm -f include/*~ *~ module.o layer.o tracker.o configparser.o TrackerGeom tkGeometry.root testConfig testObj cmsTest

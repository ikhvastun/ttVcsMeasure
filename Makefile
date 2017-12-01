CXX=g++
RM=rm -f
CPPFLAGS=-std=c++0x -g $(shell root-config --cflags)
LDFLAGS=-g $(shell root-config --ldflags)
LDLIBS=$(shell root-config --libs) -ltmva
CFLAGS = `root-config --cflags --libs`

EXECUTABLE = main

SOURCES = src/BTagCalibrationStandalone.cc src/CMS_lumi.cc src/runCode.cc src/treeReader.cc src/analysisTools.cc src/eventSelection.cc
#SOURCES = BTagCalibrationStandalone.cc CMS_lumi.cc readTreeSync_OS.cc 
OBJS=$(subst .cc,.o,$(SOURCES))

all: main

main: $(OBJS)
	$(CXX) $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDLIBS) 

%.o: %.c 
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

.PHONY: clean

clean:
	$(RM) src/*.o *.so *.d
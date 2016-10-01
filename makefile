in=main

out=$(in)
bin_name=bin/$(out)
source_name=$(in).C

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)

all: 
	g++ -ggdb -std=c++11 $(source_name) ./src/source/*.cc -I ./src/include/ $(ROOTFLAGS) $(ROOTLIBS) -o lrGlass

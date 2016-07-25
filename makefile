in=main

out=$(in)
bin_name=bin/$(out)
source_name=$(in).C

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)

all: 
	g++ -O3 -std=c++11 $(source_name) source/*.cc -I include/ $(ROOTFLAGS) $(ROOTLIBS) -o lrGlass

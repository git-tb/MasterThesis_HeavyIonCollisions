.SILENT:

CXX       := "g++"
CXX_FLAGS := -std=c++17 -g -fopenmp

BIN     := bin
SRC     := src
INCLUDE := -Iinclude -Iboost_1_82_0 -IlibInterpolate/src -Ieigen-3.4.0/Eigen #DON'T FORGET "-I"-PREFIX

DIRLIBRARIES   := #-L"/usr/lib/x86_64-linux-gnu"#DON'T FORGET "-L"-PREFIX
LIBRARIES = -l"gsl" #-l"sfml-system" -l"sfml-network" -l"sfml-audio" -l"sfml-window" -l"sfml-graphics" #DON'T FORGET "-l"-PREFIX
EXECUTABLE  := main

run: build
	./$(BIN)/$(EXECUTABLE)

build:
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $(SRC)/$(EXECUTABLE).cpp -o $(BIN)/$(EXECUTABLE) $(DIRLIBRARIES) $(LIBRARIES)


CXX       := g++
CXX_FLAGS := -std=c++2a
# CXX_FLAGS += -g
# CXX_FLAGS += -O3
CXX_FLAGS += -DNDEBUG

EXE 	:= NO_EXECUTBALE_SPECIFIED
BIN     := bin
SRC     := src

ARGS	:=

INCLUDE :=
DIRLIBRARIES := 
LIBRARIES :=

INCLUDE += -Iinclude
INCLUDE += -Iboost_1_82_0
INCLUDE += -Icubature

DIRLIBRARIES += -Lboost_1_82_0/lib
DIRLIBRARIES += -Lcubature

LIBRARIES += -lboost_program_options
LIBRARIES += -lgsl
LIBRARIES += -lhcubature -lpcubature
run:
	./$(BIN)/$(EXE) $(ARGS)

build:
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $(SRC)/$(EXE).cpp -o $(BIN)/$(EXE) $(DIRLIBRARIES) $(LIBRARIES) -Wl,-rpath=$(PWD)/boost_1_82_0/lib/ -Wl,-rpath=$(PWD)/cubature/

checklibs:
	ldd ./bin/$(EXE)

# EXAMPLE CALL IS make run EXE=spec ARGS="--initpath=data/init_20240725_145025/initialfields_piplus.csv"



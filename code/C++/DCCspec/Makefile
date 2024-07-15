CXX       := g++
CXX_FLAGS := -std=c++2a
# CXX_FLAGS += -g
# CXX_FLAGS += -O3
CXX_FLAGS += -DNDEBUG

EXE 	:= main
BIN     := bin
SRC     := src

INCLUDE :=
DIRLIBRARIES := 
LIBRARIES :=

INCLUDE += -Iinclude
INCLUDE += -Iboost_1_82_0

DIRLIBRARIES += -Lboost_1_82_0/lib

LIBRARIES += -lboost_program_options
LIBRARIES += -lgsl
run:
	./$(BIN)/$(EXE)

build:
	$(CXX) $(CXX_FLAGS) $(INCLUDE) $(SRC)/$(EXE).cpp -o $(BIN)/$(EXE) $(DIRLIBRARIES) $(LIBRARIES) -Wl,-rpath=$(PWD)/boost_1_82_0/lib/

checklibs:
	ldd ./bin/$(EXE)



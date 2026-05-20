CXX     ?= g++
CXXFLAGS ?= -std=c++14 -O3 -Wall -Wno-sign-compare -Wno-unknown-warning-option \
            -Wno-deprecated-declarations

SRC     = src/ising_tri_twisted_parallelogram.cc
INC_DIR = include

ifeq ($(OS),Windows_NT)
  EXE = bin/ising_tri_twisted_parallelogram.exe
else
  EXE = bin/ising_tri_twisted_parallelogram
endif

.PHONY: all clean

all: $(EXE)

$(EXE): $(SRC) $(wildcard $(INC_DIR)/*.h)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $< -o $@
	@echo "Built: $@"

clean:
	$(RM) bin/ising_tri_twisted_parallelogram bin/ising_tri_twisted_parallelogram.exe
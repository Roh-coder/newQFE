CXX     = g++
CXXFLAGS = -std=c++14 -O3 -Wall -Wno-sign-compare -Wno-unknown-warning-option \
           -Wno-deprecated-declarations

SRC     = src/ising_tri_twisted_parallelogram.cc
INC_DIR = include

# Detect Windows (MinGW sets OS=Windows_NT)
ifeq ($(OS),Windows_NT)
  EXE = bin/ising_tri_twisted_parallelogram.exe
else
  EXE = bin/ising_tri_twisted_parallelogram
endif

.PHONY: all clean test

all: $(EXE)

$(EXE): $(SRC) $(wildcard $(INC_DIR)/*.h)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $< -o $@
	@echo "Built: $@"

clean:
	$(RM) bin/ising_tri_twisted_parallelogram bin/ising_tri_twisted_parallelogram.exe

# Run the unit-test suite for the speedup phases (S1/S2/S3).  Tests that
# need the compiled simulator are skipped automatically when $(EXE) is
# missing — run `make all test` for the full suite.
test:
	pytest tests/ -q


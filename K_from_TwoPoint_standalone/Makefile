CXX     = g++
CXXFLAGS = -std=c++14 -O3 -Wall -Wno-sign-compare -Wno-unknown-warning-option \
           -Wno-deprecated-declarations

SRC     = src/ising_tri_twisted_parallelogram.cc
EXE     = bin/ising_tri_twisted_parallelogram
INC_DIR = include

.PHONY: all clean

all: $(EXE)

$(EXE): $(SRC) $(wildcard $(INC_DIR)/*.h)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $< -o $@
	@echo "Built: $@"

clean:
	$(RM) $(EXE)

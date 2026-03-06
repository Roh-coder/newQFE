CXX=g++
CXXFLAGS=-std=c++14 -O3 -Wall -Wextra -Wno-sign-compare

INC_DIR=../include
BIN_DIR=../bin
TARGET=potts_sw4_scan
SRC=potts_sw4_scan.cc

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) $(SRC) -o $(BIN_DIR)/$(TARGET)

clean:
	$(RM) $(BIN_DIR)/$(TARGET)

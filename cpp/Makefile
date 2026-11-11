CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Iinclude

BUILD_DIR = build
BIN_DIR = bin

TARGETS = $(BIN_DIR)/evimo2_viewer $(BIN_DIR)/imo_detect $(BIN_DIR)/test_numerics

all: $(TARGETS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/evimo2_viewer: apps/viewer/main.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN_DIR)/imo_detect: apps/imo_detect/main.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN_DIR)/test_numerics: tests/test_numerics.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

test: $(BIN_DIR)/test_numerics
	./$(BIN_DIR)/test_numerics tests/data

clean:
	rm -rf $(BIN_DIR) $(BUILD_DIR) tests/test_numerics

.PHONY: all test clean

# Default build type
BUILD_TYPE ?= release

# Compiler options
CXX         := g++

ifeq ($(BUILD_TYPE),debug)
CXXFLAGS := -Wall -Wextra -O0 -ggdb
else
CXXFLAGS    := -Wall -Wextra -O3 -march=native
endif
CXXFLAGS    += -std=c++20
# CXXFLAGS    += -Wno-unused-parameter -Wno-unused-variable
# CXXFLAGS    += -Wno-unused-but-set-parameter -Wno-comment -Wno-unused-label
# CXXFLAGS    += -Wno-sign-compare -Wno-misleading-indentation
# CXXFLAGS    += -Wno-unused-but-set-variable -Wno-changes-meaning
# CXXFLAGS    += -fpermissive
LDFLAGS     := 

# External libraries
# ifeq ($(shell uname),Darwin)
# SDK_PATH    := $(shell xcrun --show-sdk-path)
# CXXFLAGS    += -isysroot $(SDK_PATH)
# CXXFLAGS    += $(shell pkg-config --cflags mpfr)
# CXXFLAGS    += -I$(shell nix eval --raw nixpkgs#libmpc)/include
# CXXFLAGS    += $(shell pkg-config --cflags arpack)
# # CXXFLAGS    += $(shell pkg-config --cflags cln)
# CXXFLAGS    += $(shell g++ -E -x c++ - -v < /dev/null 2>&1 | \
#   sed -n '/#include <...> search starts here:/,/End of search list./p' | \
#   tail -n +2 | head -n -1 | sed 's|^ \+||' | sed 's|^|-I|')
# endif


CXXFLAGS    += -I./include
ifeq ($(shell uname),Linux)
CXXFLAGS    += -I/home/roux/mpc-compiled/include
endif

# ifeq ($(shell uname),Darwin)
# LDFLAGS     += $(shell pkg-config --libs mpfr)
# LDFLAGS     += -L$(shell nix eval --raw nixpkgs#libmpc)/lib -lmpc
# LDFLAGS     += $(shell pkg-config --libs arpack)
# # LDFLAGS     += $(shell pkg-config --libs cln)
# else
LDFLAGS     += -lmpfr -lgmp -lmpc
ifeq ($(shell uname),Linux)
LDFLAGS     += -L/home/roux/mpc-compiled/lib
endif

# Directories
EXAMPLES_DIR := examples
TESTS_DIR    := tests
PROJECTS_DIR := projects
ifeq ($(BUILD_TYPE),release)
BUILD_DIR    := build/release
else
BUILD_DIR    := build/debug
endif

# Targets
EXAMPLES := OnLoops Ising
TESTS := testVector testMpfr testMatrices

# Targets
EXAMPLES := $(EXAMPLES_DIR)/OnLoops/$(BUILD_DIR)/OnLoops \
            $(EXAMPLES_DIR)/Ising/$(BUILD_DIR)/Ising

TESTS    := $(TESTS_DIR)/$(BUILD_DIR)/testVector \
            $(TESTS_DIR)/$(BUILD_DIR)/testMpfr \
            $(TESTS_DIR)/$(BUILD_DIR)/testMatrices \
            $(TESTS_DIR)/$(BUILD_DIR)/testBigFloats

# Default target
all: examples tests projects

# Group targets
examples: $(EXAMPLES)
tests:    $(TESTS)

# ----- Specific build rules for examples -----
$(EXAMPLES_DIR)/OnLoops/$(BUILD_DIR)/OnLoops: $(wildcard $(EXAMPLES_DIR)/OnLoops/*.cpp)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(EXAMPLES_DIR)/Ising/$(BUILD_DIR)/Ising: $(wildcard $(EXAMPLES_DIR)/Ising/*.cpp)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# ----- Specific build rules for tests -----
$(TESTS_DIR)/$(BUILD_DIR)/testVector: $(TESTS_DIR)/testVector.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testMpfr: $(TESTS_DIR)/testMpfr.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testMatrices: $(TESTS_DIR)/testMatrices.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testBigFloats: $(TESTS_DIR)/testBigFloats.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Run tests after building
run-tests: tests
	@set -e; \
	for test_exec in $(TESTS); do \
	    echo "Running $$test_exec"; \
	    $$test_exec; \
	done
	@echo "✅ All tests passed."

projects: loop3pt

loop3pt: loop3pt_222 loop3pt_123 loop3pt_323 loop3pt_general

loop3pt_222: $(PROJECTS_DIR)/loop3pt_222/$(BUILD_DIR)/loop3pt_222
loop3pt_123: $(PROJECTS_DIR)/loop3pt_123/$(BUILD_DIR)/loop3pt_123
loop3pt_323: $(PROJECTS_DIR)/loop3pt_323/$(BUILD_DIR)/loop3pt_323
loop3pt_general: $(PROJECTS_DIR)/loop3pt_general/$(BUILD_DIR)/loop3pt_general

$(PROJECTS_DIR)/loop3pt_222/$(BUILD_DIR)/loop3pt_222: $(PROJECTS_DIR)/loop3pt_222/On3pt.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run-loop3pt_222: $(PROJECTS_DIR)/loop3pt_222/$(BUILD_DIR)/loop3pt_222
	$^ $(ARGS)

$(PROJECTS_DIR)/loop3pt_123/$(BUILD_DIR)/loop3pt_123: $(PROJECTS_DIR)/loop3pt_123/On3pt.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run-loop3pt_123: $(PROJECTS_DIR)/loop3pt_123/$(BUILD_DIR)/loop3pt_123
	$^ $(ARGS)

$(PROJECTS_DIR)/loop3pt_323/$(BUILD_DIR)/loop3pt_323: $(PROJECTS_DIR)/loop3pt_323/On3pt.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run-loop3pt_323: $(PROJECTS_DIR)/loop3pt_323/$(BUILD_DIR)/loop3pt_323
	$^ $(ARGS)

$(PROJECTS_DIR)/loop3pt_general/$(BUILD_DIR)/loop3pt_general: $(PROJECTS_DIR)/loop3pt_general/On3pt.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run-loop3pt_general: $(PROJECTS_DIR)/loop3pt_general/$(BUILD_DIR)/loop3pt_general
	$^ $(ARGS)

# Clean up
clean:
	rm -rf $(EXAMPLES_DIR)/*/build/*
	rm -rf $(EXAMPLES_DIR)/*/build/*
	rm -rf $(TESTS_DIR)/build/*
	rm -rf $(TESTS_DIR)/build/*
	rm -rf $(PROJECTS_DIR)/*/build/*
	rm -rf $(PROJECTS_DIR)/*/build/*

# debug mode
debug:
	$(MAKE) BUILD_TYPE=debug

.PHONY: all examples tests projects run-tests clean

BUILD_TYPE ?= release

CXX := g++

ifeq ($(BUILD_TYPE),debug)
CXXFLAGS := -Wall -Wextra -O0 -ggdb
else
CXXFLAGS := -Wall -Wextra -O3 -march=native
endif
CXXFLAGS += -std=c++20 -I./include

LDFLAGS := -lgmp -lmpfr -lmpc
# The GNU libraries gmp, mpfr and mpc are required 
# to build the program. They can be installed using
# your system's package manager, or directly downloaded
# and build from source. If they are build from source and installed
# in a non-standard location, provide
# LDFLAGS += -L/path/to/gmp/lib -L/path/to/mpfr/lib -L/path/to/mpc/lib

TESTS_DIR    := tests

ifeq ($(BUILD_TYPE),release)
BUILD_DIR    := build/release
else
BUILD_DIR    := build/debug
endif

TESTS    := $(TESTS_DIR)/$(BUILD_DIR)/testVector \
            $(TESTS_DIR)/$(BUILD_DIR)/testMpfr \
            $(TESTS_DIR)/$(BUILD_DIR)/testBigFloats

all: tests projects

tests:    $(TESTS)

$(TESTS_DIR)/$(BUILD_DIR)/testVector: $(TESTS_DIR)/testVector.cpp 
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testMpfr: $(TESTS_DIR)/testMpfr.cpp 
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testBigFloats: $(TESTS_DIR)/testBigFloats.cpp 
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run-tests: tests
	@set -e; \
	for test_exec in $(TESTS); do \
	    echo "Running $$test_exec"; \
	    $$test_exec; \
	done
	@echo "✅ All tests passed."

projects: 3pt_spin

3pt_spin: $(BUILD_DIR)/3pt_spin

$(BUILD_DIR)/3pt_spin: On3pt.cpp 
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run-3pt_spin: $(BUILD_DIR)/3pt_spin
	$^ $(ARGS)

clean:
	rm -rf $(EXAMPLES_DIR)/*/build
	rm -rf $(TESTS_DIR)/build
	rm -rf $(PROJECTS_DIR)/*/build

distclean: clean
	rm -rf $(GMP_BUILD_DIR) $(MPFR_BUILD_DIR) $(MPC_BUILD_DIR)

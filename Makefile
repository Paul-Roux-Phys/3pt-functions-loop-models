# Compiler
CXX         := clang++
CXXFLAGS    := -Wall -Wextra -O3 -stdlib=libc++
CXXFLAGS    += -std=c++23
LDFLAGS     := 

# External libraries
CXXFLAGS    += -isystem $(shell nix eval --raw nixpkgs#llvmPackages.libcxx.dev)/include/c++/v1
SDK_PATH    := $(shell xcrun --show-sdk-path)
CXXFLAGS    += -isysroot $(SDK_PATH)
CXXFLAGS    += $(shell pkg-config --cflags gmp)
CXXFLAGS    += $(shell pkg-config --cflags mpfr)
LDFLAGS     += $(shell pkg-config --libs mpfr)

# Directories
EXAMPLES_DIR := examples
TESTS_DIR    := tests
BUILD_DIR    := build

# Targets
EXAMPLES := OnLoops Ising
TESTS := testVector testMpfr testMatrices

# Default target
all: examples

# Targets
EXAMPLES := $(EXAMPLES_DIR)/OnLoops/$(BUILD_DIR)/OnLoops \
            $(EXAMPLES_DIR)/Ising/$(BUILD_DIR)/Ising

TESTS    := $(TESTS_DIR)/$(BUILD_DIR)/testVector \
            $(TESTS_DIR)/$(BUILD_DIR)/testMpfr \
            $(TESTS_DIR)/$(BUILD_DIR)/testMatrices

# Default target
all: examples tests

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
$(TESTS_DIR)/$(BUILD_DIR)/testVector: $(TESTS_DIR)/Vector.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testMpfr: $(TESTS_DIR)/Mpfr.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testMatrices: $(TESTS_DIR)/Matrices.cpp
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

# Clean up
clean:
	rm -rf $(EXAMPLES_DIR)/*/$(BUILD_DIR)
	rm -rf $(TESTS_DIR)/$(BUILD_DIR)

.PHONY: all examples tests run-tests clean

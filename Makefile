# Compiler
CXX         := g++
CXXFLAGS    := -Wall -Wextra -O3

# Directories
EXAMPLES_DIR := examples
TEST_DIR     := tests
BUILD_DIR    := build

# Find sources
EXAMPLES_SOURCES := $(wildcard $(EXAMPLES_DIR)/*.cpp)
EXAMPLES_BINARIES := $(patsubst $(EXAMPLES_DIR)/%.cpp,$(EXAMPLES_DIR)/$(BUILD_DIR)/%,$(EXAMPLES_SOURCES))

TEST_SOURCES := $(wildcard $(TEST_DIR)/*.cpp)
TEST_BINARIES := $(patsubst $(TEST_DIR)/%.cpp,$(TEST_DIR)/$(BUILD_DIR)/%,$(TEST_SOURCES))

# Default target
all: examples test

# Build examples
examples: $(EXAMPLES_BINARIES)

$(EXAMPLES_DIR)/$(BUILD_DIR)/%: $(EXAMPLES_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $< -o $@

# Build tests
tests: $(TEST_BINARIES)

$(TEST_DIR)/$(BUILD_DIR)/%: $(TEST_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $< -o $@

# Run tests after building
run-tests: tests
	@set -e; \
	for test_exec in $(TEST_BINARIES); do \
	    echo "Running $$test_exec"; \
	    $$test_exec; \
	done
	@echo "✅ All tests passed."

# Clean up
clean:
	rm -rf $(EXAMPLES_DIR)/$(BUILD_DIR)
	rm -rf $(TEST_DIR)/$(BUILD_DIR)

.PHONY: all examples tests run-tests clean

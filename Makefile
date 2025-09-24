# =================================================
# Build type
# =================================================
BUILD_TYPE ?= release

CXX := g++

ifeq ($(BUILD_TYPE),debug)
CXXFLAGS := -Wall -Wextra -O0 -ggdb
else
CXXFLAGS := -Wall -Wextra -O3 -march=native
endif
CXXFLAGS += -std=c++20 -I./include

LDFLAGS := 

# =================================================
# Dependencies
# =================================================
# Versions
GMP_VERSION  := 6.3.0
MPFR_VERSION := 4.2.2
MPC_VERSION  := 1.3.1

# Directories
THIRDPARTY     := $(abspath thirdparty)
GMP_SRC_DIR    := $(THIRDPARTY)/gmp-$(GMP_VERSION)
GMP_BUILD_DIR  := $(THIRDPARTY)/build/gmp
MPFR_SRC_DIR   := $(THIRDPARTY)/mpfr-$(MPFR_VERSION)
MPFR_BUILD_DIR := $(THIRDPARTY)/build/mpfr
MPC_SRC_DIR    := $(THIRDPARTY)/mpc-$(MPC_VERSION)
MPC_BUILD_DIR  := $(THIRDPARTY)/build/mpc

# System detection
HAVE_GMP := $(shell command -v pkg-config >/dev/null 2>&1 && pkg-config --exists gmp && echo yes || echo no)
HAVE_MPFR := $(shell command -v pkg-config >/dev/null 2>&1 && pkg-config --exists mpfr && echo yes || echo no)
HAVE_MPC := $(shell command -v pkg-config >/dev/null 2>&1 && pkg-config --exists mpc && echo yes || echo no)

# Dependency selection
ifeq ($(HAVE_GMP),yes)
  CXXFLAGS += $(shell pkg-config --cflags gmp)
  LDFLAGS  += $(shell pkg-config --libs gmp)
else
  CXXFLAGS += -I$(THIRDPARTY)/include
  LDFLAGS  += -L$(THIRDPARTY)/lib -lgmp
  NEED_GMP := $(THIRDPARTY)/lib/libgmp.a
endif

ifeq ($(HAVE_MPFR),yes)
  CXXFLAGS += $(shell pkg-config --cflags mpfr)
  LDFLAGS  += $(shell pkg-config --libs mpfr)
else
  CXXFLAGS += -I$(THIRDPARTY)/include
  LDFLAGS  += -L$(THIRDPARTY)/lib -lmpfr
  NEED_MPFR := $(THIRDPARTY)/lib/libmpfr.a
endif

ifeq ($(HAVE_MPC),yes)
  CXXFLAGS += $(shell pkg-config --cflags mpc)
  LDFLAGS  += $(shell pkg-config --libs mpc)
else
  CXXFLAGS += -I$(THIRDPARTY)/include
  LDFLAGS  += -L$(THIRDPARTY)/lib -lmpc
  NEED_MPC := $(THIRDPARTY)/lib/libmpc.a
endif

NEED_DEPS := $(NEED_GMP) $(NEED_MPFR) $(NEED_MPC)

# =================================================
# Project directories
# =================================================
EXAMPLES_DIR := examples
TESTS_DIR    := tests
PROJECTS_DIR := projects

ifeq ($(BUILD_TYPE),release)
BUILD_DIR    := build/release
else
BUILD_DIR    := build/debug
endif

# =================================================
# Targets
# =================================================
EXAMPLES := $(EXAMPLES_DIR)/OnLoops/$(BUILD_DIR)/OnLoops 
TESTS    := $(TESTS_DIR)/$(BUILD_DIR)/testVector \
            $(TESTS_DIR)/$(BUILD_DIR)/testMpfr \
            $(TESTS_DIR)/$(BUILD_DIR)/testMatrices \
            $(TESTS_DIR)/$(BUILD_DIR)/testBigFloats

all: examples tests projects

examples: $(EXAMPLES)
tests:    $(TESTS)

# =================================================
# Example build rules
# =================================================
$(EXAMPLES_DIR)/OnLoops/$(BUILD_DIR)/OnLoops: $(wildcard $(EXAMPLES_DIR)/OnLoops/*.cpp) $(NEED_DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(EXAMPLES_DIR)/Ising/$(BUILD_DIR)/Ising: $(wildcard $(EXAMPLES_DIR)/Ising/*.cpp) $(NEED_DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# =================================================
# Test build rules
# =================================================
$(TESTS_DIR)/$(BUILD_DIR)/testVector: $(TESTS_DIR)/testVector.cpp $(NEED_DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testMpfr: $(TESTS_DIR)/testMpfr.cpp $(NEED_DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testMatrices: $(TESTS_DIR)/testMatrices.cpp $(NEED_DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TESTS_DIR)/$(BUILD_DIR)/testBigFloats: $(TESTS_DIR)/testBigFloats.cpp $(NEED_DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run-tests: tests
	@set -e; \
	for test_exec in $(TESTS); do \
	    echo "Running $$test_exec"; \
	    $$test_exec; \
	done
	@echo "✅ All tests passed."

# =================================================
# Projects
# =================================================
projects: loop3pt

loop3pt: $(PROJECTS_DIR)/loop3pt/$(BUILD_DIR)/loop3pt

$(PROJECTS_DIR)/loop3pt/$(BUILD_DIR)/loop3pt: $(PROJECTS_DIR)/loop3pt/On3pt.cpp $(NEED_DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run-loop3pt: $(PROJECTS_DIR)/loop3pt/$(BUILD_DIR)/loop3pt
	$^ $(ARGS)

# =================================================
# Vendored builds
# =================================================
$(THIRDPARTY)/lib/libgmp.a: $(GMP_SRC_DIR)/configure
	mkdir -p $(GMP_BUILD_DIR)
	cd $(GMP_BUILD_DIR) && \
		$(abspath $(GMP_SRC_DIR))/configure \
			--prefix=$(THIRDPARTY) && \
		$(MAKE) && $(MAKE) install

$(THIRDPARTY)/lib/libmpfr.a: $(MPFR_SRC_DIR)/configure $(NEED_GMP)
	mkdir -p $(MPFR_BUILD_DIR)
	cd $(MPFR_BUILD_DIR) && \
		$(abspath $(MPFR_SRC_DIR))/configure \
			--prefix=$(THIRDPARTY) --with-gmp=$(THIRDPARTY) && \
		$(MAKE) && $(MAKE) install

$(THIRDPARTY)/lib/libmpc.a: $(MPC_SRC_DIR)/configure $(NEED_MPFR)
	mkdir -p $(MPC_BUILD_DIR)
	cd $(MPC_BUILD_DIR) && \
		$(abspath $(MPC_SRC_DIR))/configure \
			--prefix=$(THIRDPARTY) --with-gmp=$(THIRDPARTY) --with-mpfr=$(THIRDPARTY) && \
		$(MAKE) && $(MAKE) install

# =================================================
# Housekeeping
# =================================================
clean:
	rm -rf $(EXAMPLES_DIR)/*/build
	rm -rf $(TESTS_DIR)/build
	rm -rf $(PROJECTS_DIR)/*/build

distclean: clean
	rm -rf $(GMP_BUILD_DIR) $(MPFR_BUILD_DIR) $(MPC_BUILD_DIR)

FILENAME := On_loops

CPPC := g++
FLAGS := -std=c++20 -O3
WARNINGS := -Wall -Wno-overloaded-virtual

INCLUDES := -I/usr/local/opt/arpackpp/include \
			-I/usr/local/opt/transfer_matrices/include
LIBS := -lmpfr -lmpc -larpack -lblas

SRC := $(FILENAME).cpp
BIN_DIR := bin

ifdef bin
BIN := $(BIN_DIR)/$(bin)
else
BIN := $(BIN_DIR)/$(FILENAME)
endif

ifdef size
FLAGS += -DLATTICE_SIZE=$(size)
endif

ifdef hash_size
FLAGS += -DHASH_SIZE=$(hash_size)
endif

all: clean $(BIN)

$(BIN): | $(BIN_DIR)
	$(CPPC) $(FLAGS) $(WARNINGS) $(INCLUDES) $(SRC) -o $(BIN) $(LIBS)

clean:
	rm -f $(BIN)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
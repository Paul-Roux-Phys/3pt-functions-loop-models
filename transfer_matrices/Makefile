CPPC := g++
FLAGS := -std=c++20 -O3
WARNINGS := -Wall -Wno-overloaded-virtual

INCLUDES := -I/usr/local/opt/arpackpp/include\ 
			-I/Users/Paul/Documents/Recherche/projet_these/TransferMatricesCpp/include
LIBS := -lmpfr -lmpc -larpack -lblas

SRC := $(FILENAME).cpp

all: $(OBJ)

ifdef hash_size
FLAGS += -DHASH_SIZE=$(hash_size)
endif

$(BIN): | $(BIN_DIR)
	$(CPPC) $(FLAGS) $(WARNINGS) $(INCLUDES) $(SRC) -o $(BIN) $(LIBS)

clean:
	rm -f $(BIN)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
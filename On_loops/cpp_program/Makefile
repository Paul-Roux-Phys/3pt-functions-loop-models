BIN := bin/transfer
SRC = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp, obj/%.o, $(SRC))

INCLUDES := -I../../libs/libtransfermatrices/include

all: $(BIN)

$(BIN): $(OBJ)
	@mkdir -p bin
	g++ $(OBJ) -o $(BIN) -larpack -lblas

obj/%.o: %.cpp
	@mkdir -p obj
	g++ -std=c++20 -O3 $(INCLUDES) -c $< -o $@
	@echo "Compiled $< -> $@"

clean:
	rm -rf bin obj
	@echo "Cleaned bin and obj"

CXX = g++
CXXFLAGS = -O3 -Wall -std=c++17
LIBS = 
SRC = testRead.cc
OBJ = $(addprefix build/,$(SRC:.cc=.o))
TARGET = build/readFile

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(TARGET)

build/%.o: %.cc | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build:
	mkdir -p build

clean:
	@echo 'Cleaning up...'
	@rm -rf build

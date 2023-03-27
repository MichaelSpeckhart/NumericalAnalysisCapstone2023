# CXX: This variable specifies the C++ compiler to use. In this case, it is set to g++.
# CXXFLAGS: This variable specifies the compiler flags to use, including the optimization level (-O3), warning level (-Wall), and Pybind11 include directories.
# LIBS: This variable specifies any additional libraries to link against, in this case, the libraries needed for linking against Python.
# SRC: This variable specifies the source code file to compile.
# OBJ: This variable specifies the object file that will be generated from the source code file.
# TARGET: This variable specifies the name of the final output file, including the extension suffix for the version of Python being used.
# all: This is a phony target that specifies the default target to build when running make.
# $(TARGET): This target specifies the rules for building the target file, including compiling the object files and linking them into the final shared library.
# $(OBJ): This target specifies the rule for compiling the source code into object files.
# clean: This target specifies a rule for cleaning up the build directory, removing all object files and the final output file.

CXX = g++
CXXFLAGS = -O3 -Wall -shared -std=c++17 -fPIC
SRC = matrixMarket mmio
OBJ = $(SRC:.c=.o)
TARGET = build

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LIBS) $(OBJ) -o $(TARGET)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@echo 'Cleaning up...'
	@rm -f $(OBJ) $(TARGET)

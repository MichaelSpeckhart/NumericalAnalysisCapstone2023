CXX = g++
CXXFLAGS = -O3 -Wall -shared -Werror -fopenmp -std=c++17 -fPIC
LIBS = -lgomp
SRC = functions.cc functionsCSC.cc functionsCSR.cc functionsCOO.cc
OBJ = $(SRC:.cc=.o)
TARGET = ../../build/library.so
DEST = ../../build/

all: $(TARGET)

$(TARGET): $(OBJ)
	@echo 'Linking...'
	@echo $(OBJ) '-->' $(TARGET)
	@$(CXX) $(CXXFLAGS) -o $(TARGET) $(addprefix $(DEST), $(OBJ)) $(LIBS)

%.o: %.cc
	@echo $< '-->' $@
	@$(CXX) $(CXXFLAGS) -c $< -o $(DEST)$@

clean:
	@echo 'Cleaning up...'
	@rm -f $(OBJ) $(TARGET)
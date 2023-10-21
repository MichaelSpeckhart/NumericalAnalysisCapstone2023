CXX = arm-linux-gnueabihf-g++ -march=armv7-a -mthumb -mthumb-interwork -mfloat-abi=hard -mfpu=neon-vfpv4 -mtls-dialect=gnu  -march=armv7-a  -mthumb -mfloat-abi=hard -mfpu=neon -mvectorize-with-neon-quad 
CXXFLAGS = -O3 -Wall -shared -Werror -fopenmp -std=c++17 -fPIC
LIBS = -lgomp
SRC = functions.cc functionsCSC.cc functionsCSR.cc functionsCOO.cc
OBJ = $(SRC:.cc=.o)
TARGET = ../../build/library.so
DEST = ../../build/
SERVER_DST = ../server/

all: $(TARGET)

$(TARGET): $(OBJ)
	@echo 'Linking...'
	@echo $(OBJ) '-->' $(TARGET)
	@$(CXX) $(CXXFLAGS) -o $(TARGET) $(addprefix $(DEST), $(OBJ)) $(LIBS)
	@cp $(TARGET) $(SERVER_DST)

%.o: %.cc
	@echo $< '-->' $@
	@$(CXX) $(CXXFLAGS) -c $< -o $(DEST)$@

clean:
	@echo 'Cleaning up...'
	@rm -f $(OBJ) $(TARGET)

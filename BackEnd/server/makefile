# This defaults bits to 64, but allows it to be overridden on the command
# line
BITS = 64

# Output directory
ODIR  = build
tmp  := $(shell mkdir -p $(ODIR))

# Basic compiler configuration and flags
CXX      = g++
CXXFLAGS = -MMD -O3 -Wall -shared -std=c++17 -fopenmp -fPIC 
LDFLAGS	 = -m$(BITS) -lpthread -lrt -lcrypt -ldl -lutil -lm -lm 

# The basenames of the c++ files that this program uses
CXXFILES = http_server parse_requests server_main

# The executable we will build
TARGET = $(ODIR)/CapstoneServer.exe

# Create the .o names from the CXXFILES
OFILES = $(patsubst %, $(ODIR)/%.o, $(CXXFILES))

# Create .d files to store dependency information, so that we don't need to
# clean every time before running make
DFILES = $(patsubst %.o, %.d, $(OFILES))

# Default rule builds the executable
all: $(TARGET)

# clean up everything by clobbering the output folder
clean:
	@echo cleaning up...
	@rm -rf $(ODIR)

# build an .o file from a .cc file
$(ODIR)/%.o: src/%.cc
	@echo [CXX] $< "-->" $@
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

# Link rule for building the target from .o files
$(TARGET): $(OFILES)
	@echo [LD] $^ "-->" $@
	@$(CXX) -o $@ $^ $(LDFLAGS)

# Remember that 'all' and 'clean' aren't real targets
.PHONY: all clean

# Pull in all dependencies
-include $(DFILES)
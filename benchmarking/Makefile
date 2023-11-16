# This defaults bits to 64, but allows it to be overridden on the command
# line
BITS = 64

# Output directory
ODIR  = buildUnitTest
tmp  := $(shell mkdir -p $(ODIR))

# Basic compiler configuration and flags
CXX      = g++
CXXFLAGS = -MMD -ggdb -O3 -std=gnu++17 -m$(BITS) -mavx
LDFLAGS	 = -m$(BITS) #-lpthread -lrt -L../CSE475/oneapi-tbb-2021.8.0/lib/intel64/gcc4.8 -ltbb

# The basenames of the c++ files that this program uses
CXXFILES = functionsCSR

# The executable we will build
TARGET = $(ODIR)/unitTests

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
$(ODIR)/%.o: %.cc
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
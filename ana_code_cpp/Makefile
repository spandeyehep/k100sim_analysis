#CXX     = g++
CXX     = clang++

CXXFLAGS= $(shell root-config --cflags)
LIBS    = $(shell root-config --libs) 


SOURCES = NtupleVariables.cc Analyze_k100sim.cc
HEADERS = NtupleVariables.h Analyze_k100sim.h
OBJECTS = $(SOURCES:.cc=.o)

EXECUTABLE = k100sim_ana

all: $(SOURCES) $(EXECUTABLE)

%.o: %.cc
	@echo Compiling $<...
	$(CXX) $(CXXFLAGS) -c -o $@ $< 


$(EXECUTABLE): $(OBJECTS)
	@echo "Linking $(PROGRAM) ..."
	@echo "@$(CXX) $(LIBS) $(OBJECTS) -o $@"
	@$(CXX) -o $@ $^ $(LIBS)
	@echo "done"


# Specifying the object files as intermediates deletes them automatically after the build process.
.INTERMEDIATE:  $(OBJECTS)

# The default target, which gives instructions, can be called regardless of whether or not files need to be updated.
.PHONY : clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

###
NtupleVariables.o: ./NtupleVariables.h
Analyze_k100sim.o: ./NtupleVariables.h ./Analyze_k100sim.h

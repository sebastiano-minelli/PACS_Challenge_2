### change PACS_EXAMPLES_PATH pasting your machine path to 'pacs-examples/' git folder
PACS_EXAMPLES_PATH=$(shell echo $$PACS_ROOT)
EIGEN_INC = $(shell echo $$mkEigenInc)

CXX      ?= mpic++
CXXFLAGS ?= -std=c++20 -fopenmp

EXEC = main
OBJECT = main.cpp
COMP = main.o
LDINC = -I$(PACS_EXAMPLES_PATH)/Examples/include/muparserx -I$(EIGEN_INC) -I$(PACS_EXAMPLES_PATH)/Examples/src/Parallel/Utilities -I../lib -I.
LDLIB = $(PACS_EXAMPLES_PATH)/Examples/lib/libmuparserx.a
LDNAME = pacs

.PHONY = clean help debug standard optimized

.DEFAULT_GOAL = optimized

.SECONDARY: $(EXEC) $(COMP)

# set compilation mode
debug: CPPFLAGS = -g -Wall -fopenmp -Wno-conversion-null -Wno-deprecated-declarations
debug: all

standard: CPPFLAGS = -O0 -Wall -fopenmp -Wno-conversion-null -Wno-deprecated-declarations
standard: all

optimized: CPPFLAGS = -O3 -Wall -fopenmp -Wno-conversion-null -Wno-deprecated-declarations
optimized: all

all: clean $(EXEC)

$(EXEC): $(COMP)
	$(CXX) $(CXXFLAGS) $(COMP) $(LDLIB) -o $(EXEC)
	
$(COMP): $(OBJECT)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDINC) -c $(OBJECT)
	
clean:
	$(RM) *.o $(EXEC)

help:
	@echo "Usage: make [TARGET]"
	@echo "TARGETS:"
	@echo " - debug:        compile in debug mode"
	@echo " - standard:     compile in standard mode, i.e. no optimization"
	@echo " - optimized:    compile in -03 optimization mode"
	@echo " - clean:        remove all build artifacts"

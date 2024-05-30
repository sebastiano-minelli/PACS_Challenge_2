### change PACS_EXAMPLES_PATH pasting your machine path to 'pacs-examples/' git folder
PACS_EXAMPLES_PATH=$(shell echo $$PACS_ROOT)

CXX      ?= g++
CXXFLAGS ?= -std=c++20

EXEC = main
OBJECT = main.cpp
COMP = main.o
LDINC = $(PACS_EXAMPLES_PATH)/Examples/include/
LDLIB = $(PACS_EXAMPLES_PATH)/Examples/lib/
LDNAME = pacs

.PHONY = clean help debug standard optimized

.DEFAULT_GOAL = optimized

.SECONDARY: $(EXEC) $(COMP)

# set compilation mode
debug: CPPFLAGS = -g -I. -I../lib -Wall -Wno-conversion-null -Wno-deprecated-declarations
debug: all

standard: CPPFLAGS = -O0 -I. -I../lib -Wall -Wno-conversion-null -Wno-deprecated-declarations
standard: all

optimized: CPPFLAGS = -O3 -I. -I../lib -Wall -Wno-conversion-null -Wno-deprecated-declarations
optimized: all

all: clean $(EXEC)

$(EXEC): $(COMP)
	$(CXX) $(CXXFLAGS) -L$(LDLIB) -Wl,-rpath,$(LDLIB) -l$(LDNAME) $(COMP) -o $(EXEC)
	
$(COMP): $(OBJECT)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(LDINC) -c $(OBJECT)
	
clean:
	$(RM) *.o $(EXEC)

help:
	@echo "Usage: make [TARGET]"
	@echo "TARGETS:"
	@echo " - debug:        compile in debug mode"
	@echo " - standard:     compile in standard mode, i.e. no optimization"
	@echo " - optimized:    compile in -03 optimization mode"
	@echo " - clean:        remove all build artifacts"

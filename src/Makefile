### change PACS_EXAMPLES_PATH pasting your machine path to 'pacs-examples/' git folder
PACS_EXAMPLES_PATH=${PACS_ROOT}

CXX      ?= g++
CXXFLAGS ?= -std=c++20

EXEC = main
OBJECT = main.cpp
COMP = main.o

.PHONY = clean help debug standard optimized

.DEFAULT_GOAL = standard

.SECONDARY: $(EXEC) $(COMP)

# set compilation mode
debug: CPPFLAGS = -g -Wall -I. -I../lib -Wno-conversion-null -Wno-deprecated-declarations
debug: $(EXEC)

standard: CPPFLAGS = -O0 -Wall -I. -I../lib -Wno-conversion-null -Wno-deprecated-declarations
standard: $(EXEC)

optimized: CPPFLAGS = -O3 -Wall -I. -I../lib -Wno-conversion-null -Wno-deprecated-declarations
optimized: $(EXEC)

$(EXEC): $(COMP)
	$(CXX) $(CXXFLAGS) $(COMP) -o $(EXEC)
	
$(COMP): $(OBJECT)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(OBJECT)
	
clean:
	$(RM) *.o

help:
	@echo "Usage: make [TARGET]"
	@echo "TARGETS:"
	@echo "  debug    compile in debug mode"
	@echo "  standard compile in standard mode, i.e. no optimization"
	@echo "  optimized  compile in -03 optimization mode"
	@echo "  clean    remove all build artifacts"
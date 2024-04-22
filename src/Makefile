### change PACS_EXAMPLES_PATH pasting your machine path to 'pacs-examples/' git folder
PACS_EXAMPLES_PATH=${PACS_ROOT}

CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O0 -Wall -I. -I../lib -Wno-conversion-null -Wno-deprecated-declarations

EXEC = main
OBJECT = main.cpp
COMP = main.o

.PHONY = clean distclean

EXECUTE: COMPILE
	$(CXX) $(CXXFLAGS) $(COMP) -o $(EXEC)
	
COMPILE:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(OBJECT) -c
clean:
	$(RM) *.o
	
distclean: clean
	$(RM) *~

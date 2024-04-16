### change PACS_EXAMPLES_PATH pasting your machine path to 'pacs-examples/' git folder
PACS_EXAMPLES_PATH=../../../../pacs-examples

CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -I. -I../lib -Wno-conversion-null -Wno-deprecated-declarations -I$(PACS_EXAMPLES_PATH)/Examples/include/muparserx

EXEC = main
OBJECT = main.cpp
STATIC_LIB = $(PACS_EXAMPLES_PATH)/Examples/lib/libmuparserx.a
COMP = main.o

.PHONY = clean distclean

EXECUTE: COMPILE
	$(CXX) $(CXXFLAGS) $(COMP) $(STATIC_LIB) -o $(EXEC)
	
COMPILE:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(OBJECT) -c
clean:
	$(RM) *.o
	
distclean: clean
	$(RM) *~

all: implicitize-test

DC=../dual-contouring

INCLUDES=-I$(DC)
LIBS=-L$(DC)/build -ldualcontour

CXXFLAGS=-O3 $(INCLUDES)

implicitize-test: implicitize-test.o surf21.o surf31.o
	g++ -o $@ $^ $(LIBS)

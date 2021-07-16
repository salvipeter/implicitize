all: implicitize-test

CXXFLAGS=-O3

implicitize-test: implicitize-test.cc
	g++ -o $@ $^ $(CXXFLAGS) \
            -I../dual-contouring \
            -L../dual-contouring/build -ldualcontour

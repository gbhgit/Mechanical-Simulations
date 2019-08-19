CC = g++
CFLAGS = -g -Wall
SRCS  += src/main.cpp
PROG   = peaucellier.o

CFLAGS_CPP = -std=c++11

OPENCV = `pkg-config --cflags --libs opencv`
LIBS = $(OPENCV)

$(PROG):$(SRCS)
	$(CC) $(CFLAGS_CPP) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)

clean:
	rm -rf $(PROG)

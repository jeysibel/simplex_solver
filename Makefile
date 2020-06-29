CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-g $(shell root-config --cflags)
LDFLAGS=-g $(shell root-config --ldflags)
LDLIBS=$(shell root-config --libs)

SRCS=tool.cc support.cc
OBJS=$(subst .cc,.o,$(SRCS))

all: tool

simplex_solver: $(OBJS)
    $(CXX) $(LDFLAGS) -o tool $(OBJS) $(LDLIBS)

simplex_solver.o: tool.cc support.hh

solver.o: solver.hh solver.cc

clean:
    $(RM) $(OBJS)

distclean: clean
    $(RM) simplex_solver
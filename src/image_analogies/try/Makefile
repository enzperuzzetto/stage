ROOT=../../bcl

CC = gcc
CPPFLAGS = -I. -I$(ROOT)/include 
CFLAGS   = -Wall -Wextra -Werror -pedantic -std=c99 -g
LDFLAGS  = -L$(ROOT)/lib 
LDLIBS   = -lbcl -lm

BIN=\
	analogy\
	#synthesis

.PHONY: all
all: $(BIN)

synthesis: conversion.o features.o synthesis.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

analogy : features.o pyramid.o conversion.o BestMatch.o main.o  #	patchMatch.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
%.o:%.c
	$(CC) $(CFLAGS) -c  $< -I../../bcl/include

.PHONY: clean cleanall
clean:
	$(RM) *.o *~ *.ppm
cleanall: clean
	$(RM) $(BIN) 

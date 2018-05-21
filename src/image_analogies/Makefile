ROOT=../bcl

CC = gcc
CPPFLAGS = -I. -I$(ROOT)/include 
CFLAGS   = -Wall -Wextra -Werror -pedantic -std=c99
LDFLAGS  = -L$(ROOT)/lib 
LDLIBS   = -lbcl -lm

BIN=\
	analogy

analogy : features.o pyramid.o conversion.o patchMatch.o BestMatch.o main.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
%.o:%.c
	$(CC) $(CFLAGS) -c  $< -I../bcl/include

.PHONY: clean cleanall
clean:
	$(RM) *.o *~ *.ppm
cleanall: clean
	$(RM) $(BIN) 
